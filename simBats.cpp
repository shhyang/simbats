/*
 Copyright (C) 2014 Shenghao Yang
 
 This file is part of SimBATS.
 
 SimBATS is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 SimBATS is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with SimBATS.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <sys/times.h>

#include "Utilities.h"
#include "BatchEnc.h"
#include "BatchDec.h"
//#include "NCCoder.h"
#include "FiniteField.h"

#include "simRecoder.h"

using namespace std;

// Parameters
// code parameters
#define BATCH_SIZE 16 //M
#define PACKET_LEN 2// T
#define GFORDER 8   //m
#define PACKET_NO 64000 //K
#define DESIGN_DEC_RATIO 0.96 //1-eta
#define DEC_RATIO DESIGN_DEC_RATIO // for control only
#define ERASURE_R 0.2   //e
#define ITERATION_NO 1000
#define CHANEL_TYPE 1 // 0 tri 1 tri one 2 sq (Change channel in code!!!)


class DecoderStatus{
public:
    int nTrans;
    int nReceive;
    int nSave;
    int nError;
    int nInact;
    double * rankdist;

    DecoderStatus(int M){
        rankdist = new double[M+1];
    }
    ~DecoderStatus(){
        delete [] rankdist;
    }
};

class TimeUsed{
    private:
        struct Interval{
            private:
                struct tms start_t, end_t;
                long s, e;
                long t;
            public:
                void clear(){
                    t = 0;
                }

                void start(){
                    s = times(&start_t);
                }

                void end(){
                    e = times(&end_t);
                    t += e - s;
                }

                double time(){
                    return (double) t/60; //sysconf(_SC_CLK_TCK);
                }
        };

    public:
        Interval encoding;
        Interval decoding;
        Interval nccoding;
    public:
        void clear(){
            encoding.clear();
            decoding.clear();
            nccoding.clear();
        }
};

class BATSimulator {
private:
    int batchSize;
    int gf_order;
    int packetNum;
    int packetSize;

    int nEncLimit;
    int packetSizeWithHeader;
    int packetSizeWithHeaderAndId;

    double degreeDist[MAX_DEGREE];
    int D;

	double* rankDist;

    MTRand *psrand;
    SymbolType* input;
    SymbolType* output;

    BatsEncoder *encoder;
    BatsDecoder *decoder;
    //NCCoder *nccoder;
    SimRecoder *simcoder;

public:
    BATSimulator(int M, int m, int K, int T) : batchSize(M), gf_order(m), packetNum(K), packetSize(T) {

        nEncLimit = packetNum / batchSize * 30;
        packetSizeWithHeader = packetSize + batchSize * gf_order / 8;
        packetSizeWithHeaderAndId = packetSizeWithHeader + sizeof(KeyType);

        // Get degree distribution according to the channel
        getDegDist();

		// Get rank distribution
		getRankDist();

		//ff_init();//Tom: Added on interface change
        FF.setOrder(gf_order);

        psrand = new MTRand(); // 5344 no // 3544 all

        input = new SymbolType[packetNum * packetSize];

        output = new SymbolType[packetNum * packetSize];

        for (int i = 0; i < packetNum; i++) {
            for (int j = 0; j < packetSize; j++) {
                input[j + i * packetSize] = psrand->randInt(255);
                output[j + i * packetSize] = 255;
            }
        }

        encoder = NULL;
        decoder = NULL;
        //nccoder = NULL;
        simcoder = NULL;
    }

    ~BATSimulator() {
        delete psrand;
        delete [] input;
        delete [] output;
        delete [] rankDist;
        clearNetwork();
    }

private:
    void getDegDist() {
        stringstream iss;

        iss << "simDegreeK" << packetNum << "M" << batchSize << "m" << gf_order << ".txt";
        //    iss << "ddM" << batchSize << "m" << gf_order << "FR.txt";
        ifstream filestr;

        filestr.open(iss.str().c_str());
        if (!filestr) {
            cout << "cannot get degree distribution. File name:" << iss.str().c_str() << "\n";
            return;
        }

        // the degree file start with degree 1.
        // degree 0 has probability 0.
        degreeDist[0] = 0.0;

        double x;
        D = 1;
        while (filestr >> x && D < MAX_DEGREE) {
            degreeDist[D] = x;
            D++;
        }

        for (int i = D; i < MAX_DEGREE; i++) {
            degreeDist[i] = 0;
        }

        filestr.close();
    }

	void getRankDist() {
        rankDist = new double[batchSize+1];

		for (int i = 0; i < batchSize + 1; i++)
			rankDist[i] = 0;

        stringstream iss;

        iss << "simRankDistM" << batchSize << "m" << gf_order << ".txt";

		ifstream filestr;
		filestr.open(iss.str().c_str());
		if (!filestr) {
			cout << "cannot get simulation rank distribution. File name:" << iss.str().c_str() << endl;
			return;
		}

		double x, sum;
		sum = 0;
		for (int i = 0; i < batchSize + 1; i++) {
			if (filestr >> x) {
				rankDist[i] = x;
				sum += x;
			} else {
				cout << "Warning: file ended with only " << i << " entries out of " << batchSize + 1 << " needed." << endl;
				filestr.close();
				return;
			}
		}
		filestr.close();

		if (sum < 1.0)
			cout << "Warning: Rank distribution does not sum to 1 (" << sum << ")" << endl;
	}

    void clearNetwork(){
        if (encoder != NULL)
            delete encoder;
        if (decoder != NULL)
            delete decoder;
        //if (nccoder != NULL)
        //    delete nccoder;
        if (simcoder != NULL)
			delete simcoder;
    }

    void initNetwork() {
        // Initialize encoder
        clearNetwork();

        encoder = new BatsEncoder(batchSize, input, packetNum, packetSize);
        //encoder = new BatsEncoder(batchSize, packetNum, packetSize, input);//Tom: interface change
		encoder->setDegreeDist(degreeDist, D);

        //bool ch = encoder.verifyCheckPkg();

        // Initialize decoder
        decoder = new BatsDecoder(batchSize, output, packetNum, packetSize);
		//decoder = new BatsDecoder(batchSize, packetNum, packetSize, output);//Tom: interface change
        decoder->setDegreeDist(degreeDist, D);

        // Initialize nccoder
        //nccoder = new NCCoder(batchSize, packetSizeWithHeader);
		//nccoder = new NCCoder(batchSize, packetSize);

		// Initialize simcoder
		simcoder = new SimRecoder(batchSize, packetSizeWithHeader, (1 << gf_order) - 1, rankDist, psrand);
    }

public:
    void runOnce(TimeUsed& timeUsed, DecoderStatus& ds) {
        initNetwork();

        run(timeUsed, ds);

        decoder->rankDist(ds.rankdist);

        ds.nError = 0;
        for (int i = 0; i < packetNum; i++) {
            for (int j = 0; j < packetSize; j++) {
                if (input[j + i * packetSize] != output[j + i * packetSize]) {
                    ds.nError++;
                    break;
                }
            }
        }
    }
    void run(TimeUsed& timeUsed, DecoderStatus& ds) {

        SymbolType** batch = mallocMat<SymbolType>(batchSize, packetSizeWithHeaderAndId);
        SymbolType** batchWithoutId = (SymbolType**) malloc(batchSize * sizeof (SymbolType*));

		//Recoding batch
        SymbolType** recBatch = mallocMat<SymbolType>(batchSize, packetSizeWithHeaderAndId);
        SymbolType** recBatchWithoutId = (SymbolType**) malloc(batchSize * sizeof (SymbolType*));

        for (int i = 0; i < batchSize; i++) {
            batchWithoutId[i] = batch[i] + sizeof(KeyType);
			recBatchWithoutId[i] = recBatch[i] + sizeof(KeyType);
        }

        //SymbolType* relayOut = new SymbolType[packetSizeWithHeaderAndId];

        KeyType id;

        int nEnc = 0;
        int encEnd = false;

		int recCnt;//Number of recoded packets in current batch

        timeUsed.clear();

        while (nEnc < nEncLimit && !encEnd) {

            timeUsed.encoding.start();

            id = encoder->genBatch(batchWithoutId);

            for (int i = 0; i < batchSize; i++) {
                saveIDInPacket(batch[i], &id);
            }

            timeUsed.encoding.end();

			/* Simulation of recoding */
			recCnt = simcoder->genBatch(recBatchWithoutId, batchWithoutId);
			for (int i = 0; i < batchSize; i++) {
				saveIDInPacket(recBatch[i], &id);
			}

            // transmit output in batchSize separated packets
            //KeyType relayId;

			/* Update method for decoding */
            timeUsed.decoding.start();

            for (int i = 0; i < batchSize; i++)
            {
                //Tom: Changed interface for recieving packets and inact decoding
                //decoder->receivePacket(recBatch[i]);
                decoder->receivePacket(recBatch[i] + sizeof (KeyType), getIDFromPacket(recBatch[i]));

                if (decoder->complete(1.0))
                {
                    encEnd = true;
                    break;
                }

            }

            timeUsed.decoding.end();


            nEnc++;
        }

        ds.nTrans = nEnc * batchSize;
        ds.nReceive = decoder->nRecPkg;
        ds.nSave = decoder->nSavedPkg;
        ds.nInact = decoder->nInactVar;

        freeMat(batch, batchSize);
        free(batchWithoutId);

		freeMat(recBatch, batchSize);
		free(recBatchWithoutId);
        //delete [] relayOut;
		//delete [] recodeBatch;
    }
};

int main(int argc, char* argv[]) {

    // set parameters
    int batchSize;
    int gf_order = 8; // 1, 2, 4, 8
    int packetNum;
    int packetSize = 1; //In Bytes
    //int packetSizeInSymbol = packetSize * SymbolSize / gf_order;

    int iterationNum;
    //const float decRatio = 0.99;


    switch(argc) {
        case 1:
            batchSize = 32; // 16, 32, 64
            packetNum = 1600;
            iterationNum = 10;
            break;
        case 4:
            batchSize = atoi(argv[1]);
            packetNum = atoi(argv[2]);
            iterationNum = atoi(argv[3]);
            break;
        default:
            cout << "simbats M K numIteration" << endl;
            return 0;
    }

    BATSimulator sim(batchSize, gf_order, packetNum, packetSize);

    cout << "Simulation starts with " << "M = " << batchSize << ", q = 2^" << gf_order << ", K = " << packetNum << ", T = " << packetSize << endl;

    int iter = 0;

    int errIdx[iterationNum];
    int nSucc = 0;
    int nErrIdx = 0;
    float totalTrans = 0.0;
    TimeUsed timeUsed;
    DecoderStatus ds(batchSize);
    double accuRankDist[batchSize + 1];

    for (int i = 0; i < batchSize + 1; i++) {
        accuRankDist[i] = 0;
    }

    ofstream output;
    stringstream iss;
    time_t t = time(0);
    struct tm * now = localtime(&t);

    iss << "simK" << packetNum << "M" << batchSize << "m" << gf_order << "(" << now->tm_mon + 1 << now->tm_mday << now->tm_hour << now->tm_min << now->tm_sec << ").txt";

    output.open(iss.str().c_str());

    while (iter < iterationNum) {
        cout << "===== " << iter << " =====" << endl;
        sim.runOnce(timeUsed, ds);
        iter++;

        output << ds.nTrans << " " << ds.nReceive << " " << ds.nSave << " " << ds.nInact << " " << ds.nError << " ";

        cout << "Decoded/Total: " << packetNum - ds.nError << "/" << packetNum << " With " << ds.nTrans << " packets transmitted" << endl;

        double rate = (packetNum - ds.nError) / (float)ds.nReceive;

        cout << "Rate = " << rate << endl;

        output << rate << " ";

        cout << "Decoding time = " << timeUsed.decoding.time() << " Encoding time = " << timeUsed.encoding.time()  << endl;

        cout << "Rank Distribution: ";

        double Erk = 0.0;

        if (ds.nError == 0) {
            nSucc++;
            totalTrans += ds.nTrans;
        } else {
            errIdx[nErrIdx++] = iter-1;
        }


        for (int i = 0; i <= batchSize; i++) {
            output << ds.rankdist[i] << " ";
            cout << ds.rankdist[i] << " ";
            Erk += i * ds.rankdist[i];
            accuRankDist[i] += ds.rankdist[i];
        }

        output << "\n" << endl;

        cout << endl;

        cout << "E[rank(H)] = " << Erk << "(" << Erk / (float) batchSize << ")" << endl;
    }

    output.close();

    cout << "======== END =======" << endl;
    cout << "Simulation ends with " << nSucc << " succeeds out of " << iter << " runs" << endl;

    cout << "Average Rank Distribution: ";

    double totalRank = 0;
    double sumRank = 0;
    for(int i = 0; i <= batchSize; i++) {
        sumRank += accuRankDist[i];
    }
    for(int i = 0; i <= batchSize; i++) {
        accuRankDist[i] /= sumRank;
        totalRank += i * accuRankDist[i];
        cout << accuRankDist[i] << " ";
    }
    cout << endl;

    cout << "Average rate of succeeded runs = " <<  nSucc * packetNum / (float)totalTrans << " vs average rank = " << totalRank / (float) batchSize << endl;

    cout << "Error iterations: ";
    for (int i = 0; i < nErrIdx; i++) {
        cout << errIdx[i] << " ";
    }
    cout << endl;
}
