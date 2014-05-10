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

#ifndef BATCHDEC_H_
#define BATCHDEC_H_

#include <cstdlib>
#include <cstdio>
#include <cstring>
//#include "Utilities.h"
//#include "FiniteField.h"
#include "BatsBasic.h"

using namespace std;

struct CheckNode;
struct VariableNode;

class BEdge{
public:
    // variable node
    VariableNode* vnode;
    // check node
    CheckNode *cnode;
    int seqInCheck;
    // g is a vector of size batchSize
    // it is the row of the generator matrix of cnode corresponding vnode;
    SymbolType* g;
    // gh is a vector of the coefficients of the received packets
    // the size of the vector is equal to numRec in BatsDecoder
    // a component of gh is given by inner product of g and coding vector of a packet
    SymbolType* gh;

    BEdge *nextInVar;

    BEdge(CheckNode * c, VariableNode * v, int size): cnode(c), vnode(v){
        g = new SymbolType[size];
        gh = new SymbolType[size];
        nextInVar = NULL;
    }
    ~BEdge(){
        delete[] g;
        delete[] gh;
    }

    // set gh for the index idx
    inline void addCoef(int idx, int size, SymbolType* cv){
        gh[idx] = FF.innerprod(g,cv,size);
    }
};

class VariableNode{
private:

public:
    int id;
    int degree;
    int inactSeq;
    BEdge *edgeHead;

    bool decoded;
    SymbolType* packet;
    SymbolType* inactCoef;

    VariableNode(){
        degree = 0;
        edgeHead = NULL;
        inactSeq = -1;

        decoded = false;
    }

    ~VariableNode(){
    }

    inline void addEdge(BEdge* edge){
        edge->nextInVar = edgeHead;
        edgeHead = edge;
        degree++;
    }

    inline bool active(){
        return inactSeq < 0;
    }

};

class CheckNode {
private:

public:
    // received packets of the check node, not including the coding vector
    // the matrix Y_i
    SymbolType** packet;
    // coding vectors extracted from the packets
    // the matrix H_i
    SymbolType** codingVec;
    // rank of the current coding vectors
    int codingRank;
    // number of received packets
    int numRec;
    // number of received apckets that have been used for decoding
    int nUsedRec;

    // id of the check node
    int id;

    // the number of active and undecoded variables
    int activeUndecodeDeg;

    // the number of active variables, including decoded
    int activeDeg;

    // original degree, including inactivated
    int allDeg;

    // edge of check node
    // in edge, edges are ordered as : active undecoded, active decoded, inactivated
    BEdge **edge;

    // inactive coefficients
    SymbolType** inactCoef;

    int batchSize;

    int maxDeg;

    bool inQueue;
//    bool decoded; // replaced by decoded()

    CheckNode(int maxDeg, int batchSize, int packetSize, int maxInact) : maxDeg(maxDeg), batchSize(batchSize) {
        allDeg = 0;
        activeUndecodeDeg = 0;
        activeDeg = 0;

        numRec = 0;
        nUsedRec = 0;
        inQueue = false;

        edge = new BEdge*[maxDeg];
        for (int i = 0; i < maxDeg; i++){
            edge[i] = NULL;
        }
        inactCoef = mallocMat<SymbolType>(batchSize, maxInact);

        packet = mallocMat<SymbolType > (batchSize, packetSize);
        codingVec = mallocMat<SymbolType>(batchSize, batchSize);
        codingRank = 0;
    }

    ~CheckNode() {
        // free BiEdge
        for (int i = 0; i < maxDeg; i++) {
            if (edge[i]!=NULL)
                delete edge[i];
        }

        delete[] edge;

        freeMat(inactCoef, batchSize);
        // free packet
        freeMat(packet, batchSize);
        freeMat(codingVec, batchSize);
    }

    inline bool decoded(){
        return  (activeUndecodeDeg==0);
    }

    inline void exchangeEdge(int first, int second){
        BEdge* temp = edge[first];
        edge[first] = edge[second];
        edge[second] = temp;
        edge[first]->seqInCheck = first;
        edge[second]->seqInCheck = second;
    }

    inline BEdge* addEdge(VariableNode * v, int size){
        BEdge* e = new BEdge(this, v, size);
        v->addEdge(e);
        edge[allDeg] = e;
        e->seqInCheck = allDeg;
        allDeg++;
        if (v->active()){
            activeDeg++;
            if (activeDeg < allDeg) {
                exchangeEdge(activeDeg-1, allDeg-1);
            }
            if (!(v->decoded)){
                activeUndecodeDeg++;
                // change the order of edges   uuuaau
                if (activeUndecodeDeg < activeDeg) {
                    exchangeEdge(activeUndecodeDeg-1, activeDeg-1);
                }
            }
        }
        return e;
    }

    inline void subsInPacket(int idx, VariableNode* var, int coef, int packetSize, int maxInact){
        FF.addvvcCMP(packet[idx], var->packet, coef, packetSize);
        FF.addvvc(inactCoef[idx], var->inactCoef, coef, maxInact);
    }

    inline void subsDecodedVar(VariableNode* var, BEdge* e, int packetSize, int maxInact) {
        activeUndecodeDeg--;
        // substitute the value of pVar into the batch
        for (int i = 0; i < numRec; i++) {
            //        FF.addvvcCMP(packet[i], var->packet, e->gh[i], packetSize);
            //        FF.addvvc(inactCoef[i], var->inactCoef, e->gh[i], maxInact);
            subsInPacket(i, var, e->gh[i], packetSize, maxInact);
        }

        exchangeEdge(e->seqInCheck, activeUndecodeDeg);
    }

    inline bool codingVecIndepend(SymbolType * cv) {

        memcpy(codingVec[numRec], cv, batchSize);

        // check if the packet takes new information

        int hRank = FF.rankM(codingVec, numRec + 1, batchSize);

        if (hRank <= codingRank) {
            return false;
        } else {
            codingRank = hRank;
            return true;
        }
    }

    inline void addInact(BEdge *e, int seq){
        activeUndecodeDeg--; // since pid should be active
        activeDeg--;

        // move d to inactive part in edge in cb
        exchangeEdge(e->seqInCheck, activeUndecodeDeg);
        exchangeEdge(activeUndecodeDeg, activeDeg);

        for (int i = 0; i < numRec; i++)
            inactCoef[i][seq] = e->gh[i];
    }

    int processReceivedPacket(SymbolType *,SymbolType *, int, int);
    
    bool decode(int packetSize, int maxInact);
};


class BatsDecoder : public BatsBasic{
public:

    // current batch ID
    KeyType batchID;

    // number of received packets
    int nloss1,nloss2,nloss3;
    int nRecPkg;

    int nSavedPkg;

    int nInactVar;

    // number of total decoded packets, including check packets
    int nDecoded;

    // number of decoded original packets
    int nDecodedPkg;

private:

    // Sparse matrix decoding
    int nRecBatch;
    CheckNode* batchSet[BATSDECODER_MAXBATCH]; // save received batch

    VariableNode* var;

    ArrayQueue<CheckNode*> *decQueue;

//    int inactRedun; // the redundancy in decoding inactive
//    int receiRedun; // the redundancy in received packets
    int maxInact; // the maximum allowed number of inactivation
    int inactDecMark; // mark the position of the last run of inactivation decoding
    int nInactDec; // count the number of inact decoding has been run

    SymbolType** inactCoefs;

    int nC2; // number of C2 packets
    int maxC2; // the maximum number of C2 packets
    // temple variables
    SymbolType** tC2;
    SymbolType** Y; // inactivation decoding Y = [Y2 Y1QA]
    SymbolType** C; // inactivation decoding C = [C2 C1QA-QI]

public:
    BatsDecoder(int M, SymbolType* output, int K, int T): BatsBasic(M,K,T){

        setOutputPacket(output);
        // Sparse Matrix Decoding
        nloss1 = 0;
        nloss2 = 0;
        nloss3 = 0;

        nRecPkg = 0;
        nSavedPkg = 0;
        nInactVar = piNum;
        nRecBatch = ldpcNum;
        nDecoded = 0;
        nDecodedPkg = 0;

        decQueue = new ArrayQueue<CheckNode*>(BATSDECODER_MAXBATCH);
        decQueue->empty();

        // maximum allowed number of inactive variables
        maxInact = 3*sqrt(packetNum);
        maxInact = (piNum > maxInact)? piNum:maxInact;

//        receiRedun = 16;
//        inactRedun = maxInact*2;
//        inactRedun = (inactRedun > 8)? inactRedun:8;
        inactDecMark = 0;
        nInactDec = 0;

        inactCoefs = mallocMat<SymbolType>(totalNum, maxInact);

        // initial variable nodes
        var = new VariableNode[totalNum];

        int idx = 0;
        for(int i = 0; i<totalNum; i++){
            if (isPreInact(i)){
                var[i].inactSeq = idx;
                idx++;
            }
            var[i].id = i;
            var[i].packet = getPkgHead(i);
            var[i].inactCoef = inactCoefs[i];
        }
        // init batchSet
        for(int i=ldpcNum;i<BATSDECODER_MAXBATCH;i++){
            batchSet[i] = NULL;
        }

        // init temp variables
        nC2 = 0;
        maxC2 = maxInact;
        Y = mallocMat<SymbolType>(maxC2+hdpcNum,packetSize);
        C = mallocMat<SymbolType>(maxC2+hdpcNum,maxInact);
        tC2 = mallocMat<SymbolType>(maxC2,maxInact);
        
        // LDPC parameters init
        if (ldpcNum > 0) {
            int ldpcBlockNum = smMinLd / ldpcNum;
            int ldpcMaxCheckDeg = ldpcVarDegree * (ldpcBlockNum + 1) + 1 + 2; // including pi part

            for(int i = 0; i<ldpcNum; i++){
                batchSet[i] = new CheckNode(ldpcMaxCheckDeg,1,packetSize, maxInact);
                CheckNode* it = batchSet[i];
                it->id = i;// used for debug
                it->codingVec[0][0] = 1;
                it->codingRank = 1;
                it->numRec++;
                memset(it->packet[0], 0, packetSize);
            }

            int nb, nc, nd;
            for (int i = 0; i < smMinLd; i++) {
                nb = i / ldpcNum;
                nc = i % ldpcNum;
                for (int j = 0; j < ldpcVarDegree; j++) {
                    nd = (nc + j * nb + j) % ldpcNum;
                    BEdge* newEdge = batchSet[nd]->addEdge(&(var[i]), 1);
                    newEdge->g[0] = 1;
                    newEdge->gh[0] = 1;
                }
            }
            for (int i = 0; i < ldpcNum; i++) {
                BEdge* newEdge = batchSet[i]->addEdge(&(var[i+packetNum]), 1);
                newEdge->g[0] = 1;
                newEdge->gh[0] = 1;
            }
            // add inactive edges in the end
            for (int j = 0; j < ldpcNum; j++) {
                for (int i = 0; i < 2; i++) {
                    int k = piToExt((j + i) % piNum);
                    BEdge* newEdge = batchSet[j]->addEdge(&(var[k]), 1);
                    newEdge->g[0] = 1;
                    newEdge->gh[0] = 1;
                    batchSet[j]->inactCoef[0][var[k].inactSeq] = 1;
                }
            }
        }
    }

    ~BatsDecoder(){
        delete [] var;

        for(int i=0;i<BATSDECODER_MAXBATCH;i++){
            if(batchSet[i] != NULL)
                delete batchSet[i];
        }

        freeMat(Y,maxC2+hdpcNum);
        freeMat(C,maxC2+hdpcNum);
        freeMat(tC2, maxC2);
        freeMat(inactCoefs, totalNum);

        if (packets != NULL)
            free(packets);
        if (checkPackets != NULL)
            freeMat(checkPackets, checkNum);

        delete decQueue;
    }

    void setOutputPacket(SymbolType *output){
        packets = (SymbolType**)malloc(packetNum*sizeof(SymbolType*));

        for (int i = 0; i < packetNum; i++){
            packets[i] = output + i*packetSize;
        }

        if (checkNum <= 0)
            return;

        // precode

        if (checkPackets != NULL){
            freeMat(checkPackets,checkNum);
            checkPackets = NULL;
        }

        checkPackets = mallocMat<SymbolType>(checkNum, packetSize);
    }

    inline bool complete(double decRatio){
        return (nDecodedPkg>=packetNum * decRatio);
    }

    void rankDist(double* rd){
        int dd[batchSize + 1];
        double numRec = 0;
        for (int j = 0; j <= batchSize; j++) {
            dd[j] = 0;
        }
        CheckNode * it;
        for (int i = getSmallestBid(); i < nRecBatch; i++) {
            it = batchSet[i];
            if (it != NULL) {
                numRec += 1;
                dd[it->codingRank]++;
            }
        }

        for (int j = 0; j <= batchSize; j++) {
            rd[j] = dd[j] / (double) numRec;
        }
    }

    void receivePacket(SymbolType *packet, KeyType batchID){
        nRecPkg++;

        processReceivedPacket(packet,batchID);

        while(decQueue->isNonEmpty()){
            decodeBatch();
        }

        inactDecoding();
    }

    void inactDecoding(){

        if (nSavedPkg >= packetNum){
            bool res = true;
            int c = 0;

            while ((c<=16) && res && nInactVar < maxInact && (nDecoded + nInactVar < totalNum))
            {
                // inactivate a variable
                res = addInact();
                c++;
            }
        }


        if (nSavedPkg - inactDecMark < 16) {
//            cout << "inactDecoding " << nInactDec << " : nRecPkg = " << nRecPkg << " loss(decoded) = " << nloss1 << " loss(dependant) = " << nloss2 << endl;
            return;
        }

        if (nDecoded + nInactVar < totalNum || nC2 + hdpcNum < nInactVar) {
 //           cout << "InactDec: not ready to run inactDec(). Receive more packets!" << endl;
            return;
        }

        if(solveInactVar()){
            nDecodedPkg = packetNum;
            cout << "InactDec run "<< nInactDec+1 << " succeeds!" << endl;
        } else {
            cout << "InactDec run "<< nInactDec+1 << " fails!" << endl;
            inactDecMark = nSavedPkg;
        }
        cout << "InactDec: nRecPkg = " << nRecPkg << " nSavedPkg = " << nSavedPkg << " loss(dependant) = " << nloss2 << endl;
        cout << "InactDec: nDecoded = " << nDecoded <<" nInactVar = " << nInactVar << " nC2 = " << nC2 << " nHDPC = " << hdpcNum << endl;
        nInactDec ++;
    }

    int getDecodedPkg(){
        int n = 0;
        for(int j = 0; j < packetNum; j++){
            if(var[j].decoded){
                n++;
// same location      memcpy(getPkgHead(j),var[j].packet,packetSize);
            }
        }
        return n;
    }
private:
    // sparse matrix
    CheckNode* initNewBatch(KeyType);
    inline void tryPushDecQueue(CheckNode* it){

        if(!(it->inQueue) && (it->numRec) >= (it->activeUndecodeDeg)){
            it->inQueue = true;
            decQueue->push(it);
        }
    }
    inline void copyForInacDec(CheckNode* it){
        if (nC2 == maxC2) {
            return;
        }
        for (int i = it->nUsedRec; i < it->numRec; i++) {
            
            memcpy(tC2[nC2], it->inactCoef[i], maxInact);
            int hRank = FF.rankM(tC2, nC2+1, nInactVar);
            if (hRank <= nC2) {
                continue;
            }
            
            // Save to C
            memcpy(C[maxC2 - 1 - nC2], it->inactCoef[i], maxInact);
            // Save to Y
            memcpy(Y[maxC2 - 1 - nC2], it->packet[i], packetSize);
            nC2++;
            it->nUsedRec++;
            if (nC2==maxC2) {
                break;
            }
        }

    }
    void decodeBatch();
    bool addInact();
    bool solveCY(int);
    bool solveInactVar();
    void processReceivedPacket(SymbolType *packet, KeyType batchID);
};

#endif /* BATCHDEC_H_ */
