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

#ifndef BATCHENC_H_
#define BATCHENC_H_

#include <cstring>
#include <fstream>
#include <iostream>
#include "Utilities.h"
//#include "FiniteField.h"
#include "BatsBasic.h"

using namespace std;


class BatsEncoder : public BatsBasic{
private:

    // current batch ID
    KeyType batchID;

public:
  
    BatsEncoder(int M, SymbolType *input, int K, int T) : BatsBasic(M,K,T) {
        // init batch ID
        batchID = getSmallestBid();

        setInputPackets(input);
    }

    void setInputPackets(SymbolType *input);

    ~BatsEncoder() {
        if (packets != NULL)
            free(packets);
        if (checkPackets != NULL)
            freeMat(checkPackets, checkNum);
    }

    // generate a batch, return the degree
    void genBatchWithKey(SymbolType **batch, KeyType key);

    KeyType genBatch(SymbolType **batch){
        KeyType key = batchID;
        batchID++;
        genBatchWithKey(batch, key);
        return key;
    }
private:
    
    void genCheckPkg();
    
public:// debug
    // verify parity check
    bool verifyCheckPkg(){
        //
        SymbolType ** check = mallocMat<SymbolType>(checkNum,packetSize);
        
        // compute B_A * C_A
        for (int k = 0; k < smMinLd; k++) {
            for (int d = 0; d < ldpcVarDegree; d++) {
                FF.addvv(check[(k % ldpcNum + d * (int)(k / ldpcNum) + d) % ldpcNum], packets[k], packetSize);
            }
        }
        // compute [B_I H] * C_P
        for (int j = 0; j < ldpcNum; j++){
            for (int i = 0; i < 2; i++){
                int k = (j+i) % piNum;
                FF.addvv(check[j],getPkgHead(piToExt(k)),packetSize);
            }
        }
        // 
        
        // compute HDPC
        SymbolType** Bext = (SymbolType**)malloc(packetAndLDNum*sizeof(SymbolType*));    
        for (int i = 0; i < packetAndLDNum; i++){
            Bext[extToSM(i)] = getPkgHead(i);
        }
        
        matMulQ(& check[ldpcNum], Bext, packetSize, true);
        
        //
        int n = 0;
        
        for (int i = 0; i < checkNum; i++){
            for (int j = 0; j < packetSize; j++){
                if (check[i][j] != checkPackets[i][j]){
                    n++;
                    cout << "Wrong check pkg: " << i << endl;
                    continue;
                }
            }
        }
        
        return (n==0);
    }
};


#endif /* BATCHENC_H_ */
