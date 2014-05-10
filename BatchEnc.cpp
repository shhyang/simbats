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

#include "BatchEnc.h"

void BatsEncoder::genCheckPkg(){
//    int packetAndLDNum = packetNum + ldpcNum;
    
    SymbolType** ldpcHeader = & checkPackets[0];
    
    SymbolType** hdpcHeader = & checkPackets[ldpcNum];
    
    SymbolType** Bext = (SymbolType**)malloc(packetAndLDNum*sizeof(SymbolType*));
    
    // matrix for C_H
    SymbolType** CH = mallocMat<SymbolType>(ldpcNum, hdpcNum);
    // matrix for [I_H C_H * Q_L]
    SymbolType** AH = mallocMat<SymbolType>(hdpcNum, hdpcNum);
    // matrix for [B_A B* BI] * Q
    SymbolType** BQ = mallocMat<SymbolType>(hdpcNum, packetSize);
    
    // Compute B* = B_A C_A + B_I C_I and save in LDPC
    // compute B_A * C_A
    for (int k = 0; k < smMinLd; k++) {
        for (int d = 0; d < ldpcVarDegree; d++) {
            FF.addvv(ldpcHeader[(k % ldpcNum + d * (int)(k / ldpcNum) + d) % ldpcNum], packets[k], packetSize);
        }
    }
    // compute B_I * C_I and C_H
    for (int j = 0; j < ldpcNum; j++){
        for (int i = 0; i < 2; i++){
            int k = (j+i) % piNum;
            if (k < piMinHd){
                FF.addvv(ldpcHeader[j],packets[smMinLd+k],packetSize);
            } else { // computer CH
                CH[j][k-piMinHd] = 1;
            }
        }
    }
    
    // Compute BQ =  [B_A B* BI] * Q
        
    for (int i = 0; i < packetAndLDNum; i++){
        Bext[extToSM(i)] = getPkgHead(i);
    }
    
    matMulQ(BQ, Bext, packetSize, true);
     
    // Compute AH = I_H + C_H Q_L    

    for (int i = 0; i < packetAndLDNum; i++){
        if (i < packetNum)
            Bext[extToSM(i)] = NULL;
        else
            Bext[extToSM(i)] = CH[i-packetNum];
    }
    
    matMulQ(AH, Bext, hdpcNum, false);
    
    for(int i = 0; i < hdpcNum; i++){
        FF.incr(AH[i][i]);
    }
    
    // solve H * AH = BQ
    
    int rank_AH = FF.GaussianSolve(hdpcHeader, AH, hdpcNum, hdpcNum, BQ, packetSize);
    
    if (rank_AH < hdpcNum){
        cout << "Fatal encoding error: Cannot generate HDPC packets!" << endl;
    }
    
    // solve L = H * CH + B*
    
    for (int i = 0; i < ldpcNum; i++){
        FF.addvmcvCMP(ldpcHeader[i], hdpcHeader, CH[i], hdpcNum, packetSize);
    }
    
    // free allocated memory
    free(Bext);
    freeMat(CH, ldpcNum);
    freeMat(AH, hdpcNum);
    freeMat(BQ, hdpcNum);
}

void BatsEncoder::setInputPackets(SymbolType *input){
    packets = (SymbolType**)malloc(packetNum*sizeof(SymbolType*));
    
    for (int i = 0; i < packetNum; i++){
        packets[i] = input + i*packetSize;
    }
    
    if (checkNum <= 0)
        return;

    // precode

    if (checkPackets != NULL){
        freeMat(checkPackets,checkNum);
        checkPackets = NULL;
    }

    checkPackets = mallocMat<SymbolType>(checkNum, packetSize);
    
    genCheckPkg(); 
}


void BatsEncoder::genBatchWithKey(SymbolType **batch, KeyType key){
    int i,j;
    SymbolType* t0;

    uint16 degree = (uint16)getBatchDegree(key);

    int Tout = packetSize + nSymInHead;

    // reset
    for (i = 0 ; i< batchSize ; i ++ )
        memset(batch[i], 0 , Tout);

    // set header
    SymbolType aMask;
    for (i = 0 ; i < nSymInHead; i++){
        aMask = maskEnc;
        for (j = 0 ; j < nFFInSym; j++){
            batch[i*nFFInSym+j][i] |= aMask;
            aMask >>= fieldOrder;
        }
    }
    
    // encoding: BATS parts
    
    int *idx = new int[degree];
    int *idxI = new int[piDegree];
    SymbolType **G = mallocMat<SymbolType>(degree,batchSize);
    SymbolType **GI = mallocMat<SymbolType>(piDegree,batchSize);
    genBatchParam(degree, idx, G, idxI, GI);
    
    for(i=0;i<degree;i++){

        t0 = getPkgHead(idx[i]);
        
        for(j=0;j<batchSize;j++){
            FF.addvvcCMP(batch[j]+nSymInHead, t0, G[i][j], packetSize);
        }
    }

    // encoding: PI parts
    //degree = piDegree;

    if (piNum > 0) {
        for (i = 0; i < piDegree; i++) {
            t0 = getPkgHead(idxI[i]);

            for (j = 0; j < batchSize; j++) {
                FF.addvvcCMP(batch[j] + nSymInHead, t0, GI[i][j], packetSize);
            }
        }
    }
    
    delete [] idx;
    delete [] idxI;
    freeMat(GI, piDegree);
    freeMat(G, degree);
}
