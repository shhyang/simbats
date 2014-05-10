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

#include "BatchDec.h"

//
// Save the recieved packet to the batch. The coefficients must be sent in. 
// The coefficnets saved in codingVec of the checknode has be modified by Gaussian elimination.
//


int CheckNode::processReceivedPacket(SymbolType *p, SymbolType* coef, int packetSize, int maxInact){

    if (numRec >= batchSize) {
        return 1;
    }
    
    if (!codingVecIndepend(coef)) {
        return 2;
    }
    
    memcpy(packet[numRec],p, packetSize);
    
    // update active and inactive coefficients of the packet
    VariableNode* pVar;
    SymbolType c;
    int i;
    
    for (i = 0; i < activeUndecodeDeg; i++){ // undecoded active variables
        c = FF.innerprod(edge[i]->g, coef, batchSize);
        edge[i]->gh[numRec] = c;
        //        edge[i]->addCoef(numRec,batchSize,codingVec[numRec]);
    }
    
    for (i = activeUndecodeDeg; i < activeDeg; i++){ // decoded active variables
        c = FF.innerprod(edge[i]->g, coef, batchSize);
        //        // substitute packet
        //        FF.addvvcCMP(packet[numRec], pVar->packet, c, packetSize);
        //        // substitute inactive
        //        FF.addvvc(inactCoef[numRec], pVar->inactCoef, c, maxInact);
        subsInPacket(numRec,edge[i]->vnode,c,packetSize,maxInact);
    }
    
    for (i = activeDeg; i < allDeg; i++){ // inactive variables
        pVar = edge[i]->vnode;
        c = FF.innerprod(edge[i]->g, coef, batchSize);
        // this may not be used
        // edge[i]->gh[numRec] = c;
        FF.incr(inactCoef[numRec][pVar->inactSeq], c);
    }
    
    numRec++;
    
    return 3;
}

// definitions of BatsDecoder


// name of FUNCTION: receivePacket
// 

void BatsDecoder::processReceivedPacket(SymbolType *packet, KeyType batchID){
    
    CheckNode* it = batchSet[batchID];

    // if a new batch is received, initialize the batch
    if (it == NULL){
        it = initNewBatch(batchID);
    }
    
    // extract the packet coding vector
    SymbolType codingVec[batchSize];
    int k = 0;
    SymbolType aMask;
    for (int i = 0; i < nSymInHead; i++){
        aMask = maskDec; // maskDec = 0xFF << (8-fieldOrder);
        for (int j = 0; j < nFFInSym; j++){
            codingVec[k] = (packet[i] & aMask) >> (nFFInSym - 1 - j)*fieldOrder;
            aMask >>= fieldOrder;
            k++;
        }
    }
    

    // save the packet in the batch
    int fl = it->processReceivedPacket(packet+nSymInHead, codingVec, packetSize, maxInact);
    
    if (fl == 1){
        nloss1++;
    } else if (fl == 2){
        nloss2++;
    } else {
        // receive the new packet
        nSavedPkg++;
    }

    // try to push decoding queue
    //decQueue->empty(); // why empty the queue?
    
    tryPushDecQueue(it);
}

CheckNode* BatsDecoder::initNewBatch(KeyType batchID) {
    int i,j;
    
    nRecBatch = (batchID > nRecBatch)? batchID : nRecBatch;

    psrand->seed(batchID);
    // generate the degree
    int degree = (int) getBatchDegree(batchID);

    CheckNode *it = new CheckNode(degree+piDegree, batchSize, packetSize, maxInact);

    it->id = batchID;
    
    // encoding: BATS parts
    
    int *idx = new int[degree];
    int *idxI = new int[piDegree];
    SymbolType **G = mallocMat<SymbolType>(degree,batchSize);
    SymbolType **GI = mallocMat<SymbolType>(piDegree,batchSize);
    genBatchParam(degree, idx, G, idxI, GI);
    
    for(i=0;i<degree;i++){
        
        BEdge* newEdge = it->addEdge(&(var[idx[i]]), batchSize);
        
        for (j = 0; j < batchSize; j++) {
            newEdge->g[j] = G[i][j];
        }
    }
    
    if (piNum > 0) {
        for (i = 0; i < piDegree; i++) {
            
            BEdge* newEdge = it->addEdge(&(var[idxI[i]]), batchSize);
            
            for (j = 0; j < batchSize; j++) {
                newEdge->g[j] = GI[i][j];
            }
        }
    }
    
    delete [] idx;
    delete [] idxI;
    freeMat(GI, piDegree);
    freeMat(G, degree);   

    batchSet[batchID] = it;

    return it;
}


bool CheckNode::decode(int packetSize, int maxInact){
    int rank = 0;
    SymbolType** GH = mallocMat<SymbolType>(batchSize,batchSize);
    SymbolType** invMat = mallocMat<SymbolType>(batchSize,batchSize);
    
    for (int i = 0; i < activeUndecodeDeg; i++) { // undecoded active variables
        for (int j = 0; j < numRec; j++){
            GH[j][i] = edge[i]->gh[j];
        }
    }
    
    // check the rank of GH
    for (int i = 0 ; i < numRec; i++){
        memset(invMat[i],0,numRec);
        invMat[i][i] = 1;
    }
    
    rank = FF.GaussianElimination(GH,invMat,numRec,activeUndecodeDeg,numRec);
    
    if (rank< activeUndecodeDeg){
        freeMat(GH,batchSize);
        freeMat(invMat,batchSize);
        return false;
    }
    
    // rank == degree, solve the batch!
    
    SymbolType** tPkg = mallocMat<SymbolType>(batchSize, packetSize);
    SymbolType** tICoef = mallocMat<SymbolType>(batchSize, maxInact);
    
    for (int i = 0; i < numRec; i++) {
        // process decoded packets
        FF.mulmcvCMP(tPkg[i], packet, invMat[i], numRec, packetSize);
        // process inactive coefficients
        FF.mulmcv(tICoef[i], inactCoef, invMat[i], numRec, maxInact);
    }
    
    for (int i = 0 ; i < numRec; i++) {
        memcpy(packet[i], tPkg[i], packetSize);
        memcpy(inactCoef[i], tICoef[i], maxInact);
    }
    
    freeMat(tICoef, batchSize);
    freeMat(tPkg, batchSize);
    freeMat(GH,batchSize);
    freeMat(invMat,batchSize);
    return true;
}


void BatsDecoder::decodeBatch(){
    CheckNode* it = decQueue->pop();
    it->inQueue = false;
    
    if (it->activeUndecodeDeg == 0){ // start other procedure
        copyForInacDec(it);
        return;
    }
    
    if (!(it->decode(packetSize, maxInact))) {
        return;
    }
    
    // copy decoded packets and inact coefficients
    VariableNode* pVar;
    // process the decoded variables of the batch
    for(int i = 0; i < it->activeUndecodeDeg; i++){
        pVar = it->edge[i]->vnode;
      
        pVar->decoded = true;
        
        // copy decoded packets
        memcpy(pVar->packet, it->packet[i], packetSize);
        // copy inactive coefficients
        memcpy(pVar->inactCoef, it->inactCoef[i], maxInact);
        
        nDecoded++;
        if (pVar->id<packetNum)
            nDecodedPkg++;
//        else
//            cout << "decoded an LDPC!" << endl;

        // process the related batches
        for (BEdge* d = pVar->edgeHead; d != NULL; d = d->nextInVar) {
            if (d->cnode == it)
                continue;

            d->cnode->subsDecodedVar(pVar,d,packetSize,maxInact);
            
            tryPushDecQueue(d->cnode);
        }
    }
    
    it->nUsedRec = it->activeUndecodeDeg;
    it->activeUndecodeDeg = 0; // must be done in the last step
    
    // process remaining part
    
    copyForInacDec(it);
}

bool BatsDecoder::solveInactVar(){
    
    // try to solve the inactive variables
    int packetAndLDNum = packetNum + ldpcNum;
    
    if (nInactDec==0) {
        
        SymbolType** Bext = (SymbolType**)malloc((packetNum + ldpcNum)*sizeof(SymbolType*));
        
        // compute Y1 * Q_active
        
        for (int i = 0; i < packetAndLDNum; i++){
            if (var[i].decoded)
            Bext[extToSM(i)] = getPkgHead(i);
            else
            Bext[extToSM(i)] = NULL;
        }
        
        matMulQ(& Y[maxC2], Bext, packetSize, true);
        
        // compute [C1 I] * [Q \\ I]
        SymbolType ** IM = mallocMat<SymbolType>(nInactVar, nInactVar);
        for (int i = 0; i < nInactVar; i++)
        IM[i][i] = 1;
        
        for (int i = 0; i < packetAndLDNum; i++){
            if (var[i].decoded)
            Bext[extToSM(i)] = inactCoefs[i];
            else
            Bext[extToSM(i)] = IM[var[i].inactSeq];
        }
        
        matMulQ(& C[maxC2], Bext, nInactVar, false);
        
        for (int i = 0; i < hdpcNum; i++){
            FF.incr(C[maxC2+i][var[packetAndLDNum+i].inactSeq]);
        }
        
        freeMat(IM, nInactVar);
        
        free(Bext);
    }
    //
    // decode inactive
    //
    int nCol = nC2+hdpcNum;
    SymbolType** inactPkg = mallocMat<SymbolType>(nInactVar, packetSize);
    
    int rank = FF.GaussianSolve(inactPkg, C+(maxC2-nC2), nInactVar, nCol, Y+(maxC2-nC2), packetSize);
    
    bool res = true;
    
    if (rank < nInactVar){
        cout << "InactDec: Inactive variables cannot be decoded!" << endl;
        res = false;
    } else{
        // substitute inactive
        for (int j = 0; j < totalNum; j++) {
            if (var[j].decoded){
                FF.addvmcvCMP(var[j].packet,inactPkg,var[j].inactCoef,nInactVar,packetSize);
            } else if (!var[j].active()){
                var[j].decoded = true;
                memcpy(var[j].packet, inactPkg[var[j].inactSeq], packetSize);
            }
        }
    }
    
    // free matrices
    freeMat(inactPkg, nInactVar);
    return res;
}


bool BatsDecoder::solveCY(int nCol){
    SymbolType** invC = mallocMat<SymbolType>(nCol, nCol);
    for(int i = 0; i < nCol; i++){
        invC[i][i] = 1;
    }
    int rank = FF.GaussianElimination(C+(maxC2-nC2),invC,nCol,nInactVar,nCol);
    
    if (rank < nInactVar){
        cout << "Inactive variables cannot be decoded!";
        freeMat(invC, nCol);
        return false;
    } else{
        // substitute inactive
        SymbolType** inactPkg = mallocMat<SymbolType>(nInactVar, packetSize);
        // get the pid of the ith inactive variable
//        int inactVarID[nInactVar];
//
//        for (int i = 0; i < totalNum; i++) {
//            if (!var[i].active()) {
//                inactVarID[var[i].inactSeq] = i;
//            }
//        }
    
        for (int i = 0; i < nInactVar; i++){
            FF.mulmcvCMP(inactPkg[i],Y+(maxC2-nC2),invC[i],nCol,packetSize);
        }
        
        for (int j = 0; j < totalNum; j++) {
            if (var[j].decoded){
                FF.addvmcvCMP(var[j].packet,inactPkg,var[j].inactCoef,nInactVar,packetSize);
//                var[j].clearInact();
            } else if (!var[j].active()){
                var[j].decoded = true;
                memcpy(var[j].packet, inactPkg[var[j].inactSeq], packetSize);
                // set as active
//                var[j].inactSeq = -1;
            }
        }
        
        freeMat(inactPkg, nInactVar);
    }
    
    // free matrices
    freeMat(invC, nCol);
    return true;
}



bool BatsDecoder::addInact(){
    
    // find a variable for inactivation
    int cDeg=-1;
    int cGap=-1000;
    int cIdx=-1;
    int tmp;
    CheckNode *it;
    
    // find a batch with most likely decodable. 
    for (int i=0; i < nRecBatch; i++){
        it = batchSet[i];
        if (it!=NULL && !it->decoded()){
            tmp = it->numRec - it->activeUndecodeDeg;
            if (tmp > cGap) {
                cGap = tmp;
                cDeg = it->activeUndecodeDeg;
                cIdx = i;
            } else if (tmp == cGap && it->activeUndecodeDeg > cDeg) {
                cDeg = it->activeUndecodeDeg;
                cIdx = i;
            }   
        }
    }
    
    if (cIdx < 0 || cDeg <= 0){
        return false;
    }
    
    // find a variable node
    VariableNode* vn=NULL;
    cDeg = -1;
    
    it = batchSet[cIdx];
    for (int i=0;i<it->activeUndecodeDeg;i++){
        if(it->edge[i]->vnode->degree > cDeg){
            vn = it->edge[i]->vnode;
            cDeg = it->edge[i]->vnode->degree;
        }   
    }
    
    // set vn inactive
    vn->inactSeq = nInactVar; 
    nInactVar ++; 
    
    // substitute the inact variable
    //decQueue->empty();
    // process the related batches
    for (BEdge* d = vn->edgeHead; d != NULL; d = d->nextInVar) {
        d->cnode->addInact(d,vn->inactSeq);    
        tryPushDecQueue(d->cnode);
    }
    
    while(decQueue->isNonEmpty()){
        decodeBatch();
    }
    
    return true;
}
