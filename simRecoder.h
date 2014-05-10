/*
 Copyright (C) 2014 Tom M. Y. Lam, Shenghao Yang
 
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


#ifndef SIMRECODER_H
#define SIMRECODER_H

#include "Utilities.h"
#include "FiniteField.h"

#include <string.h>

const int shuffle_cnt = 10;

class SimRecoder {
public:
	SimRecoder(int M, int L, int Q, double* rank, MTRand* psrand) : M(M), L(L), Q(Q), rankSampler(rank, M+1), psrand(psrand) {
		H = mallocMat<SymbolType>(M, M);
        tmp = mallocMat<SymbolType>(M, M);
        psrand->seed();
	}
	~SimRecoder() {
		freeMat<SymbolType>(H, M);
        freeMat<SymbolType>(tmp, M);
	}
	int genBatch(SymbolType** dstBatchWithoutID, SymbolType** srcBatchWithoutID) {
		int curRank;

		curRank = rankSampler.sample(psrand->rand());
        
        for (int i = 0; i < curRank; i++) {
            memcpy(dstBatchWithoutID[i], srcBatchWithoutID[i], L*sizeof(SymbolType));
        }
        for (int i = curRank; i < M; i++) {
            memset(dstBatchWithoutID[i], 0, L*sizeof(SymbolType));
        }
        
        return curRank;
        
    }
private:
	//Parameters
	int M; //BatchSize
	int L; //Total packet length in bytes
	int Q; //No. of element in finite field minus one. For use in PRNG

	//Rank distribution
	DistSampler rankSampler;

	//Internal buffer for transfer matrix
	SymbolType** H;
    SymbolType** tmp;
    
	//Access to PRNG
	MTRand* psrand;
};


#endif
/* SIMRECODER_H */
