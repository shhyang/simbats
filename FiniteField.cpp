/*
 Copyright (C) 2007 James S. Plank
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

#include "FiniteField.h"

const int FiniteField::prim_poly[9] = { 0,
/*  1 */1,
/*  2 */07,   //numbers in octal
/*  3 */013,
/*  4 */023,
/*  5 */045,
/*  6 */0103,
/*  7 */0211,
/*  8 */0435 };
const int FiniteField::nw[9] = { 0, (1 << 1), (1 << 2), (1 << 3), (1 << 4), (1 << 5), (1 << 6), (1 << 7), (1 << 8) };
const int FiniteField::nwm1[9] = { 0, (1 << 1) - 1, (1 << 2) - 1, (1 << 3) - 1, (1 << 4) - 1, (1 << 5) - 1, (1 << 6) - 1, (1 << 7) - 1, (1 << 8) - 1 };



void FiniteField::init() {
	divTable = NULL;
	mulTable = NULL;
	logTable = NULL;
	expTable = NULL;
        mulTable2 = NULL;
   if (createMulintables(order) < 0) {
		cout << "create_mul_table Error: order=%d\n" << order;
	}
}

void FiniteField::finalize() {
	if (logTable != NULL) {
                delete [] logTable;
		logTable = NULL;
	}
	if (expTable != NULL) {
		delete [] expTable;
		expTable = NULL;
	}
	if (mulTable != NULL) {
		delete [] mulTable;
		mulTable = NULL;
	}
	if (divTable != NULL) {
		delete [] divTable;
		divTable = NULL;
	}
        if (order != 8 && mulTable2 != NULL) {
		delete [] mulTable2;
		mulTable2 = NULL;
	}
}

int FiniteField::createLogTables(int w) {
	int j, b;

	if (w >= 9)
		return -1;

	logTable = new int[nw[w]]; //(int *) malloc(sizeof(int) * nw[w]);
	if (logTable == NULL)
		return -1;

	expTable = new int[nw[w] * 4]; //(int *) malloc(sizeof(int) * nw[w] * 3);
	if (expTable == NULL) {
		delete[] logTable;
		logTable = NULL;
		return -1;
	}

	for (j = 0; j < nw[w]; j++) {
		logTable[j] = nwm1[w];
		expTable[j] = 0;
	}

	b = 1;
	for (j = 0; j < nwm1[w]; j++) {
		if (logTable[b] != nwm1[w]) {
			fprintf(
					stderr,
					"create_log_table Error: j=%d, b=%d, B->J[b]=%d, J->B[j]=%d (0%o)\n",
					j, b, logTable[b], expTable[j], (b << 1) ^ prim_poly[w]);
			return -1;
		}
		logTable[b] = j;
		expTable[j] = b;
		b = b << 1;
		if (b & nw[w])
			b = (b ^ prim_poly[w]) & nwm1[w];
	}

	logTable[0] = 2*nwm1[w];

	for (j = 0; j < nwm1[w]; j++) {
		expTable[j + nwm1[w]] = expTable[j];
		expTable[j + nwm1[w] * 2] = 0;
		expTable[j + nwm1[w] * 3] = 0;
	}
	return 0;
}

int FiniteField::createMulintables(int w) {
    int j, x, y, logx;

    if (w >= 9)
        return -1;

    mulTable = new int[nw[w] * nw[w]]; //(int *) malloc(sizeof(int) * nw[w] * nw[w]);
    if (mulTable == NULL)
        return -1;

    divTable = new int[nw[w] * nw[w]]; //(int *) malloc(sizeof(int) * nw[w] * nw[w]);
    if (divTable == NULL) {
        delete [] mulTable;
        mulTable = NULL;
        return -1;
    }

    if (logTable == NULL) {
        if (createLogTables(w) < 0) {
            delete [] mulTable;
            delete [] divTable;
            mulTable = NULL;
            divTable = NULL;
            return -1;
        }
    }

    /* Set mult/div tables for x = 0 */
    j = 0;
    mulTable[j] = 0; /* y = 0 */
    divTable[j] = -1;
    j++;
    for (y = 1; y < nw[w]; y++) { /* y > 0 */
        mulTable[j] = 0;
        divTable[j] = 0;
        j++;
    }

    for (x = 1; x < nw[w]; x++) { /* x > 0 */
        mulTable[j] = 0; /* y = 0 */
        divTable[j] = -1;
        j++;
        logx = logTable[x];
        for (y = 1; y < nw[w]; y++) { /* y > 0 */
            mulTable[j] = expTable[logx + logTable[y]];
            divTable[j] = expTable[logx + nwm1[w] - logTable[y]];
            j++;
        }
    }

    // set mulTable2
    if (order == 8){
        mulTable2 = mulTable;
        return 0;
    }

    mulTable2 = new int[nw[8] * nw[w]]; //(int *) malloc(sizeof(int) * nw[w] * nw[w]);
    if (mulTable2 == NULL) {
        return -1;
    }

    unsigned char mask = 0xFF;
    mask <<= (8-order);
    unsigned char aMask, a;
    int nSym = 8/order;

    j = 0;
    for (x = 0; x < nw[8]; x++){
        for (y = 0; y< nw[w]; y++){
            aMask = mask;
            mulTable2[j] = 0;
            a = x;
            for (int k = 0; k< nSym; k++){
                a = (unsigned char)x & aMask;
                a >>= (nSym-1-k)*order;
                a = mulTable[(a << order) | y];
                a <<= (nSym-1-k)*order;
                mulTable2[j] |= (int)a;
                aMask >>= order;
            }
            j++;
        }
    }

    return 0;
}


// rank of a column matrix
int FiniteField::rankM(SymbolType** a, int N, int M) {
    int r = 0;
    int i,j;
    int MN = (M<N) ? M : N;

    for (i = 0; i < MN; i++) {
        int j0 = -1;
        while (1) {
            for (j = i; j < N; j++) {
                if (a[j][r] != 0) {
                    j0 = j;
                    break;
                }
            }
            if (j0 < 0) {
                r++;
                if (r >= M)
                    return i;
                continue;
            }
            swap(a[i],a[j0]);
            break;
        }

        {
            int c = div(1, a[i][r]);
            a[i][r] = 1;
            for (j = r + 1; j < M; j++)
                a[i][j] = mul(a[i][j], c);
        }

//        for (int i2 = 0; i2 < N; i2++)
//            if (i2 != i && a[i2][r] != 0) {
//                int c = a[i2][r];
//                a[i2][r] = 0;
//                for (j = r + 1; j < M; j++) {
//                    a[i2][j] = sub(a[i2][j], mul(a[i][j], c));
//                }
//            }
        for (int i2 = i+1; i2 < N; i2++)
            if (a[i2][r] != 0) {
                int c = a[i2][r];
                a[i2][r] = 0;
                for (j = r + 1; j < M; j++) {
                    a[i2][j] = sub(a[i2][j], mul(a[i][j], c));
                }
            }
    }
    return i;
}
//
// solve linear system : Y_{rowY x col} = X * A_{rowA x col}

int FiniteField::GaussianSolve(SymbolType** X,  SymbolType** A, int rowA, int col, SymbolType** Y, int rowY, bool fu){
    SymbolType** invMat = mallocMat<SymbolType>(col,col);
    SymbolType** A2 = mallocMat<SymbolType>(col, rowA);

    for (int i = 0; i < col; i++){
        memset(invMat[i],0,col);
        invMat[i][i] = 1;
        memcpy(A2[i], A[i], rowA);
    }
    
    int rank = FF.GaussianElimination(A2,invMat,col,rowA, col);
    
    if (rank == rowA) {// rank == row
        int up = rank;
        if (fu)
            up = col;

        for (int i = 0; i < up; i++) {
            mulmcvCMP(X[i], Y, invMat[i], col, rowY);
        }
    }
    
    freeMat(invMat,col);
    freeMat(A2, col);
    return rank;
}

int FiniteField::GaussianElimination(SymbolType** A, SymbolType** Y, int col, int rowA, int rowY) {
    
    int i,j;

    for(i=0;i<rowA;i++){
        int i0 = -1;
        for (j = i; j < col; j++) {
            if (A[j][i] != 0) {
                i0 = j;
                break;
            }
        }
        if (i0 < 0) {
            break;
        }
        swap(A[i],A[i0]);
        swap(Y[i],Y[i0]);

        {
            int c = div(1, A[i][i]);
            A[i][i] = 1;
            for (j = i + 1; j < rowA; j++)
                A[i][j] = mul(A[i][j], c);

            mulvc(Y[i],c, rowY);
        }

        for (int i2 = 0; i2 < col; i2++)
            if (i2 != i && A[i2][i] != 0) {
                int c = A[i2][i];
                A[i2][i] = 0;
                for (j = i + 1; j < rowA; j++) {
                    A[i2][j] = sub(A[i2][j], mul(A[i][j], c));
                }
                addvvc(Y[i2], Y[i], c, rowY);
            }

        //               printf("A(%d)=\n",i);
        //        for(int s=0; s<A[0].size(); s++){
        //            for(int t=0;t<A.size();t++){
        //                int a = A[t][s];
        //                printf("%d\t", a);
        //            }
        //            printf("\n");
        //        }
    }

    return i;
}


FiniteField FF(1);

