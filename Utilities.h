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

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
//#include <vector>
#include "MersenneTwister.h"

//#ifndef DEBUG_BATS
//#define DEBUG_BATS
//#endif

using namespace std;


#define SymbolSize 8
#define MAX_DEGREE 2000

typedef unsigned char SymbolType;
typedef unsigned short KeyType;
typedef unsigned long uint32;
typedef unsigned short uint16;
//typedef vector<SymbolType> SymVector;

inline void saveIDInPacket(SymbolType* packet, KeyType* id){
//    SymbolType* a = (SymbolType*)id;
//    for(int i = 0; i< sizeof(KeyType); i++){
//        packet[i] = a[i];
//    }
    memcpy(packet, (SymbolType*) id, sizeof(KeyType));
}

inline KeyType getIDFromPacket(SymbolType* packet){
    KeyType id;
//    SymbolType* a = (SymbolType*)&id;
//    for(int i = 0; i< sizeof(KeyType); i++){
//        a[i] = packet[i];
//    }
    memcpy((SymbolType*)&id, packet, sizeof(KeyType));
    return id;
}

class DistSampler {
    double accu[MAX_DEGREE+1];
    int n;
public:
    DistSampler(double* dist, int D):n(D){
        double sum = 0.0;
        
        accu[0] = 0.0;

        for(int i = 0; i<D; i++){
            sum += dist[i];
            accu[i+1] = accu[i] + dist[i];
        }
    }

    inline int sample(double pos){ // pos in [0 1]
        int lb = 0, ub = n, mi;
		while(ub-lb > 1){
			mi = lb+(ub-lb)/2;
			if(pos < accu[mi])
				ub = mi;
			else
				lb = mi;
		}
		return lb;
    }
};


template<class T>
inline T** mallocMat(int row, int col){
    T** p = (T**)malloc(row*sizeof(T*));
    for(int i = 0 ; i<row; i++){
        p[i] = (T*)malloc(col*sizeof(T));
        memset(p[i],0,col*sizeof(T));
    }
    return p;
}

template<class T>
inline void freeMat(T** p, int row){
    for(int i = 0 ; i<row; i++){
        free(p[i]);
    }
    free(p);
}

template<class T>
class ArrayQueue {
private:
    T* data;
    int head;
    int tail;
    int bufsize;

public:
    ArrayQueue(int s):bufsize(s){
        data = new T[bufsize];
        tail = 0;
        head = 0;
    }
    ~ArrayQueue(){
        delete [] data;
    }

    inline int size() {
        return (tail-head);
    }
    inline void empty(){
        head = 0;
        tail = 0;
    }

    inline bool isNonEmpty(){
        return head < tail;
    }

    inline void push(T x){
        data[(tail++)%bufsize] = x;
    }

    inline T pop(){
        return data[(head++)%bufsize];
    }
};


class WrongDegree{
public:
    WrongDegree(){}
    WrongDegree(int theDegree): degree(theDegree){}
    int getDegree() const {return degree;}
private:
    int degree;
};


template <class T>
void printMat(T* a, int row, int col){
    for (int i = 0 ; i< row; i++) {
        for (int j = 0; j< col; j++)
            cout << (int)(a[i*col + j]) << " ";

        cout << endl;
    }
}


#endif /* UTILITIES_H_ */
