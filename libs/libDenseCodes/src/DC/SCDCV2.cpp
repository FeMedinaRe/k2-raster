/*  
 * Created by Fernando Silva on 12/03/17.
 *
 * Copyright (C) 2017-current-year, Fernando Silva, all rights reserved.
 *
 * 
 * Author's contact: Fernando Silva  <fernando.silva@udc.es>
 * Databases Lab, University of A Coruña. Campus de Elviña s/n. Spain
 *
 *  C++ implementation of (s,c)-Dense Code.
 *  N. Brisaboa, A. Fariña, G. Navarro, and J. Paramá. Lightweight natural language text compression. Information Retrieval. Information Retrieval(10), pp. 1-33, 2007. (online) (doi:10.1007/s10791-006-9001-9)
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include <DC/SCDCV2.h>
#include <malloc.h>
#include <cmath>
#include <map>

namespace densecode_static {

    //**********************************************************************//
    //********************* CONSTRUCTORS ***********************************//
    //**********************************************************************//
    SCDCV2::SCDCV2(uint *data, ulong size, uint s) {

        // Search S and C
        if (s == 0) {
            s = this->findBestS_Secuencial(data, size);
//            printf("Optimal S -> %u || C -> %u\n", s, MAX_SC - s);
        }

        // Set S and C
        this->s = s;
        this->c = MAX_SC - s;

        // Initialize base table
        this->initBaseTable();

        // Encode data
        this->encodeData(data, size);
    }



    SCDCV2::~SCDCV2() {
        free(this->dataEncoded);
        free(this->baseTable);
    }

    //**********************************************************************//
    //************************* ENCODE *************************************//
    //**********************************************************************//

    uint SCDCV2::encodeNumber(uint number, ulong &codeword) {
        uint len = 1;
        uint x = 0;

        if (number >= s) {
            codeword = (number % s);
            x = (number / s) - 1;

            while (x >= c) {
                codeword = codeword << NUM_BITS_MOVE;
                codeword += (ulong((x % c) + s) /*<< NUM_BITS_MOVE * len*/);
                x = (x / c) - 1;
                len++;
            }
            codeword = codeword << NUM_BITS_MOVE;
            codeword += (ulong(x + s) /*<< NUM_BITS_MOVE * len*/);
            len++;
        } else {
            codeword = number;
        }

        return len;
    }

//    void SCDCV2::encodeData(uint *data, uint size) {
//        uint x, number;
//        u_int8_t field;
//        this->dataEncoded = (unsigned char *)malloc(size * sizeof(ulong) * sizeof(unsigned char));
//        this->totalLength = 0;
//        this->size = size;
//        for (uint n = 0; n < size; n++) {
//            number = data[n];
//
//            if (number >= s) {
//                field = (number % s);
//                this->dataEncoded[this->totalLength++] = field;
//                x = (number / s) - 1;
//
//                while (x >= c) {
//                    field = (x % c) + s;
//                    this->dataEncoded[this->totalLength++] = field;
//                    x = (x / c) - 1;
//                }
//                field = (x + s);
//                this->dataEncoded[this->totalLength++] = field;
//            } else {
//                field = number;
//                this->dataEncoded[this->totalLength++] = field;
//            }
//            printf("Number %u || len -> %u\n", n , this->totalLength);
//        }
//        this->dataEncoded = (unsigned char*)realloc(this->dataEncoded, this->totalLength);
//        return;
//
//    }

    void SCDCV2::encodeData(uint *data, uint size) {
        uint len;
        ulong codeword;

        // Initialize data
        this->dataEncoded = (unsigned char *)malloc(size * sizeof(ulong) * sizeof(unsigned char));
        this->totalLength = 0;
//        this->size = size;

        for(uint n = 0; n < size; n++) {
            len = this->encodeNumber(data[n], codeword);
            for (uint b = 0; b < len; b++) {
                this->dataEncoded[this->totalLength++] = reinterpret_cast<unsigned char *>(&codeword)[b];
            }
        }
        return;
    }

    void SCDCV2::initBaseTable() {
        // Initialize base table
        this->baseTable = static_cast<ulong*>(malloc(sizeof(ulong) * sizeof(ulong)));
        this->baseTable[0] = 0;
        for (uint i = 1; i < sizeof(ulong); ++i) {
            this->baseTable[i] = this->baseTable[i - 1] + this->s * (pow(c,(i - 1)));
        }
    }


    //**********************************************************************//
    //************************* DECODE *************************************//
    //**********************************************************************//
//    uint* SCDCV2::decode() {
//        // Malloc final list
//        uint* data = (uint*)malloc(this->size * sizeof(uint));
//
//        // Decode values
//        uint j = 0, pos = 1, index = 0;
//        for (uint i = 0; i < this->totalLength; ++i) {
//            if (u_int8_t(this->dataEncoded[i]) >= s) {
//                j = j * c + (u_int8_t(this->dataEncoded[i]) - s);
//                pos++;
//            } else {
//                j = j * s + u_int8_t(this->dataEncoded[i]);
//                j = j + this->baseTable[pos - 1];
//                data[index] = j;
//                pos = 1;
//                j = 0;
//                index++;
//            }
//
//        }
//        return  data;
//    }

    uint SCDCV2::decodeNext() {
        uint j = 0, len = 0;
        uint value = u_int8_t(this->dataEncoded[this->lastPostion]);
        while (value >= s){
            j = j * c + (value - s);
            this->lastPostion++;
            len++;
            value = u_int8_t(this->dataEncoded[this->lastPostion]);
        }

        // Last chunk (stopper)
        j = j * s + value;
        j = j + this->baseTable[len];
        this->lastPostion++;
        return j;
    }

    void SCDCV2::resetNext() {
        this->lastPostion = 0;
    }

    //**********************************************************************//
    //*************************** SIZE *************************************//
    //**********************************************************************//
    size_t SCDCV2::getTotalSize() const {
        size_t totalSize =  sizeof(uint) +                              // s
//                            sizeof(ulong) +                             // size
                            sizeof(unsigned char) * this->totalLength + // data
                            sizeof(ulong);                              // totalLength

        return totalSize;
    }

    //**********************************************************************//
    //*************************** FILE *************************************//
    //**********************************************************************//
    void SCDCV2::save(std::ofstream &fp) const {
        fp.write((char*)&this->s,sizeof(ushort));                                 // s
//        fp.write((char*)&this->size,sizeof(ulong));                             // size
        fp.write((char*)&this->totalLength,sizeof(ulong));                      // totalLength
        fp.write((char*)this->dataEncoded,sizeof(unsigned char) * this->totalLength);    // data
    }

    SCDCV2* SCDCV2::load(std::ifstream &in) {
        SCDCV2 *scdc = new SCDCV2();
        in.read((char*)&scdc->s,sizeof(ushort));
        scdc->c = scdc->MAX_SC - scdc->s;
//        in.read((char*)&scdc->size,sizeof(ulong));
        in.read((char*)&scdc->totalLength,sizeof(ulong));
        scdc->dataEncoded = (unsigned char*)malloc(scdc->totalLength * sizeof(unsigned char));
        in.read((char*)scdc->dataEncoded, scdc->totalLength*sizeof(unsigned char));

        scdc->initBaseTable();
        return scdc;
    }

    //**********************************************************************//
    //********************** FIND OPTIMAL S, C *****************************//
    //**********************************************************************//
    uint SCDCV2::findBestS_Secuencial(uint* data, uint size) {
        uint maxNumber = 0;
        std::map<uint, uint> freq;
        std::map<uint, uint>::iterator it;

        // Obtain frequencies
        for (uint n = 0; n < size; n++) {
            if (data[n] > maxNumber) {
                maxNumber = data[n];
            }

            it = freq.find(data[n]);
            if (it == freq.end()) {
                // New key
                freq.insert (std::pair<uint, ulong>(data[n],1));
            } else {
                // Add +1 to frequency
                it->second++;
            }
        }

        // Create accumulator of frequencies
        uint * freqAcum = (uint*)malloc((maxNumber+1) * sizeof(uint));
        freqAcum[0] = 0;
        for (uint i = 1; i <= maxNumber; i++) {
            it = freq.find(i);
            freqAcum[i] = freqAcum[i-1];
            if (it != freq.end()) {
                freqAcum[i] += it->second;
            }
        }

        // Searching optimal S, C
        uint bestS=1;
        {
            uint  sizeBest = UINT32_MAX, newSize;
            for (uint i = 1; i < MAX_SC; i++) {
                newSize = computeSizeS(i, MAX_SC-i, freqAcum, maxNumber, sizeBest);
                if (newSize < sizeBest) {
                    sizeBest = newSize;
                    bestS = i;
                }
            }
        }
        free(freqAcum);
        return bestS;
    }

    ulong SCDCV2::computeSizeS(uint s, uint c, uint *freqAcum, ulong size, ulong sizeBest) {
        ulong  totFreq, power, left,total;

        totFreq = freqAcum[size];
        total = 0;
        power = s;
        left = 0;

        while ((left < size) && (total < sizeBest)) {
            total += totFreq - freqAcum[left];
            left += power;
            power *=c;
        }

        return total;
    }


/*------------------------------------------------------------------
Given s, c, and a list of accumulated frequencies, obtains the size of the compressed file (in bytes)
 ---------------------------------------------------------------- */
//    unsigned long int computeSizeS(int s, int c, long *freqList, unsigned long n, unsigned long SIZEBEST)
//    {
//        unsigned long int totFreq,potencia,Left,total;
//        unsigned long int i;
//
//        totFreq = freqList[n];
//        total = 0;
//        potencia =s;
//        Left = 0;
//        //      i=0;   -- DEBUGGING
//        while ((Left <n) && (total < SIZEBEST)) {
//            total += totFreq - freqList[Left];
//            Left  += potencia;
//            potencia *= c;
//            //      i++;
//            /*              if (total>=SIZEBEST){
//                        //      printf("breaking K=%d ",i);
//                                break;
//                        }
//            //              fprintf(stderr,"s=%d K=%d ",s,i);
//            */
//        }
//        return total;
//    }



} // END namespace densecode_static



