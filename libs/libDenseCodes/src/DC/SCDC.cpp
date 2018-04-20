/*
 * Created by Fernando Silva on 29/11/16.
 *
 * Copyright (C) 2016-current-year, Fernando Silva, all rights reserved.
 *
 * Author's contact: Fernando Silva  <fernando.silva@udc.es>
 * Databases Lab, University of A Coruña. Campus de Elviña s/n. Spain
 *
 *  C++ implementation of (s,c)-Dense Code. 
 *  N. Brisaboa, A. Fariña, G. Navarro, and J. Paramá. Lightweight natural language text compression. Information Retrieval. Information Retrieval(10), pp. 1-33, 2007. (online) (doi:10.1007/s10791-006-9001-9)
 *  
 *  For more information http://vios.dc.fi.udc.es/codes/semistatic.html
 *  Original code http://vios.dc.fi.udc.es/codes/files/scdc.tar.gz
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

#include <DC/SCDC.h>
#include <malloc.h>
#include <cmath>

namespace densecode_static {

    /**********************************************/
    /*              CODEWORD CLASS                */
    /**********************************************/
    SCDC::SCDC(uint s) {

        // Set S and C
        this->s = s;
        this->c = MAX_SC - s;

        // Initialize base table
        this->baseTable = static_cast<ulong*>(malloc(sizeof(ulong) * sizeof(ulong)));
        this->baseTable[0] = 0;
        for (uint i = 1; i < sizeof(ulong); ++i) {
            this->baseTable[i] = this->baseTable[i - 1] + this->s * (pow(c,(i - 1)));
//            printf("TableBase -> %lu\n", baseTable[i]);
        }
    }

    SCDC::~SCDC() {
        free(this->baseTable);
    }

    sCodeword SCDC::encode(uint number) {

        uint len = 1;
        uint x = 0;
        sCodeword codeword;
        codeword.len = 0;

        if (number >= s) {
            codeword.codeword = (number % s);
            x = (number / s) - 1;

            while (x >= c) {
                codeword.codeword = codeword.codeword << NUM_BITS_MOVE;
                codeword.codeword += (ulong((x % c) + s) /*<< NUM_BITS_MOVE * len*/);
                x = (x / c) - 1;
                len++;
            }
            codeword.codeword = codeword.codeword << NUM_BITS_MOVE;
            codeword.codeword += (ulong(x + s) /*<< NUM_BITS_MOVE * len*/);
            len++;
        } else {
            codeword.codeword = number;
        }

        codeword.len = len;
        return codeword;
    }

    void SCDC::encode(uint number, sCodeword *codeword) {

        uint len = 1;
        uint x = 0;
        codeword->len = 0;

        if (number >= s) {
            codeword->codeword = (number % s);
            x = (number / s) - 1;

            while (x >= c) {
                codeword->codeword = codeword->codeword << NUM_BITS_MOVE;
                codeword->codeword += (ulong((x % c) + s) /*<< NUM_BITS_MOVE * len*/);
                x = (x / c) - 1;
                len++;
            }
            codeword->codeword = codeword->codeword << NUM_BITS_MOVE;
            codeword->codeword += (ulong(x + s) /*<< NUM_BITS_MOVE * len*/);
            len++;
        } else {
            codeword->codeword = number;
        }

        codeword->len = len;
        return;
    }

    unsigned char* SCDC::decodeNext(unsigned char *buffer, uint *number) {
        uint i = 0, j = 0;
        uint value = static_cast<uint >(buffer[i]);
        while (value >= s){
            j = j * c + (value - s);
            i++;
            value = static_cast<uint >(buffer[i]);
        }

        // Last chunk (stopper)
        j = j * s + value;
        j = j + this->baseTable[i];
        *number = j;
        return buffer + i + 1; // TODO is valid?
    }

    uint* SCDC::decode(unsigned  char *buffer, uint bufferSize, uint listSize) {
        uint *numbers = reinterpret_cast<uint*>(malloc(listSize * sizeof(uint)));

        uint j = 0, pos = 1, index = 0;
        for (uint i = 0; i < bufferSize; ++i) {
            if (int(buffer[i]) >= s) {
                j = j * c + (int(buffer[i]) - s);
                pos++;
            } else {
                j = j * s + int(buffer[i]);
                j = j + this->baseTable[pos - 1];
                numbers[index] = j;
                pos = 1;
                j = 0;
                index++;

                if (index == listSize){
                    return  numbers;
                }
            }

        }

        return  numbers;
    }

    void SCDC::decode(unsigned  char *buffer, uint bufferSize, uint *numbers) {

        uint j = 0, pos = 1, index = 0;
        for (uint i = 0; i < bufferSize; ++i) {
            if (int(buffer[i]) >= s) {
                j = j * c + (int(buffer[i]) - s);
                pos++;
            } else {
                j = j * s + int(buffer[i]);
                j = j + this->baseTable[pos - 1];
                numbers[index] = j;
                pos = 1;
                j = 0;
                index++;
            }
        }

        return;
    }
}