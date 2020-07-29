/*
 * Created by Fernando Silva on 29/11/16.
 *
 * Copyright (C) 2016-current-year, Fernando Silva, all rights reserved.
 *
 * Author's contact: Fernando Silva  <fernando.silva@udc.es>
 * Databases Lab, University of A Coruña. Campus de Elviña s/n. Spain
 *
 *  TEST of (s,c)-Dense Code. 
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


#include <sys/types.h>
#include <malloc.h>
#include <cstdio>
#include <DC/SCDC.h>

#include <stdlib.h>     /* srand, rand */

using namespace std;

int main(int  argc, char ** argv) {

    densecode_static::sCodeword cw;
    unsigned char * cwT;
    uint nValues = 100;
    uint *values = (uint*)malloc(nValues * sizeof(uint));
    printf("STEP 0 - initializing array with %u random values\n", nValues);
    for (uint i = 0; i < nValues; i++) {
        values[i] = rand() % 10000;
    }

    unsigned char dataEncoded[nValues * sizeof(uint)];
    uint position = 0;

    densecode_static::SCDC *scdc = new densecode_static::SCDC(100);

    // 1 - Encode numbers
    printf("STEP 1 - Encoding %u numbers\n", nValues);
    for (uint i = 0; i < nValues; i++) {
        cw = scdc->encode(values[i]);
        cwT = reinterpret_cast<unsigned  char*>(&(cw.codeword));
        for (uint l = 0; l < cw.len; l++) {
            dataEncoded[position++] = cwT[l];
        }
    }
    printf("Compressed %u numbers in %u bytes\n", nValues, position);


    // 2- Decoded numbers and check values
    printf("STEP 2 - Decoding %u numbers\n", nValues);
    uint *newValues = (uint*)malloc(nValues * sizeof(uint));
    scdc->decode(dataEncoded, position, newValues);

    // 3 - Checking values
    printf("STEP 3 - Cheking values\n");
    for (uint i = 0; i < nValues; i++) {
        if (values[i] != newValues[i]) {
            printf("Error found at position %u: got %u and expected %u\n", i, newValues[i], values[i]);
            exit(-1);
        }
    }

    printf("ALL OK!!!!!\n");
}

