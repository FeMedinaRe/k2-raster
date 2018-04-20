/*
 * Created by Fernando Silva on 13/03/17.
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
#include <DC/SCDCV2.h>

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

using namespace std;

int main(int  argc, char ** argv) {

    uint nValues = 1000;
    uint *values = (uint*)malloc(nValues * sizeof(uint));
    printf("STEP 0 - initializing array with %u random values\n", nValues);
//    srand(time(NULL));
    for (uint i = 0; i < nValues; i++) {
        values[i] = rand() % 1000000;
    }

    // 1 - Encode numbers
    printf("STEP 1 - Encoding %u numbers\n", nValues);
    densecode_static::SCDCV2 *scdc = new densecode_static::SCDCV2(values, nValues);
    printf("Compressed %u numbers in %lu bytes\n", nValues, scdc->getTotalSize());


    // Save and load in a file
    ofstream outFile("./test.scdc");
    scdc->save(outFile);
    outFile.close();
    delete(scdc);
    ifstream inputFile("./test.scdc");
    scdc = densecode_static::SCDCV2::load(inputFile);
    inputFile.close();

    // 2- Decoded numbers and check values
//    printf("STEP 2 - Decoding %u numbers\n", nValues);
//    uint *newValues = scdc->decode();
//
//    // 3 - Checking values
//    printf("STEP 3 - Cheking values\n");
//    for (uint i = 0; i < nValues; i++) {
//        if (values[i] != newValues[i]) {
//            printf("Error found at position %u: got %u and expected %u\n", i, newValues[i], values[i]);
//            exit(-1);
//        }
//    }
//    free(newValues);

    // 4 - Checking decodeNext
    printf("STEP 4 - Cheking decodeNext\n");
    uint val;
    for (uint i = 0; i < nValues; i++) {
        val = scdc->decodeNext();
        if (values[i] != val) {
            printf("Error found at position %u: got %u and expected %u\n", i, val, values[i]);
            exit(-1);
        }
    }

    // 4 - Free memory
    printf("STEP 4 - Releasing memory\n");
    free(values);
    delete(scdc);

    printf("ALL OK!!!!!\n");
}

