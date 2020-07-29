/*  
 * Created by Fernando Silva on 13/03/18.
 *
 * Copyright (C) 2018-current-year, Fernando Silva, all rights reserved.
 *
 * 
 * Author's contact: Fernando Silva  <fernando.silva@udc.es>
 * Databases Lab, University of A Coruña. Campus de Elviña s/n. Spain
 *
 * DESCRIPTION
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

#include <fstream>
#include <includes/utils/cppUtils.h>
#include <RePair.h>
#include <includes/direct_access/DAC.h>

using namespace std;

int main(int  argc, char ** argv) {

    // Reads file
    ifstream inputFile("last_values.bin");
    if (!inputFile.good()) {
        printf("File %s unable to open\n", argv[1]);
        exit(-1);
    }

    uint matrixSize = 16;
    uint numOfValues = cds_utils::loadValue<uint>(inputFile);
    uint elements = numOfValues/matrixSize;
    uint* plainValues = cds_utils::loadValue<uint>(inputFile, numOfValues);


    uint maxValue = 0;
    for (uint i=0; i<numOfValues; i++) {
        if(plainValues[i] > maxValue) {
            maxValue = plainValues[i];
        }
    }

    printf("Values -> %u || elements -> %u || Max value -> %u\n", numOfValues, elements, maxValue);


    int *dict = new int[numOfValues + elements];
    uint *cdict;
    int maxRule = 0;
    RePair *rp;
    {
        uchar *strCurrent = NULL;
        uchar maxChar = 250;
        ulong lenCurrent = 0, processed = 0;

        for (uint e=0; e<elements; e++) {
            for (uint i = 0; i < matrixSize; i++) {
                dict[processed] = plainValues[e * matrixSize + i];
                processed++;
            }
            dict[processed] = 0;
            processed++;

        }
        maxChar++;
        printf("Creating RePair structure with %u elements\n", processed);
        rp = new RePair(dict, processed, maxChar);
        printf("RePair created with %u elements\n", processed);

        cdict = new uint[processed];
        ulong io = 0, ic = 0, strings = 0;
        uint maxseq = 0, currentseq = 0;
        int minRule = INT32_MAX;

        while (io<processed) {
            if (dict[io] >= 0) {
                if (dict[io] != 0) {
                    if (dict[io] > maxRule) maxRule = dict[io];
                    if (dict[io] < minRule) minRule = dict[io];
                    cdict[ic] = dict[io];
                    ic++;
                    currentseq++;
                }
                io++;
            } else {
                if (io < processed) io = -(dict[io]+1);
            }
        }
        maxRule++;
        printf("Strings -> %u || Rules -> %u  || Terminals -> %u || Rules (%u, %u)\n",
               strings, rp->getNumOfRules(), rp->getNumOfTernimals(), minRule, maxRule);
        delete []dict;

        ofstream outFile("RePair-values.rp");
        rp->save(outFile);
        DAC *dac = new cds_static::DAC(cdict, ic);
        dac->save(outFile);
        outFile.close();        // Close

    } // END BLOCK Create Repair rules


}