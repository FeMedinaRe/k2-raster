/*  
 * Created by Fernando Silva on 24/10/17.
 *
 * Copyright (C) 2017-current-year, Fernando Silva, all rights reserved.
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
#include <values.h>

using namespace std;

int main(int  argc, char ** argv) {

    if (argc != 7) {
        printf("Usage: %s <inputFile> <sizeX> <sizeY> <outputFile> <minValue> <maxValue>\n", argv[0]);
        exit(-1);
    }

    uint sizeX = atoi(argv[2]);
    uint sizeY = atoi(argv[3]);
    int minValue = atoi(argv[5]);
    int maxValue = atoi(argv[6]);

    /*********************/
    /* Read input data   */
    /*********************/
    // Reads file
    ifstream inputFile(argv[1]);
    if (!inputFile.good()) {
        printf("File %s unable to open\n", argv[1]);
        exit(-1);
    }

    // Allocs memory
    int minValueInput = MAXINT, maxValueInput = 0;
    int **data = (int **) malloc(sizeX * sizeof(int *));
    for (uint i=0; i<sizeX; i++){
        data[i] = (int *) malloc(sizeY * sizeof(int));
        for (uint j=0; j<sizeY; j++){
            data[i][j] = cds_utils::loadValue<int>(inputFile);
            if (data[i][j] < minValueInput) {
                minValueInput = data[i][j];
            }
            if (data[i][j] > maxValueInput) {
                maxValueInput = data[i][j];
            }
        }
    }
    inputFile.close();
    printf("Read %u (%ux%u) values with minValue = %i and maxValue = %i\n", sizeX*sizeY, sizeX, sizeY, minValueInput, maxValueInput);

    /*********************/
    /* Write output data */
    /*********************/
    // Open output file
    ofstream outFile(argv[4]);
    if (!outFile.good()){
        printf("File %s unable to create\n", argv[4]);
        exit(-1);
    }

    int value;
    int minValueOutput = MAXINT, maxValueOutput = 0;
    bool found= false;
    for (uint i = 0; i < sizeX; i++) {
        for (uint j = 0; j < sizeY; j++) {
            value = minValue + (((double)(data[i][j] - minValueInput) / (double)(1 + maxValueInput- minValueInput)) * (1 + maxValue - minValue));
            cds_utils::saveValue(outFile, value);
            if (value < minValueOutput) {
                minValueOutput = value;
            }
            if (value > maxValueOutput) {
                maxValueOutput = value;
            }
        }
    }
    outFile.close();
    printf("Written %u (%ux%u) values with minValue = %i (%i) and maxValue = %i (%i)\n", sizeX*sizeY, sizeX, sizeY,
           minValueOutput, minValue, maxValueOutput, maxValue);
}