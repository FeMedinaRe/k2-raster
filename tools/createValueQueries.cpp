/*  
 * Created by Fernando Silva on 7/09/17.
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

#include <iostream>
#include <queries/queriesUtils.h>
#include <k2raster/k2-raster.h>

using namespace std;

int main(int  argc, char ** argv) {

    if (argc != 5) {
        printf("Usage: %s <numQueries> <k2-rasterFile> <outFile> <valuesLimit>\n", argv[0]);
        exit(-1);
    }

    uint numQueries =  atoi(argv[1]);
    uint valuesLimit =  atoi(argv[4]);


    // Load k2-raster
    ifstream inFile(argv[2]);
    if (!inFile.good()) {
        printf("File %s unable to open\n", argv[1]);
        exit(-1);
    }

    k2raster_static::K2Raster *raster = k2raster_static::K2Raster::load(inFile);
    inFile.close();
    if (raster == nullptr) {
        printf("Error loading k2-raster from file %s\n", argv[1]);
        exit(-1);
    }


    // Open output file
    ofstream outFile(argv[3]);
    if (!outFile.good()){
        printf("File %s unable to create\n", argv[3]);
        exit(-1);
    }

    printf("Creating %u queries (MinValue -> %u || MaxValue -> %u) with limit %u\n", numQueries, raster->getMinValueNoZero(), raster->getMaxValue(), valuesLimit);
    query_utils_static::getValuesQueries(numQueries, raster->getMinValueNoZero(), raster->getMaxValue(), valuesLimit, outFile);
    outFile.close();
}


