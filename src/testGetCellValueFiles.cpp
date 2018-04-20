/*
 * Created by Fernando Silva on 23/06/16.
 *
 * Copyright (C) 2016-current-year, Fernando Silva, all rights reserved.
 *
 * Author's contact: Fernando Silva  <fernando.silva@udc.es>
 * Databases Lab, University of A Coruña. Campus de Elviña s/n. Spain
 *
 * test getCell function over some different datasets
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
#include <includes/utils/timing.h>
#include <queries/QueryBuilder.h>
#include <k2raster/k2-raster.h>
#include <include/math/statistics.h>

using namespace std;

int main(int  argc, char ** argv) {

    if (argc != 7) {
        printf("Usage: %s <k2rasterFiles> <k2rasterFilesSuffix> <maxNRow> <maxNCol> <queriesFile> <queriesFileSuffix>\n", argv[0]);
        exit(-1);
    }

    std::string k2rasterFiles = argv[1];
    std::string k2rasterFilesSuffix = argv[2];
    uint maxNRow = atoi(argv[3]);
    uint maxNCol = atoi(argv[4]);
    std::string queriesFile = argv[5];
    std::string queriesFileSuffix = argv[6];
    double totalTime = 0;
    ulong totalValue = 0;
    ulong numberOfQueries = 0;

    double *setTimes = (double *)malloc(maxNRow * maxNCol * sizeof(double));
    uint dataSet = 0;

    for (uint r = 0; r < maxNRow; ++r) {
        for(uint c = 0; c < maxNCol; ++c) {
            // Load File
            std::string filename = k2rasterFiles + "_" + to_string(r) + to_string(c) + "." + k2rasterFilesSuffix; // TODO add special separator between r and c
//            printf("Test file %s\n", filename.data());
            ifstream inFile(filename);
            if (!inFile.good()) {
                printf("File %s unable to open\n", filename.data());
                exit(-1);
            }

            // Load k2-raster
            k2raster_static::K2Raster *raster = k2raster_static::K2Raster::load(inFile);;
            inFile.close();

            // Load queries
            filename = queriesFile + "_" + to_string(r) + to_string(c) + "-" + queriesFileSuffix;
            ifstream queriesFile(filename);
            if (!queriesFile.good()) {
                printf("File %s unable to open\n", filename.data());
                exit(-1);
            }
            query_utils_static::QueryBuilder *queryBuilder = new query_utils_static::QueryBuilder(queriesFile, true, 0);
            query_utils_static::Query **queries = queryBuilder->getAllQueries(queriesFile);


            // Run queries
            start_timing();
            for (uint i = 0; i < queryBuilder->getNumberOfQueries(); i++){
                totalValue += raster->getCell(queries[i]->xini, queries[i]->yini);
            }
            double queryTime = get_timing();
            //if (queryTime != 0) {
                setTimes[dataSet] = queryTime/ 1000.;
                totalTime += queryTime;
                dataSet++;
                numberOfQueries += queryBuilder->getNumberOfQueries();
            //}

            for (uint i = 0; i < queryBuilder->getNumberOfQueries(); i++){
                free(queries[i]);
            }
            free(queries);
            delete(queryBuilder);
            delete(raster);
        }
    }
    printf("Time spent: %.6f seconds, %.6f microseconds/query, totalValue -> %lu\n",
           totalTime / 1000, (totalTime*1000) / numberOfQueries, totalValue);

    double stdev = math_utils_static::calculateStdev(dataSet, setTimes);
    double mean = math_utils_static::calculateMean(dataSet, setTimes);
    //for (uint i = 0 ; i < maxNRow * maxNCol; ++i) {
    //    printf("Set %u -> %.6f \n", i, setTimes[i] / 1000.);
    //}
    free(setTimes);

    printf("mean time = %.6f for %u datasets with standard derivation %.6f\n", mean, dataSet, stdev);
}