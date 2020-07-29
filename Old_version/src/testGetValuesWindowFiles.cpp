/*
 * Created by Fernando Silva on 23/06/16.
 *
 * Copyright (C) 2016-current-year, Fernando Silva, all rights reserved.
 *
 * Author's contact: Fernando Silva  <fernando.silva@udc.es>
 * Databases Lab, University of A Coruña. Campus de Elviña s/n. Spain
 *
 * test getValuesWindow function over several different datasets
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

    if (argc != 8) {
        printf("Usage: %s <treapFiles> <treapFilesSubfix> <maxNRow> <maxNCol> <queriesFile> <queriesFileSubfix> <check>\n", argv[0]);
        exit(-1);
    }

    std::string treapFiles = argv[1];
    std::string treapFilesSubfix = argv[2];
    uint maxNRow = atoi(argv[3]);
    uint maxNCol = atoi(argv[4]);
    std::string queriesFile = argv[5];
    std::string queriesFileSubfix = argv[6];
    bool check = atoi(argv[7]);
    double totalTime = 0;
    ulong totalCellsQuery = 0, totalCells = 0;
    ulong numberOfQueries = 0;


    if (check){
        printf("Checking activated....\n");
    }

    double *setTimes = (double *)malloc(maxNRow * maxNCol * sizeof(double));
    uint dataSet = 0;

    for (uint r = 0; r < maxNRow; ++r) {
        for (uint c = 0; c < maxNCol; ++c) {
            // Load File
            std::string filename = treapFiles + "_" + to_string(r) + to_string(c) + "." + treapFilesSubfix;
//            printf("Test file %s\n", filename.data());
            ifstream inFile(filename);
            if (!inFile.good()) {
                printf("File %s unable to open\n", filename.data());
                exit(-1);
            }

            // Load k2-raster
            k2raster_static::K2Raster *raster = k2raster_static::K2Raster::load(inFile);
            inFile.close();

            // Load queries
            filename = queriesFile + "_" + to_string(r) + to_string(c) + "-" + queriesFileSubfix;
            ifstream queriesFile(filename);
            if (!queriesFile.good()) {
                printf("File %s unable to open\n", filename.data());
                exit(-1);
            }

            query_utils_static::QueryBuilder *queryBuilder = new query_utils_static::QueryBuilder(queriesFile, 0);
            query_utils_static::Query **queries = queryBuilder->getAllQueries(queriesFile);


            int val1, val2, pos;
            int *values = (int *) malloc(raster->getNumOfPositions() * sizeof(int));

            uint cellsT = 0;
            start_timing();
            for (uint i = 0; i < queryBuilder->getNumberOfQueries(); i++) {
                totalCells = raster->getValuesWindow(*queries[i], values);
                if (check) {
                    for (uint x = queries[i]->xini; x <= queries[i]->xend; x++) {
                        for (uint y = queries[i]->yini; y <= queries[i]->yend; y++) {
                            pos = (x - queries[i]->xini) * ((queries[i]->yend - queries[i]->yini) + 1) +
                                  (y - queries[i]->yini);
                            val1 = values[pos];
                            val2 = raster->getCell(x, y);
                            if (val1 != val2) {
                                printf("ERROR -- (%u, %u [%u]) Values is %i and expect %i\n", x, y, pos, val1, val2);
                                exit(-1);
                            }
                        }
                    }
                }
                totalCellsQuery += totalCells;
                cellsT+=totalCells;
            }

            setTimes[dataSet] = get_timing();
            totalTime += setTimes[dataSet];
            dataSet++;
            numberOfQueries += queryBuilder->getNumberOfQueries();

//            printf(" TotalCells -> %u\n", cellsT);

            for (uint i = 0; i < queryBuilder->getNumberOfQueries(); i++){
                free(queries[i]);
            }
            free(queries);
            free(values);
            delete(queryBuilder);
            delete(raster);
        }
    }

    printf("time = %.6f, results = %u, us/query = %.6f, us/res=%.6f\n",
           totalTime / 1000, totalCellsQuery, (totalTime*1000) / numberOfQueries, (totalTime*1000) / totalCellsQuery );

    double stdev = math_utils_static::calculateStdev(maxNRow * maxNCol, setTimes);
    double mean = math_utils_static::calculateMean(maxNRow * maxNCol, setTimes);
    free(setTimes);
    printf("mean time = %.6f for %u datasets with standard derivation %.6f\n", mean / 1000, maxNRow * maxNCol, stdev/ 1000);
}

