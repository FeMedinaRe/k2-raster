/*
 * Created by Fernando Silva on 23/06/16.
 *
 * Copyright (C) 2016-current-year, Fernando Silva, all rights reserved.
 *
 * Author's contact: Fernando Silva  <fernando.silva@udc.es>
 * Databases Lab, University of A Coruña. Campus de Elviña s/n. Spain
 *
 * test getCellsRangeValues function over some different datasets
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
    ulong totalCellsQuery = 0,  totalCells = 0;
    ulong numberOfQueries = 0;


    if (check){
        printf("Checking activated....\n");
    }

    double *setTimes = (double *)malloc(maxNRow * maxNCol * sizeof(double));
    uint dataSet = 0;

    for (uint r = 0; r < maxNRow; ++r) {
        for(uint c = 0; c < maxNCol; ++c) {
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

            // Run queries
            uint *positions = (uint*)malloc(raster->getNumOfPositions() * sizeof(uint));
            start_timing();
            for (uint i = 0; i < queryBuilder->getNumberOfQueries(); i++){
                totalCellsQuery = raster->getCellsByValue(*queries[i], positions);
                totalCells += totalCellsQuery;

                // Check query result
                if (check) {
                    ulong totalCellsTest = 0;
                    uint value, posTest;
                    for (uint r = queries[i]->xini; r <= queries[i]->xend; r++) {
                        for (uint c = queries[i]->yini; c <= queries[i]->yend; c++) {
                            value = raster->getCell(r, c);
                            posTest = r * raster->getDimSize(2) + c;
                            if (value >= queries[i]->valini && value < queries[i]->valend) {
                                totalCellsTest++;
                                bool cellFound = false;
                                for (uint cell = 0; cell < totalCellsQuery; cell++) {
                                    if (positions[cell] == posTest) {
                                        cellFound = true;
                                        break;
                                    }
                                }
                                if (!cellFound) {
                                    printf("Cell (%u, %u) with value %i not found!!\n", r, c, value);
                                    exit(-1);
                                }
                            }
                        }
                    }

                    if (totalCellsTest != totalCellsQuery) {
                        printf("ERROR - Get %u cells, expected %u\n", totalCellsQuery, totalCellsTest);
                        exit(-1);
                    }
                }
            }

            setTimes[dataSet] = get_timing();
            totalTime += setTimes[dataSet];
            dataSet++;
            numberOfQueries += queryBuilder->getNumberOfQueries();

            for (uint i = 0; i < queryBuilder->getNumberOfQueries(); i++){
                free(queries[i]);
            }
            free(queries);
            free(positions);
            delete(queryBuilder);
            delete(raster);
        }
    }

    printf("time = %.5f, results = %u, us/query = %.5f, us/res=%.5f\n",
           totalTime / 1000, totalCells, (totalTime*1000) / numberOfQueries, (totalTime*1000) / totalCells );

    double stdev = math_utils_static::calculateStdev(maxNRow * maxNCol, setTimes);
    double mean = math_utils_static::calculateMean(maxNRow * maxNCol, setTimes);
    free(setTimes);
    printf("mean time = %.6f for %u datasets with standard derivation %.6f\n", mean / 1000, maxNRow * maxNCol, stdev/ 1000);
}