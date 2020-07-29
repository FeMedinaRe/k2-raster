/*
 * Created by Fernando Silva on 18/02/16.
 *
 * Copyright (C) 2016-current-year, Fernando Silva, all rights reserved.
 *
 * Author's contact: Fernando Silva  <fernando.silva@udc.es>
 * Databases Lab, University of A Coruña. Campus de Elviña s/n. Spain
 *
 * Execute a set of queries of getCellsRangeValue over a k2-raster
 *
 * k2-raster is a compact data structure to represent raster data that
 * uses compressed space and offers indexing capabilities.
 * It uses min/max values for indexing and improving query performance.
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

using namespace std;

int main(int  argc, char ** argv) {

    if (argc != 4) {
        printf("Usage: %s <treapFile> <queriesFile> <check>\n", argv[0]);
        exit(-1);
    }

    bool check = atoi(argv[3]);
    start_timing();
    ifstream inFile(argv[1]);
    if (!inFile.good()) {
        printf("k2-raster file %s unable to open\n", argv[1]);
        exit(-1);
    }
    k2raster_static::K2Raster *raster = k2raster_static::K2Raster::load(inFile);;
    inFile.close();
//    printf("Load RegionTreap %uX%u in %.3f seconds\n", treap->getTreapSize(), treap->getTreapSize(), get_timing()/1000);

    start_timing();
    ifstream queriesFile(argv[2]);
    if (!queriesFile.good()) {
        printf("Query file %s unable to open\n", argv[2]);
        exit(-1);
    }
    query_utils_static::QueryBuilder *queryBuilder = new query_utils_static::QueryBuilder(queriesFile, 0);
    query_utils_static::Query **queries = queryBuilder->getAllQueries(queriesFile);
//    printf("Load %u queries in %.3f seconds\n", queryBuilder->getNumberOfQueries(), get_timing()/1000);

    uint *positions = (uint*)malloc(raster->getNumOfPositions() * sizeof(uint));
    ulong totalCells = 0, totalCellsQuery, totalCellsTest, posTest;
    int value;

    if (check){
        printf("Checking activated....\n");
    }
    start_timing();
    for (uint i = 0; i < queryBuilder->getNumberOfQueries(); i++){
//        queries[i]->printQuery(i + 1);
//        printf("\n");
        totalCellsQuery = raster->getCellsByValue(*queries[i], positions);
        totalCells += totalCellsQuery;

        if (check) {
            totalCellsTest = 0;
            for (ulong r = queries[i]->xini; r <= queries[i]->xend; r++) {
                for (ulong c = queries[i]->yini; c <= queries[i]->yend; c++) {
                    value = raster->getCell(r, c);
                    posTest = r * raster->getDimSize(2) + c;
                    if (value >= queries[i]->valini && value < queries[i]->valend) {
                        totalCellsTest++;
                        bool cellFound = false;
                        for (ulong cell = 0; cell < totalCellsQuery; cell++) {
                            if (positions[cell] == posTest) {
                                cellFound = true;
                                break;
                            }
                        }
                        if (!cellFound) {
                            queries[i]->printQuery(i + 1);
                            printf("\nCell (%u, %u) with value %i not found!!\n", r, c, value);
                            exit(-1);
                        }
                    }
                }
            }

            if (totalCellsTest != totalCellsQuery) {
                queries[i]->printQuery(i + 1);
                printf("\nERROR - Get %lu cells, expected %lu\n", totalCellsQuery, totalCellsTest);
                exit(-1);
            }
        }

//        for (uint c = 0; c < cells; c++) {
//            uint row = positions[c] / treap->getDimSize(2);
//            uint col = positions[c] % treap->getDimSize(2);
//            printf("(%u, %u) -> %u\n", row, col, treap->getCell(row, col));
//        }
//        break;
    }
    double totalTime = get_timing();
    printf("time = %.5f, results = %lu, us/query = %.5f, us/res=%.5f\n",
           totalTime / 1000, totalCells, (totalTime*1000) /queryBuilder->getNumberOfQueries(), (totalTime*1000) / totalCells );

    delete(raster);
}