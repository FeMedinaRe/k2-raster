/*
 * Created by Fernando Silva on 14/04/16.
 *
 * Copyright (C) 2016-current-year, Fernando Silva, all rights reserved.
 *
 * Author's contact: Fernando Silva  <fernando.silva@udc.es>
 * Databases Lab, University of A Coruña. Campus de Elviña s/n. Spain
 *
 * Execute a set of queries of GetValuesWindow over a k2-raster
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
        printf("File %s unable to open\n", argv[2]);
        exit(-1);
    }
    k2raster_static::K2Raster *raster = k2raster_static::K2Raster::load(inFile);
    inFile.close();
//    printf("Load RegionTreap %uX%u in %.3f seconds\n", treap->getTreapSize(), treap->getTreapSize(), get_timing()/1000);

    start_timing();
    ifstream queriesFile(argv[2]);
    if (!queriesFile.good()) {
        printf("File %s unable to open\n", argv[2]);
        exit(-1);
    }

    query_utils_static::QueryBuilder *queryBuilder = new query_utils_static::QueryBuilder(queriesFile, 0);
    query_utils_static::Query **queries = queryBuilder->getAllQueries(queriesFile);
//    printf("Load %u queries in %.3f seconds\n", queryBuilder->getNumberOfQueries(), get_timing()/1000);

    int val1, val2, pos;
    int *values = (int*)malloc(raster->getNumOfPositions() * sizeof(int));
    ulong totalValues = 0, res;

    if (check) {
        printf("Checking activated....\n");
    }
    start_timing();
    for (uint i = 0; i < queryBuilder->getNumberOfQueries(); i++){
//        queries[i]->printQuery(i + 1);
        res = raster->getValuesWindow(*queries[i], values);

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
//        printf("got %u results\n", res);
        totalValues += res;
    }

    double totalTime = get_timing();
    printf("time = %.6f, results = %lu, us/query = %.6f, us/res=%.6f\n",
           totalTime / 1000, totalValues, (totalTime*1000) /queryBuilder->getNumberOfQueries(), (totalTime*1000) / totalValues );

    delete(raster);
}

