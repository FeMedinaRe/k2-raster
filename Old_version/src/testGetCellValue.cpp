/*
 * Created by Fernando Silva on 19/02/16.
 *
 * Copyright (C) 2016-current-year, Fernando Silva, all rights reserved.
 *
 * Author's contact: Fernando Silva  <fernando.silva@udc.es>
 * Databases Lab, University of A Coruña. Campus de Elviña s/n. Spain
 *
 * Execute a set of queries of GetCellValue over a k2-raster
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

    if (argc != 3) {
        printf("Usage: %s <treapFile> <queriesFile>\n", argv[0]);
        exit(-1);
    }

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
    query_utils_static::QueryBuilder *queryBuilder = new query_utils_static::QueryBuilder(queriesFile, true, 0);
    query_utils_static::Query **queries = queryBuilder->getAllQueries(queriesFile);
//    printf("Load %u queries in %.3f seconds\n", queryBuilder->getNumberOfQueries(), get_timing()/1000);

    ulong totalValue = 0/*, value*/;
    start_timing();
    for (uint i = 0; i < queryBuilder->getNumberOfQueries(); i++){
        //queries[i]->printQuery(i, true);
        totalValue += raster->getCell(queries[i]->xini, queries[i]->yini);
        //printf(" Val -> %u\n", value);
        //totalValue += value;
    }
    double totalTime = get_timing();
    printf("Time spent: %.6f seconds, %.6f microseconds/query, totalValue -> %lu\n",
           totalTime / 1000, (totalTime*1000) /queryBuilder->getNumberOfQueries(), totalValue);

    for (uint i = 0; i < queryBuilder->getNumberOfQueries(); i++){
        free(queries[i]);
    }
    free(queries);
    delete(queryBuilder);
    delete(raster);
}