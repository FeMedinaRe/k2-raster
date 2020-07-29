/*  
 * Created by Fernando Silva on 20/12/16.
 *
 * Copyright (C) 2016-current-year, Fernando Silva, all rights reserved.
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


// Own header
#include <k2raster/k2-raster.h>
#include <rtree/SpatialJoin/QueryStrategy/GetRoot.h>
#include <rtree/SpatialJoin/QueryStrategy/GetObjectsQueryStrategy.h>
#include <rtree/SpatialJoin/QueryStrategy/GetObjectsQueryStrategyQueue.h>
#include <rtree/SpatialJoin/QueryStrategy/GetObjectsQueryStrategyQueueV2.h>
#include <rtree/SpatialJoin/QueryStrategy/GetObjectsQueryStrategyTree.h>
#include <rtree/SpatialJoin/QueryStrategy/GetObjectsQueryStrategyTreeV2.h>
#include <rtree/SpatialJoin/QueryStrategy/GetObjectsQueryStrategyTreeV3.h>
#include <rtree/SpatialJoin/QueryStrategy/GetObjectsQueryStrategyLeaves.h>
#include <rtree/SpatialJoin/Visitor/GetObjectsVisitorCells.h>
// Libraries
#include <queries/QueryBuilderValues.h>
#include <includes/utils/timing.h>

int searchingByCells(ISpatialIndex *tree, k2raster_static::K2Raster *raster,
                     query_utils_static::QueryBuilder *queryBuilder, query_utils_static::Query **queries) {
    k2raster_rtree_static::GetRoot *qer = nullptr;
    qer = new k2raster_rtree_static::GetRoot();
    tree->queryStrategy(*qer);

    double sizeX, sizeY;
    double plow[2], phigh[2];
    ulong totalMBRs = 0;
    ulong totalAllCells = 0;
    uint *positions = (uint *) malloc(raster->getDimSize(1) * raster->getDimSize(2) * sizeof(uint));

    k2raster_rtree_static::GetObjectsVisitorCells vis =
            k2raster_rtree_static::GetObjectsVisitorCells(qer->getPlow(), qer->getPhigh(),
                                                          raster->getDimSize(1), raster->getDimSize(2));
//            printf("Executing %u queries.....\n", queryBuilder->getNumberOfQueries());

    start_timing();

    for (uint i = 0; i < queryBuilder->getNumberOfQueries(); i++) {
        queries[i]->xini = 0;
        queries[i]->xend = raster->getDimSize(1) - 1;
        queries[i]->yini = 0;
        queries[i]->yend = raster->getDimSize(2) - 1;
//        queries[i]->printQuery(i);
//        printf("\n");
        ulong totalCells = raster->getCellsByValue(*(queries[i]), positions);
        tree->queryStrategy(*qer);
        sizeX = qer->getPhigh()[0] - qer->getPlow()[0];
        sizeY = qer->getPhigh()[1] - qer->getPlow()[1];

//        printf("Searching %lu cells\n", totalCells);
        vis.clear();
        for (ulong c = 0; c < totalCells; c++) {
            double x = positions[c] / raster->getDimSize(2);
            double y = positions[c] % raster->getDimSize(2);
            plow[0] = phigh[0] = ((x / raster->getDimSize(1)) * sizeX) + qer->getPlow()[0];
            plow[1] = phigh[1] = ((y / raster->getDimSize(2)) * sizeY) + qer->getPlow()[1];
            Region r = Region(plow, phigh, 2);
            tree->intersectsWithQuery(r, vis);
            totalMBRs += vis.getResult().size();
            totalAllCells += vis.getAllMBRCells();
            if (c % 1000000 == 0) {
                printf("%lu of %lu\n", c, totalCells);
            }
        }
    }
    vis.clear();
    double totalTime = get_timing();
    printf("Time spent: %.6f seconds, %.6f seconds/query, MBRs -> %lu (%lu cell) %.6f microseconds/MBR, %.6f nanoseconds/cell\n",
           totalTime / 1000, (totalTime / 1000) / queryBuilder->getNumberOfQueries(),
           totalMBRs, totalAllCells,
           (totalTime * 1000) / queryBuilder->getNumberOfQueries() / totalMBRs,
           (totalTime * 1000000) / queryBuilder->getNumberOfQueries() / totalAllCells);

    return totalAllCells;
}

int main(int argc, char **argv) {

    if (argc < 5) {
        printf("Usage: %s <k2rasterFile> <rtreeFile> <queryFile> <allCells> [type]\n", argv[0]);
        exit(-1);
    }

    // Load k2-raster
    ifstream inFile(argv[1]);
    if (!inFile.good()) {
        printf("File %s unable to open\n", argv[1]);
        exit(-1);
    }
//    printf("Loading k2-rater from %s\n", argv[1]);
    k2raster_static::K2Raster *raster = k2raster_static::K2Raster::load(inFile);
    inFile.close();
    if (raster == nullptr) {
        printf("Error loading k2-raster from file %s\n", argv[1]);
        exit(-1);
    }

    // Load RTree
    string baseName = argv[2];
//    printf("Loading r-tree from %s\n", baseName.data());
    IStorageManager *diskfile = StorageManager::loadDiskStorageManager(baseName);
    StorageManager::IBuffer *file = StorageManager::createNewRandomEvictionsBuffer(*diskfile, 10, false);
    ISpatialIndex *tree = RTree::loadRTree(*file, 1);

    // Load queries
    ifstream queriesFile(argv[3]);
    if (!queriesFile.good()) {
        printf("File %s unable to open\n", argv[3]);
        exit(-1);
    }
//    printf("Loading queries from %s\n", argv[3]);
    query_utils_static::QueryBuilder *queryBuilder = new query_utils_static::QueryBuilderValues(queriesFile, 0);
    query_utils_static::Query **queries = queryBuilder->getAllQueries(queriesFile);

    // Init query strategy
    bool allCells = atoi(argv[4]);
    uint strategyType = argc == 6 ? atoi(argv[5]) : 0;
    k2raster_rtree_static::GetObjectsQueryStrategy *qe = nullptr;
    switch (strategyType) {
        case STRATEGYQUEUE:
//            printf("Using a queue....\n");
            qe = new k2raster_rtree_static::GetObjectsQueryStrategyQueue();
            break;
        case STRATEGYQUEUEV2:
//            printf("Using a queue V2....\n");
            qe = new k2raster_rtree_static::GetObjectsQueryStrategyQueueV2();
            break;
        case STRATEGYTREE:
//            printf("Using a tree.....\n");
            qe = new k2raster_rtree_static::GetObjectsQueryStrategyTree();
            break;
        case STRATEGYTREEV2:
//            printf("Using a tree V2.....\n");
            qe = new k2raster_rtree_static::GetObjectsQueryStrategyTreeV2();
            break;
        case STRATEGYTREEV3: {
//            printf("Using a tree V3.....\n");
            qe = new k2raster_rtree_static::GetObjectsQueryStrategyTreeV3();
        }
        case STRATEGYLEAVES:
//            printf("Using leaves.....\n");
            qe = new k2raster_rtree_static::GetObjectsQueryStrategyLeaves();
            break;
        case STRATEGYBYCELLS: {
//            printf("Using cells.....\n");
            return searchingByCells(tree, raster, queryBuilder, queries);
        }
        default:
//            printf("Using a queue....\n");
            qe = new k2raster_rtree_static::GetObjectsQueryStrategyQueue();
            break;
    }
    qe->setRaster(raster);
    qe->setAllCells(allCells);

    ulong totalMBRs = 0;
    ulong totalAllCells = 0;
//    printf("Executing %u queries.....\n", queryBuilder->getNumberOfQueries());
    start_timing();
    for (uint i = 0; i < queryBuilder->getNumberOfQueries(); i++) {
//        queries[i]->printQuery(i);
        qe->clear();
        qe->setRangeOfValues(queries[i]->valini, queries[i]->valend);
        tree->queryStrategy(*qe);
        totalMBRs += qe->getResult().size();
        totalAllCells += qe->getAllMBRCells();

//        printf("Query %u [%u, %u] -> %lu whites and %lu blacks || Result -> %lu\n", i, queries[i]->valini, queries[i]->valend, qe->getNumberOfWhites(), qe->getNumberOfBlacks(), qe->getResult().size());
//        break;
    }
    qe->clear();
    double totalTime = get_timing();
    printf("Time spent: %.6f seconds, %.6f seconds/query, MBRs -> %lu (%lu cell) %.6f nanoseconds/MBR, %.6f nanoseconds/cell\n",
           totalTime / 1000, (totalTime / 1000) / queryBuilder->getNumberOfQueries(),
           totalMBRs, totalAllCells,
           (totalTime * 1000000) / queryBuilder->getNumberOfQueries() / totalMBRs,
           (totalTime * 1000000) / queryBuilder->getNumberOfQueries() / totalAllCells);

    // Free memory
    // Free queries
    for (uint i = 0; i < queryBuilder->getNumberOfQueries(); i++) {
        delete (queries[i]);
    }
    free(queries);
    delete (queryBuilder);

    // Free raster
    delete (raster);

    // Free r-tree
    delete (qe);
    delete tree;
    delete file;
    delete diskfile;
}