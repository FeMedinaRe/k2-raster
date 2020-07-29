/*  
 * Created by Fernando Silva on 30/10/17.
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
// Libraries
#include <includes/utils/timing.h>
#include <rtree/Statistic/QueryStrategy/GetObjectsQueryMaxTopK.h>

int main(int argc, char **argv) {

    if (argc != 5) {
        printf("Usage: %s <k2rasterFile> <rtreeFile> <k> <withPositions>\n", argv[0]);
        exit(-1);
    }

    // ********************** //
    // Get Params             //
    // ********************** //
    uint k = atoi(argv[3]);
    if (k < 1) {
        printf("k must be greater than 0, current %i\n", k);
        exit(-1);
    }

    bool withPositions = atoi(argv[4]);

    // ********************** //
    // Load k2-raster         //
    // ********************** //
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

    // ********************** //
    // Load R-Tree            //
    // ********************** //
    string baseName = argv[2];
//    printf("Loading r-tree from %s\n", baseName.data());
    IStorageManager *diskfile = StorageManager::loadDiskStorageManager(baseName);
    StorageManager::IBuffer *file = StorageManager::createNewRandomEvictionsBuffer(*diskfile, 10, false);
    ISpatialIndex *tree = RTree::loadRTree(*file, 1);

    // ********************** //
    // Create Strategy        //
    // ********************** //
    k2raster_rtree_static::GetObjectsQueryMaxTopK *qe = new k2raster_rtree_static::GetObjectsQueryMaxTopK();
    qe->setRaster(raster);
    qe->setK(k);
    qe->setWithPositions(withPositions);
    qe->initialize();

    // ********************** //
    // Run algorithm          //
    // ********************** //
    start_timing();
    tree->queryStrategy(*qe);
    double totalTime = get_timing();
    printf("Time spent: %.6f seconds\n", totalTime / 1000);

    printf("Results (%u of %u): \n", qe->getNumResults(), k);
    for (uint mbr = 0; mbr < qe->getNumResults(); mbr++ ) {
        printf("MBR %lu: max value %i\n", qe->getMBRIds()[mbr], qe->getValues()[mbr]);
    }

    // ********************** //
    // Free memory            //
    // ********************** //
    // Strategy
    qe->clear();

    // Free raster
    delete (raster);

    // Free r-tree
    delete (qe);
    delete tree;
    delete file;
    delete diskfile;
}