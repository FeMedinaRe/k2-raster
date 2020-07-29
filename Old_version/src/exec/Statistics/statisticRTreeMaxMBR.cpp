/*  
 * Created by Fernando Silva on 25/10/17.
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
#include <rtree/Statistic/QueryStrategy/GetObjectsQueryMaxMBR.h>
// Libraries
#include <includes/utils/timing.h>

int main(int argc, char **argv) {

    if (argc != 3) {
        printf("Usage: %s <k2rasterFile> <rtreeFile>\n", argv[0]);
        exit(-1);
    }

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
    k2raster_rtree_static::GetObjectsQueryMaxMBR *qe = new k2raster_rtree_static::GetObjectsQueryMaxMBR();
    qe->setRaster(raster);

    // ********************** //
    // Run algorithm          //
    // ********************** //
    start_timing();
    qe->clear();
    tree->queryStrategy(*qe);
    double totalTime = get_timing();
    printf("Time spent: %.6f seconds, MBR %lu contains the max value %i\n", totalTime / 1000, qe->getMBRId(), qe->getMaxValue());

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