/*  
 * Created by Fernando Silva on 27/09/17.
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

// Own header
#include <k2raster/k2-raster.h>
#include <algebra/algebra.h>
#include <includes/utils/timing.h>

using namespace std;

int main(int  argc, char ** argv) {

    if (argc != 9) {
        printf("Usage: %s <k2rasterFile1> <k2rasterFile2> <outputK2rasterFile> <operation> <k1> <k2> <levelK1> <check>\n",
               argv[0]);
        exit(-1);
    }

    uint op = atoi(argv[4]);
    uint k1 = atoi(argv[5]);
    uint k2 = atoi(argv[6]);
    uint levelK1 = atoi(argv[7]);
    bool check = atoi(argv[8]);

    // ------------------------------------------------------------------- //
    // Load k2-raster 1
    // ------------------------------------------------------------------- //
    ifstream inFile1(argv[1]);
    if (!inFile1.good()) {
        printf("File %s unable to open\n", argv[1]);
        exit(-1);
    }
    k2raster_static::K2Raster *raster1 = k2raster_static::K2Raster::load(inFile1);
    inFile1.close();
    if (raster1 == nullptr) {
        printf("Error loading k2-raster from file %s\n", argv[1]);
    }

    // ------------------------------------------------------------------- //
    // Load k2-raster 2
    // ------------------------------------------------------------------- //
    ifstream inFile2(argv[2]);
    if (!inFile2.good()) {
        printf("File %s unable to open\n", argv[2]);
        exit(-1);
    }
    k2raster_static::K2Raster *raster2 = k2raster_static::K2Raster::load(inFile2);
    inFile2.close();
    if (raster2 == nullptr) {
        printf("Error loading k2-raster from file %s\n", argv[2]);
    }

    // ------------------------------------------------------------------- //
    // Run operation over both raster
    // ------------------------------------------------------------------- //
    cds_utils::start_timing();
    k2raster_static::K2Raster *rasterResult = k2raster_algebra_static::algebra(raster1, raster2,
                                                                               (OperationRaster)op,
                                                                               k1, k2, levelK1, !check);
    if (rasterResult == nullptr) {
        printf("Error!!\n");
        exit(-1);
    }
    double buildTime = cds_utils::get_timing();
    double compressionRatio = (rasterResult->getTotalSize() * 100.) / (rasterResult->getDimSize(1) * rasterResult->getDimSize(2) * sizeof(int));
    printf("Create raster \t(%u, %u) in %.3f seconds with size %lu bytes (%.2f\%)\n", rasterResult->getDimSize(1), rasterResult->getDimSize(2),
           buildTime/1000, rasterResult->getTotalSize(), compressionRatio);
    if (check) {
        double compressionRatio1 = (raster1->getTotalSize() * 100.) / (raster1->getDimSize(1) * raster1->getDimSize(2) * sizeof(int));
        double compressionRatio2 = (raster2->getTotalSize() * 100.) / (raster2->getDimSize(1) * raster2->getDimSize(2) * sizeof(int));
        printf("Originals raster: Raster1 -> %lu bytes (%.2f\%) and Raster2 -> %lu bytes (%.2f\%)\n",
               raster1->getTotalSize(), compressionRatio1, raster2->getTotalSize(), compressionRatio2);
    }

    // ------------------------------------------------------------------- //
    // Check result
    // ------------------------------------------------------------------- //
    if (check) {
        printf("Checking k2-raster........\n");
        uint sizeX = rasterResult->getDimSize(1);
        uint sizeY = rasterResult->getDimSize(2);

        for (uint x = 0; x < sizeX; x++) {
            for (uint y = 0; y < sizeY; y++) {
                int result1 = 0 ;
                switch (op) {
                    case OperationRaster::OPERATION_SUM:
                        result1 = raster1->getCell(x, y) + raster2->getCell(x, y);
                        break;
                    case OperationRaster::OPERATION_SUBT: {
                        result1 = raster1->getCell(x, y) - raster2->getCell(x, y);
                        break;
                    }
                    case OperationRaster::OPERATION_MULT:
                        result1 = raster1->getCell(x, y) * raster2->getCell(x, y);
                        break;
                    default:
                        printf("No valid operation!!\n");
                        exit(-1);
                }
                if (result1 != rasterResult->getCell(x, y)) {
                    printf("Found error at position (%u, %u), expected %i and get %i\n", x,y, result1, rasterResult->getCell(x, y));
                    exit(-1);
                }
            }
        }
        printf("ALL OK!!!!\n");
        delete(raster1);
        delete(raster2);
    }

    // ------------------------------------------------------------------- //
    // Save result
    // ------------------------------------------------------------------- //
    cds_utils::start_timing();
    ofstream outFile(argv[3]);
    rasterResult->save(outFile);
    outFile.close();        // Close file
    //printf("Save k2-tree in path %s in %.3f seconds\n",  argv[4], get_timing()/1000);
    delete(rasterResult);
}