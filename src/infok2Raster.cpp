/*
 * Created by Fernando Silva on 15/02/16.
 *
 * Copyright (C) 2016-current-year, Fernando Silva, all rights reserved.
 *
 * Author's contact: Fernando Silva  <fernando.silva@udc.es>
 * Databases Lab, University of A Coruña. Campus de Elviña s/n. Spain
 *
 * Program to execute k2-raster over a dataset (sequence of integers of 32bits)
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
#include <k2raster/k2-raster.h>
#include <k2raster/plain/k2-raster-plain.h>
#include <k2raster/plain/k2-raster-plain-DAC.h>
#include <k2raster/plain/k2-raster-plain-VByte.h>
#include <k2raster/compressLeaves/k2-raster-CompressLeaves.h>
#include <k2raster/compressLeaves/k2-raster-CompressLeavesH.h>
#include <k2raster/k2-raster-opt.h>
#include <k2raster/k2-raster-entropy.h>

using namespace std;

int main(int  argc, char ** argv) {

    if (argc < 6) {
        printf("Usage: %s <dataFile> <rows> <cols> <k2rasterFile> <check> [type] [K1 K2 levelK1 [levelK2]] [minFreqVoc] \n",
               argv[0]);
        exit(-1);
    }

    /*********************/
    /* Reads params      */
    /*********************/
    uint nrows =  atoi(argv[2]);                        // Number of rows
    uint ncols =  atoi(argv[3]);                        // Number of columns
    bool check = atoi(argv[5]);                         // If run a full check at the end of the process
    uint type = (argc < 7 ? 0 : atoi(argv[6]));         // Type of encode
    uint K1 = (argc < 10 ? 2 : atoi(argv[7]));          // Division for the first levels
    uint K2 = (argc < 10 ? 2 : atoi(argv[8]));          // Division for the next levels
    uint levelK1 = (argc < 10 ? 2 : atoi(argv[9]));     // Number of levels to apply the division K1
    uint levelK2 = (argc < 11 ? 2 : atoi(argv[10]));    // Number of levels to apply the division K2
    float minFreqVoc = (argc < 12 ? 0 : stof(argv[11]));// Minimum frequent to insert a submatrix in the vocabulary


    /*********************/
    /* Reads input data   */
    /*********************/


//    for (uint i = 0; i < nrows; i++) {
//        data[i] = (int *) malloc(ncols * sizeof(int));
//        for (uint j = 0; j < ncols; j++) {
//            data[i][j] = -1;
//        }
//    }

    // Reads file
    ifstream inputFile(argv[1]);
    if (!inputFile.good()) {
        printf("File %s unable to open\n", argv[1]);
        exit(-1);
    }

    // Allocs memory
    int **data = (int **) malloc(nrows * sizeof(int *));
    for (uint i=0; i<nrows; i++){
        data[i] = (int *) malloc(ncols * sizeof(int));
        for (uint j=0; j<ncols; j++){
            data[i][j] = cds_utils::loadValue<int>(inputFile);
        }
    }
    inputFile.close();

    /*********************/
    /* Encodes data      */
    /*********************/
    k2raster_static::K2Raster *raster;
    std::string typeName;
    cds_utils::start_timing();
    switch (type){
        case K2RASTER:
            // Normal k2-raster. It divides the first "levelK1" levels by "K1" and the rest by "K2"
            raster = new k2raster_static::K2Raster(nrows, ncols, data, K1, K2, levelK1);
            typeName = "k2-raster\t\t";
            break;
        case K2RASTER_PLAIN:
            // Plain k2-raster. It divides the first "levelK1" levels by "K1", the next "levelK2" levels by "K2"
            // and it stores the rest of the values in plain form (array of integers)
            raster = new k2raster_static::K2RasterPlain(nrows, ncols, data, K1, K2, levelK1, levelK2);
            typeName = "k2-raster-plain";
            break;
        case K2RASTER_PLAIN_DAC:
            // Plain DAC k2-raster. It divides the first "levelK1" levels by "K1", the next "levelK2" levels by "K2"
            // and it stores the rest of the values encoded with DAC
            raster = new k2raster_static::K2RasterPlainDAC(nrows, ncols, data, K1, K2, levelK1, levelK2);
            typeName = "k2-raster-DAC";
            break;
        case K2RASTER_PLAIN_VBYTE:
            // Plain VByte k2-raster. It divides the first "levelK1" levels by "K1", the next "levelK2" levels by "K2"
            // and it stores the rest of the values encoded with VByte
            raster = new k2raster_static::K2RasterPlainVByte(nrows, ncols, data, K1, K2, levelK1, levelK2);
            typeName = "k2-raster-VByte";
            break;
        case K2RASTER_VOC:
            // Vocabulary k2-raster. It divides the first "levelK1" levels by "K1", the next "levelK2" levels by "K2"
            // and it assigns a word of the vocabulary for each different final submatrix
            raster = new k2raster_static::K2RasterCompressLeaves(nrows, ncols, data, K1, K2, levelK1, levelK2);
            typeName = "k2-raster-voc";
            break;
        case K2RASTER_VOC_HYBRID:
            // Plain k2-raster. It divides the first "levelK1" levels by "K1", the next "levelK2" levels by "K2"
            // and if the sumbatrix is frequent, it is stored as a word of vocabulary, in other case, it are stored in plain form
            raster = new k2raster_static::K2RasterCompressLeavesH(nrows, ncols, data, K1, K2, levelK1, levelK2, minFreqVoc);
            typeName = "k2-raster-hybrid";
            break;
        case K2RASTER_OPT:
            // Optimized k2-raster. It divides the first "levelK1" levels by "K1" and the rest by "K2"
            // (with optimized DAC)
            raster = new k2raster_static::K2RasterOPT(nrows, ncols, data, K1, K2, levelK1);
            typeName = "k2-raster-opt";
            break;
        case K2RASTER_VOC_ENTROPY:
            // Entropy k2-raster. It divides the first "levelK1" levels by "K1", the next "levelK2" levels by "K2"
            // Implemented a entropy function to calculate, for each submatrix, whether it is better to store it in plain form
            // or as a word of vocabulary.
            raster = new k2raster_static::K2RasterEntropy(nrows, ncols, data, K1, K2, levelK1, levelK2);
            typeName = "k2-raster-entropy";
            break;
        default:
            raster = new k2raster_static::K2Raster(nrows, ncols, data, K1, K2, levelK1);
            typeName = "k2-raster\t\t";
            break;
    } // END SWITCH

    double buildTime = cds_utils::get_timing();
    double compressionRatio = (raster->getTotalSize() * 100.) / (nrows * ncols * sizeof(int));
    printf("Create %s \t(%u, %u) in %.3f seconds with size %lu bytes (%.2f\%)\n",  typeName.data(), nrows, ncols, buildTime/1000, raster->getTotalSize(), compressionRatio);

    /**********************/
    /* Check encoded data */
    /**********************/
    if (check){
        // Check cell by cell if all values are correct.
        printf("Check values........\n");
        if (!raster->checkEncode(data, nrows, ncols)){
            printf("Error found in encoded treap!!!");
            return -1;
        }
    }

    /*********************/
    /* Free data         */
    /*********************/
    for (uint i = 0; i < nrows; i++) {
        free(data[i]);
    }
    free(data);

    /****************************/
    /* Save k2-raster in a file */
    /****************************/
    cds_utils::start_timing();
    raster->saveSplit(argv[4]);
    //printf("Save RegionTreap in path %s in %.3f seconds\n",  argv[4], get_timing()/1000);
    delete(raster);
}