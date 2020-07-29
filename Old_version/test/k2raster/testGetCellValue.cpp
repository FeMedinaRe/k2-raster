/*
 * Created by Fernando Silva on 13/01/17.
 *
 * Copyright (C) 2016-current-year, Fernando Silva, all rights reserved.
 *
 * Author's contact: Fernando Silva  <fernando.silva@udc.es>
 * Databases Lab, University of A Coruña. Campus de Elviña s/n. Spain
 *
 * Test: Execute a set of queries of GetCellValue over a k2-raster
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

int main(int argc, char **argv) {
    /*********************/
    /* Set params        */
    /*********************/
    uint nrows = 5000;          // Number of rows
    uint ncols = 5000;          // Number of columns
    uint k1 = 4;                // Division for the first levels
    uint k2 = 2;                // Division for the next levels
    uint levelK1 = 4;           // Number of levels to apply the division K1
    uint plainLevels = 2;       // Number of levels to save in plain form
    float minFreqVoc = 0.10;    // Minimum frequent to insert a submatrix in the vocabulary
    uint maxValue = 486;        // Max value of a cell
    uint numOfQueries=1000;     // Number of queries

    printf("Testing Encode k2-raster with param: \n");
    printf("k1 = %u || k2 = %u\n", k1, k2);
    printf("LevelK1 = %u || plainLevels = %u\n", levelK1, plainLevels);
    printf("MinFreqVoc = %.3f\n", minFreqVoc);
    printf("\n");

    /**********************/
    /* Create random data */
    /**********************/
    printf("Create a %ux%u sample data with maximum value of %u\n", nrows, ncols, maxValue);
    // Allocs memory and set a random value
    srand(time(nullptr));
    int **data = (int **) malloc(nrows * sizeof(int *));
    for (uint i = 0; i < nrows; i++) {
        data[i] = (int *) malloc(ncols * sizeof(int));
        for (uint j = 0; j < ncols; j++) {
            data[i][j] = rand() % (maxValue + 1);
        }
    }

    /*************************/
    /* Create random queries */
    /*************************/
    query_utils_static::Query **queries = new query_utils_static::Query *[numOfQueries];
    {
        uint xini, yini;
        for (uint q = 0; q < numOfQueries; q++) {
            xini = rand() % (nrows - 1);
            yini = rand() % (ncols - 1);
            queries[q] = new query_utils_static::Query(xini, 0, yini, 0, 0, 0);
        }
    }

    /*********************/
    /* Encodes data      */
    /*********************/
    k2raster_static::K2Raster *raster;
    std::string typeName;
    printf("Encoding values........\n");
    for (uint type = K2RASTER; type <= K2RASTER_VOC_ENTROPY; type++) {
        printf("\n");
        typeName = k2raster_static::K2Raster::getTypeName(type);

        // Create k2-raster
        raster = k2raster_static::K2Raster::createK2Raster(nrows, ncols, data, k1, k2, levelK1, plainLevels,
                                                           minFreqVoc, type);
        /**********************/
        /* Check encoded data */
        /**********************/
        if (!raster->checkEncode(data, nrows, ncols)) {
            printf("%s getCellValue (encoding phase): FAIL\n", typeName.data());
            delete (raster);
            continue;
        } else {
            printf("%s getCellValue (encoding phase): OK\n", typeName.data());
        }

        /**********************/
        /* Run GetCell        */
        /**********************/
        // Check cell by cell if all values are correct.
        int val;
        bool errorFound = false;
        for (uint x = 0; x < nrows; x++) {
            for (uint y = 0; y < ncols; y++) {
                val = raster->getCell(x, y);
                if (val != data[x][y]) {
                    printf("Value in position (%u, %u) -> %i, expected -> %i\n", x, y, val, data[x][y]);
                    errorFound = true;
                    break;
                }
            }
            if (errorFound) {
                break;
            }
        }

        if (errorFound) {
            printf("%s getCellValue (cell by cell phase): FAIL\n", typeName.data());
            delete (raster);
            continue;
        } else {
            printf("%s getCellValue (cell by cell phase): OK\n", typeName.data());
        }

        // Run random queries
        errorFound = false;
        for (uint q = 0; q < numOfQueries; q++) {
            uint x = queries[q]->xini;
            uint y = queries[q]->yini;
            val = raster->getCell(x, y);
            if (val != data[x][y]) {
                printf("Value in position (%u, %u) -> %i, expected -> %i\n", x, y, val, data[x][y]);
                queries[q]->printQuery(q + 1); printf("\n");
                errorFound = true;
                break;
            }
        }
        if (errorFound) {
            printf("%s getCellValue (random queries phase): FAIL\n", typeName.data());
            delete (raster);
            continue;
        } else {
            printf("%s getCellValue (random queries phase): OK\n", typeName.data());
        }

        printf("%s getCellValue: OK\n", typeName.data());

    } // END FOR type


    /*********************/
    /* Free data         */
    /*********************/
    for (uint i = 0; i < nrows; i++) {
        free(data[i]);
    }
    free(data);
}