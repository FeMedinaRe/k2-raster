/*
 * Created by Fernando Silva on 14/03/16.
 *
 * Copyright (C) 2016-current-year, Fernando Silva, all rights reserved.
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

#include <queue>
#include <k2raster/k2-raster-opt.h>

namespace k2raster_static {

    //**********************************************************************//
    //********************* CONSTRUCTORS ***********************************//
    //**********************************************************************//
    K2RasterOPT::K2RasterOPT() {};

    K2RasterOPT::K2RasterOPT(uint sizeX, uint sizeY, int **data, uint K1, uint K2, uint levelK1)
            : K2Raster(sizeX, sizeY, data, K1, K2, levelK1){
    }

    K2RasterOPT::~K2RasterOPT() {
    }


    //**********************************************************************//
    //********************* CELL FUNCTIONS *********************************//
    //**********************************************************************//


    ulong K2RasterOPT::getCellsByValue(Query query, uint *positions) const {

        int maxValue, minValue;
        int thisbaserow, thisbasecol, thismaxrow, thismaxcol;
        int newOffsetBitmap, offsetValues;
        bool uniform, isLeaf;
        ulong totalCells = 0;
        int xini, xend, yini, yend;
        uint kLevel = 0;
        uint newK;

        query.valend--;

        // Check first level
        if (this->initialMaxValue < query.valini ||
                this->initialMinValue > query.valend){
            return 0;
        }

        //Whole matrix is uniform
        if (!this->bitmapK2Tree->access(0)){
            // All cells are valid
            for (int x = query.xini; x <= query.xend; x++){
                for (int y= query.yend; y <= query.yend; y++){
                    positions[totalCells++] = x * this->realSizeY + y;
                }
            }
            return totalCells;
        }

        // Add Root node
        queue<QNodeCounters*> queue;
        queue.push(this->createRoot());
        QNodeCounters *node;

        while (!queue.empty()) {
            node = queue.front();
            queue.pop();

            isLeaf = node->level + 1 == this->nlevels - 1;
            kLevel = this->divKLevels[node->level + 1];
            newK = this->getK(node->level);
            int posI = 0, posJ = 0;
            newOffsetBitmap = 0;

            // Search his children
            int oldNewOffsetBitmap = -1;
            int oldOffsetValues = -1;
            while (newOffsetBitmap != -1) {
                newOffsetBitmap = this->nextChild(node, posI, posJ, thisbaserow, thismaxrow, thisbasecol, thismaxcol,
                                                  query, newK, kLevel);
                if (newOffsetBitmap == -1) {
                    break;
                }
//                printf("Searching form [%u, %u] to [%u, %u]\n", thisbaserow, thisbasecol, thismaxrow, thismaxcol);

                // Get max value
                if (oldNewOffsetBitmap == -1 || oldNewOffsetBitmap != newOffsetBitmap - 1) {
                    maxValue = this->listMaxValues[node->level]->access(newOffsetBitmap - node->preValuesMax);
                } else {
                    maxValue = this->listMaxValues[node->level]->next();
                }
                maxValue = node->maxValue - maxValue;
                oldNewOffsetBitmap = newOffsetBitmap;

                if (query.valini > maxValue){
                    // Values out of range
                    continue;
                }

                uniform = isLeaf || !this->bitmapK2Tree->access(newOffsetBitmap);
                if (!uniform){
                    // If is not uniform (or leaf) get his min value
                    offsetValues = this->bitmapK2Tree->rank1(newOffsetBitmap - 1) - node->preValuesMin;
                    if (oldOffsetValues == - 1 || oldOffsetValues != offsetValues - 1) {
                        minValue = this->listMinValues[node->level]->access(offsetValues);
                    } else {
                        minValue = this->listMinValues[node->level]->next();
                    }
                    oldOffsetValues = offsetValues;
                    minValue = node->minValue + minValue;
                } else {
                    minValue = maxValue;
                }

                if (query.valend < minValue){
                    // Values out of range
                    continue;
                }

                if (query.valini <= minValue && query.valend >= maxValue){
                    // All cell values within range
                    xini = max(thisbaserow, query.xini);
                    xend = min(thismaxrow, query.xend);
                    yini = max(thisbasecol, query.yini);
                    yend = min(thismaxcol, query.yend);
                    // Add all valid positions to final result
                    for (int x = xini; x <= xend; x++){
                        for (int y= yini; y <= yend; y++){
                            positions[totalCells++] = x * this->realSizeY + y;
                        }
                    }

                } else {
                    // Search for their children and they are added to queue
                    // If it is uniform, the node has not children
                    if (!uniform) {
                        queue.push(this->createChild(node, offsetValues, thisbaserow, thisbasecol, minValue, maxValue));
                    }
                }
            }
            free(node);
        }
        return  totalCells;
    }

    ulong K2RasterOPT::getValuesWindow(Query query, int *values) const {

        int maxValue, minValue;
        int thisbaserow, thisbasecol, thismaxrow, thismaxcol;
        uint xini, xend, yini, yend, cellPos;
        int newOffsetBitmap = 0, offsetValues,oldOffsetBitmap;
        bool uniform, isLeaf;
        ulong totalValue = 0;
        uint kLevel = 0;
        uint newK;

//        query.valini = 0;
//        query.valend = this->getMaxValue();
//        query.valend--;

        // Invalid values
//        if (this->initialMinValue > query.valend ||
//            this->initialMaxValue < query.valini) {
//            return 0;
//        }

        //Whole matrix is uniform
        if (!this->bitmapK2Tree->access(0)){
            for (uint x = 0; x < (query.xend - query.xini) + 1; x++) {
                for (uint y = 0; y < (query.yend - query.yini) + 1; y++) {
                    values[x * (query.yend - query.yini) + y] = this->initialMaxValue;
                    totalValue++;
                }
            }
            return totalValue;
        }

        // Create root node
        queue<QNodeCounters*> queue;
        QNodeCounters *node;
        queue.push(this->createRoot());

        while (!queue.empty()) {
            node = queue.front();
            queue.pop();

            isLeaf = node->level + 1 == this->nlevels - 1;
            kLevel = this->divKLevels[node->level + 1];
            newK = this->getK(node->level);
            int posI = 0, posJ = 0;
            newOffsetBitmap = 0;
            oldOffsetBitmap = -1;

            // Search his children
            while (newOffsetBitmap != -1) {
                newOffsetBitmap = this->nextChild(node, posI, posJ, thisbaserow, thismaxrow, thisbasecol, thismaxcol,
                                                  query, newK, kLevel);
                if (newOffsetBitmap == -1) {
                    break;
                }

                uniform = isLeaf || !this->bitmapK2Tree->access(newOffsetBitmap);

                // Get max value
                if (oldOffsetBitmap == -1 || oldOffsetBitmap != newOffsetBitmap - 1) {
                    maxValue = this->listMaxValues[node->level]->access(newOffsetBitmap - node->preValuesMax);
                } else {
                    maxValue = this->listMaxValues[node->level]->next();
                }
                maxValue = node->maxValue - maxValue;
                oldOffsetBitmap = newOffsetBitmap;

//                if (maxValue < query.valini) {
//                    continue;
//                }

                if (uniform){
                    // All cell values within range
                    xini = max(thisbaserow, query.xini) - query.xini;
                    xend = min(thismaxrow, query.xend) - query.xini;
                    yini = max(thisbasecol, query.yini) - query.yini;
                    yend = min(thismaxcol, query.yend) - query.yini;

                    // Add all valid positions to final result
                    for (uint x = xini; x <= xend; x++){
                        cellPos = (x * ((query.yend - query.yini) + 1));
                        for (uint y= yini; y <= yend; y++){
                            values[cellPos + y] = maxValue ;
                            totalValue++;
                        }
                    }
                    continue;
                }

                offsetValues = this->bitmapK2Tree->rank1(newOffsetBitmap - 1) - node->preValuesMin;
//                minValue = node->minValue + this->listMinValues[node->level]->access(offsetValues);
//                if (minValue > query.valend) {
//                    continue;
//                }

                // Search his children
                queue.push(this->createChild(node, offsetValues, thisbaserow, thisbasecol, minValue, maxValue));
            }

            free(node);
        }

        return  totalValue;
    }

    bool K2RasterOPT::checkValues(Query query, bool allCells) const {
        int maxValue, minValue;
        int thisbaserow, thisbasecol, thismaxrow, thismaxcol;
        int newOffsetBitmap, offsetValues, oldOffsetBitmap, oldOffsetValues;
        bool uniform, isLeaf;
        uint kLevel = 0;
        uint newK;

        query.valend--;

        // Any valid value in matrix
        if (this->initialMaxValue < query.valini ||
                this->initialMinValue > query.valend){
            return 0;
        }

        //Whole matrix uniform
        if (!this->bitmapK2Tree->access(0)){
            return true;
        }

        // Size of window is equal to size of matrix
        if (query.xini == 0 && query.yini == 0
                    && query.xend == (this->realSizeX - 1)
                    && query.yend == (this->realSizeY - 1)){
            if (!allCells) {
                return true;
            } else {
                return (query.valini <= this->initialMinValue && query.valend >= this->initialMaxValue);
            }
        }

        queue<QNodeCounters*> queue;
        queue.push(this->createRoot());
        QNodeCounters *node;

        while (!queue.empty()) {
            node = queue.front();
            queue.pop();

            isLeaf = node->level + 1 == this->nlevels - 1;
            kLevel = this->divKLevels[node->level + 1];
            newK = this->getK(node->level);

            // Search his children
            int posI = 0, posJ = 0;
            newOffsetBitmap = 0;
            oldOffsetBitmap = -1;
            oldOffsetValues = - 1;

            // Search his children
            while (newOffsetBitmap != -1) {
                newOffsetBitmap = this->nextChild(node, posI, posJ, thisbaserow, thismaxrow, thisbasecol, thismaxcol,
                                                  query, newK, kLevel);
                if (newOffsetBitmap == -1) {
                    break;
                }


                // Get max value
                if (oldOffsetBitmap == -1 || oldOffsetBitmap != newOffsetBitmap - 1) {
                    maxValue = this->listMaxValues[node->level]->access(newOffsetBitmap - node->preValuesMax);
                } else {
                    maxValue = this->listMaxValues[node->level]->next();
                }
                maxValue = node->maxValue - maxValue;
                oldOffsetBitmap = newOffsetBitmap;

                // Check maxValue
                if (query.valini > maxValue){
                    // Values out of range
                    if (allCells) {
                        return false;
                    }
                    continue;
                }

                // Check it is uniform
                uniform = isLeaf || !this->bitmapK2Tree->access(newOffsetBitmap);

                // Get max value
                if (uniform){
                    if (query.valini <= maxValue && query.valend >= maxValue){
                       if (!allCells) {
                           return  true;
                       }
                    } else {
                        if (allCells) {
                            return false;
                        }
                    }
                    continue;
                }

                // If is not uniform (or leaf) get his min value
                offsetValues = this->bitmapK2Tree->rank1(newOffsetBitmap - 1) - node->preValuesMin;
                if (oldOffsetValues == - 1 || oldOffsetValues != offsetValues - 1) {
                    minValue = this->listMinValues[node->level]->access(offsetValues);
                } else {
                    minValue = this->listMinValues[node->level]->next();
                }
                oldOffsetValues = offsetValues;
                minValue = node->minValue + minValue;
//                offsetValues = this->bitmapK2Tree->rank1(newOffsetBitmap - 1) - node->preValuesMin;
//                minValue = node->minValue + this->listMinValues[node->level]->access(offsetValues);

                // Check minValue
                if (query.valend < minValue){
                    // Values out of range
                    if (allCells) {
                        return false;
                    }
                    continue;
                }

                if (query.valini <= minValue && query.valend >= maxValue){
                    // All values within range
                    if (!allCells) {
                        return true;
                    } else {
                        continue;
                    }
                } else {
                    if (thisbaserow >= query.xini && thismaxrow <= query.xend &&
                                     thisbasecol >= query.yini && thismaxcol <= query.yend) {
                        // Submatrix within of windows and not all values within range

                        if (allCells) {
                            return false;
                        } else {
                            if (minValue >= query.valini || maxValue <= query.valend) {
                                //minValue or maxValue are in range. One or more cells have a right value
                                return true;
                            }
                        }
                    }
                }

                // Search his children
                queue.push(this->createChild(node, offsetValues, thisbaserow, thisbasecol, minValue, maxValue));
            }
            free(node);
        }
        /**
         * If allcells is true -> not found any cell with a not valid value
         * if allcell is false -> not found any cell with a valid value
         */
        return allCells;
    }


    //**********************************************************************//
    //********************* FILE FUNCTIONS *********************************//
    //**********************************************************************//
    void K2RasterOPT::save(std::ofstream &of) const {
        K2Raster::save(of, K2RASTER_OPT);
    }

    /*
     * Load Ktreap from file
     */
    K2RasterOPT *K2RasterOPT::load(std::ifstream &in) {
        K2RasterOPT *ktreap = nullptr;

        try {
            ktreap = new K2RasterOPT();
            K2Raster::loadBase(in, ktreap);

        } catch (...) {
            return nullptr;
        }

        return ktreap;
    }
}
