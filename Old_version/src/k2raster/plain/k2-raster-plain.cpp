/*
 * Created by Fernando Silva on 24/05/16.
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

#include <sys/types.h>
#include <queue>
#include <k2raster/plain/k2-raster-plain.h>

namespace k2raster_static {

    //**********************************************************************//
    //********************* CONSTRUCTORS ***********************************//
    //**********************************************************************//

    K2RasterPlain::K2RasterPlain() {
        this->plainValues = nullptr;
        this->preMinValues = nullptr;
        this->preMaxValues = nullptr;
    };

    K2RasterPlain::K2RasterPlain(uint sizeX, uint sizeY, int **data, uint k1, uint k2, uint levelK1, uint plainLevels) {
        // Set K values and store real sizes
        this->K1 = k1;
        this->K2 = k2;
        this->levelK1 = levelK1;
        this->realSizeX = sizeX;
        this->realSizeY = sizeY;

        // Initialize levels  with new values of K
        this->initLevels(plainLevels);

        // Execute algorithm of construction
        uint l = this->build(sizeX, sizeY, data);

        // Insert plain values
        if (this->isLevelInPlain(l + 1)) {
            this->plainValues = (uint *)malloc(this->tmpPlainValues->size() * sizeof(uint));
            for (uint i = 0; i < tmpPlainValues->size(); i++) {
                this->plainValues[i] = tmpPlainValues->at(i);
            }
            this->numOfPlainValues = tmpPlainValues->size();
            delete(tmpPlainValues);
        }

    }

    K2RasterPlain::~K2RasterPlain() {
        if (this->plainValues != nullptr) {
            free(this->plainValues);
        }
    }

    //**********************************************************************//
    //********************* CELL FUNCTIONS *********************************//
    //**********************************************************************//

    int K2RasterPlain::getPlainCell(uint father, uint row, uint col, uint divK, uint lastValueMax) const {
        uint arrayPosition = father * divK * divK;
        row = row % divK;
        col = col % divK;
        arrayPosition +=  row * divK + col;
        return lastValueMax - this->plainValues[arrayPosition];
    }

    ulong K2RasterPlain::getPlainWindow(uint xini, uint yini, uint xend, uint yend,
                                        uint offsetValues, uint divK, uint maxValue,
                                        uint *positions, ulong totalCells, Query query) const {
        return this->getPlainWindow(xini, yini, xend, yend, query.valini, query.valend,
                                    offsetValues, divK, maxValue, positions, totalCells, false);
    }

    ulong K2RasterPlain::getPlainWindow(uint xini, uint yini, uint xend, uint yend,
                                        uint valini, uint valend,
                                        uint offsetValues, uint divK, uint maxValue,
                                        uint *positions, ulong totalCells, bool allCells) const {
        // Add all valid positions to final result
        uint value, arrayPosition;
        for (uint x = xini; x <= xend; x++){
            arrayPosition = (offsetValues * divK * divK) + (x % divK) * divK + (yini % divK);
            for (uint y= yini; y <= yend; y++){
                value = maxValue - this->plainValues[arrayPosition];
                if (value >= valini && value <= valend) {
                    positions[totalCells++] = x * this->realSizeY + y;
                } else {
                    if (allCells) {
                        return 0;
                    }
                }
                arrayPosition++;
            }
        }
        return totalCells;
    }

    bool K2RasterPlain::checkPlainCells(uint xini, uint yini, uint xend, uint yend, uint offsetValues, uint divK,
                                        uint maxValue, Query query, bool allCells) const {

        // Add all valid positions to final result
        uint value, arrayPosition;
        for (uint x = xini; x <= xend; x++){
            arrayPosition = (offsetValues * divK * divK) + (x % divK) * divK + (yini % divK);
            for (uint y= yini; y <= yend; y++){
                value = maxValue - this->plainValues[arrayPosition];
                if (value >= query.valini && value <= query.valend) {
                    if (!allCells) {
                        return true;
                    }
                } else {
                    if (allCells) {
                        return false;
                    }
                }
                arrayPosition++;
            }
        }
        return allCells;
    }

    ulong K2RasterPlain::getPlainCellValues(uint xini, uint yini, uint xend, uint yend, uint offsetValues, uint divK,
                                            uint maxValue, int *values, ulong totalCells, Query query) const {
        return getPlainCellValues(xini, yini, xend, yend, offsetValues, divK, maxValue, values,
                                  totalCells, query.xini, query.yini, query.yend);
    }

    ulong K2RasterPlain::getPlainCellValues(uint xini, uint yini, uint xend, uint yend, uint offsetValues, uint divK,
                                            uint maxValue, int *values, ulong totalCells,
                                            uint initXini, uint initYini, uint initYend) const {
        uint value, arrayPosition;
        uint cellPos;
        for (uint x = xini; x <= xend; x++){
            arrayPosition = (offsetValues * divK * divK) + (x % divK) * divK + (yini % divK);
            cellPos = ((x - initXini) * ((initYend - initYini) + 1));
            for (uint y= yini; y <= yend; y++){
                value = maxValue - this->plainValues[arrayPosition];
                values[cellPos + (y - initYini) ] = value;
                totalCells++;
                arrayPosition++;
            }
        }
        return totalCells;
    }

    //*********************//
    //****** GETCELL ******//
    //*********************//

    int K2RasterPlain::getCell(uint row, uint col) const {
        uint lastValueMax = 0;
        uint thisValueMax;
        uint posBitmap = 0, posBitmapPrev = this->getK(0) * this->getK(0);
        uint uniform, divK, offset;
        uint totalMinValues = 0, totalMinValuesPrev = 1;

        // Check firts level
        if (!this->bitmapK2Tree->access(0)){
            // All cell of the matrix are equal
            return this->initialMaxValue;
        }

        lastValueMax = this->initialMaxValue;   // Parent maximum value
        posBitmap = 1;                          // Position in k2-tree
        divK = this->divKLevels[0];             // Size of partition of submatrix
        row = row % divK;                       // Relative row in submatrix
        col = col % divK;                       // Relative column  in submatrix
        uint prevNumber = 1;                    // Number total of maximum numbers in previous levels

        // Search from level 1 to level l - 1
        uint levels = this->levelK1 + this->levelK2;
        levels = levels == this->nlevels ? levels - 1 : levels;
        for (uint l = 1; l < levels; l++) {
            divK = this->divKLevels[l];
            offset = (row / divK) * this->getK(l - 1) + col / divK;     // Relative position
            posBitmap += offset;                                        // new position in k2-tree

            // Cheks if it is uniform and get this maxValue
            uniform =! this->bitmapK2Tree->access(posBitmap);
            thisValueMax = lastValueMax - this->listMaxValues[l - 1]->access(posBitmap - prevNumber);

            if (!uniform) {
                totalMinValues = this->bitmapK2Tree->rank1(posBitmap - 1); // Number of 1s until actual position
                // Save values and search in the next level (his children)
                lastValueMax = thisValueMax;
                totalMinValues -= totalMinValuesPrev; // Number of 1s in this level until actual position
                posBitmap = posBitmapPrev  + (totalMinValues * this->getK(l) *  this->getK(l) + 1); // Children position
                row = row % divK;
                col = col % divK;
            } else {
                // It is uniform and all cells have equal value
                return thisValueMax;
            }

            /**
             * Update counters previous numbers.
             * Those counters are necessary to speed up the query and because
             * the value of "k" can be different in each level, else it is more difficult calculate the position
             * of the next submatrix
             */
            prevNumber += this->listMaxValues[l - 1]->getRealLength();
            posBitmapPrev  += this->listMinValues[l -1]->getRealLength() * this->getK(l) * this->getK(l); // Position of next level in k2-tree
            totalMinValuesPrev += this->listMinValues[l - 1]->getRealLength();  // number of minValues == number of 1s in k2-tree
        }

        if (!this->isLevelInPlain(levels)){
            // No plain values, search last level
            divK = this->divKLevels[this->nlevels - 1];
            offset = (row / divK) * this->getK(this->nlevels - 1) + col / divK;
            posBitmap += offset;
            posBitmap -= prevNumber;

            thisValueMax = this->listMaxValues[this->nlevels - 2]->access(posBitmap);
            thisValueMax = lastValueMax - thisValueMax;
            return thisValueMax;
        } else {
            // Values of last levels are in plain form
            return this->getPlainCell(totalMinValues, row, col, divK, lastValueMax);
        }
    }

    //*****************************//
    //****** GETCELLBYVALUE ******//
    //****************************//
    ulong K2RasterPlain::getCellsByValue(Query query, uint *positions) const {

        int maxValue, minValue;
        int thisbaserow, thisbasecol, thismaxrow, thismaxcol;
        int newOffsetBitmap, offsetValues, oldNewOffsetBitmap, oldOffsetValues;
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
            oldNewOffsetBitmap = -1;
            oldOffsetValues = -1;

            // Search his children
            while (newOffsetBitmap != -1) {
                newOffsetBitmap = this->nextChild(node, posI, posJ, thisbaserow, thismaxrow, thisbasecol, thismaxcol,
                                                  query, newK, kLevel);
                if (newOffsetBitmap == -1) {
                    break;
                }

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

                if (!uniform && this->isLevelInPlain(node->level + 2)) {
                    // Search in plain values
                    xini = max(thisbaserow, query.xini);
                    xend = min(thismaxrow, query.xend);
                    yini = max(thisbasecol, query.yini);
                    yend = min(thismaxcol, query.yend);
//                    uint prevCells = totalCells;
                    totalCells = this->getPlainWindow(xini, yini, xend, yend,
                                                      offsetValues, kLevel, maxValue,
                                                      positions, totalCells, query);
//                    printf("Added %u cells between [%u, %u] - [%u, %u]\n", totalCells - prevCells,
//                           xini, xend, yini, yend);
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
                    if (!uniform) {
                        queue.push(this->createChild(node, offsetValues, thisbaserow, thisbasecol, minValue, maxValue));
                    }
                }
            }
            free(node);
        }
        return  totalCells;
    }

    ulong K2RasterPlain::getCellsByValue(uint xini, uint xend, uint yini, uint yend, uint valini, uint valend,
                                         MinBasicSubMatrix *minSubMatrix, uint *positions, bool allCells) const {
        if (minSubMatrix->minValue >= valini && minSubMatrix->maxValue <= valend) {
            // All cells are equal
            ulong totalCells = 0;
            for (ulong x = xini; x <= xend; x++) {
                for (ulong y = yini; y <= yend; y++) {
                    positions[totalCells++] = x * this->realSizeY + y;
                }
            }
            return totalCells;
        }

        uint divKLevel = this->divKLevels[minSubMatrix->level];

        if (this->isLevelInPlain(minSubMatrix->level + 1)) {
            // Submatrix represented by plain values
            uint countOnes = this->bitmapK2Tree->rank1(minSubMatrix->bitmapPosition - 1);
            return this->getPlainWindow(xini, yini, xend, yend, valini, valend,
                                        countOnes - this->preMinValues[minSubMatrix->level - 1],
                                        divKLevel, minSubMatrix->maxValue, positions, 0, allCells);

        } else {
            // Calculate initial position
            uint baseRow = (xini / divKLevel) * divKLevel;
            uint baseCol = (yini / divKLevel) * divKLevel;

            // Check if the MBR overlap completely with the submatrix
            if (allCells && xini == baseRow && xend == baseRow + divKLevel
                && yini == baseCol && yend == baseCol + divKLevel) {
                // Some value is without of the range of values
                return 0;
            }

            return this->getCellsByValue(xini, xend, yini, yend, minSubMatrix->minValue, minSubMatrix->maxValue,
                                         minSubMatrix->level, minSubMatrix->bitmapChildren, valini, valend,
                                         baseRow, baseCol, 0, positions, allCells);
        }
    }

    ulong
    K2RasterPlain::getCellsByValue(uint xini, uint xend, uint yini, uint yend, uint fatherMinValue, uint fatherMaxValue,
                                   uint l, uint bitmapChildren, uint valini, uint valend,
                                   uint baseRow, uint baseCol,
                                   ulong totalCells, uint *positions, bool allCells) const {

        uint divKLevel = this->divKLevels[l + 1];      // Size of its children
        uint k = this->getK(l);
        uint level = l + 1;
        bool uniform, isLeaf = level == this->nlevels - 1;

        uint bitmapPosition;
        int minValue, maxValue;
        uint countOnes;
        uint thisBaseRow, thisBaseCol;

        // Set limits
        uint limit1X = (xini - baseRow) / divKLevel;
        uint limit2X = (xend - baseRow) / divKLevel;
        uint limit1Y = (yini - baseCol) / divKLevel;
        uint limit2Y = (yend - baseCol) / divKLevel;

        bool firstMin = false;
        for (uint x = limit1X; x <= limit2X; x++) {
            thisBaseRow = baseRow + x * divKLevel;
            firstMin = firstMin && limit1Y == 0 && limit2Y == (k - 1);
            for (uint y = limit1Y; y <= limit2Y; y++) {
                // Get position at k2-tree
                bitmapPosition = bitmapChildren + x * k + y;

                // Calculate base position (row and column)
                thisBaseCol = baseCol + y * divKLevel;

                // Get max value
                if (y == limit1Y && !firstMin) {
                    maxValue = fatherMaxValue -
                               this->listMaxValues[level - 1]->access(
                                       bitmapPosition - this->preMaxValues[level - 1]);
                } else {
                    maxValue = fatherMaxValue - this->listMaxValues[level - 1]->next();
                }

                if (valini > maxValue) {
                    // Values out of range
                    if (allCells) {
                        return 0;
                    }
                    firstMin = false;
                    continue;
                }

                uniform = level == isLeaf || !this->bitmapK2Tree->access(bitmapPosition);
                if (!uniform) {
                    // If is not uniform (or leaf) get his min value
                    if (!firstMin) {
                        countOnes = this->bitmapK2Tree->rank1(bitmapPosition - 1);
                        minValue = fatherMinValue +
                                   this->listMinValues[level - 1]->access(countOnes - this->preMinValues[level - 1]);
                        firstMin = true;
                    } else {
                        countOnes++;
                        minValue = fatherMinValue + this->listMinValues[level - 1]->next();
                    }
                } else {
                    minValue = maxValue;
                }

                if (valend < minValue) {
                    // Values out of range
                    if (allCells) {
                        return 0;
                    }
                    continue;
                }

                uint newxini = max(thisBaseRow, xini);
                uint newxend = min(thisBaseRow + divKLevel - 1, xend);
                uint newyini = max(thisBaseCol, yini);
                uint newyend = min(thisBaseCol + divKLevel - 1, yend);

                if (valini <= minValue && valend >= maxValue) {
                    // Add all valid positions to final result
                    for (ulong x = newxini; x <= newxend; x++) {
                        for (ulong y = newyini; y <= newyend; y++) {
                            positions[totalCells++] = x * this->realSizeY + y;
                        }
                    }
                } else {
                    // Search for their children and they are added to queue
                    // If it is uniform, the node has not children
                    if (!uniform) {
                        if (this->isLevelInPlain(level + 1)) {
                            // Submatrix represented by plain values
                            totalCells = this->getPlainWindow(newxini, newyini, newxend, newyend, valini, valend,
                                                              countOnes - this->preMinValues[level - 1], divKLevel,
                                                              maxValue, positions, totalCells, allCells);
                        } else {
                            // Check if the MBR overlap completely with the submatrix
                            if (allCells && newxini == thisBaseRow && newxend == thisBaseRow + divKLevel - 1
                                && newyini == thisBaseCol && newyend == thisBaseCol + divKLevel - 1) {
                                // Some value is without of the range of values
                                return 0;
                            }

                            // Continue with the division
                            uint nextBitmapChildren = 0;
                            // Calculate children position at bitmap
                            if (level < this->levelK1) {
                                nextBitmapChildren = countOnes * this->K1 * this->K1 + 1;
                            } else {
                                nextBitmapChildren = (this->countOnesk1 * this->K1 * this->K1) +
                                                     ((countOnes - this->countOnesk1) * this->K2 * this->K2 + 1);
                            }

                            totalCells = this->getCellsByValue(newxini, newxend, newyini, newyend, minValue, maxValue,
                                                               level, nextBitmapChildren, valini, valend,
                                                               thisBaseRow, thisBaseCol, totalCells, positions,
                                                               allCells);
                            if (allCells && totalCells == 0) {
                                return 0;
                            }
                        } // END IF plainLevels

                        if (allCells && totalCells == 0) {
                            return 0;
                        }
                    } // END IF uniform
                } // END IF continue search
            } // ENF FOR y
        } // END FOR x
        return totalCells;
    }

    //*************************//
    //****** CHECKVALUES ******//
    //*************************//
    bool K2RasterPlain::checkValues(Query query, bool allCells) const {
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
            oldOffsetValues = -1;

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

                if (this->isLevelInPlain(node->level + 2)) {
                    // Search in plain values
                    uint xini = max(thisbaserow, query.xini);
                    uint xend = min(thismaxrow, query.xend);
                    uint yini = max(thisbasecol, query.yini);
                    uint yend = min(thismaxcol, query.yend);
                    bool value = this->checkPlainCells(xini, yini, xend, yend,
                                                      offsetValues, kLevel, maxValue,
                                                       query, allCells);
                    if (value && !allCells) {
                        return true;
                    } else {
                        if (!value && allCells) {
                            return false;
                        }
                    }
                    continue;
                } else {
                    // Search his children
                    queue.push(this->createChild(node, offsetValues, thisbaserow, thisbasecol, minValue, maxValue));
                }
            }
            free(node);
        }
        /**
         * If allcells is true -> not found any cell with a not valid value
         * if allcell is false -> not found any cell with a valid value
         */
        return allCells;
    }

    //******************************//
    //****** GETVALUESWINDOW ******//
    //*****************************//
    ulong K2RasterPlain::getValuesWindow(uint xini, uint xend, uint yini, uint yend, int *values) const {

        int maxValue, minValue;
        int thisbaserow, thisbasecol, thismaxrow, thismaxcol;
        int newOffsetBitmap = 0, offsetValues, oldOffsetBitmap;
        bool uniform, isLeaf;
        uint totalValue = 0, cellPos;
        uint kLevel = 0;
        uint newK;

        //Whole matrix is uniform
        if (!this->bitmapK2Tree->access(0)){
            for (uint x = 0; x < (xend - xini) + 1; x++) {
                for (uint y = 0; y < (yend - yini) + 1; y++) {
                    values[x * (yend - yini) + y] = this->initialMaxValue;
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
                                                  xini, xend, yini, yend, newK, kLevel);
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
                    uint localXini = max((uint)thisbaserow, xini) - xini;
                    uint localXend = min((uint)thismaxrow, xend) - xini;
                    uint localYini = max((uint)thisbasecol, yini) - yini;
                    uint localYend = min((uint)thismaxcol, yend) - yini;

                    // Add all valid positions to final result
                    for (uint x = localXini; x <= localXend; x++){
                        cellPos = (x * ((yend - yini) + 1));
                        for (uint y = localYini; y <= localYend; y++){
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

                if (this->isLevelInPlain(node->level + 2)) {
                    // Search in plain values
                    uint localXini = max((uint)thisbaserow, xini);
                    uint localXend = min((uint)thismaxrow, xend);
                    uint localYini = max((uint)thisbasecol, yini);
                    uint localYend = min((uint)thismaxcol, yend);

                    totalValue += this->getPlainCellValues(localXini, localYini, localXend, localYend,
                                                          offsetValues, kLevel, maxValue,
                                                          values, totalValue,
                                                           xini, yini, yend);
                    continue;
                } else {
                    // Search his children
                    queue.push(this->createChild(node, offsetValues, thisbaserow, thisbasecol, minValue, maxValue));
                }
            }

            free(node);
        }

        return  totalValue;
    }


    //**********************************************************************//
    //********************* FILE FUNCTIONS *********************************//
    //**********************************************************************//
    void K2RasterPlain::save(std::ofstream &of) const {
        K2Raster::save(of, K2RASTER_PLAIN);

        // Levelsk2
        saveValue(of, this->levelK2);

        /// Save plain values
        saveValue<uint>(of, this->numOfPlainValues);
        saveValue<uint>(of, this->plainValues, this->numOfPlainValues);


        std::string outFileName = "last_values.bin";
        ofstream outFile(outFileName.data());
        if (!outFile.good()) {
            printf("File %s unable to open\n", outFileName.data());
            return;
        }

        saveValue<uint>(outFile, this->numOfPlainValues);
        saveValue<uint>(outFile, this->plainValues, this->numOfPlainValues);
        outFile.close();


    }

    void K2RasterPlain::saveSplit(char *filename) const {
        K2Raster::saveSplit(filename, K2RASTER_PLAIN);

        // Args *************************
        {
            std::string outFileName = filename;
            outFileName += ".argsPlain.k2raster";
            ofstream outFile(outFileName.data());
            if (!outFile.good()) {
                printf("File %s unable to open\n", outFileName.data());
                return;
            }

            saveValue(outFile, this->levelK2);

            outFile.close();
        } // END BLOCK args

        // Plain Value *************************
        {
            std::string outFileName = filename;
            outFileName += ".plainValue.k2raster";
            ofstream outFile(outFileName.data());
            if (!outFile.good()) {
                printf("File %s unable to open\n", outFileName.data());
                return;
            }

            saveValue<uint>(outFile, this->numOfPlainValues);
            saveValue<uint>(outFile, this->plainValues, this->numOfPlainValues);

            outFile.close();
        } // END BLOCK plain values
    }


    /*
     * Load Ktreap from file
     */
    K2RasterPlain *K2RasterPlain::load(std::ifstream &in) {
        K2RasterPlain *ktreap = nullptr;

        try {
            ktreap = new K2RasterPlain();
            K2Raster::loadBase(in, ktreap);

            // Load k2 value
            ktreap->levelK2 = loadValue<uint>(in);

            // load plain values
            ktreap->numOfPlainValues = loadValue<uint>(in);
            ktreap->plainValues = loadValue<uint>(in, ktreap->numOfPlainValues);
        } catch (...) {
            return nullptr;
        }

        return ktreap;
    }


    //**********************************************************************//
    //********************* SIZE FUNCTIONS *********************************//
    //**********************************************************************//
    size_t K2RasterPlain::getTotalSize() const {
        size_t totalSize = sizeof(uint)             // Type
                           + sizeof(uint)           // realSize
                           + (sizeof(uint) * 4)     // K1, K2, levelK1, levelK2
                           + (sizeof(uint) * 2);    // initialMaxValue and initialMinValue

        // k2-tree
        totalSize += this->bitmapK2Tree->getSize(); // bitMapsK2tree

        // DACs
        totalSize += this->getSizeDACs(this->listMaxValues, this->numOfMaxDACs); // Max DACs
        totalSize += this->getSizeDACs(this->listMinValues, this->numOfMinDACs); // Min DACs

        //Plain values
        totalSize += sizeof(uint); //numOfPlainValues;
        totalSize += sizeof(uint) * this->numOfPlainValues;
        return totalSize;
    }

    size_t K2RasterPlain::printSize() const {
        printf("RealSize \t-> %lu\n", sizeof(uint));
        printf("K1, K2, levelK1, levelK2 \t-> %lu\n", sizeof(uint) * 4);
        printf("InitialValues (Max and min) \t-> %lu\n", sizeof(uint) * 2);
        printf("BitMapK2Tree \t-> %lu\n", this->bitmapK2Tree->getSize());


        this->printSizeDACs(this->listMaxValues, this->numOfMaxDACs, false);
        this->printSizeDACs(this->listMinValues, this->numOfMinDACs, false);

        printf("PlainValue -> %lu (+ %lu)\n", this->numOfPlainValues * sizeof(uint), sizeof(uint));
        return this->getTotalSize();
    }

    //**********************************************************************//
    //***************** Auxiliary Functions  *******************************//
    //**********************************************************************//
    void K2RasterPlain::initLevels(uint plainLevels) {
        K2Raster::initLevels();

        // Recalculate levelK1 and levelK2
        if (this->nlevels - plainLevels >= this->levelK1) {
            this->levelK2 = this->nlevels - this->levelK1 - plainLevels;
        } else {
            this->levelK2 = 0;
            this->levelK1 = this->nlevels - plainLevels;
        }
        // Calculate the num of DACs for minValue and maxValues
        this->numOfMaxDACs = this->levelK1 + this->levelK2 - 1;
        this->numOfMinDACs = this->levelK1 + this->levelK2 == this->nlevels ? this->nlevels - 2 : this->levelK1 + this->levelK2 - 1;
    }


    bool K2RasterPlain::isLevelInPlain(uint level) const {
        return (level == (this->levelK1 + this->levelK2) && this->levelK1 + this->levelK2 < this->nlevels );
    }

}

