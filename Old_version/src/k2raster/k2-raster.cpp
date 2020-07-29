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
#include <stdexcept>      // std::out_of_range

#include <k2raster/k2-raster.h>
#include <k2raster/plain/k2-raster-plain.h>
#include <k2raster/k2-raster-opt.h>
#include <k2raster/k2-raster-entropy.h>
#include <include/spatialindex/SpatialIndex.h>
#include <k2raster/plain/k2-raster-plain-DAC.h>
#include <k2raster/plain/k2-raster-plain-VByte.h>
#include <k2raster/compressLeaves/k2-raster-CompressLeaves.h>
#include <k2raster/compressLeaves/k2-raster-CompressLeavesH.h>

namespace k2raster_static {

    //**********************************************************************//
    //********************* CONSTRUCTORS ***********************************//
    //**********************************************************************//
    K2Raster::K2Raster() {
        // Initialize
        this->preMinValues = nullptr;
        this->preMaxValues = nullptr;
        this->minValueNoZero = MAXINT;
        this->displacement = 0;
    };

    K2Raster::K2Raster(uint sizeX, uint sizeY, int **data, uint k1, uint k2, uint levelK1) : K2Raster() {
        // Set K values and store real sizes
        this->K1 = k1;
        this->K2 = k2;
        this->levelK1 = levelK1;
        this->realSizeX = sizeX;
        this->realSizeY = sizeY;

        // Initialize levels with new values of K
        this->initLevels();

        // Execute algorithm of construction
        this->build(sizeX, sizeY, data);
    }

    //*********************//
    //****** ALGEBRA ******//
    //*********************//
    K2Raster::K2Raster(uint k1, uint k2, uint levelK1, K2Raster *k2raster1, K2Raster *k2raster2,
                       OperationRaster operation, bool deleteInputs) : K2Raster() {

        // Set K values and store real sizes
        this->K1 = k1;
        this->K2 = k2;
        this->levelK1 = levelK1;
        // TODO displacement is a BETA Version
        if (operation == OperationRaster::OPERATION_SUBT &&
                k2raster1->initialMinValue < k2raster2->initialMaxValue) {
            // The displacement is used to avoid negative numbers.
            this->displacement = (uint)abs((int)(k2raster1->initialMinValue - k2raster2->initialMaxValue));
        }

        // Size of k2raster1 and size of k2raster2 must be the same
        this->realSizeX = min(k2raster1->realSizeX, k2raster2->realSizeX);
        this->realSizeY = min(k2raster1->realSizeY, k2raster2->realSizeY);

        // Initialize levels  with new values of K
        this->initLevels();

        //buildOperation
        this->buildOperation(k2raster1, k2raster2, operation, deleteInputs);

    }

    // Equal K's
    K2Raster::K2Raster(K2Raster *k2raster1, K2Raster *k2raster2, OperationRaster operation, bool deleteInputs) {
        // Initialize
        this->preMinValues = nullptr;
        this->preMaxValues = nullptr;
        this->minValueNoZero = MAXINT;

        // Set K values and store real sizes
        this->K1 = k2raster1->K1;
        this->K2 = k2raster1->K2;
        this->levelK1 = k2raster1->levelK1;
        // TODO displacement is a BETA Version
        if (operation == OperationRaster::OPERATION_SUBT &&
            k2raster1->initialMinValue < k2raster2->initialMaxValue) {
            // The displacement is used to avoid negative numbers.
            this->displacement = (uint)abs((int)(k2raster1->initialMinValue - k2raster2->initialMaxValue));
        }

        // Size of k2raster1 and size of k2raster2 must be the same
        this->realSizeX = k2raster1->realSizeX;
        this->realSizeY = k2raster1->realSizeY;

        // Initialize levels  with new values of K
        this->initLevels();

        this->buildOperationEqualKs(k2raster1, k2raster2, operation, deleteInputs);
    }




    //************************//
    //****** DESTRUCTOR ******//
    //************************//
    K2Raster::~K2Raster() {
        free(this->divKLevels);
        delete(this->bitmapK2Tree);
        if (this->listMaxValues != nullptr) {
            for (uint i = 0; i < this->numOfMaxDACs; i++){
                delete(this->listMaxValues[i]);
            }
            free(this->listMaxValues);
        }
        if (this->preMaxValues != nullptr) {
            free(this->preMaxValues);
        }
        if (this->listMinValues != nullptr) {
            for (uint i = 0; i < this->numOfMinDACs; i++){
                delete(this->listMinValues[i]);
            }
            free(this->listMinValues);
        }
        if (this->preMinValues != nullptr) {
            free(this->preMinValues);
        }
    }


    //**********************************************************************//
    //******************** DECOMPRESSION ***********************************//
    //**********************************************************************//
    int* K2Raster::decompress() {
//         Malloc memory
        int *data = new int[this->realSizeX * this->realSizeY];
//        this->getValuesWindow(0, this->realSizeX-1, 0, this->realSizeY-1, data);
//        return data;

        //Whole matrix is uniform
        if (!this->bitmapK2Tree->access(0)){
            for (uint x = 0; x < this->realSizeX; x++) {
                for (uint y = 0; y < this->realSizeY; y++) {
                    data[x * this->realSizeY + y] = this->initialMaxValue;
                }
            }
            return data;
        }

        uint *countOnes = new uint[this->numOfMaxDACs]();
        this->decompress(1, 1, this->initialMaxValue, 0, 0, countOnes, data);
        delete[] countOnes;
        return data;
    }

    void K2Raster::decompress(uint currentLevel, uint childrenPosition, uint maxValueParent, uint baseRow, uint baseCol, uint *countOnes, int *data){
        uint thisBaseRow = baseRow, thisBaseCol = baseCol;  // Base position of the parent submatrix
        uint size = this->divKLevels[currentLevel];     //Size of the following submatrices (next level)
        uint newK = this->getK(currentLevel-1);               // K for next level

        uint maxValue;
        ulong ones = countOnes[currentLevel - 1] + this->preMinValues[currentLevel - 1];
        ulong nextChildrenPosition;
        ulong jumpNextChildren;
        if (currentLevel < this->levelK1) {
            jumpNextChildren = this->K1 * this->K1;
            nextChildrenPosition = (ones) * jumpNextChildren + 1;
        } else {
            jumpNextChildren = this->K2 * this->K2;
            nextChildrenPosition = (this->countOnesk1 * this->K1 * this->K1) +
                                   ((ones - this->countOnesk1) * jumpNextChildren + 1);
        }

        for (uint i = 0; i < newK; i++) {
            for (uint j = 0; j < newK; j++) {
                if (childrenPosition == this->preMaxValues[currentLevel-1]){
                    maxValue = maxValueParent - this->listMaxValues[currentLevel - 1]->access(0);
                } else {
                    maxValue = maxValueParent - this->listMaxValues[currentLevel - 1]->next();
                }

                if (currentLevel == this->nlevels-1 || this->bitmapK2Tree->access(childrenPosition) == 0) {
                    // Reach an uniform node, store the corresponding values.
                    uint xend = min(thisBaseRow + size, this->realSizeX);
                    uint yend = min(thisBaseCol + size, this->realSizeY);
                    for (uint x = thisBaseRow; x < xend; x++) {
                        for (uint y = thisBaseCol; y < yend; y++) {
                            data[x * this->realSizeY + y] = maxValue - this->displacement;
                        } // END FOR y
                    } // END FOR x

                } else {
                    // Need to go down one level
                    countOnes[currentLevel-1]++;
                    this->decompress(currentLevel +1, nextChildrenPosition, maxValue, thisBaseRow, thisBaseCol, countOnes, data);
                    nextChildrenPosition += jumpNextChildren;
                }
                // Next child
                childrenPosition++;
                thisBaseCol += size;
            }
            thisBaseCol = baseCol;
            thisBaseRow += size;
        }
    }

    //**********************************************************************//
    //********************* CELL FUNCTIONS *********************************//
    //**********************************************************************//

    //*********************//
    //****** GETCELL ******//
    //*********************//

    int K2Raster::getCell(uint row, uint col) const {
        uint lastValueMax = 0;
        uint thisValueMax;
        uint posBitmap = 0, posBitmapPrev = this->getK(0) * this->getK(0);
        uint divK, offset;
        bool uniform;
        uint totalMinValues, totalMinValuesPrev = 1;

        if (row > this->realSizeX || col > this->realSizeY) {
            return 0;
        }

        // Check firts level
        if (!this->bitmapK2Tree->access(0)){
            // All cell of the matrix are equal
            return this->initialMaxValue - this->displacement;
        }

        lastValueMax = this->initialMaxValue;   // Parent maximum value
        posBitmap = 1;                          // Position in k2-tree
        divK = this->divKLevels[0];             // Size of partition of submatrix
        row = row % divK;                       // Relative row in submatrix
        col = col % divK;                       // Relative column  in submatrix
        uint prevNumber = 1;                    // Number total of maximum numbers in previous levels


        // Search from level 1 to level l - 1
        for (uint l = 1; l < this->nlevels - 1; l++) {
            divK = this->divKLevels[l];
            offset = (row / divK) * this->getK(l - 1) + col / divK;     // Relative position
            posBitmap += offset;                                        // new position in k2-tree

            // Checks if it is uniform and get this maxValue
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
                return thisValueMax- this->displacement;
            }

            // -------------------------------------------------------------------------------------------------- //
            // Update counters previous numbers.
            // Those counters are necessary to speed up the query and because
            // the value of "k" can be different in each level, else it is more difficult calculate the position
            // of the next submatrix
            // -------------------------------------------------------------------------------------------------- //
            prevNumber += this->listMaxValues[l - 1]->getRealLength();
            posBitmapPrev  += this->listMinValues[l -1]->getRealLength() * this->getK(l) * this->getK(l); // Position of next level in k2-tree
            totalMinValuesPrev += this->listMinValues[l - 1]->getRealLength();  // number of minValues == number of 1s in k2-tree
        }

        // Search last level (leaves)
        divK = this->divKLevels[this->nlevels - 1];
        offset = (row / divK) * this->getK(this->nlevels - 2) + col / divK;
        posBitmap += offset;
        posBitmap -= prevNumber;

        thisValueMax = this->listMaxValues[this->nlevels - 2]->access(posBitmap);
        thisValueMax = lastValueMax - thisValueMax;
       return thisValueMax - this->displacement;
    }

    //*****************************//
    //****** GETCELLBYVALUE ******//
    //****************************//

    ulong K2Raster::getCellsByValue(Query query, uint *positions) const {

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
            while (newOffsetBitmap != -1) {
                newOffsetBitmap = this->nextChild(node, posI, posJ, thisbaserow, thismaxrow, thisbasecol, thismaxcol,
                                                  query, newK, kLevel);
                if (newOffsetBitmap == -1) {
                    break;
                }
//                printf("Searching form [%u, %u] to [%u, %u]\n", thisbaserow, thisbasecol, thismaxrow, thismaxcol);

                // Get max value
                maxValue = node->maxValue -
                           this->listMaxValues[node->level]->access(newOffsetBitmap - node->preValuesMax);
                if (query.valini > maxValue){
                    // Values out of range
                    continue;
                }

                uniform = isLeaf || !this->bitmapK2Tree->access(newOffsetBitmap);
                if (!uniform){
                    // If is not uniform (or leaf) get his min value
                    offsetValues = this->bitmapK2Tree->rank1(newOffsetBitmap - 1) - node->preValuesMin;
                    minValue = node->minValue + this->listMinValues[node->level]->access(offsetValues);
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

    ulong K2Raster::getCellsByValue(uint xini, uint xend, uint yini, uint yend, uint valini, uint valend,
                                    BasicSubMatrix *minSubMatrix, uint *positions, bool allCells) const {
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

        // Calculate initial base position
        uint baseRow = (xini / divKLevel) * divKLevel;
        uint baseCol = (yini / divKLevel) * divKLevel;

        if (minSubMatrix->bitmapChildren == 0) {
            ulong countOnes = this->bitmapK2Tree->rank1(minSubMatrix->bitmapPosition - 1);
            if (minSubMatrix->level < this->levelK1) {
                minSubMatrix->bitmapChildren = (countOnes) * this->K1 * this->K1 + 1;
            } else {
                minSubMatrix->bitmapChildren = (this->countOnesk1 * this->K1 * this->K1) +
                                            ((countOnes - this->countOnesk1) * this->K2 * this->K2 + 1);
            }
        }

        return this->getCellsByValue(xini, xend, yini, yend, minSubMatrix->minValue, minSubMatrix->maxValue,
                                     minSubMatrix->level, minSubMatrix->bitmapChildren, valini, valend,
                                     baseRow, baseCol, 0, positions, allCells);
    }


    ulong
    K2Raster::getCellsByValue(uint xini, uint xend, uint yini, uint yend, uint fatherMinValue, uint fatherMaxValue,
                              uint l, uint bitmapChildren, uint valini, uint valend,
                              uint baseRow, uint baseCol,
                              ulong totalCells, uint *positions, bool allCells) const {

        uint divKLevel = this->divKLevels[l + 1];      // Size of its children
        uint k = this->getK(l);
        uint bitmapPosition;
        int minValue, maxValue;
        bool uniform;
        ulong countOnes;
        uint thisBaseRow, thisBaseCol;
        uint level = l + 1;

        for (uint x = (xini - baseRow) / divKLevel; x <= (xend - baseRow) / divKLevel; x++) {
            thisBaseRow = baseRow + x * divKLevel;
            for (uint y = (yini - baseCol) / divKLevel; y <= (yend - baseCol) / divKLevel; y++) {
                // Get position at k2-tree
                bitmapPosition = bitmapChildren + x * k + y;
                thisBaseCol = baseCol + y * divKLevel;

                // Get max value
                if (y == (yini - baseCol) / divKLevel) {
                    maxValue = fatherMaxValue -
                               this->listMaxValues[level - 1]->access(bitmapPosition - this->preMaxValues[level - 1]);
                } else {
                    maxValue = fatherMaxValue - this->listMaxValues[level - 1]->next();
                }

                if (valini > maxValue) {
                    // Values out of range
                    if (allCells) {
                        return 0;
                    }
                    continue;
                }

                uniform = level == this->nlevels - 1 || !this->bitmapK2Tree->access(bitmapPosition);
                if (!uniform) {
                    // If is not uniform (nor leaf) get his min value
                    countOnes = this->bitmapK2Tree->rank1(bitmapPosition - 1);
                    minValue = fatherMinValue +
                               this->listMinValues[level - 1]->access(countOnes - this->preMinValues[level - 1]);
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
                        uint nextBitmapChildren = 0;
                        // Calculate children position at bitmap
                        if (level < this->levelK1) {
                            nextBitmapChildren = countOnes * this->K1 * this->K1 + 1;
                        } else {
                            nextBitmapChildren = (this->countOnesk1 * this->K1 * this->K1) +
                                                 ((countOnes - this->countOnesk1) * this->K2 * this->K2 +
                                                  1);
                        }

                        totalCells = this->getCellsByValue(newxini, newxend, newyini, newyend, minValue, maxValue,
                                                           level, nextBitmapChildren, valini, valend,
                                                           thisBaseRow, thisBaseCol, totalCells, positions, allCells);
                        if (allCells && totalCells == 0) {
                            return 0;
                        }
                    } // END IF uniform
                } // END IF continue search
            } // ENF FOR y
        } // END FOR x
        return totalCells;
    }


    //*****************************//
    //****** GETVALUESWINDOW ******//
    //*****************************//
    ulong K2Raster::getValuesWindow(Query query, int *values) const {
        return this->getValuesWindow(query.xini, query.xend, query.yini, query.yend, values);
    }

    ulong K2Raster::getValuesWindow(uint xini, uint xend, uint yini, uint yend, int *values) const {
        int maxValue, minValue;
        int thisbaserow, thisbasecol, thismaxrow, thismaxcol;
        int newOffsetBitmap = 0, offsetValues,oldOffsetBitmap;
        bool uniform, isLeaf;
        ulong totalValue = 0;
        uint kLevel = 0, newK, cellPos;

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

                if (uniform){
                    // All cell values within range
                    uint localXini = max((uint)thisbaserow, xini) - xini;
                    uint localXend = min((uint)thismaxrow, xend) - xini;
                    uint localYini = max((uint)thisbasecol, yini) - yini;
                    uint localYend = min((uint)thismaxcol, yend) - yini;

                    // Add all valid positions to final result
                    for (uint x = localXini; x <= localXend; x++){
                        cellPos = (x * ((yend - yini) + 1));
                        for (uint y= localYini; y <= localYend; y++){
                            values[cellPos + y] = maxValue ;
                            totalValue++;
                        }
                    }
                    continue;
                }

                offsetValues = this->bitmapK2Tree->rank1(newOffsetBitmap - 1) - node->preValuesMin;

                // Search his children
                queue.push(this->createChild(node, offsetValues, thisbaserow, thisbasecol, minValue, maxValue));
            }
            delete(node);
        }
        return  totalValue;
    }

    //*************************//
    //****** CHECKVALUES ******//
    //*************************//
    bool K2Raster::checkValues(Query query, bool allCells) const {
        int maxValue, minValue;
        int thisbaserow, thisbasecol, thismaxrow, thismaxcol;
        int newOffsetBitmap, offsetValues;
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

            // Search his children
            while (newOffsetBitmap != -1) {
                newOffsetBitmap = this->nextChild(node, posI, posJ, thisbaserow, thismaxrow, thisbasecol, thismaxcol,
                                                  query, newK, kLevel);
                if (newOffsetBitmap == -1) {
                    break;
                }
                maxValue = node->maxValue -
                           this->listMaxValues[node->level]->access(newOffsetBitmap - node->preValuesMax);

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

                offsetValues = this->bitmapK2Tree->rank1(newOffsetBitmap - 1) - node->preValuesMin;
                minValue = node->minValue + this->listMinValues[node->level]->access(offsetValues);

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
        return allCells;
    }

    //*****************************//
    //***** GETMAXVALUEWINDOW *****//
    //*****************************//
    int K2Raster::getMaxValueWindow(Query query) {
        return this->getMaxValueWindow(query.xini, query.xend, query.yini, query.yend);
    }

    int K2Raster::getMaxValueWindow(uint xini, uint xend, uint yini, uint yend) {
        BasicSubMatrix *submatrix = new BasicSubMatrix; 
        this->initSubMatrix(submatrix, this->initialMinValue, this->initialMaxValue);
        int value = this->getMaxValueWindow(xini, xend, yini, yend, 0, 0, submatrix);
        delete submatrix;
        return value;
    }

    int K2Raster::getMaxValueWindow(uint xini, uint xend, uint yini, uint yend, BasicSubMatrix *BasicSubMatrix) {


        uint size = this->divKLevels[BasicSubMatrix->level];

        // Calculate initial base position
        uint baseRow = (xini / size) * size;
        uint baseCol = (yini / size) * size;
        return this->getMaxValueWindow(xini, xend, yini, yend, baseRow, baseCol, BasicSubMatrix) - this->displacement;
    }

    int K2Raster::getMaxValueWindow(uint xini, uint xend, uint yini, uint yend,
                                    uint baseRow, uint baseCol, BasicSubMatrix *submatrix) {


        if (submatrix->level == nlevels - 1) {
            return submatrix->maxValue;
        }

        // Check if the nodes is uniform
        if (!this->bitmapK2Tree->access(submatrix->bitmapPosition)) {
            return submatrix->maxValue;
        }

        // Submatrix size
        uint sizeParent = this->divKLevels[submatrix->level];
        // The current submatrix completely overlaps with the window
        if (xini <= baseRow && yini <= baseCol
            && xend >= (baseRow + sizeParent)
            && yend >= (baseCol + sizeParent)) {
            return submatrix->maxValue;
        }

        // Adjust the window
        uint localXini = max(baseRow, xini);
        uint localXend = min(baseRow + sizeParent - 1, xend);
        uint localYini = max(baseCol, yini);
        uint localYend = min(baseCol + sizeParent - 1, yend);

        // Go down one level
        uint divKSize = this->divKLevels[submatrix->level+1];
        uint childXini = (localXini - baseRow) / divKSize;
        uint childXend = (localXend - baseRow) / divKSize;
        uint childYini = (localYini - baseCol) / divKSize;
        uint childYend = (localYend - baseCol) / divKSize;


        // Child params
        uint thisBaseRow = baseRow + childXini * divKSize;
        uint thisBaseCol = baseCol + childYini * divKSize;

        if (this->isLevelInPlain(submatrix->level + 1)) {
            // TODO implement
        }

        // Create a copy with child's information
        BasicSubMatrix *childSubmatrix = new BasicSubMatrix; 
        childSubmatrix->level = submatrix->level +1;
        childSubmatrix->bitmapChildren = 0;
        ulong countOnes = this->bitmapK2Tree->rank1(submatrix->bitmapPosition - 1);
        uint k = submatrix->level < this->levelK1 ? this->K1 : this->K2;

//        if (submatrix->bitmapChildren == 0) {
            if (submatrix->level < this->levelK1) {
                submatrix->bitmapChildren = (countOnes) * this->K1 * this->K1 + 1;
            } else {
                submatrix->bitmapChildren = (this->countOnesk1 * this->K1 * this->K1) +
                                            ((countOnes - this->countOnesk1) * this->K2 * this->K2 + 1);
            }
//        }

        int maxValueWindow = INT_MIN;
        int value;
        for (uint x = childXini; x <= childXend; x++) {
            for (uint y = childYini; y <= childYend; y++) {
                // Get position at k2-tree
                childSubmatrix->bitmapPosition = submatrix->bitmapChildren + x * k + y;

                // Get max value
                if (y == childYini) {
                    childSubmatrix->maxValue = submatrix->maxValue -
                               this->listMaxValues[childSubmatrix->level - 1]->access(childSubmatrix->bitmapPosition - this->preMaxValues[childSubmatrix->level - 1]);
                } else {
                    childSubmatrix->maxValue = submatrix->maxValue - this->listMaxValues[childSubmatrix->level - 1]->next();
                }

                if (childSubmatrix->level != (this->nlevels - 1)) {
                    value = this->getMaxValueWindow(xini, xend, yini, yend, thisBaseRow, thisBaseCol, childSubmatrix);
                } else {
                    // Leaf node
                    value = childSubmatrix->maxValue - this->displacement;
                }

                if (value == submatrix->maxValue) {
                    // There is no higher value than "value" in the rest of the matrix
                    delete(childSubmatrix);
                    return value;
                }

                if (value > maxValueWindow) {
                    maxValueWindow = value;
                }

                thisBaseCol += divKSize;
            } // END FOR y
            thisBaseCol = baseCol + childYini * divKSize;
            thisBaseRow += divKSize;
        } // END FOR x

        delete(childSubmatrix);
        return maxValueWindow;
    }

    //*****************************//
    //***** GETMINVALUEWINDOW *****//
    //*****************************//

    int K2Raster::getMinValueWindow(Query query) {
        return this->getMinValueWindow(query.xini, query.xend, query.yini, query.yend);
    }

    int K2Raster::getMinValueWindow(uint xini, uint xend, uint yini, uint yend) {
        BasicSubMatrix *submatrix = new BasicSubMatrix;
        this->initSubMatrix(submatrix, this->initialMinValue, this->initialMaxValue);
        int value = this->getMaxValueWindow(xini, xend, yini, yend, 0, 0, submatrix);
        delete submatrix;
        return value;
    }

    int K2Raster::getMinValueWindow(uint xini, uint xend, uint yini, uint yend, BasicSubMatrix *BasicSubMatrix) {


        uint size = this->divKLevels[BasicSubMatrix->level];

        // Calculate initial base position
        uint baseRow = (xini / size) * size;
        uint baseCol = (yini / size) * size;
        return this->getMaxValueWindow(xini, xend, yini, yend, baseRow, baseCol, BasicSubMatrix) - this->displacement;
    }

    int K2Raster::getMinValueWindow(uint xini, uint xend, uint yini, uint yend,
                                    uint baseRow, uint baseCol, BasicSubMatrix *submatrix) {


        if (submatrix->level == nlevels - 1) {
            return submatrix->minValue; // maxValue == minValue
        }

        // Check if the nodes is uniform
        if (!this->bitmapK2Tree->access(submatrix->bitmapPosition)) {
            return submatrix->minValue; // maxValue == minValue
        }

        // Submatrix size
        uint sizeParent = this->divKLevels[submatrix->level];
        // The current submatrix completely overlaps with the window
        if (xini <= baseRow && yini <= baseCol
            && xend >= (baseRow + sizeParent)
            && yend >= (baseCol + sizeParent)) {
            return submatrix->minValue;
        }


        uint childXini, childXend, childYini, childYend;
        uint divKSize = this->divKLevels[submatrix->level+1];

        // Adjust the window
        {
            uint localXini = max(baseRow, xini);
            uint localXend = min(baseRow + sizeParent - 1, xend);
            uint localYini = max(baseCol, yini);
            uint localYend = min(baseCol + sizeParent - 1, yend);

            // Go down one level
            childXini = (localXini - baseRow) / divKSize;
            childXend = (localXend - baseRow) / divKSize;
            childYini = (localYini - baseCol) / divKSize;
            childYend = (localYend - baseCol) / divKSize;
        } // END BLOCK Adjust the window


        // Child params
        uint thisBaseRow = baseRow + childXini * divKSize;
        uint thisBaseCol = baseCol + childYini * divKSize;

        // Create a copy with child's information
        BasicSubMatrix *childSubmatrix = new BasicSubMatrix;
        childSubmatrix->level = submatrix->level +1;
        childSubmatrix->bitmapChildren = 0;
        ulong countOnes = this->bitmapK2Tree->rank1(submatrix->bitmapPosition - 1);
        uint k = submatrix->level < this->levelK1 ? this->K1 : this->K2;

//        if (submatrix->bitmapChildren == 0) {
        if (submatrix->level < this->levelK1) {
            submatrix->bitmapChildren = (countOnes) * this->K1 * this->K1 + 1;
        } else {
            submatrix->bitmapChildren = (this->countOnesk1 * this->K1 * this->K1) +
                                        ((countOnes - this->countOnesk1) * this->K2 * this->K2 + 1);
        }
//        }

        int minValueWindow = INT_MAX;
        int value;
        for (uint x = childXini; x <= childXend; x++) {
            for (uint y = childYini; y <= childYend; y++) {
                // Get position at k2-tree
                childSubmatrix->bitmapPosition = submatrix->bitmapChildren + x * k + y;

                // Get min value
                if (y == childYini) {
                    ulong countOnes = this->bitmapK2Tree->rank1(childSubmatrix->bitmapPosition - 1);
                    childSubmatrix->minValue = submatrix->minValue -
                                               this->listMinValues[childSubmatrix->level - 1]->access(countOnes - this->preMinValues[childSubmatrix->level - 1]);
                } else {
                    childSubmatrix->minValue = submatrix->minValue - this->listMinValues[childSubmatrix->level - 1]->next();
                }

                if (childSubmatrix->level != (this->nlevels - 1)) {
                    value = this->getMinValueWindow(xini, xend, yini, yend, thisBaseRow, thisBaseCol, childSubmatrix);
                } else {
                    // Leaf node
                    value = childSubmatrix->minValue;
                }

                if (value == submatrix->minValue) {
                    // There is no higher value than "value" in the rest of the matrix
                    delete(childSubmatrix);
                    return value;
                }

                if (value < minValueWindow) {
                    minValueWindow = value;
                }

                thisBaseCol += divKSize;
            } // END FOR y
            thisBaseCol = baseCol + childYini * divKSize;
            thisBaseRow += divKSize;
        } // END FOR x

        delete(childSubmatrix);
        return minValueWindow;
    }


    //**********************************************************************//
    //******************** CELL AUX FUNCTIONS ******************************//
    //**********************************************************************//
    BasicSubMatrixWithUsed *K2Raster::getMinMatrix(uint xini, uint xend, uint yini, uint yend, BasicSubMatrixWithUsed *minSubMatrix) const {

        uint divKLevel, k;
        uint subXini, subXend, subYini, subYend;

        if (this->isLevelInPlain(minSubMatrix->level + 1)) {
            // We don't need to go down more levels
            return minSubMatrix;
        }

        if (!this->bitmapK2Tree->access(minSubMatrix->bitmapPosition)){
            // Uniform matrix
            return minSubMatrix;
        }

        // Calculate local positions
        {
            divKLevel = this->divKLevels[minSubMatrix->level];
            xini = xini % divKLevel;
            xend = xend % divKLevel;
            yini = yini % divKLevel;
            yend = yend % divKLevel;
        } // END BLOCK  Calculate local positions

        // Calculate which children have some cell within the region of interest
        divKLevel = this->divKLevels[minSubMatrix->level + 1];      // Size of its children
        subXini = xini / divKLevel;
        subXend = xend / divKLevel;
        subYini = yini / divKLevel;
        subYend = yend / divKLevel;
        if (subXini != subXend || subYini != subYend) {
            // The region of interest belongs to several children, return parent
            return minSubMatrix;
        }

        // Copy structure
        BasicSubMatrixWithUsed *newMinMatrix = new BasicSubMatrixWithUsed;
        this->initSubMatrix(newMinMatrix, minSubMatrix->minValue, minSubMatrix->maxValue);
        newMinMatrix->bitmapPosition = minSubMatrix->bitmapPosition;
        newMinMatrix->bitmapChildren = minSubMatrix->bitmapChildren;
        newMinMatrix->usedBy = 0;
        uint initLevel = minSubMatrix->level + 1;

        // Continue with the rest of the levels
        for (uint l = initLevel; l < this->nlevels; l++) {

            newMinMatrix->level = l;
            k = this->getK(l - 1);

            // Go down one level, calculate next child
            xini = xini % divKLevel;
            xend = xend % divKLevel;
            yini = yini % divKLevel;
            yend = yend % divKLevel;

            // Calculate bitmap position and check if the submatrix is uniform or not
            newMinMatrix->bitmapPosition = newMinMatrix->bitmapChildren + subXini * k + subYini;

            // Calculate max value
            newMinMatrix->maxValue -= this->listMaxValues[l - 1]->access(
                    newMinMatrix->bitmapPosition - this->preMaxValues[l - 1]);

//            uniform = (l == this->nlevels - 1) || !this->bitmapK2Tree->access(newMinMatrix->bitmapPosition);
            if ((l == this->nlevels - 1) || !this->bitmapK2Tree->access(newMinMatrix->bitmapPosition)) {  // True if belong to last level or it is an uniform submatrix
                newMinMatrix->minValue = newMinMatrix->maxValue;
                newMinMatrix->bitmapChildren = 0;
                return  newMinMatrix;
            } else {
                // Get min value from DAC
                ulong countOnes = this->bitmapK2Tree->rank1(newMinMatrix->bitmapPosition - 1);
                newMinMatrix->minValue = minSubMatrix->minValue +
                                         this->listMinValues[l - 1]->access(countOnes - this->preMinValues[l - 1]);

                // Get position of its children
                if (l < this->levelK1) {
                    newMinMatrix->bitmapChildren = (countOnes) * this->K1 * this->K1 + 1;
                } else {
                    newMinMatrix->bitmapChildren = (this->countOnesk1 * this->K1 * this->K1) +
                                                   ((countOnes - this->countOnesk1) * this->K2 * this->K2 + 1);
                }
            }

            if (this->isLevelInPlain(l + 1)) {
                // Values of the next level are in plain form
                return newMinMatrix;
            }

            // Check next children
            divKLevel = this->divKLevels[l + 1];      // Size of its children
            subXini = xini / divKLevel;
            subXend = xend / divKLevel;
            subYini = yini / divKLevel;
            subYend = yend / divKLevel;
            if (subXini != subXend || subYini != subYend) {
                // Region belongs to several children, return parent
                return newMinMatrix;
            }
        }
        return newMinMatrix;
    }



    MinSubMatrix *K2Raster::getMinMatrix(uint xini, uint xend, uint yini, uint yend,
                                         uint valini, uint valend, MinSubMatrix *minSubMatrix) const {

        uint divKLevel, k;
        uint subXini, subXend, subYini, subYend;
        bool uniform;

        if (this->isLevelInPlain(minSubMatrix->level + 1) ||                     // Values of the next level are in plain form
            minSubMatrix->color == SubMatrixColor::COLOR_SUBMATRIX_WHITE ||
            minSubMatrix->color == SubMatrixColor::COLOR_SUBMATRIX_BLACK) {
            // We don't need to go down more levels
            return minSubMatrix;
        }

        // Calculate local positions
        {
            divKLevel = this->divKLevels[minSubMatrix->level];
            xini = xini % divKLevel;
            xend = xend % divKLevel;
            yini = yini % divKLevel;
            yend = yend % divKLevel;
        } // END BLOCK  Calculate local positions

        // Calculate which children have some cell within the region of interest
        divKLevel = this->divKLevels[minSubMatrix->level + 1];      // Size of its children
        subXini = xini / divKLevel;
        subXend = xend / divKLevel;
        subYini = yini / divKLevel;
        subYend = yend / divKLevel;
        if (subXini != subXend || subYini != subYend) {
            // The region of interest belongs to several children, return parent
            return minSubMatrix;
        }

        // Copy structure
        MinSubMatrix *newMinMatrix = (MinSubMatrix *) malloc(sizeof(MinSubMatrix));
        newMinMatrix->bitmapPosition = minSubMatrix->bitmapPosition;
        newMinMatrix->bitmapChildren = minSubMatrix->bitmapChildren;
        newMinMatrix->minValue = minSubMatrix->minValue;
        newMinMatrix->maxValue = minSubMatrix->maxValue;
        newMinMatrix->usedBy = 0;
        uint initLevel = minSubMatrix->level + 1;

        if (minSubMatrix->usedBy == 0) {
            // The father matrix is not used anymore, we free it.
            free(minSubMatrix);
        }

        // Continue with the rest of the levels
        for (uint l = initLevel; l < this->nlevels; l++) {

            newMinMatrix->level = l;
            k = this->getK(l - 1);

            // Go down one level, calculate next child
            xini = xini % divKLevel;
            xend = xend % divKLevel;
            yini = yini % divKLevel;
            yend = yend % divKLevel;

            // Calculate bitmap position and check if the submatrix is uniform or not
            newMinMatrix->bitmapPosition = newMinMatrix->bitmapChildren + subXini * k + subYini;
            uniform = (l == this->nlevels - 1) ||
                      (this->bitmapK2Tree->access(newMinMatrix->bitmapPosition) == 0);  // True if belong to last level or it is an uniform submatrix

            // Calculate max value
//            newMinMatrix->preMaxValues += (l - 1 == 0 ? 1 : this->listMaxValues[l - 2]->getRealLength());
            newMinMatrix->maxValue -= this->listMaxValues[l - 1]->access(
                    newMinMatrix->bitmapPosition - this->preMaxValues[l - 1]);


            uint countOnes = 0;
            if (uniform) {
                // Uniform matrix, it has no children
                newMinMatrix->minValue = newMinMatrix->maxValue;
            } else {
                countOnes = this->bitmapK2Tree->rank1(newMinMatrix->bitmapPosition - 1);

                // Calculate min Value
//                newMinMatrix->preMinValues += (l - 1 == 0 ? 1 : this->listMinValues[l - 2]->getRealLength());
                newMinMatrix->minValue += this->listMinValues[l - 1]->access(
                        countOnes - this->preMinValues[l - 1]);
            } // END IF UNIFORM

            // Set submaxtrix color
            if (newMinMatrix->minValue > valend || newMinMatrix->maxValue < valini) {
                // It is an invalid node, value out of range
                newMinMatrix->color = SubMatrixColor::COLOR_SUBMATRIX_WHITE;
                return newMinMatrix;
            } else {
                if (newMinMatrix->minValue >= valini && newMinMatrix->maxValue <= valend) {
                    // All cells lie within the range of value
                    newMinMatrix->color = SubMatrixColor::COLOR_SUBMATRIX_BLACK;
                    return newMinMatrix;
                } else {
                    // We can not determinate if all cells lie within the range of values
                    newMinMatrix->color = SubMatrixColor::COLOR_SUBMATRIX_GREY;
                }
            } // END IF SET COLOR

            // Calculate children position at bitmap
            if (l < this->levelK1) {
                newMinMatrix->bitmapChildren = (countOnes) * this->K1 * this->K1 + 1;
            } else {
                newMinMatrix->bitmapChildren = (this->countOnesk1 * this->K1 * this->K1) +
                                               ((countOnes - this->countOnesk1) * this->K2 * this->K2 + 1);
            }

            if (this->isLevelInPlain(l + 1)) {
                // Values of the next level are in plain form
                return newMinMatrix;
            }

            // Check next children
            divKLevel = this->divKLevels[l + 1];      // Size of its children
            subXini = xini / divKLevel;
            subXend = xend / divKLevel;
            subYini = yini / divKLevel;
            subYend = yend / divKLevel;
            if (subXini != subXend || subYini != subYend) {
                // Region belongs to several children, return parent
                return newMinMatrix;
            }
        }
        return newMinMatrix;
    }

    MinNodeTreeNode *K2Raster::getMinMatrix(uint xini, uint xend, uint yini, uint yend, uint valini, uint valend,
                                              MinNodeTreeNode *initialNode,
                                              map<ulong, MinNodeTreeNode *> &nextMBRs) const {

        // Get initial tree node (where the MBR was stored)
        MinNodeTreeNode *node = initialNode;
        if (node->color == SubMatrixColor::COLOR_SUBMATRIX_WHITE ||
            node->color == SubMatrixColor::COLOR_SUBMATRIX_BLACK ||
            this->isLevelInPlain(node->level + 1)) {
            // We don't need go down more levels
            return node;
        }

        uint divKLevel, k, bitmapPosition;
        uint subXini, subXend, subYini, subYend;

        // Local positions
        {
            divKLevel = this->divKLevels[node->level];
            xini = xini % divKLevel;
            xend = xend % divKLevel;
            yini = yini % divKLevel;
            yend = yend % divKLevel;
        }

        // Go down the tree
        for (uint l = node->level + 1; l < this->nlevels; l++) {

            // Calculate its children that have some cell of the window
            divKLevel = this->divKLevels[l];      // Size of its children
            subXini = xini / divKLevel;
            subXend = xend / divKLevel;
            subYini = yini / divKLevel;
            subYend = yend / divKLevel;
            if (subXini != subXend || subYini != subYend) {
                // it belongs to several children, return parent
                return node;
            }

            // Calculate next child
            xini = xini % divKLevel;
            xend = xend % divKLevel;
            yini = yini % divKLevel;
            yend = yend % divKLevel;

            // Calculate bitmap position
            k = this->getK(l - 1);
            bitmapPosition = node->bitmapChildren + subXini * k + subYini;

            // Search the position in the map
            MinNodeTreeNode *nextNode;
            try {
                nextNode = nextMBRs.at(bitmapPosition);
            } catch (const std::out_of_range &oor) {
                // There is no node for position 'bitmapPosition'
                // Create a new one
                if (l < this->levelK1 - 1) {
                    this->createChildrenNodeMinSubMatrix(node, k, valini, valend, &nextMBRs);
                    nextNode = nextMBRs[bitmapPosition];
                } else {
                    nextNode = this->createNexMinNodeTreeNode(bitmapPosition, l, node, valini, valend);
                    nextMBRs[nextNode->bitmapPosition] = nextNode;
                }

            } // END CATCH at
//            nextNode->used = true;

            // Check color
            if (nextNode->color == SubMatrixColor::COLOR_SUBMATRIX_WHITE ||
                nextNode->color == SubMatrixColor::COLOR_SUBMATRIX_BLACK ||
                this->isLevelInPlain(nextNode->level + 1)) {
                // We don't need go down more levels
                return nextNode;
            }
            node = nextNode;
        } // END FOR level

        return nullptr;
    }

    MinNodeTreeNode *K2Raster::getMinMatrixV2(uint xini, uint xend, uint yini, uint yend, uint valini, uint valend,
                                                MinNodeTreeNode *node) const {

        if (node->color == SubMatrixColor::COLOR_SUBMATRIX_WHITE ||
            node->color == SubMatrixColor::COLOR_SUBMATRIX_BLACK ||
            this->isLevelInPlain(node->level + 1)) {
            // We don't need to go down more levels
            return node;
        }

        uint divKLevel, k;
        ulong bitmapPosition;
        uint subXini, subXend, subYini, subYend;

        // Local positions
        {
            divKLevel = this->divKLevels[node->level];
            xini = xini % divKLevel;
            xend = xend % divKLevel;
            yini = yini % divKLevel;
            yend = yend % divKLevel;
        }

        // Go down the tree
        uint initLevel = node->level + 1;
        for (uint l = initLevel; l < this->nlevels; l++) {

            // Calculate the children that have some cell of the window
            divKLevel = this->divKLevels[l];      // Size of its children
            subXini = xini / divKLevel;
            subXend = xend / divKLevel;
            subYini = yini / divKLevel;
            subYend = yend / divKLevel;
            if (subXini != subXend || subYini != subYend) {
                // it belongs to several children, return parent
                return node;
            }

            // Calculate next child
            xini = xini % divKLevel;
            xend = xend % divKLevel;
            yini = yini % divKLevel;
            yend = yend % divKLevel;

            // Calculate bitmap position
            k = this->getK(l - 1);
            bitmapPosition = node->bitmapChildren + subXini * k + subYini;

            // Search the position
            MinNodeTreeNode *nextNode;
            if (node->children == nullptr) {
                // First access to its children, alloc necessary memory to store the pointers of its children.
                node->children = (MinNodeTreeNode **) malloc(node->numberOfChildren * sizeof(MinNodeTreeNode *));
                for (uint c = 0; c < node->numberOfChildren; c++) {
                    node->children[c] = nullptr;
                }
            }
            if (node->children[subXini * k + subYini] != nullptr) {
                // Child was created before. Get it from the children array
                nextNode = node->children[subXini * k + subYini];
            } else {
                // Create and get the child
                if (l < this->levelK1 - 1) {
                    // Create all children
                    this->createChildrenNodeMinSubMatrix(node, k, valini, valend, nullptr, true);
                    nextNode = node->children[subXini * k + subYini];
                } else {
                    // Only create the request node
                    nextNode = this->createNexMinNodeTreeNode(bitmapPosition, l, node, valini, valend, true);
                    node->children[subXini * k + subYini] = nextNode;
                }
            }
//            nextNode->used = true;

            // Check color
            if (nextNode->color == SubMatrixColor::COLOR_SUBMATRIX_WHITE ||
                nextNode->color == SubMatrixColor::COLOR_SUBMATRIX_BLACK ||
                this->isLevelInPlain(nextNode->level + 1)) {
                // We don't need go down more levels
                return nextNode;
            }
            node = nextNode;
        } // END FOR level

        return nullptr;
    }

    BasicSubMatrix* K2Raster::createRootBasicSubMatrix() const {
        // First node (root node)
        BasicSubMatrix *subMatrix = new BasicSubMatrix;
        this->initSubMatrix(subMatrix, this->initialMinValue, this->initialMaxValue);
        return subMatrix;
    }

    BasicSubMatrixWithUsed* K2Raster::createRootBasicSubMatrixWithUsed() const {
        // First node (root node)
        BasicSubMatrixWithUsed *subMatrix = new BasicSubMatrixWithUsed;
        this->initSubMatrix(subMatrix, this->initialMinValue, this->initialMaxValue);
        subMatrix->usedBy = 0;
        return subMatrix;
    }


    MinSubMatrix *K2Raster::createRootMinSubMatrix(uint valini, uint valend) const {
        // First node (root node)
        MinSubMatrix *subMatrix = (MinSubMatrix *) malloc(sizeof(MinSubMatrix));
        this->initMinSubMatrix(subMatrix, valini, valend);
        subMatrix->usedBy = 0;
        return subMatrix;
    }

    MinNodeTreeNode *K2Raster::createRootNodeMinSubMatrix(uint valini, uint valend, bool addChildren) const {
        // First node (root node)
        MinNodeTreeNode *node = (MinNodeTreeNode *) malloc(sizeof(MinNodeTreeNode));
        this->initMinSubMatrix(node, valini, valend);
        node->queue = nullptr;

        if (addChildren) {
            node->numberOfChildren = this->getK(0) * this->getK(0);
            node->children = nullptr;
        }

        // Create queue
        node->queue = nullptr;
        return node;
    }

    //**********************************************************************//
    //********************* INFO FUNCTIONS *********************************//
    //**********************************************************************//

    uint K2Raster::getMinValue() const {
        return this->initialMinValue;
    }

    uint K2Raster::getMaxValue() const {
        return this->initialMaxValue;
    }

    uint K2Raster::getMinValueNoZero() const {
        return this->minValueNoZero;
    }

    uint K2Raster::getNumOfPositions() const {
        return this->realSizeX * this->realSizeY;
    }

    uint K2Raster::getDimSize(uint dim) const {
        switch (dim) {
            case 1:
                return this->realSizeX;
            case 2:
                return this->realSizeY;
            default:
                return 0;
        }
    }

    std::string K2Raster::getTypeName() const {
        return "k2-raster";
    }

    std::string K2Raster::getTypeName(uint type) {
        switch (type) {
            case K2RASTER:
                return "k2-raster\t\t";
            case K2RASTER_PLAIN:
                return "k2-raster-plain";
            case K2RASTER_PLAIN_DAC:
                return "k2-raster-DAC";
            case K2RASTER_PLAIN_VBYTE:
                return "k2-raster-VByte";
            case K2RASTER_VOC:
                return "k2-raster-voc";
            case K2RASTER_VOC_HYBRID:
                return "k2-raster-hybrid";
            case K2RASTER_OPT:
                return "k2-raster-opt";
            case K2RASTER_VOC_ENTROPY:
                return "k2-raster-entropy";
            default:
                return "k2-raster\t\t";
        } // END SWITCH
    }

    //**********************************************************************//
    //********************* FILE FUNCTIONS *********************************//
    //**********************************************************************//
    void K2Raster::save(std::ofstream &of) const {
        this->save(of, K2RASTER);
    }

    void K2Raster::save(std::ofstream &of, uint type) const {
        saveValue(of, type);
        saveValue(of, this->realSizeX);
        saveValue(of, this->realSizeY);
        saveValue(of, this->K1);
        saveValue(of, this->K2);
        saveValue(of, this->levelK1);
        saveValue(of, this->initialMaxValue);
        saveValue(of, this->initialMinValue);
        saveValue(of, this->minValueNoZero);
        saveValue(of, this->displacement);
        saveValue(of, this->numOfMaxDACs);
        saveValue(of, this->numOfMinDACs);

        // BitMaps
        this->bitmapK2Tree->save(of);

        // DACs
        // Max values (nlevels - 1) (not stored level 0)
        for (uint i = 0; i < this->numOfMaxDACs; i++){
            this->listMaxValues[i]->save(of);
        }
        // Min values (nlevels - 2) (not stores level 0 and last level)
        for (uint i = 0; i < this->numOfMinDACs; i++){
            this->listMinValues[i]->save(of);
        }
    }

    void K2Raster::saveSplit(char *filename) const {
        this->saveSplit(filename, K2RASTER);
    }

    void K2Raster::saveSplit(char *filename, uint type) const {

        // Args *************************
        {
            std::string outFileName = filename;
            outFileName += ".args.k2raster";
            ofstream outFile(outFileName.data());
            if (!outFile.good()) {
                printf("File %s unable to open\n", outFileName.data());
                return;
            }

            saveValue(outFile, type);
            saveValue(outFile, this->realSizeX);
            saveValue(outFile, this->realSizeY);
            saveValue(outFile, this->K1);
            saveValue(outFile, this->K2);
            saveValue(outFile, this->levelK1);
            saveValue(outFile, this->initialMaxValue);
            saveValue(outFile, this->initialMinValue);
            saveValue(outFile, this->minValueNoZero);
            saveValue(outFile, this->displacement);
            saveValue(outFile, this->numOfMaxDACs);
            saveValue(outFile, this->numOfMinDACs);

            outFile.close();
        } // END BLOCK args


        // k2-tree *************************
        {
            std::string outFileName = filename;
            outFileName += ".k2-tree.k2raster";
            ofstream outFile(outFileName.data());
            if (!outFile.good()) {
                printf("File %s unable to open\n", outFileName.data());
                return;
            }

            this->bitmapK2Tree->save(outFile);

            outFile.close();
        } // END BLOCK k2-tree


        // DACs *************************
        // Max DACs
        for (uint i = 0; i < this->numOfMaxDACs; i++) {
            std::string outFileName = filename;
            outFileName += ".DAC-max-";
            outFileName += std::to_string(i);
            outFileName += ".k2raster";
            ofstream outFile(outFileName.data());
            if (!outFile.good()) {
                printf("File %s unable to open\n", outFileName.data());
                return;
            }
            this->listMaxValues[i]->save(outFile);

            outFile.close();
        } // END MAX DACs

        // Min DACs
        for (uint i = 0; i < this->numOfMinDACs; i++) {
            std::string outFileName = filename;
            outFileName += ".DAC-min-";
            outFileName += std::to_string(i);
            outFileName += ".k2raster";
            ofstream outFile(outFileName.data());
            if (!outFile.good()) {
                printf("File %s unable to open\n", outFileName.data());
                return;
            }
            this->listMinValues[i]->save(outFile);

            outFile.close();
        } // END MIN DACs
    }

    /*
     * Load Ktreap from file
     */
    K2Raster *K2Raster::load(std::ifstream &in) {
        uint type = loadValue<uint>(in);
        size_t pos = in.tellg();
        in.seekg(pos-sizeof(uint));

        switch (type) {
            case K2RASTER: {
                K2Raster *raster = new K2Raster();
                K2Raster::loadBase(in, raster);
                return raster;
//                return nullptr;
            }
            case K2RASTER_PLAIN:
                return K2RasterPlain::load(in);
            case K2RASTER_PLAIN_DAC:
                return K2RasterPlainDAC::load(in);
            case K2RASTER_PLAIN_VBYTE:
                return K2RasterPlainVByte::load(in);
            case K2RASTER_VOC:
                return K2RasterCompressLeaves::load(in);
            case K2RASTER_VOC_HYBRID:
                return K2RasterCompressLeavesH::load(in);
            case K2RASTER_OPT:
                return K2RasterOPT::load(in);
            case K2RASTER_VOC_ENTROPY:
                return K2RasterEntropy::load(in);
            default:
                return nullptr;
        }
        return nullptr;
    }

    K2Raster *K2Raster::createK2Raster(uint sizeX, uint sizeY, int **data, uint k1, uint k2, uint levelK1,
                                       uint plainLevels, float minFreqVoc, uint type) {
        K2Raster *raster = nullptr;
        switch (type) {
            case K2RASTER:
                // Normal k2-raster. It divides the first "levelK1" levels by "K1" and the rest by "K2"
                raster = new k2raster_static::K2Raster(sizeX, sizeY, data, k1, k2, levelK1);
                break;
            case K2RASTER_PLAIN:
                // Plain k2-raster. It divides the first "levelK1" levels by "K1", the next "levelK2" levels by "K2"
                // and it stores the rest of the values in plain form (array of integers)
                raster = new k2raster_static::K2RasterPlain(sizeX, sizeY, data, k1, k2, levelK1, plainLevels);
                break;
            case K2RASTER_PLAIN_DAC:
                // Plain DAC k2-raster. It divides the first "levelK1" levels by "K1", the next "levelK2" levels by "K2"
                // and it stores the rest of the values encoded with DAC
                raster = new k2raster_static::K2RasterPlainDAC(sizeX, sizeY, data, k1, k2, levelK1, plainLevels);
                break;
            case K2RASTER_PLAIN_VBYTE:
                // Plain VByte k2-raster. It divides the first "levelK1" levels by "K1", the next "levelK2" levels by "K2"
                // and it stores the rest of the values encoded with VByte
                raster = new k2raster_static::K2RasterPlainVByte(sizeX, sizeY, data, k1, k2, levelK1, plainLevels);
                break;
            case K2RASTER_VOC:
                // Vocabulary k2-raster. It divides the first "levelK1" levels by "K1", the next "levelK2" levels by "K2"
                // and it assigns a word of the vocabulary for each different final submatrix
                raster = new k2raster_static::K2RasterCompressLeaves(sizeX, sizeY, data, k1, k2, levelK1, plainLevels);
                break;
            case K2RASTER_VOC_HYBRID:
                // Plain k2-raster. It divides the first "levelK1" levels by "K1", the next "levelK2" levels by "K2"
                // and if the sumbatrix is frequent, it is stored as a word of vocabulary, in other case, it are stored in plain form
                raster = new k2raster_static::K2RasterCompressLeavesH(sizeX, sizeY, data, k1, k2, levelK1, plainLevels,
                                                                      minFreqVoc);
                break;
            case K2RASTER_OPT:
                // Optimized k2-raster. It divides the first "levelK1" levels by "K1" and the rest by "K2"
                // (with optimized DAC)
                raster = new k2raster_static::K2RasterOPT(sizeX, sizeY, data, k1, k2, levelK1);
                break;
            case K2RASTER_VOC_ENTROPY:
                // Entropy k2-raster. It divides the first "levelK1" levels by "K1", the next "levelK2" levels by "K2"
                // Implemented a entropy function to calculate, for each submatrix, whether it is better to store it in plain form
                // or as a word of vocabulary.
                raster = new k2raster_static::K2RasterEntropy(sizeX, sizeY, data, k1, k2, levelK1, plainLevels);
                break;
            default:
                raster = nullptr;
                break;
        } // END SWITCH
        return raster;
    }

    /*
     * Load Ktreap from file
     */
    K2Raster *K2Raster::loadBase(std::ifstream &in, K2Raster *k2raster) {

        try {
            loadValue<uint>(in); //Ktreap type
            k2raster->realSizeX = loadValue<uint>(in);
            k2raster->realSizeY = loadValue<uint>(in);
            k2raster->K1 = loadValue<uint>(in);
            k2raster->K2 = loadValue<uint>(in);
            k2raster->levelK1 = loadValue<uint>(in);
            k2raster->initialMaxValue = loadValue<uint>(in);
            k2raster->initialMinValue = loadValue<uint>(in);
            k2raster->minValueNoZero = loadValue<uint>(in);
            k2raster->displacement = loadValue<uint>(in);
            k2raster->initLevels();
            k2raster->numOfMaxDACs = loadValue<uint>(in);
            k2raster->numOfMinDACs = loadValue<uint>(in);

            // BitMaps
            k2raster->bitmapK2Tree = cds_static::BitSequence::load(in);

            // DACs
            k2raster->listMaxValues = (cds_static::DirectAccess **) malloc(
                    k2raster->numOfMaxDACs * sizeof(cds_static::DirectAccess *));
            k2raster->preMaxValues = (ulong *) malloc(k2raster->numOfMaxDACs * sizeof(ulong));
            for (uint i = 0; i < k2raster->numOfMaxDACs; i++) {
                k2raster->listMaxValues[i] = cds_static::DAC::load(in);
                k2raster->preMaxValues[i] =
                        i == 0 ? 1 : k2raster->preMaxValues[i - 1] + k2raster->listMaxValues[i - 1]->getRealLength();
            }

            k2raster->listMinValues = (cds_static::DirectAccess **) malloc(
                    k2raster->numOfMinDACs * sizeof(cds_static::DirectAccess *));
            k2raster->countOnesk1 = 1;
            k2raster->preMinValues = (ulong *) malloc(k2raster->numOfMinDACs * sizeof(ulong));
            for (uint i = 0; i < k2raster->numOfMinDACs; i++) {
                k2raster->listMinValues[i] = cds_static::DAC::load(in);
                k2raster->preMinValues[i] =
                        i == 0 ? 1 : k2raster->preMinValues[i - 1] + k2raster->listMinValues[i - 1]->getRealLength();
                if (i + 1 < k2raster->levelK1) {
                    k2raster->countOnesk1 += k2raster->listMinValues[i]->getRealLength();
                }
            }
        } catch (...) {
            return nullptr;
        }

        return k2raster;
    }


    //**********************************************************************//
    //********************* SIZE FUNCTIONS *********************************//
    //**********************************************************************//
    size_t K2Raster::getTotalSize() const {
        size_t totalSize = sizeof(uint)             // Type
                           + sizeof(uint)  * 2      // realSize
                           + sizeof(uint)           // Displacement
                           + (sizeof(uint) * 3)     // K1, K2, levelK1
                           + (sizeof(uint) * 3);    // initialMaxValue, initialMinValue and minValueNoZero

        // k2-tree
        totalSize += this->bitmapK2Tree->getSize(); // bitMapsK2tree

        // Min and Max DACs
        for (uint i = 0; i < this->numOfMaxDACs; i++) {
            totalSize += this->listMaxValues[i]->getSize();
        }
        for (uint i = 0; i < this->numOfMinDACs; i++) {
            totalSize += this->listMinValues[i]->getSize();
        }
        return totalSize;
    }

    size_t K2Raster::printSize() const {
        printf("RealSize \t-> %lu\n", sizeof(uint));
        printf("K1, K2, levelK1 \t-> %lu\n", sizeof(uint) * 3);
        printf("InitialValues (Max, min and minNoZero) \t-> %lu\n", sizeof(uint) * 3);
        printf("BitMapK2Tree \t-> %lu\n", this->bitmapK2Tree->getSize());
        size_t totalDACs = 0;
        for (uint i = 0; i < this->numOfMaxDACs; i++){
//            printf("DAC max level %u -> %lu\n", i + 1, this->listMaxValues[i]->getSize());
            totalDACs += this->listMaxValues[i]->getSize();
        }
        printf("Total DACs max \t -> %lu\n", totalDACs);

        totalDACs = 0;
        for (uint i = 0; i < this->numOfMinDACs; i++){
//            printf("DAC min level %u -> %lu\n", i + 1, this->listMinValues[i]->getSize());
            totalDACs += this->listMinValues[i]->getSize();
        }
        printf("Total DACs min \t -> %lu\n", totalDACs);
        return this->getTotalSize();
    }

    //**********************************************************************//
    //******************** CHECK FUNCTIONS *********************************//
    //**********************************************************************//
    void K2Raster::checkData(uint sizeX, uint sizeY, int **data) {
//        int replacedNegativeNumbers = 0;
        int maxNegative = 0;
        for (uint i = 0; i < sizeX; i++) {
            for (uint j = 0; j < sizeY; j++) {
                if (data[i][j] < 0) {
                    if (data[i][j] < maxNegative){
                        maxNegative = data[i][j];
                    }
                }
            }
        }
        this->displacement = -maxNegative;
//        if (replacedNegativeNumbers) {
//            printf("Warning: found %d negative numbers.\n", replacedNegativeNumbers);
//        }
    }

    bool K2Raster::checkEncode(int **data, uint nrows, uint ncols) const {

        int val;
        for (uint x = 0; x < nrows; x++){
//            printf("Searching position (%u, %u) to (%u, %u)\n", x, 0, x, ncols-1);
            for (uint y = 0; y < ncols; y++){
                if (data[x][y] != -1) {
                    val = this->getCell(x, y);
                    if (val != data[x][y]) {
                        printf("Value in position (%u, %u) -> %i, expected -> %i\n", x, y, val, data[x][y]);
                        return false;
                    }
                }
            }
        }
        return true;
    }

    //**********************************************************************//
    //********************** DACs FUNCTIONS ********************************//
    //**********************************************************************//
    size_t K2Raster::getSizeDACs(cds_static::DirectAccess **dacs, uint levels) const {
        size_t totalDACs = 0;
        for (uint i = 0; i < levels; i++) {
            totalDACs += dacs[i]->getSize();
        }
        return totalDACs;
    }

    size_t K2Raster::printSizeDACs(cds_static::DirectAccess **dacs, uint levels, bool printEach) const {
        size_t totalDACs = 0;
        for (uint i = 0; i < levels; i++){
            if (printEach){
                printf("DAC level %u -> %lu\n", i + 1, dacs[i]->getSize());
            }
            totalDACs += dacs[i]->getSize();
        }
        printf("Total DACs \t -> %lu\n", totalDACs);
        return totalDACs;
    }

    //**********************************************************************//
    //********************* AUX FUNCTIONS **********************************//
    //**********************************************************************//

    void K2Raster::initLevels() {
        // Calculate size, number of levels and other important params
        uint size = max(this->realSizeX, this->realSizeY);
        this->size = 1 << (bits(size - 1)); //FIXME: Only works for K=2
        this->nlevels = bits(this->size);

        // Calculate the number of elements for submatrix and level
        //TODO (store real levelK1) proof with a few levels
        this->divKLevels = (uint *) malloc(this->nlevels * sizeof(uint));
        uint divK = this->size;
        for (uint i = 0; i < this->nlevels; i++){
            this->divKLevels[i] = divK;
//            printf("Level %u -> %u and k=%u\n", i, divK, this->getK(i));
            if (divK == 1){
                this->nlevels = i + 1;
                if (this->nlevels < this->levelK1){
                    this->levelK1 = this->nlevels;
                }
                break;
            }
            if (divK < this->getK(i)){
                this->levelK1 = i;
            }
            divK /= this->getK(i);
        }

        // Calculate the num of DACs for minValue and maxValues
        this->numOfMaxDACs = this->nlevels - 1;
        this->numOfMinDACs = this->nlevels - 2;

        return;
    }

    uint K2Raster::build(uint sizeX, uint sizeY, int **data) {

        // Check source data
        this->checkData(sizeX, sizeY, data);

        // ------------------------------------------------------------------- //
        // Initialize nodes where store max a min values and their children
        // ------------------------------------------------------------------- /
        uint totalNodes = 1; // Total number of nodes in the conceptual tree
        std::vector<uint> **maxValues = (std::vector<uint> **)malloc(this->numOfMaxDACs * sizeof(std::vector<uint> *));
        std::vector<uint> **minValues = (std::vector<uint> **)malloc(this->numOfMinDACs * sizeof(std::vector<uint> *));
        this->tmpPlainValues = new std::vector<uint>();
        uint **tmpBitmap = (uint**)malloc(this->numOfMaxDACs * sizeof(uint*));
        uint *levelPositions = (uint*)malloc(this->numOfMaxDACs * sizeof(uint));

        // Malloc memory to store the temporal k2-tree (for levels)
        {
            uint nBits = 1;
            for (uint l = 0; l < this->numOfMaxDACs; l++) {
                maxValues[l] = new std::vector<uint>();
                nBits *= (this->getK(l) * this->getK(l));
                tmpBitmap[l] = (uint *) calloc(uint_len(1, nBits), sizeof(uint));
                levelPositions[l] = 0;
                if (l < this->numOfMinDACs) {
                    minValues[l] = new std::vector<uint>();
                }
            }
        }

        // ------------------------------------------------------------------- //
        // Run the algorithm of construction over the matrix
        // Return the conceptual tree that represent the k2-raster
        // ------------------------------------------------------------------- //
        NodeMatrix *root = checkSubmatrix(0, 0, this->divKLevels[0], 0, data, totalNodes,
                                          maxValues, minValues, this->tmpPlainValues,
                                          tmpBitmap, levelPositions);

        // ------------------------------------------------------------------- //
        // Transform the conceptual tree in a compact tree
        // ------------------------------------------------------------------- //
        return this->conceptualToCompact(root, totalNodes, maxValues, minValues, tmpBitmap, levelPositions);

    }

    uint K2Raster::buildOperation(K2Raster *k2Raster1, K2Raster *k2Raster2, OperationRaster operation, bool deleteInputs) {
        // ------------------------------------------------------------------- //
        // Initialize nodes where store max a min values and their children
        // ------------------------------------------------------------------- /
        uint totalNodes = 1; // Total number of nodes in the conceptual tree
        std::vector<uint> **maxValues = (std::vector<uint> **)malloc(this->numOfMaxDACs * sizeof(std::vector<uint> *));
        std::vector<uint> **minValues = (std::vector<uint> **)malloc(this->numOfMinDACs * sizeof(std::vector<uint> *));
        this->tmpPlainValues = new std::vector<uint>();
        uint **tmpBitmap = (uint**)malloc(this->numOfMaxDACs * sizeof(uint*));
        uint *levelPositions = (uint*)malloc(this->numOfMaxDACs * sizeof(uint));

        // Malloc memory to store the temporal k2-tree (for levels)
        {
            uint nBits = 1;
            for (uint l = 0; l < this->numOfMaxDACs; l++) {
                maxValues[l] = new std::vector<uint>();
                nBits *= (this->getK(l) * this->getK(l));
                if (l + 1 != this->nlevels - 1) {
                    tmpBitmap[l] = (uint *) calloc(uint_len(1, nBits), sizeof(uint));
                }
                levelPositions[l] = 0;
                if (l < this->numOfMinDACs) {
                    minValues[l] = new std::vector<uint>();
                }
            }
        }

        // ------------------------------------------------------------------- //
        // Run the algorithm of construction over the matrix
        // Return the conceptual tree that represent the k2-raster
        // ------------------------------------------------------------------- //
        BasicSubMatrixWithChildren **ptrK2raster1  = new BasicSubMatrixWithChildren*[1];
        BasicSubMatrixWithChildren *node1 = new BasicSubMatrixWithChildren;
        this->initSubMatrix(node1, k2Raster1->initialMinValue, k2Raster1->initialMaxValue);
        node1->children = nullptr;
        ptrK2raster1[0] = node1;

        BasicSubMatrixWithChildren **ptrK2raster2  = new BasicSubMatrixWithChildren*[1];
        BasicSubMatrixWithChildren *node2 = new BasicSubMatrixWithChildren;
        this->initSubMatrix(node2, k2Raster2->initialMinValue, k2Raster2->initialMaxValue);
        node2->children = nullptr;
        ptrK2raster2[0] = node2;

        NodeMatrix *root = checkSubmatrixOperation(0, 0, 0, totalNodes,
                                                   maxValues, minValues, this->tmpPlainValues, tmpBitmap, levelPositions,
                                                   ptrK2raster1, 1, ptrK2raster2, 1,
                                                   k2Raster1, k2Raster2, operation);
        removeChildren(node1);
        delete(node1);
        delete[](ptrK2raster1);
        removeChildren(node2);
        delete(node2);
        delete[](ptrK2raster2);

        if (deleteInputs) {
            delete k2Raster1;
            delete k2Raster2;
        }
        // ------------------------------------------------------------------- //
        // Transform the conceptual tree in a compact tree
        // ------------------------------------------------------------------- //
        return this->conceptualToCompact(root, totalNodes, maxValues, minValues, tmpBitmap, levelPositions);
    }

    // Equal K's
    uint K2Raster::buildOperationEqualKs(K2Raster *k2Raster1, K2Raster *k2Raster2, OperationRaster operation, bool deleteInputs) {
//        cds_utils::start_timing();

        // ------------------------------------------------------------------- //
        // Initialize nodes where store max a min values and their children
        // ------------------------------------------------------------------- /
        uint totalNodes = 1; // Total number of nodes in the conceptual tree
        std::vector<uint> **maxValues = (std::vector<uint> **)malloc(this->numOfMaxDACs * sizeof(std::vector<uint> *));
        std::vector<uint> **minValues = (std::vector<uint> **)malloc(this->numOfMinDACs * sizeof(std::vector<uint> *));
        this->tmpPlainValues = new std::vector<uint>();
        uint **tmpBitmap = (uint**)malloc(this->numOfMaxDACs * sizeof(uint*));
        uint *levelPositions = (uint*)malloc(this->numOfMaxDACs * sizeof(uint));

        // Malloc memory to store the temporal k2-tree (for levels)
        {
            uint nBits = 1;
            for (uint l = 0; l < this->numOfMaxDACs; l++) {
                maxValues[l] = new std::vector<uint>();
                nBits *= (this->getK(l) * this->getK(l));
                if (l + 1 != this->nlevels - 1) {
                    tmpBitmap[l] = (uint *) calloc(uint_len(1, nBits), sizeof(uint));
                }
                levelPositions[l] = 0;
                if (l < this->numOfMinDACs) {
                    minValues[l] = new std::vector<uint>();
                }
            }
        }

        // ------------------------------------------------------------------- //
        // Run the algorithm of construction over the matrix
        // Return the conceptual tree that represent the k2-raster
        // ------------------------------------------------------------------- //
        BasicSubMatrix *ptrK2raster1  = new BasicSubMatrix; 
        this->initSubMatrix(ptrK2raster1,  k2Raster1->getMinValue(), k2Raster1->getMaxValue());
        uint *countOnes1 = (uint*)calloc(k2Raster1->numOfMaxDACs, sizeof(uint));
        BasicSubMatrix *ptrK2raster2  = new BasicSubMatrix; 
        this->initSubMatrix(ptrK2raster2,  k2Raster2->getMinValue(), k2Raster2->getMaxValue());
        uint *countOnes2 = (uint*)calloc(k2Raster2->numOfMaxDACs, sizeof(uint));

//        double timeInitialization = cds_utils::get_timing();
//        cds_utils::start_timing();
        NodeMatrix *root = checkSubmatrixOperation(0, 0, 0, totalNodes,
                                                   maxValues, minValues, this->tmpPlainValues, tmpBitmap, levelPositions,
                                                   ptrK2raster1, ptrK2raster2, k2Raster1, k2Raster2, countOnes1, countOnes2,
                                                   operation);

        delete(ptrK2raster1);
        free(countOnes1);
        delete(ptrK2raster2);
        free(countOnes2);

        if (deleteInputs) {
            delete k2Raster1;
            delete k2Raster2;
        }
//        double timeConstruction = cds_utils::get_timing();
        // ------------------------------------------------------------------- //
        // Transform the conceptual tree in a compact tree
        // ------------------------------------------------------------------- //
//        cds_utils::start_timing();
        uint result = this->conceptualToCompact(root, totalNodes, maxValues, minValues, tmpBitmap, levelPositions);
//        double timeConceptualToCompact = cds_utils::get_timing();
//        printf("Time: Init = %.3f seconds || Construction =  %.3f seconds || Compact =  %.3f seconds\n",
//               timeInitialization/1000.0, timeConstruction/1000.0, timeConceptualToCompact/1000.0);
        return result;
    }

    NodeMatrix* K2Raster::checkSubmatrix(uint baseRow, uint baseCol, uint subMatrixSize,
                                         uint level, int **data, uint &totalNodes,
                                         std::vector<uint> **maxValues, std::vector<uint> **minValues,
                                         std::vector<uint> *plainValues, uint **tmpBitmap, uint *levelPosition) {
        if (subMatrixSize == 1){
            // --------------------------------------------------- //
            // One cell, last level.
            // MinValue is equal to MaxValue and has not children
            // --------------------------------------------------- //
            NodeMatrix *childNode = new NodeMatrix();
            childNode->minValue = data[baseRow][baseCol] + this->displacement;
            childNode->maxValue = childNode->minValue;

            if (childNode->minValue != 0 && childNode->minValue < this->minValueNoZero) {
                this->minValueNoZero = childNode->minValue;
            }
            return childNode;
        } else {
            // Search submatrix
            uint pos = 0;
            uint thisBaseRow, thisBaseCol;
            NodeMatrix *node = new NodeMatrix();
            node->minValue = MAXINT;                    // MinValue of this submatrix
            node->maxValue = -1;                        // MaxValue of this submatrix
            node->children = nullptr;
            node->plainChildren = nullptr;

            if (this->isLevelInPlain(level + 1)) {
                // Only enter here if the type of k2-raster uses plain values
                // Return the values of submatrix in plain
                node->plainChildren = (int *)malloc(subMatrixSize * subMatrixSize * sizeof(int));

                // Get all values and seek the max and min value
                for (uint r = baseRow; r < baseRow + subMatrixSize; r++) {
                    for (uint c = baseCol; c < baseCol + subMatrixSize; c++) {
                        if (r < this->realSizeX && c < this->realSizeY) {
                            node->plainChildren[pos] = data[r][c];
                            if (node->plainChildren[pos] > node->maxValue) {
                                node->maxValue = node->plainChildren[pos];
                            }
                            if (node->plainChildren[pos] < node->minValue) {
                                node->minValue = node->plainChildren[pos];
                            }
                            if ( node->plainChildren[pos] != 0 &&  node->plainChildren[pos] < this->minValueNoZero) {
                                this->minValueNoZero =  node->plainChildren[pos];
                            }
                        } else {
                            node->plainChildren[pos] = -1;
                        }
                        pos++;
                    }
                }

                // Empty or uniform matrix, free all children
                if (node->maxValue == -1 || node->minValue == node->maxValue) {
                    free(node->plainChildren);
                    node->plainChildren = nullptr;
                } else {
                    for (uint c = 0; c < pos; c++) {
                        if (node->plainChildren[c] == -1){
                            plainValues->push_back(0);
                        } else {
                            plainValues->push_back(node->maxValue - node->plainChildren[c]);
                        }

                    }
                    free(node->plainChildren);
                    node->plainChildren = nullptr;
                }
            } else {
                // Follow with the partition
                uint size = this->divKLevels[level + 1];    //Size of the following submatrices (next level)
                uint newK = this->getK(level);              // K for next level
                uint numberOfChildren = newK * newK;
                node->children = (NodeMatrix**)malloc(numberOfChildren * sizeof(NodeMatrix*));

                // Check children and get min and max values of each one.
                for (uint i = 0; i < newK; i++){
                    for (uint j = 0; j < newK; j++){
                        thisBaseRow = baseRow + (size * i);
                        thisBaseCol = baseCol + (size * j);
                        node->children[pos] = nullptr;
                        if (thisBaseRow < this->realSizeX && thisBaseCol < this->realSizeY) {
                            // Recursive search
                            node->children[pos] =  checkSubmatrix(thisBaseRow, thisBaseCol, size, level + 1, data,totalNodes,
                                                                  maxValues, minValues, plainValues, tmpBitmap, levelPosition);
                            if (node->children[pos]->minValue < node->minValue) {
                                node->minValue = node->children[pos]->minValue;
                            }
                            if (node->children[pos]->maxValue > node->maxValue) {
                                node->maxValue = node->children[pos]->maxValue;
                            }
                        }
                        pos++;
                    }
                }

                // Check if matrix is empty or uniform
                if (node->maxValue != -1 && node->minValue != node->maxValue) {
                    // --------------------------------------------------------------------- //
                    // Apply delta to min and max values of children.
                    // All values ​​are always equal to or greater than 0
                    // newMaxValue = parentMaxValue - maxValue (parentMaxValue >= maxValue)
                    // newMinValue = minValue - parentMinValue (parentMinValue <= minValue)
                    // --------------------------------------------------------------------- //
                    for (uint c = 0; c < (newK * newK); c++){
                        if (node->children[c] != nullptr){
                            //node->children[c]->maxValue = node->maxValue - node->children[c]->maxValue;
                            //node->children[c]->minValue -= node->minValue;
                            maxValues[level]->push_back(node->maxValue - node->children[c]->maxValue);
                            if (node->children[c]->maxValue != node->children[c]->minValue) {
                                minValues[level]->push_back(node->children[c]->minValue - node->minValue);
                                bit_set(tmpBitmap[level], levelPosition[level]);
                            } else {
                                bitclean(tmpBitmap[level], levelPosition[level]);
                            }
                            delete(node->children[c]);
                        } else {
                            maxValues[level]->push_back(0);
                            bitclean(tmpBitmap[level], levelPosition[level]);
                        }
                        levelPosition[level]++;
                        totalNodes++;
                    }
                } else {
                    // All children with the same value (we do nothing)
                    for (uint c = 0; c < (newK * newK); c++){
                        if (node->children[c] != nullptr){
                            delete(node->children[c]);
                        }
                    }
                }
                free(node->children);
                node->children = nullptr;
            }
            return node;
        }
    }

    NodeMatrix* K2Raster::checkSubmatrixOperation(uint baseRow, uint baseCol, uint currentLevel, uint &totalNodes,
                                                  std::vector<uint> **maxValues, std::vector<uint> **minValues,
                                                  std::vector<uint> *plainValues, uint **tmpBitmap, uint *levelPosition,
                                                  BasicSubMatrixWithChildren **ptrK2raster1, uint numOfNodes1,
                                                  BasicSubMatrixWithChildren **ptrK2raster2, uint numOfNodes2,
                                                  K2Raster *k2Raster1, K2Raster *k2Raster2, OperationRaster operation) {
        // Check if the current node is an uniform node or leaf node
        bool uniformRaster1 = false;
        BasicSubMatrix* subMatrix1;
        if (numOfNodes1 == 1) {
            subMatrix1 = ptrK2raster1[0];
            uniformRaster1 = subMatrix1->minValue == subMatrix1->maxValue /*||
                    subMatrix1->level == (k2Raster1->nlevels-1)*/;
        }
        bool uniformRaster2 = false;
        BasicSubMatrix* subMatrix2;
        if (numOfNodes2 == 1) {
            subMatrix2 = ptrK2raster2[0];
            uniformRaster2 = subMatrix2->minValue == subMatrix2->maxValue /*||
                    subMatrix2->level == (k2Raster2->nlevels - 1)*/;
        }

        // ------------------------------------------------------------------- //
        // If both nodes are uniforms, we only need their values
        // ------------------------------------------------------------------- //
        if (uniformRaster1 && uniformRaster2) {
            NodeMatrix *childNode = new NodeMatrix();
            switch (operation) {
                case OPERATION_SUM:
                    childNode->maxValue = subMatrix1->maxValue + subMatrix2->maxValue;
                    break;
                case OPERATION_SUBT:
                    childNode->maxValue = this->displacement + (subMatrix1->maxValue - subMatrix2->maxValue);
                    break;
                case OPERATION_MULT:
                    childNode->maxValue = subMatrix1->maxValue * subMatrix2->maxValue;
                    break;
                default:
                    // No valid operation. The process never reaches this point
                    exit(-1);
            }

            childNode->minValue = childNode->maxValue;

            // Store the minimum value of the rater which is different to zero
            if (childNode->minValue != 0 && childNode->minValue < this->minValueNoZero) {
                this->minValueNoZero = childNode->minValue;
            }
            return childNode;
        } else {
            // Follow both trees
            uint pos = 0;
            uint thisBaseRow, thisBaseCol;
            NodeMatrix *node = new NodeMatrix();

            if (this->isLevelInPlain(currentLevel + 1)) {
                exit(-1); //TODO change this
                return nullptr;
            } else {
                // Follow with the partition
                uint size = this->divKLevels[currentLevel + 1];    //Size of the following submatrices (next level)
                uint newK = this->getK(currentLevel);              // K for next level
                uint numberOfChildren = newK * newK;
                node->children = (NodeMatrix **) malloc(numberOfChildren * sizeof(NodeMatrix *));

                BasicSubMatrixWithChildren **ptrChildren1;
                uint numOfNodesContained1=1;
                BasicSubMatrixWithChildren **ptrChildren2;
                uint numOfNodesContained2=1;

                // Check children and get min and max values of each one.
                for (uint i = 0; i < newK; i++) {
                    for (uint j = 0; j < newK; j++) {
                        thisBaseRow = baseRow + (size * i);
                        thisBaseCol = baseCol + (size * j);
                        node->children[pos] = nullptr;

                        // Recursive search
                        if (thisBaseRow < this->realSizeX && thisBaseCol < this->realSizeY) {
                            ptrChildren1 = this->getNodesContained(i, j, newK, size, ptrK2raster1, numOfNodes1,
                                                    k2Raster1, numOfNodesContained1);
                            ptrChildren2 = this->getNodesContained(i, j, newK, size, ptrK2raster2, numOfNodes2,
                                                    k2Raster2, numOfNodesContained2);
//                            ptrChildren1 = !uniformRaster1 ? this->getNodesContained(i, j, newK, size, ptrK2raster1, numOfNodes1,
//                                                                       k2Raster1, numOfNodesContained1) : ptrK2raster1;
//
//                            ptrChildren2 = !uniformRaster2 ? this->getNodesContained(i, j, newK, size, ptrK2raster2, numOfNodes2,
//                                                                                     k2Raster2, numOfNodesContained2) : ptrK2raster2;

                            // Create a new node
                            node->children[pos] = checkSubmatrixOperation(thisBaseRow, thisBaseCol, currentLevel + 1, totalNodes,
                                                                          maxValues, minValues, plainValues, tmpBitmap, levelPosition,
                                                                          ptrChildren1, numOfNodesContained1, ptrChildren2, numOfNodesContained2,
                                                                          k2Raster1, k2Raster2, operation);

                            // Check the new node
                            if (node->children[pos]->minValue < node->minValue) {
                                node->minValue = node->children[pos]->minValue;
                            }
                            if (node->children[pos]->maxValue > node->maxValue) {
                                node->maxValue = node->children[pos]->maxValue;
                            }

                            // Free memory
                            {
                                for (uint c = 0; c < numOfNodesContained1; c++) {
                                    removeChildren(ptrChildren1[c]);
                                }
                                delete[](ptrChildren1);
                                for (uint c = 0; c < numOfNodesContained2; c++) {
                                    removeChildren(ptrChildren2[c]);
                                }
                                delete[](ptrChildren2);
                            } // END BLOCK Free memory
                        }
                        pos++;
                    } // END FOR  j < newk
                } // END FOR  i < newk

                // Check if matrix is empty or uniform
                if (node->maxValue != -1 && node->minValue != node->maxValue) {
                    // ------------------------------------------------------------------------ //
                    // Apply delta to min and max values of the children.
                    // All values ​​are always equal to or greater than 0
                    // newMaxValue = parentMaxValue - maxValue (parentMaxValue >= maxValue)
                    // newMinValue = minValue - parentMinValue (parentMinValue <= minValue)
                    // ------------------------------------------------------------------------ //
                    for (uint c = 0; c < numberOfChildren; c++){
                        if (node->children[c] != nullptr){
                            //node->children[c]->maxValue = node->maxValue - node->children[c]->maxValue;
                            //node->children[c]->minValue -= node->minValue;
                            maxValues[currentLevel]->push_back(node->maxValue - node->children[c]->maxValue);
                            if (currentLevel + 1 != this->nlevels - 1) {
                                if (node->children[c]->maxValue != node->children[c]->minValue) {
                                    minValues[currentLevel]->push_back(node->children[c]->minValue - node->minValue);
                                    bit_set(tmpBitmap[currentLevel], levelPosition[currentLevel]);
                                } else {
                                    bitclean(tmpBitmap[currentLevel], levelPosition[currentLevel]);
                                }
                            }
                            delete(node->children[c]);
                        } else {
                            maxValues[currentLevel]->push_back(0);
                            if (currentLevel + 1 != this->nlevels - 1) {
                                bitclean(tmpBitmap[currentLevel], levelPosition[currentLevel]);
                            }
                        }
                        levelPosition[currentLevel]++;
                        totalNodes++;
                    }
                } else {
                    // All children with the same value (we do nothing)
                    for (uint c = 0; c < numberOfChildren; c++){
                        if (node->children[c] != nullptr){
                            delete(node->children[c]);
                        }
                    }
                }
                free(node->children);
                node->children = nullptr;
            } // END IF isLevelInPlain(currentLevel + 1)
            return node;
        } // END IF uniformRaster1 && uniformRaster2
    }

    // Equal K's
    NodeMatrix* K2Raster::checkSubmatrixOperation(uint baseRow, uint baseCol, uint currentLevel, uint &totalNodes,
                                                  std::vector<uint> **maxValues, std::vector<uint> **minValues,
                                                  std::vector<uint> *plainValues, uint **tmpBitmap, uint *levelPosition,
                                                  BasicSubMatrix *ptrK2raster1, BasicSubMatrix *ptrK2raster2,
                                                  K2Raster *k2Raster1, K2Raster *k2Raster2, uint *countOnes1, uint *countOnes2,
                                                  OperationRaster operation) {

        // Check if the current node is an uniform node or leaf node
        bool uniformRaster1 = ptrK2raster1->maxValue == ptrK2raster1->minValue;
        bool uniformRaster2 = ptrK2raster2->maxValue == ptrK2raster2->minValue;

        // ------------------------------------------------------------------- //
        // If both nodes are uniforms, we only need their values
        // ------------------------------------------------------------------- //
        if (uniformRaster1 && uniformRaster2) {
            NodeMatrix *childNode = new NodeMatrix();
            switch (operation) {
                case OPERATION_SUM:
                    childNode->maxValue = ptrK2raster1->maxValue + ptrK2raster2->maxValue;
                    break;
                case OPERATION_SUBT:
                    childNode->maxValue = this->displacement + (ptrK2raster1->maxValue - ptrK2raster2->maxValue);
                    break;
                case OPERATION_MULT:
                    childNode->maxValue = ptrK2raster1->maxValue * ptrK2raster2->maxValue;
                    break;
                default:
                    // No valid operation. The process never reaches this point
                    exit(-1);
            }

            childNode->minValue = childNode->maxValue;

            // Store the minimum value of the rater which is different to zero
            if (childNode->minValue != 0 && childNode->minValue < this->minValueNoZero) {
                this->minValueNoZero = childNode->minValue;
            }
            return childNode;
        } else {
            // Follow both trees
            uint pos = 0;
            uint thisBaseRow = baseRow, thisBaseCol = baseCol;
            NodeMatrix *node = new NodeMatrix();

//            if (this->isLevelInPlain(currentLevel + 1)) {
//                exit(-1); //TODO change this
//                return nullptr;
//            } else {
                // Follow with the partition
                uint size = this->divKLevels[currentLevel + 1];    //Size of the following submatrices (next level)
                uint newK = this->getK(currentLevel);              // K for next level
                uint numberOfChildren = newK * newK;
                node->children = (NodeMatrix**)malloc(numberOfChildren * sizeof(NodeMatrix*));

                // ------------------------------------------------------------------- //
                // Child of k2-Raster 1
                // ------------------------------------------------------------------- //
                BasicSubMatrix *ptrChild1 = uniformRaster1 ? ptrK2raster1 : this->getFirstChild(k2Raster1, ptrK2raster1, uniformRaster1,
                                                                                                countOnes1);

                // ------------------------------------------------------------------- //
                // Child of k2-Raster 2
                // ------------------------------------------------------------------- //
//                BasicSubMatrix *ptrChild2 = this->getFirstChild(k2Raster2, ptrK2raster2, uniformRaster2, countOnes2);
                BasicSubMatrix *ptrChild2 = uniformRaster2 ? ptrK2raster2 : this->getFirstChild(k2Raster2, ptrK2raster2, uniformRaster2,
                                                                                                countOnes2);

                // Check children and get min and max values of each one.
                for (uint i = 0; i < newK; i++) {
                    for (uint j = 0; j < newK; j++) {
                        node->children[pos] = nullptr;

                        if (i != 0 || j != 0) { // TODO update this when thisBaseRow < this->realSizeX && thisBaseCol < this->realSizeY
                            // ------------------------------------------------------------------- //
                            // Update children values. Values of the first child are  already updated
                            // ------------------------------------------------------------------- //
                            if (!uniformRaster1) {
                                this->getNextChild(k2Raster1, ptrK2raster1, ptrChild1, countOnes1);
                            }

                            if (!uniformRaster2) {
                                this->getNextChild(k2Raster2, ptrK2raster2, ptrChild2, countOnes2);
                            }
                        }

                        if (thisBaseRow < this->realSizeX && thisBaseCol < this->realSizeY) {
                            // Recursive search
                            node->children[pos] = checkSubmatrixOperation(thisBaseRow, thisBaseCol, currentLevel + 1,
                                                                          totalNodes,
                                                                          maxValues, minValues, plainValues, tmpBitmap,
                                                                          levelPosition,
                                                                          ptrChild1, ptrChild2, k2Raster1, k2Raster2,
                                                                          countOnes1, countOnes2,
                                                                          operation);
                            if (node->children[pos]->minValue < node->minValue) {
                                node->minValue = node->children[pos]->minValue;
                            }
                            if (node->children[pos]->maxValue > node->maxValue) {
                                node->maxValue = node->children[pos]->maxValue;
                            }
                        }
                        pos++;
                        thisBaseCol +=size;
                    } // END FOR  j < newk
                    thisBaseCol = baseCol;
                    thisBaseRow +=size;
                } // END FOR  i < newk

                // Free ptrChild
                if (!uniformRaster1){
                    delete(ptrChild1);
                }
                if (!uniformRaster2){
                    delete(ptrChild2);
                };

                // Check if matrix is empty or uniform
                if (node->maxValue != -1 && node->minValue != node->maxValue) {
                    /**
                     * Apply delta to min and max values of children.
                     * All values ​​are always equal to or greater than 0
                     * newMaxValue = parentMaxValue - maxValue (parentMaxValue >= maxValue)
                     * newMinValue = minValue - parentMinValue (parentMinValue <= minValue)
                     */
                    for (uint c = 0; c < numberOfChildren; c++){
                        if (node->children[c] != nullptr){
                            //node->children[c]->maxValue = node->maxValue - node->children[c]->maxValue;
                            //node->children[c]->minValue -= node->minValue;
                            maxValues[currentLevel]->push_back(node->maxValue - node->children[c]->maxValue);
                            if (currentLevel + 1 != this->nlevels - 1) {
                                if (node->children[c]->maxValue != node->children[c]->minValue) {
                                    minValues[currentLevel]->push_back(node->children[c]->minValue - node->minValue);
                                    bit_set(tmpBitmap[currentLevel], levelPosition[currentLevel]);
                                } else {
                                    bitclean(tmpBitmap[currentLevel], levelPosition[currentLevel]);
                                }
                            }
                            delete(node->children[c]);
                        } else {
                            maxValues[currentLevel]->push_back(0);
                            bitclean(tmpBitmap[currentLevel], levelPosition[currentLevel]);
                        }
                        levelPosition[currentLevel]++;
                        totalNodes++;
                    }
                } else {
                    // All children with the same value (we do nothing)
                    for (uint c = 0; c < numberOfChildren; c++){
                        if (node->children[c] != nullptr){
                            delete(node->children[c]);
                        }
                    }
                }
                free(node->children);
                node->children = nullptr;
//            } // ENF IF (this->isLevelInPlain(level + 1))
//            if (node->minValue == 0) printf("Node at level %u has a 0 value", currentLevel);
            return node;
        }
    }

    BasicSubMatrix* K2Raster::getFirstChild(K2Raster *k2Raster, BasicSubMatrix * ptrK2raster, bool uniformNode, uint *countOnes){
        BasicSubMatrix *ptrChild = new BasicSubMatrix; 
        this->initSubMatrix(ptrChild, ptrK2raster->minValue, ptrK2raster->maxValue);
        ptrChild->level = ptrK2raster->level + 1;

        // Calculate the position of the first child and update its values
        if (!uniformNode) {
            ptrChild->bitmapPosition = ptrK2raster->bitmapChildren;

            // Update max value
            if (ptrChild->bitmapPosition == k2Raster->preMaxValues[ptrChild->level - 1]) {
                ptrChild->maxValue -= k2Raster->listMaxValues[ptrChild->level - 1]->access(0);
            } else {
                ptrChild->maxValue -= k2Raster->listMaxValues[ptrChild->level - 1]->next();
            }

            // Update min value
            if (ptrChild->level != (k2Raster->nlevels - 1) &&
                k2Raster->bitmapK2Tree->access(ptrChild->bitmapPosition)) {
                ulong ones = countOnes[ptrChild->level - 1] + k2Raster->preMinValues[ptrChild->level - 1];
                ptrChild->minValue = ptrChild->maxValue + 1;
                if (ptrChild->level < k2Raster->levelK1) {
                    ptrChild->bitmapChildren = (ones) * k2Raster->K1 * k2Raster->K1 + 1;
                } else {
                    ptrChild->bitmapChildren = (k2Raster->countOnesk1 * k2Raster->K1 * k2Raster->K1) +
                                               ((ones - k2Raster->countOnesk1) * k2Raster->K2 * k2Raster->K2 + 1);
                }
                countOnes[ptrChild->level - 1]++;
            } else {
                ptrChild->minValue = ptrChild->maxValue;
            }
        }
        return ptrChild;
    }

    void K2Raster::getNextChild(K2Raster *k2Raster, BasicSubMatrix *ptrK2raster, BasicSubMatrix *ptrChild, uint *countOnes){
        // update position
        ptrChild->bitmapPosition++;

        // Update max value
        ptrChild->maxValue = ptrK2raster->maxValue - k2Raster->listMaxValues[ptrChild->level - 1]->next();

        // Update min value
        if (ptrChild->level != (k2Raster->nlevels-1) && k2Raster->bitmapK2Tree->access(ptrChild->bitmapPosition)) {
            ptrChild->minValue = ptrChild->maxValue+1;
            ulong ones = countOnes[ptrChild->level-1] + k2Raster->preMinValues[ptrChild->level - 1];
            if (ptrChild->level < k2Raster->levelK1) {
                ptrChild->bitmapChildren = (ones) * k2Raster->K1 * k2Raster->K1 + 1;
            } else {
                ptrChild->bitmapChildren = (k2Raster->countOnesk1 * k2Raster->K1 * k2Raster->K1) +
                                           ((ones - k2Raster->countOnesk1) * k2Raster->K2 * k2Raster->K2 + 1);
            }
            countOnes[ptrChild->level-1]++;
        } else {
            ptrChild->minValue = ptrChild->maxValue;
        }
    }

    BasicSubMatrixWithChildren** K2Raster::getNodesContained(uint i, uint j, uint k, uint size,
                                                                           BasicSubMatrixWithChildren **ptrK2raster, uint numOfNodesPtr,
                                                                           K2Raster *raster, uint &numOfNodesContained) {

        BasicSubMatrixWithChildren **ptrChildren;
        numOfNodesContained = 1;
        uint sizeSubmatrix = raster->divKLevels[ptrK2raster[0]->level];
        BasicSubMatrixWithChildren **currentNodes = ptrK2raster;
        uint numOfChildren = numOfNodesPtr;

        // It is an uniform node, return the same node
        if (numOfNodesPtr == 1 && ptrK2raster[0]->minValue == ptrK2raster[0]->maxValue) {
            BasicSubMatrixWithChildren *node = ptrK2raster[0];
            ptrChildren = new BasicSubMatrixWithChildren*[1];
            ptrChildren[0] = node;
            return ptrChildren;
        } // END Uniform node

        while (1) {
            // CASE 1
            if (size == sizeSubmatrix) {
                // The result node and the node of tree have same size
                BasicSubMatrixWithChildren *node = currentNodes[i * k + j];
                ptrChildren = new BasicSubMatrixWithChildren*[1];
                ptrChildren[0] = node;
                break;
            } // END IF size == sizeSubmatrix1

            // CASE 2
            if (size > sizeSubmatrix) {
                // The node result is bigger than nodes of the tree.
                // Get several nodes.
                uint numOfNodes = size / sizeSubmatrix;
                ptrChildren = new BasicSubMatrixWithChildren*[numOfNodes * numOfNodes];
                numOfNodesContained = numOfNodes * numOfNodes;
                uint c = 0;
                for (uint x = i * numOfNodes; x <  i * numOfNodes + numOfNodes; x++) {
                    for (uint y = j * numOfNodes; y < j * numOfNodes + numOfNodes; y++) {
                        BasicSubMatrixWithChildren *node = currentNodes[x * numOfNodes * k + y];
                        ptrChildren[c] = node;
                        c++;
                    }
                }
                break;
            } // END IF size > sizeSubmatrix1

            // CASE 3
            if (size < sizeSubmatrix) {
                // Get the node where it is contained.
                uint numOfNodes = sizeSubmatrix / size;
                k = sqrt(numOfChildren);
                uint x = i / numOfNodes;
                uint y = j / numOfNodes;
                BasicSubMatrixWithChildren *node = currentNodes[x * k + y];

                if (node->maxValue == node->minValue) {
                    ptrChildren = new BasicSubMatrixWithChildren*[1];
                    ptrChildren[0] = node;
                    break;
                } else {
                    currentNodes = this->getChildrenNode(node, raster, numOfChildren);
                    sizeSubmatrix = raster->divKLevels[node->level + 1];
                    i = i % numOfNodes;
                    j = j % numOfNodes;
                    k = size < sizeSubmatrix ? raster->getK(node->level) : numOfNodes;
                }
            }
        }
        return ptrChildren;
    }

    BasicSubMatrixWithChildren ** K2Raster::getChildrenNode(BasicSubMatrixWithChildren *node, K2Raster *raster, uint &numOfChildren) {
        if (node->children == nullptr) {
            // Get position of the children of the currentNode
            ulong countOnes = node->bitmapPosition == 0 ? 0 : raster->bitmapK2Tree->rank1(
                    node->bitmapPosition - 1);
            ulong childrenPosition = 0;
            if (node->level < raster->levelK1) {
                childrenPosition = (countOnes) * raster->K1 * raster->K1 + 1;
            } else {
                childrenPosition = (raster->countOnesk1 * raster->K1 * raster->K1) +
                                   ((countOnes - raster->countOnesk1) * raster->K2 *
                                            raster->K2 + 1);
            } // END IF get childrenPosition

            // Create children
            numOfChildren = raster->getK(node->level) * raster->getK(node->level);
            node->children = new BasicSubMatrixWithChildren*[numOfChildren];
            node->numOfChildren = numOfChildren;
            bool accessMinDac = false;
            for (uint c = 0; c < numOfChildren; c++) {
                BasicSubMatrixWithChildren *child = new BasicSubMatrixWithChildren;
                child->bitmapPosition = childrenPosition + c;
                child->level = node->level + 1;
                child->children = nullptr;
                child->numOfChildren = 0;

                // Get max value
                if (c == 0) {
                    child->maxValue = node->maxValue -
                            raster->listMaxValues[child->level - 1]->access(
                                              child->bitmapPosition -
                                                      raster->preMaxValues[child->level - 1]);
                } else {
                    child->maxValue = node->maxValue - raster->listMaxValues[child->level -1]->next();
                } // END IF get max value

                // Get min value
                if (child->level != (raster->nlevels-1) // If it is not a leaf
                    && raster->bitmapK2Tree->access(child->bitmapPosition)) {
                    child->minValue = child->maxValue+1; // All values of the node are equal
//                    if (!accessMinDac) {
//                        countOnes = raster->bitmapK2Tree->rank1(child->bitmapPosition - 1);
//                        child->minValue = node->minValue +
//                                raster->listMinValues[child->level - 1]->access(
//                                                  countOnes - raster->preMinValues[child->level - 1]);
//                        accessMinDac = true;
//                    } else {
//                        child->minValue = node->minValue + raster->listMinValues[child->level -1]->next();
//                    }
                } else {
                    child->minValue = child->maxValue; // All values of the node are equal
                } // END IF get min value

                node->children[c] = child;

            } // END FOR numOfChildren
        } // END IF  node->children == nullptr
        numOfChildren = node->numOfChildren;
        return node->children;
    }

    uint K2Raster::conceptualToCompact(NodeMatrix *root, uint totalNodes, std::vector<uint> **maxValues, std::vector<uint> **minValues,
                                       uint **tmpBitmap, uint *levelPositions) {
        // ------------------------------------------------------------------- //
        // Prepare all structures of k2-raster to store the conceptual tree
        // ------------------------------------------------------------------- //

        // k2-tree
        uint *tmpK2Tree = (uint * )calloc(uint_len(1, totalNodes), sizeof(uint));

        // Prepare DACs (first level and levels from levelK1 + levelK2 to nlevels are stored in plain form)
        this->listMinValues = (cds_static::DirectAccess **)malloc(this->numOfMinDACs * sizeof(cds_static::DirectAccess *));
        this->preMinValues = (ulong *) malloc(this->numOfMinDACs * sizeof(ulong));
        this->listMaxValues = (cds_static::DirectAccess **)malloc(this->numOfMaxDACs * sizeof(cds_static::DirectAccess *));
        this->preMaxValues = (ulong *) malloc(this->numOfMaxDACs * sizeof(ulong));

        // ------------------------------------------------------------------- //
        // The algorithm traverses the conceptual tree
        // and adds the values to their respective structures
        // ------------------------------------------------------------------- //

        // Start with the root node
        totalNodes = 0;
        this->initialMaxValue = root->maxValue;
        this->initialMinValue = root->minValue;
        if (root->maxValue != root->minValue) {
            bit_set(tmpK2Tree, totalNodes++);
        } else {
            bitclean(tmpK2Tree, totalNodes++);
        }
        delete(root);
        this->countOnesk1 = 1;       // Add root node

        // Run the rest of the levels
        uint l = 0;
        for (l = 0; l < this->numOfMaxDACs; l++) {

            if (levelPositions[l] == 0) {
                // Empty level. No more nodes in the tree
                // "l" is the last level of the tree
                break;
            }

            // If the current level is not the last level
            // set the corresponding bit in the k2-tree
            if (l + 1 != (this->nlevels - 1)) {
                // Add bit to k2-tree
                for (uint b = 0; b < levelPositions[l]; b++) {
                    if (bitget(tmpBitmap[l], b)) {
                        bit_set(tmpK2Tree, totalNodes++);
                    } else {
                        bitclean(tmpK2Tree, totalNodes++);
                    }
                }
                free(tmpBitmap[l]);
            } // END IF not lastLevel

            // Create a DAC of maximum values for the current level "l"
            this->listMaxValues[l] = new cds_static::DAC(maxValues[l]->data(), maxValues[l]->size(), 3);
            this->preMaxValues[l] =  l == 0 ? 1 : this->preMaxValues[l - 1] + this->listMaxValues[l - 1]->getRealLength();
            delete(maxValues[l]);

            // Create a DAC of minimum values for the current level "l", if it is necessary
            if (l < this->numOfMinDACs) {
                this->listMinValues[l] = new cds_static::DAC(minValues[l]->data(), minValues[l]->size(), 3);
                this->preMinValues[l] = l == 0 ? 1 : this->preMinValues[l - 1] + this->listMinValues[l - 1]->getRealLength();
                if (l + 1 < this->levelK1) {
                    // Count the number of "1s" in the k2-tree for partitions of "k1"
                    this->countOnesk1 += this->listMinValues[l]->getRealLength();
                }
                delete(minValues[l]);
            } // END IF minimum values
        } // END FOR run the rest of the levels

        // The current level of the conceptual tree is "l".
        // Free the memory for the rest of the levels and
        // update numOfMaxDACs and numOfMinDACs with the new level
        if (l != this->numOfMaxDACs) {
            for (uint i = l; i < this->numOfMaxDACs; i++) {
                delete(maxValues[i]);
                if ( i < this->numOfMinDACs) {
                    delete(minValues[i]);
                }
                if (i + 1 != this->nlevels - 1) {
                    free(tmpBitmap[i]);
                }
            }

            this->numOfMaxDACs = l;
            this->numOfMinDACs = l;
        }

        // ------------------------------------------------------------------- //
        // Free the memory for the temporal estructures
        // ------------------------------------------------------------------- //
        free(levelPositions);
        free(tmpBitmap);
        free(maxValues);
        free(minValues);

        // ------------------------------------------------------------------- //
        // Create k2-tree with factor 20 (overhead 5%).
        // It is necessary  a rank function
        // ------------------------------------------------------------------- //
        cds_static::BitSequenceBuilder *builder = new cds_static::BitSequenceBuilderRG(20);
        this->bitmapK2Tree = builder->build(tmpK2Tree, totalNodes);
        free(tmpK2Tree);
        delete(builder);

        if (!this->isLevelInPlain(l + 1)) {
            // Not need tmpPlanValues vector anymore
            delete (this->tmpPlainValues);
        }
        return l;
    }

    uint K2Raster::getK(uint level) const {
        if (level < this->levelK1){
            return this->K1;
        }
        return this->K2;
    }

    bool K2Raster::isLevelInPlain(uint level) const {
        return false;
    }


    QNodeCounters* K2Raster::createRoot() const {
        QNodeCounters *node = new QNodeCounters(1, 0, 0, 0, this->initialMinValue, this->initialMaxValue);
        node->preValuesMax = 1;
        node->preValuesMin = 1;
        node->offsetPrev = this->getK(0) * this->getK(0);
        return node;
    }

    QNodeCounters * K2Raster::createChild(QNodeCounters *parentNode, uint offsetValues,
                                      uint thisbaserow, uint thisbasecol, uint minValue, uint maxValue) const {
        uint vk = this->getK(parentNode->level + 1);
        vk *= vk;
        QNodeCounters *newNode = new QNodeCounters(parentNode->offsetPrev + (offsetValues * vk) + 1, parentNode->level + 1,
                                   thisbaserow, thisbasecol, minValue, maxValue);

        // *********************************************************************************************************
        //  Update counters previous numbers.
        //  Those counters are necessary to speed up the query and because
        //  the value of "k" can be different in each level, else it is more difficult calculate the position
        //  of the next submatrix
        // *********************************************************************************************************
        newNode->preValuesMax = parentNode->preValuesMax + this->listMaxValues[parentNode->level]->getRealLength();     // Number of maximum values in previous levels
        newNode->preValuesMin = parentNode->preValuesMin + this->listMinValues[parentNode->level]->getRealLength();     // Number of minimum values in previous levels
        newNode->offsetPrev = parentNode->offsetPrev + this->listMinValues[parentNode->level]->getRealLength() * vk;    // Position in k2-tree of next level
        return newNode;
    }

    int K2Raster::nextChild(QNodeCounters *parent,int &posI, int &posJ,
                         int &thisbaserow, int &thismaxrow, int &thisbasecol, int &thismaxcol,
                         Query query, uint newK, uint kLevel) const {
        return this->nextChild(parent, posI, posJ, thisbaserow, thismaxrow, thisbasecol, thismaxcol,
                               query.xini, query.xend, query.yini, query.yend, newK, kLevel);
    }

    int K2Raster::nextChild(QNodeCounters *parent, int &posI, int &posJ, int &thisbaserow, int &thismaxrow,
                            int &thisbasecol, int &thismaxcol, uint xini, uint xend, uint yini, uint yend, uint newK,
                            uint kLevel) const {

        for (uint i = posI; i < newK; i++) {
            thisbaserow = parent->baseRow + i * kLevel;
            thismaxrow = thisbaserow + kLevel - 1;

            if (thismaxrow < xini) {
                posJ = 0;
                continue;
            }
            if (thisbaserow > xend || thisbaserow >= this->realSizeX) {
                //Positions out of range
                posI = 0;
                break;
            }

            for (uint j = posJ; j < newK; j++) {
                thisbasecol = parent->baseCol + j * kLevel;
                thismaxcol = thisbasecol + kLevel - 1;

                if (thismaxcol < yini) {
                    continue;
                }
                if (thisbasecol > yend || thisbasecol >= this->realSizeY) {
                    //Positions out of range
                    posJ = 0;
                    break;
                }

                posI = i;
                posJ = j + 1;
                return parent->position + i * newK + j; // Tree position of the node
            }
            posJ = 0;
        }
        return -1;
    }

    void K2Raster::initSubMatrix(BasicSubMatrix *subMatrix, uint valini, uint valend) const{
        subMatrix->bitmapPosition = 0;
        subMatrix->bitmapChildren = 1;
        subMatrix->level = 0;
        subMatrix->minValue = valini;
        subMatrix->maxValue = valend;
        return;
    }

    void K2Raster::initMinSubMatrix(MinBasicSubMatrix *subMatrix, uint valini, uint valend) const{
        this->initSubMatrix(subMatrix, valini, valend);
        subMatrix->bitmapChildren = 1;

        // Set submaxtrix color
        if (subMatrix->minValue > valend || subMatrix->maxValue < valini) {
            // It is an invalid node, value out of range
            subMatrix->color = SubMatrixColor::COLOR_SUBMATRIX_WHITE;
        } else {
            if (subMatrix->minValue >= valini && subMatrix->maxValue <= valend) {
                // All cells lie within the range of value
                subMatrix->color = SubMatrixColor::COLOR_SUBMATRIX_BLACK;
            } else {
                // We can not determinate if all cells lie within the range of values
                subMatrix->color = SubMatrixColor::COLOR_SUBMATRIX_GREY;
            }
        } // END IF SET COLOR

        return;
    }

    MinNodeTreeNode *K2Raster::createNexMinNodeTreeNode(uint bitmapPosition, uint level, MinNodeTreeNode *node,
                                                            uint valini, uint valend, bool addChildren) const {
        MinNodeTreeNode *nextNode;
        nextNode = (MinNodeTreeNode *) malloc(sizeof(MinNodeTreeNode));
        nextNode->bitmapPosition = bitmapPosition;
        nextNode->level = level;
        nextNode->queue = nullptr;
//        nextNode->used = false;
        if (addChildren) {
            nextNode->numberOfChildren = this->getK(level) * this->getK(level);
            nextNode->children = nullptr;
//            nextNode->children = (MinNodeTreeNode **) malloc(
//                    nextNode->numberOfChildren * sizeof(MinNodeTreeNode *));
//            for (uint c = 0; c < nextNode->numberOfChildren; c++) {
//                nextNode->children[c] = nullptr;
//            }
        }

        // Create queue
//        nextNode->queue = new std::queue<int64_t>();

        // Check whether it is uniform or not
        bool uniform = (level == this->nlevels - 1) ||
                       (this->bitmapK2Tree->access(nextNode->bitmapPosition) == 0);  // Last level or an uniform node

        // Calculate max value
        nextNode->maxValue = node->maxValue - this->listMaxValues[level - 1]->access(
                nextNode->bitmapPosition - this->preMaxValues[level - 1]);

        if (uniform) {
            // Uniform matrix, it has no children (minValue == maxValue)
            nextNode->minValue = nextNode->maxValue;
        } else {
            uint countOnes = this->bitmapK2Tree->rank1(nextNode->bitmapPosition - 1);

            // Calculate min Value
            nextNode->minValue = node->minValue + this->listMinValues[level - 1]->access(
                    countOnes - this->preMinValues[level - 1]);

            // Calculate children position at bitmap
            if (level < this->levelK1) {
                nextNode->bitmapChildren = (countOnes) * this->K1 * this->K1 + 1;
            } else {
                nextNode->bitmapChildren = (this->countOnesk1 * this->K1 * this->K1) +
                                           ((countOnes - this->countOnesk1) * this->K2 * this->K2 + 1);
            }
        } // End IF uniform

        // Set submaxtrix color
        if (nextNode->minValue > valend || nextNode->maxValue < valini) {
            // It is an invalid node, value out of range
            nextNode->color = SubMatrixColor::COLOR_SUBMATRIX_WHITE;
        } else {
            if (nextNode->minValue >= valini && nextNode->maxValue <= valend) {
                // nextNode cells lie within the range of value
                nextNode->color = SubMatrixColor::COLOR_SUBMATRIX_BLACK;
            } else {
                // We can not determinate if all cells lie within the range of values
                nextNode->color = SubMatrixColor::COLOR_SUBMATRIX_GREY;
            }
        } // END IF set color

        return nextNode;
    }

    void K2Raster::createChildrenNodeMinSubMatrix(MinNodeTreeNode *node, uint k, uint valini, uint valend,
                                                  std::map<ulong, MinNodeTreeNode *> *nextMBRs,
                                                  bool addChildren) const {
        uint bitmapPosition = node->bitmapChildren;
        uint countOnes = this->bitmapK2Tree->rank1(bitmapPosition - 1);
        uint level = node->level + 1;
        MinNodeTreeNode *nextNode;
        bool accessMin = false;

        for (uint x = 0; x < k; x++) {
            for (uint y = 0; y < k; y++) {
                nextNode = (MinNodeTreeNode *) malloc(sizeof(MinNodeTreeNode));
                nextNode->bitmapPosition = bitmapPosition;
                nextNode->level = level;
                nextNode->queue = nullptr;
//                nextNode->used = false;

                if (addChildren) {
                    nextNode->numberOfChildren = this->getK(level) * this->getK(level);
                    nextNode->children = nullptr;
//                    nextNode->children = (MinNodeTreeNode **) malloc(
//                            nextNode->numberOfChildren * sizeof(MinNodeTreeNode *));
//                    for (uint c = 0; c < nextNode->numberOfChildren; c++) {
//                        nextNode->children[c] = nullptr;
//                    }
                }

                // Check whether it is uniform or not
                bool uniform = (nextNode->level == this->nlevels - 1) ||
                               (this->bitmapK2Tree->access(nextNode->bitmapPosition) ==
                                0);  // Last level or an uniform node

                // Calculate max value
                if (x == 0 && y == 0) {
                    nextNode->maxValue = node->maxValue - this->listMaxValues[level - 1]->access(
                            nextNode->bitmapPosition - this->preMaxValues[level - 1]);
                } else {
                    nextNode->maxValue = node->maxValue - this->listMaxValues[level - 1]->next();
                }

                if (uniform) {
                    // Uniform matrix, it has no children (minValue == maxValue)
                    nextNode->minValue = nextNode->maxValue;
                } else {
                    // Calculate min Value
                    if (!accessMin) {
                        // First child
                        nextNode->minValue = node->minValue + this->listMinValues[level - 1]->access(
                                countOnes - this->preMinValues[level - 1]);
                        accessMin = true;
                    } else {
                        nextNode->minValue = node->minValue + this->listMinValues[level - 1]->next();
                    }
                    // Calculate children position at bitmap
                    if (level < this->levelK1) {
                        nextNode->bitmapChildren = (countOnes) * this->K1 * this->K1 + 1; // + 1 -> root node
                    } else {
                        nextNode->bitmapChildren = (this->countOnesk1 * this->K1 * this->K1) +
                                                   ((countOnes - this->countOnesk1) * this->K2 * this->K2 +
                                                    1); // + 1 -> root node
                    }
                    countOnes++;
                } // End IF uniform

                // Set submaxtrix color
                if (nextNode->minValue > valend || nextNode->maxValue < valini) {
                    // It is an invalid node, value out of range
                    nextNode->color = SubMatrixColor::COLOR_SUBMATRIX_WHITE;
                } else {
                    if (nextNode->minValue >= valini && nextNode->maxValue <= valend) {
                        // nextNode cells lie within the range of value
                        nextNode->color = SubMatrixColor::COLOR_SUBMATRIX_BLACK;
                    } else {
                        // We can not determinate if all cells lie within the range of values
                        nextNode->color = SubMatrixColor::COLOR_SUBMATRIX_GREY;
                    }
                } // END IF set color

                // Add new node to map of MBRs
                if (addChildren) {
                    node->children[x * k + y] = nextNode;
                } else {
                    (*nextMBRs)[nextNode->bitmapPosition] = nextNode;
                }
                bitmapPosition++;
            } // END FOR y
        } // END FOR x
        return;
    }

}
