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
#include <k2raster/plain/k2-raster-plain-VByte.h>

namespace k2raster_static {

    //**********************************************************************//
    //********************* CONSTRUCTORS ***********************************//
    //**********************************************************************//

    K2RasterPlainVByte::K2RasterPlainVByte() { };

    K2RasterPlainVByte::K2RasterPlainVByte(uint sizeX, uint sizeY, int **data, uint K1, uint K2, uint levelK1, uint plainLevels) :
            K2RasterPlain(sizeX, sizeY, data, K1, K2, levelK1, plainLevels) {
        if (this->levelK1 + this->levelK2 < this->nlevels ) {
            // Save values of the last levels
           /// this->lastValues = new cds_static::DAC(this->plainValues, this->numOfPlainValues, 3);
            this->vByteEncoder = new libencoders_static::VByte(this->plainValues, this->numOfPlainValues, 3, 1000);
            free(this->plainValues);
            this->plainValues = nullptr;
//            printf("Saved %u plain values from level %u to %u\n", this->numOfPlainValues, this->levelK1 + this->levelK2, this->nlevels);
        }
    }

    K2RasterPlainVByte::~K2RasterPlainVByte() {
        if (this->vByteEncoder != nullptr) {
            delete(this->vByteEncoder);
        }
    }

    //**********************************************************************//
    //********************* CELL FUNCTIONS *********************************//
    //**********************************************************************//

    int K2RasterPlainVByte::getPlainCell(uint father, uint row, uint col, uint divK, uint lastValueMax) const {
        uint arrayPosition = father * divK * divK;
        row = row % divK;
        col = col % divK;
        arrayPosition +=  row * divK + col;
        return lastValueMax - this->vByteEncoder->access(arrayPosition);
    }

    ulong K2RasterPlainVByte::getPlainWindow(uint xini, uint yini, uint xend, uint yend,
                                             uint offsetValues, uint divK, uint maxValue,
                                             uint *positions, ulong totalCells, Query query) const {

        // Add all valid positions to final result
        uint value, arrayPosition;
        arrayPosition = (offsetValues * divK * divK) + (xini % divK) * divK + (yini % divK);
        uint p = this->vByteEncoder->searchNextPosition(0, arrayPosition);
        for (uint x = xini; x <= xend; x++){
            for (uint y= yini; y <= yend; y++){
                p = this->vByteEncoder->decodeNext(p, value);
                value = maxValue - value;
                if (value >= query.valini && value <= query.valend) {
                    positions[totalCells++] = x * this->realSizeY + y;
                }
                arrayPosition++;
            }
            p = this->vByteEncoder->searchNextPosition(p, divK - (1 + yend - yini));
        }
        return totalCells;
    }

    bool K2RasterPlainVByte::checkPlainCells(uint xini, uint yini, uint xend, uint yend, uint offsetValues, uint divK,
                                        uint maxValue, Query query, bool allCells) const {

        // Add all valid positions to final result
        uint value, arrayPosition;
        arrayPosition = (offsetValues * divK * divK) + (xini % divK) * divK + (yini % divK);
        uint p = this->vByteEncoder->searchNextPosition(0, arrayPosition);
        for (uint x = xini; x <= xend; x++){
            for (uint y= yini; y <= yend; y++){
                p = this->vByteEncoder->decodeNext(p, value);
                value = maxValue - value;
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
            p = this->vByteEncoder->searchNextPosition(p, divK - (1 + yend - yini));
        }
        return allCells;
    }

    ulong
    K2RasterPlainVByte::getPlainCellValues(uint xini, uint yini, uint xend, uint yend, uint offsetValues, uint divK,
                                           uint maxValue, int *values, ulong totalCells,
                                           uint initXini, uint initYini, uint initYend) const {
        uint value, arrayPosition, cellPos;
        arrayPosition = (offsetValues * divK * divK) + (xini % divK) * divK + (yini % divK);
        uint p = this->vByteEncoder->searchNextPosition(0, arrayPosition);
        for (uint x = xini; x <= xend; x++){
            cellPos = ((x - initXini) * ((initYend - initYini) + 1));
            for (uint y= yini; y <= yend; y++){
                p = this->vByteEncoder->decodeNext(p, value);
                value = maxValue - value;
                values[cellPos + (y - initYini) ] = value;
                totalCells++;
                arrayPosition++;
            }
            p = this->vByteEncoder->searchNextPosition(p, divK - (1 + yend - yini));
        }
        return totalCells;
    }


    //**********************************************************************//
    //********************* FILE FUNCTIONS *********************************//
    //**********************************************************************//
    void K2RasterPlainVByte::save(std::ofstream &of) const {
        K2Raster::save(of, K2RASTER_PLAIN_VBYTE);

        // Levelsk2
        saveValue(of, this->levelK2);

        // Save last values
        this->vByteEncoder->save(of);
    }

    /*
     * Load Ktreap from file
     */
    K2RasterPlainVByte *K2RasterPlainVByte::load(std::ifstream &in) {
        K2RasterPlainVByte *ktreap = nullptr;
        try {
            ktreap = new K2RasterPlainVByte();
            K2Raster::loadBase(in, ktreap);

            // Load k2 value
            ktreap->levelK2 = loadValue<uint>(in);

            // load last values
            ktreap->vByteEncoder = libencoders_static::VByte::load(in);
        } catch (...) {
            return nullptr;
        }
        return ktreap;
    }

    //**********************************************************************//
    //********************* SIZE FUNCTIONS *********************************//
    //**********************************************************************//
    size_t K2RasterPlainVByte::getTotalSize() const {
        size_t totalSize = sizeof(uint)             // Type
                           + sizeof(uint)           // realSize
                           + (sizeof(uint) * 4)     // K1, K2, levelK1, levelK2
                           + (sizeof(uint) * 2);    // initialMaxValue and initialMinValue

        // k2-tree
        totalSize += this->bitmapK2Tree->getSize(); // bitMapsK2tree

        // DACs
        totalSize += this->getSizeDACs(this->listMaxValues, this->numOfMaxDACs); // Max DACs
        totalSize += this->getSizeDACs(this->listMinValues, this->numOfMinDACs); // Min DACs

        //VByte last values
        totalSize += this->vByteEncoder->getSize();
        return totalSize;
    }

    size_t K2RasterPlainVByte::printSize() const {
        printf("RealSize \t-> %lu\n", sizeof(uint));
        printf("K1, K2, levelK1, levelK2 \t-> %lu\n", sizeof(uint) * 4);
        printf("InitialValues (Max and min) \t-> %lu\n", sizeof(uint) * 2);
        printf("BitMapK2Tree \t-> %lu\n", this->bitmapK2Tree->getSize());


        this->printSizeDACs(this->listMaxValues, this->numOfMaxDACs, false);
        this->printSizeDACs(this->listMinValues, this->numOfMinDACs, false);

        // Last values
        printf("VByte last values -> %lu\n", this->vByteEncoder->getSize());
        return this->getTotalSize();
    }
}

