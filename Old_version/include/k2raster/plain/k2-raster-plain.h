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

#ifndef K2RASTER_PLAIN_H_
#define K2RASTER_PLAIN_H_

#include <k2raster/k2-raster.h>

namespace k2raster_static {

class K2RasterPlain : public K2Raster {
  public:
    /****** Constructors ******/
    K2RasterPlain(uint sizeX, uint size, int **data, uint K1, uint K2, uint levelK1, uint plainLevels);
    ~K2RasterPlain();

    /****** Cell functions ******/
    // GetCell
    virtual int getCell(uint row, uint col) const;
    virtual int getPlainCell(uint father, uint row, uint col, uint divK, uint lastValueMax) const;

    // GetCellsByValue
    virtual ulong getCellsByValue(Query query, uint *positions) const;
    virtual ulong getCellsByValue(uint xini, uint xend, uint yini, uint yend, uint valini, uint valend,
                                  MinBasicSubMatrix *minSubMatrix, uint *positions, bool allCells) const;

    virtual ulong getPlainWindow(uint xini, uint yini, uint xend, uint yend,
                                 uint offsetValues, uint divK, uint maxValue,
                                 uint *positions, ulong totalCells, Query query) const;

    virtual ulong getPlainWindow(uint xini, uint yini, uint xend, uint yend,
                                 uint valini, uint valend,
                                 uint offsetValues, uint divK, uint maxValue,
                                 uint *positions, ulong totalCells, bool allCells) const;

    // checkValues
    virtual bool checkValues(Query query, bool allCells) const;
    virtual bool checkPlainCells(uint xini, uint yini, uint xend, uint yend,
                                 uint offsetValues, uint divK, uint maxValue,
                                 Query query, bool allCells) const;

    // getValuesWindow
    virtual ulong getValuesWindow(uint xini, uint xend, uint yini, uint yend, int *values) const;

    virtual ulong getPlainCellValues(uint xini, uint yini, uint xend, uint yend,
                                     uint offsetValues, uint divK, uint maxValue,
                                     int *values, ulong totalCells, Query query) const;
    virtual ulong getPlainCellValues(uint xini, uint yini, uint xend, uint yend,
                                     uint offsetValues, uint divK, uint maxValue,
                                     int *values, ulong totalCells,
                                     uint initXini, uint initYini, uint initYend) const;

    // getMinMatrix
//    virtual minSubMatrix *getMinMatrix(uint xini, uint xend, uint yini, uint yend,
//                                       uint valini, uint valend, minSubMatrix *minSubMatrix) const;
//
//    // V2 (using a tree of precalculated submatrices)
//    virtual tNodeMinSubMatrix *getMinMatrix(uint xini, uint xend, uint yini, uint yend,
//                                            uint valini, uint valend,
//                                            ulong bitmapPosition,  std::map<ulong, tNodeMinSubMatrix*> &nextMBRs) const;

    /****** File functions ******/
    virtual void save(std::ofstream &of) const;
    virtual void saveSplit(char *filename) const;
    static K2RasterPlain *load(std::ifstream &in);

    /****** Size functions ******/
    virtual size_t getTotalSize() const;
    virtual size_t printSize() const;

  protected:
    K2RasterPlain();

    // Attributes
    uint levelK2;
    uint *plainValues;
    uint numOfPlainValues;

    /****** Cell aux functions ******/
    virtual ulong getCellsByValue(uint xini, uint xend, uint yini, uint yend, uint fatherMinValue, uint fatherMaxValue,
                                  uint l, uint bitmapChildren, uint valini, uint valend,
                                  uint baseRow, uint baseCol,
                                  ulong totalCells, uint *positions, bool allCells) const;

    /****** Auxiliary functions ******/
    void initLevels(uint plainLevels);
    virtual bool isLevelInPlain(uint level) const;
};
}


#endif // K2RASTER_PLAIN_H_
