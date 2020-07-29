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

#ifndef K2RASTER_COMPRESS_LEAVES_H_
#define K2RASTER_COMPRESS_LEAVES_H_

#include <k2raster/plain/k2-raster-plain.h>
#include <include/hash/HashTableWord.h>

namespace k2raster_static {

/****** Aux hash functions ******/
int compareFreqListDesc(const void *a, const void *b);

class K2RasterCompressLeaves : public K2RasterPlain {
  public:
    /****** Constructors ******/
    K2RasterCompressLeaves(uint sizeX, uint size, int **data, uint K1, uint K2, uint levelK1, uint plainLevels);
    ~K2RasterCompressLeaves();

    /****** Cell functions ******/
    virtual int getPlainCell(uint father, uint row, uint col, uint divK, uint lastValueMax) const;

    virtual ulong getPlainWindow(uint xini, uint yini, uint xend, uint yend,
                                 uint offsetValues, uint divK, uint maxValue,
                                 uint *positions, ulong totalCells, Query query) const;
    virtual bool checkPlainCells(uint xini, uint yini, uint xend, uint yend,
                                 uint offsetValues, uint divK,uint maxValue,
                                 Query query, bool allCells) const;

    virtual ulong getPlainCellValues(uint xini, uint yini, uint xend, uint yend,
                                     uint offsetValues, uint divK, uint maxValue,
                                     int *values, ulong totalCells, Query query) const;

    /****** File functions ******/
    virtual void save(std::ofstream &of) const;
    static K2RasterCompressLeaves *load(std::ifstream &in);

    /****** Size functions ******/
    virtual size_t getTotalSize() const;
    virtual size_t printSize() const;

  protected:
    K2RasterCompressLeaves();

    // Attributes
    uint lenWord;
    uint numOfWords;                      // Number of words of the vocabulary
    uint totalLengthOfWords;                    // Length of the vocabulary
    unsigned char* words;               // Vocabulaby
    cds_static::DirectAccess *valuesLeaves;   // Encoded values of leaves


};
}


#endif // K2RASTER_COMPRESS_LEAVES_H_
