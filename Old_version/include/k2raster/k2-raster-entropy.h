/*
 * Created by Fernando Silva on 13/06/16.
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

#ifndef K2RASTER_ENTROPY_H_H_
#define K2RASTER_ENTROPY_H_

#include <k2raster/plain/k2-raster-plain.h>
#include <include/hash/HashTableWord.h>

namespace k2raster_static {

/****** Aux hash functions ******/
int compareFreqListDescHEntropy(const void *a, const void *b);
bool compareFreqListDescHEntropy2(ulong left, ulong right);

class K2RasterEntropy : public K2RasterPlain {
  public:
    /****** Constructors ******/
    K2RasterEntropy(uint sizeX, uint size, int **data, uint K1, uint K2, uint levelK1, uint plainLevels);
    ~K2RasterEntropy();

    /****** Cell functions ******/
    // GetCell
    virtual int getPlainCell(uint father, uint row, uint col, uint divK, uint lastValueMax) const;

    // GetCellsByValue
    virtual ulong getPlainWindow(uint xini, uint yini, uint xend, uint yend,
                                 uint valini, uint valend,
                                 uint offsetValues, uint divK, uint maxValue,
                                 uint *positions, ulong totalCells, bool allCells) const;

    // checkValues
    virtual bool checkPlainCells(uint xini, uint yini, uint xend, uint yend,
                                 uint offsetValues, uint divK,uint maxValue,
                                 Query query, bool allCells) const;

    // getValuesWindow
    virtual ulong getPlainCellValues(uint xini, uint yini, uint xend, uint yend,
                                     uint offsetValues, uint divK, uint maxValue,
                                     int *values, ulong totalCells,
                                     uint initXini, uint initYini, uint initYend) const;

    /****** File functions ******/
    virtual void save(std::ofstream &of) const;
    virtual void saveSplit(char *filename) const;
    static K2RasterEntropy *load(std::ifstream &in);

    /****** Size functions ******/
    virtual size_t getTotalSize() const;
    virtual size_t printSize() const;

  protected:
    K2RasterEntropy();

    // Attributes
    uint lenWord;
    uint numOfWords;                            // Number of words of the vocabulary
    uint totalLengthOfWords;                    // Length of the vocabulary
    unsigned char* words;                       // Vocabulary
    cds_static::BitSequence *bitmapVocLeaves;  // k2-tree representation of the data: single bitmap for T:L
    cds_static::DirectAccess *vocLeaves;     // Encoded vocabulary of leaves
    cds_static::DirectAccess *valuesLeaves;       // Encode values of leaves that they are not in the vocabulary

};
}


#endif // KTREAP_RASTER_COMPRESS_LEAVES_H_
