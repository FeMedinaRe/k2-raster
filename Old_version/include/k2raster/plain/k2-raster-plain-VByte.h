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

#ifndef K2RASTER_PLAIN_VBYTE_H_
#define K2RASTER_PLAIN_VBYTE_H_

#include <k2raster/plain/k2-raster-plain.h>
#include <include/VByte/VByte.h>

namespace k2raster_static {

class K2RasterPlainVByte : public K2RasterPlain {
  public:
    /****** Constructors ******/
    K2RasterPlainVByte(uint sizeX, uint size, int **data, uint K1, uint K2, uint levelK1, uint plainLevels);
    ~K2RasterPlainVByte();

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
                                     int *values, ulong totalCells,
                                     uint initXini, uint initYini, uint initYend) const;

    /****** File functions ******/
    virtual void save(std::ofstream &of) const;
    static K2RasterPlainVByte *load(std::ifstream &in);

    /****** Size functions ******/
    virtual size_t getTotalSize() const;
    virtual size_t printSize() const;

  protected:
    K2RasterPlainVByte();

    // Attributes
    libencoders_static::VByte * vByteEncoder;
};
}


#endif // K2RASTER_PLAIN_VBYTE_H_
