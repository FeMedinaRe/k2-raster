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

#ifndef K2RASTER_OPT_H_
#define K2RASTER_OPT_H_

#include <k2raster/k2-raster.h>
#include <k2raster/util/NodeMatrix.h>

namespace k2raster_static {

class K2RasterOPT : public K2Raster {
public:
    // Constructors
    K2RasterOPT(uint sizeX, uint size, int **data, uint K1, uint k2, uint levelK1);

    ~K2RasterOPT();

    // Cell functions
//    virtual int getCell(uint row, uint col) const;

    virtual ulong getCellsByValue(Query query, uint *positions) const;
    virtual bool checkValues(Query query, bool allCells) const;

    virtual ulong getValuesWindow(Query query, int *values) const;

    // File functions
    virtual void save(std::ofstream &of) const;

    static K2RasterOPT *load(std::ifstream &in);

protected:
    // Constructor
    K2RasterOPT();
};
}
#endif // KTREAP_REGION_H_
