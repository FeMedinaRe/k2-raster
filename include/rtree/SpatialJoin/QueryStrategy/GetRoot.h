/*  
 * Created by Fernando Silva on 17/01/17.
 *
 * Copyright (C) 2017-current-year, Fernando Silva, all rights reserved.
 *
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

#ifndef K2_RASTER_GETROOT_H_
#define K2_RASTER_GETROOT_H_

// Libraries
#include <include/spatialindex/SpatialIndex.h>

using namespace SpatialIndex;

namespace k2raster_rtree_static {

    class GetRoot : public SpatialIndex::IQueryStrategy {
    public:
        /****** Entries ******/
        virtual void getNextEntry(const IEntry &entry, id_type &nextEntry, bool &hasNext);

        /****** Getters and setters ******/
        double *getPlow();

        double *getPhigh();

    private:
        double plow[2];
        double phigh[2];

    };
}


#endif // K2_RASTER_GETROOT_H_
