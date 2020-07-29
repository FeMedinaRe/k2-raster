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

#include <rtree/SpatialJoin/QueryStrategy/GetRoot.h>

namespace k2raster_rtree_static {

    //**********************************************************************//
    //*************************** ENTRIES **********************************//
    //**********************************************************************//

    /****** Process MBR ******/
    void GetRoot::getNextEntry(const IEntry &entry, id_type &nextEntry, bool &hasNext) {
        IShape *ps;
        entry.getShape(&ps);
        Region *pr = dynamic_cast<Region *>(ps);
        this->plow[0] = pr->m_pLow[0];
        this->plow[1] = pr->m_pLow[1];
        this->phigh[0] = pr->m_pHigh[0];
        this->phigh[1] = pr->m_pHigh[1];

        hasNext = false;
    }

    //**********************************************************************//
    //******************* SETTERS AND GETTERS ******************************//
    //**********************************************************************//
    double *GetRoot::getPlow() {
        return this->plow;
    }

    double *GetRoot::getPhigh() {
        return this->phigh;
    }
}

