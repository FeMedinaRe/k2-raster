/*  
 * Created by Fernando Silva on 30/10/17.
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



#ifndef K2_RASTER_GETOBJECTSQUERYMINTOP_K_H_
#define K2_RASTER_GETOBJECTSQUERYMINTOP_K_H_

// System libraries
// Own header
#include <rtree/Statistic/QueryStrategy/GetObjectsQueryTopK.h>
// Libraries


namespace k2raster_rtree_static {

    class GetObjectsQueryMinTopK : public GetObjectsQueryTopK<cmpPriorityPointerMBRMin> {
    public:


    protected:
        /****** AUX ******/
        virtual int getMBRLeafValue(uint xini, uint xend, uint yini, uint yend, BasicSubMatrix *basicSubMatrix);
        virtual int getMBRValue(BasicSubMatrix *basicSubMatrix);
        virtual ulong getCells(priorityPointerMBR *ptr);
    };
} // END NAMESPACE k2raster_rtree_static


#endif // K2_RASTER_GETOBJECTSQUERYMINTOP_K_H_
