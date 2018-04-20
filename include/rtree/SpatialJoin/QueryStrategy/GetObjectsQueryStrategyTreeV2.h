/*  
 * Created by Fernando Silva on 2/01/17.
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

#ifndef K2_RASTER_GETOBJECTSQUERYSTRATEGY_TREEV2_H_
#define K2_RASTER_GETOBJECTSQUERYSTRATEGY_TREEV2_H_

// System libraries
#include <stack>
#include <list>
#include <tuple>
// Own header
#include <k2raster/k2-raster.h>
#include <k2raster/util/InfoMinSubMatrix.h>
#include <rtree/SpatialJoin/QueryStrategy/GetObjectsQueryStrategy.h>
// Libraries


using namespace SpatialIndex;
using namespace k2raster_static;
using namespace std;

namespace k2raster_rtree_static {

    class GetObjectsQueryStrategyTreeV2 : public GetObjectsQueryStrategy {

    public:
        /****** Process MBR ******/
        virtual void getNextEntry(const IEntry &entry, id_type &nextEntry, bool &hasNext);

        /****** Getters and setters ******/
        virtual void clear();

    private:
        // Queues to store the next MBR
        stack<MinNodeTreeNode *> nextMBRs;
        MinNodeTreeNode *next;

        // Temp experimental
//        ulong totalMBRsCreated;
//        ulong totalMBRsNotUsed;

    }; // END CLASS GetObjectsQueryStrategy

} // END NAMESPACE k2raster_rtree_static

#endif // K2_RASTER_GETOBJECTSQUERYSTRATEGY_TREEV2_H_
