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

#ifndef K2_RASTER_GETOBJECTSQUERYSTRATEGY_H_
#define K2_RASTER_GETOBJECTSQUERYSTRATEGY_H_

// Types
#define STRATEGYQUEUE 0
#define STRATEGYQUEUEV2 6
#define STRATEGYTREE 1
#define STRATEGYTREEV2 2
#define STRATEGYTREEV3 3
#define STRATEGYLEAVES 4
#define STRATEGYBYCELLS 5
#define STRATEGYVALUESQUEUE 10


// System libraries
#include <queue>
#include <stack>
#include <list>
#include <tuple>
// Own header
#include <k2raster/k2-raster.h>
// Libraries
#include <include/spatialindex/SpatialIndex.h>

using namespace SpatialIndex;
using namespace k2raster_static;
using namespace std;

namespace k2raster_rtree_static {

    class GetObjectsQueryStrategy : public SpatialIndex::IQueryStrategy {
    public:
        /****** Getters and setters ******/
//        ulong getNumberOfWhites(){return this->numberOfWhites;};
//        ulong getNumberOfBlacks(){return this->numberOfBlacks;};
        void setRaster(K2Raster *raster);
        void setRangeOfValues(uint valini, uint valend);
        void setAllCells(bool allCells);
        virtual void clear();

        /****** Result ******/
        list<tuple<id_type, uint *, ulong> *> getResult() const;

        ulong getAllMBRCells() const;

    protected:
        // k2-raster where search the values of the MBRs
        K2Raster *raster;
        uint sizeX;
        uint sizeY;

        // valini and valend of query
        uint valini;
        uint valend;
        bool allCell;

        // r-tree params to transform r-tree position into raster position
        double offsetX;
        double offsetY;
        double treeSizeX;
        double treeSizeY;

        /****** Temporal (only for testing) ******/
//        ulong numberOfWhites;
//        ulong numberOfBlacks;

        /****** Results ******/
        list<tuple<id_type, uint *, ulong> *> result;
        ulong allMBRCells;

        /****** AUX ******/
        void calculatePositionParams(double *m_pLow, double *m_pHigh);

        void calculatePositions(double *m_pLow, double *m_pHigh,
                                uint &xini, uint &xend, uint &yini, uint &yend);

        void processChildren(const INode *n, MinBasicSubMatrix *minMatrix);

    }; // END CLASS GetObjectsQueryStrategy

} // END NAMESPACE k2raster_rtree_static

#endif // K2_RASTER_GETOBJECTSQUERYSTRATEGY_H_
