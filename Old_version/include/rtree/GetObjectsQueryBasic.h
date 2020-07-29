/*  
 * Created by Fernando Silva on 25/10/17.
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



#ifndef K2_RASTER_GETOBJECTSQUERYBASIC_H_
#define K2_RASTER_GETOBJECTSQUERYBASIC_H_

// System libraries
// Own header
#include <k2raster/k2-raster.h>
// Libraries
#include <include/spatialindex/SpatialIndex.h>


using namespace SpatialIndex;
using namespace k2raster_static;
using namespace std;

namespace k2raster_rtree_static {
    class GetObjectsQueryBasic : public SpatialIndex::IQueryStrategy {
    public:

        /****** Getters and setters ******/
        void setRaster(K2Raster *raster);
        virtual void initialize();
        virtual void clear();

        /****** Result ******/

    protected:
        // k2-raster where search the values of the MBRs
        K2Raster *raster;
        uint sizeX;
        uint sizeY;

        // R-tree params to transform r-tree position into raster position
        double offsetX;
        double offsetY;
        double treeSizeX;
        double treeSizeY;

        // Others
        bool first; // Initially true (we will process the root node)

        /****** AUX ******/
        void calculatePositionParams(double *m_pLow, double *m_pHigh);
        void calculatePositions(double *m_pLow, double *m_pHigh,
                                uint &xini, uint &xend, uint &yini, uint &yend);
    };
} // END NAMESPACE k2raster_rtree_static


#endif // K2_RASTER_GETOBJECTSQUERYBASIC_H_
