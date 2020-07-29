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



#ifndef K2_RASTER_GETOBJECTSQUERYMAXTOP_K_H_OLD
#define K2_RASTER_GETOBJECTSQUERYMAXTOP_K_H_OLD

// System libraries
// Own header
#include <rtree/GetObjectsQueryBasic.h>
#include <k2raster/util/InfoMinSubMatrix.h>
// Libraries
#include <queue>


namespace k2raster_rtree_static {

    struct cmpPriorityPointerMBR {
        bool operator()(const priorityPointerMBR* lhs, const priorityPointerMBR* rhs) const {
            return lhs->value < rhs->value;
        }
    };


    class GetObjectsQueryMaxTopKOLD : public GetObjectsQueryBasic {
    public:

        /****** Entries ******/
        virtual void getNextEntry(const IEntry &entry, id_type &nextEntry, bool &hasNext);
        virtual void setNextEntry(id_type &nextEntry, bool &hasNext);

        /****** Getters and setters ******/
        void setK(uint k) {this->k = k;};
        void setWithPositions(bool withPositions) {this->withPositions = withPositions;};
        virtual void initialize();
        virtual void clear();

        // Result
        uint getNumResults() {return this->numResults;};
        int* getMaxValues() {return this->maxValues;};
        id_type* getMBRIds() {return this->MBRIds;};

        /****** Result ******/

    protected:
        // ****** Params ****** //
        uint k;
        bool withPositions;

        // ****** MBR Ids ****** //
        std::priority_queue<priorityPointerMBR*, vector<priorityPointerMBR*>, cmpPriorityPointerMBR> pointers;
        BasicSubMatrixWithUsed *nextPointer;

        // ****** Result ****** //
        uint numResults;
        int *maxValues;
        id_type *MBRIds;
        uint **positions;
        uint *numOfPositions;


        /****** AUX ******/
    };
} // END NAMESPACE k2raster_rtree_static


#endif // K2_RASTER_GETOBJECTSQUERYMAXTOP_K_H_OLD
