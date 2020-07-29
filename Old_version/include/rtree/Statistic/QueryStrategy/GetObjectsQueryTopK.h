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



#ifndef K2_RASTER_GETOBJECTSQUERYTOP_K_H_
#define K2_RASTER_GETOBJECTSQUERYTOP_K_H_

// System libraries
// Own header
#include <rtree/GetObjectsQueryBasic.h>
#include <k2raster/util/InfoMinSubMatrix.h>
// Libraries
#include <queue>


namespace k2raster_rtree_static {
    struct cmpPriorityPointerMBRMax {
        bool operator()(const priorityPointerMBR* lhs, const priorityPointerMBR* rhs) const {
            return lhs->value < rhs->value;
        }
    };

    struct cmpPriorityPointerMBRMin {
        bool operator()(const priorityPointerMBR* lhs, const priorityPointerMBR* rhs) const {
            return lhs->value > rhs->value;
        }
    };

    template <class T>
    class GetObjectsQueryTopK : public GetObjectsQueryBasic {
    public:

        /****** Entries ******/
        virtual void getNextEntry(const IEntry &entry, id_type &nextEntry, bool &hasNext);

        /****** Getters and setters ******/
        void setK(uint k) {this->k = k;};
        void setWithPositions(bool withPositions) {this->withPositions = withPositions;};
        virtual void initialize();
        virtual void clear();

        // Result
        uint getNumResults() {return this->numResults;};
        int* getValues() {return this->values;};
        id_type* getMBRIds() {return this->MBRIds;};

        /****** Result ******/

    protected:
        // ****** Params ****** //
        uint k;
        bool withPositions;

        // ****** MBR Ids ****** //
        std::priority_queue<priorityPointerMBR*, vector<priorityPointerMBR*>, T> pointers;
        BasicSubMatrixWithUsed *nextPointer;

        // ****** Result ****** //
        uint numResults;
        int *values;
        id_type *MBRIds;
        uint **positions;
        ulong *numOfPositions;


        /****** AUX ******/
        virtual void setNextEntry(id_type &nextEntry, bool &hasNext);
        virtual int getMBRLeafValue(uint xini, uint xend, uint yini, uint yend, BasicSubMatrix *basicSubMatrix)=0;
        virtual int getMBRValue(BasicSubMatrix *basicSubMatrix)=0;
        virtual ulong getCells(priorityPointerMBR *ptr)=0;
    };
} // END NAMESPACE k2raster_rtree_static


#endif // K2_RASTER_GETOBJECTSQUERYTOP_K_H_
