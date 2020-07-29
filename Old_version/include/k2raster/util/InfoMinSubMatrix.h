/*  
 * Created by Fernando Silva on 5/01/17.
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

#ifndef K2_RASTER_INFOMINSUBMATRIX_H_
#define K2_RASTER_INFOMINSUBMATRIX_H_

#include <queue>
#include <vector>
#include <sys/types.h>

namespace k2raster_static {



    // ************************************************* //
    // ********************* BASIC ********************* //
    // ************************************************* //

    struct BasicSubMatrix {
        ulong bitmapPosition;
        ulong bitmapChildren;
        uint level;
        uint minValue; // TODO remove this one
        uint maxValue;
    };

    struct BasicSubMatrixWithChildren : public  BasicSubMatrix {
        BasicSubMatrixWithChildren **children;
        uint numOfChildren;
    };

    struct BasicSubMatrixWithUsed : public  BasicSubMatrix {
        uint usedBy;
    };

    // Functions
    void removeChildren(BasicSubMatrixWithChildren *matrix);


    // ************************************************* //
    // ***************** SPATIAL JOIN ****************** //
    // ************************************************* //

    // COLOR_SUBMATRIX_WHITE -> All cells without of range of values
    // COLOR_SUBMATRIX_BLACK -> All cells within of range of values
    // COLOR_SUBMATRIX_GREY  -> We cannot determinate it
    enum SubMatrixColor {
        COLOR_SUBMATRIX_WHITE, COLOR_SUBMATRIX_BLACK, COLOR_SUBMATRIX_GREY
    };

    struct MinBasicSubMatrix : public BasicSubMatrix {
        ushort color;
    };

    struct MinSubMatrix : public MinBasicSubMatrix {
        uint usedBy;
    };

    struct MinNodeTreeNode : public MinBasicSubMatrix {
        uint numberOfChildren;
        MinNodeTreeNode **children;

        // Id to process
        std::queue<int64_t> *queue;
    };


    // ************************************************* //
    // ******************* STATISTIC ******************* //
    // ************************************************* //
    struct PositionsStruct {
        uint xini;
        uint xend;
        uint yini;
        uint yend;
    };

    struct priorityPointerMBR {
        int64_t id;
        bool isRealValue;                   // If the maxValue (pointer.maxValue) is the max value of the MBR (true) or it is a possible value (false)
        int value;
        BasicSubMatrixWithUsed *pointer;
        PositionsStruct *positions;
    };

} // END NAMESPACE k2raster_static


#endif // K2_RASTER_INFOMINSUBMATRIX_H_
