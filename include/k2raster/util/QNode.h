/*
 * Created by Fernando Silva on 11/02/16.
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

#ifndef UTIL_QNODE_H_
#define UTIL_QNODE_H_

#include <sys/types.h>

namespace k2raster_static    {

    class QNode {

    public:
        QNode() {};
        QNode(uint minValue, uint maxValue);
        QNode(uint position, uint level, uint baseRow, uint baseCol, uint minValue, uint maxValue);

        uint position;
        uint level;
        uint baseRow;
        uint baseCol;
        uint maxValue;
        uint minValue;
    };


    class QNodePositions {

    public:
        QNodePositions() {};
        QNodePositions(uint position, uint level, uint baseRow, uint baseCol,
                       uint maxValue, uint row, uint col);
        uint position;
        uint level;
        uint baseRow;
        uint baseCol;
        uint maxValue;
        uint row;
        uint col;
    };

    class QNodeCounters {

    public:
        QNodeCounters() {};
        QNodeCounters(uint position, uint level, uint baseRow, uint baseCol, uint minValue, uint maxValue);

        uint position;
        uint level;
        uint baseRow;
        uint baseCol;
        uint maxValue;
        uint minValue;
        uint preValuesMin;
        uint preValuesMax;
        uint offsetPrev;
    };
}

#endif // UTIL_QNODE_H_
