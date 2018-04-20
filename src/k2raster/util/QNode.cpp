/*
 * Created by Fernando Silva on 15/02/16.
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

#include <k2raster/util/QNode.h>

namespace k2raster_static {

    QNode::QNode(uint minValue, uint maxValue) {
            this->position = 1;
            this->minValue = minValue;
            this->maxValue = maxValue;
            this->level = 0;
            this->baseRow = 0;
            this->baseCol = 0;
    }

    QNode::QNode(uint position, uint level, uint baseRow, uint baseCol, uint minValue, uint maxValue) {
        this->position = position;
        this->level = level;
        this->baseRow = baseRow;
        this->baseCol = baseCol;
        this->minValue = minValue;
        this->maxValue = maxValue;
    }

    QNodePositions::QNodePositions(uint position, uint level, uint baseRow, uint baseCol, uint maxValue,
                                   uint row, uint col) {
        this->position = position;
        this->level = level;
        this->baseRow = baseRow;
        this->baseCol = baseCol;
        this->maxValue = maxValue;
        this->row = row;
        this->col = col;
    }

    QNodeCounters::QNodeCounters(uint position, uint level, uint baseRow, uint baseCol, uint minValue,
                                 uint maxValue) {
        this->position = position;
        this->level = level;
        this->baseRow = baseRow;
        this->baseCol = baseCol;
        this->minValue = minValue;
        this->maxValue = maxValue;
    }

}

