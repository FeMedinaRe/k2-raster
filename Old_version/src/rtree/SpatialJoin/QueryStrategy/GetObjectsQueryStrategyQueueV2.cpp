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

#include <rtree/SpatialJoin/QueryStrategy/GetObjectsQueryStrategyQueueV2.h>

namespace k2raster_rtree_static {

    void GetObjectsQueryStrategyQueueV2::getNextEntry(const IEntry &entry, id_type &nextEntry, bool &hasNext) {
        IShape *ps;
        entry.getShape(&ps);
        Region *pr = dynamic_cast<Region *>(ps);
        MinSubMatrix *minMatrix;

        if (this->minSubMatrices.empty()) {
            // Calculate r-tree params (only once)
            this->calculatePositionParams(pr->m_pLow, pr->m_pHigh);
            // Create root node
            minMatrix = this->raster->createRootMinSubMatrix(this->valini, this->valend);
        } else {
            // Get next MBR id to process
            minMatrix = minSubMatrices.top();
            minMatrix->usedBy--;
            if (this->nChildren.top() == 1){
                this->minSubMatrices.pop();
                this->nChildren.pop();
            } else {
                this->nChildren.top()--;
            }
        }

        uint xini, xend, yini, yend;
        this->calculatePositions(pr->m_pLow, pr->m_pHigh, xini, xend, yini, yend);
        // Search the minimum matrix containing the region of the MBR completely
        if (minMatrix->color == SubMatrixColor::COLOR_SUBMATRIX_GREY) {
            minMatrix = raster->getMinMatrix(xini, xend, yini, yend, this->valini, this->valend, minMatrix);
        }
        delete ps;

        // Check the color of the submatrix
        switch (minMatrix->color) {
            case SubMatrixColor::COLOR_SUBMATRIX_WHITE:     // The MBR is invalid (its values do not meet the requirements)
                // Do not process its children, continue with the next MBR
//                this->numberOfWhites++;
                break;
            default:
                const INode *n = dynamic_cast<const INode *>(&entry);
                if (n != 0) {
                    if (n->isLeaf()) {
                        // Search children
                        this->processChildren(n, minMatrix);
                    } else {
                        // Internal node
                        // Traverse the MBR and add its children to be processed
                        if (n->getChildrenCount() > 0) {
                            this->minSubMatrices.push(minMatrix);
                            this->nChildren.push(n->getChildrenCount());
                            for (uint32_t cChild = 0; cChild < n->getChildrenCount(); cChild++) {
                                this->ids.push(n->getChildIdentifier(cChild));
                                minMatrix->usedBy++;
                            }
                        }
                    } // ENF IF leaf
                }
                break;
        }

        // Check if we can free the initial minMatrix
        if (minMatrix->usedBy == 0) {
            free(minMatrix);
        }

        // Select the next MBR to be processed
        if (!ids.empty()) {
            nextEntry = ids.top();
            ids.pop();
            hasNext = true;
        } else {
            // No more MBR in the list, finish the process
            hasNext = false;
        }
    }

}
