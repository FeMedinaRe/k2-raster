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

#include <rtree/SpatialJoin/QueryStrategy/GetObjectsQueryStrategyTree.h>

namespace k2raster_rtree_static {

    void GetObjectsQueryStrategyTree::getNextEntry(const IEntry &entry, id_type &nextEntry, bool &hasNext) {
        IShape *ps;
        entry.getShape(&ps);
        Region *pr = dynamic_cast<Region *>(ps);

        MinNodeTreeNode *currentNode;
        if (this->nextMBRs.empty()) {
            // Calculate r-tree params (only once)
            this->calculatePositionParams(pr->m_pLow, pr->m_pHigh);
            // Create root node
            currentNode = this->raster->createRootNodeMinSubMatrix(this->valini, this->valend);
            this->nextMBRs[currentNode->bitmapPosition] = currentNode;
            this->next = this->nextMBRs.begin();
        } else {
            currentNode = this->next->second;
        }

        uint xini, xend, yini, yend;
        this->calculatePositions(pr->m_pLow, pr->m_pHigh, xini, xend, yini, yend);
        if (currentNode->color == SubMatrixColor::COLOR_SUBMATRIX_GREY) {
            currentNode = raster->getMinMatrix(xini, xend, yini, yend, this->valini, this->valend,
                                               currentNode,
                                               nextMBRs);
        }
        delete ps;

        // Check the color of the submatrix
        switch (currentNode->color) {
            case SubMatrixColor::COLOR_SUBMATRIX_WHITE:     // The MBR is invalid (its values do not meet the requirements)
                // Do not process its children, continue with the next MBR
                break;
            default:
                const INode *n = dynamic_cast<const INode *>(&entry);
                if (n != 0) {
                    if (n->isLeaf()) {
                        // Search children
                        this->processChildren(n, currentNode);
                    } else {
                        // Internal node
                        // Traverse the MBR and add its children to be processed
                        if (currentNode->queue == nullptr) {
                            currentNode->queue = new std::queue<int64_t>();
                        }
                        for (uint32_t cChild = 0; cChild < n->getChildrenCount(); cChild++) {
                            currentNode->queue->push(n->getChildIdentifier(cChild));
                        }
                    } // ENF IF leaf
                }
                break;
        }

        hasNext = false;
        map<ulong, MinNodeTreeNode *>::iterator initPos = this->next;
        for (map<ulong, MinNodeTreeNode *>::iterator it = initPos; it != nextMBRs.end(); ++it) {
            if (it->second->queue == nullptr || it->second->queue->empty()) {
                // Free node. It has no any MBR to process
                if (it->second->queue != nullptr) {
                    delete (it->second->queue);
                }
//                this->totalMBRsCreated++;
//                if (!it->second->used) {
//                    this->totalMBRsNotUsed++;
//                }
                free(it->second);
            } else {
                this->next = it;
                hasNext = true;
                nextEntry = it->second->queue->front();
                it->second->queue->pop();
                break;
            }
        }
//        if (!hasNext) {
//            printf("Used %lu MBRs of %lu (%.3f\%)\n", this->totalMBRsCreated - this->totalMBRsNotUsed, this->totalMBRsCreated,
//                   (this->totalMBRsCreated - this->totalMBRsNotUsed) * 100.0 / this->totalMBRsCreated);
//        }
    }


    void GetObjectsQueryStrategyTree::clear() {
        GetObjectsQueryStrategy::clear();
        this->nextMBRs.clear();
    }
}
