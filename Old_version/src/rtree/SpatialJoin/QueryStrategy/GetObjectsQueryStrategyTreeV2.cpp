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

#include <rtree/SpatialJoin/QueryStrategy/GetObjectsQueryStrategyTreeV2.h>

namespace k2raster_rtree_static {

    void GetObjectsQueryStrategyTreeV2::getNextEntry(const IEntry &entry, id_type &nextEntry, bool &hasNext) {
        IShape *ps;
        entry.getShape(&ps);
        Region *pr = dynamic_cast<Region *>(ps);

        if (this->next == nullptr) {
            // Calculate r-tree params (only once)
            this->calculatePositionParams(pr->m_pLow, pr->m_pHigh);
            this->next = this->raster->createRootNodeMinSubMatrix(this->valini, this->valend, true);
            this->next->queue = new queue<id_type>();
        }

        uint xini, xend, yini, yend;
        this->calculatePositions(pr->m_pLow, pr->m_pHigh, xini, xend, yini, yend);
        MinNodeTreeNode *currentNode = raster->getMinMatrixV2(xini, xend, yini, yend,
                                                                this->valini, this->valend, this->next);
        delete ps;

        // Check the color of the submatrix
        switch (currentNode->color) {
            case SubMatrixColor::COLOR_SUBMATRIX_WHITE:     // The MBR is invalid (its values do not meet the requirements)
                // Do not process its children, continue with the next MBR
                break;
            default:
                const auto *n = dynamic_cast<const INode *>(&entry);
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

        // Check current submatrix
        if (!this->next->queue->empty()) {
            // Get next MBR of this node
            hasNext = true;
            nextEntry = this->next->queue->front();
            this->next->queue->pop();
        } else {
            // No more MBR for current node.
            // Add its children to queue
            if (this->next->children != nullptr) {
                for (uint c = 0; c < this->next->numberOfChildren; c++) {
                    if (this->next->children[c] != nullptr) {
                        this->nextMBRs.push(this->next->children[c]);
                    }
                }
                free(this->next->children);
            }
            // Delete it
            delete (this->next->queue);
            free(this->next);

            // Get next node of the queue
            while (!this->nextMBRs.empty()) {
                this->next = this->nextMBRs.top();
                this->nextMBRs.pop();

                if (this->next->queue == nullptr) {
                    // Free node. It has no any MBR to process
                    if (this->next->children != nullptr) {
                        for (uint c = 0; c < this->next->numberOfChildren; c++) {
                            if (this->next->children[c] != nullptr) {
                                this->nextMBRs.push(this->next->children[c]);
                            }
                        }
                        free(this->next->children);
                    }
                    free(this->next);
                } else {
                    hasNext = true;
                    nextEntry = this->next->queue->front();
                    this->next->queue->pop();
                    break;
                }
            }
        }
    }


    void GetObjectsQueryStrategyTreeV2::clear() {
        GetObjectsQueryStrategy::clear();
        this->next = nullptr;
    }
}
