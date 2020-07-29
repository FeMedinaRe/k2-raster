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

#include <rtree/SpatialJoin/QueryStrategy/GetObjectsQueryStrategyTreeV3.h>

namespace k2raster_rtree_static {

    void GetObjectsQueryStrategyTreeV3::getNextEntry(const IEntry &entry, id_type &nextEntry, bool &hasNext) {
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
                const INode *n = dynamic_cast<const INode *>(&entry);
                if (n != 0) {
                    if (n->isLeaf()) {
                        // Search children
                        for (uint32_t cChild = 0; cChild < n->getChildrenCount(); cChild++) {
                            // Search individual object
                            IShape *ps2;
                            Region re2;
                            n->getChildShape(cChild, &ps2);
//                            re2 = dynamic_cast<Region *>(ps2);
                            ps2->getMBR(re2);
                            uint xini, xend, yini, yend;
                            this->calculatePositions(re2.m_pLow, re2.m_pHigh, xini, xend, yini, yend);
                            delete(ps2);
//                            printf("Positions -> [%u, %u][%u, %u]\n", xini, yini, xend, yend);

                            // Alloc array os positions
                            ulong totalCells = 0;
                            uint *positionsFinal = nullptr;
                            ulong positionsSize = ((xend - xini + 1) * (yend - yini + 1));
                            uint *positions = (uint *) malloc( positionsSize * sizeof(uint));

                            // Search cells
                            if (currentNode->color ==
                                SubMatrixColor::COLOR_SUBMATRIX_BLACK) {  // The MBR is valid (and all objects within of the MBR)
                                for (uint x = xini; x <= xend; x++) {
                                    for (uint y = yini; y <= yend; y++) {
                                        positions[totalCells++] = x * raster->getDimSize(2) + y;
                                    } // END FOR y
                                } // ENF FOR x
                            } else {
                                MinNodeTreeNode *leafNode = raster->getMinMatrixV2(xini, xend, yini, yend,
                                                                                        this->valini, this->valend, currentNode);
                                switch (leafNode->color) {
                                    case SubMatrixColor::COLOR_SUBMATRIX_WHITE: break;
                                    case SubMatrixColor::COLOR_SUBMATRIX_BLACK:
                                        for (uint x = xini; x <= xend; x++) {
                                            for (uint y = yini; y <= yend; y++) {
                                                positions[totalCells++] = x * raster->getDimSize(2) + y;
                                            } // END FOR y
                                        } // ENF FOR x
                                        break;
                                    case SubMatrixColor::COLOR_SUBMATRIX_GREY:
                                        totalCells = raster->getCellsByValue(xini, xend, yini, yend, this->valini, this->valend,
                                                                             leafNode, positions, this->allCell);
                                        break;
                                    default:
                                        break;
                                }

                                if (totalCells == 0) {
                                    free(positions); // No valid cells
                                    continue;
                                } else {
                                    if (positionsSize != totalCells) {
                                        positionsFinal = (uint *) realloc(positions, totalCells);
                                        if (!positionsFinal) {
                                            printf("Error realloc -> MBR ID: %lu, cells %lu\n", n->getIdentifier(),
                                                   totalCells);
                                            free(positions);
                                            exit(-1);
                                        };
                                    } else {
                                        positionsFinal = positions;
                                    } // END IF positionsSize != totalCells
                                } // END IF totalCells == 0
                            } // END IF check color black

                            // Add result
                            this->result.push_back(
                                    new tuple<id_type, uint *, ulong>(n->getIdentifier(), positions, totalCells));
                            this->allMBRCells += totalCells;
                        } // END FOR children
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
//        if (!hasNext && this->totalMBRsCreated > 0) {
//            printf("Used %lu MBRs of %lu (%.3f\%)\n", this->totalMBRsCreated - this->totalMBRsNotUsed, this->totalMBRsCreated,
//                   (this->totalMBRsCreated - this->totalMBRsNotUsed) * 100.0 / this->totalMBRsCreated);
//        }
    }


    void GetObjectsQueryStrategyTreeV3::clear() {
        GetObjectsQueryStrategy::clear();
        this->next = nullptr;
    }
}
