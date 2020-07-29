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

#include <rtree/Statistic/QueryStrategy/GetObjectsQueryMaxMBR.h>

namespace k2raster_rtree_static {

    //**********************************************************************//
    //*************************** ENTRIES **********************************//
    //**********************************************************************//

    /****** Process MBR ******/
    void GetObjectsQueryMaxMBR::getNextEntry(const IEntry &entry, id_type &nextEntry, bool &hasNext) {
        BasicSubMatrixWithUsed *minSubMatrix;

        if (this->first) {
            this->first = false;

            // Get object
            IShape *ps;
            entry.getShape(&ps);
            Region *pr = dynamic_cast<Region *>(ps);

            // First iteration of getNextEntry
            this->calculatePositionParams(pr->m_pLow, pr->m_pHigh);

            // Pointer to root node of the k2-raster
            this->nextPointer = this->raster->createRootBasicSubMatrixWithUsed();
            this->nextPointer->usedBy = 1;
            uint xini, xend, yini, yend;
            this->calculatePositions(pr->m_pLow, pr->m_pHigh, xini, xend, yini, yend);
            minSubMatrix = this->raster->getMinMatrix(xini, xend, yini, yend, this->nextPointer);
            minSubMatrix->usedBy++;
            this->nextPointer->usedBy--;

            if (this->nextPointer->usedBy == 0) {
                delete(this->nextPointer);
            }
            this->nextPointer = minSubMatrix;

            // Delete MBR
            delete(ps);
        }
        this->nextPointer->usedBy--;

        const INode *n = dynamic_cast<const INode *>(&entry);
        if (n != 0) {
            uint xini, xend, yini, yend;
            for (uint32_t cChild = 0; cChild < n->getChildrenCount(); cChild++) {
                // Search individual object
                IShape *ps2;
                Region re2;
                n->getChildShape(cChild, &ps2);
                ps2->getMBR(re2);

                // ----------------------------------------------------------------------------- //
                // Get the corresponding positions of the R-Tree (pr->m_pLow, pr->m_pHigh)       //
                // in the k2-raster (xini, xend, yini, yend)                                     //
                // ----------------------------------------------------------------------------- //
                this->calculatePositions(re2.m_pLow, re2.m_pHigh, xini, xend, yini, yend);
                delete (ps2);

                // ----------------------------------------------------------------------------- //
                // Get minimum submatrix (of k2-raster) that completely contains the current MBR //
                // ----------------------------------------------------------------------------- //
                minSubMatrix = this->raster->getMinMatrix(xini, xend, yini, yend, this->nextPointer);

                if (n->isLeaf()) {
                    // Search maxValue of the children
                    int value = this->raster->getMaxValueWindow(xini, xend, yini, yend, minSubMatrix);

                    if (minSubMatrix != this->nextPointer) {
                        delete (minSubMatrix);
                    }

                    if (value > this->maxValue) {
                        this->maxValue = value;
                        this->MBRId = n->getChildIdentifier(cChild);
                    }
                    if (value == this->nextPointer->maxValue) { // Higher value of parent .There is no other value higher than the current one
                        hasNext = false;
                        return;
                    }
                } else {
                    minSubMatrix->usedBy++;

                    priorityPointerMBR *childPointer = new priorityPointerMBR;
                    childPointer->id = n->getChildIdentifier(cChild);
                    childPointer->pointer = minSubMatrix;
                    this->pointers.push(childPointer);
                } // END IF isLeaf
            } // END FOR children

            if (this->nextPointer->usedBy == 0){
                delete(this->nextPointer);
                this->nextPointer = nullptr;
            }

            if (!this->pointers.empty() && (int)this->pointers.top()->pointer->maxValue < this->maxValue) { // The next node has a value less than the current one
                hasNext = false;
                return;
            }
        }

        // Select the next MBR to be processed
        if (!this->pointers.empty()) {
            nextEntry = this->pointers.top()->id;
            nextPointer = this->pointers.top()->pointer;
            delete(this->pointers.top());
            this->pointers.pop();
            hasNext = true;
        } else {
            // No more MBR in the list, finish the process
            hasNext = false;
        }
    }


    //**********************************************************************//
    //******************* SETTERS AND GETTERS ******************************//
    //**********************************************************************//
    void GetObjectsQueryMaxMBR::clear() {
        GetObjectsQueryBasic::clear();

        // Remove submatrix
        if (this->nextPointer != nullptr && this->nextPointer->usedBy == 0) {
            delete(this->nextPointer);
        }

        while (!this->pointers.empty()) {
            this->pointers.top()->pointer->usedBy--;
            if (this->pointers.top()->pointer->usedBy == 0) {
                delete(this->pointers.top()->pointer);
            }

            delete(this->pointers.top());
            this->pointers.pop();
        }

        // Clean params
        this->maxValue = INT_MIN;
        this->nextPointer = nullptr;
    }

} // END NAMESPACE k2raster_rtree_static