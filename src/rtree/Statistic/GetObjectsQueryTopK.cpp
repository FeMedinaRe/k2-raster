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

#include <rtree/Statistic/QueryStrategy/GetObjectsQueryTopK.h>

namespace k2raster_rtree_static {

    //**********************************************************************//
    //*************************** ENTRIES **********************************//
    //**********************************************************************//

    /****** Process MBR ******/
    template <class T>
    void GetObjectsQueryTopK<T>::getNextEntry(const IEntry &entry, id_type &nextEntry, bool &hasNext) {
        BasicSubMatrixWithUsed *minSubMatrix;

        if (this->first) {
            // ----------------------- //
            // Process Root Node       //
            // ----------------------- //
            this->first = false;

            // Get object
            IShape *ps;
            entry.getShape(&ps);
            Region *pr = dynamic_cast<Region *>(ps);

            // First iteration of getNextEntry
            this->calculatePositionParams(pr->m_pLow, pr->m_pHigh);

            // Pointer to root node of the k2-raster
            this->nextPointer = this->raster->createRootBasicSubMatrixWithUsed(); // Root node
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
        } else {
            // ----------------------- //
            // Process node            //
            // ----------------------- //
            // Get the current node from queue
            this->nextPointer = this->pointers.top()->pointer;
            delete (this->pointers.top());
            this->pointers.pop();

        }
        this->nextPointer->usedBy--;

        // ----------------------- //
        // Check current node      //
        // ----------------------- //
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
                minSubMatrix->usedBy++;
                priorityPointerMBR *childPointer = new priorityPointerMBR;
                childPointer->positions = nullptr;
                childPointer->id = n->getChildIdentifier(cChild);

                if (n->isLeaf()) {
                    // Search maxValue/minValue of each child
                    int value = this->getMBRLeafValue(xini, xend, yini, yend, minSubMatrix);

                    if (this->withPositions) {
                        PositionsStruct *pos = new PositionsStruct;
                        pos->xini = xini;
                        pos->xend = xend;
                        pos->yini = yini;
                        pos->yend = yend;

                        childPointer->positions = pos;
                    } else {
                        if (minSubMatrix != this->nextPointer) {
                            delete (minSubMatrix);
                        }
                        minSubMatrix = nullptr;
                    }
                    childPointer->isRealValue = true;               // It is the real max/min value of the MBR
                    childPointer->value = value;
                } else {
                    childPointer->isRealValue = false;              // pointer.value is the possible max/min value of the MBR
                    childPointer->value = this->getMBRValue(minSubMatrix);
                } // END IF isLeaf
                childPointer->pointer = minSubMatrix;
                this->pointers.push(childPointer); // Add node to queue
            } // END FOR children

            if (this->nextPointer->usedBy == 0){
                delete(this->nextPointer);
                this->nextPointer = nullptr;
            }
        }

        this->setNextEntry(nextEntry, hasNext);
    }

    template <class T>
    void GetObjectsQueryTopK<T>::setNextEntry(id_type &nextEntry, bool &hasNext) {

        hasNext = false; // if hasNext = false; No more MBR in the list, finish the process
        priorityPointerMBR *ptr;
        while (!this->pointers.empty()) {
            ptr = this->pointers.top();
            if (this->pointers.top()->isRealValue) {

                // ----------------------- //
                // Add Result              //
                // ----------------------- //
                this->values[this->numResults] = ptr->value;
                this->MBRIds[this->numResults] = ptr->id;

                // Add positions if withPositions is true
                if (this->withPositions) {

                    // Allocate memory
                    ulong num = (ptr->positions->xend - ptr->positions->xini + 1) * (ptr->positions->yend - ptr->positions->yini + 1);
                    this->positions[this->numResults] = (uint*)malloc(num * sizeof(uint));

                    ptr->pointer->bitmapChildren = 0;
                    ulong  num2 = this->getCells(ptr);

                    if (num2 == 0) {
                        printf("Error getting position for MBR %lu and value %i\n", ptr->id, ptr->value);
                        exit(-1);
                    }


                    if (num != num2) {
                        uint* result = (uint*)realloc(this->positions[this->numResults], num2);
                        if (result != nullptr) {
                            this->positions[this->numResults] = result;
                        }
                    }
                    this->numOfPositions[this->numResults] = num2;

                    ptr->pointer->usedBy--;
                    if (ptr->pointer->usedBy == 0) {
                        delete (ptr->pointer);
                    }
                    delete(ptr->positions);
                }

                this->numResults++;
                delete(ptr);
                this->pointers.pop();

                if (this->numResults == k) {
                    // No more MBRs to collect
                    hasNext = false;
                    break;
                }
            } else {
                // ----------------------- //
                // Set next MBR                //
                // ----------------------- //
                nextEntry = ptr->id;
                hasNext = true;
                break;
            }
        }
    }


    //**********************************************************************//
    //******************* SETTERS AND GETTERS ******************************//
    //**********************************************************************//
    template <class T>
    void GetObjectsQueryTopK<T>::initialize() {
        GetObjectsQueryBasic::initialize();
        this->nextPointer = nullptr;
        this->numResults = 0;
        this->values = new int[this->k];
        this->MBRIds = new id_type[this->k];

        if (this->withPositions) {
            this->positions = new uint*[this->k];
            this->numOfPositions = new ulong[this->k];
        } else {
            this->positions = nullptr;
            this->numOfPositions = nullptr;
        }
    }

    template <class T>
    void GetObjectsQueryTopK<T>::clear() {
        GetObjectsQueryBasic::clear();

        // Remove submatrix
        if (this->nextPointer != nullptr && this->nextPointer->usedBy == 0) {
            delete(this->nextPointer);
        }

        priorityPointerMBR *ptr;
        while (!this->pointers.empty()) {
            ptr = this->pointers.top();
            if (ptr->pointer != nullptr) {
                ptr->pointer->usedBy--;
                if (ptr->pointer->usedBy == 0) {
                    delete (ptr->pointer);
                }
            }
            if (ptr->positions != nullptr) {
                delete(ptr->positions);
            }

            delete(this->pointers.top());
            this->pointers.pop();
        }

        // Clean params
        delete[](this->values);
        delete[](this->MBRIds);

        if (this->withPositions) {
            for (uint mbr = 0; mbr < this->numResults; mbr++) {
//                delete[](this->positions[mbr]);
                free(this->positions[mbr]);
            }
            delete[](this->positions);
            delete[](this->numOfPositions);
        }
    }


    template class GetObjectsQueryTopK<cmpPriorityPointerMBRMin>;
    template class GetObjectsQueryTopK<cmpPriorityPointerMBRMax>;
} // END NAMESPACE k2raster_rtree_static