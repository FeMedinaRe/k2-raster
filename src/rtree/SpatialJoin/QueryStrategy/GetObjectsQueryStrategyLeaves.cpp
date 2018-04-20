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

#include <rtree/SpatialJoin/QueryStrategy/GetObjectsQueryStrategyLeaves.h>

namespace k2raster_rtree_static {

    //**********************************************************************//
    //*************************** ENTRIES **********************************//
    //**********************************************************************//

    /****** Process MBR ******/
    void GetObjectsQueryStrategyLeaves::getNextEntry(const IEntry &entry, id_type &nextEntry, bool &hasNext) {
        IShape *ps;
        entry.getShape(&ps);
        Region *pr = dynamic_cast<Region *>(ps);

        if (this->first) {
            // Calculate r-tree params (only once)
            this->offsetX = pr->m_pLow[0] * -1;
            this->offsetY = pr->m_pLow[1] * -1;

            this->treeSizeX = pr->m_pHigh[0] + this->offsetX;
            this->treeSizeY = pr->m_pHigh[1] + this->offsetY;

            this->first = false;
            this->rootNode = this->raster->createRootMinSubMatrix(this->valini, this->valend);
        }

        uint xini, xend, yini, yend;
        this->calculatePositions(pr->m_pLow, pr->m_pHigh, xini, xend, yini, yend);
        delete ps;

        const INode *n = dynamic_cast<const INode *>(&entry);
        if (n != 0) {
            if (n->isLeaf()) {
                // Search children
                this->processChildren(n, this->rootNode);
            } else {
                // Internal node
                // Traverse the MBR and add its children to be processed
                for (uint32_t cChild = 0; cChild < n->getChildrenCount(); cChild++) {
                    ids.push(n->getChildIdentifier(cChild));
                }
            } // ENF IF leaf
        }

        // Select the next MBR to be processed
        if (!ids.empty()) {
            nextEntry = ids.top();
            ids.pop();
            hasNext = true;
        } else {
            // No more MBR in the list, finish the process
            hasNext = false;
            free(this->rootNode);
        }
    }

    //**********************************************************************//
    //******************* SETTERS AND GETTERS ******************************//
    //**********************************************************************//
    void GetObjectsQueryStrategyLeaves::clear() {
        uint *a;
        std::list<tuple<id_type, uint *, ulong> *>::iterator it = this->result.begin();
        while (it != this->result.end()) {
            a = std::get<1>(*(*it));
            free(a);
            delete (*it);
            this->result.erase(it++);
        }

        this->result.clear();
        this->allMBRCells = 0;
        this->first = true;
    }
}
