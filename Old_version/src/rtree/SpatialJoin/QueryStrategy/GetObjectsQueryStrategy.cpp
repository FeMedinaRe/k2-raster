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

#include <rtree/SpatialJoin/QueryStrategy/GetObjectsQueryStrategy.h>


namespace k2raster_rtree_static {

    //**********************************************************************//
    //******************* SETTERS AND GETTERS ******************************//
    //**********************************************************************//
    void GetObjectsQueryStrategy::setRaster(K2Raster *raster) {
        this->raster = raster;
        this->sizeX = raster->getDimSize(1);
        this->sizeY = raster->getDimSize(2);
    }

    void GetObjectsQueryStrategy::clear() {
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

//        this->numberOfWhites = 0;
//        this->numberOfBlacks = 0;
    }

    void GetObjectsQueryStrategy::setRangeOfValues(uint valini, uint valend) {
        this->valini = valini;
        this->valend = valend;
    }

    void GetObjectsQueryStrategy::setAllCells(bool allCells) {
        this->allCell = allCells;
    }

    //**********************************************************************//
    //***************************** RESULTS ********************************//
    //**********************************************************************//
    list<tuple<id_type, uint *, ulong> *> GetObjectsQueryStrategy::getResult() const {
        return this->result;
    }

    ulong GetObjectsQueryStrategy::getAllMBRCells() const {
        return this->allMBRCells;
    }

    //**********************************************************************//
    //*************************** POSITIONS ********************************//
    //**********************************************************************//
    void GetObjectsQueryStrategy::calculatePositionParams(double *m_pLow, double *m_pHigh) {
        this->offsetX = m_pLow[0] * -1;
        this->offsetY = m_pLow[1] * -1;

        this->treeSizeX = m_pHigh[0] + this->offsetX;
        this->treeSizeY = m_pHigh[1] + this->offsetY;
    }

    void GetObjectsQueryStrategy::calculatePositions(double *m_pLow, double *m_pHigh,
                                                     uint &xini, uint &xend, uint &yini, uint &yend) {
        xini = ((m_pLow[0] + this->offsetX) / this->treeSizeX) * (this->sizeX - 1);
        xend = ((m_pHigh[0] + this->offsetX) / this->treeSizeX) * (this->sizeX - 1);
        yini = ((m_pLow[1] + this->offsetY) / this->treeSizeY) * (this->sizeY - 1);
        yend = ((m_pHigh[1] + this->offsetY) / this->treeSizeY) * (this->sizeY - 1);

    }

    //**********************************************************************//
    //************************* MBR PROCESSING *****************************//
    //**********************************************************************//
    void GetObjectsQueryStrategy::processChildren(const INode *n, MinBasicSubMatrix *minMatrix) {
        for (uint32_t cChild = 0; cChild < n->getChildrenCount(); cChild++) {
            // Search individual object
            IShape *ps2;
            Region re2;
            n->getChildShape(cChild, &ps2);
//            re2 = dynamic_cast<Region *>(ps2);
            ps2->getMBR(re2);
            uint xinic, xendc, yinic, yendc;
            this->calculatePositions(re2.m_pLow, re2.m_pHigh, xinic, xendc, yinic, yendc);
            delete(ps2);;

            // Alloc array os positions
            ulong totalCells = 0;
            uint *positionsFinal = nullptr;
            ulong positionsSize = ((xendc - xinic + 1) * (yendc - yinic + 1));
            uint *positions = (uint *) malloc( positionsSize * sizeof(uint));

            // Search cells
            if (minMatrix->color == SubMatrixColor::COLOR_SUBMATRIX_BLACK) {
                // The MBR is valid (and all objects within of the MBR)
                for (uint x = xinic; x <= xendc; x++) {
                    for (uint y = yinic; y <= yendc; y++) {
                        positions[totalCells++] = x * this->sizeY + y;
                    } // END FOR y
                } // ENF FOR x
                positionsFinal = positions;
            } else {
                totalCells = raster->getCellsByValue(xinic, xendc, yinic, yendc, this->valini, this->valend,
                                                     minMatrix, positions, this->allCell);

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
                    new tuple<id_type, uint *, ulong>(n->getIdentifier(), positionsFinal,
                                                      totalCells));
            this->allMBRCells += totalCells;

        } // END FOR children
    }
}
