/*  
 * Created by Fernando Silva on 17/01/17.
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

#include <rtree/SpatialJoin/Visitor/GetObjectsVisitorCells.h>

namespace k2raster_rtree_static {
    //**********************************************************************//
    //************************** Constructor *******************************//
    //**********************************************************************//
    GetObjectsVisitorCells::GetObjectsVisitorCells(double *m_pLow, double *m_pHigh, uint rasterSizeX,
                                                   uint rasterSizeY) {
        this->calculatePositionParams(m_pLow, m_pHigh, rasterSizeX, rasterSizeY);
    }

    //**********************************************************************//
    //*************************** VISIT ************************************//
    //**********************************************************************//
    void GetObjectsVisitorCells::visitNode(const INode &n) {
        if (n.isLeaf()) {
            try {
                std::tuple<uint *, ulong> *data = this->result.at(n.getIdentifier());
                std::get<1>(*data) = std::get<1>(*data) + 1;
                std::get<0>(*data)[std::get<1>(*data)] = this->point;
//                printf("MBR %lu has %lu cells\n", n.getIdentifier(), std::get<1>(*data));
            } catch (const std::out_of_range &oor) {
                IShape *ps;
                n.getShape(&ps);
                Region *pr = dynamic_cast<Region *>(ps);

                uint xini, xend, yini, yend;
                this->calculatePositions(pr->m_pLow, pr->m_pHigh, xini, xend, yini, yend);
                delete (ps);
                uint *positions = (uint *) malloc((xend - xini + 1) * (yend - yini + 1) * sizeof(uint));
                std::tuple<uint *, ulong> *data = new std::tuple<uint *, ulong>(positions, 0);
                this->result[n.getIdentifier()] = data;
//                printf("Added MBR %lu\n", n.getIdentifier());
            }
        }
    }

    //**********************************************************************//
    //*************************** POSITIONS ********************************//
    //**********************************************************************//
    void GetObjectsVisitorCells::calculatePositionParams(double *m_pLow, double *m_pHigh, uint rasterSizeX,
                                                         uint rasterSizeY) {
        this->offsetX = m_pLow[0] * -1;
        this->offsetY = m_pLow[1] * -1;

        this->treeSizeX = m_pHigh[0] + this->offsetX;
        this->treeSizeY = m_pHigh[1] + this->offsetY;

        this->rasterSizeX = rasterSizeX;
        this->rasterSizeY = rasterSizeY;
    }

    void GetObjectsVisitorCells::calculatePositions(double *m_pLow, double *m_pHigh,
                                                    uint &xini, uint &xend, uint &yini, uint &yend) {
        xini = ((m_pLow[0] + this->offsetX) / this->treeSizeX) * (this->rasterSizeX - 1);
        xend = ((m_pHigh[0] + this->offsetX) / this->treeSizeX) * (this->rasterSizeX - 1);
        yini = ((m_pLow[1] + this->offsetY) / this->treeSizeY) * (this->rasterSizeY - 1);
        yend = ((m_pHigh[1] + this->offsetY) / this->treeSizeY) * (this->rasterSizeY - 1);
    }


    //**********************************************************************//
    //***************************** RESULT *********************************//
    //**********************************************************************//
    std::map<id_type, std::tuple<uint *, ulong> *> GetObjectsVisitorCells::getResult() const {
        return this->result;
    }

    ulong GetObjectsVisitorCells::getAllMBRCells() const {
        return this->allMBRCells;
    }

    //**********************************************************************//
    //*********************** GETTERs AND SETTERs **************************//
    //**********************************************************************//
    void GetObjectsVisitorCells::clear() {
        uint *a;
        std::map<id_type, std::tuple<uint *, ulong> *>::iterator it = this->result.begin();
        while (it != this->result.end()) {
            a = std::get<0>(*(*it).second);
            free(a);
            delete (it->second);
            this->result.erase(it++);
        }

        this->result.clear();
        this->allMBRCells = 0;
    }
}

