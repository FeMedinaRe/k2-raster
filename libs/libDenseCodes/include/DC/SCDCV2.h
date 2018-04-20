/*
 * Created by Fernando Silva on 11/03/17.
 *
 * Copyright (C) 2016-current-year, Fernando Silva, all rights reserved.
 *
 * Author's contact: Fernando Silva  <fernando.silva@udc.es>
 * Databases Lab, University of A Coruña. Campus de Elviña s/n. Spain
 *
 *  C++ implementation of (s,c)-Dense Code. 
 *  N. Brisaboa, A. Fariña, G. Navarro, and J. Paramá. Lightweight natural language text compression. Information Retrieval. Information Retrieval(10), pp. 1-33, 2007. (online) (doi:10.1007/s10791-006-9001-9)
 *  
 *  For more information http://vios.dc.fi.udc.es/codes/semistatic.html
 *  Original code http://vios.dc.fi.udc.es/codes/files/scdc.tar.gz
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


#ifndef INCLUDE_DC_SCDC_V2_H_
#define INCLUDE_DC_SCDC_V2_H_

#include <sys/types.h>
#include <fstream>

namespace densecode_static {

class SCDCV2 {

uint const MAX_SC = 256;
uint const NUM_BITS_MOVE = 8;

public:
    /****** Constructors ******/
    SCDCV2(uint* data, ulong size, uint s=0);
    virtual ~SCDCV2();

    /****** Encode ******/
    uint encodeNumber(uint number, ulong &codeword);

    /****** Decode ******/
//    uint * decode();
    uint decodeNext();
    void resetNext();

    /****** Size ******/
    virtual size_t getTotalSize() const;

    /****** FILE ******/
    void save(std::ofstream &fp) const;
    static SCDCV2 *load(std::ifstream &in);

protected:
    /*** Params ***/
    ulong *baseTable;
    ushort s;
    ushort c;

    /*** Params data ***/
//    ulong size;
    unsigned char *dataEncoded;
    ulong totalLength;

    /*** Params aux access ***/
    ulong lastPostion=0;

    /****** Constructors ******/
    SCDCV2(){};

    /****** Aux functions ******/
    void encodeData(uint *data, uint size);
    void initBaseTable();
//    void encodeData2(uint *data, uint size);

    /*** Compute S and C ***/
    uint findBestS_Secuencial(uint* data, uint size);
    ulong computeSizeS(uint s, uint c, uint *freqAcum, ulong size, ulong sizeBest);

    };

} // namespace densecode_static



#endif  // INCLUDE_DC_SCDC_H_
