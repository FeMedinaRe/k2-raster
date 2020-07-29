/*
 * Created by Fernando Silva on 24/05/16.
 *
 * Copyright (C) 2016-current-year, Fernando Silva, all rights reserved.
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

#include <sys/types.h>
#include <queue>
#include <k2raster/compressLeaves/k2-raster-CompressLeaves.h>

namespace k2raster_static {

    // Tmp attributes
    utils_encoder_static::HashTableWord *hash;

    //**********************************************************************//
    //********************* CONSTRUCTORS ***********************************//
    //**********************************************************************//

    K2RasterCompressLeaves::K2RasterCompressLeaves() {
        this->plainValues = nullptr;
    };

    K2RasterCompressLeaves::K2RasterCompressLeaves(uint sizeX, uint sizeY, int **data, uint K1, uint K2, uint levelK1, uint plainLevels) :
            K2RasterPlain(sizeX, sizeY, data, K1, K2, levelK1, plainLevels) {

        if (this->levelK1 + this->levelK2 < this->nlevels ) {
            // Transform uints into char values. This is because we need to manage values ​​as words,
            // making it easier to use the hash table.
            unsigned char * plainWords = (unsigned char *)this->plainValues;

            // The word consists of the fusion of all values ​​of the submatrix
            uint submatrixSize = this->divKLevels[this->levelK1 + this->levelK2 - 1] * this->divKLevels[this->levelK1 + this->levelK2 - 1];
            this->lenWord = sizeof(uint) * submatrixSize;
            uint numOfSubmatrices = this->numOfPlainValues / submatrixSize;
//            printf("Number of Values- > %u, Number of Submatrices -> %u (elements %u), wordLen -> %u\n", this->numOfPlainValues, numOfSubmatrices, submatrixSize, this->lenWord);

            /********************** Beginning of the first pass **********************/

            // FIST STEP - Calculate frequency of each "word"
            unsigned char * word;
            hash = new utils_encoder_static::HashTableWord(42049060);
            uint *positionInTH = (uint*) malloc (numOfSubmatrices * sizeof(uint));
            uint j;
            ulong addrInTH;
            this->numOfWords = 0;
            for (uint i = 0; i < numOfSubmatrices; i ++) {
                word = &plainWords[i * this->lenWord];

//                if (i < 100) {
//                    printf("addrInTH -> %u \t", addrInTH);
//                }
                j = hash->searchElement(word, this->lenWord, &addrInTH);
//                if (i < 100) {
//                    printf("addrInTH -> %u (j -> %u)\n", addrInTH, j);
//                }

                if (j == this->numOfWords) {
                    hash->insertElement(word, this->lenWord, &addrInTH);
                    positionInTH[this->numOfWords++] = addrInTH;
                }
                hash->addWeight(addrInTH, 1);
            }

//            printf("Found %u words in leaves values\n", this->numOfWords);

            // SECOND STEP - Sort vocabulary by frequency
            std::qsort(positionInTH, this->numOfWords, sizeof(uint), compareFreqListDesc);

            // THIRD STEP - Assign codewords to vocabulary
            for (uint i = 0; i < this->numOfWords; i++) {
                hash->setCodeword(positionInTH[i], i);
            }

            // FOURTH STEP - Copy words
            this->words = (unsigned char *)malloc(this->numOfWords * this->lenWord * sizeof(unsigned char));
            this->totalLengthOfWords = 0;
            for (uint w = 0; w < this->numOfWords; w++) {
                for (uint l = 0; l < this->lenWord; l++) {
                    this->words[this->totalLengthOfWords++] = hash->getWord(positionInTH[w])[l];
                }
            }
            free(positionInTH);

            /********************** Beginning of the second pass **********************/

            // FIFTH STEP - Compact leaves
            uint *compactedValues = (uint*)malloc(numOfSubmatrices * sizeof(uint));
            for (uint i = 0; i < numOfSubmatrices; i++) {
                word = &plainWords[i * this->lenWord];
                j = hash->searchElement(word, this->lenWord, &addrInTH);
                compactedValues[i] = hash->getCodeword(addrInTH);
            }

            // SIXTH STEP - Encode wit DAC
            free(this->plainValues);
            this->plainValues = nullptr;
            delete(hash);
            this->valuesLeaves = new cds_static::DAC(compactedValues, numOfSubmatrices, 3);
            free(compactedValues);
        }
    }

    K2RasterCompressLeaves::~K2RasterCompressLeaves() {
        if (this->words != nullptr) {
            free(this->words);
        }
        if (this->valuesLeaves) {
            delete (valuesLeaves);
        }
    }

    //**********************************************************************//
    //********************* CELL FUNCTIONS *********************************//
    //**********************************************************************//

    int K2RasterCompressLeaves::getPlainCell(uint father, uint row, uint col, uint divK, uint lastValueMax) const {
        row = row % divK;
        col = col % divK;
        uint  arrayPosition =  row * divK + col;
        uint compressedValue = this->valuesLeaves->access(father);
        unsigned char * word = &this->words[compressedValue * this->lenWord];
        uint *values = (uint*)word;

        return lastValueMax - values[arrayPosition];
    }

    ulong K2RasterCompressLeaves::getPlainWindow(uint xini, uint yini, uint xend, uint yend,
                                                 uint offsetValues, uint divK, uint maxValue,
                                                 uint *positions, ulong totalCells, Query query) const {

        // Add all valid positions to final result
        uint value, arrayPosition;
        uint compressedValue = this->valuesLeaves->access(offsetValues);
        unsigned char * word = &this->words[compressedValue * this->lenWord];
        uint *values = (uint*)word;
        for (uint x = xini; x <= xend; x++){
            arrayPosition = (x % divK) * divK + (yini % divK);
            for (uint y= yini; y <= yend; y++){
                value = maxValue - values[arrayPosition];
                if (value >= query.valini && value <= query.valend) {
                    positions[totalCells++] = x * this->realSizeY + y;
                }
                arrayPosition++;
            }
        }
        return totalCells;
    }

    bool K2RasterCompressLeaves::checkPlainCells(uint xini, uint yini, uint xend, uint yend, uint offsetValues, uint divK,
                                        uint maxValue, Query query, bool allCells) const {

        // Add all valid positions to final result
        uint value, arrayPosition;
        uint compressedValue = this->valuesLeaves->access(offsetValues);
        unsigned char * word = &this->words[compressedValue * this->lenWord];
        uint *values = (uint*)word;
        for (uint x = xini; x <= xend; x++){
            arrayPosition = (x % divK) * divK + (yini % divK);
            for (uint y= yini; y <= yend; y++){
                value = maxValue - values[arrayPosition];
                if (value >= query.valini && value <= query.valend) {
                    if (!allCells) {
                        return true;
                    }
                } else {
                    if (allCells) {
                        return false;
                    }
                }
                arrayPosition++;
            }
        }
        return allCells;
    }

    ulong
    K2RasterCompressLeaves::getPlainCellValues(uint xini, uint yini, uint xend, uint yend, uint offsetValues, uint divK,
                                               uint maxValue, int *values, ulong totalCells, Query query) const {
        uint value, arrayPosition;
        uint cellPos;
        uint compressedValue = this->valuesLeaves->access(offsetValues);
        unsigned char * word = &this->words[compressedValue * this->lenWord];
        uint *wordValues = (uint*)word;
        for (uint x = xini; x <= xend; x++){
            arrayPosition = (x % divK) * divK + (yini % divK);
            cellPos = ((x - query.xini) * ((query.yend - query.yini) + 1));
            for (uint y= yini; y <= yend; y++){
                value = maxValue - wordValues[arrayPosition];
               // if (value >= query.valini && value <= query.valend) {
                    values[cellPos + (y - query.yini) ] = value;
                    totalCells++;
                //}
                arrayPosition++;
            }
        }
        return totalCells;
    }


    //**********************************************************************//
    //********************* FILE FUNCTIONS *********************************//
    //**********************************************************************//
    void K2RasterCompressLeaves::save(std::ofstream &of) const {
        K2Raster::save(of, K2RASTER_VOC);

        // Levelsk2
        saveValue(of, this->levelK2);

        // Save vocabulary
        saveValue(of, this->numOfWords);
        saveValue(of, this->totalLengthOfWords);
        saveValue(of, this->words, this->totalLengthOfWords);

        // Save DAC values of leaves
        this->valuesLeaves->save(of);
    }

    /*
     * Load Ktreap from file
     */
    K2RasterCompressLeaves *K2RasterCompressLeaves::load(std::ifstream &in) {
        K2RasterCompressLeaves *ktreap = nullptr;

        try {
            ktreap = new K2RasterCompressLeaves();
            K2Raster::loadBase(in, ktreap);

            // Load k2 value
            ktreap->levelK2 = loadValue<uint>(in);

            // Load vocabulary
            uint submatrixSize = ktreap->divKLevels[ktreap->levelK1 + ktreap->levelK2 - 1] * ktreap->divKLevels[ktreap->levelK1 + ktreap->levelK2 - 1];
            ktreap->lenWord = sizeof(uint) * submatrixSize;
            ktreap->numOfWords = loadValue<uint>(in);
            ktreap->totalLengthOfWords = loadValue<uint>(in);
            ktreap->words = loadValue<unsigned char>(in, ktreap->totalLengthOfWords);

            // Load DAC values of leaves
            ktreap->valuesLeaves = cds_static::DAC::load(in);
        } catch (...) {
            return nullptr;
        }

        return ktreap;
    }

    //**********************************************************************//
    //********************* SIZE FUNCTIONS *********************************//
    //**********************************************************************//
    size_t K2RasterCompressLeaves::getTotalSize() const {
        size_t totalSize = sizeof(uint)             // Type
                           + sizeof(uint)           // realSize
                           + (sizeof(uint) * 4)     // K1, K2, levelK1, levelK2
                           + (sizeof(uint) * 2);    // initialMaxValue and initialMinValue

        // k2-tree
        totalSize += this->bitmapK2Tree->getSize(); // bitMapsK2tree

        // DACs
        totalSize += this->getSizeDACs(this->listMaxValues, this->numOfMaxDACs); // Max DACs
        totalSize += this->getSizeDACs(this->listMinValues, this->numOfMinDACs); // Min DACs

        // Vocabulary
        totalSize += sizeof(uint)                                           // Number of words
                     + sizeof(uint)                                         // Total length of words
                     + sizeof(unsigned char) * this->totalLengthOfWords;    // Words

        // DAC values of leaves
        totalSize += this->valuesLeaves->getSize();
        return totalSize;
    }

    size_t K2RasterCompressLeaves::printSize() const {
        printf("RealSize \t-> %lu\n", sizeof(uint));
        printf("K1, K2, levelK1, levelK2 \t-> %lu\n", sizeof(uint) * 4);
        printf("InitialValues (Max and min) \t-> %lu\n", sizeof(uint) * 2);
        printf("BitMapK2Tree \t-> %lu\n", this->bitmapK2Tree->getSize());


        this->printSizeDACs(this->listMaxValues, this->numOfMaxDACs, false);
        this->printSizeDACs(this->listMinValues, this->numOfMinDACs, false);

        // Vocabulary
        printf("NumOfWords, TotalLengthOfWords \t-> %lu\n", sizeof(uint) * 2);
        printf("Words \t-> %lu\n", this->totalLengthOfWords * sizeof(unsigned char));

        printf("DAC values of leaves -> %lu\n", this->valuesLeaves->getSize());
        return this->getTotalSize();
    }


    //**********************************************************************//
    //******************* AUX HASH FUNCTIONS *******************************//
    //**********************************************************************//
    int compareFreqListDesc(const void *a, const void *b) {
        unsigned int *left,*right;
        left =  (uint*) a;
        right = (uint*) b;
//        printf("Compare %u with %u\n", *left, *right);
        if (hash->getWeight(*left) < hash->getWeight(*right))
            return 1;
        if (hash->getWeight(*left) > hash->getWeight(*right))
            return -1;

        return 0;
    }
}

