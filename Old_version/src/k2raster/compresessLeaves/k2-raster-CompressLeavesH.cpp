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
#include <k2raster/compressLeaves/k2-raster-CompressLeavesH.h>

namespace k2raster_static {

    // Tmp attributes
    utils_encoder_static::HashTableWord *hashH;

    //**********************************************************************//
    //********************* CONSTRUCTORS ***********************************//
    //**********************************************************************//

    K2RasterCompressLeavesH::K2RasterCompressLeavesH() {
        this->plainValues = nullptr;
    };

    // TODO optimize memory
    K2RasterCompressLeavesH::K2RasterCompressLeavesH(uint sizeX, uint sizeY, int **data, uint K1, uint K2, uint levelK1, uint plainLevels, float minFreqVoc) :
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
            hashH = new utils_encoder_static::HashTableWord(42049060);
            uint *positionInTH = (uint*) malloc (numOfSubmatrices * sizeof(uint));
            uint j;
            ulong addrInTH;
            this->numOfWords = 0;
//            printf("Adding submatrices to hash table....\n");
            for (uint i = 0; i < numOfSubmatrices; i ++) {
                word = &plainWords[i * this->lenWord];
                j = hashH->searchElement(word, this->lenWord, &addrInTH);
                if (j == this->numOfWords) {
                    hashH->insertElement(word, this->lenWord, &addrInTH);
                    positionInTH[this->numOfWords++] = addrInTH;
                }
                hashH->addWeight(addrInTH, 1);
            }
//            printf("Found %u words in leaves values\n", this->numOfWords);
//            positionInTH = (uint*)realloc(positionInTH, this->numOfWords * sizeof(uint));

            // SECOND STEP - Sort vocabulary by frequency
            std::qsort(positionInTH, this->numOfWords, sizeof(uint), compareFreqListDescH);

            // THIRD STEP - Assign codewords to vocabulary (only when weight if bigger than minFreqVoc
            int codeWord = 0;
            // Recalcule minFreqVoc, from percent to value
            uint nminFreqVoc = (minFreqVoc * numOfSubmatrices / 100.0);
            for (uint i = 0; i < this->numOfWords; i++) {
                if (hashH->getWeight(positionInTH[i]) > nminFreqVoc) {
                    hashH->setCodeword(positionInTH[i], codeWord);
                    codeWord++;
                } else {
                    hashH->setCodeword(positionInTH[i], -1);
                }
            }
            this->numOfWords = codeWord; // Only store words with frequency bigger than minFreqVoc
//            printf("Found %u words with frequency bigger than %0.4f\% (%u)\n", this->numOfWords, minFreqVoc, nminFreqVoc);

            /********************** Beginning of the second pass **********************/

            // FOURTH STEP - Compact leaves
            uint *compactedValues = (uint*)malloc(numOfSubmatrices * sizeof(uint));
            uint *noCompactedValues = (uint *)malloc((numOfSubmatrices - this->numOfWords) * submatrixSize  * sizeof(uint));
            uint *tmpBitmapVoc = (uint*)malloc(cds_static::uint_len(numOfSubmatrices, 1) * sizeof(uint));
            uint numOfCompacted = 0, numOfNoCompacted = 0;
            for (uint i = 0; i < numOfSubmatrices; i++) {
                word = &plainWords[i * this->lenWord];
                j = hashH->searchElement(word, this->lenWord, &addrInTH);
                codeWord = hashH->getCodeword(addrInTH);
                if (codeWord != -1) {
                    // Use vocabulary
                    compactedValues[numOfCompacted++] = codeWord;
                    cds_static::bit_set(tmpBitmapVoc, i);
                } else {
                    // Store in basic form
                    for (uint c = 0; c < submatrixSize; c++) {
                        noCompactedValues[numOfNoCompacted++] = this->plainValues[i * submatrixSize + c];
                    }
                    cds_static::bitclean(tmpBitmapVoc, i);
                }
            }
//            compactedValues = (uint*)realloc(compactedValues, numOfCompacted* sizeof(uint));
//            noCompactedValues = (uint*)realloc(noCompactedValues, numOfNoCompacted * sizeof(uint));
            free(this->plainValues);
            this->plainValues = nullptr;

            // SIXTH STEP - Copy words
//            printf("Coping %u words ..........\n", this->numOfWords);
            this->words = (unsigned char *)malloc(this->numOfWords * this->lenWord * sizeof(unsigned char));
            this->totalLengthOfWords = 0;
            for (uint w = 0; w < this->numOfWords; w++) {
                for (uint l = 0; l < this->lenWord; l++) {
                    this->words[this->totalLengthOfWords++] = hashH->getWord(positionInTH[w])[l];
                }
            }
            free(positionInTH);
            delete(hashH);

            // FIFTH STEP - Encode wit DAC and bitsequence
//            printf("Creating bitmapVocLeaves with %u submatrices\n", numOfSubmatrices);
            cds_static::BitSequenceBuilder *builder = new cds_static::BitSequenceBuilderRG(20);
            this->bitmapVocLeaves = builder->build(tmpBitmapVoc, numOfSubmatrices);
            free(tmpBitmapVoc);
            delete(builder);
//            printf("Creating DAC for %u compacted values\n", numOfCompacted);
            this->vocLeaves = new cds_static::DAC(compactedValues, numOfCompacted, 3);
            free(compactedValues);
//            printf("Creating DAC for %u non compacted values\n", numOfNoCompacted);
            this->valuesLeaves = new cds_static::DAC(noCompactedValues, numOfNoCompacted, 3);
            free(noCompactedValues);
        }
    }

    K2RasterCompressLeavesH::~K2RasterCompressLeavesH() {
        if (this->words != nullptr) {
            free(this->words);
        }
        if (this->vocLeaves != nullptr) {
            delete (vocLeaves);
        }
        if (this->valuesLeaves != nullptr) {
            delete(valuesLeaves);
        }
        if (this->bitmapVocLeaves != nullptr) {
            delete(bitmapVocLeaves);
        }
    }

    //**********************************************************************//
    //********************* CELL FUNCTIONS *********************************//
    //**********************************************************************//

    int K2RasterCompressLeavesH::getPlainCell(uint father, uint row, uint col, uint divK, uint lastValueMax) const {
        row = row % divK;
        col = col % divK;
        uint  arrayPosition =  row * divK + col;

        // Check in bitmap if submatrix if in vocabulary or not
        if (this->bitmapVocLeaves->access(father)) {
            uint pos = this->bitmapVocLeaves->rank1(father) - 1;
            uint compressedValue = this->vocLeaves->access(pos);
            unsigned char * word = &this->words[compressedValue * this->lenWord];
            uint *values = (uint*)word;
            return lastValueMax - values[arrayPosition];
        } else {
            uint pos = (this->bitmapVocLeaves->rank0(father) - 1) * divK * divK;
            return lastValueMax - this->valuesLeaves->access(pos + arrayPosition);
        }
    }

    ulong K2RasterCompressLeavesH::getPlainWindow(uint xini, uint yini, uint xend, uint yend,
                                                  uint offsetValues, uint divK, uint maxValue,
                                                  uint *positions, ulong totalCells, Query query) const {
        // Add all valid positions to final result
        uint value, arrayPosition;

        // Check in bitmap if submatrix if in vocabulary or not
        if (this->bitmapVocLeaves->access(offsetValues)) {
            uint pos = this->bitmapVocLeaves->rank1(offsetValues) - 1;
            uint compressedValue = this->vocLeaves->access(pos);
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
        } else {
            uint pos = (this->bitmapVocLeaves->rank0(offsetValues) - 1) * divK * divK;
            for (uint x = xini; x <= xend; x++){
                arrayPosition = (x % divK) * divK + (yini % divK);
                // First position in y
                value = maxValue - this->valuesLeaves->access(pos + arrayPosition);
//                printf("value\n");
                if (value >= query.valini && value <= query.valend) {
                    positions[totalCells++] = x * this->realSizeY + yini;
                }
                for (uint y= yini+1; y <= yend; y++){
                    value = maxValue - this->valuesLeaves->next();
                    if (value >= query.valini && value <= query.valend) {
                        positions[totalCells++] = x * this->realSizeY + y;
                    }
                }
            }
        }
        return totalCells;
    }

    bool K2RasterCompressLeavesH::checkPlainCells(uint xini, uint yini, uint xend, uint yend, uint offsetValues, uint divK,
                                        uint maxValue, Query query, bool allCells) const {

        // Add all valid positions to final result
        uint value, arrayPosition;

        if (this->bitmapVocLeaves->access(offsetValues)) {
            uint pos = this->bitmapVocLeaves->rank1(offsetValues) - 1;
            uint compressedValue = this->vocLeaves->access(pos);
            unsigned char *word = &this->words[compressedValue * this->lenWord];
            uint *values = (uint *) word;
            for (uint x = xini; x <= xend; x++) {
                arrayPosition = (x % divK) * divK + (yini % divK);
                for (uint y = yini; y <= yend; y++) {
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
        } else {
            uint pos = (this->bitmapVocLeaves->rank0(offsetValues) - 1) * divK * divK;
            for (uint x = xini; x <= xend; x++){
                arrayPosition = (x % divK) * divK + (yini % divK);
                // First position in y
                value = maxValue - this->valuesLeaves->access(pos + arrayPosition);
                if (value >= query.valini && value <= query.valend) {
                    if (!allCells) {
                        return true;
                    }
                } else {
                    if (allCells) {
                        return false;
                    }
                }
                for (uint y= yini+1; y <= yend; y++){
                    value = maxValue - this->valuesLeaves->next();
                    if (value >= query.valini && value <= query.valend) {
                        if (!allCells) {
                            return true;
                        }
                    } else {
                        if (allCells) {
                            return false;
                        }
                    }
                }
            }
        }
        return allCells;
    }

    ulong K2RasterCompressLeavesH::getPlainCellValues(uint xini, uint yini, uint xend, uint yend, uint offsetValues,
                                                      uint divK,
                                                      uint maxValue, int *values, ulong totalCells, Query query) const {
        uint value, arrayPosition;
        uint cellPos;

        if (this->bitmapVocLeaves->access(offsetValues)) {
            uint pos = this->bitmapVocLeaves->rank1(offsetValues) - 1;
            uint compressedValue = this->vocLeaves->access(pos);
            unsigned char *word = &this->words[compressedValue * this->lenWord];
            uint *wordValues = (uint *) word;
            for (uint x = xini; x <= xend; x++) {
                arrayPosition = (x % divK) * divK + (yini % divK);
                cellPos = ((x - query.xini) * ((query.yend - query.yini) + 1));
                for (uint y = yini; y <= yend; y++) {
                    value = maxValue - wordValues[arrayPosition];
                    // if (value >= query.valini && value <= query.valend) {
                    values[cellPos + (y - query.yini)] = value;
                    totalCells++;
                    //}
                    arrayPosition++;
                }
            }
        } else {
            uint pos = (this->bitmapVocLeaves->rank0(offsetValues) - 1) * divK * divK;
            for (uint x = xini; x <= xend; x++){
                arrayPosition = (x % divK) * divK + (yini % divK);
                cellPos = ((x - query.xini) * ((query.yend - query.yini) + 1));
                value = maxValue - this->valuesLeaves->access(pos + arrayPosition);
                values[cellPos + (yini - query.yini)] = value;
                totalCells++;
                for (uint y= yini+1; y <= yend; y++){
                    value = maxValue - this->valuesLeaves->next();
                    values[cellPos + (y - query.yini)] = value;
                    totalCells++;
                }
            }
        }
        return totalCells;
    }


    //**********************************************************************//
    //********************* FILE FUNCTIONS *********************************//
    //**********************************************************************//
    void K2RasterCompressLeavesH::save(std::ofstream &of) const {
        K2Raster::save(of, K2RASTER_VOC_HYBRID);

        // Levelsk2
        saveValue(of, this->levelK2);

        // Save vocabulary
        saveValue(of, this->numOfWords);
        saveValue(of, this->totalLengthOfWords);
        saveValue(of, this->words, this->totalLengthOfWords);

        // Save DAC values of leaves
        this->vocLeaves->save(of);
        this->valuesLeaves->save(of);

        // Save bitmap vocabulary of leaves
        this->bitmapVocLeaves->save(of);
    }

    /*
     * Load Ktreap from file
     */
    K2RasterCompressLeavesH *K2RasterCompressLeavesH::load(std::ifstream &in) {
        K2RasterCompressLeavesH *ktreap = nullptr;

        try {
            ktreap = new K2RasterCompressLeavesH();
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
            ktreap->vocLeaves = cds_static::DAC::load(in);
            ktreap->valuesLeaves = cds_static::DAC::load(in);

            // Load bitmap of vocabulary of leaves
            ktreap->bitmapVocLeaves = cds_static::BitSequence::load(in);
        } catch (...) {
            return nullptr;
        }

        return ktreap;
    }

    //**********************************************************************//
    //********************* SIZE FUNCTIONS *********************************//
    //**********************************************************************//
    size_t K2RasterCompressLeavesH::getTotalSize() const {
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
        totalSize += this->vocLeaves->getSize();
        totalSize += this->valuesLeaves->getSize();

        // Bitmap vocabulary of leaves
        totalSize += this->bitmapVocLeaves->getSize();
        return totalSize;
    }

    size_t K2RasterCompressLeavesH::printSize() const {
        printf("RealSize \t-> %lu\n", sizeof(uint));
        printf("K1, K2, levelK1, levelK2 \t-> %lu\n", sizeof(uint) * 4);
        printf("InitialValues (Max and min) \t-> %lu\n", sizeof(uint) * 2);
        printf("BitMapK2Tree \t-> %lu\n", this->bitmapK2Tree->getSize());


        this->printSizeDACs(this->listMaxValues, this->numOfMaxDACs, false);
        this->printSizeDACs(this->listMinValues, this->numOfMinDACs, false);

        // Vocabulary
        printf("NumOfWords, TotalLengthOfWords \t-> %lu\n", sizeof(uint) * 2);
        printf("Words \t-> %lu\n", this->totalLengthOfWords * sizeof(unsigned char));

        printf("DAC vocabulary of leaves -> %lu\n", this->vocLeaves->getSize());
        printf("DAC values of leaves -> %lu\n", this->vocLeaves->getSize());
        printf("Bitmap vocabulary of leaves -> %lu\n", this->bitmapVocLeaves->getSize());
        return this->getTotalSize();
    }


    //**********************************************************************//
    //******************* AUX HASH FUNCTIONS *******************************//
    //**********************************************************************//
    int compareFreqListDescH(const void *a, const void *b) {
        unsigned int *left,*right;
        left =  (uint*) a;
        right = (uint*) b;
//        printf("Compare %u with %u\n", *left, *right);
        if (hashH->getWeight(*left) < hashH->getWeight(*right))
            return 1;
        if (hashH->getWeight(*left) > hashH->getWeight(*right))
            return -1;

        return 0;
    }
}

