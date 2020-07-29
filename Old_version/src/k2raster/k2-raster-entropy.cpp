/*
 * Created by Fernando Silva on 13/06/16.
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
#include <k2raster/k2-raster-entropy.h>

namespace k2raster_static {


    // TODO Remove this, Only it is declared in this way to be used in the sort function. We need a better solution.
    utils_encoder_static::HashTableWord *hashHEntropy;

    //**********************************************************************//
    //********************* CONSTRUCTORS ***********************************//
    //**********************************************************************//

    K2RasterEntropy::K2RasterEntropy() {
        this->plainValues = nullptr;
        this->bitmapVocLeaves = nullptr;
        this->vocLeaves = nullptr;
        this->valuesLeaves = nullptr;
    };

    // TODO optimize memory
    K2RasterEntropy::K2RasterEntropy(uint sizeX, uint sizeY, int **data, uint K1, uint K2, uint levelK1, uint plainLevels) :
            K2RasterPlain(sizeX, sizeY, data, K1, K2, levelK1, plainLevels) {

        if (this->levelK1 + this->levelK2 < this->nlevels ) {
            // Transforms uints into char values. This is because we need to manage values ​​as words (strings),
            // making it easier to use the hash table.
            // The word consists of the fusion of all values ​​of the submatrix
            unsigned char * plainWords = (unsigned char *)this->plainValues;

            uint submatrixSize = this->divKLevels[this->levelK1 + this->levelK2 - 1]
                                 * this->divKLevels[this->levelK1 + this->levelK2 - 1];     // Size of each submatrix
            this->lenWord = sizeof(uint) * submatrixSize;                                   // Size, in bytes, of each submatrix
            uint numOfSubmatrices = this->numOfPlainValues / submatrixSize;                 // Total number of submatrices
//            printf("Number of Values- > %u, Number of Submatrices -> %u (elements %u), wordLen -> %u\n", this->numOfPlainValues, numOfSubmatrices, submatrixSize, this->lenWord);


            /******************************************************************************/
            // FIST STEP - Calculate frequency of each "word" and each different "value"
            /******************************************************************************/
            unsigned char * word, *wordValue;
            uint numOfValues = 0, minValue=MAXINT, maxValue=0, value, j;
            ulong addrInTH;

            // Creates hash tables for words and values
            //TODO Calculate the initial memory of the hash table dynamically
            hashHEntropy = new utils_encoder_static::HashTableWord(numOfSubmatrices);
            utils_encoder_static::HashTableWord *hashValues = new utils_encoder_static::HashTableWord(numOfSubmatrices);
            ulong *positionInTH = (ulong*) malloc (numOfSubmatrices * sizeof(ulong)); // TODO change array by vector (or another std type)
//            std::vector<ulong> positionInTH;

            this->numOfWords = 0;

            // Adds submatrices and their values into the hash tables
            //printf("Adding submatrices to hash table....\n");
            for (uint i = 0; i < numOfSubmatrices; ++i) {
                word = &plainWords[i * this->lenWord];                              // Reads the "word"
                j = hashHEntropy->searchElement(word, this->lenWord, &addrInTH);    // Check if it exists
                if (j == this->numOfWords) {
                    // If the word does not exist, it is inserted into the hash table
                    hashHEntropy->insertElement(word, this->lenWord, &addrInTH);
                    positionInTH[this->numOfWords++] = addrInTH;
//                    positionInTH.push_back(addrInTH);
//                    this->numOfWords++;
                }
                // Goes up the frequncy
                hashHEntropy->addWeight(addrInTH, 1);

                // Check all values of the submatrix
                for (uint c = 0; c < submatrixSize; ++c) {
                    value = this->plainValues[i * submatrixSize + c];                   // Reads the value
                    wordValue = (unsigned char*)&value;                                 // It is converted in a string
                    j = hashValues->searchElement(wordValue, sizeof(uint), &addrInTH);   // Check if it exists
                    if (j == numOfValues) {
                        // If the value does not exist, it is inserted into the hash table
                        hashValues->insertElement(wordValue, sizeof(uint), &addrInTH);
                        if (value < minValue) {
                            minValue = value;
                        }
                        if (value > maxValue) {
                            maxValue = value;
                        }
                        numOfValues++;
                    }
                    hashValues->addWeight(addrInTH, 1);
                }
            }

            // printf("Found %u words in leaves values and %u values\n", this->numOfWords, numOfValues);

            /******************************************************************************/
            // SECOND STEP - Calculate entropy (vocabulary and values)
            /******************************************************************************/

            // Calculate entropy of values
            double entropyValues = 0;
            {
                uint totalNumOfValue = numOfSubmatrices * submatrixSize;
                double p;
                for (uint i = minValue; i <= maxValue; i++) {
                    wordValue = (unsigned  char*)&i;
                    j = hashValues->searchElement(wordValue, sizeof(uint), &addrInTH);
                    if (j != numOfValues) {
                        p = (double)hashValues->getWeight(addrInTH) / totalNumOfValue;
                        entropyValues += (p * log2(p));
                    }

                }
                entropyValues = -entropyValues;
                delete(hashValues);
            } // END Calculate entropy of values

            // Calculate entropy of vocabulary
            double entropyVoc = 0;
            {
                double p;
                for (uint i = 0; i < this->numOfWords; ++i) {
                    p  = hashHEntropy->getWeight(positionInTH[i]/*.at(i)*/) / (double)numOfSubmatrices;
                    entropyVoc += (p * log2(p));
                }
                entropyVoc = -entropyVoc;
            } // END Calculate entropy of vocabulary


            /******************************************************************************/
            // THIRD STEP - Sort vocabulary by frequency
            /******************************************************************************/

//            std::sort (positionInTH.begin(), positionInTH.end(), compareFreqListDescHEntropy2);
            std::qsort(positionInTH, this->numOfWords, sizeof(ulong), compareFreqListDescHEntropy);
//
            /******************************************************************************/
            // FOURTH STEP - Assign codewords to vocabulary
            /******************************************************************************/
            int codeWord = 0;
            uint maxCompactedValues=0, maxNoCompactedValues=0;
            {
                double pPlain, pVoc;
                uint p;
                for (uint i = 0; i < this->numOfWords; i++) {
                    p = hashHEntropy->getWeight(positionInTH[i]/*.at(i)*/);
                    pPlain = p * (submatrixSize * entropyValues);
                    pVoc = submatrixSize * sizeof(uint) * 8 + (p * entropyVoc);
                    if (pVoc <= pPlain) {
                        hashHEntropy->setCodeword(positionInTH[i]/*.at(i)*/, codeWord);
                        codeWord++;
                        maxCompactedValues += p;
                    } else {
                        hashHEntropy->setCodeword(positionInTH[i]/*.at(i)*/, -1);
                        maxNoCompactedValues += p;
                    }
                }

//                printf("Created vocabulary with %u words (of %u words). Value entropy -> %.03f and Voc entropy -> %.03f\n",
//                       codeWord, this->numOfWords, entropyValues, entropyVoc);
                this->numOfWords = codeWord;
            }

            /******************************************************************************/
            // FIFTH STEP - Compact
            /******************************************************************************/
            uint *compactedValues = (uint*)malloc(maxCompactedValues * sizeof(uint));
            uint *noCompactedValues = (uint *)malloc(maxNoCompactedValues * submatrixSize * sizeof(uint));
            uint *tmpBitmapVoc = (uint*)malloc(cds_static::uint_len(numOfSubmatrices, 1) * sizeof(uint));
            uint numOfCompacted = 0, numOfNoCompacted = 0;

            for (uint i = 0; i < numOfSubmatrices; i++) {
                word = &plainWords[i * this->lenWord];
                j = hashHEntropy->searchElement(word, this->lenWord, &addrInTH);
                codeWord = hashHEntropy->getCodeword(addrInTH);
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
//            printf("Reserved %u cells, used %u\n", numOfSubmatrices * submatrixSize, numOfNoCompacted);
            free(this->plainValues);
            this->plainValues = nullptr;

            /******************************************************************************/
            // SIXTH STEP - Copy words
            /******************************************************************************/

            //printf("Coping %u words ..........\n", this->numOfWords);
//            this->words = (unsigned char *)malloc(this->numOfWords * this->lenWord * sizeof(unsigned char));
            this->words = new uchar [this->numOfWords * this->lenWord];
            this->totalLengthOfWords = 0;
            for (uint w = 0; w < this->numOfWords; w++) {
                for (uint l = 0; l < this->lenWord; l++) {
                    this->words[this->totalLengthOfWords++] = hashHEntropy->getWord(positionInTH[w]/*.at(w)*/)[l];
                }
            }
            free(positionInTH);
//            positionInTH.clear();
//            positionInTH.shrink_to_fit();
            delete(hashHEntropy);

            /******************************************************************************/
            // SEVENTH STEP - Encode wit DAC and bitsequence
            /******************************************************************************/

            //printf("Creating bitmapVocLeaves with %u submatrices\n", numOfSubmatrices);
            cds_static::BitSequenceBuilder *builder = new cds_static::BitSequenceBuilderRG(20);
            this->bitmapVocLeaves = builder->build(tmpBitmapVoc, numOfSubmatrices);
            free(tmpBitmapVoc);
            delete(builder);
//            printf("Creating DAC for %u compacted values\n", numOfCompacted);
            this->vocLeaves = new cds_static::DAC(compactedValues, numOfCompacted, 3);
            free(compactedValues);
//            printf("Creating DAC for %u no compacted values\n", numOfNoCompacted);
            this->valuesLeaves = new cds_static::DAC(noCompactedValues, numOfNoCompacted, 3);
//            printf("Created DAC for %u no compacted values\n", numOfNoCompacted);
//            printf("free %u no compacted values\n", numOfNoCompacted);
            free(noCompactedValues);
        } else {
            // No plain levels
            uint * emptyArray = new uint[1];
            emptyArray[0] = 0;

            cds_static::BitSequenceBuilder *builder = new cds_static::BitSequenceBuilderRG(20);
            this->bitmapVocLeaves = builder->build(emptyArray, 1);
            delete(builder);
            this->vocLeaves = new cds_static::DAC(emptyArray, 1, 3);
            this->valuesLeaves = new cds_static::DAC(emptyArray, 1, 3);
            delete [] emptyArray;
        }// ENF IF
    }

    K2RasterEntropy::~K2RasterEntropy() {
        if (this->words != nullptr) {
            delete[](this->words);
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

    int K2RasterEntropy::getPlainCell(uint father, uint row, uint col, uint divK, uint lastValueMax) const {
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

    ulong K2RasterEntropy::getPlainWindow(uint xini, uint yini, uint xend, uint yend, uint valini, uint valend,
                                          uint offsetValues, uint divK, uint maxValue, uint *positions,
                                          ulong totalCells, bool allCells) const {
        // Add all valid positions to final result
        uint value, arrayPosition;

        // Check in bitmap if submatrix is in vocabulary or not
        if (this->bitmapVocLeaves->access(offsetValues)) {
            uint pos = this->bitmapVocLeaves->rank1(offsetValues) - 1;
            uint compressedValue = this->vocLeaves->access(pos);
            unsigned char * word = &this->words[compressedValue * this->lenWord];
            uint *values = (uint*)word;
            for (uint x = xini; x <= xend; x++){
                arrayPosition = (x % divK) * divK + (yini % divK);
                for (uint y= yini; y <= yend; y++){
                    value = maxValue - values[arrayPosition];
                    if (value >= valini && value <= valend) {
                        positions[totalCells++] = x * this->realSizeY + y;
                    } else {
                        if (allCells) {
                            return 0;
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
                if (value >= valini && value <= valend) {
                    positions[totalCells++] = x * this->realSizeY + yini;
                } else {
                    if (allCells) {
                        return 0;
                    }
                }
                for (uint y= yini+1; y <= yend; y++){
                    value = maxValue - this->valuesLeaves->next();
                    if (value >= valini && value <= valend) {
                        positions[totalCells++] = x * this->realSizeY + y;
                    } else {
                        if (allCells) {
                            return 0;
                        }
                    }
                }
            }
        }
        return totalCells;
    }

    bool K2RasterEntropy::checkPlainCells(uint xini, uint yini, uint xend, uint yend, uint offsetValues, uint divK,
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

    ulong K2RasterEntropy::getPlainCellValues(uint xini, uint yini, uint xend, uint yend,
                                     uint offsetValues, uint divK, uint maxValue,
                                     int *values, ulong totalCells,
                                     uint initXini, uint initYini, uint initYend) const {
        uint value, arrayPosition;
        uint cellPos;

        if (this->bitmapVocLeaves->access(offsetValues)) {
            uint pos = this->bitmapVocLeaves->rank1(offsetValues) - 1;
            uint compressedValue = this->vocLeaves->access(pos);
            unsigned char *word = &this->words[compressedValue * this->lenWord];
            uint *wordValues = (uint *) word;
            for (uint x = xini; x <= xend; x++) {
                arrayPosition = (x % divK) * divK + (yini % divK);
                cellPos = ((x - initXini) * ((initYend - initYini) + 1));
                for (uint y = yini; y <= yend; y++) {
                    value = maxValue - wordValues[arrayPosition];
                    // if (value >= query.valini && value <= query.valend) {
                    values[cellPos + (y - initYini)] = value;
                    totalCells++;
                    //}
                    arrayPosition++;
                }
            }
        } else {
            uint pos = (this->bitmapVocLeaves->rank0(offsetValues) - 1) * divK * divK;
            for (uint x = xini; x <= xend; x++){
                arrayPosition = (x % divK) * divK + (yini % divK);
                cellPos = ((x - initXini) * ((initYend - initYini) + 1));
                value = maxValue - this->valuesLeaves->access(pos + arrayPosition);
                values[cellPos + (yini - initYini)] = value;
                totalCells++;
                for (uint y= yini+1; y <= yend; y++){
                    value = maxValue - this->valuesLeaves->next();
                    values[cellPos + (y - initYini)] = value;
                    totalCells++;
                }
            }
        }
        return totalCells;
    }


    //**********************************************************************//
    //********************* FILE FUNCTIONS *********************************//
    //**********************************************************************//
    void K2RasterEntropy::save(std::ofstream &of) const {
        K2Raster::save(of, K2RASTER_VOC_ENTROPY);

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
        this->saveSplit("proba");
    }

    void K2RasterEntropy::saveSplit(char *filename) const {
//        K2Raster::saveSplit(filename, K2RASTER_VOC_ENTROPY);

        // Args *************************
        {
            std::string outFileName = filename;
            outFileName += ".argsEntropy.k2raster";
            ofstream outFile(outFileName.data());
            if (!outFile.good()) {
                printf("File %s unable to open\n", outFileName.data());
                return;
            }

            saveValue(outFile, this->levelK2);

            outFile.close();
        } // END BLOCK args

        // Vocabulary *************************
        {
            std::string outFileName = filename;
            outFileName += ".vocabulary.k2raster";
            ofstream outFile(outFileName.data());
            if (!outFile.good()) {
                printf("File %s unable to open\n", outFileName.data());
                return;
            }

            saveValue(outFile, this->numOfWords);
            saveValue(outFile, this->totalLengthOfWords);
            saveValue(outFile, this->words, this->totalLengthOfWords);

            outFile.close();
        } // END BLOCK Vocabulary

        // DAC encodedValues *************************
        {
            std::string outFileName = filename;
            outFileName += ".encodedValues.k2raster";
            ofstream outFile(outFileName.data());
            if (!outFile.good()) {
                printf("File %s unable to open\n", outFileName.data());
                return;
            }

            this->vocLeaves->save(outFile);

            outFile.close();
        } // END BLOCK DAC values of leaves

        // DAC plainValue *************************
        {
            std::string outFileName = filename;
            outFileName += ".plainValues.k2raster";
            ofstream outFile(outFileName.data());
            if (!outFile.good()) {
                printf("File %s unable to open\n", outFileName.data());
                return;
            }
            this->valuesLeaves->save(outFile);

            outFile.close();
        } // END BLOCK DAC values of leaves

        // bitmap vocabulary of leaves *************************
        {
            std::string outFileName = filename;
            outFileName += ".bitmapIsInVoc.k2raster";
            ofstream outFile(outFileName.data());
            if (!outFile.good()) {
                printf("File %s unable to open\n", outFileName.data());
                return;
            }

            this->bitmapVocLeaves->save(outFile);

            outFile.close();
        } // END BLOCK DAC values of leaves
    }

    /*
     * Load Ktreap from file
     */
    K2RasterEntropy *K2RasterEntropy::load(std::ifstream &in) {
        K2RasterEntropy *ktreap = nullptr;

        try {
            ktreap = new K2RasterEntropy();
            K2Raster::loadBase(in, ktreap);

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
    size_t K2RasterEntropy::getTotalSize() const {
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

    size_t K2RasterEntropy::printSize() const {
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
    int compareFreqListDescHEntropy(const void *a, const void *b) {
        unsigned int *left,*right;
        left =  (uint*) a;
        right = (uint*) b;
        if (hashHEntropy->getWeight(*left) < hashHEntropy->getWeight(*right))
            return 1;
        if (hashHEntropy->getWeight(*left) > hashHEntropy->getWeight(*right))
            return -1;

        return 0;
    }

    bool compareFreqListDescHEntropy2(ulong left, ulong right) {
        return hashHEntropy->getWeight(left) < hashHEntropy->getWeight(right);
    }
}

