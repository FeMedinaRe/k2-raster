/*
 * Created by Fernando Silva on 14/03/16.
 *
 * Copyright (C) 2016-current-year, Fernando Silva, all rights reserved.
 *
 * Author's contact: Fernando Silva  <fernando.silva@udc.es>
 * Databases Lab, University of A Coruña. Campus de Elviña s/n. Spain
 *
 * Base class with all general params and functions to use the k2-raster
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

#ifndef K2RASTER_H_
#define K2RASTER_H_

// Type of k2-rasters
#define K2RASTER 2
#define K2RASTER_PLAIN 3
#define K2RASTER_PLAIN_DAC 4
#define K2RASTER_PLAIN_VBYTE 5
#define K2RASTER_VOC 6
#define K2RASTER_VOC_HYBRID 7
#define K2RASTER_OPT 8
#define K2RASTER_VOC_ENTROPY 9


// System libraries
#include <map>
#include <vector>

// Own headers
#include <k2raster/util/NodeMatrix.h>
#include <k2raster/util/QNode.h>
#include <k2raster/util/InfoMinSubMatrix.h>

// Libraries
#include <bitsequence/BitSequenceBuilder.h>
#include <direct_access/DirectAccess.h>
#include <queries/Query.h>


using namespace query_utils_static;

namespace k2raster_static {

/****** Operations ******/
enum OperationRaster {
    OPERATION_SUM,
    OPERATION_SUBT,
    OPERATION_MULT
};

class K2Raster {
 public:
    /****** Constructors ******/
    K2Raster(uint sizeX, uint sizeY, int **data, uint K1, uint k2, uint levelK1);
    K2Raster(uint K1, uint k2, uint levelK1,
             K2Raster *k2raster1, K2Raster *k2raster2, OperationRaster operation,
             bool deleteInputs=false);
    K2Raster(K2Raster *k2raster1, K2Raster *k2raster2, OperationRaster operation,
             bool deleteInputs=false); // equal K's


    virtual ~K2Raster();

    /****** Decompression ******/
    int *decompress();

    /****** Cell functions ******/
    // GetCell
    virtual int getCell(uint row, uint col) const;

    // GetCellsByValue
    virtual ulong getCellsByValue(Query query, uint *positions) const;
    virtual ulong getCellsByValue(uint xini, uint xend, uint yini, uint yend, uint valini, uint valend,
                                  BasicSubMatrix *minSubMatrix, uint *positions, bool allCells) const;

    // checkValues
    virtual bool checkValues(Query query, bool allCells) const;

    // getValuesWindow
    virtual ulong getValuesWindow(Query query, int *values) const;
    virtual ulong getValuesWindow(uint xini, uint xend, uint yini, uint yend, int *values) const;

    // GetMaxValueWindow
    virtual int getMaxValueWindow(Query query);
    virtual int getMaxValueWindow(uint xini, uint xend, uint yini, uint yend);
    virtual int getMaxValueWindow(uint xini, uint xend, uint yini, uint yend,
                                  BasicSubMatrix *basicSubMatrix);
    virtual int getMaxValueWindow(uint xini, uint xend, uint yini, uint yend,
                                  uint baseRow, uint baseCol,
                                  BasicSubMatrix *basicSubMatrix); // TODO change to protected

    // GetMinValueWindow
    virtual int getMinValueWindow(Query query);
    virtual int getMinValueWindow(uint xini, uint xend, uint yini, uint yend);
    virtual int getMinValueWindow(uint xini, uint xend, uint yini, uint yend,
                                  BasicSubMatrix *basicSubMatrix);
    virtual int getMinValueWindow(uint xini, uint xend, uint yini, uint yend,
                                  uint baseRow, uint baseCol,
                                  BasicSubMatrix *basicSubMatrix); // TODO change to protected


    /****** Node functions ******/

    // getMinMatrix
    virtual BasicSubMatrixWithUsed *getMinMatrix(uint xini, uint xend, uint yini, uint yend,
                                                 BasicSubMatrixWithUsed *minSubMatrix) const;

    // With color
    // V1 (using a queue)
    virtual MinSubMatrix *getMinMatrix(uint xini, uint xend, uint yini, uint yend,
                                       uint valini, uint valend,
                                       MinSubMatrix *minSubMatrix) const;

    // Using a tree of precalculated submatrices
    virtual MinNodeTreeNode *getMinMatrix(uint xini, uint xend, uint yini, uint yend,
                                            uint valini, uint valend,
                                            MinNodeTreeNode *initialNode,
                                            std::map<ulong, MinNodeTreeNode *> &nextMBRs) const;

    // V2 Using a tree of precalculated submatrices
    virtual MinNodeTreeNode *getMinMatrixV2(uint xini, uint xend, uint yini, uint yend,
                                              uint valini, uint valend,
                                                MinNodeTreeNode *minSubMatrix) const;

    BasicSubMatrix *createRootBasicSubMatrix() const;
    BasicSubMatrixWithUsed *createRootBasicSubMatrixWithUsed() const;
    MinSubMatrix *createRootMinSubMatrix(uint valini, uint valend) const;
    MinNodeTreeNode *createRootNodeMinSubMatrix(uint valini, uint valend, bool addChildren = false) const;


    /****** INFO functions ******/
    virtual uint getMinValue() const;
    virtual uint getMaxValue() const;
    virtual uint getMinValueNoZero() const;
    virtual uint getNumOfPositions () const;
    virtual uint getDimSize(uint dim) const;

    virtual std::string getTypeName() const;

    static std::string getTypeName(uint type);

    /****** File functions ******/
    virtual void save(std::ofstream &of) const;
    virtual void save(std::ofstream &of, uint type) const;
    virtual void saveSplit(char *filename) const;
    virtual void saveSplit(char *filename, uint type) const;
    static K2Raster *load(std::ifstream &in);

    static K2Raster *createK2Raster(uint sizeX, uint sizeY, int **data, uint k1, uint k2, uint levelK1,
                                    uint plainLevels = 2, float minFreqVoc = 0.01, uint type = 1);

    /****** Size functions ******/
    virtual size_t getTotalSize() const;
    virtual size_t printSize() const;

    /****** Check functions ******/
    bool checkEncode(int **data, uint nrows, uint ncols) const;

protected:
    /****** Basic params ******/
    uint realSizeX;                             // Real row size of matrix
    uint realSizeY;                             // Real column size of matrix
    uint size;                                  // Virtual size of the matrix (power of K)
    uint K1;                                    // Division for first levels
    uint K2;                                    // Division for last levels
    uint levelK1;                               // Number of level that apply K1 subdivision
    uint *divKLevels;                           // Partition point at each level
    uint nlevels;                               // Height of the tree (+1)
    uint countOnesk1;                           // Number of "ones" that use 'k1'
    uint displacement;                          // Displacement to avoid negative values

    /****** Structures ******/
    cds_static::BitSequence *bitmapK2Tree;      // k2-tree representation of the data: single bitmap for T:L

    /****** Values ******/
    uint numOfMaxDACs;
    uint numOfMinDACs;
    cds_static::DirectAccess **listMaxValues;   // List of max values by level
    cds_static::DirectAccess **listMinValues;   // List of min values by level
    uint initialMaxValue;                        // Max value of whole matrix
    uint initialMinValue;                        // Min value of whole matrix
    uint minValueNoZero;                        // Min value of whole matrix (and it is not a zero)

    ulong *preMinValues;                        // List with the total number of previous min values at level "l"
    ulong *preMaxValues;                        // List with the total number of previous max values at level "l"

    /****** Auxiliary params (only construction) ******/
    std::vector<uint> *tmpPlainValues;

    /****** Constructor ******/
    K2Raster();

    /****** Decompress ******/
    void decompress(uint currentLevel, uint childrenPosition, uint maxValueParent, uint baseRow, uint baseCol,
                    uint *countOnes, int *data);

    /****** Aux functions ******/
    uint build(uint sizeX, uint sizeY, int **data);
    uint buildOperation(K2Raster *k2Raster1, K2Raster *k2Raster2, OperationRaster operation, bool deleteInputs);
    uint buildOperationEqualKs(K2Raster *k2Raster1, K2Raster *k2Raster2, OperationRaster operation, bool deleteInputs); // Equal K's

    NodeMatrix* checkSubmatrix(uint baseRow, uint baseCol, uint subMatrixSize,
                               uint level, int **data, uint &totalNodes,
                               std::vector<uint> **maxValues, std::vector<uint> **minValues,
                               std::vector<uint> *plainValues, uint **tmpBitmap, uint *levelPosition);
    NodeMatrix* checkSubmatrixOperation(uint baseRow, uint baseCol, uint currentLevel, uint &totalNodes,
                                        std::vector<uint> **maxValues, std::vector<uint> **minValues,
                                        std::vector<uint> *plainValues, uint **tmpBitmap, uint *levelPosition,
                                        BasicSubMatrixWithChildren **ptrK2raster1, uint numOfNodes1,
                                        BasicSubMatrixWithChildren **ptrK2raster2, uint numOfNodes2,
                                        K2Raster *k2Raster1, K2Raster *k2Raster2, OperationRaster operation);
    NodeMatrix* checkSubmatrixOperation(uint baseRow, uint baseCol, uint currentLevel, uint &totalNodes,
                                        std::vector<uint> **maxValues, std::vector<uint> **minValues,
                                        std::vector<uint> *plainValues, uint **tmpBitmap, uint *levelPosition,
                                        BasicSubMatrix *ptrK2raster1, BasicSubMatrix *ptrK2raster2,
                                        K2Raster *k2Raster1, K2Raster *k2Raster2, uint *countOnes1, uint *countOnes2,
                                        OperationRaster operation); // Equal K's

    // Map algebra - Nodes functions
    inline BasicSubMatrix *getFirstChild(K2Raster *k2Raster, BasicSubMatrix * ptrK2raster, bool uniformNode, uint *countOnes);
    inline void getNextChild(K2Raster *k2Raster, BasicSubMatrix * ptrK2raster, BasicSubMatrix * ptrChild, uint *countOnes);
    BasicSubMatrixWithChildren **getNodesContained(uint i, uint j, uint k, uint size,
                                                                  BasicSubMatrixWithChildren **ptrK2raster, uint numOfNodes,
                                                                  K2Raster *raster, uint &numOfNodesContained);
    BasicSubMatrixWithChildren ** getChildrenNode(BasicSubMatrixWithChildren *node, K2Raster *raster, uint &numOfChildren);

    uint conceptualToCompact(NodeMatrix *root, uint totalNodes, std::vector<uint> **maxValues, std::vector<uint> **minValues,
                             uint **tmpBitmap, uint *levelPositions);
    void initLevels();
    uint getK(uint level) const;
    virtual bool isLevelInPlain(uint level) const;

    /****** Values ******/
    void checkData(uint sizeX, uint sizeY, int **data);

    /****** Cell aux functions ******/
    virtual ulong getCellsByValue(uint xini, uint xend, uint yini, uint yend, uint fatherMinValue, uint fatherMaxValue,
                                  uint l, uint bitmapChildren, uint valini, uint valend,
                                  uint baseRow, uint baseCol,
                                  ulong totalCells, uint *positions, bool allCells) const;

    /****** Files ******/
    static K2Raster *loadBase(std::ifstream &in, K2Raster *ktreap);


    /****** DACs functions ******/
    size_t getSizeDACs(cds_static::DirectAccess **dacs, uint levels) const;
    size_t printSizeDACs(cds_static::DirectAccess **dacs, uint levels, bool printEach) const;


    /****** Nodes functions ******/
    QNodeCounters * createRoot() const;
    QNodeCounters * createChild(QNodeCounters *parentNode, uint offsetValues,
                        uint thisbaserow, uint thisbasecol, uint minValue, uint maxValue) const;

    int nextChild(QNodeCounters *parent, int &posI, int &posJ,
                  int &thisbaserow, int &thismaxrow, int &thisbasecol, int &thismaxcol,
                  Query query, uint newK, uint kLevel) const;
    int nextChild(QNodeCounters *parent, int &posI, int &posJ,
                  int &thisbaserow, int &thismaxrow, int &thisbasecol, int &thismaxcol,
                  uint xini, uint xend, uint yini, uint yend,
                  uint newK, uint kLevel) const;


    /****** Submatrix functions ******/
    void initSubMatrix(BasicSubMatrix *subMatrix, uint valini, uint valend) const;
    void initMinSubMatrix(MinBasicSubMatrix *subMatrix, uint valini, uint valend) const;
    MinNodeTreeNode *createNexMinNodeTreeNode(uint bitmapPosition, uint level, MinNodeTreeNode *node,
                                                      uint valini, uint valend, bool addChildren = false) const;

    void createChildrenNodeMinSubMatrix(MinNodeTreeNode *node, uint k, uint valini, uint valend,
                                            std::map<ulong, MinNodeTreeNode *> *nextMBRs,
                                            bool addChildren = false) const;
};
}
#endif // K2RASTER_H_
