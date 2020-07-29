k<sup>2</sup>-raster -  A compacta data strucutre for raster data
=========================

We present the k2-raster,  a new technique for representing raster datasets that uses compressed space and offers indexing capabilities, thus, improving query times over the raster data.

Let M be a raster matrix of size n×n, being n a power of k, containing values v ≥ 0 for each cell M<sub>i</sub><sub>j</sub>. 
The k<sup>2</sup>-raster recursively partitions the matrix M into k<sup>2</sup> submatrices, analogous to the original k<sup>2</sup>-tree, 
and builds a tree representing this recursive subdivision. 
In addition, the representation is coupled with an efficient representation of the maximum and minimum values of each submatrix, 
which are needed for the representation of the raster data, but also provide the indexing functionality.
 
 
## Compilation ##
 
 To compile the code you need to execute the following commands:
 ```bash
 cd k2-raster
 mkdir build
 cd build
 cmake ..
 make DEBUG=FALSE
 ```
 
 Executables files are stored in
 * ./build/bin/
 
 Library files are stored in
 * ./build/lib/

 
## Encode data ##
 
To compress a raster, you have to run the command:

 ```bash
 ./Encodek2Raster <dataFile> <rows> <cols> <k2rasterFile> <check> [type] [k1 k2 levelK1 [plainLevels]] [minFreqVoc]  
 ```
 
 Where:
 * **\<dataFile>** input file (list of intergers of 32 bits)
 * **\<rows>** number of rows
 * **\<cols>** number of columns
 * **\<treapFile>** path of the file where store the k<sup>2</sup>-raster
 * **\<check>** A 1 if after constructing the structure performs a check, cell by cell, of the values. Compare the original value with the value returned by the structure. Or a 0 otherwise
 * **\<type>** type of  k<sup>2</sup>-raster which want to use to encode 
 * **\<k1>** k for the first levels. It partitions a matrix into k1 equal-size submatrices.
 * **\<k2>** k for the rest of the levels. It partitions a matrix into k1 equal-size submatrices
 * **\<levelk1>** number of levels we use the partition k1
 * **\<plainLevels>** number of levels represented as plain values (no used with type 2 and 8)
 * **\<minFreqVoc>** (only for type 7) minimum frequency of a word in the vocabulary.
 

#### k<sup>2</sup>-raster types ####
 
 * **Type 1**: -- NOT USED --
 * [**Type 2**](src/k2raster/k2-raster.cpp): Standard k2-raster. Version used on [__*SSDBM16*__](http://doi.org/10.1145/2949689.2949710)
 * [**Type 3**](src/k2raster/plain/k2-raster-plain.cpp): Stored the lastest \<plainLevels> as plain values (array of integers of 32 bits)
 * [**Type 4**](src/k2raster/plain/k2-raster-plain-DAC.cpp): Stored the lastest \<plainLevels> using a DAC encode. 
 * [**Type 5**](src/k2raster/plain/k2-raster-plain-VByte.cpp): Stored the lastest \<plainLevels> using VByte encode.
 * [**Type 6**](src/k2raster/compresessLeaves/k2-raster-CompressLeaves.cpp): Creates a vocabulary to represent all submatrices defined by \<plainLevels>
 * [**Type 7**](src/k2raster/compresessLeaves/k2-raster-CompressLeavesH.cpp): Creates a vocabulary to represent all sumbatrices defined by \<plainLevels> than have a frequency greater than \<minFreVoc>. The rest of submatrices are stored as plain values (array of interges of 32 bits) 
 * [**Type 8**](src/k2raster/k2-raster-opt.cpp): Optimization of **Type 2**. It executes the function *next()* instead of *access()* when we want to obtain a value from a DAC.
 * [**Type 9**](src/k2raster/k2-raster-entropy.cpp): It uses a heuristic function to decide if a submatrix is stored in the vocabulary or in plain. The heuristic function determines which of the two options uses less space.
 
 

#### Example ####
 
For example, a 5x5 raster data (5 rows and 5 columns) is stored in _./inputData.bin_ and each cell has an integer value of 32 bits.

\[ 4 5 2 6 8 ... 33 45 23 25 ] (Total 25 numbers, the first 5 numbers correspond to the first row, the next five numbers correspond to the second row and so on) 
 
 To use **Type 2** and *k1* equal to 4 during the 4 first levels and *k2* equal to 2 for the rest of the levels.
 
```bash
 ./Encodek2Raster ./inputData.bin 5 5 outputFile.k2r 0 2 4 2 4
```

## Queries ##

 We have implemented 4 useful queries to retrieve information using indexing techniques to improve the time consumption:

#### Obtaining a cell value ####

 Given a position in the raster matrix, this query obtains its cell value.

 ```bash
 ./GetCellValuek2Raster <k2rasterFile> <queriesFile>
 ```
 Where:
 
 * **\<k2rasterFile>** Path of the file where the k<sup>2</sup>-raster is stored 
 * **\<queriesFile>**  File with a set of queries. Each line is a query with format "posX posY"

#### Obtaining all the values of a region ####

 Given a range of values, this query retrieves all raster positions whose value lies within that range.

 ```bash
 ./GetValuesWindowk2Raster <k2rasterFile> <queriesFile> <check>
 ```
 Where:
 
 * **\<k2rasterFile>** Path of the file where the k<sup>2</sup>-raster is stored 
 * **\<queriesFile>**  File with a set of queries. Each line is a query with format "posX1 posX2 posY1 posY2 val1 val2" (val1 and val2 are not used in this version)
 * **\<check>** 1 to check if the final result is correct. 0 in another case.
 

#### Retrieving cells with a given value or range of values ####

Given a range of values, this query retrieves all raster positions whose value lies within that range.

 ```bash
 ./GetCellsRangeValuek2Raster <k2rasterFile> <queriesFile> <check>
 ```
 Where:
 
 * **\<k2rasterFile>** Path of the file where the k<sup>2</sup>-raster is stored 
 * **\<queriesFile>**  File with a set of queries. Each line is a query with format "posX1 posX2 posY1 posY2 val1 val2"
 * **\<check>** 1 to check if the final result is correct. 0 in another case.

#### Checking the existence of a given value or range of values ####

Given a region and a range, this query checks if all cell values of the region are within the range of values or if there exists at least one cell value in the region whose value lies within the range of values.

 ```bash
 ./CheckValuesk2Raster <k2rasterFile> <queriesFile> <check> [allCells]
 ```
 Where:
 
 * **\<k2rasterFile>** Path of the file where the k<sup>2</sup>-raster is stored 
 * **\<queriesFile>**  File with a set of queries. Each line is a query with format "posX1 posX2 posY1 posY2 val1 val2"
 * **\<check>** 1 to check if the final result is correct. 0 in another case.
 * **\<allCells>** sets whether all cells or at least one must meet query conditions

#### Spatial join - k<sup>2</sup>-raster  and R-tree ####

```bash
 ./joinRTreek2Raster <k2rasterFile> <rtreeFile> <queryFile> <allCells> [strategyType]
 ```
 
 Where:
  
 * **\<k2rasterFile>** Path of the file where the k<sup>2</sup>-raster is stored 
 * **\<rtreeFile>** Path of the file where the R-tree is stored 
 * **\<queriesFile>**  File with a set of queries. Each line is a query with format "posX1 posX2 posY1 posY2 val1 val2"
 * **\<allCells>** sets whether all cells or at least one must meet query conditions
 * **\<strategyType>** Select the strategy used during the spatial join:
    * [**Type 0**](src/rtree/SpatialJoin/QueryStrategy/GetObjectsQueryStrategyQueue.cpp). Using a queue: Submatrices to be processed are added to a queue.
    * [**Type 1**](src/rtree/SpatialJoin/QueryStrategy/GetObjectsQueryStrategyTree.cpp). Using a tree: A tree is created where it stores by levels the submatrices to be processed.
    * [**Type 2**](src/rtree/SpatialJoin/QueryStrategy/GetObjectsQueryStrategyTree2.cpp). Using a tree V2: A tree is created where it stores by levels the submatrices to be processed (version 2).
    * [**Type 3**](src/rtree/SpatialJoin/QueryStrategy/GetObjectsQueryStrategyLeaves.cpp). Using leaves: First, the leaves of the R-tree are obtained and afterward their cells are searched in the k<sup>2</sup>-raster.
    * [**Type 4**](src/joinRTree.cpp). Using cells: First, the cells of the k2-raster are obtained and then they are searched in the R-tree if they match in a leaf.

## Publications ##

Ladra, S., Paramá, J. R., & Silva-Coira, F. (2016). Compact and queryable representation of raster datasets. In Proceedings of the 28th International Conference on Scientific and Statistical Database Management - SSDBM ’16 (pp. 1–12). New York, New York, USA: ACM Press. http://doi.org/10.1145/2949689.2949710

Ladra, S., Paramá, J. R., & Silva-Coira, F. (2017). Scalable and queryable compressed storage structure for raster data, Information System, 72 (2017). doi:10.1016/j.is.2017.10.007.

## Authors ##

Susana Ladra - <susana.ladra@udc.es>     
José R. Paramá - <jose.parama@udc.es>    
Fernando Silva-Coira - <fernando.silva@udc.es>

#### Implemented By ####
Fernando Silva-Coira - <fernando.silva@udc.es>


