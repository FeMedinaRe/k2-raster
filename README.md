k<sup>2</sup>-raster -  A compacta data strucutre for raster data
=========================

The k2-raster is structure for representing raster datasets that uses compressed space and offers indexing capabilities, thus, improving query times over the raster data.

Let M be a raster matrix of size n×n, being n a power of k, containing values v ≥ 0 for each cell M<sub>i</sub><sub>j</sub>. 
The k<sup>2</sup>-raster recursively partitions the matrix M into k<sup>2</sup> submatrices, analogous to the original k<sup>2</sup>-tree, 
and builds a tree representing this recursive subdivision. 
In addition, the representation is coupled with an efficient representation of the maximum and minimum values of each submatrix, 
which are needed for the representation of the raster data, but also provide the indexing functionality.

## Compilation ##
To compile the code you need to execute the following commands:
```bash
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

or 

```bash
sh compile.sh
```

Executables are stored in
 * ./build/bin/
 
 Libraries are stored in
 * ./build/lib/
 
 ## Encode data ##
  
 To compress a raster, you run the  next command:
 
  ```bash
  ./k2_raster_encode <input_data> <rows> <cols> <output_data> <set_check> <type> <k1> [k2 level_k1 [plain_levels]]  
  ```
 Where:
 * **\<input_data>** input file (list of values)
 * **\<rows>** number of rows
 * **\<cols>** number of columns
 * **\<output_data>** path of the file where store the k<sup>2</sup>-raster
 * **\<set_check>** Testing the output k<sup>2</sup>-raster. A 1 if after building the structure performs a check, cell by cell, of the values. Compare the original value with the value returned by the structure. Or a 0 otherwise
 * **\<type>** type code of the k<sup>2</sup>-raster.
 * **\<k1>** k for the first levels. It partitions a matrix into k1 equal-size submatrices.
 * **\<k2>** k for the rest of the levels. It partitions a matrix into k1 equal-size submatrices
 * **\<level_k1>** number of levels we use the partition k1
 * **\<plain_levels>** number of levels represented as plain values
 
 #### k<sup>2</sup>-raster types ####
 * [**k<sup>2</sup>-raster**](include/k2_raster.hpp): (Type code 1). Standard k2-raster. Version used on [__*SSDBM16*__](http://doi.org/10.1145/2949689.2949710)
 * [**k<sup>2</sup>-raster_plain**](include/k2_raster_plain.hpp): (Type code 10). Stored the lastest \<plain_levels> as plain values
 * [**k<sup>2</sup><sub>H</sub>-raster (*heuristic k<sup>2</sup>-raster*)**](include/k2_raster_heuristic.hpp): (Type code 11). It uses a heuristic function to decide if a submatrix is stored in the vocabulary or in plain. The heuristic function determines which of the two options uses less space. Version used on [__*Information System (2017)*__](http://doi.org/10.1016/j.is.2017.10.007)


 #### Example ####
  
 For example, a 5x5 raster data (5 rows and 5 columns) is stored in _./input_data.bin_ and each cell is represented with an integer value of 32 bits.
 
 \[ 4 5 2 6 8 ... 33 45 23 25 ] (Total 25 numbers, the first 5 numbers correspond to the first row, the next five numbers correspond to the second row and so on) 
  
  To use the standard **k2-raster** *(Type code 1)* and *k1* equal to 4 for the 4 first levels and *k2* equal to 2 for the rest of the levels. Testing is enabled.
  
 ```bash
  ./k2_raster_encode ./input_ata.bin 5 5 output_file.k2r 1 1 4 2 4
 ```
 
 ## Queries (In progress.........) ##
 
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
 * **\<queriesFile>**  File with a set of queries. Each line is a query with format "posX1 posX2 posY1 posY2 val1 val2" (val1 and val2 are not used on this version)
 * **\<check>** 1 to check if the final result is correct. 0 in other case.
 

#### Retrieving cells with a given value or range of values ####

Given a range of values, this query retrieves all raster positions whose value lies within that range.

 ```bash
 ./GetCellsRangeValuek2Raster <k2rasterFile> <queriesFile> <check>
 ```
 Where:
 
 * **\<k2rasterFile>** Path of the file where the k<sup>2</sup>-raster is stored 
 * **\<queriesFile>**  File with a set of queries. Each line is a query with format "posX1 posX2 posY1 posY2 val1 val2"
 * **\<check>** 1 to check if the final result is correct. 0 in other case.

#### Checking the existence of a given value or range of values ####

Given a region and a range, this query checks if all cell values of the region are within the range of values or if there exists at least one cell value in the region whose value lies within the range of values.

 ```bash
 ./CheckValuesk2Raster <k2rasterFile> <queriesFile> <check> [allCells]
 ```
 Where:
 
 * **\<k2rasterFile>** Path of the file where the k<sup>2</sup>-raster is stored 
 * **\<queriesFile>**  File with a set of queries. Each line is a query with format "posX1 posX2 posY1 posY2 val1 val2"
 * **\<check>** 1 to check if the final result is correct. 0 in other case.
 * **\<allCells>** sets whether all cells or at least one must meet query conditions


## Publications ##

Ladra, S., Paramá, J. R., & Silva-Coira, F. (2016). Compact and queryable representation of raster datasets. In Proceedings of the 28th International Conference on Scientific and Statistical Database Management - SSDBM ’16 (pp. 1–12). New York, New York, USA: ACM Press. [doi:10.1145/2949689.2949710](http://doi.org/10.1145/2949689.2949710)

Ladra, S., Paramá, J. R., & Silva-Coira, F. (2017). Scalable and queryable compressed storage structure for raster data, Information System, 72 (2017). [doi:10.1016/j.is.2017.10.007](http://doi.org/10.1016/j.is.2017.10.007).

## Web Sites ##

* [Official Web Site](https://lbd.udc.es/research/k2-raster/)
* [Code](https://gitlab.lbd.org.es/fsilva/k2-raster)


## Authors ##
Susana Ladra - <susana.ladra@udc.es>     
José R. Paramá - <jose.parama@udc.es>    
Fernando Silva-Coira - <fernando.silva@udc.es>

### Implemented By ###
Fernando Silva-Coira - <fernando.silva@udc.es>

### Acknowledgements ###
This work was supported by Ministerio de Economía y Competitividad (PGE and FEDER) [grant numbers TIN2016-78011-C4-1-R ; TIN2016-77158-C4-3-R ; TIN2013-46238-C4-3-R ; TIN2013-46801-C4-3-R ], Centro para el desarrollo Tecnológico e Industrial [grant numbers IDI-20141259; ITC-20151247, ITC-20161074]; Xunta de Galicia (co-founded with FEDER) [grant numbers ED431C 2017/58; ED431G/01].


  
  