k<sup>2</sup>-raster -  A compact data strucutre for raster data
=========================

The k<sup>2</sup>-raster is structure for representing raster datasets that uses compressed space and offers indexing capabilities, thus, improving query times over the raster data.

Let M be a raster matrix of size n×n, being n a power of k, containing values v ≥ 0 for each cell M<sub>i</sub><sub>j</sub>. 
The k<sup>2</sup>-raster recursively partitions the matrix M into k<sup>2</sup> submatrices, analogous to the original k<sup>2</sup>-tree, 
and builds a tree representing this recursive subdivision. 
In addition, the representation is coupled with an efficient representation of the maximum and minimum values of each submatrix, 
which are needed for the representation of the raster data, but also provide the indexing functionality.

## Table of Contents

[[_TOC_]]

## 1. Compilation

### 1.1. Dependencies
It is necessary **cmake** (minimum version 2.8.7), a c++ compiler (**g++**), and **git**.

* cmake (minimum version 2.8.7)
* A c++ compiler (g++)
* git


 ```bash
sudo apt install cmake g++ git
 ```

### 1.2. Compilation

1. Download the project from the repository.
    ```bash
    git clone https://gitlab.lbd.org.es/fsilva/k2-raster
    ```
2. Enter the project folder.
    ```bash
    cd k2-lidar
    ```
3. Run *compile.sh* script:
    ```bash
    sh compile.sh
    ```
3. (alternative) Or the following lines:
    ```bash
    sh ./clean.sh
    git submodule init
    git submodule update --init --recursive
    mkdir -p build
    cd build
    cmake ..
    make
    ```

These files can be easily removed with:
 ```bash
 sh clean.sh
 ```

### 1.3. Files

Executables files are stored in
* ./build/bin/

Library files are stored in
* ./build/lib/


# Usage
## 2.1 Raster data

The input is a raster file represented as a file of integers of 32-bits in row-major order. 

### 2.1. Encode
 
 The next command compress a raster file:

Usage:
 ```bash
./encode_k2r <input_data> <rows> <cols> <output_data> [-c] [-t <type>] [-k <k1>] [-q <k2> -l <level_k1> [-p <plain_levels>]] 
 ```

Where:

* **\<input_data>**: [string] - path to the input file (list of binary values).
* **\<rows>**: [number] - number of rows.
* **\<cols>** [number] - number of columns.
* **\<output>**: [string] - path to the output file where store the raster result.
* **-k \<k1>**: [number] - *k* value for the first levels (def. 2).
* **-q \<k2>**: [number] - *k* value for the last levels (def. 0, same value as k1).
* **-l \<level_k1>**: [number] - Number of level for k1. (def. 0).
* **-p \<plain_levels>**: [number] - Number of levels encoded with heuristic (11) or as plain (10) type.
* **-c**: [boolean] - Enable checking process. (def. false).
* **-t \<type>**: [number] - raster type (def. 10):\n"
    * 1 ->  hybrid k2-raster.
    * 10 ->  plain k2-raster.
    * 11 ->  heuristic k2-raster.
    
 More info: [k2_raster_encode](src/k2_raster_encode.cpp).

#### Example
We have a raster file of size 5x5 called "input_data.bin". Each cell is represented with an integer values of 32 bits, and they are stored following the row-major order.  For example, [4 5 2 6 8 ... 33 45 23 25 ], a total of 25 values, where the first five numbers correspond to the first row, the next five numbers correspond to the second row and so on.

We want to encode the raster using the standard **k<sup>2</sup>-raster** *(Type code 1)* with *k<sub>1</sub>*=4 for the 4 first levels and *k<sub>2</sub>*=2 for the rest of the levels. The result will be stored in a file called "output_file.k2r". Testing is enabled.

The corresponding command is:
 ```bash
  ./encode_k2r ./input_data.bin 5 5 output_file.k2r -c  -t 1 -k 4 -q 2 -l 4
 ```

 ### 2.2. k<sup>2</sup>-raster types
 * [**k<sup>2</sup>-raster**](include/k2_raster.hpp): (Type code 1). Standard k2-raster. Version used on [__*SSDBM16*__](http://doi.org/10.1145/2949689.2949710)
 * [**k<sup>2</sup>-raster_plain**](include/k2_raster_plain.hpp): (Type code 10). Stored the lasts level as plain values.
 * [**k<sup>2</sup><sub>H</sub>-raster (*heuristic k<sup>2</sup>-raster*)**](include/k2_raster_heuristic.hpp): (Type code 11). It uses a heuristic function to decide if a submatrix is stored in the vocabulary or in plain. The heuristic function determines which of the two options uses less space. Version used on [__*Information System (2017)*__](http://doi.org/10.1016/j.is.2017.10.007)



### 2.3. Queries
 
  We have implemented 4 useful queries to retrieve information using indexing techniques to improve the time consumption:
  
#### 2.3.1. Obtaining a cell value ([get_cell](src/k2_raster_get_cell.cpp)) ####

 Given a position in the raster matrix, this query obtains its cell value.
 
 ```bash
 ./get_cell_k2r <k2_raster_file> <queries_file> [-n <n_reps>]
 ```
 Where:
 
 * **\<k2_raster_file>**: [string] -  Path of the file where the k<sup>2</sup>-raster is stored. 
 * **\<queries_file>**: [string] -   File with a set of queries. Each line is a query with format *"posX posY"* (without ").
 * **-n \<n_reps>**: [number] -  Number of repetitions of the set of queries (def. 1).

#### 2.3.2. Obtaining all the values of a region ([get_values_window](src/k2_raster_get_values_window.cpp)) ####

Given a window range,. This query retrieves all values whose value lies within that window.

 ```bash
 ./get_values_window_k2r <k2_raster_file> <queries_file> [-c] [-n <n_reps>]>
 ```
 Where:
 
 * **\<k2_raster_file>** Path of the file where the k<sup>2</sup>-raster is stored 
 * **\<queries_file>**  File with a set of queries. Each line is a query with format:
   * *"pos_x_ini pos_x_end pos_y_ini pos_y_end "* (without ").  For a region [pos_x_ini, pos_y_ini] X [pos_x_end, pos_y_end]
 * **-c**: [boolean] - Enable checking process. (def. false).
 * **-n \<n_reps>**: [number] -  Number of repetitions of the set of queries (def. 1).

#### 2.3.3. Retrieving cells with a given value or range of values ([get_cells_by_value](src/k2_raster_get_cells_by_value.cpp)) ####

Given a range of values, this query retrieves all raster positions whose value lies within that range.

 ```bash
 ./get_cells_by_value <k2_raster_file> <queries_file> [-c] [-n <n_reps>]>
 ```
 Where:
 
 * **\<k2_raster_file>** Path of the file where the k<sup>2</sup>-raster is stored 
 * **\<queries_file>**  File with a set of queries. Each line is a query with format:
   * *"pos_x_ini pos_x_end pos_y_ini pos_y_end val_ini val_end"* (without ").  For a region [pos_x_ini, pos_y_ini] X [pos_x_end, pos_y_end] and range value [val_ini-val_end]
 * **-c**: [boolean] - Enable checking process. (def. false).
 * **-n \<n_reps>**: [number] -  Number of repetitions of the set of queries (def. 1).


#### 2.3.4. Checking the existence of a given value or range of values ([check_values_window](src/k2_raster_check_values_window.cpp)) ####

Given a region and a range, this query checks if all cell values of the region are within the range of values or if there exists at least one cell value in the region whose value lies within the range of values.

 ```bash
 ./CheckValuesk2Raster <k2_raster_file> <queries_file> [-c] [-w] [-n <n_reps>]>
 ```
 Where:
 
 * **\<k2_raster_file>** Path of the file where the k<sup>2</sup>-raster is stored 
 * **\<queries_file>**  File with a set of queries. Each line is a query with format:
     * *"pos_x_ini pos_x_end pos_y_ini pos_y_end val_ini val_end"* (without ").  For a region [pos_x_ini, pos_y_ini] X [pos_x_end, pos_y_end] and range value [val_ini-val_end]
 * **-c**: [boolean] - Enable checking process. (def. false).
 * **-w**: [boolean] -  (weak) Sets whether least one must meet query conditions (def. (strong) All cells e must meet query conditions).
 * **-n \<n_reps>**: [number] -  Number of repetitions of the set of queries (def. 1).



## 3 Temporal raster data

Soon
 
### 3.1 Encode


#### 3.1.2 Example

### 3.2 T-k<sup>2</sup>-raster types

### 3.2 Queries

### 3.2.1 Obtaining a cell value ([get_cell](src/temporal/k2_raster_temporal_get_cell.cpp))



## 4. Map algebra

### 4.1. Point Wise operation

Run a point wise operation (sum, subtraction or multiplication) between two k<sup>2</sup>-raster. 

Usage:
 ```bash
 ./algebra_k2r <raster1> <raster2> [-o <operation>] [-s <output>] [-c] [-m] [-t <type>]
 ```

Where:

* **\<raster1>**: [string] - path to the first k2-raster file.
* **\<scalar_value>**: [number] - Number to operate with each cell.
* **-o \<operation>**: [number] - map algebra operation (def. 0). Possible operations are:
  * 0 ->  Sum.
  * 1 ->  Subtraction.
  * 3 ->  Multiplication.
* **-s \<output>**: [string] - path to the output file where store the raster result.
* **-c**: [boolean] - Enable checking process. (def. false)\.
* **-m**: [boolean] - Enable the reduction of consumed memory. (def. false)\.
* **-t \<type>**: [number] - raster type (def. 10):\n"
  * 10 ->  hybrid k2-raster.
  
### 4.2. Scalar operation

Runa a scalar operation (sum, subtraction or multiplication) to one k<sup>2</sup>-raster.

Usage:
 ```bash
 ./algebra_sc_k2r <raster1> <scalar_value>[-o <operation>] [-s <output>] [-c] [-m] [-t <type>] [-r nreps]
 ```

Where:

* **\<raster1>**: [string] - path to the first k2-raster file.
* **\<raster2>**: [string] - path to the second k2-raster file.
* **-o \<operation>**: [number] - map algebra operation (def. 0). Possible operations are:
    * 0 ->  Sum.
    * 1 ->  Subtraction.
    * 3 ->  Multiplication.
* **-s \<output>**: [string] - path to the output file where store the raster result.
* **-c**: [boolean] - Enable checking process. (def. false)\.
* **-m**: [boolean] - Enable the reduction of consumed memory. (def. false)\.
* **-t \<type>**: [number] - raster type (def. 10):
    * 10 ->  hybrid k2-raster.
* **-r**: [number] - number of executions. (def. 1)\.


### 4.3. Thresholding operation

Run a thresholding operation to one k<sup>2</sup>-raster.

Usage:
 ```bash
 ./algebra_th_k2r <raster1> <thr_value>[-o <operation>] [-s <output>] [-c] [-m] [-t <type>] [-r nreps]
 ```

Where:

* **\<raster1>**: [string] - path to the first k2-raster file.
* **\<thr_value>**: [number] - Number to apply the thresholding.
* **-s \<output>**: [string] - path to the output file where store the raster result.
* **-c**: [boolean] - Enable checking process. (def. false)\.
* **-m**: [boolean] - Enable the reduction of consumed memory. (def. false)\.
* **-t \<type>**: [number] - raster type (def. 10):
    * 10 ->  hybrid k2-raster.
* **-r**: [number] - number of executions. (def. 1)\.

### 4.4. Zonal operation

Run zonal map algebra (sum, subtraction or multiplication) between two k<sup>2</sup>-raster.

Usage:
 ```bash
 ./algebra_k2r <raster1> <raster_zonal> [-o <operation>] [-s <output>] [-c] [-m] [-t <type>]
 ```

Where:

* **\<raster1>**: [string] - path to the first k2-raster file.
* **\<raster_zonal>**: [string] - path to the zonal k2-raster file.
* **-o \<operation>**: [number] - map algebra operation (def. 0). Possible operations are:
    * 0 ->  Sum.
    * 1 ->  Subtraction.
    * 3 ->  Multiplication.
* **-s \<output>**: [string] - path to the output file where store the raster result.
* **-c**: [boolean] - Enable checking process. (def. false)\.
* **-m**: [boolean] - Enable the reduction of consumed memory. (def. false)\.
* **-f**: [boolean] - Enable version 2. (def. false)\.
* **-t \<type>**: [number] - raster type (def. 10):
    * 10 ->  hybrid k2-raster.


# Publications ##

Ladra, S., Paramá, J. R., & Silva-Coira, F. (2016). Compact and queryable representation of raster datasets. In Proceedings of the 28th International Conference on Scientific and Statistical Database Management - SSDBM ’16 (pp. 1–12). New York, New York, USA: ACM Press. [doi:10.1145/2949689.2949710](http://doi.org/10.1145/2949689.2949710)

Ladra, S., Paramá, J. R., & Silva-Coira, F. (2017). Scalable and queryable compressed storage structure for raster data, Information System, 72. [doi:10.1016/j.is.2017.10.007](http://doi.org/10.1016/j.is.2017.10.007).

Cerdeira-Pena, A., de Bernardo, G., Fariña, A., Paramá, J. R., & Silva-Coira, F. (2018). Towards a compact representation of temporal rasters. In Proceedings of the 25th String Processing and Information Retrieval. SPIRE 2018 (Vol. 11147 LNCS, pp. 117–130). Springer Verlag. [doi:978-3-030-00479-8_10](https://doi.org/10.1007/978-3-030-00479-8_10)

Silva-Coira F., Paramá J. R., Ladra S., López J. R., Gutiérrez G. (2020). Efficient processing of raster and vector data. PLOS ONE 15(1): e0226943. [doi:10.1371/journal.pone.0226943](https://doi.org/10.1371/journal.pone.0226943)

Silva-Coira, F., Paramá, J. R., de Bernardo, G., & Seco, D. (2021). Space-efficient representations of raster time series. Information Sciences, 566, 300–325. [doi:10.1016/j.ins.2021.03.035](https://doi.org/10.1016/j.ins.2021.03.035)

# Web Sites ##

* [Official Web Site](https://lbd.udc.es/research/k2-raster/)
* [Code](https://gitlab.lbd.org.es/fsilva/k2-raster)


# Authors ##
* Susana Ladra - <susana.ladra@udc.es>     
* José R. Paramá - <jose.parama@udc.es>    
* Fernando Silva-Coira - <fernando.silva@udc.es>
* Guillermo de Bernardo - <guillermo.debernardo@udc.es> (T-k<sup>2</sup>-raster)

### Implemented By ###
* Fernando Silva-Coira - <fernando.silva@udc.es> (k<sup>2</sup>-raster and T-k<sup>2</sup>-raster)
* Guillermo de Bernardo - <guillermo.debernardo@udc.es> (T-k<sup>2</sup>-raster)

## License
This project is licensed under the GNU License - see the [LICENSE.md](LICENSE.md) file for details.

## Acknowledgements ##
This work was supported by European Union's Horizon 2020 research and innovation programme under the Marie Sk\l odowska-Curie [grant number 690941], Ministerio de Economía y Competitividad (PGE and FEDER) [grant numbers TIN2016-78011-C4-1-R ; TIN2016-77158-C4-3-R ; TIN2013-46238-C4-3-R ; TIN2013-46801-C4-3-R ], Centro para el desarrollo Tecnológico e Industrial [grant numbers IDI-20141259; ITC-20151247, ITC-20161074]; Xunta de Galicia (co-founded with FEDER) [grant numbers ED431C 2017/58; ED431G/01].

=========================

## ¿Cómo desplegar este proyecto en Docker?

### 1. Software necesario para desplegar el proyecto

* Git
* Docker Engine o Docker Desktop según el sistema operativo que se utilice

### 2. Clonar el repositorio

* Ejecutar el siguiente comando:
```bash
git clone https://github.com/FeMedinaRe/k2-raster
```

### 2. ¿Cómo crear la imagen Docker?

* Ejecutar el siguiente comando:
```bash
docker build -t k2raster .
```

### 3. Ejecutar el contenedor

* Ejecutar el siguiente comando:
```bash
docker run -d -p 8080:80 k2raster
```

### 4. Acceder a la visualización web

* Acceder a la IP del dispositivo donde se está ejecutando el contenedor y agregar el puerto asignado previamente:

localhost:8080