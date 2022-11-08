/*
 * Created by Fernando Silva on 15/02/22.
 *
 * Copyright (C) 2016-current-year, Fernando Silva, all rights reserved.
 *
 * Author's contact: Fernando Silva  <fernando.silva@udc.es>
 * Databases Lab, University of A Coruña. Campus de Elviña s/n. Spain
 *
 * Library to set and check command line arguments
 *
 * k2-raster is a compact data structure to represent raster data that
 * uses compressed space and offers indexing capabilities.
 * It uses min/max values for indexing and improving query performance.
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

#ifndef INCLUDE_UTILS_UTILS_ARGS
#define INCLUDE_UTILS_UTILS_ARGS

// System libraries
#include <unistd.h>
// Own libraries
// Third libraries


//**********************************************************************//
//************************** ALGEBRA ***********************************//
//**********************************************************************//
struct args_algebra
{
    std::string raster1;                                    // Path to raster 1
    std::string raster2;                                    // Path to raster 2
    std::string output_data;                                // Path to the output file where store the raster result
    ushort      operation=0;                                // Operation (number)
    bool        set_check = false;                          // Enable checking process
    ushort      k2_raster_type = k2raster::K2_RASTER_TYPE;  // Raster type (remove)
};

// Print command line options
void print_usage_algebra( char *const argv[]) {
    std::cout << "usage: " + std::string(argv[0]) + " <raster1> <raster2> [-o <operation>] [-s <output>] [-c] [-t <type>]  \n\n" +
                 "Run map algebra between two k2-raster.\n" +
                 "   raster1: [string] - path to the first k2-raster file\n" +
                 "   raster2: [string] - path to the second k2-raster file\n" +
                 "   -o operation: [number] - map algebra operation (def. 0):\n" +
                 "                           " + std::to_string(k2raster::OperationRaster::OPERATION_SUM) + " ->  Sum\n" +
                 "                           " + std::to_string(k2raster::OperationRaster::OPERATION_SUBT) + " ->  Subtraction\n" +
                 "                           " + std::to_string(k2raster::OperationRaster::OPERATION_MULT) + " ->  Multiplication\n" +
                 "   -s output: [string] - path to the output file where store the raster result\n" +
                 "   -c: [boolean] - enable checking process. (def. false)\n" +
                 "   -t type: [number] - raster type (def. 0):\n"
                 "                           " + std::to_string(k2raster::K2_RASTER_TYPE) + " ->  hybrid k2-raster.\n";
}

// Parse arguments
void parse_args_algebra(int argc, char *const argv[], args_algebra &arg)
{
    int c;
    extern char *optarg;
    extern int optind;

    std::string sarg;
    while ((c = getopt(argc, argv, "o:s:ct:h")) != -1)
    {
        switch (c)
        {
            case 'o':
                sarg.assign(optarg);
                arg.operation = stoi(sarg);
                break;
            case 's':
                arg.output_data.assign(optarg);
                break;
            case 'c':
                arg.set_check = true;
                break;
            case 't':
                sarg.assign(optarg);
                arg.k2_raster_type = stoi(sarg);
                break;
            case 'h':
                print_usage_algebra(argv);
                exit(1);
            case '?':
                std::cout << "Unknown option. Use -h for help" << std::endl;
                print_usage_algebra(argv);
                exit(1);
        }

    }
    // the only input parameter is the file name
    if (argc == optind + 2) {
        arg.raster1.assign(argv[optind]);
        arg.raster2.assign(argv[optind+1]);
    } else {
        std::cout << "Invalid number of arguments" << std::endl;
        print_usage_algebra(argv);
        exit(-1);
    }

    /***********************/
    /* Check params values */
    /***********************/

    // Check if there is an operation with that code
    switch (arg.operation) {
        case k2raster::OperationRaster::OPERATION_SUM:
        case k2raster::OperationRaster::OPERATION_SUBT:
        case k2raster::OperationRaster::OPERATION_MULT:
            break;
        default:
            std::cout << "[Error] - No valid operation " << arg.operation << std::endl;
            print_usage_algebra(argv);
            exit(-1);
    }
}

//**********************************************************************//
//*************************** SCALAR ***********************************//
//**********************************************************************//
struct args_algebra_scalar
{
    std::string raster1;                                    // Path to raster 1
    int         scalar_value;                               // Integer value
    std::string output_data;                                // Path to the output file where store the raster result
    ushort      operation=0;                                // Operation (number)
    bool        set_check = false;                          // Enable checking process
    ushort      k2_raster_type = k2raster::K2_RASTER_TYPE;  // Raster type (remove)
    uint        n_reps=1;                                   // Number of executions
};

// Print command line options
void print_usage_algebra_scalar( char *const argv[]) {
    std::cout << "usage: " + std::string(argv[0]) + " raster1 scalar_value [-o operation] [-s output] [-c] [-t type] [-r nreps] \n\n" +
                 "Run a scalar operation to one k2-raster.\n" +
                 "   raster1: [string] - path to the first k2-raster file\n" +
                 "   scalar_value: [number] - Number to operate with each cell\n" +
                 "   -o operation: [number] - map algebra operation (def. 0):\n" +
                 "                           " + std::to_string(k2raster::OperationRaster::OPERATION_SUM) + " ->  Sum\n" +
                 "                           " + std::to_string(k2raster::OperationRaster::OPERATION_SUBT) + " ->  Subtraction\n" +
                 "                           " + std::to_string(k2raster::OperationRaster::OPERATION_MULT) + " ->  Multiplication\n" +
                 "   -s output: [string] - path to the output file where store the raster result\n" +
                 "   -c: [boolean] - enable checking process. (def. false)\n" +
                 "   -t type: [number] - raster type (def. 0):\n"
                 "                           " + std::to_string(k2raster::K2_RASTER_TYPE) + " ->  hybrid k2-raster.\n" +
                 "   -r: [number] - number of executions. (def. 1)\n";
}

// Parse arguments
void parse_args_algebra_scalar(int argc, char *const argv[], args_algebra_scalar &arg)
{
    int c;
    extern char *optarg;
    extern int optind;

    std::string sarg;
    while ((c = getopt(argc, argv, "o:s:ct:r:h")) != -1)
    {
        switch (c)
        {
            case 'o':
                sarg.assign(optarg);
                arg.operation = stoi(sarg);
                break;
            case 's':
                arg.output_data.assign(optarg);
                break;
            case 'c':
                arg.set_check = true;
                break;
            case 't':
                sarg.assign(optarg);
                arg.k2_raster_type = stoi(sarg);
                break;
            case 'r':
                sarg.assign(optarg);
                arg.n_reps = stoi(sarg);
                break;
            case 'h':
                print_usage_algebra_scalar(argv);
                exit(1);
            case '?':
                std::cout << "Unknown option. Use -h for help" << std::endl;
                print_usage_algebra_scalar(argv);
                exit(1);
        }

    }
    // the only input parameter is the file name
    if (argc == optind + 2) {
        arg.raster1.assign(argv[optind]);
        arg.scalar_value = (atoi(argv[optind+1]));
    } else {
        std::cout << "Invalid number of arguments" << std::endl;
        print_usage_algebra_scalar(argv);
        exit(-1);
    }

    /***********************/
    /* Check params values */
    /***********************/

    // Check if there is an operation with that code
    switch (arg.operation) {
        case k2raster::OperationRaster::OPERATION_SUM:
        case k2raster::OperationRaster::OPERATION_SUBT:
        case k2raster::OperationRaster::OPERATION_MULT:
            break;
        default:
            std::cout << "[Error] - No valid operation " << arg.operation << std::endl;
            print_usage_algebra_scalar(argv);
            exit(-1);
    }

    if (arg.n_reps < 1) {
        std::cout << "Invalid number of executions (-r). Should be 1 or greater" << std::endl;
        print_usage_algebra_scalar(argv);
        exit(-1);
    }
}

//**********************************************************************//
//*********************** Thresholding *********************************//
//**********************************************************************//
struct args_algebra_thr
{
    std::string raster1;                                    // Path to raster 1
    int         thr_value;                                  // Integer value
    std::string output_data;                                // Path to the output file where store the raster result
    bool        set_check = false;                          // Enable checking process
    ushort      k2_raster_type = k2raster::K2_RASTER_TYPE;  // Raster type (remove)
    uint        n_reps=1;                                   // Number of executions
};

// Print command line options
void print_usage_algebra_thr( char *const argv[]) {
    std::cout << "usage: " + std::string(argv[0]) + " raster1 thr_value [-o operation] [-s output] [-c] [-t type] [-r nreps] \n\n" +
                 "Run a thresholding operation to one k2-raster.\n" +
                 "   raster1: [string] - path to the first k2-raster file\n" +
                 "   thr_value: [number] - Number to apply the thresholding\n" +
                 "   -s output: [string] - path to the output file where store the raster result\n" +
                 "   -c: [boolean] - enable checking process. (def. false)\n" +
                 "   -t type: [number] - raster type (def. 0):\n"
                 "                           " + std::to_string(k2raster::K2_RASTER_TYPE) + " ->  hybrid k2-raster.\n" +
                 "   -r: [number] - number of executions. (def. 1)\n";

}

// Parse arguments
void parse_args_algebra_thr(int argc, char *const argv[], args_algebra_thr &arg)
{
    int c;
    extern char *optarg;
    extern int optind;

    std::string sarg;
    while ((c = getopt(argc, argv, "s:ct:r:h")) != -1)
    {
        switch (c)
        {
            case 's':
                arg.output_data.assign(optarg);
                break;
            case 'c':
                arg.set_check = true;
                break;
            case 't':
                sarg.assign(optarg);
                arg.k2_raster_type = stoi(sarg);
                break;
            case 'r':
                sarg.assign(optarg);
                arg.n_reps = stoi(sarg);
                break;
            case 'h':
                print_usage_algebra_thr(argv);
                exit(1);
            case '?':
                std::cout << "Unknown option. Use -h for help" << std::endl;
                print_usage_algebra_thr(argv);
                exit(1);
        }

    }
    // the only input parameter is the file name
    if (argc == optind + 2) {
        arg.raster1.assign(argv[optind]);
        arg.thr_value = (atoi(argv[optind+1]));
    }
    else {
        std::cout << "Invalid number of arguments" << std::endl;
        print_usage_algebra_thr(argv);
        exit(-1);
    }

    /*********************/
    /* Check params      */
    /*********************/
    if (arg.n_reps < 1) {
        std::cout << "Invalid number of executions (-r). Should be 1 or greater" << std::endl;
        print_usage_algebra_thr(argv);
        exit(-1);
    }
}

//**********************************************************************//
//*************************** ZONAL ************************************//
//**********************************************************************//
struct args_algebra_zonal
{
    std::string raster1;                                    // Path to raster 1
    std::string raster_zonal;                               // Path to zonal raster
    std::string output_data;                                // Path to the output file where store the raster result
    ushort      operation=0;                                // Operation (number)
    bool        set_check = false;                          // Enable checking process
    ushort      k2_raster_type = k2raster::K2_RASTER_TYPE;  // Raster type (remove)
};

// Print command line options
void print_usage_algebra_zonal( char *const argv[]) {
    std::cout << "usage: " + std::string(argv[0]) + " raster1 raster_zonal [-o operation] [-s output] [-c] [-t type]  \n\n" +
                 "Run zonal map algebra between two k2-raster.\n" +
                 "   raster1: [string] - path to the first k2-raster file\n" +
                 "   raster2: [string] - path to the second k2-raster file\n" +
                 "   -o operation: [number] - map algebra operation (def. 0):\n" +
                 "                           " + std::to_string(k2raster::OperationRaster::OPERATION_SUM) + " ->  Sum\n" +
                 "                           " + std::to_string(k2raster::OperationRaster::OPERATION_SUBT) + " ->  Subtraction\n" +
                 "                           " + std::to_string(k2raster::OperationRaster::OPERATION_MULT) + " ->  Multiplication\n" +
                 "   -s output: [string] - path to the output file where store the raster result\n" +
                 "   -c: [boolean] - enable checking process. (def. false)\n" +
                 "   -t type: [number] - raster type (def. 0):\n"
                 "                           " + std::to_string(k2raster::K2_RASTER_TYPE) + " ->  hybrid k2-raster.\n";
}

// Parse arguments
void parse_args_algebra_zonal(int argc, char *const argv[], args_algebra_zonal &arg)
{
    int c;
    extern char *optarg;
    extern int optind;

    std::string sarg;
    while ((c = getopt(argc, argv, "o:s:ct:h")) != -1)
    {
        switch (c)
        {
            case 'o':
                sarg.assign(optarg);
                arg.operation = stoi(sarg);
                break;
            case 's':
                arg.output_data.assign(optarg);
                break;
            case 'c':
                arg.set_check = true;
                break;
            case 't':
                sarg.assign(optarg);
                arg.k2_raster_type = stoi(sarg);
                break;
            case 'h':
                print_usage_algebra_zonal(argv);
                exit(1);
            case '?':
                std::cout << "Unknown option. Use -h for help" << std::endl;
                print_usage_algebra_zonal(argv);
                exit(1);
        }

    }
    // the only input parameter is the file name
    if (argc == optind + 2) {
        arg.raster1.assign(argv[optind]);
        arg.raster_zonal.assign(argv[optind+1]);
    } else {
        std::cout << "Invalid number of arguments" << std::endl;
        print_usage_algebra_zonal(argv);
        exit(-1);
    }

    /***********************/
    /* Check params values */
    /***********************/

    // Check if there is an operation with that code
    switch (arg.operation) {
        case k2raster::OperationZonalRaster::OPERATION_ZONAL_SUM:
            break;
        default:
            std::cout << "[Error] - No valid operation " << arg.operation << std::endl;
            print_usage_algebra_zonal(argv);
            exit(-1);
    }
}
#endif //INCLUDE_UTILS_UTILS_ARGS
