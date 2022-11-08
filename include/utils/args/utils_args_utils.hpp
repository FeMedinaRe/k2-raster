/*
 * Created by Fernando Silva on 14/06/22.
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

#ifndef INCLUDE_UTILS_ARGS_UTILS
#define INCLUDE_UTILS_ARGS_UTILS

// System libraries
#include <unistd.h>
#include <string>
#include <iostream>

// Own libraries

// Third libraries


//**********************************************************************//
//********************* Extract last level *****************************//
//**********************************************************************//
struct args_extract_last
{
    std::string raster;                                    // Path to raster 1
    std::string output_data;                                // Path to the output file where store the last level (binary)
};

// Print command line options
void print_usage_extract_last( char *const argv[]) {
    std::cout << "usage: " + std::string(argv[0]) + " <raster> [-s <output>]\n\n" +
                 "Extract the values of the last level of the k2-raster and they are stored in a binary file.\n" +
                 "   raster: [string] - Path to the first k2-raster file.\n" +
                 "   -s output: [string] - Path to the output file where store the last level (binary).\n";
}

// Parse arguments
void parse_args_extract_last(int argc, char *const argv[], args_extract_last &arg)
{
    int c;
    extern char *optarg;
    extern int optind;

    std::string sarg;
    while ((c = getopt(argc, argv, "s:h")) != -1)
    {
        switch (c)
        {
            case 's':
                arg.output_data.assign(optarg);
                break;
            case 'h':
                print_usage_extract_last(argv);
                exit(1);
            case '?':
                std::cout << "Unknown option. Use -h for help" << std::endl;
                print_usage_extract_last(argv);
                exit(1);
        }

    }
    // the only input parameter is the file name
    if (argc == optind + 1) {
        arg.raster.assign(argv[optind]);
    } else {
        std::cout << "Invalid number of arguments" << std::endl;
        print_usage_extract_last(argv);
        exit(-1);
    }

    /***********************/
    /* Check params values */
    /***********************/
    if (arg.raster.empty()) {
        std::cout << "Please, set the output path with -s." << std::endl;
        print_usage_extract_last(argv);
        exit(-1);
    }
}
#endif // INCLUDE_UTILS_ARGS_UTILS
