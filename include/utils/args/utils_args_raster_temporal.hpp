/*
 * Created by Fernando Silva on 20/06/22.
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

#ifndef INCLUDE_UTILS_ARGS_RASTER_TEMPORAL
#define INCLUDE_UTILS_ARGS_RASTER_TEMPORAL

// System libraries
#include <unistd.h>
#include <string>
// Own libraries
#include <temporal/k2_raster_temporal_base.hpp>
// Third libraries


//**********************************************************************//
//*************************** ENCODE ***********************************//
//**********************************************************************//

struct args_encode_temporal
{
    std::string input_file;                                 // This file contains the parameter of configurations (first row) and the name of each raster.
    std::string input_folder;                               // This folder contains the set of raster
    std::string output_data;                                // Path to the output file where store the t-k2-raster
    ushort      type = k2raster::ATH_K2_RASTER_TYPE;        // Raster type
    uint        snap_freq=10;                               // Frequency of a snapshot in the temporal serie.
    ushort      scale_factor=0;                             // Input values are converted as: value *  10^scale_factor
    bool        set_check = false;                          // Enable checking processs
};

// Print command line options
void print_usage_encode_temporal( char *const argv[]) {
    std::cout << "usage: " + std::string(argv[0]) + " <input_file> <input_folder> <output_data> [-c] [-t <type>] [-f <snap_freq>] [-s scale_factor] \n\n" +
                 "Construct a tk2-raster structure.\n" +
                 "   <input_data>: [string] - Path to the file that contains the parameter of configurations (first row) and the name of each raster. \n" +
                 "   <input_folder>: [string] - Path to the folder that contains the set of raster. \n" +
                 "   <output_data>: [string] - Path to the output file where store the tk2-raster.\n" +
                 "   -c: [boolean] - Enable checking process. (def. false)\n" +
                 "   -t <type>: [number] - Raster type (def. 112 - Adaptive temporal k2-raster):\n"
                 "                           " + std::to_string(k2raster::K2_RASTER_TEMPORAL_TYPE) + " ->  Each raster is k2-raster.\n" +
                 "                           " + std::to_string(k2raster::K2_RASTER_TEMPORAL_TYPE_H) + " ->  Each raster is a Heuristic k2-raster.\n" +
                 "                           " + std::to_string(k2raster::T_K2_RASTER_TYPE) + " ->  Temporal k2-raster.\n" +
                 "                           " + std::to_string(k2raster::AT_K2_RASTER_TYPE) + " ->  Adaptive Heuristic k2-raster.\n" +
                 "                           " + std::to_string(k2raster::ATH_K2_RASTER_TYPE) + " ->  Adaptive Heuristic k2-raster.\n" +
                 "   -f <snap_freq>: [number] - Frequency of a snapshot in the temporal serie. (def. 10). Only for non-adaptive versions.\n" +
                 "                           Only for types [" << std::to_string(k2raster::K2_RASTER_TEMPORAL_TYPE) +
                                             "-" << std::to_string(k2raster::T_K2_RASTER_TYPE) + "].\n" +
                 "   -s <scale_factor>: [number] - Input values are converted as: value *  10^scale_factor. (def. 0).\n";
}

// Parse arguments
void parse_args_encode_temporal(int argc, char *const argv[], args_encode_temporal &arg)
{
    int c;
    extern char *optarg;
    extern int optind;

    std::string sarg;
    while ((c = getopt(argc, argv, "ct:f:s:h")) != -1)
    {
        switch (c)
        {
            case 'c':
                arg.set_check = true;
                break;
            case 't':
                sarg.assign(optarg);
                arg.type = stoi(sarg);
                break;
            case 'f':
                sarg.assign(optarg);
                arg.snap_freq = stoi(sarg);
                break;
            case 's':
                sarg.assign(optarg);
                arg.scale_factor = stoi(sarg);
                break;
            case 'h':
                print_usage_encode_temporal(argv);
                exit(1);
            case '?':
            default:
                std::cout << "Unknown option. Use -h for help" << std::endl;
                print_usage_encode_temporal(argv);
                exit(1);
        }

    }
    // the only input parameter is the file name
    if (argc == optind + 3) {
        arg.input_file.assign(argv[optind]);
        arg.input_folder.assign(argv[optind]+1);
        arg.output_data.assign(argv[optind+2]);
    } else {
        std::cout << "Invalid number of arguments" << std::endl;
        print_usage_encode_temporal(argv);
        exit(-1);
    }

    /***********************/
    /* Check params values */
    /***********************/

    switch (arg.type) {
        case k2raster::K2_RASTER_TEMPORAL_TYPE:
        case k2raster::K2_RASTER_TEMPORAL_TYPE_H:
        case k2raster::T_K2_RASTER_TYPE:
            // Check snapshot frequency
            if (arg.snap_freq < 1) {
                std::cout << "[Error] - No valid value snap_freq =  " << arg.snap_freq << ". Should be greater or equal to 1" << std::endl;
            }
            print_usage_encode_temporal(argv);
            exit(-1);
        case k2raster::AT_K2_RASTER_TYPE:
        case k2raster::ATH_K2_RASTER_TYPE:
        case k2raster::K2_RASTER_TEMPORAL_GLOBAL_TYPE:
            break;
        default:
            std::cout << "[Error] - No valid type " << arg.type << std::endl;
            print_usage_encode_temporal(argv);
            exit(-1);
    }
}

//**********************************************************************//
//************************** GET CELL **********************************//
//**********************************************************************//
struct args_get_cell
{
    std::string input_file;     // Path to k2-raster.
    std::string query_file;     // Path to the query file
    uint        n_reps = 1;     // Number of repetitions of the query set
};

// Print command line options
void print_usage_get_cell( char *const argv[]) {
    std::cout << "usage: " + std::string(argv[0]) + " <k2_raster_file> <query_file> [-n <n_reps>]  \n\n" +
                 "Query a k2-raster.\n" +
                 "   <k2_raster_file>: [string] - path to the k2-raster file.\n" +
                 "   <query_file>: [string] - path to the query file with the set of queries. Each row is a query with format 'x y'.\n" +
                 "   -n <n_reps>: [boolean] - Number of repetitions of the set of queries (def. 1)\n";
}

// Parse arguments
void parse_args_get_cell(int argc, char *const argv[], args_get_cell &arg)
{
    int c;
    extern char *optarg;
    extern int optind;

    std::string sarg;
    while ((c = getopt(argc, argv, "n:h")) != -1)
    {
        switch (c)
        {
            case 'n':
                sarg.assign(optarg);
                arg.n_reps = stoi(sarg);
                break;
            case 'h':
                print_usage_get_cell(argv);
                exit(1);
            case '?':
                std::cout << "Unknown option. Use -h for help" << std::endl;
                print_usage_get_cell(argv);
                exit(1);
        }

    }
    // the only input parameter is the file name
    if (argc == optind + 2) {
        arg.input_file.assign(argv[optind]);
        arg.query_file.assign(argv[optind+1]);
    } else {
        std::cout << "Invalid number of arguments" << std::endl;
        print_usage_get_cell(argv);
        exit(-1);
    }

    /***********************/
    /* Check params values */
    /***********************/

    if (arg.n_reps < 1) {
        std::cout << "[Error] - No valid value n_reps =  " << arg.n_reps << ". Should be greater or equal to 1" << std::endl;
        print_usage_get_cell(argv);
        exit(-1);
    }
}

//**********************************************************************//
//*********************** GET CELLS/WINDOW *****************************//
//**********************************************************************//
struct args_get_cells
{
    std::string input_file;         // Path to k2-raster.
    std::string query_file;         // Path to the query file
    bool        set_check = false;  // Enable checking process
    uint        n_reps = 1;         // Number of repetitions of the query set
};

// Print command line options
void print_usage_get_cells_or_window( char *const argv[], bool is_get_window=false) {
    std::cout << "usage: " + std::string(argv[0]) + " <k2_raster_file> <query_file> [-c] [-n <n_reps>]  \n\n";
    if (is_get_window) {
        std::cout << "Given a window range,. This query retrieves all values whose value lies within that window.\n";
    } else {
        std::cout << "Given a range of values, this query retrieves all raster positions whose value lies within that range.\n";
    }
    std::cout << "   <k2_raster_file>: [string] - path to the k2-raster file. \n";
    std::cout << "   <query_file>: [string] - path to the query file with the set of queries. Each row is a query with format 'x y'.\n";
    std::cout << "   -c: [boolean] - Enable checking process. (def. false)\n";
    std::cout << "   -n <n_reps>: [boolean] - Number of repetitions of the set of queries (def. 1)\n";
}

// Parse arguments
void parse_args_get_cells_or_window(int argc, char *const argv[], args_get_cells &arg, bool is_get_window=false)
{
    int c;
    extern char *optarg;
    extern int optind;

    std::string sarg;
    while ((c = getopt(argc, argv, "cn:h")) != -1)
    {
        switch (c)
        {
            case 'c':
                arg.set_check = true;
                break;
            case 'n':
                sarg.assign(optarg);
                arg.n_reps = stoi(sarg);
                break;
            case 'h':
                print_usage_get_cells_or_window(argv, is_get_window);
                exit(1);
            case '?':
                std::cout << "Unknown option. Use -h for help" << std::endl;
                print_usage_get_cells_or_window(argv, is_get_window);
                exit(1);
        }

    }
    // the only input parameter is the file name
    if (argc == optind + 2) {
        arg.input_file.assign(argv[optind]);
        arg.query_file.assign(argv[optind+1]);
    } else {
        std::cout << "Invalid number of arguments" << std::endl;
        print_usage_get_cells_or_window(argv, is_get_window);
        exit(-1);
    }

    /***********************/
    /* Check params values */
    /***********************/

    if (arg.n_reps < 1) {
        std::cout << "[Error] - No valid value n_reps =  " << arg.n_reps << ". Should be greater or equal to 1" << std::endl;
        print_usage_get_cells_or_window(argv, is_get_window);
        exit(-1);
    }
}


void print_help(char * argv0) {
    printf("Usage: %s <k2_raster_file> <query_file> <strong_check> <check>\n", argv0);
}

//**********************************************************************//
//*********************** GET CELLS/WINDOW *****************************//
//**********************************************************************//
struct args_check_values
{
    std::string input_file;             // Path to k2-raster.
    std::string query_file;             // Path to the query file
    bool        strong_check = true;   // Set a strong  or weak check
    bool        set_check = false;      // Enable checking process
    uint        n_reps = 1;             // Number of repetitions of the query set
};

// Print command line options
void print_usage_check_values( char *const argv[]) {
    std::cout << "usage: " + std::string(argv[0]) + " <k2_raster_file> <query_file> [-c] [-w] [-n <n_reps>]  \n\n" +
                 "Given a region and a range, this query checks if all cell values of the region are within the range of values (strong) or if there exists at least one cell value (weak) in the region whose value lies within the range of values.\n" +
                 "   <k2_raster_file>: [string] - path to the k2-raster file. \n" +
                 "   <query_file>: [string] - path to the query file with the set of queries. Each row is a query with format 'x y'.\n" +
                 "   -c: [boolean] - Enable checking process. (def. false)\n" +
                 "   -w: [boolean] - Enable weak check. (def. 'strong')\n" +
                 "   -n <n_reps>: [boolean] - Number of repetitions of the set of queries (def. 1)\n";
}

// Parse arguments
void parse_args_check_values(int argc, char *const argv[], args_check_values &arg)
{
    int c;
    extern char *optarg;
    extern int optind;

    std::string sarg;
    while ((c = getopt(argc, argv, "cwn:h")) != -1)
    {
        switch (c)
        {
            case 'c':
                arg.set_check = true;
                break;
            case '2':
                arg.strong_check = false;
                break;
            case 'n':
                sarg.assign(optarg);
                arg.n_reps = stoi(sarg);
                break;
            case 'h':
                print_usage_check_values(argv);
                exit(1);
            case '?':
                std::cout << "Unknown option. Use -h for help" << std::endl;
                print_usage_check_values(argv);
                exit(1);
        }

    }
    // the only input parameter is the file name
    if (argc == optind + 2) {
        arg.input_file.assign(argv[optind]);
        arg.query_file.assign(argv[optind+1]);
    } else {
        std::cout << "Invalid number of arguments" << std::endl;
        print_usage_check_values(argv);
        exit(-1);
    }

    /***********************/
    /* Check params values */
    /***********************/
    if (arg.n_reps < 1) {
        std::cout << "[Error] - No valid value n_reps =  " << arg.n_reps << ". Should be greater or equal to 1" << std::endl;
        print_usage_check_values(argv);
        exit(-1);
    }
}



#endif // INCLUDE_UTILS_ARGS_RASTER_TEMPORAL
