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

#ifndef INCLUDE_UTILS_ARGS_RASTER
#define INCLUDE_UTILS_ARGS_RASTER

// System libraries
#include <unistd.h>
#include <string>
#include "k2_raster_base.hpp"
// Own libraries
// Third libraries


//**********************************************************************//
//*************************** ENCODE ***********************************//
//**********************************************************************//
struct args_encode
{
    std::string input_data;                                 // Path to raster (binary file)
    size_t      rows;
    size_t      cols;
    std::string output_data;                                // Path to the output file where store the k2-raster
    uint        k1=2;
    uint        k2=0;
    uint        level_k1=0;
    uint        plain_levels = 0;
    bool        set_check = false;                          // Enable checking process
    ushort      type = k2raster::K2_RASTER_TYPE;            // Raster type

};

// Print command line options
void print_usage_encode( char *const argv[]) {
    std::cout << "usage: " + std::string(argv[0]) + " <input_data> <rows> <cols> <output_data> [-c] [-t <type>] [-k <k1>] [-q <k2> -l <level_k1> [-p <plain_levels>]]  \n\n" +
                 "Construct a k2-raster structure.\n" +
                 "   <input_data>: [string] - path to the raster file.\n" +
                 "   <rows>: [number] - Number of rows.\n" +
                 "   <cols>: [number] - Number of cols.\n" +
                 "   <output_data>: [string] - Path to the output file where store the k2-raster.\n" +
                 "   -c: [boolean] - Enable checking process. (def. false)\n" +
                 "   -t <type>: [number] - raster type (def. 0):\n"
                 "                           " + std::to_string(k2raster::K2_RASTER_TYPE) + " ->  hybrid k2-raster.\n" +
                 "                           " + std::to_string(k2raster::K2_RASTER_TYPE_PLAIN) + " ->  Plain k2-raster (last level).\n" +
                 "                           " + std::to_string(k2raster::K2_RASTER_TYPE_HEURISTIC) + " ->  Heuristic k2-raster.\n" +
                 "   -k <k1>: [number] - k1 value. (def. 2)\n" +
                 "   -q <k2>: [number] - k2 value . (def. 0, same value as k1).\n" +
                 "   -l <level_k1>: [number] - Number of level for k1. (def. 0).\n" +
                 "   -p <plain_levels>: [number] - Number of levels encoded  with heuristic or as plain.\n" +
                 "                           Only for types [" << std::to_string(k2raster::K2_RASTER_TYPE_PLAIN) +
                       "-" << std::to_string(k2raster::K2_RASTER_TYPE_HEURISTIC) + "].\n";


}

// Parse arguments
void parse_args_encode(int argc, char *const argv[], args_encode &arg)
{
    int c;
    extern char *optarg;
    extern int optind;

    std::string sarg;
    while ((c = getopt(argc, argv, "ct:k:q:l:p:h")) != -1)
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
            case 'k':
                sarg.assign(optarg);
                arg.k1 = stoi(sarg);
                break;
            case 'q':
                sarg.assign(optarg);
                arg.k2 = stoi(sarg);
                break;
            case 'l':
                sarg.assign(optarg);
                arg.level_k1 = stoi(sarg);
                break;
            case 'p':
                sarg.assign(optarg);
                arg.plain_levels = stoi(sarg);
                break;
            case 'h':
                print_usage_encode(argv);
                exit(1);
            case '?':
                std::cout << "Unknown option. Use -h for help" << std::endl;
                print_usage_encode(argv);
                exit(1);
        }

    }
    // the only input parameter is the file name
    if (argc == optind + 4) {
        arg.input_data.assign(argv[optind]);
        arg.rows = atol(argv[optind+1]);
        arg.cols = atol(argv[optind+2]);
        arg.output_data.assign(argv[optind+3]);
    } else {
        std::cout << "Invalid number of arguments" << std::endl;
        print_usage_encode(argv);
        exit(-1);
    }

    if (arg.k2 == 0) {
        arg.k2 = arg.k1;
    }

    /***********************/
    /* Check params values */
    /***********************/

    if (arg.k1 < 2) {
        std::cout << "[Error] - No valid value k1 =  " << arg.k2 << ". Should be greater or equal to 2" << std::endl;
        print_usage_encode(argv);
        exit(-1);
    }

    // Check if there is an operation with that code
    switch (arg.type) {
        case k2raster::K2_RASTER_TYPE:
        case k2raster::K2_RASTER_TYPE_PLAIN:
        case k2raster::K2_RASTER_TYPE_HEURISTIC:
            break;
        default:
            std::cout << "[Error] - No valid type " << arg.type << std::endl;
            print_usage_encode(argv);
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

//**********************************************************************//
//************************** SHOW INFO *********************************//
//**********************************************************************//
struct args_show_info
{
    std::string input_file;     // Path to k2-raster.
};

// Print command line options
void print_usage_show_info( char *const argv[]) {
    std::cout << "usage: " + std::string(argv[0]) + " <k2_raster_file>  \n\n" +
                 "Show info about k2-raster.\n" +
                 "   <k2_raster_file>: [string] - path to the k2-raster file.\n";
}

// Parse arguments
void parse_args_show_info(int argc, char *const argv[], args_show_info &arg)
{
    int c;
    extern char *optarg;
    extern int optind;

    std::string sarg;
    while ((c = getopt(argc, argv, "n:h")) != -1)
    {
        switch (c)
        {
            case 'h':
                print_usage_show_info(argv);
                exit(1);
            case '?':
                std::cout << "Unknown option. Use -h for help" << std::endl;
                print_usage_show_info(argv);
                exit(1);
        }

    }
    // the only input parameter is the file name
    if (argc == optind + 1) {
        arg.input_file.assign(argv[optind]);
    } else {
        std::cout << "Invalid number of arguments" << std::endl;
        print_usage_show_info(argv);
        exit(-1);
    }
}



#endif //INCLUDE_UTILS_ARGS_RASTER
