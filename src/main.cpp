// Copyright (c) 2023 Anton Obukhov
// All rights reserved.
// This file is part of the Planar Isotropic Remesher project.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <boost/program_options.hpp>
#include <string>
#include <filesystem>

#include "remesh_isotropic_planar/utils.h"
#include "remesh_isotropic_planar/remesh_isotropic_planar.h"

namespace fs = std::filesystem;
namespace PO = boost::program_options;
namespace PMP = CGAL::Polygon_mesh_processing;


void parse_cli_args(
        int argc,
        const char *argv[],
        std::string &path_in,
        std::string &path_out,
        size_t &resolution,
        bool &verbose,
        bool &dump_intermediates
) {
    try {
        PO::options_description desc("Allowed options");
        desc.add_options()
                ("help,h", "Produce help message")
                ("input", PO::value<std::string>()->required(), "Input file *.off (required)")
                ("output", PO::value<std::string>()->required(), "Output file *.off (required)")
                ("resolution,r", PO::value<size_t>()->default_value(128), "Resolution")
                ("verbose,v", PO::value<bool>()->default_value(true), "Verbose")
                ("dump,d", PO::value<bool>()->default_value(false), "Dump intermediates")
                ;

        PO::positional_options_description pos;
        pos.add("input", 1);
        pos.add("output", 1);

        PO::variables_map vm;
        PO::store(PO::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
        PO::notify(vm);

        if (vm.count("help") > 0) {
            std::cout << "Planar Isotropic Remesher" << std::endl <<
                         "Command line syntax: <input.off> <output.off> [--resolution <int>]" << std::endl;
            exit(0);
        }

        path_in = vm["input"].as<std::string>();
        path_out = vm["output"].as<std::string>();
        resolution = vm["resolution"].as<size_t>();
        verbose = vm["verbose"].as<bool>();
        dump_intermediates = vm["dump"].as<bool>();
    } catch (const PO::error &e) {
        std::cerr << "Error: " << e.what() << "\n";
        exit(0);
    } catch (const std::exception &e) {
        std::cerr << "Unhandled Exception: " << e.what() << "\n";
        exit(0);
    }
}


int main(int argc, const char *argv[]) {
    std::string path_in, path_out;
    size_t resolution;
    bool verbose, dump_intermediates;
    parse_cli_args(argc, argv, path_in, path_out, resolution, verbose, dump_intermediates);

    if (verbose) {
        std::cout << "Loading from " << path_in << "..." << std::endl;
    }
    Surface_mesh mesh = read_and_repair_input_or_exit(path_in, verbose);

    double mesh_extent = aabb_extent(mesh);
    double max_edge_len = mesh_extent / double(resolution);
    double tolerance_meshlab_default = std::sqrt(3) * mesh_extent / 100000;
    double tolerance = std::min(std::max(tolerance_meshlab_default, max_edge_len / 100), max_edge_len / 3);
    if (verbose) {
        std::cout << "Inferred max_edge_len=" << max_edge_len << " and tolerance=" << tolerance <<
                " from resolution=" << resolution << "..." << std::endl;
    }

    Surface_mesh mesh_out = remesh_isotropic_planar(mesh, max_edge_len, tolerance, dump_intermediates, verbose);

    if (verbose) {
        std::cout << "Storing to " << path_out << "..." << std::endl;
    }
    fs::path path_dirs = fs::path(path_out).parent_path();
    fs::create_directories(path_dirs);
    CGAL::IO::write_polygon_mesh(path_out, mesh_out);

    if (verbose) {
        std::cout << "Done" << std::endl;
    }
}
