// Copyright (c) 2023 Anton Obukhov
// All rights reserved.
// This file is part of the Planar Isotropic Remesher project.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/cluster_point_set.h>

#include "remesh_isotropic_planar/utils.h"

namespace PMP = CGAL::Polygon_mesh_processing;


void tqdm(int counter, int total) {
    double progress = 1. * counter / (total - 1);
    int bar_width = 60;
    std::cout << "[";
    int position = int(bar_width * progress);
    for (int i = 0; i < bar_width; i++) {
        if (i < position) std::cout << "=";
        else if (i == position) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << "%  ";
    if (counter == total - 1) {
        std::cout << std::endl;
    } else {
        std::cout << "\r";
    }
    std::cout.flush();
}


void merge_close_vertices(
        std::vector<Point_3> &points,
        std::vector<std::vector<std::size_t>> &polygons,
        double tolerance
) {
    using Point_set = CGAL::Point_set_3<Point_3>;

    Point_set point_set;
    for (const auto &p : points) {
        point_set.insert(p);
    }

    size_t idx_vec = 0;
    for (const Point_set::Index &idx_set : point_set) {
        const Point_3 &p_set = point_set.point(idx_set);
        const Point_3 &p_vec = points[idx_vec++];
        if (p_set != p_vec) {
            throw std::runtime_error("Indexing check failed");
        }
    }

    auto cluster_map = point_set.add_property_map<int>("cluster", -1).first;
    std::vector<std::pair<std::size_t, std::size_t> > adjacencies;
    std::size_t nb_clusters = CGAL::cluster_point_set(
            point_set,
            cluster_map,
            point_set.parameters()
                .neighbor_radius(tolerance)
                .adjacencies(std::back_inserter(adjacencies))
    );

    std::vector<Point_3> cluster_aabb_min(nb_clusters, Point_3(INFINITY, INFINITY, INFINITY));
    std::vector<Point_3> cluster_aabb_max(nb_clusters, Point_3(-INFINITY, -INFINITY, -INFINITY));
    std::vector<Point_3> cluster_centroids(nb_clusters, Point_3(0, 0, 0));
    std::vector<int> cluster_sizes(nb_clusters, 0);
    for (const Point_set::Index &idx : point_set) {
        int cluster_id = cluster_map[idx];
        const Point_3 &p = point_set.point(idx);
        cluster_aabb_min[cluster_id] = CGAL::min(cluster_aabb_min[cluster_id], p);
        cluster_aabb_max[cluster_id] = CGAL::max(cluster_aabb_max[cluster_id], p);
        cluster_centroids[cluster_id] = cluster_centroids[cluster_id] + (p - CGAL::ORIGIN);
        cluster_sizes[cluster_id]++;
    }

    std::vector<Point_3> points_new;
    std::vector<std::vector<std::size_t>> polygons_new;

    std::vector<size_t> map_cluster_to_new_idx(nb_clusters, SIZE_MAX);
    std::vector<bool> cluster_collapsible(nb_clusters, false);
    for (std::size_t i = 0; i < nb_clusters; i++) {
        cluster_centroids[i] = CGAL::ORIGIN + (cluster_centroids[i] - CGAL::ORIGIN) / cluster_sizes[i];
        auto sz = cluster_aabb_max[i] - cluster_aabb_min[i];
        if (cluster_sizes[i] >= 1 && sz.x() < tolerance && sz.y() < tolerance && sz.z() < tolerance) {
            cluster_collapsible[i] = true;
            map_cluster_to_new_idx[i] = points_new.size();
            points_new.push_back(cluster_centroids[i]);
        }
        if (cluster_sizes[i] == 0) {
            throw std::runtime_error("Empty cluster");
        }
    }

    for (const auto &polygon : polygons) {
        std::vector<std::size_t> poly_new;
        for (std::size_t idx_old : polygon) {
            int cluster_id = cluster_map[idx_old];
            size_t idx_new = map_cluster_to_new_idx[cluster_id];
            poly_new.push_back(idx_new);
        }
        polygons_new.push_back(poly_new);
    }

    points = std::move(points_new);
    polygons = std::move(polygons_new);
}


void clean_up_polygon_soup(
        std::vector<Point_3> &points,
        std::vector<std::vector<std::size_t>> &polygons,
        double tolerance,
        bool verbose
) {
    size_t num_points_before = points.size();
    size_t num_polygons_before = polygons.size();
    size_t num_points_after;
    size_t num_polygons_after;

    if (tolerance > 0) {
        merge_close_vertices(points, polygons, tolerance);

        if (verbose) {
            num_points_after = points.size();
            num_polygons_after = polygons.size();

            std::cout << "Merged: vertices (" << num_points_before << "->" << num_points_after << ") " <<
                      "faces (" << num_polygons_before << "->" << num_polygons_after << ") with tol=" <<
                      tolerance << "..." << std::endl;

            num_points_before = num_points_after;
            num_polygons_before = num_polygons_after;
        }
    }

    PMP::repair_polygon_soup(points, polygons);
    PMP::orient_polygon_soup(points, polygons);

    if (verbose) {
        num_points_after = points.size();
        num_polygons_after = polygons.size();

        std::cout << "Repaired: vertices (" << num_points_before << "->" << num_points_after << ") " <<
                     "faces (" << num_polygons_before << "->" << num_polygons_after << ")..." << std::endl;
    }
}


void read_and_repair_polygon_soup(
        const std::string &path_in,
        std::vector<Point_3> &points,
        std::vector<std::vector<std::size_t>> &polygons,
        bool verbose
) {
    points.resize(0);
    polygons.resize(0);

    if (!CGAL::IO::read_polygon_soup(path_in, points, polygons) || points.empty()) {
        std::cerr << "Error: cannot read file " << path_in << std::endl;
        exit(0);
    }

    clean_up_polygon_soup(points, polygons, 0.0, verbose);
}


Surface_mesh read_and_repair_input_or_exit(const std::string &path_in, bool verbose) {
    std::vector<Point_3> points;
    std::vector<std::vector<std::size_t> > polygons;
    read_and_repair_polygon_soup(path_in, points, polygons, verbose);

    Surface_mesh mesh;
    PMP::polygon_soup_to_polygon_mesh(points, polygons, mesh);
    return mesh;
}


Surface_mesh clean_up_mesh(const Surface_mesh &mesh, double tolerance, bool verbose) {
    std::vector<Point_3> points;
    std::vector<std::vector<std::size_t>> polygons;
    PMP::polygon_mesh_to_polygon_soup(mesh, points, polygons);

    clean_up_polygon_soup(points, polygons, tolerance, verbose);

    Surface_mesh mesh_out;
    PMP::polygon_soup_to_polygon_mesh(points, polygons, mesh_out);

    return mesh_out;
}


CGAL::Bbox_3 aabb_surface_mesh(const Surface_mesh &mesh) {
    std::vector<Point_3> points;
    points.reserve(mesh.number_of_vertices());

    for (const auto &v : mesh.vertices()) {
        points.push_back(mesh.point(v));
    }

    auto bbox = CGAL::bounding_box(points.begin(), points.end());
    return {
            bbox.xmin(),
            bbox.ymin(),
            bbox.zmin(),
            bbox.xmax(),
            bbox.ymax(),
            bbox.zmax()
    };
}


double aabb_extent(const CGAL::Bbox_3 &bbox) {
    double x_extent = bbox.xmax() - bbox.xmin();
    double y_extent = bbox.ymax() - bbox.ymin();
    double z_extent = bbox.zmax() - bbox.zmin();
    return std::max({x_extent, y_extent, z_extent});
}


double aabb_extent(const Surface_mesh &mesh) {
    return aabb_extent(aabb_surface_mesh(mesh));
}
