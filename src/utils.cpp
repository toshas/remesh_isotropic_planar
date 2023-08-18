// Copyright (c) 2023 Anton Obukhov
// All rights reserved.
// This file is part of the Planar Isotropic Remesher project.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup_extension.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
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


void check_soup(const std::vector<Point_3> &points, const std::vector<std::vector<std::size_t>> &polygons) {
    std::size_t num_pt = points.size();
    for (const auto &poly : polygons) {
        if (poly.size() != 3) {
            throw std::runtime_error("Not a triangle");
        }
        for (std::size_t idx : poly) {
            if (idx == (std::size_t)(-1)) {
                throw std::runtime_error("Bad index (-1)");
            } else if (idx >= num_pt) {
                throw std::runtime_error(std::string("Bad index greater than ") + std::to_string(num_pt));
            }
        }
    }
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
    std::map<int, std::set<std::size_t>> map_cluster_to_idx_old;
    for (const Point_set::Index &idx : point_set) {
        int cluster_id = cluster_map[idx];
        if (cluster_id == -1) {
            throw std::runtime_error("Found a point without a cluster assigned");
        }

        map_cluster_to_idx_old[cluster_id].insert(idx);
        cluster_sizes[cluster_id]++;

        const Point_3 &p = point_set.point(idx);
        cluster_aabb_min[cluster_id] = CGAL::min(cluster_aabb_min[cluster_id], p);
        cluster_aabb_max[cluster_id] = CGAL::max(cluster_aabb_max[cluster_id], p);
        cluster_centroids[cluster_id] = cluster_centroids[cluster_id] + (p - CGAL::ORIGIN);
    }

    std::vector<Point_3> points_new;
    std::vector<size_t> map_vertex_idx_old_to_idx_new(points.size(), SIZE_MAX);
    for (std::size_t i = 0; i < nb_clusters; i++) {
        if (cluster_sizes[i] == 0) {
            throw std::runtime_error("Empty cluster");
        }
        auto sz = cluster_aabb_max[i] - cluster_aabb_min[i];
        auto all_cluster_idx_old = map_cluster_to_idx_old[int(i)];
        if (cluster_sizes[i] >= 1 && sz.x() < tolerance && sz.y() < tolerance && sz.z() < tolerance) {
            size_t idx_new = points_new.size();
            auto centroid = CGAL::ORIGIN + (cluster_centroids[i] - CGAL::ORIGIN) / cluster_sizes[i];
            points_new.push_back(centroid);
            for (const auto &idx_old : all_cluster_idx_old) {
                map_vertex_idx_old_to_idx_new[idx_old] = idx_new;
            }
        } else {
            for (const auto &idx_old : all_cluster_idx_old) {
                size_t idx_new = points_new.size();
                const Point_3 &p = point_set.point(idx_old);
                points_new.push_back(p);
                map_vertex_idx_old_to_idx_new[idx_old] = idx_new;
            }
        }
    }

    std::vector<std::vector<std::size_t>> polygons_new;
    for (const auto &polygon : polygons) {
        if (polygon.size() != 3) {
            throw std::runtime_error("Not a triangle");
        }

        std::vector<std::size_t> poly_new(polygon.size(), SIZE_MAX);
        std::vector<Point_3> verts(3);

        std::size_t i = 0;
        for (const auto &idx_old : polygon) {
            const auto &idx_new = map_vertex_idx_old_to_idx_new[idx_old];
            verts[i] = points_new[idx_new];
            poly_new[i] = idx_new;
            i++;
        }

        K_EPICK::Triangle_3 triangle(verts[0], verts[1], verts[2]);
        if (triangle.squared_area() == 0.0) {
            continue;
        }

        polygons_new.push_back(poly_new);
    }

    points = std::move(points_new);
    polygons = std::move(polygons_new);

    check_soup(points, polygons);
}


void clean_up_polygon_soup(
        std::vector<Point_3> &points,
        std::vector<std::vector<std::size_t>> &polygons,
        const Surface_mesh *reference,
        double tolerance,
        bool verbose
) {
    size_t num_points_before = points.size();
    size_t num_polygons_before = polygons.size();
    size_t num_points_after;
    size_t num_polygons_after;

    merge_close_vertices(points, polygons, tolerance);

    PMP::repair_polygon_soup(points, polygons);
    if (reference == nullptr) {
        PMP::orient_polygon_soup(points, polygons);
    } else {
        PMP::orient_triangle_soup_with_reference_triangle_mesh(*reference, points, polygons);
        PMP::duplicate_non_manifold_edges_in_polygon_soup(points, polygons);
    }

    if (verbose) {
        num_points_after = points.size();
        num_polygons_after = polygons.size();

        std::cout << "Mesh repair: vertices (" << num_points_before << "->" << num_points_after << ") " <<
                     "faces (" << num_polygons_before << "->" << num_polygons_after << ") with tolerance=" << tolerance <<
                     "..." << std::endl;
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

    clean_up_polygon_soup(points, polygons, nullptr, 0.0, verbose);
}


Surface_mesh read_and_repair_input_or_exit(const std::string &path_in, bool verbose) {
    std::vector<Point_3> points;
    std::vector<std::vector<std::size_t> > polygons;
    read_and_repair_polygon_soup(path_in, points, polygons, verbose);

    Surface_mesh mesh;
    PMP::polygon_soup_to_polygon_mesh(points, polygons, mesh);
    return mesh;
}


Surface_mesh clean_up_mesh(const Surface_mesh &mesh, const Surface_mesh *reference, double tolerance, bool verbose) {
    std::vector<Point_3> points;
    std::vector<std::vector<std::size_t>> polygons;
    PMP::polygon_mesh_to_polygon_soup(mesh, points, polygons);

    clean_up_polygon_soup(points, polygons, reference, tolerance, verbose);

    Surface_mesh mesh_out;
    PMP::polygon_soup_to_polygon_mesh(points, polygons, mesh_out);
    PMP::remove_isolated_vertices(mesh_out);

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
