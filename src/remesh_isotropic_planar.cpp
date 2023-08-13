// Copyright (c) 2023 Anton Obukhov
// All rights reserved.
// This file is part of the Planar Isotropic Remesher project.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <CGAL/exceptions.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/region_growing.h>
#include <CGAL/Polygon_mesh_processing/surface_Delaunay_remeshing.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/mark_domain_in_triangulation.h>
#include <boost/property_map/vector_property_map.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>

#include "remesh_isotropic_planar/remesh_isotropic_planar.h"
#include "remesh_isotropic_planar/utils.h"

typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor halfedge_descriptor;

typedef CGAL::Delaunay_mesh_vertex_base_2<K_EPICK> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K_EPICK> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K_EPICK, Tds> CDT;

namespace PMP = CGAL::Polygon_mesh_processing;


void detect_planar_regions(
        const Surface_mesh &mesh,
        std::vector<std::size_t> &region_ids,
        boost::vector_property_map<CGAL::Epick::Vector_3> &region_normals,
        std::vector<Point_3> &region_offsets,
        std::size_t &nb_regions,
        double tolerance = 1e-5
) {
    region_ids.resize(num_faces(mesh), -1);

    // detect planar regions in the mesh
    nb_regions = PMP::region_growing_of_planes_on_faces(
            mesh,
            CGAL::make_random_access_property_map(region_ids),
            CGAL::parameters::
                    cosine_of_maximum_angle(1 - tolerance).
                    region_primitive_map(region_normals).
                    maximum_distance(tolerance)
    );

    // Compute the centroid for each face and add it to the region's centroid
    region_offsets.resize(nb_regions, CGAL::ORIGIN);
    std::vector<int> face_counts(nb_regions, 0);
    for (auto face : faces(mesh)) {
        std::size_t region_id = region_ids[face];
        auto vertices = vertices_around_face(halfedge(face, mesh), mesh);
        std::vector<Point_3> points;
        for (auto v_it = vertices.first; v_it != vertices.second; ++v_it) {
            points.push_back(mesh.point(*v_it));
        }
        Point_3 centroid = CGAL::centroid(points.begin(), points.end());
        region_offsets[region_id] = region_offsets[region_id] + (centroid - CGAL::ORIGIN);
        face_counts[region_id]++;
    }
    for (size_t region_id=0; region_id<nb_regions; region_id++) {
        Point_3 &p = region_offsets[region_id];
        double denom = face_counts[region_id];
        region_offsets[region_id] = Point_3(p.x() / denom, p.y() / denom, p.z() / denom);
    }
}


void extract_region_mesh(
        const Surface_mesh &mesh,
        /*const */std::vector<std::size_t> &region_ids,
        size_t region_id,
        Surface_mesh &region_mesh
) {
    std::map<Surface_mesh::Vertex_index, Surface_mesh::Vertex_index> vertex_map;

    for (auto f: faces(mesh)) {
        if (region_id == get(CGAL::make_random_access_property_map(region_ids), f)) {
            std::vector<Surface_mesh::Vertex_index> face_vertices;

            for (auto h: halfedges_around_face(halfedge(f, mesh), mesh)) {
                auto v = source(h, mesh);
                Surface_mesh::Vertex_index new_vertex;

                if (vertex_map.find(v) == vertex_map.end()) {
                    new_vertex = region_mesh.add_vertex(mesh.point(v));
                    vertex_map[v] = new_vertex;
                } else {
                    new_vertex = vertex_map[v];
                }

                face_vertices.push_back(new_vertex);
            }

            region_mesh.add_face(face_vertices);
        }
    }
}


void dump_region_mesh(
        const Surface_mesh &mesh,
        /*const */std::vector<std::size_t> &region_ids,
        size_t region_id,
        const std::string &path
) {
    Surface_mesh region_mesh;
    extract_region_mesh(mesh, region_ids, region_id, region_mesh);
    CGAL::IO::write_polygon_mesh(path, region_mesh);
}


void dump_boundaries_2d(const std::vector<std::vector<Eigen::Vector2d>> &data, const std::string &path) {
    double minX = std::numeric_limits<double>::max();
    double minY = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::min();
    double maxY = std::numeric_limits<double>::min();
    for (const auto &polyline : data) {
        for (const auto &point : polyline) {
            minX = std::min(minX, point[0]);
            minY = std::min(minY, point[1]);
            maxX = std::max(maxX, point[0]);
            maxY = std::max(maxY, point[1]);
        }
    }
    double scale = 512 * std::min(1.0 / (maxX - minX), 1.0 / (maxY - minY));
    std::ofstream outFile(path);
    outFile << "<svg viewBox=\"0 0 512 512\" xmlns='http://www.w3.org/2000/svg'>\n";
    for (const auto &polyline : data) {
        outFile << "<polyline points='";
        for (const auto &point : polyline) {
            outFile << (point[0] - minX) * scale << "," << 512 - (point[1] - minY) * scale << " ";
        }
        outFile << (polyline[0][0] - minX) * scale << "," << 512 - (polyline[0][1] - minY) * scale << " ";
        outFile << "' style='fill:none;stroke:black;stroke-width:1' />\n";
    }
    outFile << "</svg>\n";
    outFile.close();
}


void dump_boundaries_3d(
        const std::vector<std::vector<Eigen::Vector3d>> &boundaries_3d,
        const std::string &path
) {
    int num_verts = 0;
    for (const auto &polyline : boundaries_3d) {
        for (const auto &point : polyline) {
            num_verts++;
        }
    }
    std::ofstream outFile(path);
    outFile << "ply\n";
    outFile << "format ascii 1.0\n";
    outFile << "element vertex " << num_verts << "\n";
    outFile << "property float x\n";
    outFile << "property float y\n";
    outFile << "property float z\n";
    outFile << "end_header\n";
    for (const auto &polyline : boundaries_3d) {
        for (const auto &p : polyline) {
            outFile << p[0] << " " << p[1] << " " << p[2] << "\n";
        }
    }
}


void compute_region_transformations(
        const Eigen::Vector3d &normal,
        const Eigen::Vector3d &offset,
        Eigen::Matrix4d &mat_to_2d,
        Eigen::Matrix4d &mat_to_3d
) {
    Eigen::Vector3d v_target(0, 0, 1);
    Eigen::Quaterniond q = Eigen::Quaterniond::FromTwoVectors(normal, v_target);
    Eigen::Matrix3d rot_to_2d = q.toRotationMatrix();
    Eigen::Matrix3d rot_to_3d = rot_to_2d.transpose();
    mat_to_2d = Eigen::Matrix4d::Identity();
    mat_to_2d.block<3, 3>(0, 0) = rot_to_2d;
    mat_to_2d.block<3, 1>(0, 3) = -rot_to_2d * offset;
    mat_to_3d = Eigen::Matrix4d::Identity();
    mat_to_3d.block<3, 3>(0, 0) = rot_to_3d;
    mat_to_3d.block<3, 1>(0, 3) = offset;
}


Eigen::Vector3d transform_2d_to_3d(const Eigen::Matrix4d &m, const Eigen::Vector2d &p) {
    Eigen::Vector4d p4(p.x(), p.y(), 0.0, 1.0);
    Eigen::Vector4d out = m * p4;
    return out.head<3>();
}


Eigen::Vector2d transform_3d_to_2d(const Eigen::Matrix4d &m, const Eigen::Vector3d &p) {
    Eigen::Vector4d p4(p.x(), p.y(), p.z(), 1.0);
    Eigen::Vector4d out = m * p4;
    return out.head<2>();
}


void extract_boundaries_3d(
        const Surface_mesh &region_mesh,
        std::vector<std::vector<halfedge_descriptor>> &boundaries_3d
) {
    std::vector<halfedge_descriptor> border_cycles;
    PMP::extract_boundary_cycles(region_mesh, std::back_inserter(border_cycles));
    boundaries_3d.resize(0);
    for (halfedge_descriptor h : border_cycles) {
        std::vector<halfedge_descriptor> b;
        for (halfedge_descriptor hc : CGAL::halfedges_around_face(h, region_mesh)) {
            b.push_back(hc);
        }
        boundaries_3d.push_back(b);
    }
}


void convert_boundaries_3d_representations(
        const Surface_mesh &region_mesh,
        const std::vector<std::vector<halfedge_descriptor>> &boundaries_3d_halfedge,
        std::vector<std::vector<Eigen::Vector3d>> &boundaries_3d_eigen
) {
    boundaries_3d_eigen.resize(0);
    for (const auto &cycle_3d_halfedge : boundaries_3d_halfedge) {
        std::vector<Eigen::Vector3d> cycle_3d_eigen;
        for (halfedge_descriptor hc : cycle_3d_halfedge) {
            const Point_3 &p = region_mesh.point(target(hc, region_mesh));
            Eigen::Vector3d p3(p.x(), p.y(), p.z());
            cycle_3d_eigen.push_back(p3);
        }
        boundaries_3d_eigen.push_back(cycle_3d_eigen);
    }
}


void convert_boundaries_3d_to_2d(
        const std::vector<std::vector<Eigen::Vector3d>> &boundaries_3d,
        const Eigen::Matrix4d &transform,
        std::vector<std::vector<Eigen::Vector2d>> &boundaries_2d
) {
    boundaries_2d.resize(0);
    for (const auto &cycle_3d : boundaries_3d) {
        std::vector<Eigen::Vector2d> cycle_2d;
        for (const auto &p3 : cycle_3d) {
            Eigen::Vector2d p2 = transform_3d_to_2d(transform, p3);
            cycle_2d.push_back(p2);
        }
        boundaries_2d.push_back(cycle_2d);
    }
}


Eigen::Vector4d cdt_compute_aabb(const CDT& cdt) {
    double xmin = std::numeric_limits<double>::max();
    double ymin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::lowest();
    double ymax = std::numeric_limits<double>::lowest();

    for (auto vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit) {
        double x = vit->point().x();
        double y = vit->point().y();
        xmin = std::min(xmin, x);
        ymin = std::min(ymin, y);
        xmax = std::max(xmax, x);
        ymax = std::max(ymax, y);
    }

    return {xmin, ymin, xmax, ymax};
}


bool is_point_close_to_boundaries_2d(
        const std::vector<std::vector<Eigen::Vector2d>> &boundaries_2d,
        const Eigen::Vector2d &point,
        double delta
) {
    for (const auto &boundary : boundaries_2d) {
        for (std::size_t i = 0; i < boundary.size(); i++) {
            const Eigen::Vector2d &p1 = boundary[i];
            const Eigen::Vector2d &p2 = boundary[(i+1) % boundary.size()];

            double line_length_squared = (p2 - p1).squaredNorm();
            if (line_length_squared == 0.0) {
                continue;
            }

            double t = (point - p1).dot(p2 - p1) / line_length_squared;
            if (t < 0.0 || t > 1.0) {
                if ((point - p1).norm() <= delta || (point - p2).norm() <= delta) {
                    return true;
                }
            } else {
                Eigen::Vector2d projection = p1 + t * (p2 - p1);
                if ((point - projection).norm() <= delta) {
                    return true;
                }
            }
        }
    }

    return false;
}


void cdt_add_subdivide_edge(
        CDT &cdt,
        const Eigen::Vector2d &v1,
        const Eigen::Vector2d &v2,
        double diameter,
        bool add_v1_constraint,
        bool add_v2_constraint,
        CDT::Vertex_handle &vh1,
        CDT::Vertex_handle &vh2
) {
    Eigen::Vector2d v = v2 - v1;
    double len = v.squaredNorm();
    if (len == 0.0) {
        return;
    }
    len = std::sqrt(len);
    v = v / len;

    int nb_chunks = int(ceil(len / diameter));
    double chunk_len = len / nb_chunks;

    if (add_v1_constraint) {
        vh1 = cdt.insert(CDT::Point(v1.x(), v1.y()));
    }
    CDT::Vertex_handle vh_last = vh1;

    for (int c=0; c<nb_chunks; c++) {
        CDT::Vertex_handle vh_next;
        if (c == nb_chunks-1) {
            if (add_v2_constraint) {
                vh_next = cdt.insert(CDT::Point(v2.x(), v2.y()));
                vh2 = vh_next;
            } else {
                vh_next = vh2;
            }
        } else {
            Eigen::Vector2d p_next = v1 + v * ((c + 1) * chunk_len);
            vh_next = cdt.insert(CDT::Point(p_next.x(), p_next.y()));
        }
        cdt.insert_constraint(vh_last, vh_next);
        vh_last = vh_next;
    }
}


void cdt_add_subdivide_cycle(CDT &cdt, const std::vector<Eigen::Vector2d> &cycle_2d, double diameter) {
    std::vector<CDT::Vertex_handle> boundary_cdt;
    CDT::Vertex_handle vh_first, vh_last;

    for (size_t i=0; i<cycle_2d.size(); i++) {
        const Eigen::Vector2d &v1 = cycle_2d[i];
        const Eigen::Vector2d &v2 = cycle_2d[(i+1) % cycle_2d.size()];

        bool is_first = i == 0;
        bool is_last = i + 1 == cycle_2d.size();
        CDT::Vertex_handle vh1, vh2;

        if (!is_first) {
            vh1 = vh_last;
        }
        if (is_last) {
            vh2 = vh_first;
        }

        cdt_add_subdivide_edge(
                cdt,
                v1,
                v2,
                diameter,
                /*add_v1_constraint=*/is_first,
                /*add_v2_constraint=*/!is_last,
                vh1,
                vh2
        );

        if (is_first) {
            vh_first = vh1;
        }
        vh_last = vh2;
    }
}


void cdt_add_boundaries_and_meshgrid(
        CDT &cdt,
        const std::vector<std::vector<Eigen::Vector2d>> &boundaries_2d,
        double diameter,
        double rel_margin=0.3
) {
    for (const auto &boundary_2d : boundaries_2d) {
        cdt_add_subdivide_cycle(cdt, boundary_2d, diameter);
    }

    Eigen::Vector4d aabb = cdt_compute_aabb(cdt);
    double xmin = aabb[0], ymin = aabb[1], xmax = aabb[2], ymax = aabb[3];

    double l = diameter;
    double h = sqrt(3) * l / 2;
    double margin = diameter * rel_margin;
    bool offset = false;
    for (double y = ymin; y <= ymax; y += h, offset = !offset) {
        for (double x = xmin + (offset ? l / 2 : 0); x <= xmax; x += l) {
            Eigen::Vector2d p(x, y);
            if (is_point_close_to_boundaries_2d(boundaries_2d, p, margin)) {
                continue;
            }
            cdt.insert(CDT::Point(x, y));
        }
    }
}


void triangulate_region(
        const std::vector<std::vector<Eigen::Vector2d>> &boundaries_2d,
        const Eigen::Matrix4d &mat_to_3d,
        Surface_mesh &mesh_triangulated,
        double diameter
) {
    if (boundaries_2d.empty()) {
        throw std::runtime_error("Region does not have a boundary");
    }

    CDT cdt;
    cdt_add_boundaries_and_meshgrid(cdt, boundaries_2d, diameter);
    CGAL::mark_domain_in_triangulation(cdt);

    std::map<CDT::Vertex_handle, Surface_mesh::Vertex_index> vertex_map;
    for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
        if (!fit->is_in_domain()) {
            continue;
        }

        for (int i=0; i<3; i++) {
            auto vit = fit->vertex(i);
            if (vertex_map.find(vit) == vertex_map.end()) {
                auto p2 = vit->point();
                auto p3 = transform_2d_to_3d(mat_to_3d, {p2.x(), p2.y()});
                vertex_map[vit] = mesh_triangulated.add_vertex(Point_3 (p3[0], p3[1], p3[2]));
            }
        }

        auto v1 = vertex_map[fit->vertex(0)];
        auto v2 = vertex_map[fit->vertex(1)];
        auto v3 = vertex_map[fit->vertex(2)];
        mesh_triangulated.add_face(v1, v2, v3);
    }
}


void process_one_region(
        const Surface_mesh &mesh,
        /*const */std::vector<std::size_t> &region_ids,
        const boost::vector_property_map<CGAL::Epick::Vector_3> &region_normals,
        const std::vector<Point_3> &region_offsets,
        Surface_mesh &region_mesh_triangulated,
        size_t region_id,
        double max_edge_len,
        bool dump_intermediates=false,
        bool debug=false
) {
    if (debug) {
        dump_intermediates = true;
    }

    std::string path_base;
    if (dump_intermediates) {
        std::ostringstream oss;
        oss << "region." << std::setw(5) << std::setfill('0') << region_id;
        path_base = oss.str();
    }

    Surface_mesh region_mesh;
    extract_region_mesh(mesh, region_ids, region_id, region_mesh);
    if (dump_intermediates) {
        dump_region_mesh(mesh, region_ids, region_id, path_base + ".0.input.off");
    }

    std::vector<std::vector<halfedge_descriptor>> region_mesh_boundaries_3d_halfedge;
    extract_boundaries_3d(region_mesh, region_mesh_boundaries_3d_halfedge);

    std::vector<std::vector<Eigen::Vector3d>> region_mesh_boundaries_3d_eigen;
    convert_boundaries_3d_representations(
            region_mesh, region_mesh_boundaries_3d_halfedge, region_mesh_boundaries_3d_eigen);
    if (dump_intermediates) {
        dump_boundaries_3d(region_mesh_boundaries_3d_eigen, path_base + ".1.boundaries.ply");
    }

    Eigen::Vector3d region_mesh_normal(
            region_normals[region_id].x(), region_normals[region_id].y(), region_normals[region_id].z());
    Eigen::Vector3d region_mesh_offset(
            region_offsets[region_id].x(), region_offsets[region_id].y(), region_offsets[region_id].z());
    Eigen::Matrix4d region_mesh_mat_to_2d, region_mesh_mat_to_3d;
    compute_region_transformations(
            region_mesh_normal,
            region_mesh_offset,
            region_mesh_mat_to_2d,
            region_mesh_mat_to_3d
    );

    std::vector<std::vector<Eigen::Vector2d>> region_mesh_boundaries_2d;
    convert_boundaries_3d_to_2d(
            region_mesh_boundaries_3d_eigen,
            region_mesh_mat_to_2d,
            region_mesh_boundaries_2d
    );
    if (dump_intermediates) {
        dump_boundaries_2d(region_mesh_boundaries_2d, path_base + ".2.boundaries.svg");
    }

    if (debug) {
        return;
    }

    triangulate_region(
            region_mesh_boundaries_2d,
            region_mesh_mat_to_3d,
            region_mesh_triangulated,
            /*diameter=*/ max_edge_len
    );
    if (dump_intermediates) {
        CGAL::IO::write_polygon_mesh(path_base + ".3.remeshed.off", region_mesh_triangulated);
    }
}


Surface_mesh remesh_isotropic_planar(
        const Surface_mesh &mesh,
        double max_edge_len,
        double tolerance,
        bool dump_intermediates,
        bool verbose
) {
    std::vector<std::size_t> region_ids;
    boost::vector_property_map<CGAL::Epick::Vector_3> region_normals;
    std::vector<Point_3> region_offsets;
    size_t nb_regions;

    if (verbose) {
        std::cout << "Extracting planar regions..." << std::endl;
    }
    detect_planar_regions(
            mesh,
            region_ids,
            region_normals,
            region_offsets,
            nb_regions
    );

    Surface_mesh mesh_out;
    int failed_regions = 0;

    if (verbose) {
        std::cout << "Remeshing " << nb_regions << " regions..." << std::endl;
    }
    for (size_t r=0; r<nb_regions; r++) {
        if (verbose) {
            tqdm(int(r), int(nb_regions));
        }

        Surface_mesh region_mesh_triangulated;
        auto process_fn = [&] (bool debug) {
            process_one_region(
                    mesh,
                    region_ids,
                    region_normals,
                    region_offsets,
                    region_mesh_triangulated,
                    /*region_id=*/ r,
                    max_edge_len,
                    dump_intermediates,
                    debug
            );
        };

        try {
            process_fn(false);

            std::map<Surface_mesh::Vertex_index, Surface_mesh::Vertex_index> vertex_map;

            for (auto v : region_mesh_triangulated.vertices()) {
                vertex_map[v] = mesh_out.add_vertex(region_mesh_triangulated.point(v));
            }

            for (auto f : region_mesh_triangulated.faces()) {
                auto hc = region_mesh_triangulated.halfedge(f);
                std::vector<Surface_mesh::Vertex_index> face_vertices;
                for (int i = 0; i < 3; i++) {
                    face_vertices.push_back(vertex_map[region_mesh_triangulated.target(hc)]);
                    hc = region_mesh_triangulated.next(hc);
                }
                if (face_vertices.size() == 3) {
                    mesh_out.add_face(face_vertices[0], face_vertices[1], face_vertices[2]);
                }
            }
        } catch (const CDT::Intersection_of_constraints_exception &e) {
            std::cerr << "Region " << r << " could not be remeshed" << std::endl;
            failed_regions++;
#ifndef NDEBUG
            if (!dump_intermediates) {
                process_fn(true);
            }
#endif
        }
    }

    if (verbose) {
        std::cout << "Cleaning up the final mesh..." << std::endl;
    }
    mesh_out = clean_up_mesh(mesh_out, &mesh, tolerance, verbose);

    if (failed_regions > 0) {
        std::cerr << "Output mesh is missing " << failed_regions << " surfaces" << std::endl;
    }

    return mesh_out;
}
