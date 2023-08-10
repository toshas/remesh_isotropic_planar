// Copyright (c) 2023 Anton Obukhov
// All rights reserved.
// This file is part of the Planar Isotropic Remesher project.
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K_EPICK;
typedef K_EPICK::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;


Surface_mesh remesh_isotropic_planar(
        const Surface_mesh &mesh,
        double max_edge_len,
        double tolerance=1e-5,
        bool dump_intermediates=false,
        bool verbose=false
);
