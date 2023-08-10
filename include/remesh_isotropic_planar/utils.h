// Copyright (c) 2023 Anton Obukhov
// All rights reserved.
// This file is part of the Planar Isotropic Remesher project.
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "remesh_isotropic_planar/remesh_isotropic_planar.h"


void tqdm(int counter, int total);
double aabb_extent(const Surface_mesh &mesh);
Surface_mesh read_and_repair_input_or_exit(const std::string &path_in, bool verbose=false);
Surface_mesh clean_up_mesh(const Surface_mesh &mesh, double tolerance, bool verbose=false);
