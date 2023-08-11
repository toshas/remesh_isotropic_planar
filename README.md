# Planar Isotropic Remesher

This code implements _almost_ isotropic remeshing of 3D meshes with planarity constraints.
It fully retains the input geometry and refines flat areas specifically.
All output _mesh edges_ will not exceed a predetermined length, and _mesh vertices_ will be evenly spaced in flat regions, preserving the original shape.

![teaser](docs/teaser.png)

## Features

- **Near-Isotropic Output**: _Mesh edges_ are subdivided up to a defined size, with input geometry constraints given precedence.
- **Single-Pass Algorithm**: Contrary to iterative pure isotropic methods, ours delivers superior performance.
- **Geometry Integrity**: Input geometry is retained without alteration; reference our [Gallery](#gallery) for examples.
- **Needle Face Elimination**: Check the teaser image, especially near the stairs, for evident reduction.
- **High-Resolution Capability**: Our method excels where pure isotropic remeshers lag.
- **Robust Handling**: Efficiently processes non-manifold and imperfect input meshes.
- **CAD Compatibility**: Especially tailored for engineered and other CAD models.

## Usage

### Building

Ensure that the dependencies are installed:
- [CMake](https://cmake.org/) build system and a C++17 toolchain
- [CGAL](https://www.cgal.org/) 3D processing library
- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) matrix computation library
- [Boost](https://www.boost.org/) C++ extensions library

Then build as follows:

```shell
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

### Running

The code is compatible with the Object File Format (`*.off`) file format.  
To test:

```shell
./build/remesh_isotropic_planar data/jacuzzi.off out.off --resolution 128 
```

Here, `resolution` indicates the ratio of the model's diameter to the output _mesh edge_ length. 
For the given test command, the edge in the output will be specifically 1/128th of the input model's dimensions.
Opting for higher resolutions will yield finer meshes with increased faces and vertices.

## Gallery

The first row presents a sample input model alongside its wireframe. 
We compare two techniques: MeshLab and our method, demonstrated at resolutions 32 and 128, respectively. 
For MeshLab, settings of the Isotropic Explicit Remeshing filter included 10 iterations, adaptive remeshing, a crease angle of 0.01, among standard options. 
The animation clearly highlights our method's superior consistency and absence of artifacts.

<img src="docs/comparison.gif"/>

## License and Citation

Copyright (c) 2023, Anton Obukhov.

The code is provided under GPL-3.0-or-later, refer to [LICENSE](LICENSE) for complete terms.

