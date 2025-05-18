# Pose Transfer for 3D Objects

## Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should find and build the dependencies and create a `example` binary.

## Run

From within the `build` directory just issue:

    ./example <sourceMesh> <targetMesh>

A glfw app should launch displaying 2 meshes - 2 of the source and 1 for the target.

## Instructions

There are 3 major steps in Pose transfer.

### 1. Correspondence establishment

Map points on 1st mesh to points on 3rd mesh.

### 2. Correspondence verification

Click points on the 1st mesh to view corresponding points on 3rd mesh.

### 3. Deformation transfer

Deform 2nd mesh and select option 5 in the menu to transfer deformation to target.

## Dependencies

The only dependencies are STL, Eigen, [libigl](http://libigl.github.io/libigl/) and the dependencies
of the `igl::opengl::glfw::Viewer` (OpenGL, glad and GLFW).

The CMake build system will automatically download libigl and its dependencies using
[CMake FetchContent](https://cmake.org/cmake/help/latest/module/FetchContent.html),
thus requiring no setup on your part.
