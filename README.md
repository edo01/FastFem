# FastFem

FastFem is a finite element method (FEM) library developed as part of the course "Implementation of Finite Element Methods" (MU5MAM30) at Sorbonne University. The project was created by Edoardo Carrà, Lorenzo Gentile, and Pierpaolo Marzo.

## Table of Contents

- [FastFem](#fastfem)
  - [Table of Contents](#table-of-contents)
  - [Overview](#overview)
  - [Features](#features)
  - [Installation](#installation)
    - [Requirements](#requirements)
    - [Build Instructions](#build-instructions)
  - [FastFem Usage](#fastfem-usage)
    - [FastFem linear algebra](#fastfem-linear-algebra)
      - [**FastFem matrices**](#fastfem-matrices)
      - [**MatrixTools Utility Class**](#matrixtools-utility-class)
      - [Example: Applying Dirichlet Boundary Conditions](#example-applying-dirichlet-boundary-conditions)
      - [**Iterative Solvers**](#iterative-solvers)
      - [Example: Solving a System with CGSolver](#example-solving-a-system-with-cgsolver)
  - [Testing](#testing)

## Overview

FastFem is designed for efficient FEM computations, providing support for various matrix representations, solvers, and mesh handling. It supports P1, P2, and P3 elements, and it has been used to solve linear and nonlinear problems. It includes tools for performance analysis and testing.

## Features

- Support for **Compressed Sparse Row (CSR)**, **Skylyne**, **COO** matrices.
- Multiple test modules for performance and correctness verification.
- Various finite element operations including affine transformations and mesh handling.
- Benchmarking tools for convergence analysis.

## Installation

### Requirements

Ensure you have the following dependencies installed:
- CMake
- C++17
- Python 3

### Build Instructions

1. Create the folders needed for building the project and move into it:
```sh
mkdir build && cd build
```
2. Run CMake to configure the project:
```sh
cmake ..
```
3. Build the project:
```sh
make
```
4. Execute the tests:
```sh
./[test_name]
```
NB: some tests require an additional argument to specify the number of elements in the mesh.

## FastFem Usage

An example of how to use FastFem can be found in the `test_p2`. The example demonstrates how to solve a Poisson problem using the library.

A brief overview of the usage of the main components is provided below:

### FastFem linear algebra

#### **FastFem matrices**

The FastFem library provides various sparse matrix formats optimized for finite element computations:

- **SparseMatrix**: Abstract base class for sparse matrices.

- **COOMatrix**: Stores (row, column, value) triplets for flexible assembly.
- **CSRMatrix**: Uses compressed row storage for efficient operations.
- **SymCSRMatrix**: CSR variant storing only the upper triangular part.
- **SkylineMatrix**: Optimized for banded symmetric matrices.
Each format balances memory efficiency and computational performance.

These matrices exploit some components representing their pattern:

- **CSRPattern**: Represents the sparsity pattern of a CSRMatrix.
- **SkylinePattern**: Represents the sparsity pattern of a SkylineMatrix.

This patterns have private constructors and can be created using the following static `create_from_dof_handler` method:
```cpp
// Initialize CSR pattern from DoFHandler
linalg::CSRPattern csr_pattern = linalg::CSRPattern::create_from_dof_handler(dof_handler);

// Create a CSRMatrix
linalg::CSRMatrix A(n_dofs, csr_pattern);
```

#### **MatrixTools Utility Class**

`MatrixTools` provides static utility functions for assembling and modifying matrices in finite element computations. Key functionalities include:

- **Boundary Condition Application**: Modifies system matrices to enforce Dirichlet boundary conditions.
- **Matrix Assembly**: Adds local element contributions to the global sparse matrix.
- **Vector Assembly**: Accumulates local contributions to the right-hand side vector.
- **Interpolation**: Projects functions onto the finite element space.

This utility streamlines matrix operations, ensuring efficiency in assembling and solving FEM systems.

#### Example: Applying Dirichlet Boundary Conditions

```cpp
linalg::MatrixTools::apply_homogeneous_dirichlet(A, rhs, dof_handler, 0);
```
#### **Iterative Solvers**

- **IterativeSolver**: Base class defining the solver interface.

- **CGSolver**: Implements the Conjugate Gradient method for symmetric positive-definite systems.

#### Example: Solving a System with CGSolver

```cpp
linalg::CGSolver solver(1000, 1e-12);
linalg::Vector sol = solver.solve(A, rhs);
```

## Testing

FastFem includes multiple test modules to validate performance and accuracy. These tests can be found in the `test/` directory and include:
- **Module test** for individual components.
- **Performance tests** for benchmarking execution times.
- **Poisson solver tests** to validate FEM implementations.

