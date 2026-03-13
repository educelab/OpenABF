# Product Definition

## Project Name
OpenABF

## Description
A header-only C++ library implementing angle-based flattening algorithms for mesh parameterization.

## Problem Statement
Developers lack a lightweight, easy-to-integrate C++ implementation of ABF++ and related mesh parameterization algorithms. Existing mesh flattening tools are too heavy or require complex dependencies to embed in projects.

## Target Users
- C++ developers building geometry processing, 3D graphics, or mesh analysis applications
- Researchers in computational geometry, UV unwrapping, or virtual unwrapping of physical artifacts
- Game developers needing lightweight mesh parameterization without heavy dependencies

## Key Goals
- Provide a zero-to-low dependency, header-only ABF++ implementation that is trivial to integrate
- Support multiple flattening algorithms (ABF, ABF++, LSCM)
- Provide interfaces for multi-chart flattening and packing
- Maintain high numerical accuracy for research use
- Provide memory- and time-efficient implementations for production use
