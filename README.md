# SWE_basics
Basic routines for Spherical Wave modelling in the Far-Field

# Main functions: 

## q = myVc_swe(F,J,varargin)
Vector field Series expansion in terms of Vector spherical harmonics
Inputs
   * F: Vector Field to expand. (Struct with fields theta and phi)
   * J: Number of harmonics, J = 2N(N+2)
   * varargin: theta, phi (rad)
Outputs
   * q: Vector with odd (s=1) and even (s=2) expansion coefficients in alternate fashion

## F = myVc_sws(q,varargin)
Vectorial farfield reconstruction from a Spherical wave expansion
Synthesis, given the coefficients of a series expansion in spherical harmonics
Inputs
   * q: Vector of coefficients
   * varargin: theta, phi
Outputs:
   * Reconstructed field. (Struct with fields theta and phi)
