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

## [F,Pnm,mP,dP] = my_SHgen(m,n,s,varargin)
Vectorial spherical harmonics generator, 
Inputs
   * n: (scalar) Degree 
   * m: (scalar) Order
   * s: (scalar) Function (1: odd <-> TE, 2: even <-> TM)
   * varargin: 
       - theta: vector
       - phi: vector
       - test: boolean, enables the comparison between the calculated Functions and direct and (selected) theoretical expressions taken from Hansen book p.322-324, for n == 2 | n == 3 | n == 4 
Outputs 
   * F: Vectorial Spherical Harmonic (struct with fields theta and phi)
   * Pnm: normalized associated Legendre function
   * mP: normalized associated Legendre function, times m, over sin\theta
   * dP: derivative of normalized associated Legendre function

## plotPattern(E,theta,phi,varargin)
Flexible routine for Radiation pattern plot, in a given coordinate system (1D: Rectangular, polar, 2D: Rectangular grid, Projected, 3D: Polar)
Inputs:
   * E:                (Array), Supported sizes: N_theta x N_Phi x N_Pats, N_theta x N_Phi, N_theta x N_Pats, N_theta x 1, N_phi x 1
   * theta:            (Array), 1D or 2D
   * phi:              (Array), 1D or 2D
   * varargin:     
       - norma:        (Logic), true if normalized plot is desired
       - scale:        (String), indicate the field axis scaling, 'dB', 'lin'
       - plStyle:      (String), representation of the input field, 'asis' -> E, 'fld' -> abs(E), 'pwr' -> abs(E).^2, 'dir' -> dir*abs(E)/max(abs(E))
       - coords:       (String), format of the axes, 1D fields: 'rect','polar', 2D fields: 'rect','proj2D','pol3D'
                           'rect' for 3D patterns creates a projection on a rectangular (theta, phi) grid
                           'proj2D' for 3D patterns creates a projection on a Direction cosine coordinate system (u=sin(th)cos(ph),v=sin(th)sin(ph)).
       - limsScale:    (Vector), [minScale,maxScale]
       - thetaCut:     (scalar), -1 for sweep, other value for an specific angle cut (3D patterns)
       - phiCut:       (scalar), -1 for sweep, other value for an specific angle cut (3D patterns)
       - multPatts:	(logic), true to indicate that multiple patterns are provided, false by default.
       - aggrPatts:    (logic), true to indicate that patterns should be added before plotting, false by default. 
Outputs:
   * No variable neither file output is generated, just the pattern plot
