PROGRAM cdfvrmp
  !!======================================================================
  !!                     ***  PROGRAM  cdfvrmp  ***
  !!=====================================================================
  !!  ** Purpose : Remap variables from one vertical grid to another.
  !!
  !!  ** Method  : Uses the the PPR pakcage of Engwirda & Kelley 2016
  !!               (https://doi.org/10.48550/arXiv.1606.08188), a 
  !!               Fortran-90 library designed to compute high-order 
  !!               piecewise polynomial reconstructions and conservative 
  !!               integral re-mappings on structured grids.
  !!
  !!
  !! History : 4.0  : 02/2025  : D. Bruciaferri    : Original code
  !!
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------
  USE cdfio
  USE vremap
  USE cdftools
  USE modcdfnames
  USE modutils    ! for heading

  IMPLICIT NONE

END PROGRAM cdfvrmp
  
