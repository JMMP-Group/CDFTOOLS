
#   include "../PPR/src/ppr_1d.F90"

MODULE vremap
  !!======================================================================
  !!                     ***  MODULE  vremap  ***
  !! Implement Engwirda & Kelley 2016 conservative vertical remapping 
  !! operators in CDFTOOLS.
  !!=====================================================================
  !! History : 4.0 : 09/2019 : J. Chanut      : Original code in NEMO
  !!               : 02/2025 : D. Bruciaferri : Porting to CDFTOOLS
  !!----------------------------------------------------------------------

  USE ppr_1d

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: reconstructandremap, remap_linear

CONTAINS

   SUBROUTINE reconstructandremap(ptin, phin, ptout, phout, kjpk_in, kjpk_out, kn_var)
      !!----------------------------------------------------------------------
      !!                    *** ROUTINE  reconstructandremap ***
      !!
      !! ** Purpose :   Conservative remapping of a vertical column 
      !!                from one set of layers to an other one.
      !!
      !! ** Method  :   Uses D. Engwirda Piecewise Polynomial Reconstruction library.
      !!                https://github.com/dengwirda/PPR
      !!                
      !!
      !! References :   Engwirda, Darren & Kelley, Maxwell. (2015). A WENO-type 
      !!                slope-limiter for a family of piecewise polynomial methods. 
      !!                https://arxiv.org/abs/1606.08188
      !!-----------------------------------------------------------------------
      INTEGER(KIND=4), INTENT(in   )                              ::   kjpk_in    ! Number of input levels
      INTEGER(KIND=4), INTENT(in   )                              ::   kjpk_out   ! Number of output levels
      INTEGER(KIND=4), INTENT(in   )                              ::   kn_var     ! Number of variables
      REAL(KIND=4)   , INTENT(in   ), DIMENSION(kjpk_in)          ::   phin       ! Input thicknesses
      REAL(KIND=4)   , INTENT(in   ), DIMENSION(kjpk_out)         ::   phout      ! Output thicknesses
      REAL(KIND=4)   , INTENT(in   ), DIMENSION(kjpk_in , kn_var) ::   ptin       ! Input data
      REAL(KIND=4)   , INTENT(inout), DIMENSION(kjpk_out, kn_var) ::   ptout      ! Remapped data
      !
      INTEGER(KIND=4), PARAMETER :: ndof = 1
      INTEGER(KIND=4)            :: jk, jn
      REAL(KIND=8)               ::  zwin(kjpk_in+1) ,  ztin(ndof, kn_var, kjpk_in)    ! rmap1d uses dp
      REAL(KIND=8)               :: zwout(kjpk_out+1), ztout(ndof, kn_var, kjpk_out)   ! rmap1d uses dp
      TYPE(rmap_work) :: work
      TYPE(rmap_opts) :: opts
      TYPE(rcon_ends) :: bc_l(kn_var)
      TYPE(rcon_ends) :: bc_r(kn_var)
      !!--------------------------------------------------------------------

      ! Set interfaces and input data:
      zwin(1) = 0.0
      DO jk = 2, kjpk_in + 1
         zwin(jk) = zwin(jk-1) + phin(jk-1)
      END DO

      DO jn = 1, kn_var
         DO jk = 1, kjpk_in
            ztin(ndof, jn, jk) =  ptin(jk, jn)
         END DO
      END DO

      zwout(1) = 0.0
      DO jk = 2, kjpk_out + 1
         zwout(jk) = zwout(jk-1) + phout(jk-1)
      END DO

      ! specify methods
!      opts%edge_meth = p1e_method     ! 1st-order edge interp.
!      opts%cell_meth = pcm_method
!      opts%cell_meth = plm_method     ! PLM method in cells
      opts%edge_meth = p3e_method     ! 3rd-order edge interp.
      opts%cell_meth = ppm_method     ! PPM method in cells    
!      opts%edge_meth = p5e_method     ! 5th-order edge interp.
!      opts%cell_meth = pqm_method     ! PQM method in cells

      ! limiter
!      opts%cell_lims = null_limit     ! no lim.
!      opts%cell_lims = weno_limit
      opts%cell_lims = mono_limit     ! monotone limiter   

      ! set boundary conditions
      bc_l%bcopt = bcon_loose         ! "loose" = extrapolate
      bc_r%bcopt = bcon_loose
!      bc_l%bcopt = bcon_slope        
!      bc_r%bcopt = bcon_slope

      ! init. method workspace
      CALL work%init(kjpk_in+1, kn_var, opts)

      ! remap
      CALL rmap1d(kjpk_in+1, kjpk_out+1, kn_var, ndof, &
      &           zwin, zwout, ztin, ztout,            &
      &           bc_l, bc_r, work, opts)

      ! clear method workspace
      CALL work%free()

      DO jn = 1, kn_var
         DO jk = 1, kjpk_out
            ptout(jk, jn) = ztout(1, jn, jk)
         END DO
      END DO

   END SUBROUTINE reconstructandremap
 
   SUBROUTINE remap_linear(ptin, pzin, ptout, pzout, kjpk_in, kjpk_out, kn_var)
      !!----------------------------------------------------------------------
      !!                    *** ROUTINE  remap_linear ***
      !!
      !! ** Purpose :   Linear interpolation based on input/ouputs depths
      !!
      !!-----------------------------------------------------------------------
      INTEGER(KIND=4), INTENT(in   )                              ::   kjpk_in    ! Number of input levels
      INTEGER(KIND=4), INTENT(in   )                              ::   kjpk_out   ! Number of output levels
      INTEGER(KIND=4), INTENT(in   )                              ::   kn_var     ! Number of variables
      REAL(KIND=4)   , INTENT(in   ), DIMENSION(kjpk_in)          ::   pzin       ! Input depths
      REAL(KIND=4)   , INTENT(in   ), DIMENSION(kjpk_out)         ::   pzout      ! Output depths
      REAL(KIND=4)   , INTENT(in   ), DIMENSION(kjpk_in , kn_var) ::   ptin       ! Input data
      REAL(KIND=4)   , INTENT(inout), DIMENSION(kjpk_out, kn_var) ::   ptout      ! Interpolated data
      !
      INTEGER(KIND=4)                                             :: jkin, jkout, jn
      !!--------------------------------------------------------------------
      !      
      DO jkout = 1, kjpk_out !  Loop over destination grid
         !
         IF     ( pzout(jkout) <=  pzin(  1    ) ) THEN ! Surface extrapolation 
            DO jn = 1, kn_var
! linear
!               ptout(jkout,jn) = ptin(1 ,jn) + &
!                               & (pzout(jkout) - pzin(1)) / (pzin(2)    - pzin(1)) &
!                               &                          * (ptin(2,jn) - ptin(1,jn))
               ptout(jkout,jn) = ptin(1,jn)
            END DO
         ELSEIF ( pzout(jkout) >= pzin(kjpk_in) ) THEN ! Bottom extrapolation 
            DO jn = 1, kn_var
! linear
!               ptout(jkout,jn) = ptin(kjpk_in ,jn) + &
!                               & (pzout(jkout) - pzin(kjpk_in)) / (pzin(kjpk_in)    - pzin(kjpk_in-1)) &
!                               &                                * (ptin(kjpk_in,jn) - ptin(kjpk_in-1,jn))
               ptout(jkout,jn) = ptin(kjpk_in ,jn)
            END DO
         ELSEIF ( ( pzout(jkout) > pzin(1) ).AND.( pzout(jkout) < pzin(kjpk_in) )) THEN
            DO jkin = 1, kjpk_in - 1 !  Loop over source grid
               IF ( pzout(jkout) < pzin(jkin+1) ) THEN
                  DO jn = 1, kn_var
                     ptout(jkout,jn) =  ptin(jkin,jn) + &
                                     & (pzout(jkout) - pzin(jkin)) / (pzin(jkin+1)    - pzin(jkin)) &
                                     &                             * (ptin(jkin+1,jn) - ptin(jkin,jn))
                  END DO
                  EXIT
               ENDIF
            END DO
         ENDIF
         !
      END DO

   END SUBROUTINE remap_linear

END MODULE vremap
