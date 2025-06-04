PROGRAM cdfsaltc
  !!======================================================================
  !!                     ***  PROGRAM  cdfsaltc  ***
  !!=====================================================================
  !!  ** Purpose : Compute the salt content of the ocean : 1 single value
  !!
  !!  ** Method  : compute the sum ( rho Abs_Sal  * e1t *e2t * e3t * tmask )
  !!               N.B: the correct way of computing the salt content is
  !!                    using the Absolute Salinity, NOT the Practical Salinity.
  !!
  !! History : 4.0  : 06/2025  : D. Bruciaferri : Original code
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class integration
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4), PARAMETER                :: jp_sc3d=1, jp_sc2d=2 , jp_scvol=3
  INTEGER(KIND=4)                           :: jk, jt              ! dummy loop index
  INTEGER(KIND=4)                           :: ik, it              ! working integer
  INTEGER(KIND=4)                           :: ierr                ! working integer
  INTEGER(KIND=4)                           :: iimin=0, iimax=0    ! domain limitation for computation
  INTEGER(KIND=4)                           :: ijmin=0, ijmax=0    ! domain limitation for computation
  INTEGER(KIND=4)                           :: ikmin=0, ikmax=0    ! domain limitation for computation
  INTEGER(KIND=4)                           :: mxloption=0         ! mixed layer option    
  INTEGER(KIND=4)                           :: narg, iargc, ijarg  ! command line 
  INTEGER(KIND=4)                           :: npiglo, npjglo      ! size of the domain
  INTEGER(KIND=4)                           :: npk, npkk,npt       ! size of the domain
  INTEGER(KIND=4)                           :: nvpk                ! vertical levels in working variable
  INTEGER(KIND=4)                           :: ncout
  INTEGER(KIND=4), DIMENSION(:),  ALLOCATABLE  :: ipk, id_varout   ! for output variables

  REAL(KIND=4), PARAMETER                   :: pprho0=1020.        ! water density (kg/m3)

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e1t, e2t            ! horizontal metrics
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e3t                 ! vertical metric
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: salin               ! absolute salinity
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: tmask               ! tmask
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rmxldep             ! mixed layer depth
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: gdepw               ! depth
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: e31d                ! vertical metrics in case of full step

  REAL(KIND=8)                              :: dvol                ! 3D volume of the ocean
  REAL(KIND=8)                              :: dsum                ! weighted sum 3D
  REAL(KIND=8)                              :: dvol2d              ! volume of a layer
  REAL(KIND=8)                              :: dsum2d              ! weigthed sum per layer
  REAL(KIND=8)                              :: dsurf               ! surface of a layer
  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dtim                ! time counter
  REAL(KIND=8), DIMENSION(1,1)              :: dl_dum              ! working pseudo array for nc output

  TYPE(variable), DIMENSION(:),    ALLOCATABLE :: stypvar          ! structure for attributes

  CHARACTER(LEN=256)                        :: cf_tfil             ! input gridT file
  CHARACTER(LEN=256)                        :: cf_mfil             ! input MXL file
  CHARACTER(LEN=256)                        :: cf_out='saltc.nc'   ! netcdf output file
  CHARACTER(LEN=256)                        :: cldum               ! dummy character variable
  CHARACTER(LEN=256)                        :: cv_msk='tmask'      ! variable for masking

  LOGICAL                                   :: lfull=.FALSE.       ! flag for full step computation
  LOGICAL                                   :: lchk                ! flag for missing files
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage :  cdfsaltc  -f T-file [-mxloption option] [-mxlf MXL-file] ...'
     PRINT *,'     [-zoom imin imax jmin jmax kmin kmax] [-full] [-o OUT-file]'
     PRINT *,'     [-M MSK-file VAR-mask ] [-vvl ]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Compute the salt content in the specified 3D area (kg)'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f T-file : name of the input file with absolute salinity (and MLD if needed).'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-zoom imin imax jmin jmax kmin kmax] : limit of a sub domain where'
     PRINT *,'              the salt content will be calculated.'
     PRINT *,'                - if imin = 0 then ALL i are taken'
     PRINT *,'                - if jmin = 0 then ALL j are taken'
     PRINT *,'                - if kmin = 0 then ALL k are taken'
     PRINT *,'       [-full ] : assume full step model output instead of default'
     PRINT *,'               partial steps.'
     PRINT *,'       [-mxloption option]: option= 1 : compute only in the mixed layer,'
     PRINT *,'                            option=-1 : exclude mixed layer in the computation'
     PRINT *,'                            option= 0 : [Default], do not take care of mxl.'
     PRINT *,'       [-mxlf MXL-file ]: give name where MLD is found, if not in T-file.'
     PRINT *,'                  Only usefull if mxl option is not 0.'
     PRINT *,'       [-o OUT-file ] : specify netcdf output filename instead of ',TRIM(cf_out)
     PRINT *,'       [-M MSK-file VAR-mask] : Allow the use of a non standard mask file '
     PRINT *,'              with VAR-mask, instead of ',TRIM(cn_fmsk),' and ',TRIM(cv_msk) 
     PRINT *,'              This option is a usefull alternative to -zoom option, when the '
     PRINT *,'              area of interest is not ''box-like'' '
     PRINT *,'       [ -vvl ] : use time-varying  e3t for integration'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       Files ',TRIM(cn_fhgr),', ',TRIM(cn_fzgr),' and ',TRIM(cn_fmsk) 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : saltc.nc unless -o option is used.'
     PRINT *,'              variables: saltc3d (kg)'
     PRINT *,'                       : saltc2d(dep) (kg) '
     PRINT *,'                       : saltc3dpervol (kg/m3) '
     PRINT *,'       Standard output'
     PRINT *,'       '
     PRINT *,'      SEE ALSO: '
     PRINT *,'          cdfpolymask '
     STOP 
  ENDIF

  cf_mfil='none'
  ijarg = 1 
  DO WHILE ( ijarg <= narg ) 
     CALL getarg ( ijarg, cldum) ; ijarg = ijarg + 1
     SELECT CASE ( cldum )
     CASE ( '-f'       ) ; CALL getarg ( ijarg, cf_tfil) ; ijarg = ijarg + 1
     CASE ( '-full'    ) ; lfull = .true.
     CASE ( '-vvl'     ) ; lg_vvl = .true.
     CASE ('-mxloption') ; CALL getarg ( ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) mxloption
     CASE ( '-mxlf'    ) ; CALL getarg ( ijarg, cf_mfil) ; ijarg = ijarg + 1
     CASE ( '-o   '    ) ; CALL getarg ( ijarg, cf_out)  ; ijarg = ijarg + 1 
     CASE ( '-zoom'    ) ; CALL getarg ( ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) iimin
        ;                  CALL getarg ( ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) iimax
        ;                  CALL getarg ( ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) ijmin
        ;                  CALL getarg ( ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) ijmax
        ;                  CALL getarg ( ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) ikmin
        ;                  CALL getarg ( ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) ikmax
     CASE ( '-M'       ) ; CALL getarg ( ijarg, cn_fmsk) ; ijarg = ijarg + 1 
        ;                    CALL getarg ( ijarg, cv_msk ) ; ijarg = ijarg + 1
     CASE DEFAULT        ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  END DO

  IF ( cf_mfil == 'none' ) cf_mfil = cf_tfil

  lchk = chkfile(cn_fhgr)
  lchk = chkfile(cn_fzgr) .OR. lchk
  lchk = chkfile(cn_fmsk) .OR. lchk
  lchk = chkfile(cf_tfil) .OR. lchk
  lchk = chkfile(cf_mfil) .OR. lchk
  IF ( lchk ) STOP 99 ! missing files

  IF ( lg_vvl ) THEN 
     cn_fe3t = cf_tfil
     cn_ve3t = cn_ve3tvvl
  ENDIF

  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)
  npt    = getdim (cf_tfil,cn_t)

  npkk=npk

  IF (iimin /= 0 ) THEN ; npiglo = iimax - iimin + 1;  ELSE ; iimin=1             ; ENDIF
  IF (ijmin /= 0 ) THEN ; npjglo = ijmax - ijmin + 1;  ELSE ; ijmin=1             ; ENDIF
  IF (ikmin /= 0 ) THEN ; npkk   = ikmax - ikmin + 1;  ELSE ; ikmin=1 ; ikmax=npk ; ENDIF

  nvpk   = getvdim(cf_tfil,cn_vosaline)
  IF (nvpk == 2 ) nvpk = 1
  IF (nvpk == 3 ) nvpk = npkk

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt
  PRINT *, 'nvpk   = ', nvpk

  ! Allocate arrays
  ALLOCATE ( tmask(npiglo,npjglo))
  ALLOCATE ( salin(npiglo,npjglo))
  ALLOCATE ( e1t  (npiglo,npjglo), e2t(npiglo,npjglo), e3t(npiglo,npjglo))

  IF (mxloption /= 0) THEN
     ALLOCATE ( rmxldep(npiglo,npjglo))
  ENDIF

  ALLOCATE ( gdepw(npk), dtim(npt))
  IF ( lfull ) ALLOCATE ( e31d(npk))

  e1t(:,:) = getvar(cn_fhgr, cn_ve1t, 1, npiglo, npjglo, kimin=iimin, kjmin=ijmin)
  e2t(:,:) = getvar(cn_fhgr, cn_ve2t, 1, npiglo, npjglo, kimin=iimin, kjmin=ijmin)
  gdepw(:) = getvare3(cn_fzgr, cn_gdepw,  npk)

  IF ( lfull ) e31d(:) = getvare3(cn_fzgr, cn_ve3t1d, npk)

  CALL CreateOutput

  DO jt=1,npt
     IF ( lg_vvl ) THEN ; it = jt
     ELSE ;               it = 1
     ENDIF
     dvol = 0.d0
     dsum = 0.d0
     PRINT * ,'TIME : ', dtim(jt)
     IF (mxloption /= 0) rmxldep(:,:) = getvar(cf_mfil, cn_somxl010, 1, npiglo, npjglo, kimin=iimin, kjmin=ijmin, ktime=jt)

     DO jk = 1,nvpk
        ik = jk + ikmin -1
        ! Get absolute salinity at ik
        salin( :,:)   = getvar(cf_tfil, cn_vosaline, ik, npiglo, npjglo, kimin=iimin, kjmin=ijmin, ktime=jt)
        tmask(:,:)   = getvar(cn_fmsk, cv_msk,      ik, npiglo, npjglo, kimin=iimin, kjmin=ijmin          )           

        ! get e3t at level ik ( ps...)
        IF ( lfull ) THEN
           e3t(:,:) = e31d(ik)
        ELSE
           e3t(:,:) = getvar(cn_fe3t, cn_ve3t, ik, npiglo, npjglo, kimin=iimin, kjmin=ijmin, ktime=it, ldiom=.NOT.lg_vvl)
        ENDIF

        SELECT CASE ( mxloption ) 
        CASE ( 1 ) 
           e3t(:,:) = MAX ( 0., MIN( e3t,rmxldep-gdepw(ik) ) )
        CASE ( -1 )
           e3t(:,:) = MIN ( e3t, MAX( 0.,gdepw(ik)+e3t(:,:)-rmxldep ) )
        END SELECT

        dsurf  = sum(e1t * e2t       * tmask)
        dvol2d = sum(e1t * e2t * e3t * tmask)
        dvol   = dvol + dvol2d

        dsum2d = sum(e1t * e2t * e3t * salin * tmask)
        dsum   = dsum + dsum2d

        IF (dvol2d /= 0 )THEN
           PRINT *, ' Salt Content  at level ',ik,'(',gdepw(ik),' m) ',pprho0*dsum2d, 'surface = ',dsurf/1.e6,' km^2'
        ELSE
           PRINT *, ' No points in the water at level ',ik,'(',gdepw(ik),' m) '
        ENDIF
        dl_dum(1,1) = pprho0*dsum2d
        ierr = putvar(ncout, id_varout(jp_sc2d), dl_dum(:,:),jk, 1, 1, ktime=jt )

     END DO

     PRINT * ,' Total Salt content        : ', pprho0*dsum ,' kg'
     PRINT * ,' Total Salt content/volume : ', pprho0*dsum/dvol ,' kg/m3 '
     dl_dum(1,1)=pprho0*dsum
     ierr = putvar(ncout, id_varout(jp_sc3d), dl_dum(:,:),1, 1, 1, ktime=jt )
     dl_dum(1,1)=dl_dum(1,1)/dvol
     ierr = putvar(ncout, id_varout(jp_scvol), dl_dum(:,:),1, 1, 1, ktime=jt )
  END DO
  ierr = closeout(ncout )
CONTAINS
  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create the netcdf outputfile 
    !!
    !!----------------------------------------------------------------------
    ! so far in cdfsaltc, only 4 variables willbe output.
    ! indeed 4 scalar but that will be considered as (x,y,t) ie (1,1,t)
    INTEGER(KIND=4) :: ivar=3, ierr
    REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zdumlon, zdumlat
    !!----------------------------------------------------------------------
    ALLOCATE(stypvar(ivar) )
    ALLOCATE(    ipk(ivar), id_varout(ivar) )
    ALLOCATE( zdumlon(1,1), zdumlat(1,1) )
    zdumlon(:,:) = 0.
    zdumlat(:,:) = 0.
    ipk(:)= 1 
    ! define new variables for output 

    stypvar%scale_factor      = 1.
    stypvar%add_offset        = 0.
    stypvar%savelog10         = 0.
    stypvar%conline_operation = 'N/A'
    stypvar%caxis             = 'T'
    stypvar%cprecision        = 'r8'

    stypvar(jp_sc3d)%cname          = 'saltc3d'
    stypvar(jp_sc3d)%cunits         = 'kg'
    stypvar(jp_sc3d)%clong_name     = 'Total Salt Content'
    stypvar(jp_sc3d)%cshort_name    = 'saltc3d'

    stypvar(jp_scvol)%cname          = 'saltc3dpervol'
    stypvar(jp_scvol)%cunits         = 'kg/m3'
    stypvar(jp_scvol)%clong_name     = 'Total Salt Content per unit volume'
    stypvar(jp_scvol)%cshort_name    = 'saltc3dpervol'

    ipk(jp_sc2d) = npkk
    stypvar(jp_sc2d)%cname          = 'saltc2d'
    stypvar(jp_sc2d)%cunits         = 'kg'
    stypvar(jp_sc2d)%clong_name     = 'Salt Content at each selected level'
    stypvar(jp_sc2d)%cshort_name    = 'saltc2d'

    ncout =  create     (cf_out, 'none',  1,      1, npkk,   cdep='depthw' )
    ierr  = createvar   (ncout,  stypvar, ivar, ipk, id_varout             )
    ierr  = putheadervar(ncout,  cf_tfil, 1,      1, npkk,                 &
          &                 pnavlon=zdumlon, pnavlat=zdumlat,              &
          &                 pdep=gdepw(ikmin:ikmax),                       &
          &                 cdep='depthw'                                  )

    dtim = getvar1d(cf_tfil, cn_vtimec, npt)
    ierr = putvar1d(ncout,  dtim,  npt, 'T')

    DEALLOCATE( zdumlon, zdumlat)

  END SUBROUTINE CreateOutput

END PROGRAM cdfsaltc
