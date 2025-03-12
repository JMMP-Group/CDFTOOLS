PROGRAM cdfvrmp
  !!======================================================================
  !!                     ***  PROGRAM  cdfvrmp  ***
  !!=====================================================================
  !!  ** Purpose : Remap variables from one vertical grid to another.
  !!
  !!  ** Method  : Two methods are available:
  !!
  !!               1) the Piecewise Parabolic Method (PPM) via the PPR 
  !!                  pakcage of Engwirda & Kelley 2016
  !!                  (https://doi.org/10.48550/arXiv.1606.08188), a 
  !!                  Fortran-90 library designed to compute high-order 
  !!                  piecewise polynomial reconstructions and conservative 
  !!                  integral re-mappings on structured grids.
  !!               2) simple linear interpolation
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

  INTEGER(KIND=4)                                     :: ierr, ipos              ! working integer
  INTEGER(KIND=4)                                     :: narg, iargc, ijarg      ! command line
  INTEGER(KIND=4)                                     :: ji, jj, jk, jt, jvar    ! dummy loop index
  INTEGER(KIND=4)                                     :: npiglo, npjglo, npkinp  ! size of the INPUT mesh
  INTEGER(KIND=4)                                     :: npiout, npjout, npkout  ! size of the TARGET mesh
  INTEGER(KIND=4)                                     :: npt                     ! time-records of the INPUT file
  INTEGER(KIND=4)                                     :: nk_inp, nk_out          ! number of local wet levels
  INTEGER(KIND=4)                                     :: iimin=0, iimax=0        ! domain i-limits for computation
  INTEGER(KIND=4)                                     :: ijmin=0, ijmax=0        ! domain j-limits for computation
  INTEGER(KIND=4)                                     :: idep, idep_max          ! possible depth index, maximum
  INTEGER(KIND=4)                                     :: ncout                   ! ncid of output file
  INTEGER(KIND=4)   , DIMENSION(1)                    :: ipk, id_varout          ! only one output variable
  INTEGER(KIND=4)   , DIMENSION(:,:)    , ALLOCATABLE :: navloninp, navlatinp
  INTEGER(KIND=4)   , DIMENSION(:,:)    , ALLOCATABLE :: mbathyinp, mbathyout, msk_out, msk_inp !tmaskutil
  INTEGER(KIND=4)   , DIMENSION(:,:)    , ALLOCATABLE :: ssrmask, mbk_inp, mbk_out, wndmask

  REAL(KIND=4)                                        :: zspval                  ! missing value
  REAL(KIND=8)      , PARAMETER                       :: eps = 1.e-15            ! accuracy param
  REAL(KIND=8)                                        :: depinp, depout          ! ocean depth
  REAL(KIND=8)                                        :: depdiff                 ! to check accuracy
  REAL(KIND=8)      , DIMENSION(:)      , ALLOCATABLE :: dtim                    ! time counter
  REAL(KIND=8)      , DIMENSION(:)      , ALLOCATABLE :: ddep                    ! depth variable
  REAL(KIND=8)      , DIMENSION(:,:)    , ALLOCATABLE :: var_inp
  REAL(KIND=8)      , DIMENSION(:,:)    , ALLOCATABLE :: e3p_inp, e3p_out
  REAL(KIND=8)      , DIMENSION(:,:)    , ALLOCATABLE :: gdepp_inp, gdepp_out
  REAL(KIND=8)      , DIMENSION(:,:,:)  , ALLOCATABLE :: var_out
  !REAL(KIND=4)                                        :: depinp, depout          ! ocean depth
  !REAL(KIND=4)                                        :: depdiff                 ! to check accuracy
  !REAL(KIND=4)      , DIMENSION(:)      , ALLOCATABLE :: dtim                    ! time counter
  !REAL(KIND=4)      , DIMENSION(:)      , ALLOCATABLE :: ddep                    ! depth variable
  !REAL(KIND=4)      , DIMENSION(:,:)    , ALLOCATABLE :: var_inp
  !REAL(KIND=4)      , DIMENSION(:,:)    , ALLOCATABLE :: e3p_inp, e3p_out
  !REAL(KIND=4)      , DIMENSION(:,:)    , ALLOCATABLE :: gdepp_inp, gdepp_out
  !REAL(KIND=4)      , DIMENSION(:,:,:)  , ALLOCATABLE :: var_out


  CHARACTER(LEN=256)                                  :: cf_msh             ! input TARGET_MESH-file
  CHARACTER(LEN=256)                                  :: cf_inp             ! input T U V W files
  CHARACTER(LEN=256)                                  :: cldum              
  CHARACTER(LEN=256)                                  :: scheme             ! Type of remapping scheme
  CHARACTER(LEN=256)                                  :: cv_inp             ! input variable name
  CHARACTER(LEN=256)                                  :: cvtype             ! type of C-grid point to work with
  CHARACTER(LEN=256)                                  :: cv_e3              ! e3 scale factor name
  CHARACTER(LEN=256)                                  :: cv_msk             ! mask name
  CHARACTER(LEN=256)                                  :: cf_out=''          ! output file name
  CHARACTER(LEN=256)                                  :: clunits            ! attribute of output file : units
  CHARACTER(LEN=256)                                  :: cllong_name        ! attribute of output file : long name
  CHARACTER(LEN=256)                                  :: clshort_name       ! attribute of output file : short name
  CHARACTER(LEN=256)                                  :: cglobal            ! attribute of output file : global
  CHARACTER(LEN=256)                                  :: cdep, cv_dep       ! deptht name for dim and var
  CHARACTER(LEN=256), DIMENSION(:)      , ALLOCATABLE :: clv_dep            ! array of possible depth name 
                                                                            ! (or 3rd dimension)

  TYPE (variable)   , DIMENSION(:)      , ALLOCATABLE :: stypvar            ! Type variable is defined in cdfio.
  LOGICAL                                             :: lchk     = .FALSE. ! flag for missing files
  LOGICAL                                             :: lerror   = .FALSE. ! flag for missing arguments
  LOGICAL                                             :: lspwnd   = .FALSE. ! flag for spatial window
  LOGICAL                                             :: lverbose = .FALSE. ! flag for verbosity
  ! -------------------------------------------------------------------------------------------------------------
  !
  ! -----------------------
  ! 1. : Initialization
  ! -----------------------
  CALL ReadCdfNames()
  
  ! check argument number and show usage if necessary
  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage :  cdfvrmp -m TRG_MSH-file -f IN-file -v IN-var -p C-point ...'
     PRINT *,'       ... [-w imin imax jmin jmax kmin kmax] [-o OUT-file] [-verbose]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'      This tool uses conservative remapping operators to remap all the '
     PRINT *,'      variables included in the IN-file onto the TARGET vertical grid defined'
     PRINT *,'      in the TRG_MSH-file. '
     PRINT *,'      The INPUT  model geometry is defined in the files '
     PRINT *,'      ',   TRIM(cn_fhgr),' and ',TRIM(cn_fzgr),' in the current directory.'
     PRINT *,'      '
     PRINT *,'      This tool works only for the case where the INPUT model geometry and'
     PRINT *,'      the TARGET model geometry share the same horizontal grid. '
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -m TRG_MSH-file :  target mesh_mask.nc file.'
     PRINT *,'       -f IN-file      :  input netCDF file that needs to be vertically remapped'
     PRINT *,'                          onto the TARGET model geometry'
     PRINT *,'       -v IN-var       :  IN-var = name of netCDF variable to work with; '
     PRINT *,'       -p C-point      :  one of T|U|V|F|W indicating the position of IN-var on the'
     PRINT *,'                          C-grid.'
     PRINT *,'       -s REMAP-SCHEME :  one of ppm|linear indicating the scheme used for remapping.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-w imin imax jmin jmax] : spatial window where remapping is conducted'
     PRINT *,'       [-o OUT-file ]           : output filename instead of ''vrmp_<IN-var>.nc'' '
     PRINT *,'       [-verbose]               : increase verbosity  '
     PRINT *,'     '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ', TRIM(cn_fhgr),' and ',TRIM(cn_fzgr),' in the current directory '
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : cdfrmp.nc (default).'
     PRINT *,'         variables : all the variables included in the IN-file.'
     STOP
  ENDIF

  ! Mandatory arguments are set to none by default for further check
  cf_msh = 'none' ; cf_inp='none' ; cv_inp='none'; cvtype='none'
  ijarg = 1
  ! Parse command line
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1
     SELECT CASE  ( cldum   )
     CASE ( '-m'       ) ; CALL getarg(ijarg, cf_msh ) ; ijarg = ijarg + 1
     CASE ( '-f'       ) ; CALL getarg(ijarg, cf_inp ) ; ijarg = ijarg + 1
     CASE ( '-v'       ) ; CALL getarg(ijarg, cv_inp ) ; ijarg = ijarg + 1
     CASE ( '-p'       ) ; CALL getarg(ijarg, cvtype ) ; ijarg = ijarg + 1
     CASE ( '-s'       ) ; CALL getarg(ijarg, scheme ) ; ijarg = ijarg + 1
     ! options
     CASE ( '-w'       ) ; CALL getarg(ijarg, cldum  ) ; ijarg = ijarg + 1 ;  READ(cldum,*) iimin
                           CALL getarg(ijarg, cldum  ) ; ijarg = ijarg + 1 ;  READ(cldum,*) iimax
                           CALL getarg(ijarg, cldum  ) ; ijarg = ijarg + 1 ;  READ(cldum,*) ijmin
                           CALL getarg(ijarg, cldum  ) ; ijarg = ijarg + 1 ;  READ(cldum,*) ijmax
     CASE ( '-o'       ) ; CALL getarg(ijarg, cf_out ) ; ijarg = ijarg + 1 
     CASE ( '-verbose' ) ; lverbose=.TRUE.
     CASE DEFAULT        ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  IF ( cf_inp == 'none' )  THEN
     PRINT *,' You must specify an input file with -f option'
     lerror= lerror .OR. .TRUE.
  ENDIF
  IF ( cv_inp == 'none' )  THEN
     PRINT *,' You must specify an input variable with -v option'
     lerror= lerror .OR. .TRUE.
  ENDIF
  IF ( cvtype == 'none' )  THEN
     PRINT *,' You must specify a point type with -p option'
     lerror= lerror .OR. .TRUE.
  ENDIF
  IF (lerror ) STOP 99

  ! Checking all the needed file exists
  lchk = chkfile(cn_fhgr)
  lchk = chkfile(cn_fzgr) .OR. lchk
  lchk = chkfile(cn_fmsk) .OR. lchk
  lchk = chkfile(cf_msh ) .OR. lchk
  lchk = chkfile(cf_inp ) .OR. lchk
  IF ( lchk ) STOP 99 ! missing file

  ! Checking if spatial window is requested
  IF ( (iimin /= 0).AND.(iimax /= 0).AND.(iimin /= 0).AND.(iimax /= 0) ) lspwnd = .TRUE.

  ! Checking remapping scheme
  IF ( TRIM(scheme) == 'ppm') THEN
     PRINT*, "Using Piecewise Parabolic Method (PPM)."
  ELSEIF ( TRIM(scheme) == 'linear') THEN
     PRINT*, "Using linear interpolation (LINEAR)."
  ELSE
     PRINT*, "ERROR: we couldn't recognise the remapping method you specified."
     STOP 99
  ENDIF

  ! Get dimensions of input geometry and file
  npiglo = getdim (cf_inp, cn_x)
  npjglo = getdim (cf_inp, cn_y)
  npkinp = getdim (cf_inp, cn_z)
  npt    = getdim (cf_inp, cn_t)
  
  PRINT *, ' INPUT FILE:'
  PRINT *, ' ---------- '
  PRINT *, '   NPIGLO = ', npiglo
  PRINT *, '   NPJGLO = ', npjglo
  PRINT *, '   NPKinp = ', npkinp
  PRINT *, '   NPT    = ', npt

  ! Get dimensions of target geometry
  npiout = getdim (cf_msh, cn_x)
  npjout = getdim (cf_msh, cn_y)

  ! looking for npk among various possible name
  idep_max=4
  ALLOCATE ( clv_dep(idep_max) )
  clv_dep(:) = (/cn_z,'z','nav_lev','levels'/)
  idep=1  ; ierr=1000
  DO WHILE ( ierr /= 0 .AND. idep <= idep_max )
     npkout  = getdim (cf_msh, clv_dep(idep), cdtrue=cv_dep, kstatus=ierr)
     idep = idep + 1
  ENDDO

  IF ( ierr /= 0 ) THEN  ! none of the dim name was found
     PRINT *,' ERROR: we could not find the npk dimension of the target mesh!'
     STOP 99
  ENDIF

  PRINT *, ' TARGET MESH:'
  PRINT *, ' ---------- '
  PRINT *, '   NPIout = ', npiout
  PRINT *, '   NPJout = ', npjout
  PRINT *, '   NPKout = ', npkout

  ! Check that lateral dimensions of input and 
  ! target model geometry agree
  IF ( (  npiglo /= npiout ) .OR. ( npjglo /= npjout ) ) THEN
     PRINT *,'  ERROR: The input and target meshes MUST have the same lateral dimensions!'
     STOP 99
  ENDIF

  ! Allocate working arrays
  ALLOCATE ( stypvar  ( 1              ) )
  ALLOCATE ( e3p_inp  ( npiglo, npkinp ) )
  ALLOCATE ( e3p_out  ( npiglo, npkout ) )
  ALLOCATE ( msk_inp  ( npiglo, npkinp ) )
  ALLOCATE ( msk_out  ( npiglo, npkout ) )
  ALLOCATE ( var_inp  ( npiglo, npkinp ) )
  ALLOCATE ( mbk_inp  ( npiglo, npjglo ) )
  ALLOCATE ( mbk_out  ( npiglo, npjglo ) )
  ALLOCATE ( ssrmask  ( npiglo, npjglo ) )
  ALLOCATE ( wndmask  ( npiglo, npjglo ) )
  ALLOCATE ( navloninp( npiglo, npjglo ) )
  ALLOCATE ( navlatinp( npiglo, npjglo ) )
  ALLOCATE ( mbathyinp( npiglo, npjglo ) )
  ALLOCATE ( mbathyout( npiglo, npjglo ) )
  ALLOCATE ( var_out  ( npiglo, npjglo, npkout ) )

  IF ( scheme == 'linear' ) THEN
     ALLOCATE ( gdepp_inp  ( npiglo, npkinp ) )
     ALLOCATE ( gdepp_out  ( npiglo, npkout ) )
  ENDIF

  ! Workout the name and values of some mesh related variables
  mbathyinp(:,:) = getvar(cn_fzgr, cn_mbathy, 1, npiglo, npjglo)
  navloninp(:,:) = getvar(cf_inp , cn_vlon2d, 1, npiglo, npjglo)
  navlatinp(:,:) = getvar(cf_inp , cn_vlat2d, 1, npiglo, npjglo)

  SELECT CASE (TRIM(cvtype))
  CASE ( 'T' )
     IF ( scheme == 'ppm' ) THEN
        cv_e3 = cn_ve3t
     ELSE
        cv_e3 = cn_ve3w
     ENDIF                   
     cv_dep       = cn_gdept
     cdep         = cn_vdeptht
     cv_msk       = cn_tmask
     mbk_inp(:,:) = mbathyinp(:,:)
     mbk_out(:,:) = getvar(cf_msh , "mbkt", 1, npiglo, npjglo)
  CASE ( 'U' )
     IF ( scheme == 'ppm' ) THEN
        cv_e3 = cn_ve3u
     ELSE
        cv_e3 = 'e3uw_0'            
     ENDIF
     cv_dep       = cn_gdept
     cdep         = cn_vdepthu
     cv_msk       = cn_umask
     mbk_out(:,:) = getvar(cf_msh , "mbku", 1, npiglo, npjglo)
     DO jj=1,npjglo
        DO ji=1,npiglo
           mbk_inp(ji,jj) = MIN( mbathyinp(ji+1,jj), mbathyinp(ji,jj) )
        END DO
     END DO 
  CASE ( 'V' )
     IF ( scheme == 'ppm' ) THEN
        cv_e3 = cn_ve3v
     ELSE
        cv_e3 = 'e3vw_0'
     ENDIF
     cv_dep       = cn_gdept
     cdep         = cn_vdepthv
     cv_msk       = cn_vmask
     mbk_out(:,:) = getvar(cf_msh , "mbkv", 1, npiglo, npjglo)
     DO jj=1,npjglo
        DO ji=1,npiglo
           mbk_inp(ji,jj) = MIN( mbathyinp(ji,jj+1), mbathyinp(ji,jj) )
        END DO
     END DO
  !CASE ( 'W' )
  !   cv_e3        = cn_ve3w
  !   cv_dep       = cn_gdepw
  !   cdep         = cn_vdepthw
  !   cv_msk       = cn_vmask
  !   mbk_inp(:,:) = mbathyinp(:,:)
  !   mbk_out(:,:) = getvar(cf_msh , "mbkt", 1, npiglo, npjglo)
  CASE DEFAULT
     PRINT *, 'this type of variable is not known :', TRIM(cvtype)
     STOP 99
  END SELECT

  CALL CreateOutputFile

  ! Main loop
  DO jt=1,npt
     ! Initisalise output variable
     var_out(:,:,:) = 0.0
     IF ( jt == 1 ) THEN
        ssrmask(:,:) = getvar(cn_fmsk, cv_msk , 1, npiglo, npjglo)
        IF ( lspwnd ) THEN
           wndmask(:,:) = 0
           wndmask(iimin:iimax,ijmin:ijmax) = 1 
        ELSE
           wndmask(:,:) = 1
        END IF  
     ENDIF
     DO jj=1,npjglo
        ! get vertical scale factors
        IF ( jt == 1 ) THEN
           msk_inp(:,:) = getvarxz(cn_fzgr, cv_msk, jj, npiglo, npkinp)
           e3p_inp(:,:) = getvarxz_dp(cn_fzgr, cv_e3 , jj, npiglo, npkinp)
           e3p_out(:,:) = getvarxz_dp(cf_msh , cv_e3 , jj, npiglo, npkout)
           IF ( TRIM(scheme) == 'linear') THEN
              gdepp_inp(:,1) = 0.5 * e3p_inp(:,1)
              gdepp_out(:,1) = 0.5 * e3p_out(:,1)
              DO jk = 2, npkinp
                 gdepp_inp(:,jk) = gdepp_inp(:,jk-1) + e3p_inp(:,jk)
              END DO
              DO jk = 2, npkout
                 gdepp_out(:,jk) = gdepp_out(:,jk-1) + e3p_out(:,jk)
              END DO
           ENDIF
        END IF
        ! get variable
        var_inp(:,:) = getvarxz_dp(cf_inp, cv_inp, jj, npiglo, npkinp, ktime=jt) * msk_inp(:,:)
        DO ji=1,npiglo
           IF ( ssrmask(ji,jj)*wndmask(ji,jj) == 1 ) THEN
              nk_inp = mbk_inp(ji,jj)
              nk_out = mbk_out(ji,jj)
              IF ( lverbose .AND. scheme == "ppm") THEN 
                 depinp  = SUM(e3p_inp(ji,1:nk_inp))
                 depout  = SUM(e3p_out(ji,1:nk_out))
                 depdiff = ABS(depinp-depout)
                 IF ( depdiff > eps ) THEN
                    PRINT*, 'i=', ji, ', j=', jj, ':'
                    PRINT*, '     depdiff=', depdiff
                    PRINT*, '     INP depth: sum(e3)=', depinp, ', nk_inp=', nk_inp
                    PRINT*, '     OUT depth: sum(e3)=', depout, ', nk_out=', nk_out
                 ENDIF
              ENDIF
              IF ( TRIM(scheme) == 'ppm') THEN
                 CALL reconstructandremap( var_inp(ji,   1:nk_inp), e3p_inp(ji,1:nk_inp), & 
                              &            var_out(ji,jj,1:nk_out), e3p_out(ji,1:nk_out), &
                              &            nk_inp                 , nk_out              , 1 )
                 !IF ( (ji==66) .AND. (jj==171) ) THEN
                 !   PRINT*, 'i=', ji, ', j=', jj, ':'
                 !   DO jk=1,nk_inp
                 !      PRINT*, 'jk=', jk, ' e3p_inp=', e3p_inp(ji,jk), ' var_inp=', var_inp(ji,jk)
                 !   ENDDO
                 !   DO jk=1,nk_out
                 !      PRINT*, 'jk=', jk, ' e3p_out=', e3p_out(ji,jk), ' var_out=', var_out(ji,jj,jk)
                 !   ENDDO
                 !ENDIF
              ELSE
                 CALL remap_linear(var_inp(ji,   1:nk_inp), gdepp_inp(ji,1:nk_inp), &
                              &    var_out(ji,jj,1:nk_out), gdepp_out(ji,1:nk_out), & 
                              &    nk_inp                 , nk_out                , 1)
              ENDIF   
           ELSE
              var_out(ji,jj,:) = var_inp(ji,:)
           END IF
        END DO
     END DO
     DO jk=1,npkout
        ierr = putvar(ncout, id_varout(1), var_out(:,:,jk), jk, npiglo, npjglo, ktime=jt)
     END DO         
  END DO 

  ierr = closeout(ncout)

CONTAINS

  SUBROUTINE CreateOutputFile
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutputFile  ***
    !!
    !! ** Purpose :  Perform output file creation with all the variables 
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
    ipk(:)                       = npkout
    ierr = getvaratt(cf_inp, cv_inp, clunits, zspval, cllong_name, clshort_name)
    !
    stypvar(1)%rmissing_value    = 1.e+20
    stypvar(1)%scale_factor      = 1.
    stypvar(1)%add_offset        = 0.
    stypvar(1)%savelog10         = 0.
    stypvar(1)%conline_operation = 'N/A'

    stypvar(1)%cunits            = TRIM(clunits)
    stypvar(1)%valid_min         = -50.
    stypvar(1)%valid_max         =  50.
    stypvar(1)%cname             = TRIM(cv_inp)
    stypvar(1)%clong_name        = TRIM(cllong_name)
    stypvar(1)%cshort_name       = TRIM(clshort_name)

    ! create output fileset
    IF ( cf_out == '' ) cf_out = 'vrmp_'//TRIM(cv_inp)//'_'//TRIM(scheme)//'.nc'  ! set default name if not 
                                                               ! specified on command line
    !ncout = create      (cf_out, cf_inp, npiglo, npjglo, npkout, cdep=cdep, cdepvar=cv_dep)
    ncout = create      (cf_out, cf_inp, npiglo, npjglo, npkout, cdep=cdep, cdepvar=cdep)
    ierr  = createvar   (ncout ,      stypvar , 1 , ipk,   id_varout)
    !ierr  = putheadervar(ncout ,      cf_inp  ,  ikx, iky, npk, pnavlon=rdumlon, pnavlat=rdumlat, pdep=gdep, cdep=cv_dep)
    dtim  = getvar1d(cf_inp, cn_vtimec, npt )
    ddep  = getvar1d(cf_msh, cv_dep, npkout )
    ierr  = putvar1d(ncout, dtim, npt, 'T')
    ierr  = putvar1d(ncout, ddep, npkout, 'D')

  END SUBROUTINE CreateOutputFile

END PROGRAM cdfvrmp

  
