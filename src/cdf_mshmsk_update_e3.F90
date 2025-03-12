PROGRAM cdf_mshmsk_update_e3
  !!======================================================================
  !!                     ***  PROGRAM cdf_mshmsk_update_e3  ***
  !!======================================================================
  !!  ** Purpose : Udate the vertical scale factors of a mesh_mask.nc file 
  !!               to match the volume of another mesh_mask.nc file.
  !!               Needed when conducting conservative vertical remapping
  !!               with cdfvrmp.f90.
  !!
  !!  ** Method  : In each grid point, update the bottom e3t to match the 
  !!               ocean volume of another model domain passed in input.
  !!
  !! History : 4.0  : 02/2025  : D. Bruciaferri    : Original code  
  !!=====================================================================
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------
  USE netcdf
  USE cdfio
  USE modcdfnames

  IMPLICIT NONE

  INTEGER(KIND=4)                                     :: ierr                    ! working integer
  INTEGER(KIND=4)                                     :: narg, iargc, ijarg      ! command line
  INTEGER(KIND=4)                                     :: ji, jj, jk, jt          ! dummy loop index
  INTEGER(KIND=4)                                     :: npiglo, npjglo, npkinp  ! size of the INPUT mesh
  INTEGER(KIND=4)                                     :: npiout, npjout, npkout  ! size of the TARGET mesh
  INTEGER(KIND=4)                                     :: npt                     ! time-records of the INPUT file
  INTEGER(KIND=4)                                     :: idep, idep_max          ! possible depth index, maximum
  INTEGER(KIND=4)                                     :: nczgr                   ! ncid of output file
  INTEGER                                             :: id_dept1d, id_depw1d, id_gdept, id_gdepw
  INTEGER                                             :: id_e3t, id_e3u , id_e3v
  INTEGER                                             :: id_e3w, id_e3uw, id_e3vw
  INTEGER                                             :: id_mbkt, id_mbku, id_mbkv, id_bathy
  INTEGER                                             :: id_navlat, id_navlon, id_navlev, id_time
  INTEGER                                             :: idx, idy, idz, idt
  INTEGER                                             :: id_tmsk, id_umsk, id_vmsk, id_fmsk
 
  INTEGER(KIND=4)   , DIMENSION(:,:)    , ALLOCATABLE :: mbkt_trg, mbk_trg, mskup
  INTEGER(KIND=4)   , DIMENSION(:,:)    , ALLOCATABLE :: mbkt, mbku, mbkv
  INTEGER(KIND=4)   , DIMENSION(:,:)    , ALLOCATABLE :: tmask, umask, vmask, fmask
  INTEGER           ,                       PARAMETER :: dp = SELECTED_REAL_KIND(15,307) ! double precision (real 8)

  REAL(dp)          , DIMENSION(:)      , ALLOCATABLE :: dtim                    ! time counter
  REAL(dp)          , DIMENSION(:)      , ALLOCATABLE :: gdept_1d, gdepw_1d      ! depth variable
  REAL(dp)          , DIMENSION(:,:)    , ALLOCATABLE :: hdep_trg, e3_trg 
  REAL(dp)          , DIMENSION(:,:)    , ALLOCATABLE :: navlon, navlat
  REAL(dp)          , DIMENSION(:,:)    , ALLOCATABLE :: e3w
  REAL(dp)          , DIMENSION(:,:)    , ALLOCATABLE :: gdepw !, gdept
  REAL(dp)          , DIMENSION(:,:)    , ALLOCATABLE :: zwtmp, zttmp
  REAL(dp)          , DIMENSION(:,:,:)  , ALLOCATABLE :: e3, gdept

  CHARACTER(LEN=256)                                  :: cf_trg             ! input TARGET_MESH-file
  CHARACTER(LEN=256)                                  :: cf_inp             ! input T U V W files
  CHARACTER(LEN=256)                                  :: cv_lev             ! type of vertical coordinates
  CHARACTER(LEN=256)                                  :: cldum              ! can handle a long list of section files
  CHARACTER(LEN=256)                                  :: cf_zgr='zgr.nc'    ! output file name
  CHARACTER(LEN=256)                                  :: cdep, cv_dep       ! deptht name for dim and var
  CHARACTER(LEN=256), DIMENSION(:)      , ALLOCATABLE :: clv_dep            ! array of possible depth name 
                                                                            ! (or 3rd dimension)

  TYPE (variable)   , DIMENSION(:)      , ALLOCATABLE :: stypvar            ! Type variable is defined in cdfio.
  LOGICAL                                             :: lchk     = .FALSE. ! flag for missing files
  LOGICAL                                             :: lerror   = .FALSE. ! flag for missing arguments
  ! -------------------------------------------------------------------------------------------------------------
  !
  ! -----------------------
  ! 1. : Initialization
  ! -----------------------
  CALL ReadCdfNames()
  
  ! check argument number and show usage if necessary
  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage :  cdf_mshmsk_update_e3 -i INP-file -lev COORD -t TRG-file [-o OUT-file] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       This tool updates the vertical scale factors of an INPUT mesh_mask.nc'
     PRINT *,'       file (i.e., INP-file) to match the volume of a TARGET mesh_mask.nc file' 
     PRINT *,'       (i.e., TRG-file). Needed when conducting conservative vertical remapping'
     PRINT *,'       with cdfvrmp.f90.'
     PRINT *,'       This tool works only for the case where the INPUT model geometry and'
     PRINT *,'       the TARGET model geometry share the same horizontal grid. '
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -i INP-file     :  input mesh_mask.nc file that will be updated.'
     PRINT *,'       -lev COORD      :  one of z|zps|gvc '
     PRINT *,'       -t TRG-file     :  target mesh_mask.nc file.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-o OUT-file ]  :  output filename instead of ''vrmp_<IN-var>.nc'' '
     PRINT *,'     '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       None.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : updated_mesh_mask.nc (default).'
     STOP
  ENDIF

  ! Mandatory arguments are set to none by default for further check
  cf_inp='none' ; cv_lev='none'; cf_trg='none';
  ijarg = 1
  ! Parse command line
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1
     SELECT CASE  ( cldum   )
     CASE ( '-i'       ) ; CALL getarg(ijarg, cf_inp ) ; ijarg = ijarg + 1
     CASE ( '-lev'     ) ; CALL getarg(ijarg, cv_lev ) ; ijarg = ijarg + 1
     CASE ( '-t'       ) ; CALL getarg(ijarg, cf_trg ) ; ijarg = ijarg + 1
     ! options
     CASE ( '-o'       ) ; CALL getarg(ijarg, cf_zgr ) ; ijarg = ijarg + 1 
     CASE DEFAULT        ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  IF ( cf_inp == 'none' )  THEN
     PRINT *,' You must specify an input mesh_mask.nc file with the -i option'
     lerror= lerror .OR. .TRUE.
  ENDIF
  IF ( cv_lev == 'none' )  THEN
     PRINT *,' You must specify the type of vertical coordinates of the input mesh_mask.nc file with the -lev option'
     lerror= lerror .OR. .TRUE.
  ENDIF
  IF ( ( cv_lev /= 'z' ) .AND. ( cv_lev /= 'zps' ) .AND. ( cv_lev /= 'gvc' ) ) THEN
     PRINT *,' The type of vertical coordinates is unkown - please chosee between "z", "zps", or "gvc"!'
     lerror= lerror .OR. .TRUE.
  ENDIF

  IF ( cf_trg == 'none' )  THEN
     PRINT *,' You must specify a target mesh_mask.nc file with the -t option'
     lerror= lerror .OR. .TRUE.
  ENDIF
  IF (lerror ) STOP 99

  ! Checking all the needed file exists
  lchk = chkfile(cf_inp ) .OR. lchk
  lchk = chkfile(cf_trg ) .OR. lchk
  IF ( lchk ) STOP 99 ! missing file

  ! Get dimensions of input geometry and file
  npiglo = getdim (cf_inp, cn_x)
  npjglo = getdim (cf_inp, cn_y)

  ! looking for npk among various possible name
  idep_max=4
  ALLOCATE ( clv_dep(idep_max) )
  clv_dep(:) = (/cn_z,'z','nav_lev','levels'/)
  idep=1  ; ierr=1000
  DO WHILE ( ierr /= 0 .AND. idep <= idep_max )
     npkinp  = getdim (cf_inp, clv_dep(idep), cdtrue=cv_dep, kstatus=ierr)
     idep = idep + 1
  ENDDO

  IF ( ierr /= 0 ) THEN  ! none of the dim name was found
     PRINT *,' ERROR: we could not find the npk dimension of the input mesh!'
     STOP 99
  ENDIF

  PRINT *, ' INPUT MESH_MASK.NC:'
  PRINT *, ' ---------- '
  PRINT *, '   NPIGLO = ', npiglo
  PRINT *, '   NPJGLO = ', npjglo
  PRINT *, '   NPKinp = ', npkinp

  ! Get dimensions of target geometry
  npiout = getdim(cf_trg, cn_x)
  npjout = getdim(cf_trg, cn_y)

  ! looking for npk among various possible name
  idep=1  ; ierr=1000
  DO WHILE ( ierr /= 0 .AND. idep <= idep_max )
     npkout  = getdim (cf_trg, clv_dep(idep), cdtrue=cv_dep, kstatus=ierr)
     idep = idep + 1
  ENDDO

  IF ( ierr /= 0 ) THEN  ! none of the dim name was found
     PRINT *,' ERROR: we could not find the npk dimension of the target mesh!'
     STOP 99
  ENDIF

  PRINT *, ' TARGET MESH_MASK.NC:'
  PRINT *, ' ---------- '
  PRINT *, '   NPItrg = ', npiout
  PRINT *, '   NPJtrg = ', npjout
  PRINT *, '   NPKtrg = ', npkout

  ! Check that lateral dimensions of input and 
  ! target model geometry agree
  IF ( (  npiglo /= npiout ) .OR. ( npjglo /= npjout ) ) THEN
     PRINT *,'  ERROR: The input and target meshes MUST have the same lateral dimensions!'
     STOP 99
  ENDIF

  ! Allocate working arrays
  ALLOCATE ( gdept_1d   ( npkout ) )
  ALLOCATE ( gdepw_1d   ( npkout ) )
  ALLOCATE ( mbkt       ( npiglo, npjglo ) )
  ALLOCATE ( mbku       ( npiglo, npjglo ) )
  ALLOCATE ( mbkv       ( npiglo, npjglo ) )
  ALLOCATE ( mskup      ( npiglo, npjglo ) )
  ALLOCATE ( navlon     ( npiglo, npjglo ) )
  ALLOCATE ( navlat     ( npiglo, npjglo ) )
  ALLOCATE ( hdep_trg   ( npiglo, npjglo ) )
  ALLOCATE ( mbkt_trg   ( npiglo, npjglo ) )
  ALLOCATE ( mbk_trg    ( npiglo, npjglo ) )
  ALLOCATE ( e3_trg     ( npiglo, npjglo ) )
  ALLOCATE ( e3w        ( npiglo, npjglo) )
  ALLOCATE ( gdepw      ( npiglo, npjglo ) )
  ALLOCATE ( zttmp      ( npiglo, npjglo ) )
  ALLOCATE ( zwtmp      ( npiglo, npjglo ) ) 
  ALLOCATE ( tmask      ( npiglo, npjglo ) )
  ALLOCATE ( umask      ( npiglo, npjglo ) )
  ALLOCATE ( vmask      ( npiglo, npjglo ) )
  ALLOCATE ( fmask      ( npiglo, npjglo ) )
  ALLOCATE ( e3         ( npiglo, npjglo, npkinp ) )
  ALLOCATE ( gdept      ( npiglo, npjglo, npkinp ) )
 
  mbkt_trg(:,:) = getvar  (cf_trg, cn_mbathy, 1  , npiglo, npjglo)
  tmask(:,:)    = getvar  (cf_trg, cn_tmask , 1  , npiglo, npjglo)   
  umask(:,:)    = getvar  (cf_trg, cn_umask , 1  , npiglo, npjglo)
  vmask(:,:)    = getvar  (cf_trg, cn_vmask , 1  , npiglo, npjglo)
  gdept_1d(:)   = getvar1d(cf_inp, cn_gdept , npkinp             )
  gdepw_1d(:)   = getvar1d(cf_inp, cn_gdepw , npkinp             )
  navlon(:,:)   = getvar  (cf_inp, cn_vlon2d, 1  , npiglo, npjglo)
  navlat(:,:)   = getvar  (cf_inp, cn_vlat2d, 1  , npiglo, npjglo)

  CALL CreateMeshZgrFile

  ! Initisalise output variable
  mskup(:,:) = 1
  e3(:,:,:)  = 0.0d0
  gdepw(:,:) = 0.0d0
  gdept(:,:,:) = 0.0d0

  ! ----------------------------------------------------------------------------------
  ! T-GRID
  ! ----------------------------------------------------------------------------------

  ! Compute the bathymetry of the target mesh_mask @ T-grid
  hdep_trg(:,:) = 0.0d0
  DO jj = 1, npjglo
     e3_trg(:,:) = getvarxz_dp(cf_trg, cn_ve3t, jj, npiglo, npkout)
     DO ji = 1, npiglo
        hdep_trg(ji,jj) = SUM( e3_trg(ji, 1:mbkt_trg(ji,jj) ) ) * tmask(ji,jj)
     END DO
  END DO

  ! Initialising land points
  WHERE ( hdep_trg(:,:) == 0 ) mskup(:,:)  = 0

  ! Updating ocean e3t and related variables to match volume of target mesh_mask
  DO jk = 1, npkinp
     ! Store useful arrays
     zwtmp(:,:) = gdepw(:,:)    ! W-level @ jk
     IF ( jk > 1 ) THEN
        zttmp(:,:) = gdept(:,:,jk-1) ! T-level @ jk-1
     ELSE
        zttmp(:,:) = 0.0d0
     ENDIF
     ! Update e3t @ jk if needed
     e3(:,:,jk) = getvar(cf_inp, cn_ve3t, jk, npiglo, npjglo) 
     gdepw(:,:) = zwtmp(:,:) + e3(:,:,jk) ! W-level @ jk+1, first guess
     WHERE ( ( gdepw(:,:) >= hdep_trg(:,:) ) .AND. ( mskup(:,:) == 1 ) )
        e3(:,:,jk) = e3(:,:,jk) - (gdepw(:,:) - hdep_trg(:,:))
        mskup(:,:) = 0
     END WHERE
     ! Compute gdepw, gdept and e3w
     gdepw(:,:) = zwtmp(:,:) + e3(:,:,jk)               ! W-level @ jk+1, correct
     gdept(:,:,jk) = 0.5d0 * ( zwtmp(:,:) +  gdepw(:,:) ) ! T-level @ jk, as mean value
     e3w(:,:)   = gdept(:,:,jk) - zttmp(:,:)            ! e3w @ jk
     IF ( jk == 1 ) e3w(:,:) = 2.0d0 * e3w(:,:)           ! @ the surface
     ! e3w @ jk
     ierr = NF90_PUT_VAR( nczgr, id_e3w, e3w, start=(/1,1,jk,1/), count=(/npiglo,npjglo,1,1/) )
     IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
     ENDIF
  END DO

  ! e3t
  ierr = NF90_PUT_VAR( nczgr, id_e3t, e3, start=(/1,1,1,1/), count=(/npiglo,npjglo,npkinp,1/) )
  IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
  ENDIF
  
  ! bathy_meter
  ierr = NF90_PUT_VAR( nczgr, id_bathy, hdep_trg, start=(/1,1/), count=(/npiglo,npjglo/) )
  IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
  ENDIF
  
  ! mbkt
  WHERE( hdep_trg(:,:) == 0. ) ; mbkt(:,:) = 0        ! land
  ELSE WHERE                   ; mbkt(:,:) = npkinp-1 ! ocean 
  END WHERE

  DO jj = 1, npjglo
     DO ji = 1, npiglo
        DO jk = 1, npkinp-1
           IF ( hdep_trg(ji,jj) >= gdept(ji,jj,jk) ) mbkt(ji,jj) = MAX( 2, jk )
        END DO
     END DO
  END DO
  WHERE (hdep_trg(:,:)<=0) mbkt(:,:) = 0

  ierr = NF90_PUT_VAR( nczgr, id_mbkt, mbkt, start=(/1,1/), count=(/npiglo,npjglo/) )
  IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
  ENDIF

  ! ----------------------------------------------------------------------------------
  ! U-GRID
  ! ----------------------------------------------------------------------------------

  ! Compute the bathymetry of the target mesh_mask @ U-grid
  hdep_trg(:,:) = 0.0d0
  e3(:,:,:)     = 0.0d0
  mskup(:,:)    = 1
  mbk_trg(:,:)  = npkout-1

  ! Compute number of wet levels at U-points
  DO jj=1,npjglo
     DO ji=1,npiglo
        mbk_trg(ji,jj) = MIN( mbkt_trg(ji+1,jj), mbkt_trg(ji,jj) )
     END DO
  END DO

  ! Compute ocean depth at U-points 
  DO jj = 1, npjglo
     e3_trg(:,:) = getvarxz_dp(cf_trg, cn_ve3u, jj, npiglo, npkout)
     DO ji = 1, npiglo
        hdep_trg(ji,jj) = SUM( e3_trg(ji, 1:mbk_trg(ji,jj) ) ) * umask(ji,jj)
     END DO
  END DO

  ! Initialising land points
  WHERE ( hdep_trg(:,:) == 0 ) mskup(:,:)  = 0

  ! Updating ocean e3u to match volume of target mesh_mask
  gdepw(:,:) = 0.0d0
  gdept(:,:,:) = 0.0d0 ! Depth of T-levels @ U-points
  DO jk = 1, npkinp
     ! Store useful arrays
     zwtmp(:,:) = gdepw(:,:) ! W-level @ jk
     IF ( jk > 1 ) THEN
        zttmp(:,:) = gdept(:,:,jk-1) ! T-level @ jk-1
     ELSE
        zttmp(:,:) = 0.0d0
     ENDIF
     ! Update e3u @ jk if needed
     e3(:,:,jk) = getvar(cf_inp, cn_ve3u, jk, npiglo, npjglo)
     gdepw(:,:) = zwtmp(:,:) + e3(:,:,jk)
     WHERE ( ( gdepw(:,:) >= hdep_trg(:,:) ) .AND. ( mskup(:,:) == 1 ) )
        e3(:,:,jk)    = e3(:,:,jk) - (gdepw(:,:) - hdep_trg(:,:))
        mskup(:,:) = 0
     END WHERE
     ! Compute gdepw, gdept
     gdepw(:,:) = zwtmp(:,:) + e3(:,:,jk)               ! W-level @ jk+1, correct
     gdept(:,:,jk) = 0.5d0 * ( zwtmp(:,:) +  gdepw(:,:) ) ! T-level @ jk, as mean value
     e3w(:,:)   = gdept(:,:,jk) - zttmp(:,:)            ! e3uw @ jk
     ! e3uw @ jk
     ierr = NF90_PUT_VAR( nczgr, id_e3uw, e3w, start=(/1,1,jk,1/), count=(/npiglo,npjglo,1,1/) )
     IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
     ENDIF
  END DO

  ! Saving e3u
  ierr = NF90_PUT_VAR( nczgr, id_e3u, e3, start=(/1,1,1,1/), count=(/npiglo,npjglo,npkinp,1/) )
  IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
  ENDIF

  ! mbku
  WHERE( hdep_trg(:,:) == 0.0d0 ) ; mbku(:,:) = 0        ! land
  ELSE WHERE                   ; mbku(:,:) = npkinp-1 ! ocean 
  END WHERE

  DO jj = 1, npjglo
     DO ji = 1, npiglo
        DO jk = 1, npkinp-1
           IF ( hdep_trg(ji,jj) >= gdept(ji,jj,jk) ) mbku(ji,jj) = MAX( 2, jk )
        END DO
     END DO
  END DO
  WHERE (hdep_trg(:,:)<=0.0d0) mbku(:,:) = 0

  ierr = NF90_PUT_VAR( nczgr, id_mbku, mbku, start=(/1,1/), count=(/npiglo,npjglo/) )
  IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
  ENDIF

  ! ----------------------------------------------------------------------------------
  ! V-GRID
  ! ----------------------------------------------------------------------------------

  ! Compute the bathymetry of the target mesh_mask @ V-grid
  hdep_trg(:,:) = 0.0d0
  e3(:,:,:)     = 0.0d0
  mskup(:,:)    = 1
  mbk_trg(:,:)  = npkout-1

  ! Compute number of wet levels at V-points
  DO jj=1,npjglo
     DO ji=1,npiglo
        mbk_trg(ji,jj) = MIN( mbkt_trg(ji,jj+1), mbkt_trg(ji,jj) )
     END DO
  END DO

  ! Compute ocean depth at V-points 
  DO jj = 1, npjglo-1
     e3_trg(:,:) = getvarxz_dp(cf_trg, cn_ve3v, jj, npiglo, npkout)
     DO ji = 1, npiglo
        hdep_trg(ji,jj) = SUM( e3_trg(ji, 1:mbk_trg(ji,jj) ) ) * vmask(ji,jj)
     END DO
  END DO

  ! Initialising land points
  WHERE ( hdep_trg(:,:) == 0.0d0 ) mskup(:,:)  = 0

  ! Updating ocean e3v to match volume of target mesh_mask
  gdepw(:,:) = 0.0d0
  gdept(:,:,:) = 0.0d0 ! Depth of T-levels @ V-points
  DO jk = 1, npkinp
     ! Store useful arrays
     zwtmp(:,:) = gdepw(:,:) ! W-level @ jk
     IF ( jk > 1 ) THEN
        zttmp(:,:) = gdept(:,:,jk-1) ! T-level @ jk-1
     ELSE
        zttmp(:,:) = 0.0d0
     ENDIF
     ! Update e3v @ jk if needed
     e3(:,:,jk) = getvar(cf_inp, cn_ve3v, jk, npiglo, npjglo)
     gdepw(:,:) = zwtmp(:,:) + e3(:,:,jk)
     WHERE ( ( gdepw(:,:) >= hdep_trg(:,:) ) .AND. ( mskup(:,:) == 1 ) )
        e3(:,:,jk)    = e3(:,:,jk) - (gdepw(:,:) - hdep_trg(:,:))
        mskup(:,:)  = 0
     END WHERE
     ! Compute gdepw, gdept
     gdepw(:,:) = zwtmp(:,:) + e3(:,:,jk)               ! W-level @ jk+1, correct
     gdept(:,:,jk) = 0.5d0 * ( zwtmp(:,:) +  gdepw(:,:) ) ! T-level @ jk, as mean value
     e3w(:,:)   = gdept(:,:,jk) - zttmp(:,:)            ! e3uw @ jk
     ! e3vw @ jk
     ierr = NF90_PUT_VAR( nczgr, id_e3vw, e3w, start=(/1,1,jk,1/), count=(/npiglo,npjglo,1,1/) )
     IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
     ENDIF
  END DO

  ! Saving e3v
  ierr = NF90_PUT_VAR( nczgr, id_e3v, e3, start=(/1,1,1,1/), count=(/npiglo,npjglo,npkinp,1/) )
  IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
  ENDIF
 
  ! mbkv
  WHERE( hdep_trg(:,:) == 0.0d0 ) ; mbkv(:,:) = 0        ! land
  ELSE WHERE                   ; mbkv(:,:) = npkinp-1 ! ocean 
  END WHERE

  DO jj = 1, npjglo
     DO ji = 1, npiglo
        DO jk = 1, npkinp-1
           IF ( hdep_trg(ji,jj) >= gdept(ji,jj,jk) ) mbkv(ji,jj) = MAX( 2, jk )
        END DO
     END DO
  END DO
  WHERE (hdep_trg(:,:)<=0.0d0) mbkv(:,:) = 0
 
  ierr = NF90_PUT_VAR( nczgr, id_mbkv, mbkv, start=(/1,1/), count=(/npiglo,npjglo/) )
  IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
  ENDIF

  ! ----------------------------------------------------------------------------------
  ! MASKS
  ! ----------------------------------------------------------------------------------

  DO jk = 1, npkinp
     tmask(:,:) = 0.
     umask(:,:) = 0.
     vmask(:,:) = 0.
     WHERE ( mbkt(:,:) >= jk ) tmask(:,:) = 1.
     WHERE ( mbku(:,:) >= jk ) umask(:,:) = 1.
     WHERE ( mbkv(:,:) >= jk ) vmask(:,:) = 1.
     ! Saving variables @ jk
     ierr = NF90_PUT_VAR( nczgr, id_tmsk, tmask, start=(/1,1,jk,1/), count=(/npiglo,npjglo,1,1/) )
     IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
     ENDIF
     ierr = NF90_PUT_VAR( nczgr, id_umsk, umask, start=(/1,1,jk,1/), count=(/npiglo,npjglo,1,1/) )
     IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
     ENDIF
     ierr = NF90_PUT_VAR( nczgr, id_vmsk, vmask, start=(/1,1,jk,1/), count=(/npiglo,npjglo,1,1/) )
     IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
     ENDIF

  END DO
     
  ierr = NF90_CLOSE(nczgr)

CONTAINS

  SUBROUTINE CreateMeshZgrFile

    ierr= NF90_CREATE(cf_zgr, or(NF90_CLOBBER,NF90_NETCDF4), nczgr)

    ierr= NF90_DEF_DIM(nczgr, 'x', npiglo, idx)
    ierr= NF90_DEF_DIM(nczgr, 'y', npjglo, idy)
    ierr= NF90_DEF_DIM(nczgr, 'nav_lev', npkinp, idz)
    ierr= NF90_DEF_DIM(nczgr, 'time_counter', NF90_UNLIMITED,  idt)
 
    ierr=NF90_DEF_VAR(nczgr, 'nav_lon',       NF90_FLOAT, (/idx,idy/), id_navlon )
    ierr=NF90_DEF_VAR(nczgr, 'nav_lat',       NF90_FLOAT, (/idx,idy/), id_navlat )
    ierr=NF90_DEF_VAR(nczgr, 'nav_lev',       NF90_FLOAT, (/idz/)    , id_navlev )
    ierr=NF90_DEF_VAR(nczgr, 'time_counter',  NF90_FLOAT, (/idt/)    , id_time   )

    ierr=NF90_DEF_VAR(nczgr, 'gdept_1d',       NF90_FLOAT, (/idz,idt/), id_dept1d )
    ierr=NF90_DEF_VAR(nczgr, 'gdepw_1d',       NF90_FLOAT, (/idz,idt/), id_depw1d )

    ierr=NF90_DEF_VAR(nczgr, 'mbkt',          NF90_FLOAT, (/idx,idy,idt/), id_mbkt )
    ierr=NF90_DEF_VAR(nczgr, 'mbku',          NF90_FLOAT, (/idx,idy,idt/), id_mbku )
    ierr=NF90_DEF_VAR(nczgr, 'mbkv',          NF90_FLOAT, (/idx,idy,idt/), id_mbkv )
    ierr=NF90_DEF_VAR(nczgr, 'bathy_metry',   NF90_FLOAT, (/idx,idy,idt/), id_bathy )

    ierr=NF90_DEF_VAR(nczgr, 'e3t_0',         NF90_DOUBLE, (/idx,idy,idz,idt/), id_e3t )
    ierr=NF90_DEF_VAR(nczgr, 'e3w_0',         NF90_DOUBLE, (/idx,idy,idz,idt/), id_e3w )
    ierr=NF90_DEF_VAR(nczgr, 'e3u_0',         NF90_DOUBLE, (/idx,idy,idz,idt/), id_e3u )
    ierr=NF90_DEF_VAR(nczgr, 'e3v_0',         NF90_DOUBLE, (/idx,idy,idz,idt/), id_e3v )
    ierr=NF90_DEF_VAR(nczgr, 'e3uw_0',        NF90_DOUBLE, (/idx,idy,idz,idt/), id_e3uw )
    ierr=NF90_DEF_VAR(nczgr, 'e3vw_0',        NF90_DOUBLE, (/idx,idy,idz,idt/), id_e3vw )

    ierr=NF90_DEF_VAR(nczgr, 'tmask',         NF90_BYTE, (/idx,idy,idz,idt/), id_tmsk )
    ierr=NF90_DEF_VAR(nczgr, 'umask',         NF90_BYTE, (/idx,idy,idz,idt/), id_umsk )
    ierr=NF90_DEF_VAR(nczgr, 'vmask',         NF90_BYTE, (/idx,idy,idz,idt/), id_vmsk )

    ierr = NF90_ENDDEF(nczgr)

    ! put dimension related variables
    ierr = NF90_PUT_VAR(nczgr, id_navlon, navlon)
    ierr = NF90_PUT_VAR(nczgr, id_navlat, navlat)
    ierr = NF90_PUT_VAR(nczgr, id_navlev, gdept_1d )
    ierr = NF90_PUT_VAR(nczgr, id_time  , (/0./) )

    ! put 1D vertical variables (reference depth)
    ierr = NF90_PUT_VAR(nczgr, id_dept1d, gdept_1d, start=(/1,1/), count=(/npkinp,1/) )
    ierr = NF90_PUT_VAR(nczgr, id_depw1d, gdepw_1d, start=(/1,1/), count=(/npkinp,1/) )

  END SUBROUTINE CreateMeshZgrFile

END PROGRAM cdf_mshmsk_update_e3 
