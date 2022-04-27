!###############################################################################
!#                                                                             #
!# aed_pesticides.F90                                                          #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2022 -  The University of Western Australia                      #
!#                                                                             #
!#   AED is free software: you can redistribute it and/or modify               #
!#   it under the terms of the GNU General Public License as published by      #
!#   the Free Software Foundation, either version 3 of the License, or         #
!#   (at your option) any later version.                                       #
!#                                                                             #
!#   AED is distributed in the hope that it will be useful,                    #
!#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
!#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
!#   GNU General Public License for more details.                              #
!#                                                                             #
!#   You should have received a copy of the GNU General Public License         #
!#   along with this program.  If not, see <http://www.gnu.org/licenses/>.     #
!#   For use in a commercial seeting please contact the authors.               #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Originally created Dec 2021 by Matthew Hipsey and Thanh Hoang, UWA          #
!# Follow updates @ https://github.com/AquaticEcoDynamics/libaed-water         #
!#                                                                             #
!# For information about the module, please refer to the documentation         #
!# published here:                                                             #
!#                                                                             #
!# https://aquaticecodynamics.github.io/aed-science/index.html                 #
!#                                                                             #
!###############################################################################

#include "aed.h"


MODULE aed_pesticides
!-------------------------------------------------------------------------------
!  aed_pesticides --- pesticide contaminant model
!-------------------------------------------------------------------------------
   USE aed_core
   USE aed_util

   IMPLICIT NONE

   PRIVATE   ! By default make everything private
!
   PUBLIC aed_pesticides_data_t
!
!-------------------------------------------------------------------------------
!  %% NAMELIST   %%  pesticide_param_t
   TYPE pest_sorb_t
      CHARACTER(64)     :: pest_sorbent ! State variable name for sorbents
      AED_REAL          :: Kpst_sorb    ! Sorption factors for sorbents
   END TYPE pest_sorb_t
   TYPE pesticide_param_t
      CHARACTER(64)     :: name
      AED_REAL          :: Rhydrol
      AED_REAL          :: Rphoto
      AED_REAL          :: Ruptake
      AED_REAL          :: theta_hydrol
      AED_REAL          :: K_gpp
      AED_REAL          :: Fsed_pst
      AED_REAL          :: coef_light_kb_vis, coef_light_kb_uva, coef_light_kb_uvb  !-- Light inactivation
      INTEGER           :: sorption_model
      INTEGER           :: num_sorb
      TYPE(pest_sorb_t) :: sorbents(MAX_PSTC_SORB)
   END TYPE pesticide_param_t
   TYPE,extends(pesticide_param_t) :: pesticide_data_t
      INTEGER  :: id_sorb(MAX_PSTC_SORB)
      INTEGER  :: id_sorbv(MAX_PSTC_SORB)
   END TYPE
   !  %% END NAMELIST   %%  pesticide_param_t


   TYPE,extends(aed_model_data_t) :: aed_pesticides_data_t
      !# Variable identifiers
      INTEGER,ALLOCATABLE :: id_pstd(:)                    ! Column ID of pesticides
      INTEGER,ALLOCATABLE :: id_psta(:,:)                  ! Column ID of pesticides attached to components
      INTEGER,ALLOCATABLE :: id_psts(:)                    ! Column ID of pesticides in sediment
      INTEGER,ALLOCATABLE :: id_ss(:)                      ! Column ID of ss if chosen

      ! Diagnostic IDs for processes
      INTEGER,ALLOCATABLE :: id_atmvolat(:), id_sedflux(:), id_sorption(:),    &
                             id_photolysis(:), id_hydrolysis(:), id_uptake(:), &
                             id_settling(:), id_resus(:), id_total(:)

      INTEGER  :: id_oxy, id_pH,  id_doc, id_tss           ! Dependency ID
      INTEGER  :: id_tem, id_sal, id_gpp                   ! Environmental IDs (3D)
      INTEGER  :: id_par, id_nir, id_uva, id_uvb           ! Environmental IDs (3D)
      INTEGER  :: id_I_0                                   ! Environmental ID (2D)
      INTEGER  :: resuspension
      INTEGER  :: id_epsilon, id_taub

      !# Model parameters
      INTEGER  :: num_pesticides !, num_sorp
      TYPE(pesticide_data_t),DIMENSION(:),ALLOCATABLE :: pesticides

      INTEGER  :: pst_piston_model,pst_sorption_model
     !INTEGER  :: num_ss
     !AED_REAL,ALLOCATABLE  :: Fsed_pst(:)
     !AED_REAL,ALLOCATABLE  :: Rhydrol(:),Rphoto(:),Ruptake(:),theta_hydrol(:),K_gpp(:)
     !AED_REAL :: tau_0_min, kTau_0
     !AED_REAL,ALLOCATABLE :: epsilon(:), tau_0(:), tau_r(:), Ke_ss(:)
     !AED_REAL,ALLOCATABLE :: epsilonP(:), tauP_0(:)
     !AED_REAL,DIMENSION(:),ALLOCATABLE :: ss_set, ss_tau, ss_ke

      LOGICAL :: simSediment, simResuspension, simSorption, &
                 simVolatilisation, simPhotolysis, simUptake
      AED_REAL :: att_ts

     CONTAINS
         PROCEDURE :: define            => aed_define_pesticides
         PROCEDURE :: calculate         => aed_calculate_pesticides
         PROCEDURE :: calculate_surface => aed_calculate_surface_pesticides
         PROCEDURE :: calculate_benthic => aed_calculate_benthic_pesticides
         PROCEDURE :: equilibrate       => aed_equilibrate_pesticides
         PROCEDURE :: mobility          => aed_mobility_pesticides
         PROCEDURE :: light_extinction  => aed_light_extinction_pesticides
        !PROCEDURE :: delete            => aed_delete_pesticides
   END TYPE

!MODULE GLOBALS
   INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
                                              ! 1 = basic diagnostic outputs
                                              ! 2 = flux rates, and supporitng
                                              ! 3 = other metrics
                                              !10 = all debug & checking outputs

!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed_define_pesticides(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the pesticide model
!
!  Here, the aed_p_m namelist is read and the variables exported
!  by the model are registered with AED core.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed_pesticides_data_t),INTENT(inout) :: data

!
!LOCALS
   INTEGER  :: status
   INTEGER  :: i, pst_i, sorp_i
   CHARACTER(4) :: trac_name

!  %% NAMELIST   %% /aed_pesticides/
!  %% Last Checked 26/12/2021
   INTEGER  :: num_pesticides = 1
   INTEGER  :: the_pesticides(MAX_PATHO_TYPES)
   INTEGER  :: num_sorp = 0
   INTEGER  :: pst_piston_model = 1
   INTEGER  :: pst_sorption_model = 1
   INTEGER  :: resuspension = 0

   LOGICAL  :: simSediment       = .false.
   LOGICAL  :: simResuspension   = .false.
   LOGICAL  :: simSorption       = .false.
   LOGICAL  :: simVolatilisation = .false.
   LOGICAL  :: simPhotolysis     = .false.
   LOGICAL  :: simUptake         = .false.
   CHARACTER(len=64)  :: oxy_variable = ''
   CHARACTER(len=64)  :: gpp_variable = ''
   CHARACTER(len=128) :: dbase='aed_pesticide_pars.csv'

! From Module Globals
   LOGICAL  :: extra_diag = .FALSE.      !## Obsolete Use diag_level = 10
!  INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
!                                             ! 1 = basic diagnostic outputs
!                                             ! 2 = flux rates, and supporitng
!                                             ! 3 = other metrics
!                                             !10 = all debug & checking outputs
!  %% END NAMELIST   %% /aed_pesticides/

   NAMELIST /aed_pesticides/ num_pesticides, the_pesticides,       &
                             oxy_variable,                         &
                             simSediment,                          &
                             simResuspension, resuspension,        &
                             simVolatilisation, pst_piston_model,  &
                             simSorption,                          &
                             simPhotolysis,                        &
                             simUptake, gpp_variable,              &
                             dbase, diag_level

!-----------------------------------------------------------------------
!BEGIN
   print *,"        aed_pesticides configuration"
   print *,"          NOTE : UNDER DEVELOPMENT ... "

   !stop "Please disable the pesticide model in your configuration"


   ! Read the namelist
   read(namlst,nml=aed_pesticides,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_pesticides'

   IF ( extra_diag ) diag_level = 10

   ! Confirm configuration; defaults, requested config options and checks
   data%simSediment       = simSediment
   data%simPhotolysis     = simPhotolysis
   data%simResuspension   = simResuspension
   data%simSorption       = simSorption
   data%simVolatilisation = simVolatilisation
   data%simUptake         = simUptake

   IF(resuspension == 0) simResuspension = .false.
   IF(pst_piston_model == 0) simVolatilisation = .false.
   IF(pst_sorption_model == 0 .or. num_sorp == 0) simSorption = .false.
   IF(gpp_variable .EQ. '') simUptake = .false.

   ! Set module values to user provided numbers in the namelist
   data%resuspension = resuspension

   ! Store pesticide specific parameter values in module data type
   !   NB: all rates must be provided in values per day,
   !   and are converted here to values per second.
   CALL aed_pesticides_load_params(data, dbase, num_pesticides, the_pesticides)


   ! Register link to sorbent state variables
   DO pst_i = 1,num_pesticides
      !phy_i = 0
      DO sorp_i = 1,data%pesticides(pst_i)%num_sorb
          ! Find the sorbent concentration
          data%pesticides(pst_i)%id_sorb(sorp_i) = &
              aed_locate_variable(data%pesticides(pst_i)%sorbents(sorp_i)%pest_sorbent)
          ! Find the sorbent vvel
          data%pesticides(pst_i)%id_sorbv(sorp_i) = &
              aed_locate_variable(TRIM(data%pesticides(pst_i)%sorbents(sorp_i)%pest_sorbent)//'_vvel')

      ENDDO
   ENDDO

   ! Register state dependencies
   data%id_tss=-1 ; data%id_doc=-1 ; data%id_pH=-1 ; data%id_oxy=-1 ; data%id_gpp=-1
   IF (oxy_variable .NE. '') data%id_oxy = aed_locate_variable(oxy_variable)
   IF (gpp_variable .NE. '') data%id_gpp = aed_locate_variable(gpp_variable)

   ! Register environmental dependencies
   data%id_I_0 = aed_locate_sheet_global('par_sf')
   data%id_tem = aed_locate_global('temperature')
   data%id_sal = aed_locate_global('salinity')
   data%id_par = aed_locate_global('par')
   data%id_nir = aed_locate_global('nir')
   data%id_uva = aed_locate_global('uva')
   data%id_uvb = aed_locate_global('uvb')
   IF ( resuspension > 0 ) &
      data%id_taub = aed_locate_sheet_global('taub')


END SUBROUTINE aed_define_pesticides
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
INTEGER FUNCTION load_csv(dbase, pd)
!-------------------------------------------------------------------------------
   USE aed_csv_reader
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(len=*),INTENT(in) :: dbase
   TYPE(pesticide_data_t) :: pd(MAX_PSTC_TYPES)
!
!LOCALS
   INTEGER :: unit, nccols, ccol, dcol
   CHARACTER(len=32),POINTER,DIMENSION(:) :: csvnames
   CHARACTER(len=32) :: name_
   TYPE(AED_SYMBOL),DIMENSION(:),ALLOCATABLE :: values
   INTEGER :: idx_col = 0
   LOGICAL :: meh
   INTEGER :: ret = 0
!
!BEGIN
!-------------------------------------------------------------------------------
   unit = aed_csv_read_header(dbase, csvnames, nccols)
   IF (unit <= 0) THEN
      load_csv = -1
      RETURN !# No file found
   ENDIF

   ALLOCATE(values(nccols))

   DO WHILE ( aed_csv_read_row(unit, values) )
      DO ccol=2,nccols
         dcol = ccol - 1
         pd(dcol)%name = csvnames(ccol)
         CALL copy_name(values(1), name_)

         SELECT CASE (name_)
            CASE ('Rhydrol')           ; pd(dcol)%Rhydrol           = extract_double(values(ccol))
            CASE ('Rphoto')            ; pd(dcol)%Rphoto            = extract_double(values(ccol))
            CASE ('Ruptake')           ; pd(dcol)%Ruptake           = extract_double(values(ccol))
            CASE ('theta_hydrol')      ; pd(dcol)%theta_hydrol      = extract_double(values(ccol))
            CASE ('K_gpp')             ; pd(dcol)%K_gpp             = extract_double(values(ccol))
            CASE ('Fsed_pst')          ; pd(dcol)%Fsed_pst          = extract_double(values(ccol))
            CASE ('coef_light_kb_vis') ; pd(dcol)%coef_light_kb_vis = extract_double(values(ccol))
            CASE ('coef_light_kb_uva') ; pd(dcol)%coef_light_kb_uva = extract_double(values(ccol))
            CASE ('coef_light_kb_uvb') ; pd(dcol)%coef_light_kb_uvb = extract_double(values(ccol))
            CASE ('sorption_model')    ; pd(dcol)%sorption_model    = extract_integer(values(ccol))
            CASE ('num_sorb')          ; pd(dcol)%num_sorb          = extract_integer(values(ccol))

            CASE ('sorb(1)%pest_sorbent') ; CALL copy_name(values(ccol), pd(dcol)%sorbents(1)%pest_sorbent)
            CASE ('sorb(1)%Kpst_sorb') ; pd(dcol)%sorbents(1)%Kpst_sorb = extract_double(values(ccol))

            CASE ('sorb(2)%pest_sorbent') ; CALL copy_name(values(ccol), pd(dcol)%sorbents(2)%pest_sorbent)
            CASE ('sorb(2)%Kpst_sorb') ; pd(dcol)%sorbents(2)%Kpst_sorb = extract_double(values(ccol))

            CASE ('sorb(3)%pest_sorbent') ; CALL copy_name(values(ccol), pd(dcol)%sorbents(3)%pest_sorbent)
            CASE ('sorb(3)%Kpst_sorb') ; pd(dcol)%sorbents(3)%Kpst_sorb = extract_double(values(ccol))

            CASE DEFAULT ; print *, 'Unknown pesticide CSV parameter row "', TRIM(name_), '"'
         END SELECT
      ENDDO
   ENDDO

   meh = aed_csv_close(unit)
   !# don't care if close fails

   IF (ASSOCIATED(csvnames)) DEALLOCATE(csvnames)
   IF (ALLOCATED(values))    DEALLOCATE(values)

   load_csv = ret
END FUNCTION load_csv
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_pesticides_load_params(data, dbase, count, list)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_pesticides_data_t),INTENT(inout) :: data
   CHARACTER(len=*),INTENT(in) :: dbase
   INTEGER,INTENT(in) :: count
   INTEGER,INTENT(in) :: list(*)
!
!LOCALS
   INTEGER  :: status

   INTEGER  :: i,tfil,ns
   AED_REAL :: min_conc

   CHARACTER(4) :: pst_name

   TYPE(pesticide_data_t) :: pd(MAX_PSTC_TYPES)
   NAMELIST /pesticide_data/ pd   ! %% pesticide_param_t - see above
!-------------------------------------------------------------------------------
!BEGIN
    min_conc = 1e-8
    SELECT CASE (param_file_type(dbase))
       CASE (CSV_TYPE)
           status = load_csv(dbase, pd)
       CASE (NML_TYPE)
           print *,'NML file type for pesticides is not supported, please convert to CSV'
!          tfil = find_free_lun()
!          open(tfil,file=dbase, status='OLD',iostat=status)
!          IF (status /= 0) STOP 'Error opening namelist pesticide_data'
!          read(tfil,nml=pesticide_data,iostat=status)
!          close(tfil)
       CASE DEFAULT
           print *,'Unknown file type "',TRIM(dbase),'"'; status=1
    END SELECT
    IF (status /= 0) STOP 'Error reading namelist pesticide_data'

    data%num_pesticides = count
    ALLOCATE(data%pesticides(count))
    ALLOCATE(data%id_pstd(count))
    ALLOCATE(data%id_psta(count,10)) !need to get max num_sorb from pd
    IF (data%simSediment) THEN
       ALLOCATE(data%id_psts(count))
    ENDIF

    IF ( diag_level >= 2 ) THEN
       ALLOCATE(data%id_atmvolat(count))
       ALLOCATE(data%id_sedflux(count))
       ALLOCATE(data%id_settling(count))
       ALLOCATE(data%id_resus(count))
       ALLOCATE(data%id_sorption(count))
       ALLOCATE(data%id_photolysis(count))
       ALLOCATE(data%id_hydrolysis(count))
       ALLOCATE(data%id_uptake(count))
       ALLOCATE(data%id_total(count))
    ENDIF

    DO i=1,count
       ! Assign parameters from database to simulated groups
       data%pesticides(i) = pd(list(i))

       ! Unit adjustments from read-in parameters
       data%pesticides(i)%Fsed_pst = data%pesticides(i)%Fsed_pst / secs_per_day
       data%pesticides(i)%Rhydrol = data%pesticides(i)%Rhydrol / secs_per_day
       data%pesticides(i)%Rphoto = data%pesticides(i)%Rphoto / secs_per_day
       data%pesticides(i)%Ruptake = data%pesticides(i)%Ruptake / secs_per_day

       ! Register group as a state variable
       data%id_pstd(i) = aed_define_variable(                                  &
                             TRIM(data%pesticides(i)%name)//'_d',              &
                             'mmol/m3', 'pesticide dissolved concentration',   &
                             min_conc,                                         &
                             minimum=min_conc)

       ! Check if we need to registrer a variable for the sorbed fraction(s)
       IF (data%pesticides(i)%num_sorb > 0) THEN
         pst_name = '0'
         DO ns = 1, data%pesticides(i)%num_sorb
           pst_name(1:1) = CHAR(ICHAR('0') + ns)
           data%id_psta(i,ns) = aed_define_variable(                             &
                             TRIM(data%pesticides(i)%name)//'_'//TRIM(pst_name), &
                             'mmol/m3', 'pesticide sorbed concentration',        &
                             min_conc,                                           &
                             minimum=min_conc,                                   &
                             mobility = zero_)
          ENDDO
       ENDIF

       IF (data%simSediment) THEN
          data%id_psts(i) = aed_define_sheet_variable( TRIM(data%pesticides(i)%name)//'_s', 'mmol/m2', 'pesticides in sediment')
          PRINT *,'WARNING: simSediment is not complete'
       ENDIF

       !data%id_total(i) = aed_define_diag_variable( TRIM(data%pesticides(i)%name)//'_t', 'orgs/m3', 'total')
       IF ( diag_level >= 2 ) THEN
         data%id_atmvolat(i)  = &
                aed_define_sheet_diag_variable( TRIM(data%pesticides(i)%name)//'_atm', 'mmol/m2/d', 'volatilisation')
         data%id_sedflux(i)   = &
                aed_define_sheet_diag_variable( TRIM(data%pesticides(i)%name)//'_dsf', 'mmol/m2/d', 'dissolved sediment flux')
         data%id_resus(i)     = &
                aed_define_sheet_diag_variable( TRIM(data%pesticides(i)%name)//'_res', 'mmol/m2/d', 'resuspension flux')
         data%id_settling(i)  = aed_define_diag_variable( TRIM(data%pesticides(i)%name)//'_set', 'mmol/m3/d', 'settling rate')
         data%id_sorption(i)  = aed_define_diag_variable( TRIM(data%pesticides(i)%name)//'_srp', 'mmol/m3/d', 'sorption rate')
         data%id_photolysis(i)= aed_define_diag_variable( TRIM(data%pesticides(i)%name)//'_pht', 'mmol/m3/d', 'photolysis rate')
         data%id_hydrolysis(i)= aed_define_diag_variable( TRIM(data%pesticides(i)%name)//'_hyd', 'mmol/m3/d', 'hydrolysis rate')
         data%id_uptake(i)    = aed_define_diag_variable( TRIM(data%pesticides(i)%name)//'_upt', 'mmol/m3/d', 'uptake rate')
         data%id_total(i)     = &
                aed_define_diag_variable( TRIM(data%pesticides(i)%name)//'_tot', 'mmol/m3'  , 'total pesticide concentration')
       ENDIF
   ENDDO
END SUBROUTINE aed_pesticides_load_params
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_surface_pesticides(data,column,layer_idx)
!------------------------------------------------------------------------------+
! Air-water exchange (volatilisation) for the aed pesticide model
!------------------------------------------------------------------------------+
!ARGUMENTS
   CLASS (aed_pesticides_data_t),INTENT(in) :: data
   TYPE  (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp, salt, windHt, wind, vel, depth
   ! State
   AED_REAL :: volat, pst
   ! Temporary variables
   INTEGER  :: pst_i
   AED_REAL :: k600
!
!------------------------------------------------------------------------------+
!BEGIN

  temp = _STATE_VAR_(data%id_tem)      ! local temperature
  salt = _STATE_VAR_(data%id_sal)       ! local salinity

  windHt = 1
  wind = 1
  vel = 0.1
  depth = 1

  DO pst_i=1,data%num_pesticides

    volat = zero_

    !-----------------------------------------------
    ! Compute necessary piston velocity and air-sea flux
    IF( data%simVolatilisation ) THEN

      pst = _STATE_VAR_(data%id_pstd(pst_i))

      k600 = aed_gas_piston_velocity(windHt,wind,temp,salt,                    &
        vel=vel,depth=depth,schmidt_model=2,piston_model=data%pst_piston_model)

      volat = k600 * pst
    ENDIF

    !-----------------------------------------------
    ! Set surface exchange value (mmmol/m2/s) for AED ODE solution
    _FLUX_VAR_T_(data%id_pstd(pst_i)) = volat

    !-----------------------------------------------
    ! Set surface exchange value (mmmol/m2/d) as a diagnostic
    _DIAG_VAR_S_(data%id_atmvolat(pst_i)) = volat * secs_per_day

  ENDDO

END SUBROUTINE aed_calculate_surface_pesticides
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_pesticides(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of pesticide biogeochemical model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_pesticides_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL :: pth_f, pth_a, pest_d, pest_a
   AED_REAL :: temp,salinity,oxy,pH,doc
   AED_REAL :: Io,par,uva,uvb
   AED_REAL :: hydrolysis, photolysis, uptake
   AED_REAL :: f_AOC,f_pH,f_DO,phi,lightBW,phstar,att_frac
   AED_REAL :: f_upt, f_pht

   INTEGER  :: pst_i, sorp_i, pst_s

!-------------------------------------------------------------------------------
!BEGIN


   ! Retrieve current environmental conditions
   temp = _STATE_VAR_(data%id_tem)     ! local temperature
   salinity = _STATE_VAR_(data%id_sal) ! local salinity
   IF (data%id_oxy>0) THEN  ! & use_oxy
      oxy = _STATE_VAR_(data%id_oxy)   ! local oxygen
   ELSE
      oxy = 10.0 !mg/L
   ENDIF

   !doc = _STATE_VAR_(data%id_doc)    ! local DOC
   !ph = _STATE_VAR_(data%id_ph)      ! local pH
   phstar = 0.0                       ! abs(ph-7.)

   ! Get light bandwidth intensities
   Io = _STATE_VAR_S_(data%id_I_0)    ! surface short wave radiation
   par = _STATE_VAR_(data%id_par)     ! local photosynthetically active radiation (45% of sw)
   IF ( data%id_uva > 0 ) THEN
      uva = _STATE_VAR_(data%id_uva)
   ELSE
      uva = (par/0.45)*0.03           ! uva is 3% of sw (Kirk 1994)
   ENDIF
   IF ( data%id_uvb > 0 ) THEN
      uvb = _STATE_VAR_(data%id_uvb)
   ELSE
      uvb = (par/0.45)*0.003          ! uvb is 0.3% of sw
   ENDIF


   DO pst_i=1,data%num_pesticides

      !-----------------------------------------------------------------
      ! RETREIVE THIS PESTICIDE GROUP

      pest_d = _STATE_VAR_(data%id_pstd(pst_i))
    !  IF ( data%num_sorp > 0 ) THEN
    !    pest_a(pst_s) = _STATE_VAR_(data%id_psta(pst_i,sorp_i))
    !  END IF

      !-----------------------------------------------------------------
      ! COMPUTE PESTICIDE FLUX RATES

      hydrolysis = zero_
      photolysis = zero_
      uptake     = zero_

      _DIAG_VAR_(data%id_hydrolysis(pst_i)) = zero_
      _DIAG_VAR_(data%id_photolysis(pst_i)) = zero_
      _DIAG_VAR_(data%id_uptake(pst_i)) = zero_
      _DIAG_VAR_(data%id_total(pst_i))  = zero_

      ! 1. Breakdown under ambient conditions

      hydrolysis = data%pesticides(pst_i)%Rhydrol * (data%pesticides(pst_i)%theta_hydrol**(temp-20.0))

      ! 2. Sunlight breakdown

      IF ( data%simPhotolysis ) THEN
        f_pht = par / ( par + 500. )
        photolysis = data%pesticides(pst_i)%coef_light_kb_vis * f_pht

        f_pht = uva / ( uva + 50. )
        photolysis = photolysis+ data%pesticides(pst_i)%coef_light_kb_uva * f_pht

        f_pht = uvb / ( uvb + 5. )
        photolysis = photolysis+ data%pesticides(pst_i)%coef_light_kb_uvb * f_pht

       !       photolysis = photo(vis,cdom,1) + photo(uva,cdom,2) + photo(uvb,cdom,3)
       !       !# Limit photolysis to 90% of doc pool within 1 hour
       !       IF(photolysis > 0.9*docr/3.6e3) photolysis = 0.9*docr/3.6e3
      ENDIF


      ! 3. Biological uptake
      IF ( data%simUptake ) THEN
        f_upt = _DIAG_VAR_(data%id_GPP) / ( _DIAG_VAR_(data%id_GPP) + data%pesticides(pst_i)%K_gpp )
        uptake = data%pesticides(pst_i)%Ruptake * f_upt
      ENDIF


      !-----------------------------------------------------------------
      ! SET TEMPORAL DERIVATIVES FOR ODE SOLVER

      ! Pesticide breakdown and uptake
      _FLUX_VAR_(data%id_pstd(pst_i)) = _FLUX_VAR_(data%id_pstd(pst_i))        &
                                    -  (hydrolysis + photolysis + uptake)*pest_d

      DO sorp_i=1,data%pesticides(pst_i)%num_sorb
        pest_a = _STATE_VAR_(data%id_psta(pst_i,sorp_i))
        _FLUX_VAR_(data%id_psta(pst_i,sorp_i)) = _FLUX_VAR_(data%id_psta(pst_i,sorp_i))        &
                                      -  (hydrolysis + photolysis/2.)*pest_a
        _DIAG_VAR_(data%id_total(pst_i)) = _DIAG_VAR_(data%id_total(pst_i)) + pest_a
      ENDDO
      _DIAG_VAR_(data%id_total(pst_i)) = _DIAG_VAR_(data%id_total(pst_i)) + pest_d


      !-----------------------------------------------------------------
      ! SET DIAGNOSTICS

      IF ( diag_level >= 2 ) THEN
        _DIAG_VAR_(data%id_hydrolysis(pst_i)) = &
            hydrolysis * (_DIAG_VAR_(data%id_total(pst_i))+pest_d) * secs_per_day
        _DIAG_VAR_(data%id_photolysis(pst_i)) = &
            photolysis * (_DIAG_VAR_(data%id_total(pst_i))/2.+pest_d) * secs_per_day
        _DIAG_VAR_(data%id_uptake(pst_i)) = uptake * pest_d * secs_per_day
      ENDIF

   ENDDO

END SUBROUTINE aed_calculate_pesticides
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_benthic_pesticides(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic sedimentation of pesticide.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_pesticides_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL :: ss       ! State
   INTEGER  :: pst_i, ss_i
   AED_REAL :: diss_flux, sett_flux, kin_flux
   AED_REAL :: pest_d, pest_sed_d

   AED_REAL :: bottom_stress, dummy_tau

   ! Parameters
!
!-------------------------------------------------------------------------------
!BEGIN
   diss_flux = zero_

   IF (data%simSediment) THEN
      ! Dynamic sediment pool of pesticide that increase and decrease

      DO pst_i=1,data%num_pesticides
         ! Retrieve current (local) state variable values.
         pest_d = _STATE_VAR_(data%id_pstd(pst_i))

         pest_sed_d = _STATE_VAR_S_(data%id_psts(pst_i)) ! pesticide (sediment pool)

         diss_flux = data%pesticides(pst_i)%Fsed_pst

         ! Org flux to / from the sediment (orgs/m2/s)
  !       _FLUX_VAR_B_(data%id_ps(pst_i)) = _FLUX_VAR_B_(data%id_ps(pst_i)) - diss_flux

         ! Add to respective pools in water (free/attached)
         _FLUX_VAR_(data%id_pstd(pst_i)) = _FLUX_VAR_(data%id_pstd(pst_i)) + diss_flux

         IF ( diag_level >= 2 ) &
           _DIAG_VAR_S_ (data%id_sedflux(pst_i)) = diss_flux * secs_per_day

         IF ( diag_level >= 2 ) &
           _DIAG_VAR_S_ (data%id_resus(pst_i)) = zero_ * secs_per_day
      ENDDO
   ELSE
      ! No sediment pool is resolved, but we will still predict generic diss flux
      DO pst_i=1,data%num_pesticides
        diss_flux = data%pesticides(pst_i)%Fsed_pst

        ! Flux from the sediment
        _FLUX_VAR_B_(data%id_pstd(pst_i)) = _FLUX_VAR_B_(data%id_pstd(pst_i)) + diss_flux

        IF ( diag_level >= 2 ) &
          _DIAG_VAR_S_ (data%id_sedflux(pst_i)) = diss_flux * secs_per_day

        IF ( diag_level >= 2 ) &
          _DIAG_VAR_S_ (data%id_resus(pst_i)) = zero_ * secs_per_day
      ENDDO
   ENDIF
END SUBROUTINE aed_calculate_benthic_pesticides
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_equilibrate_pesticides(data,column,layer_idx)
!------------------------------------------------------------------------------+
! Update partitioning of pesticides between dissolved and particulate pools;
! updated after kinetic transformations are applied
!------------------------------------------------------------------------------+
!ARGUMENTS
   CLASS (aed_pesticides_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp, tss

   ! State
   AED_REAL :: pest_d_preequil, pest_t, pest_d, pest_s(10)

   ! Temporary variables
   AED_REAL :: sorbents(10), Kpstp(10)
   INTEGER  :: pst_i, sorp_i

   AED_REAL, PARAMETER :: dt = 900  ! Needs to be linked to aed_core/common

!-------------------------------------------------------------------------------
!BEGIN
   IF(.NOT. data%simSorption) RETURN

   ! Retrieve current environmental conditions for the cell.
   temp = _STATE_VAR_(data%id_tem)    ! local temperature

   pest_d = zero_ ; pest_t = zero_ ; pest_s(:) = zero_

   DO pst_i=1,data%num_pesticides

     ! Record dissolved conc, before the sorption algorithm
     pest_d_preequil = _STATE_VAR_(data%id_pstd(pst_i))

     ! Find the total pesticide across all forms
     pest_t = _STATE_VAR_(data%id_pstd(pst_i))
     DO sorp_i=1,data%pesticides(pst_i)%num_sorb

       sorbents(sorp_i) = _STATE_VAR_(data%pesticides(pst_i)%id_sorb(sorp_i))

       pest_t = pest_t + _STATE_VAR_(data%id_psta(pst_i,sorp_i))

       Kpstp(sorp_i) = data%pesticides(pst_i)%sorbents(sorp_i)%Kpst_sorb
     ENDDO

     ! Re-distribute based on component concentrations
     pest_d = pesticide_sorption(temp,sorbents,pest_t,pest_s,Kpstp,data%pesticides(pst_i)%sorption_model)

     ! Update core data arrays
     _STATE_VAR_(data%id_pstd(pst_i))    = pest_d              ! Dissolved
     DO sorp_i=1,data%pesticides(pst_i)%num_sorb
      _STATE_VAR_(data%id_psta(pst_i,sorp_i)) = pest_s(sorp_i) ! Adsorped to particle group
     ENDDO

     ! Update diagnostic
     _DIAG_VAR_(data%id_sorption(pst_i)) = ((pest_d - pest_d_preequil)/dt)* secs_per_day
   ENDDO


END SUBROUTINE aed_equilibrate_pesticides
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_mobility_pesticides(data,column,layer_idx,mobility)
!-------------------------------------------------------------------------------
! Get the vertical movement values based on linked particulates (sorbents)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_pesticides_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: mobility(:)
!
!LOCALS
   INTEGER  :: sorp_i,pst_i
   AED_REAL :: temp
   AED_REAL :: sorbent_vvel(10)
!
!-------------------------------------------------------------------------------
!BEGIN
   temp = _STATE_VAR_(data%id_tem)

   ! Set velocity of sorped pesticides, if simulated.
   DO pst_i=1,data%num_pesticides

     IF ( diag_level >= 2 ) _DIAG_VAR_(data%id_settling(pst_i)) = zero_

     DO sorp_i=1,data%pesticides(pst_i)%num_sorb
       sorbent_vvel(sorp_i) = _DIAG_VAR_(data%pesticides(pst_i)%id_sorbv(sorp_i)) / secs_per_day
       mobility(data%id_psta(pst_i,sorp_i)) = sorbent_vvel(sorp_i)

       IF ( diag_level >= 2 ) THEN
         _DIAG_VAR_(data%id_settling(pst_i)) = _DIAG_VAR_(data%id_settling(pst_i)) &
            + _STATE_VAR_(data%id_psta(pst_i,sorp_i)) * sorbent_vvel(sorp_i)  * secs_per_day
       ENDIF
     ENDDO
   ENDDO
END SUBROUTINE aed_mobility_pesticides
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_light_extinction_pesticides(data,column,layer_idx,extinction)
!-------------------------------------------------------------------------------
! Get the light extinction coefficient due to ss variables in this module
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_pesticides_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: extinction
!
!LOCALS
   AED_REAL :: ss
   INTEGER  :: ss_i, pst_i
!
!-----------------------------------------------------------------------
!BEGIN
   RETURN
END SUBROUTINE aed_light_extinction_pesticides
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION pesticide_sorption(temp,sorbents,pest_t,pest_s,Kpstp_,sorption_model) &
                                                                  RESULT(pest_d)
!-------------------------------------------------------------------------------
! Function to partition pesticide concentration amongst several sorbents
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER, INTENT(in)    :: sorption_model
   AED_REAL,INTENT(in)    :: temp,pest_t,Kpstp_(:),sorbents(:)
   AED_REAL,INTENT(inout) :: pest_s(:)
   AED_REAL :: pest_d
!
!LOCALS
   INTEGER   :: sorp_i
   AED_REAL  :: PSTtot, SorbentConc, pstpar , pstdis, Kpstp
!
!-------------------------------------------------------------------------------
!BEGIN

   ! Default return variables to all dissolved
   pest_s(:) = zero_
   pest_d = pest_t

   ! Set local values of total PST and SS etc
   PSTtot      = MAX(1e-10, pest_t )             ! Co in Chao (mg)
   SorbentConc = MAX(1e-10, sum(sorbents(:)) )   ! s in Chao  (mg = mol/L * g/mol * mg/g)

   IF(SorbentConc < 1e-5) RETURN

   IF(sorption_model == 1) THEN
     !-----------------------------------------------------
     ! This is the model for sorption from Ji 2008:
     !
     ! Ji, Z-G. 2008. Hydrodynamics and Water Quality. Wiley Press.

     Kpstp = Kpstp_(1)  ! simplification for now, 1st value is used

     pstpar = (Kpstp*SorbentConc) / (one_+Kpstp*SorbentConc) * PSTtot
     pstdis = one_ / (one_+Kpstp*SorbentConc) * PSTtot

     ! Set return variables for each sorbent and dissolved
     DO sorp_i = 1,SIZE(sorbents)
       pest_s(sorp_i) = pstpar * (sorbents(sorp_i)/SorbentConc)
     ENDDO
     pest_d = pstdis

   ENDIF


END FUNCTION pesticide_sorption
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE aed_pesticides
