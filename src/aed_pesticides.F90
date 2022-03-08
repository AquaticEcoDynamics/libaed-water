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
!# https://aquaticecodynamics.github.io/aed-science/index.html                 #                                     #
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
TYPE pest_sorb_t
   !State variable name for pesticide pest
   CHARACTER(64) :: pest_sorbent
   !Sorption factors for pesticide onto sorbents
   AED_REAL      :: Kpst_sorb
END TYPE pest_sorb_t

!----------------------------------------------------------------------------
!  %% NAMELIST   %%  pesticide_param_t
   TYPE pesticide_param_t
      CHARACTER(64) :: name
      AED_REAL          :: Rhydrol
      AED_REAL          :: Rphoto
      AED_REAL          :: Ruptake
      AED_REAL          :: theta_hydrol
      AED_REAL          :: K_gpp
      AED_REAL          :: coef_light_kb_vis, coef_light_kb_uva, coef_light_kb_uvb  !-- Light inactivation
      INTEGER           :: num_sorb
      TYPE(pest_sorb_t) :: sorbents(MAX_ZOOP_PREY)
   END TYPE pesticide_param_t
!  %% END NAMELIST   %%  pesticide_param_t

TYPE,extends(pesticide_param_t) :: pesticide_data_t
   INTEGER  :: id_sorb(MAX_ZOOP_PREY)
  !INTEGER  :: id_phyIN(MAX_ZOOP_PREY),id_phyIP(MAX_ZOOP_PREY)
END TYPE




   TYPE,extends(aed_model_data_t) :: aed_pesticides_data_t
      !# Variable identifiers
      INTEGER,ALLOCATABLE :: id_pstd(:)                    ! Column ID of pesticides
      INTEGER,ALLOCATABLE :: id_psta(:,:)                  ! Column ID of pesticides attached to components
      INTEGER,ALLOCATABLE :: id_psts(:)                    ! Column ID of pesticides in sediment
      INTEGER,ALLOCATABLE :: id_ss(:)                      ! Column ID of ss if chosen

      ! Diagnostic IDs for processes
      INTEGER,ALLOCATABLE :: id_atmvolat(:), id_sedflux(:), id_sorption(:), &
                             id_photolysis(:), id_hydrolysis(:), id_uptake(:)

      INTEGER  :: id_oxy, id_pH,  id_doc, id_tss           ! Dependency ID
      INTEGER  :: id_tem, id_sal, id_gpp                   ! Environmental IDs (3D)
      INTEGER  :: id_par, id_nir, id_uva, id_uvb           ! Environmental IDs (3D)
      INTEGER  :: id_I_0                                   ! Environmental ID (2D)
      INTEGER  :: resuspension
      INTEGER  :: id_epsilon, id_taub

      !# Model parameters
      INTEGER  :: num_pesticides, num_sorp
      TYPE(pesticide_param_t),DIMENSION(:),ALLOCATABLE :: pesticides

      INTEGER  :: num_ss
      INTEGER  :: pst_piston_model
      AED_REAL,ALLOCATABLE  :: Fsed_pst(:), sssss(:)
      AED_REAL,ALLOCATABLE  :: Rhydrol(:),Rphoto(:),Ruptake(:),theta_hydrol(:),K_gpp(:)
      AED_REAL :: tau_0_min, kTau_0
      AED_REAL,ALLOCATABLE :: epsilon(:), tau_0(:), tau_r(:), Ke_ss(:)
      AED_REAL,ALLOCATABLE :: epsilonP(:), tauP_0(:)
      AED_REAL,DIMENSION(:),ALLOCATABLE :: ss_set, ss_tau, ss_ke

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
   INTEGER  :: i
   CHARACTER(4) :: trac_name

!  %% NAMELIST   %% /aed_pesticides/
!  %% Last Checked 26/12/2021
   INTEGER  :: num_pesticides, num_sorp
   INTEGER  :: the_pesticides(MAX_PATHO_TYPES)
   INTEGER  :: num_ss = 0
   AED_REAL :: ss_set(MAX_PATHO_TYPES)=zero_
   AED_REAL :: ss_tau(MAX_PATHO_TYPES)=one_
   AED_REAL :: ss_ke(MAX_PATHO_TYPES) =zero_
   AED_REAL :: ss_initial = zero_
   AED_REAL :: epsilon(100)
   AED_REAL :: tau_0(100)
   AED_REAL :: tau_r(100)
   AED_REAL :: tau_0_min
   AED_REAL :: Ktau_0
   AED_REAL :: att_ts = (1./(secs_per_day*7.)) ! attachment fraction re-equilibrates over a week
   INTEGER  :: resuspension
   LOGICAL  :: simSediment = .FALSE.
   CHARACTER(len=64)  :: oxy_variable = ''
   CHARACTER(len=64)  :: gpp_variable = ''
   CHARACTER(len=128) :: dbase='aed_pesticide_pars.nml'

! From Module Globals
   LOGICAL  :: extra_diag = .FALSE.      !## Obsolete Use diag_level = 10
!  INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
!                                             ! 1 = basic diagnostic outputs
!                                             ! 2 = flux rates, and supporitng
!                                             ! 3 = other metrics
!                                             !10 = all debug & checking outputs
!  %% END NAMELIST   %% /aed_pesticides/

   NAMELIST /aed_pesticides/ num_pesticides, the_pesticides, &
                             num_sorp, oxy_variable,         &
            num_ss, ss_set, ss_tau, ss_ke, simSediment,      &
            resuspension, epsilon, tau_0, tau_0_min, Ktau_0, &
            dbase, extra_diag, att_ts,diag_level, &
            gpp_variable
!-----------------------------------------------------------------------
!BEGIN
   print *,"        aed_pesticides configuration"
   print *,"          NOTE : UNDER DEVELOPMENT ... STOPPING"

   !stop "Please disable the pesticide model in your configuration"

   ! Read the namelist
   read(namlst,nml=aed_pesticides,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_pesticides'

   IF ( extra_diag ) diag_level = 10

   data%tau_0_min = tau_0_min
   data%Ktau_0 = Ktau_0
   data%att_ts = att_ts

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   CALL aed_pesticides_load_params(data, dbase, num_pesticides, the_pesticides)

   IF ( num_ss > 0 ) THEN
      ALLOCATE(data%id_ss(num_ss))
      ALLOCATE(data%ss_set(num_ss)) ; data%ss_set(1:num_ss) = ss_set(1:num_ss)/secs_per_day
      ALLOCATE(data%ss_ke(num_ss))  ; data%ss_ke(1:num_ss)  = ss_ke(1:num_ss)
      ALLOCATE(data%ss_tau(num_ss)) ; data%ss_tau(1:num_ss) = ss_tau(1:num_ss)

      ALLOCATE(data%epsilon(num_ss)); data%epsilon(1:num_ss) = epsilon(1:num_ss)
      ALLOCATE(data%tau_0(num_ss))  ; data%tau_0(1:num_ss)   = tau_0(1:num_ss)

      trac_name = 'ss0'
      ! Register state variables
      DO i=1,num_ss
         trac_name(3:3) = CHAR(ICHAR('0') + i)
         data%id_ss(i) = aed_define_variable(TRIM(trac_name),'g/m**3','pesticide ss', &
                                                  ss_initial,minimum=zero_)
      ENDDO
   ENDIF
   IF ( resuspension == 2 ) THEN
      data%id_epsilon =  aed_define_sheet_diag_variable('epsilon','g/m**2/s', 'Resuspension rate')
   ENDIF

   ! Register state dependancies
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

   data%resuspension = resuspension

END SUBROUTINE aed_define_pesticides
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
INTEGER FUNCTION load_csv(dbase, pd)
!-------------------------------------------------------------------------------
   USE aed_csv_reader
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(len=*),INTENT(in) :: dbase
   TYPE(pesticide_param_t) :: pd(MAX_PHYTO_TYPES)
!
!LOCALS
   INTEGER :: unit, nccols, ccol
   CHARACTER(len=32),POINTER,DIMENSION(:) :: csvnames
   CHARACTER(len=32) :: name
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
         pd(ccol)%name = csvnames(ccol)

         CALL copy_name(values(1), name)
         SELECT CASE (name)



            CASE ('coef_light_kb_vis')   ; pd(ccol)%coef_light_kb_vis   = extract_double(values(ccol))
            CASE ('coef_light_kb_uva')   ; pd(ccol)%coef_light_kb_uva   = extract_double(values(ccol))
            CASE ('coef_light_kb_uvb')   ; pd(ccol)%coef_light_kb_uvb   = extract_double(values(ccol))

            CASE DEFAULT ; print *, 'Unknown row "', TRIM(name), '"'
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

   TYPE(pesticide_param_t) :: pd(MAX_PATHO_TYPES)
   NAMELIST /pesticide_data/ pd   ! %% pesticide_param_t - see above
!-------------------------------------------------------------------------------
!BEGIN
    min_conc = 1e-8
    SELECT CASE (param_file_type(dbase))
       CASE (CSV_TYPE)
           status = load_csv(dbase, pd)
       CASE (NML_TYPE)
           tfil = find_free_lun()
           open(tfil,file=dbase, status='OLD',iostat=status)
           IF (status /= 0) STOP 'Error opening namelist pesticide_data'
           read(tfil,nml=pesticide_data,iostat=status)
           close(tfil)
       CASE DEFAULT
           print *,'Unknown file type "',TRIM(dbase),'"'; status=1
    END SELECT
    IF (status /= 0) STOP 'Error reading namelist pesticide_data'


    data%num_pesticides = count
    ALLOCATE(data%pesticides(count))
    ALLOCATE(data%id_pstd(count))
    ALLOCATE(data%id_psta(count,data%num_sorp))
    IF (data%simSediment) THEN
       ALLOCATE(data%id_psts(count))
    ENDIF
    IF ( diag_level >= 2 ) THEN
       ALLOCATE(data%id_atmvolat(count))
       ALLOCATE(data%id_sedflux(count))
       ALLOCATE(data%id_sorption(count))
       ALLOCATE(data%id_photolysis(count))
       ALLOCATE(data%id_hydrolysis(count))
       ALLOCATE(data%id_uptake(count))
    ENDIF

    DO i=1,count
       ! Assign parameters from database to simulated groups
       data%pesticides(i) = pd(list(i))
       !data%pesticides(i)%coef_sett_w_path = data%pesticides(i)%coef_sett_w_path/secs_per_day
       !data%pesticides(i)%coef_mort_kd20 = data%pesticides(i)%coef_mort_kd20/secs_per_day

       ! Register group as a state variable
       data%id_pstd(i) = aed_define_variable(                                  &
                             TRIM(data%pesticides(i)%name)//'_d',            &
                             'mmol/m**3', 'pesticide dissolved concetration',  &
                             min_conc,                                         &
                            ! pd(list(i))%p_initial,                           &
                             minimum=min_conc)

       ! Check if we need to registrer a variable for the attached fraction
       IF (data%num_sorp > zero_) THEN
         pst_name = '0'
         DO ns = 1, data%num_sorp
           pst_name(1:1) = CHAR(ICHAR('0') + ns)
           data%id_psta(i,ns) = aed_define_variable(                                &
                             TRIM(data%pesticides(i)%name)//'_'//TRIM(pst_name),            &
                             'mmol/m**3', 'pesticide sorbed concentration',                &
                             min_conc,                                         &
                            ! pd(list(i))%p_initial,                          &
                             minimum=min_conc,                                 &
                             !minimum=pd(list(i))%p0,                         &
                             mobility = zero_)
          ENDDO
       ENDIF


       IF (data%simSediment) THEN
          data%id_psts(i) = aed_define_sheet_variable( TRIM(data%pesticides(i)%name)//'_s', 'orgs/m2', 'pesticides in sediment')
          PRINT *,'WARNING: simSediment is not complete'
       ENDIF

       !data%id_total(i) = aed_define_diag_variable( TRIM(data%pesticides(i)%name)//'_t', 'orgs/m3', 'total')
       IF ( diag_level >= 2 ) THEN
         data%id_atmvolat(i)  = aed_define_sheet_diag_variable( TRIM(data%pesticides(i)%name)//'_vol', 'mmol/m2/day', 'growth')
         data%id_sedflux(i)   = aed_define_sheet_diag_variable( TRIM(data%pesticides(i)%name)//'_dsf', 'mmol/m2/day', 'growth')
         data%id_sorption(i)  = aed_define_diag_variable( TRIM(data%pesticides(i)%name)//'_srp', 'mmol/m3/day', 'growth')
         data%id_photolysis(i)= aed_define_diag_variable( TRIM(data%pesticides(i)%name)//'_pht', 'mmol/m3/day', 'growth')
         data%id_hydrolysis(i)= aed_define_diag_variable( TRIM(data%pesticides(i)%name)//'_hyd', 'mmol/m3/day', 'growth')
         data%id_uptake(i)    = aed_define_diag_variable( TRIM(data%pesticides(i)%name)//'_upt', 'mmol/m3/day', 'growth')

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
   AED_REAL :: pth_f, pth_a, pest_d, pest_a(data%num_sorp)
   AED_REAL :: temp,salinity,oxy,pH,doc
   AED_REAL :: Io,par,uva,uvb
   AED_REAL :: hydrolysis, photolysis, uptake
   AED_REAL :: f_AOC,f_pH,f_DO,phi,lightBW,phstar,att_frac
   AED_REAL :: f_upt

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
      IF ( data%num_sorp > 0 ) THEN
        pest_a(pst_s) = _STATE_VAR_(data%id_psta(pst_i,sorp_i))
      END IF

      !-----------------------------------------------------------------
      ! COMPUTE PESTICIDE FLUX RATES

      hydrolysis = zero_
      photolysis = zero_
      uptake     = zero_

      ! 1. Breakdown under ambient conditions

      hydrolysis = data%Rhydrol(pst_i) * (data%theta_hydrol(pst_i)**(temp-20.0))

      ! 2. Sunlight breakdown

      IF ( data%simPhotolysis ) THEN
        photolysis = data%Rphoto(pst_i) ! something about light here
       !       photolysis = photo(vis,cdom,1) + photo(uva,cdom,2) + photo(uvb,cdom,3)
       !       !# Limit photolysis to 90% of doc pool within 1 hour
       !       IF(photolysis > 0.9*docr/3.6e3) photolysis = 0.9*docr/3.6e3
      ENDIF


      ! 3. Biological uptake
      IF ( data%simUptake ) THEN
        f_upt = _DIAG_VAR_(data%id_GPP) / ( _DIAG_VAR_(data%id_GPP) + data%K_gpp(pst_i) )
        uptake = data%Ruptake(pst_i) * f_upt
      ENDIF


      !-----------------------------------------------------------------
      ! SET TEMPORAL DERIVATIVES FOR ODE SOLVER

      ! Pesticide breakdown and uptake
      _FLUX_VAR_(data%id_pstd(pst_i)) = _FLUX_VAR_(data%id_pstd(pst_i))        &
                                    -  (hydrolysis + photolysis + uptake)*pest_d


      !-----------------------------------------------------------------
      ! SET DIAGNOSTICS
    !  _DIAG_VAR_(data%id_total(pst_i)) =  pth_f + pth_a + pth_d
      IF ( diag_level >= 10 ) THEN
    !     _DIAG_VAR_(data%id_growth(pst_i)) =
    !     _DIAG_VAR_(data%id_sunlight(pst_i)) =
    !     _DIAG_VAR_(data%id_mortality(pst_i)) =
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

      diss_flux = data%Fsed_pst(pst_i)


      ! Org flux to / from the sediment (orgs/m2/s)
  !    _FLUX_VAR_B_(data%id_ps(pst_i)) = _FLUX_VAR_B_(data%id_ps(pst_i)) - diss_flux

      ! Add to respective pools in water (free/attached)
      _FLUX_VAR_(data%id_pstd(pst_i)) = _FLUX_VAR_(data%id_pstd(pst_i)) + diss_flux

      IF ( diag_level >= 2 ) THEN
        _DIAG_VAR_S_ (data%id_sedflux(pst_i)) = diss_flux * secs_per_day
      ENDIF
   ENDDO

ELSE

   ! No sediment pool is resolved, but we will still predict generic diss flux
   DO pst_i=1,data%num_pesticides

     diss_flux = data%Fsed_pst(pst_i)

     ! Flux from the sediment
     _FLUX_VAR_B_(data%id_pstd(pst_i)) = _FLUX_VAR_B_(data%id_pstd(pst_i)) + diss_flux

     IF ( diag_level >= 2 ) THEN
       _DIAG_VAR_S_ (data%id_sedflux(pst_i)) = diss_flux * secs_per_day
     ENDIF

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
   AED_REAL :: frp,frpads,pH

   ! Temporary variables
   AED_REAL :: pest_t, pest_d, pest_s(data%num_sorp)
   AED_REAL :: sorbents(data%num_sorp)
   INTEGER  :: pst_i, sorp_i

!-------------------------------------------------------------------------------
!BEGIN
   IF(.NOT. data%simSorption) RETURN

   ! Retrieve current environmental conditions for the cell.
   temp = _STATE_VAR_(data%id_tem)    ! local temperature

   pest_d = zero_ ; pest_t = zero_ ; pest_s(:) = zero_

   DO pst_i=1,data%num_pesticides

     pest_t = _STATE_VAR_(data%id_pstd(pst_i))
     DO sorp_i=1,data%num_sorp
       sorbents(sorp_i) = 0.

       pest_t = pest_t + _STATE_VAR_(data%id_psta(pst_i,sorp_i))
     ENDDO

     ! distribute based on component concentrations
     pest_d = pesticide_sorption(temp,sorbents,pest_t,pest_s,data%sssss(pst_i))

     ! Update core data arrays
     _STATE_VAR_(data%id_pstd(pst_i))    = pest_d    ! Dissolved
     DO sorp_i=1,data%num_sorp
      _STATE_VAR_(data%id_psta(pst_i,sorp_i)) = pest_s(sorp_i)    ! Adsorped to particle group
     ENDDO
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
   AED_REAL :: sorbent_vvel(data%num_sorp)
!
!-------------------------------------------------------------------------------
!BEGIN
   temp = _STATE_VAR_(data%id_tem)

   ! Set velocity of sorped pesticides, if simulated.
   DO pst_i=1,data%num_pesticides

     DO sorp_i=1,data%num_sorp
       sorbent_vvel(sorp_i) = 0. ! INSERT LINKS HERE
       mobility(data%id_psta(pst_i,sorp_i)) = sorbent_vvel(sorp_i)
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
   DO ss_i=1,data%num_ss
      ! Retrieve current (local) state variable values.
      ss = _STATE_VAR_(data%id_ss(ss_i))
      ! Self-shading with contribution from this phytoplankton concentration.
      extinction = extinction + (data%Ke_SS(ss_i)*ss)
   ENDDO
END SUBROUTINE aed_light_extinction_pesticides
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
FUNCTION pesticide_sorption(temp,sorbents,pest_t,pest_s,kkkk) RESULT(pest_d)
!-------------------------------------------------------------------------------
! Function to partition pesticide concentration amongst several sorbents
!-------------------------------------------------------------------------------
!ARGUMENTS
!   INTEGER,INTENT(in)  :: method
   AED_REAL,INTENT(in) :: temp,sorbents(:),pest_t,kkkk
   AED_REAL,INTENT(inout) :: pest_s(:)
   AED_REAL :: pest_d
!
!LOCALS
   AED_REAL  :: fT        !-- Value of the temperature function
   AED_REAL,PARAMETER  :: tp = 20.0
!
!-------------------------------------------------------------------------------
!BEGIN

   pest_s(:) = zero_
   pest_d = pest_t

END FUNCTION pesticide_sorption
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_pesticides
