!###############################################################################
!#                                                                             #
!# aed_phosphorus.F90                                                          #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2013 - 2025 -  The University of Western Australia               #
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
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Created 24 August 2011                                                      #
!#                                                                             #
!###############################################################################
!                                                                              !
!          .----------------.  .----------------.  .----------------.          !
!          | .--------------. || .--------------. || .--------------. |        !
!          | |   ______     | || |  ____  ____  | || |    _______   | |        !
!          | |  |_   __ \   | || | |_   ||   _| | || |   /  ___  |  | |        !
!          | |    | |__) |  | || |   | |__| |   | || |  |  (__ \_|  | |        !
!          | |    |  ___/   | || |   |  __  |   | || |   '.___`-.   | |        !
!          | |   _| |_      | || |  _| |  | |_  | || |  |`\____) |  | |        !
!          | |  |_____|     | || | |____||____| | || |  |_______.'  | |        !
!          | |              | || |              | || |              | |        !
!          | '--------------' || '--------------' || '--------------' |        !
!          '----------------'  '----------------'  '----------------'          !
!                                                                              !
!###############################################################################

#include "aed.h"


MODULE aed_phosphorus
!------------------------------------------------------------------------------+
! aed_phosphorus --- phosphorus biogeochemical model
!
! The AED module phosphorus contains equations that describe exchange of
! soluable reactive phosphorus across the air/water interface, sediment flux
! and sorption/desorption
!------------------------------------------------------------------------------+
   USE aed_core
   USE aed_util, ONLY: PO4AdsorptionFraction

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed_phosphorus_data_t
!
   TYPE,extends(aed_model_data_t) :: aed_phosphorus_data_t
      !# Variable identifiers
      INTEGER  :: id_frp, id_frpads, id_oxy,  id_tss, id_pH
      INTEGER  :: id_Fsed_frp, id_l_resus
      INTEGER  :: id_e_temp, id_e_salt, id_e_rain, id_tssext, id_dz
      INTEGER  :: id_sed_frp, id_frpads_vvel, id_atm_dep, &
                  id_frpads_set, id_frp_srp, id_frpads_res, id_frpads_swi

      !# Model parameters
      AED_REAL :: Fsed_frp,Ksed_frp,theta_sed_frp             ! Benthic
      AED_REAL :: atm_pip_dd, atm_frp_conc                    ! Deposition
      AED_REAL :: Kpo4p, Kadsratio, Qmax, theta_Kpo4, K_sal   ! Adsoprtion
      AED_REAL :: w_po4ads                                    ! Sedimentation
      LOGICAL  :: simDryDeposition,simWetDeposition
      LOGICAL  :: ben_use_oxy,ben_use_aedsed
      INTEGER  :: PO4AdsorptionModel, resuspension
      LOGICAL  :: simPO4Adsorption, ads_use_pH, ads_use_external_tss

     CONTAINS
         PROCEDURE :: define            => aed_define_phosphorus
         PROCEDURE :: calculate_benthic => aed_calculate_benthic_phosphorus
         PROCEDURE :: calculate_surface => aed_calculate_surface_phosphorus
         PROCEDURE :: equilibrate       => aed_equilibrate_phosphorus
         PROCEDURE :: mobility          => aed_mobility_phosphorus
!        PROCEDURE :: light_extinction  => aed_light_extinction_phosphorus
!        PROCEDURE :: delete            => aed_delete_phosphorus

   END TYPE

! MODULE GLOBALS
   INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
                                              ! 1 = basic diagnostic outputs
                                              ! 2 = flux rates, and supporitng
                                              ! 3 = other metrics
                                              !10 = all debug & checking outputs

!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed_define_phosphorus(data, namlst)
!------------------------------------------------------------------------------+
! Initialise the AED model
!
!  Here, the aed namelist is read and the variables to be simulated
!  by the model are registered with AED.
!------------------------------------------------------------------------------+
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed_phosphorus_data_t),INTENT(inout) :: data

!
!LOCALS
   INTEGER  :: status

!  %% NAMELIST   %%  /aed_phosphorus/
!  %% Last Checked 20/08/2021
   ! Initial
   AED_REAL          :: frp_initial   = 4.5
   AED_REAL          :: frp_min       = zero_
   AED_REAL          :: frp_max       = nan_
   ! Benthic
   AED_REAL          :: Fsed_frp      = 3.5
   AED_REAL          :: Ksed_frp      = 30.0
   AED_REAL          :: theta_sed_frp = 1.05
   CHARACTER(len=64) :: phosphorus_reactant_variable=''
   CHARACTER(len=64) :: Fsed_frp_variable=''
   ! Adsorption
   LOGICAL           :: simPO4Adsorption = .FALSE.
   LOGICAL           :: ads_use_external_tss = .FALSE.
   INTEGER           :: PO4AdsorptionModel = 1
   LOGICAL           :: ads_use_pH   = .FALSE.
   AED_REAL          :: Kpo4p        = 1.05
   AED_REAL          :: theta_Kpo4   = 1.02
   AED_REAL          :: K_sal        = 1000
   AED_REAL          :: Kadsratio    = 1.05
   AED_REAL          :: Qmax         = 1.05
   AED_REAL          :: w_po4ads     = zero_
   CHARACTER(len=64) :: po4sorption_target_variable=''
   CHARACTER(len=64) :: frp_ads_particle_link=''            !   For FV API 2.0 (To be implemented)
   CHARACTER(len=64) :: pH_variable  ='CAR_pH'
   CHARACTER(len=64) :: resus_link   =''
   INTEGER           :: resuspension = 0
   ! Atmospheric deposition
   LOGICAL           :: simDryDeposition = .false.
   LOGICAL           :: simWetDeposition = .false.
   AED_REAL          :: atm_pip_dd   = zero_
   AED_REAL          :: atm_frp_conc = zero_
! %% From Module Globals
!  INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
!                                             ! 1 = basic diagnostic outputs
!                                             ! 2 = flux rates, and supporitng
!                                             ! 3 = other metrics
!                                             !10 = all debug & checking outputs
!  %% END NAMELIST   %%  /aed_phosphorus/

   NAMELIST /aed_phosphorus/ frp_initial,frp_min,frp_max,                      &
                            Fsed_frp,Ksed_frp,theta_sed_frp,Fsed_frp_variable, &
                            phosphorus_reactant_variable,                      &
                            simPO4Adsorption,PO4AdsorptionModel,               &
                            ads_use_external_tss,frp_ads_particle_link,        &
                            po4sorption_target_variable,                       &
                            Kpo4p,theta_Kpo4,K_sal,Kadsratio,Qmax,             &
                            ads_use_pH,pH_variable, w_po4ads,                  &
                            simDryDeposition, simWetDeposition,                &
                            atm_pip_dd, atm_frp_conc,                          &
                            resuspension, resus_link, diag_level
!
!------------------------------------------------------------------------------+
!BEGIN
   print *,"        aed_phosphorus configuration"

   ! Read the namelist
   read(namlst,nml=aed_phosphorus,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist for &aed_phosphorus'

   ! Store parameter values in the module level data object
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   data%Fsed_frp             = Fsed_frp/secs_per_day
   data%Ksed_frp             = Ksed_frp
   data%theta_sed_frp        = theta_sed_frp
   data%simPO4Adsorption     = simPO4Adsorption
   data%ads_use_external_tss = ads_use_external_tss
   data%PO4AdsorptionModel   = PO4AdsorptionModel
   data%ads_use_pH           = ads_use_pH
   data%Kpo4p                = Kpo4p
   data%theta_Kpo4           = theta_Kpo4
   data%K_sal                = K_sal
   data%Kadsratio            = Kadsratio
   data%Qmax                 = Qmax
   data%w_po4ads             = w_po4ads/secs_per_day
   data%atm_pip_dd           = MAX(zero_,atm_pip_dd/secs_per_day)
   data%atm_frp_conc         = MAX(zero_,atm_frp_conc)
   data%simDryDeposition     = simDryDeposition
   data%simWetDeposition     = simWetDeposition
   data%resuspension         = resuspension

   ! Register main state variable
   data%id_frp = aed_define_variable( 'frp', 'mmol P/m3',                      &
                                    'filterable reactive phosphate',           &
                                    frp_initial,minimum=frp_min,maximum=frp_max)

   ! Register external state variable dependencies (for benthic flux)
   data%ben_use_oxy = phosphorus_reactant_variable .NE. '' !This means oxygen module switched on
   IF (data%ben_use_oxy) &
     data%id_oxy = aed_locate_variable(phosphorus_reactant_variable)

   data%ben_use_aedsed = Fsed_frp_variable .NE. '' !This means aed sediment module switched on
   IF (data%ben_use_aedsed) &
     data%id_Fsed_frp = aed_locate_sheet_variable(Fsed_frp_variable)

   data%id_frpads = -1
   data%id_frpads_vvel = -1
   ! Check if particles and PO4 adsorption are simulated
   IF (data%simPO4Adsorption) THEN
     IF (data%ads_use_external_tss) THEN
         PRINT *,'        PO4 adsorption is configured to use external TSS var'
         data%id_tssext = aed_locate_global('tss')
     ELSE
       IF (po4sorption_target_variable .NE. '' ) THEN
          print *,'          PO4 is adsorbing to ',TRIM(po4sorption_target_variable)
          print *,'          ... found'
          data%id_tss = aed_locate_variable(po4sorption_target_variable)
          IF(w_po4ads<-999.) THEN
            print *,'          Checking for associated _vvel link array ',TRIM(po4sorption_target_variable)//'_vvel'
            data%id_frpads_vvel = aed_locate_variable(TRIM(po4sorption_target_variable)//'_vvel')
            print *,'          ... found'
            data%w_po4ads = zero_
          ELSE
            PRINT *,'  ERROR PO4 adsorption vvel link variable not found even though w_po4ads specifies link'
          ENDIF
       ELSE
          PRINT *,'  ERROR PO4 adsorption is configured but no internal or external target variable is set'
          STOP
       ENDIF
     ENDIF

     data%id_frpads = aed_define_variable('frp_ads','mmol P/m3','adsorbed phosphate', &
                      zero_,minimum=zero_,mobility=data%w_po4ads)

     IF (data%ads_use_pH) THEN
       data%id_pH = aed_locate_variable(pH_variable)
     ENDIF

     !-- resuspension link variable
     IF ( resuspension>0 .AND. .NOT.resus_link .EQ. '' ) THEN
        data%id_l_resus  = aed_locate_sheet_variable(TRIM(resus_link)) ! ('NCS_resus')
       !IF ( resuspension>1 ) THEN
       !    data%id_sedpipfr  = aed_define_sheet_diag_variable('pip_sed_frac','w/w',  'Sediment PIP fraction')
       !  ENDIF
     ELSE
        data%id_l_resus = 0
        data%resuspension = 0
     ENDIF

     !-- Check diagnostics specific for adsorbed phosphate
     data%id_frpads_set = aed_define_diag_variable('frp_ads_set','mmol P/m3/d',&
                                           'adsobed PO4 sedimentation flux')
     data%id_frpads_res = aed_define_sheet_diag_variable('frp_ads_res','mmol P/m2/d',&
                                           'adsobed PO4 resuspension flux')
     data%id_frpads_swi = aed_define_sheet_diag_variable('frp_ads_swi','mmol P/m2/d',&
                                           'adsobed PO4 net flux at the swi')
     data%id_frp_srp = aed_define_diag_variable('frp_srp','mmol P/m3/d',       &
                                           'PO4 adsorption rate')

   ENDIF

   ! Register diagnostic variables
   data%id_sed_frp = aed_define_sheet_diag_variable('frp_dsf','mmol P/m2/d', &
                                         'PO4 exchange across sed/water interface')
   IF( simWetDeposition .OR. simDryDeposition ) THEN
    data%id_atm_dep = aed_define_sheet_diag_variable('dip_atm','mmol P/m2/d', &
                                         'DIP atmospheric deposition flux')
   ENDIF

   ! Register environmental dependencies
   data%id_e_temp = aed_locate_global('temperature')
   data%id_e_salt = aed_locate_global('salinity')
   data%id_dz     = aed_locate_global('layer_ht')
   IF( simWetDeposition ) data%id_E_rain = aed_locate_sheet_global('rain')

END SUBROUTINE aed_define_phosphorus
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_equilibrate_phosphorus(data,column,layer_idx)
!------------------------------------------------------------------------------+
! Update partitioning of phosphate between dissolved and particulate pools
! after kinetic transformations are applied
!------------------------------------------------------------------------------+
!ARGUMENTS
   CLASS (aed_phosphorus_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp, tss, salt

   ! State
   AED_REAL :: frp,frpads,pH

   ! Temporary variables
   AED_REAL :: PO4dis, PO4par, PO4tot, Kpo4p

!-------------------------------------------------------------------------------
!BEGIN

   IF(.NOT. data%simPO4Adsorption) RETURN

   tss = zero_

   ! Retrieve current environmental conditions for the cell.
   temp = _STATE_VAR_(data%id_E_temp)    ! local temperature
   salt = _STATE_VAR_(data%id_E_salt)    ! local salinity
   IF(data%ads_use_external_tss) THEN
     tss = _STATE_VAR_(data%id_tssext) ! externally supplied total susp solids
   END IF

   ! Retrieve current (local) state variable values.
   frp = _STATE_VAR_(data%id_frp)            ! dissolved PO4
   frpads = _STATE_VAR_(data%id_frpads)      ! adsorped PO4
   IF (.NOT.data%ads_use_external_tss) &
      tss = _STATE_VAR_(data%id_tss)         ! local total susp solids


   ! Adjust local sorption coefficients for temperature or salinity  (PO4AdsorptionModel = 1 only)
   Kpo4p = data%Kpo4p * Kpo4p_fT_fSal(data%theta_Kpo4, data%K_sal, salt, temp)

   ! Compute sorption
   IF(data%ads_use_pH) THEN
     pH = _STATE_VAR_(data%id_pH)

     CALL PO4AdsorptionFraction(data%PO4AdsorptionModel,              &  ! Dependencies
                                 frp+frpads,                          &
                                 tss,                                 &
                                 Kpo4p,data%Kadsratio,data%Qmax,      &
                                 PO4dis,PO4par,                       &  ! Returning variables
                                 thepH=pH)

   ELSE
     CALL PO4AdsorptionFraction(data%PO4AdsorptionModel,              &  ! Dependecies
                                 frp+frpads,                          &
                                 tss,                                 &
                                 Kpo4p,data%Kadsratio,data%Qmax,      &
                                 PO4dis,PO4par,                       &  ! Returning variables
                                 temp_=temp,salt_=salt)
   ENDIF

   ! Set back to core variables
   _STATE_VAR_(data%id_frp)    = PO4dis    ! Dissolved PO4 (FRP)
   _STATE_VAR_(data%id_frpads) = PO4par    ! Adsorped PO4  (PIP)

END SUBROUTINE aed_equilibrate_phosphorus
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_surface_phosphorus(data,column,layer_idx)
!------------------------------------------------------------------------------+
! Air-water exchange for the aed phosphorus model. Includes wet/dry deposition,
! depending on the configuration.
!------------------------------------------------------------------------------+
!ARGUMENTS
   CLASS (aed_phosphorus_data_t),INTENT(in) :: data
   TYPE  (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp, salt, wind, vel, depth, rain
   ! State
   AED_REAL :: n2o, nox, no2, nh4
   ! Temporary variables
   AED_REAL :: n2o_atm_flux = zero_ ! N2O atmos exchange with water
   AED_REAL :: Cn2o_air = zero_     ! N2O in the air phase
   AED_REAL :: kn2o_trans = zero_   ! N2O piston velocity
   AED_REAL :: wind_hgt             ! Height of wind data, for correction
   AED_REAL :: f_pres = 1.0         ! Pressure correction function (currently
                                    !  unsued but applicable at high altitudes)
!
!------------------------------------------------------------------------------+
!BEGIN

  !----------------------------------------------------------------------------+
  !# Atmosphere loading of DIP to the water, due to dry or wet deposition
  IF( data%simDryDeposition ) THEN
    !-----------------------------------------------
    ! Set surface exchange value (mmmol/m2/s) for AED ODE solution.
   IF (data%simPO4Adsorption) & !# id_frpads is not set unless simPO4Adsorption is true
    _FLUX_VAR_T_(data%id_frpads) = data%atm_pip_dd
  ENDIF

  IF( data%simWetDeposition ) THEN
    !-----------------------------------------------
    ! Get the necessary environmental variables (from physical driver)
    rain = _STATE_VAR_S_(data%id_E_rain) / secs_per_day   ! Rain (m/s)

    !-----------------------------------------------
    ! Set surface exchange value (mmmol/m2/s) for AED ODE solution.
    _FLUX_VAR_T_(data%id_frp) = _FLUX_VAR_T_(data%id_frp) &
                              + rain * data%atm_frp_conc
  ENDIF

  IF( data%simDryDeposition .OR. data%simWetDeposition ) THEN
    !-----------------------------------------------
    ! Also store deposition across the atm/water interface as a
    ! diagnostic variable (mmmol/m2/day).
   IF (data%simPO4Adsorption) THEN !# id_frpads is not set unless simPO4Adsorption is true
        _DIAG_VAR_S_(data%id_atm_dep) = (_FLUX_VAR_T_(data%id_frp) + _FLUX_VAR_T_(data%id_frpads)) * secs_per_day
   ELSE
        _DIAG_VAR_S_(data%id_atm_dep) = _FLUX_VAR_T_(data%id_frp) * secs_per_day
   ENDIF
  ENDIF

END SUBROUTINE aed_calculate_surface_phosphorus
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_benthic_phosphorus(data,column,layer_idx)
!------------------------------------------------------------------------------+
! Calculate the bottom fluxes and benthic sink & source terms of AED phosphorus.
! Everything in units per surface area (not volume!) per time.
!------------------------------------------------------------------------------+
!ARGUMENTS
   CLASS (aed_phosphorus_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp,dz

   ! State
   AED_REAL :: frp,oxy

   ! Temporary variables
   AED_REAL :: frp_flux, Fsed_frp, fT, fDO, sedpipfr, Fsed_pip

!-------------------------------------------------------------------------------
!BEGIN

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   temp = _STATE_VAR_(data%id_E_temp)       ! local temperature
   dz = _STATE_VAR_(data%id_dz)

   ! Retrieve current (local) state variable values.
   frp = _STATE_VAR_(data%id_frp)           ! phosphorus

   ! Compute the sediment flux dependent on overlying oxygen & temperature
   fT = data%theta_sed_frp**(temp-20.0)
   IF (data%ben_use_oxy) THEN
      ! Sediment flux dependent on oxygen and temperature
      oxy = _STATE_VAR_(data%id_oxy)
      fDO = data%Ksed_frp/(data%Ksed_frp+oxy)
   ELSE
      ! Sediment flux dependent on temperature only
      fDO = one_
   ENDIF

   ! Set the flux dependent on if and how it is linked
   IF (data%ben_use_aedsed) THEN
     ! Linked to aed_sedflux, check if its constant or dynamically set
     IF ( aed_is_const_var(data%id_Fsed_frp) ) THEN
        Fsed_frp = _DIAG_VAR_S_(data%id_Fsed_frp) * MIN(3.,fDO * fT) / secs_per_day
     ELSE
        Fsed_frp = _DIAG_VAR_S_(data%id_Fsed_frp) / secs_per_day
     ENDIF
   ELSE
     Fsed_frp = data%Fsed_frp * MIN(3.,fDO * fT)
   ENDIF
   !IF (data%ben_use_aedsed) THEN
   !    Fsed_frp = _STATE_VAR_S_(data%id_Fsed_frp)
   !ELSE
   !    Fsed_frp = data%Fsed_frp
   !ENDIF
   frp_flux = Fsed_frp

   ! Set bottom fluxes for the pelagic (change per surface area per second)
   _FLUX_VAR_(data%id_frp) = _FLUX_VAR_(data%id_frp) + frp_flux

   ! Set sink and source terms for the benthos (change per surface area per second)
   ! Note that this should include the fluxes to and from the pelagic.
   !_FLUX_VAR_B_(data%id_ben_frp) = _FLUX_VAR_B_(data%id_ben_frp) + (-frp_flux)

   ! Also store sediment flux as diagnostic variable.
   _DIAG_VAR_S_(data%id_sed_frp) = frp_flux * secs_per_day

   IF (data%simPO4Adsorption) THEN
     Fsed_pip = zero_
     IF( data%resuspension > 0 ) THEN
       sedpipfr = 100./1e6 ! 100mg/kg
       Fsed_pip = _DIAG_VAR_S_(data%id_l_resus) * sedpipfr * 1e3/30.91   ! g/m2/s * gP/gSed * mmol/mol / g/mol = mmol/m2/s
       _DIAG_VAR_S_(data%id_frpads_res) = Fsed_pip
     ENDIF
     _DIAG_VAR_S_(data%id_frpads_swi) = _DIAG_VAR_(data%id_frpads_set)*dz + Fsed_pip*secs_per_day
   ENDIF


END SUBROUTINE aed_calculate_benthic_phosphorus
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_mobility_phosphorus(data,column,layer_idx,mobility)
!-------------------------------------------------------------------------------
! Get the vertical movement velocities of frp_ads (+ve up; -ve down)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_phosphorus_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: mobility(:)
!
!LOCALS
   AED_REAL :: vvel, dz
!-------------------------------------------------------------------------------
!BEGIN

   IF(.NOT. data%simPO4Adsorption) RETURN

   vvel = zero_
   dz = _STATE_VAR_(data%id_dz)

   IF( data%id_frpads_vvel>0 ) THEN
     ! adopt vertical velocity of host particle
     vvel = _DIAG_VAR_(data%id_frpads_vvel) / secs_per_day
   ELSE
     ! adopt constant value read in from nml
     vvel = data%w_po4ads
   ENDIF
   mobility(data%id_frpads) = vvel

   _DIAG_VAR_(data%id_frpads_set) = (vvel/dz)*_STATE_VAR_(data%id_frpads)*secs_per_day

END SUBROUTINE aed_mobility_phosphorus
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
PURE AED_REAL FUNCTION Kpo4p_fT_fSal(theta_Kpo4, K_sal, sal, temp)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for P-adsorption sensitivity to Sal & T
! Ref: Zhang and Huang 2011, in Floriday bay, suggesting increasing P-adsorption
! capability with increasing temperature and decreasing salinity
!-------------------------------------------------------------------------------
!ARGUMENTS
AED_REAL,INTENT(in) :: theta_Kpo4 ! theta for P-adsorp on limitation, default = 1.02
AED_REAL,INTENT(in) :: K_sal ! half-saturation of salinity, default = 60
AED_REAL,INTENT(in) :: temp
AED_REAL,INTENT(in) :: sal
AED_REAL,PARAMETER  :: Topt = 45. ! optimum temperature for P-adsorp, default = 45 as the
                                  !       ref showed increasing P-adsorp with temperature

AED_REAL :: fT, fSal
!
!-------------------------------------------------------------------------------
!BEGIN

fT = (theta_Kpo4**(temp-Topt))

IF(K_sal==zero_)THEN
fSal = one_
ELSE
fSal = K_sal/(K_sal+sal)
ENDIF

Kpo4p_fT_fSal = fT * fSal

END FUNCTION Kpo4p_fT_fSal
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_phosphorus
