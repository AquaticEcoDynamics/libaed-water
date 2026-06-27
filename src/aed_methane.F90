!###############################################################################
!#                                                                             #
!# aed_methane.F90                                                              #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2013 - 2022 -  The University of Western Australia               #
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
!# Created March 2012                                                          #
!#  Track changes on GitHub @ https://github.com/AquaticEcoDynamics/libaed-water
!#                                                                             #
!#                                                                             #
!###############################################################################
!                                                                              !
!         .----------------.  .----------------.  .----------------.           !
!         | .--------------. || .--------------. || .--------------. |         !
!         | |     ______   | || |      __      | || |  _______     | |         !
!         | |   .' ___  |  | || |     /  \     | || | |_   __ \    | |         !
!         | |  / .'   \_|  | || |    / /\ \    | || |   | |__) |   | |         !
!         | |  | |         | || |   / ____ \   | || |   |  __ /    | |         !
!         | |  \ `.___.'\  | || | _/ /    \ \_ | || |  _| |  \ \_  | |         !
!         | |   `._____.'  | || ||____|  |____|| || | |____| |___| | |         !
!         | |              | || |              | || |              | |         !
!         | '--------------' || '--------------' || '--------------' |         !
!         '----------------'  '----------------'  '----------------'           !
!                                                                              !
!###############################################################################

#include "aed.h"

!
MODULE aed_methane
!-------------------------------------------------------------------------------
! aed_methane --- Methane model
!
! The AED methane module contains equations that describe the dissolved and bubble 
! methane dynamics in the sediments, water column and water-air interface. 
!-------------------------------------------------------------------------------
   USE aed_core

   USE aed_util,  ONLY: aed_gas_piston_velocity

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed_methane_data_t

   TYPE,extends(aed_model_data_t) :: aed_methane_data_t
!
      !# Variable identifiers
      INTEGER  :: id_ch4, id_oxy, id_ch4_bub
      INTEGER  :: id_Fsed_ch4, id_Fsed_ch4_ebb
      INTEGER  :: id_temp, id_salt
      INTEGER  :: id_wind, id_vel, id_depth, id_depth_c, id_dens, id_zone
      INTEGER  :: id_par, id_extc, id_dz, id_tau, id_air_pres, id_layer_ht, id_area, id_volume, id_q_cub, id_u_star, id_Q_net
      INTEGER  :: id_ch4ox
      INTEGER  :: id_sed_ch4, id_sed_ch4_ebb, id_sed_ch4_ebb_3d
      INTEGER  :: id_atm_ch4, id_atm_ch4_ebb, id_ch4_ebb_df, id_kCH4, id_beta, id_epsilon
      INTEGER  :: id_ch4_s, id_poc_s, id_ch4s_prod, id_ch4s_oxid, id_ch4s_diff, id_ch4s_ebb, id_ch4_crit, id_poc_resp, id_poc_set, id_phy_set
      INTEGER  :: id_poc_swi, id_phy_swi, id_ch4_oxid_target, id_gpp, id_ch4_oxic_prod, id_ch4_pw
      INTEGER  :: id_tempzon
      INTEGER  :: id_ch4_bub_s, id_ice_w, id_ice_b
      

      !# Model parameters
      AED_REAL :: Fsed_ch4, Ksed_ch4, theta_sed_ch4, Fsed_ch4_ebb
      AED_REAL :: Rch4ox, Kch4ox, vTch4ox, Kch4ox_ch4, beta_ch4, alpha_oxy
      AED_REAL :: atm_ch4
      AED_REAL :: ch4_bub_aLL, ch4_bub_cLL, ch4_bub_kLL
      AED_REAL :: ch4_bub_disf1, ch4_bub_disf2, ch4_bub_disdp, ch4_oxic_prod
      AED_REAL :: poc_s_initial, poc_s_unr, ch4_bub_s_initial
      AED_REAL :: Rch4ox_sed, Kch4ox_sed, vTch4ox_sed
      AED_REAL :: Rpoc_resp, theta_poc_resp, tau_poc_burial
      AED_REAL :: Kpocox_sed
      AED_REAL :: ch4_s_diff_coef, z_dbl, Phi, sed_z, Ksed_beta, oxy_sed_p, gamma_max, fredox
      AED_REAL :: alpha, LA, Kch4ox_diff
      


      !# Model options
      LOGICAL  :: use_oxy, use_sed_model_ch4, use_sed_model_ebb
      LOGICAL  :: simCH4, simCH4ebb
      INTEGER  :: ebb_model, MOx_model
      INTEGER  :: ch4_piston_model, ch4_schmidt_model

     CONTAINS
         PROCEDURE :: define            => aed_define_methane
         PROCEDURE :: calculate_surface => aed_calculate_surface_methane
         PROCEDURE :: calculate         => aed_calculate_methane
         PROCEDURE :: calculate_benthic => aed_calculate_benthic_methane
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
SUBROUTINE aed_define_methane(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and the variables exported
!  by the model are registered with AED.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_methane_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst
!
!LOCALS
   INTEGER           :: status

!  %% NAMELIST   %%  /aed_methane/
!  %% Last Checked 26/07/2024
   AED_REAL          :: ch4_initial      = 4.5
   LOGICAL           :: simCH4
   AED_REAL          :: Fsed_ch4         = 0.0
   AED_REAL          :: Ksed_ch4         = 30.0
   AED_REAL          :: theta_sed_ch4    = 1.0
   CHARACTER(len=64) :: Fsed_ch4_variable=''
   AED_REAL          :: atm_ch4          = 1.76e-6
   AED_REAL          :: Rch4ox           = 0.01
   AED_REAL          :: Kch4ox           = 0.01
   AED_REAL          :: KCH4ox_ch4       = 20.0
   AED_REAL          :: vTch4ox          = 1.05
   AED_REAL          :: beta_ch4         = 1.2
   AED_REAL          :: alpha_oxy        = 10.0
   CHARACTER(len=64) :: methane_reactant_variable=''
   CHARACTER(len=64)  :: ch4_oxid_target_variable=''

   INTEGER           :: ebb_model        = 0
   INTEGER           :: MOx_model        = 1
   INTEGER           :: ch4_piston_model = 1
   INTEGER           :: ch4_schmidt_model = 5 !freshwater

   LOGICAL           :: simCH4ebb
   AED_REAL          :: Fsed_ch4_ebb     = zero_
   CHARACTER(len=64) :: Fsed_ebb_variable=''
   AED_REAL          :: ch4_bub_aLL      = 42.9512677
   AED_REAL          :: ch4_bub_cLL      = 0.634
   AED_REAL          :: ch4_bub_kLL      = -0.8247
   AED_REAL          :: ch4_bub_disf1    = 0.07
   AED_REAL          :: ch4_bub_disf2    = 0.33
   AED_REAL          :: ch4_bub_disdp    = 20.0
   AED_REAL          :: ch4_s_initial    = 100.0
   AED_REAL          :: ch4_bub_s_initial = 0.0 !Initial methane bubble storage in the sediment
   AED_REAL          :: poc_s_initial    = 600000 !Initial poc in the sediment
   AED_REAL          :: poc_s_unr        = 500000 !Initial unreactive POC in the sediment
   AED_REAL          :: Rch4ox_sed       = 0.01
   AED_REAL          :: Kch4ox_sed       = 0.02 !this is mmol/m2/second
   AED_REAL          :: vTch4ox_sed      = 1.05
   AED_REAL          :: Rpoc_resp        = 0.00011
   AED_REAL          :: tau_poc_burial   = 1.5e8
   AED_REAL          :: oxy_sed_p        = 0.03
   AED_REAL          :: fredox           = 0.3
   AED_REAL          :: gamma_max        = 0.7
   AED_REAL          :: theta_poc_resp   = 1.07
   AED_REAL          :: Kpocox_sed       = 30
   AED_REAL          :: ch4_s_diff_coef  = 0.000001
   AED_REAL          :: z_dbl            = 0.1
   AED_REAL          :: Phi              = 0.7
   AED_REAL          :: sed_z            = 0.5
   AED_REAL          :: alpha            = 0.5
   AED_REAL          :: LA               = 0.0
   AED_REAL          :: Ksed_beta        = 0.08
   AED_REAL          :: ch4_oxic_prod    = 0.05
   AED_REAL          :: Kch4ox_diff      = 30
   !CHARACTER(len=64) :: Fsed_poc_variable=''
   

! %% From Module Globals
!  INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
!                                             ! 1 = basic diagnostic outputs
!                                             ! 2 = flux rates, and supporitng
!                                             ! 3 = other metrics
!                                             !10 = all debug & checking outputs

!  %% END NAMELIST   %%  /aed_carbon/

   NAMELIST /aed_methane/ch4_initial,                                                &
                         Rch4ox,Kch4ox,vTch4ox,Kch4ox_ch4,                           &
                         beta_ch4, alpha_oxy,                                        &
                         methane_reactant_variable, simCH4,                          &
                         Fsed_ch4,Ksed_ch4,theta_sed_ch4,Fsed_ch4_variable,          &
                         atm_ch4, ch4_oxid_target_variable,                          &
                         ch4_piston_model, ch4_schmidt_model,                        &
                         ebb_model,MOx_model,                                        &
                         Fsed_ch4_ebb, Fsed_ebb_variable,                            &
                         ch4_bub_aLL,ch4_bub_cLL, ch4_bub_kLL,                       &
                         ch4_bub_disf1, ch4_bub_disf2, ch4_bub_disdp,                &
                         ch4_s_initial, poc_s_initial, poc_s_unr, ch4_bub_s_initial, &
                         Rch4ox_sed, Kch4ox_sed, vTch4ox_sed,                        &
                         Rpoc_resp, theta_poc_resp, tau_poc_burial,                  &
                         Kpocox_sed, gamma_max, oxy_sed_p,                           &
                         ch4_s_diff_coef, z_dbl, Phi, sed_z,                         & 
                         LA, alpha, Ksed_beta,fredox,                                &
                         ch4_oxic_prod, Kch4ox_diff,                                 &
                         diag_level


!-------------------------------------------------------------------------------
!BEGIN

   print *,"        aed_methane configuration"

   !# Set defaults
   data%simCH4        = .true.
   data%simCH4ebb     = .false.

   !# Read the namelist
   read(namlst,nml=aed_methane,iostat=status)
   IF (status /= 0) THEN
      print *,'Error reading namelist for &aed_methane'
      STOP
   ENDIF

   IF (ebb_model>0) THEN 
      data%simCH4ebb = .true.
   ELSE
      data%simCH4ebb = .false.
   ENDIF

   !# Store parameter values in modules own derived type
   !  Note: rates are provided in values per day, and
   !        are converted here to values per second.
   data%simCH4           = simCH4
   data%Fsed_ch4         = Fsed_ch4/secs_per_day
   data%Ksed_ch4         = Ksed_ch4
   data%theta_sed_ch4    = theta_sed_ch4
   data%Rch4ox           = Rch4ox/secs_per_day
   data%Kch4ox           = Kch4ox
   data%Kch4ox_ch4       = Kch4ox_ch4
   data%vTch4ox          = vTch4ox
   data%atm_ch4          = atm_ch4
   data%ch4_piston_model = ch4_piston_model
   data%ch4_schmidt_model = ch4_schmidt_model
   data%ebb_model        = ebb_model
   data%MOx_model        = MOx_model
   data%Fsed_ch4_ebb     = Fsed_ch4_ebb
   data%ch4_bub_aLL      = ch4_bub_aLL
   data%ch4_bub_cLL      = ch4_bub_cLL
   data%ch4_bub_kLL      = ch4_bub_kLL
   data%ch4_bub_disf1    = ch4_bub_disf1
   data%ch4_bub_disf2    = ch4_bub_disf2
   data%ch4_bub_disdp    = ch4_bub_disdp
   data%Rch4ox_sed       = Rch4ox_sed /secs_per_day
   data%Kch4ox_sed       = Kch4ox_sed
   data%vTch4ox_sed      = vTch4ox_sed
   data%Rpoc_resp        = Rpoc_resp/secs_per_day
   data%tau_poc_burial   = tau_poc_burial
   data%oxy_sed_p        = oxy_sed_p
   data%fredox           = fredox
   data%gamma_max        = gamma_max
   data%theta_poc_resp   = theta_poc_resp
   data%Kpocox_sed       = Kpocox_sed
   data%ch4_s_diff_coef  = ch4_s_diff_coef
   data%poc_s_initial    = poc_s_initial
   data%poc_s_unr        = poc_s_unr
   data%ch4_bub_s_initial = ch4_bub_s_initial
   data%z_dbl            = z_dbl
   data%Phi              = Phi
   data%sed_z            = sed_z
   data%alpha            = alpha
   data%LA               = LA
   data%Ksed_beta        = Ksed_beta
   data%ch4_oxic_prod    = ch4_oxic_prod
   data%Kch4ox_diff      = Kch4ox_diff
   data%beta_ch4         = beta_ch4
   data%alpha_oxy        = alpha_oxy
  

   !# Register state variables
      data%id_ch4 = aed_define_variable('ch4','mmol C/m3','methane',    &
                                     ch4_initial,minimum=zero_)
      IF( data%simCH4ebb ) THEN
           data%id_ch4_bub = aed_define_variable('ch4_bub','mmol C/m3', &
                                     'methane bubbles',zero_,minimum=zero_)
          IF( data%ebb_model==2 ) THEN
           data%id_ch4_s = aed_define_sheet_variable('ch4_sed', 'mmol C/m2', &
                                     'sediment methane concentration', ch4_s_initial, minimum=zero_)
           data%id_poc_s = aed_define_sheet_variable('poc_sed', 'mmol C/m2', &
                                     'sediment POC concentration', poc_s_initial - poc_s_unr, minimum=zero_)
           data%id_ch4_bub_s = aed_define_sheet_variable('ch4_bub_s', 'mmol C/m2', &
                                     'sediment methane bubble storage', ch4_bub_s_initial, minimum=zero_)
          ENDIF
      ENDIF

   !# Register external state variable dependencies
   data%use_oxy = methane_reactant_variable .NE. '' !i.e., oxygen module engaged
   IF (data%use_oxy) THEN
      data%id_oxy = aed_locate_variable(methane_reactant_variable)
   ENDIF

   data%use_sed_model_ch4 = Fsed_ch4_variable .NE. ''
   IF (data%use_sed_model_ch4) &
      data%id_Fsed_ch4 = aed_locate_sheet_variable(Fsed_ch4_variable)

   data%use_sed_model_ebb = Fsed_ebb_variable .NE. ''
   IF (data%use_sed_model_ebb) &
      data%id_Fsed_ch4_ebb = aed_locate_sheet_variable(Fsed_ebb_variable)

   !IF (data%ebb_model == 2) data%id_poc_swi = aed_locate_sheet_variable("OGM_poc_swi")
   IF (data%ebb_model == 2) data%id_poc_swi = aed_locate_variable("OGM_poc_set")
   IF (data%ebb_model == 2) data%id_phy_swi = aed_locate_sheet_variable("PHY_phy_swi_c")

   data%id_ch4_oxid_target = aed_locate_variable(ch4_oxid_target_variable)

   data%id_gpp = aed_locate_variable("PHY_gpp")

   !surface renewal model — locate GLM-provided globals
   data%id_u_star = aed_locate_sheet_global('u_star')
   data%id_Q_net  = aed_locate_sheet_global('Q_net')

   !# Register diagnostic variables
   IF (diag_level>0) THEN
     !IF( data%simCH4 ) THEN
       data%id_ch4ox   = aed_define_diag_variable('ch4ox','mmol C/m3/d', 'methane oxidation rate')
       data%id_sed_ch4 = aed_define_sheet_diag_variable('ch4_dsf','mmol C/m2/d', &
                            'CH4 exchange across sed/water interface')
       data%id_atm_ch4 = aed_define_sheet_diag_variable('ch4_atm',        &
                            'mmol C/m2/d', 'CH4 exchange across atm/water interface', surf=.TRUE.)
       data%id_ch4_oxic_prod = aed_define_diag_variable('ch4_oxic_prod',  &
                           'mmol/m3/day', 'methane oxic production')
       data%id_sed_ch4_ebb_3d = aed_define_diag_variable('ch4_ebb_dsfv','mmol C/m2/s', & !used for bubble dissolution
                            'CH4 ebullition release rate', zavg=.TRUE., rezero=.FALSE.)
       data%id_kCH4 = aed_define_sheet_diag_variable('kCH4',        &
                            'm/s', 'Piston velocity', surf=.TRUE.)
       data%id_beta = aed_define_sheet_diag_variable('beta',        &
                            'm2/s3', 'Buoyancy flux', surf=.TRUE.)
       data%id_epsilon = aed_define_sheet_diag_variable('epsilon',   &
                            'm2/s3', 'Turbulent kinetic energy dissipation', surf=.TRUE.)

       IF( data%simCH4ebb ) THEN
         IF( data%ebb_model==1) THEN
            data%id_ch4_ebb_df = aed_define_sheet_diag_variable('ch4_ebb_dis','mmol C/m2/d', &
                            'CH4 bubble dissolution rate')
            data%id_sed_ch4_ebb = aed_define_sheet_diag_variable('ch4_ebb_dsf','mmol C/m2/d', &
                            'CH4 ebullition across sed/water interface')
            data%id_atm_ch4_ebb = aed_define_sheet_diag_variable('ch4_ebb_atm', &
                            'mmol C/m2/d', 'CH4 ebullition across atm/water interface', surf=.TRUE.)
         ENDIF
         IF( data%ebb_model==2 ) THEN
            data%id_ch4s_prod = aed_define_sheet_diag_variable('ch4s_prod','mmol C/m2/day', 'sediment ch4 production')
            data%id_ch4s_oxid = aed_define_sheet_diag_variable('ch4s_oxid', 'mmol C/m2/day', 'sediment ch4 oxidation')
            data%id_ch4s_diff = aed_define_sheet_diag_variable('ch4s_diff', 'mmol C/m2/day', 'sediment ch4 diffusion')
            data%id_ch4s_ebb = aed_define_sheet_diag_variable('ch4s_ebb', 'mmol C/m2/day', 'sediment ch4 ebullition')
            data%id_ch4_crit = aed_define_sheet_diag_variable('ch4_crit', 'mmol/m2/day', 'critical sediment ch4 concentration for ebullition')
            data%id_poc_resp = aed_define_sheet_diag_variable('poc_resp', 'mmol/m2/day' , 'poc resperation in the sediment')
            data%id_poc_set = aed_define_sheet_diag_variable('poc_set', 'mmol/m2/day', 'poc settling')
            data%id_phy_set = aed_define_sheet_diag_variable('phy_set', 'mmol/m2/day', 'phy settling')
            data%id_tempzon = aed_define_sheet_diag_variable('temp_zon', 'C degree', 'temperature of different sediment zones')
            data%id_ch4_pw = aed_define_sheet_diag_variable('ch4_pw', 'mmol/m3', 'methane porewater concentration')
       ENDIF
     ENDIF
    ENDIF
  !ENDIF

   !# Register environmental dependencies
   data%id_temp = aed_locate_global('temperature')
   data%id_salt = aed_locate_global('salinity')
   data%id_extc = aed_locate_global('extc_coef')
   data%id_par  = aed_locate_global('par')
   data%id_dz   = aed_locate_global('layer_ht')
   data%id_vel  = aed_locate_global('cell_vel')           ! needed for k600
   data%id_depth= aed_locate_global('depth')
   data%id_layer_ht= aed_locate_global('layer_ht')
   data%id_wind = aed_locate_sheet_global('wind_speed')
   data%id_tau  = aed_locate_sheet_global('taub')
   data%id_depth_c = aed_locate_sheet_global('col_depth')
   data%id_dens = aed_locate_global('density')
   data%id_air_pres = aed_locate_sheet_global('air_pres')
   data%id_area = aed_locate_global('layer_area')


END SUBROUTINE aed_define_methane
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_methane(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of aed_methane model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_methane_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL :: ch4,oxy,temp, dic, gpp
   AED_REAL :: ch4oxidation, prodch4oxic, fch4ox, Tabs
!
!-------------------------------------------------------------------------------
!BEGIN
      ! DIC to be linked to methane oxidation
      dic = _STATE_VAR_(data%id_ch4_oxid_target) 
      ! Retrieve current (local) state variable values.
      ch4 = _STATE_VAR_(data%id_ch4)    ! CH4
      gpp = _DIAG_VAR_(data%id_gpp)

      !# Retrieve current dependent state variable values.
      IF (data%use_oxy) THEN
         oxy = _STATE_VAR_(data%id_oxy) ! O2
      ELSE
         oxy = zero_
      ENDIF

      !# Retrieve current environmental conditions.
      temp = _STATE_VAR_(data%id_temp)  ! temperature
      
      !# Compute rates of change (mmol C/m3/day)
      !ch4oxidation = aed_methane_fch4ox(data%use_oxy,data%Rch4ox,data%Kch4ox,data%vTch4ox,oxy,temp)

      !# Set temporal derivatives
      !_FLUX_VAR_(data%id_ch4) = _FLUX_VAR_(data%id_ch4) + (-ch4*ch4oxidation)

      !# If a linked oxygen pool is present, take oxidation from it assume 1:1 stoichometry
      
      !_FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) - ch4*ch4oxidation  ! *(32./12.) not mass
      !_FLUX_VAR_(data%id_ch4_oxid_target) = _FLUX_VAR_(data%id_ch4_oxid_target) + (ch4*ch4oxidation)

      !New oxidation code
      Tabs = temp + 273.15 ! Kelvin
      !# Compute rates of change (mmol C/m3/day)
      IF (data%MOx_model == 1) THEN !option one is Michealis Menten
         ch4oxidation = aed_methane_fch4ox(data%use_oxy,data%Rch4ox,data%Kch4ox,data%vTch4ox,oxy,temp,data%Kch4ox_ch4,ch4) ! /s
      ELSEIF (data%MOx_model == 2) THEN
         ch4oxidation = aed_methane_fch4ox_Tott(data%use_oxy,ch4,oxy,Tabs)
         !print*, "ch4oxidation", ch4oxidation*secs_per_day
      ELSEIF (data%MOx_model == 3) THEN
         ch4oxidation = aed_methane_fch4ox_NL(data%use_oxy,ch4,oxy,temp,data%Rch4ox,data%vTch4ox,data%beta_ch4,data%alpha_oxy)
      ENDIF

      _FLUX_VAR_(data%id_ch4) = _FLUX_VAR_(data%id_ch4) - ch4oxidation !mmol/m3/s
      !# If a linked oxygen pool is present, take oxidation from it assume 1:1 stoichometry
      _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) - 2.0 * ch4oxidation  ! *(32./12.) not mass
      _FLUX_VAR_(data%id_ch4_oxid_target) = _FLUX_VAR_(data%id_ch4_oxid_target) + ch4oxidation

      !# Export diagnostic variables
      IF (diag_level>0) _DIAG_VAR_(data%id_ch4ox) =  ch4oxidation*secs_per_day

      !Oxic methane production [mmol/m3/d]
      prodch4oxic = gpp * data%ch4_oxic_prod
      !print*, "prodch4oxic", prodch4oxic
      _DIAG_VAR_(data%id_ch4_oxic_prod) = prodch4oxic

      _FLUX_VAR_(data%id_ch4) = _FLUX_VAR_(data%id_ch4) + (prodch4oxic/secs_per_day)

END SUBROUTINE aed_calculate_methane
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_surface_methane(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Air-sea exchange for the aed carbon model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_methane_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp, salt, wind, S, T, depth, vel

   ! State
   AED_REAL :: ch4

   ! Temporary variables

   AED_REAL :: FCH4
   AED_REAL :: Ko, kCH4, CH4solub
   AED_REAL :: Tabs,windHt,atm
   AED_REAL :: q_cub, kin_vsc, dens, u_star, Qnet, denom, numer, thermexp, beta, epsilon, zsurf
   AED_REAL :: A1,A2,A3,A4,B1,B2,B3,logC
   INTEGER  :: iter
   INTEGER  :: ii
   AED_REAL :: karman=0.41, cpl = 4186.0, g=9.8

!-------------------------------------------------------------------------------
!BEGIN

   IF(.NOT.data%simCH4) RETURN

   Ko = 0.

   !----------------------------------------------------------------------------
   !# Get dependent state variables from physical driver
   windHt = 10.
   wind   = _STATE_VAR_S_(data%id_wind) ! Wind speed at 10 m above surface (m/s)
   temp   = _STATE_VAR_(data%id_temp)   ! Temperature (degrees Celsius)
   salt   = _STATE_VAR_(data%id_salt)   ! Salinity (psu)
   depth  = MAX( _STATE_VAR_(data%id_depth), one_ )
   dens = _STATE_VAR_(data%id_dens)
   ch4 = _STATE_VAR_(data%id_ch4)
   Tabs = temp + 273.15
   zsurf = _STATE_VAR_(data%id_layer_ht)

   IF (data%id_vel > 0 ) THEN
     vel = _STATE_VAR_(data%id_vel)
   ELSE
    ! vel  = SQRT(_STATE_VAR_(data%id_E_tau)/_STATE_VAR_(data%id_E_dens))
    ! vel = vel/0.41 * log(depth/0.01)
     vel = 0.0001
   ENDIF
   
   !surface renewal model
   IF (data%ch4_piston_model == 10) THEN
      u_star = _STATE_VAR_S_(data%id_u_star)
      Qnet = _STATE_VAR_S_(data%id_Q_net)
      kin_vsc = (2.414d-5 * (10.0**(247.8/(Tabs-140.0))))/dens !m2/s
      denom = 9.99868d+2 + 1.0d-3*(65.185*temp - 8.4878*temp**2 + 0.05607*temp**3)
      numer = 1.0d-3 * (65.185 - 16.9756*temp + 0.16821*temp**2)
      thermexp = -numer / denom 
      beta = (g * thermexp * Qnet) / (cpl * dens) !m2/s3
      IF (beta > 0.0) THEN
         epsilon = (0.6 * u_star**3)/(karman * zsurf) !m2/s3
      ELSE 
         epsilon = 0.77 * abs(beta) + ((0.56 * u_star**3)/(karman * zsurf)) !m2/s3
      ENDIF
   ENDIF
   !----------------------------------------------------------------------------
   !# CH4 flux
     ! Algorithm from Arianto Santoso <abs11@students.waikato.ac.nz>
   
     ! Piston velocity for CH4
     kCH4 = aed_gas_piston_velocity(wshgt=windHt,wind=wind,tem=temp,sal=salt, &
         vel=vel,depth=depth,epsilon=epsilon,kin_vsc=kin_vsc,LA=data%LA, &
         schmidt_model=data%ch4_schmidt_model,piston_model=data%ch4_piston_model) !need to change it to schmidt model 5?

     _DIAG_VAR_S_(data%id_kCH4) = kCH4
     _DIAG_VAR_S_(data%id_epsilon) = epsilon
     _DIAG_VAR_S_(data%id_beta) = beta
     

     ! Solubility, Ko (mol/L/atm)
     atm = data%atm_ch4   ! 1.76 e-6 !## recent atmospheric CH4 from NOAA (in atm)

     A1 = -415.2807
     A2 =  596.8104
     A3 =  379.2599
     A4 =  -62.0757
     B1 =   -0.05916
     B2 =    0.032174
     B3 =   -0.0048198

     logC = (log(atm)) + A1                                                   &
          + (A2 * (100./Tabs)) + (A3 * log (Tabs/100.)) + (A4 * (Tabs/100.))  &
          + salt * (B1 + (B2  * (Tabs/100.)) + (B3 * (Tabs/100.)*(Tabs/100.)))

     CH4solub = exp(logC) * 1e-3
     
     !# Now compute methane atm flux  (mmol/m2/s = m/s * mmol/m3)
     FCH4 = kCH4 *  (ch4 - CH4solub)

     !----------------------------------------------------------------------------
     !# Transfer surface exchange value to AED (mmmol/m2) converted by driver.
     _FLUX_VAR_T_(data%id_ch4) = -FCH4

     _DIAG_VAR_S_(data%id_atm_ch4) = FCH4*secs_per_day
     
   !----------------------------------------------------------------------------

END SUBROUTINE aed_calculate_surface_methane
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_benthic_methane(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED carbon.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_methane_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp, par, extc, dz, depth, yearday, depth_c, dens, tstart, Tabs, Hch4, air_pres, matz
   ! State
   AED_REAL :: oxy, pocs, ch4, ch4s

   AED_REAL :: balance, q_cub

   ! Temporary variables
   AED_REAL :: ch4_flux, Fsed_ch4, ebb_flux, Fsed_ch4_ebb, ch4_bub_disf, poc_resp, ch4_prod, gamma, ch4_diff, ch4_ox_sed, ptot, ch4_ebb, ch4_cr, k_HL, ch4_pw, oxy_sed, layer_ht, depth_hydro_p, poc_swi, phy_swi, poc_set, poc_burial, mu, mu_temp, beta, layer_area
   AED_REAL :: poc_resp_potential, pN2, ptotal
   !AED_REAL, PARAMETER :: maxMPBProdn = 40.     ! mmolC/m2/day                     !
   !AED_REAL, PARAMETER :: IkMPB       = 180.0   ! Light sensitivity of MPB  !
   AED_REAL :: ch4_bub_s, excess_mass, excess, k_release, ch4_eff !new bubble storage variables
   AED_REAL :: ox_SWI, ch4_diff_max 
!-------------------------------------------------------------------------------
!BEGIN

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   temp = _STATE_VAR_(data%id_temp) ! local temperature
   par  = _STATE_VAR_(data%id_par)  ! local par
   dz   = _STATE_VAR_(data%id_dz)   ! local layer thickness
   depth= _STATE_VAR_(data%id_depth)! local layer depth
   extc = _STATE_VAR_(data%id_extc) ! local extinction
   oxy = _STATE_VAR_(data%id_oxy) ! local oxygen
    ! Retrieve current (local) state variable values.
   layer_ht = _STATE_VAR_(data%id_layer_ht)
   depth_c = _STATE_VAR_S_(data%id_depth_c)
   layer_area = _STATE_VAR_(data%id_area)
   matz = _STATE_VAR_S_(data%id_zone)
   _DIAG_VAR_S_(data%id_tempzon) = temp
   

   IF ( data%use_sed_model_ch4 ) THEN
      Fsed_ch4 = _STATE_VAR_S_(data%id_Fsed_ch4) / secs_per_day
      IF( data%simCH4ebb ) Fsed_ch4_ebb = _STATE_VAR_S_(data%id_Fsed_ch4_ebb) / secs_per_day
   ELSE
      Fsed_ch4 = data%Fsed_ch4
      IF( data%simCH4ebb ) Fsed_ch4_ebb = data%Fsed_ch4_ebb
   ENDIF

   IF( data%simCH4ebb ) THEN
     IF( data%ebb_model == 2 ) THEN
      ! POC Settling
         poc_swi = _DIAG_VAR_(data%id_poc_swi)/secs_per_day !POC settling
         phy_swi = _DIAG_VAR_S_(data%id_phy_swi)/secs_per_day !PHY settling
         pocs = _STATE_VAR_S_(data%id_poc_s)

         poc_set = poc_swi + phy_swi
         _DIAG_VAR_S_(data%id_poc_set) = poc_swi * secs_per_day
         _DIAG_VAR_S_(data%id_phy_set) = phy_swi * secs_per_day
         _FLUX_VAR_B_(data%id_poc_s) = _FLUX_VAR_B_(data%id_poc_s) - (phy_swi) !5.8e-04 need actual settling model average 2.581018e-06
         _FLUX_VAR_B_(data%id_poc_s) = _FLUX_VAR_B_(data%id_poc_s) - (poc_swi)

         !POC respiration
         poc_resp = data%Rpoc_resp * pocs * (data%theta_poc_resp**(temp-20.0))
         _DIAG_VAR_S_(data%id_poc_resp) = poc_resp * secs_per_day
         _FLUX_VAR_B_(data%id_poc_s) = _FLUX_VAR_B_(data%id_poc_s) - (poc_resp)

         !POC burial
         poc_burial = pocs / data%tau_poc_burial
         _FLUX_VAR_B_(data%id_poc_s) = _FLUX_VAR_B_(data%id_poc_s) - (poc_burial)
         
         !dissolved methane and bubble state variables
         ch4s = _STATE_VAR_S_(data%id_ch4_s)
         ch4_bub_s = _STATE_VAR_S_(data%id_ch4_bub_s)
         !methane porewater diagnostic
         ch4_pw = ch4s / (data%Phi * data%sed_z) 
         
         !Critical methane concentration for bubble formation
         dens = _STATE_VAR_(data%id_dens) !m3/kg
         air_pres = _STATE_VAR_S_(data%id_air_pres) * 100 !convert from hPa to Pa
         Tabs = temp + 273.15 ! Kelvin
         Hch4 = 1.4e-5 * exp(1600.0*((1/Tabs)-(1/298.0))) * 8.3144621 * Tabs 
         depth_hydro_p = depth_c - depth + (layer_ht/2)
         ptotal = air_pres + (dens * 9.81 * depth_hydro_p)
         pN2 = 0.79 * air_pres
         ch4_cr = Hch4 * ((ptotal-pN2)/(8.3144621 * Tabs))
         _DIAG_VAR_S_(data%id_ch4_crit) = ch4_cr * 1000! mol to mmol

         !ch4 production
         oxy = _STATE_VAR_(data%id_oxy)
         oxy_sed = data%oxy_sed_p * _STATE_VAR_(data%id_oxy) ! sediment oxy ~ X% of bottom layer oxy
         gamma = min(data%Kpocox_sed/(data%Kpocox_sed+oxy_sed), data%gamma_max)
         ch4_prod = gamma * poc_resp
         _DIAG_VAR_S_(data%id_ch4s_prod) = ch4_prod * secs_per_day
         _FLUX_VAR_B_(data%id_ch4_s) = _FLUX_VAR_B_(data%id_ch4_s) + (ch4_prod)
         
         !ch4 oxidation
         !ch4_ox_sed = data%Rch4ox_sed * oxy_sed/(oxy_sed+data%Kch4ox_sed) * (data%vTch4ox_sed**(temp-20.0)) * ch4s !mmol/m2/s
         !_DIAG_VAR_S_(data%id_ch4s_oxid) = ch4_ox_sed * secs_per_day !convert to /day
         !_FLUX_VAR_B_(data%id_ch4_s) = _FLUX_VAR_B_(data%id_ch4_s) - (ch4_ox_sed)

         !ch4 diffusion
         !ch4_pw = ch4s / (data%Phi * data%sed_z) 
         
         !ch4_diff = max(0.0, (data%ch4_s_diff_coef/data%z_dbl)*(ch4_pw - ch4)) !* (data%Kch4ox_diff/(oxy+data%Kch4ox_diff))
         !_DIAG_VAR_S_(data%id_ch4s_diff) = ch4_diff * secs_per_day ! convert to /day
         !_FLUX_VAR_B_(data%id_ch4_s) = _FLUX_VAR_B_(data%id_ch4_s) - (ch4_diff)


         !original ebullition formulation bubbles get released immadiately when poreater CH4 exceeds critical conc.
         k_HL = log(2.0) / 1800 ! ebullition half life per second !k_HL = 0.000001
         ch4_ebb = MAX(0.0, k_HL * (ch4_pw - 1000*ch4_cr) * data%sed_z) ! mmol/m2/s?
         _DIAG_VAR_S_(data%id_ch4s_ebb) = ch4_ebb * secs_per_day
         _FLUX_VAR_B_(data%id_ch4_s) = _FLUX_VAR_B_(data%id_ch4_s) - (ch4_ebb)

         
         !bubble formation
         !excess_mass = MAX(0.0, (ch4_pw - 1000*ch4_cr) * data%sed_z)  
         !excess = excess_mass / 3600.0
         !_FLUX_VAR_B_(data%id_ch4_s) = _FLUX_VAR_B_(data%id_ch4_s) - (excess)
         !_FLUX_VAR_B_(data%id_ch4_bub_s) = _FLUX_VAR_B_(data%id_ch4_bub_s) + (excess)
         
         
         !bubble release
         !k_release = log(2.0) / (2.0 * 86400.0)   ! ~2 days (start here)
         !ch4_ebb = k_release * ch4_bub_s    ! mmol m⁻2 s⁻1
         !_DIAG_VAR_S_(data%id_ch4s_ebb) = ch4_ebb * secs_per_day
         !_FLUX_VAR_B_(data%id_ch4_bub_s) = _FLUX_VAR_B_(data%id_ch4_bub_s) - (ch4_ebb)
         
         !Diffusion
         !ch4_diff = max(0.0, (data%ch4_s_diff_coef/data%z_dbl)*(ch4_pwd - ch4)) !* (data%Kch4ox_diff/(oxy+data%Kch4ox_diff))
         ch4_eff = min(ch4_pw, 1000*ch4_cr)
         ch4_diff = max(0.0, (data%ch4_s_diff_coef/data%z_dbl)*(ch4_eff - ch4)) * (data%Kch4ox_diff/(oxy+data%Kch4ox_diff))
         _DIAG_VAR_S_(data%id_ch4s_diff) = ch4_diff * secs_per_day ! convert to /day
         _FLUX_VAR_B_(data%id_ch4_s) = _FLUX_VAR_B_(data%id_ch4_s) - (ch4_diff)
         _DIAG_VAR_S_(data%id_ch4_pw) = ch4_pw
      ENDIF
   ENDIF

   
     IF (data%use_oxy) THEN
       ! Sediment flux dependent on oxygen and temperature
       oxy = _STATE_VAR_(data%id_oxy)
       IF (data%ebb_model < 2) THEN 
         ch4_flux = Fsed_ch4 * data%Ksed_ch4/(data%Ksed_ch4+oxy) * (data%theta_sed_ch4**(temp-20.0))
       ELSE
         ch4_flux = ch4_diff
       ENDIF
       IF( data%simCH4ebb ) THEN
          IF( data%ebb_model == 1 ) THEN 
             ebb_flux = Fsed_ch4_ebb * (data%theta_sed_ch4**(temp-20.0))
          ELSEIF (data%ebb_model == 2) THEN
             ebb_flux = ch4_ebb
          ENDIF
       ENDIF
       
       !IF( data%simCH4ebb ) THEN
        ! ebb_flux = Fsed_ch4_ebb * (data%theta_sed_ch4**(temp-20.0))
         ! Kinneret special lake level equations
         !ebb_flux = ebb_flux * data%ch4_bub_cLL * exp(data%ch4_bub_kLL*(data%ch4_bub_aLL-depth))
      ENDIF


   ! Re-distribute bubbles to the water or atmosphere, or dissolve
   IF( data%simCH4ebb ) THEN
      IF ( data%ebb_model == 1 ) THEN 
         ! Add bubbles to layer
         _FLUX_VAR_(data%id_ch4_bub) = _FLUX_VAR_(data%id_ch4_bub) + (ebb_flux)
         ! Dissolve bubbles in this bottom layer, depending on depth
         ch4_bub_disf = data%ch4_bub_disf1
         IF( depth > data%ch4_bub_disdp) ch4_bub_disf = data%ch4_bub_disf2
         _FLUX_VAR_(data%id_ch4) = _FLUX_VAR_(data%id_ch4) + ebb_flux*ch4_bub_disf !/ dz
      ENDIF

      IF (diag_level>0) THEN
         IF ( data%ebb_model == 1) THEN
            _DIAG_VAR_S_(data%id_ch4_ebb_df) = (ebb_flux*ch4_bub_disf) * secs_per_day ! MH Something wrong with dz here?
               ! Release the remainder to the atmosphere (mmol/m2/day)
            _DIAG_VAR_S_(data%id_atm_ch4_ebb) = _DIAG_VAR_S_(data%id_atm_ch4_ebb)+ ebb_flux * (1.-ch4_bub_disf) * secs_per_day
               ! Note the bubble flux, as the zone sees it  (mmol/m2/day)
            _DIAG_VAR_S_(data%id_sed_ch4_ebb) = ebb_flux * secs_per_day
               !print*, 'sed_ch4_ebb', _DIAG_VAR_S_(data%id_sed_ch4_ebb)
         ENDIF
            ! Note the bubble flux, as the water sees it  (mmol/m2/s)
            _DIAG_VAR_(data%id_sed_ch4_ebb_3d) = ebb_flux 
      ENDIF
   ENDIF

    ! Set bottom fluxes for the pelagic (flux per surface area, per second)
   ! Increment sediment flux value into derivative of water column variable
   _FLUX_VAR_(data%id_ch4) = _FLUX_VAR_(data%id_ch4) + (ch4_flux)
   
   ! Store dissolved sediment fluxes as diagnostic variables (flux per surface area, per day
   _DIAG_VAR_S_(data%id_sed_ch4) = ch4_flux * secs_per_day


END SUBROUTINE aed_calculate_benthic_methane
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
PURE AED_REAL FUNCTION aed_methane_fch4ox(use_oxy,Rch4ox,Kch4ox,vTch4ox,oxy,temp,Kch4ox_ch4,ch4)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for methane oxidation
!
! Here, the classical Michaelis-Menten formulation for nitrification
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in)  :: use_oxy
   AED_REAL,INTENT(in) :: Rch4ox,Kch4ox,vTch4ox,oxy,temp,Kch4ox_ch4,ch4
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (use_oxy) THEN
      aed_methane_fch4ox = Rch4ox * oxy/(Kch4ox+oxy) * ch4/(Kch4ox_ch4+ch4) * (vTch4ox**(temp-20.0))
   ELSE
      aed_methane_fch4ox = Rch4ox * (vTch4ox**(temp-20.0))
   ENDIF

END FUNCTION aed_methane_fch4ox
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
PURE AED_REAL FUNCTION aed_methane_fch4ox_Tott(use_oxy,ch4,oxy,Tabs)
!-------------------------------------------------------------------------------
! Non-linear oxygen function
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in)  :: use_oxy
   AED_REAL,INTENT(in) :: ch4,oxy,Tabs
   AED_REAL :: b0,b1,b2,a,k,ln_rate, oxyk
!-------------------------------------------------------------------------------
b0 = 20.8
b1 = 0.79
b2 = -5669.61
!a = 0.01
!k = 0.18

IF (use_oxy) THEN
   ln_rate = b0 + b1 * log(ch4) + b2 * (1.0 / Tabs) &
            + log(exp(-0.01 * oxy) - exp(-0.19 * oxy))
   
   aed_methane_fch4ox_Tott = exp(ln_rate)

   aed_methane_fch4ox_Tott = aed_methane_fch4ox_Tott/secs_per_day
ELSE
   aed_methane_fch4ox_Tott = zero_
ENDIF

END FUNCTION aed_methane_fch4ox_Tott
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
PURE AED_REAL FUNCTION aed_methane_fch4ox_NL(use_oxy,ch4,oxy,temp,Rch4ox,vTch4ox,beta_ch4,alpha_oxy)
!-------------------------------------------------------------------------------
! Non-linear oxygen function
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in)  :: use_oxy
   AED_REAL,INTENT(in) :: ch4,oxy,temp,Rch4ox,vTch4ox,beta_ch4,alpha_oxy
!-------------------------------------------------------------------------------
IF (use_oxy) THEN
   aed_methane_fch4ox_NL = Rch4ox * (ch4**beta_ch4) * exp(-alpha_oxy * oxy)*(1-exp(-18.0* alpha_oxy * oxy)) * (vTch4ox**(temp-20.0))
ELSE
   aed_methane_fch4ox_NL = zero_
ENDIF

END FUNCTION aed_methane_fch4ox_NL
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE aed_methane
