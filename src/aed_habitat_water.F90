!###############################################################################
!#                                                                             #
!# aed_habitat.F90                                                             #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2016 - 2020 -  The University of Western Australia               #
!#                                                                             #
!#   GLM is free software: you can redistribute it and/or modify               #
!#   it under the terms of the GNU General Public License as published by      #
!#   the Free Software Foundation, either version 3 of the License, or         #
!#   (at your option) any later version.                                       #
!#                                                                             #
!#   GLM is distributed in the hope that it will be useful,                    #
!#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
!#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
!#   GNU General Public License for more details.                              #
!#                                                                             #
!#   You should have received a copy of the GNU General Public License         #
!#   along with this program.  If not, see <http://www.gnu.org/licenses/>.     #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Created March 2016                                                          #
!#                                                                             #
!###############################################################################

#include "aed+.h"

!
MODULE aed_habitat_water
!-------------------------------------------------------------------------------
! aed_habitat --- habitat model
!
!-------------------------------------------------------------------------------
   USE aed_core
   USE aed_util

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed_habitat_water_data_t
!
   TYPE,extends(aed_model_data_t) :: aed_habitat_water_data_t
      INTEGER :: num_habitats
      !# Variable identifiers
      INTEGER :: id_bird, id_mtox
      INTEGER :: id_chsi, id_chpl, id_chfl, id_chsd
      INTEGER :: id_fssi, id_fsdp, id_fssd, id_fsst, id_fsls, id_fsdw, id_fsmt
      INTEGER :: id_rhsi, id_rhpl, id_rhfl, id_rhsd, id_rhtr, id_rhsp
      INTEGER :: id_chhsi, id_chhpl, id_chhfl, id_chhsd, id_chhtr, id_chhsp
      INTEGER :: id_crhsi, id_crhpl, id_crhfl, id_crhsd, id_crhtr
      INTEGER :: id_pshsi, id_pshsd, id_pshpl, id_pshfl
      INTEGER :: id_wettime, id_drytime
      INTEGER, ALLOCATABLE :: id_d_rupfs(:),id_d_rupft(:),id_d_rupfl(:),id_d_rupfa(:),id_d_rupfd(:)
      INTEGER, ALLOCATABLE :: id_d_chafs(:),id_d_chaft(:),id_d_chafl(:),id_d_chafa(:),id_d_chafd(:),id_d_chafv(:)
      INTEGER, ALLOCATABLE :: id_d_crcfm(:),id_d_crcft(:),id_d_crcfv(:),id_d_crcfs(:),id_d_crcfd(:),id_d_crcfp(:)
      INTEGER, ALLOCATABLE :: id_d_pssfm(:), id_d_pssft(:),id_d_pssfv(:),id_d_pssfs(:),id_d_pssfd(:), id_d_pssfu(:)

      !# Dependencies
      INTEGER :: id_l_ph, id_l_hab, id_l_aass, id_l_rveg, id_l_bveg
      INTEGER :: id_l_salg, id_l_falg, id_d_turb, id_l_ncs1, id_l_ncs2, id_l_tau0
      INTEGER :: id_l_otrc, id_l_oxy, id_l_sav
      INTEGER :: id_l_svwc, id_l_stmp25, id_l_stmp, id_l_veg1, id_l_veg2, id_l_pass
      INTEGER, ALLOCATABLE :: id_l_mtox(:)

      !# Environment variables
      INTEGER :: id_E_temp, id_E_salt, id_E_bathy, id_E_matz, id_E_depth
      INTEGER :: id_E_nearlevel, id_E_extc, id_E_Io, id_E_stress, id_E_airtemp

      !# Model switches
      LOGICAL :: simBirdForaging,simBenthicProd,simFishTolerance,simGalaxiidSpawning
      LOGICAL :: simCrabHabitat,simRuppiaHabitat,simCharaHabitat
      LOGICAL :: simMosquitoRisk,simCyanoRisk
      LOGICAL :: simMetalTox,simClearWater
      LOGICAL :: simCrocEggs,simPassiflora

      !# Model parameters
      AED_REAL, ALLOCATABLE :: mtox_lims(:)
      INTEGER :: num_mtox, n_zones_chara, n_zones_fishspawn
      INTEGER :: n_zones_crocs, n_zones_pass
      AED_REAL,ALLOCATABLE :: active_zones_chara(:), active_zones_fishspawn(:)
      AED_REAL,ALLOCATABLE :: active_zones_crocs(:), active_zones_pass(:)


     CONTAINS
         PROCEDURE :: define             => aed_define_habitat_water
!        PROCEDURE :: calculate          => aed_calculate_habitat_water
!        PROCEDURE :: calculate_benthic  => aed_calculate_benthic_habitat_water
         PROCEDURE :: calculate_riparian => aed_calculate_riparian_habitat_water
!        PROCEDURE :: mobility           => aed_mobility_habitat_water
!        PROCEDURE :: light_extinction   => aed_light_extinction_habitat_water
!        PROCEDURE :: delete             => aed_delete_habitat_water

   END TYPE

!-------------------------------------------------------------------------------
!MODULE VARIABLES
   AED_REAL, PARAMETER :: DDT = 0.25/24.    ! Currently assuming 15 min timestep
   LOGICAL :: extra_diag
   INTEGER :: diag_level = 10

!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed_define_habitat_water(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED HABITAT (WATER) module
!
!  Here, the aed namelist is read and the variables exported
!  are registered with AED.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed_habitat_water_data_t),INTENT(inout) :: data
!
!LOCALS
   INTEGER :: i, z, status, num_mtox
   INTEGER :: n_zones_chara = 0, active_zones_chara(MAX_ZONES)
   INTEGER :: n_zones_fishspawn = 0,active_zones_fishspawn(MAX_ZONES)
   LOGICAL :: simFishTolerance, simGalaxiidSpawning, simBenthicProd, &
              simCyanoRisk, simMosquitoRisk, simCrabHabitat, simClearWater,      &
              simRuppiaHabitat, simCharaHabitat, simMetalTox
   AED_REAL          :: mtox_lims(10)
   CHARACTER(len=64) :: bird_acid_link, bird_habs_link, bird_aass_link, bird_rveg_link, bird_bveg_link
   CHARACTER(len=64) :: fshsi_veg_link, fshsi_oxy_link, fshsi_otrc_link
   CHARACTER(len=64) :: chsi_otrc_link, chsi_oxy_link, chsi_veg_link
   CHARACTER(len=64) :: chhsi_salg_link,chhsi_falg_link
   CHARACTER(len=64) :: chhsi_ncs1_link,chhsi_ncs2_link,chhsi_tau0_link
   CHARACTER(len=64) :: crhsi_ncs1_link,crhsi_ncs2_link
   CHARACTER(len=64) :: crhsi_stmp_link,crhsi_svwc_link,crhsi_pass_link
   CHARACTER(len=64) :: pshsi_stmp_link,pshsi_svwc_link,pshsi_veg1_link,pshsi_veg2_link
   CHARACTER(len=64) :: pshsi_ncs1_link, pshsi_ncs2_link
   CHARACTER(len=64) :: rhsi_salg_link, rhsi_falg_link
   CHARACTER(len=40) :: mtox_acid_link, mtox_aass_link
   CHARACTER(len=40) :: mtox_vars(10)

   NAMELIST /aed_habitat_water/                                                &
                               simFishTolerance,                               &
                               simClearWater,                                  &
                               simCyanoRisk,                                   &
                               simMetalTox, mtox_vars, mtox_lims,              &
                               diag_level
!
!-------------------------------------------------------------------------------
!BEGIN
   print *,"        aed_habitat_water initialization"
   print *," WARNING! aed_habitat model is currently under development"

   ! Default
   simFishTolerance = .false.
   simClearWater = .false.
   simCyanoRisk = .false.
   simMetalTox = .false.


   ! Read the namelist
   read(namlst,nml=aed_habitat_water,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_habitat'

   ! Update module level switches
   data%num_habitats = 0
   data%simFishTolerance = simFishTolerance ; IF(simFishTolerance) data%num_habitats=data%num_habitats+1
   data%simMetalTox      = simMetalTox      ; IF(simMetalTox) data%num_habitats=data%num_habitats+1
   data%simCyanoRisk     = simCyanoRisk     ; IF(simCyanoRisk) data%num_habitats=data%num_habitats+1
   data%simClearWater    = simClearWater    ; IF(simClearWater) data%num_habitats=data%num_habitats+1

   print *,"          ... # habitat templates simulated: ",data%num_habitats

   !----------------------------------------------------------------------------
   ! Define variables and dependencies


   !-- CONTAMINATION
   IF( simMetalTox ) THEN
     data%id_mtox =  aed_define_sheet_diag_variable('toxicity','-', 'Suitability')

     mtox_acid_link = 'CAR_pH'
     mtox_aass_link = 'ASS_uzaass'

     mtox_vars = '' ;  mtox_lims = 1.0
     DO i=1,10 ; IF (mtox_vars(i)  .EQ. '' ) THEN ; num_mtox = i-1 ; EXIT ; ENDIF ; ENDDO
     ALLOCATE(data%id_l_mtox(num_mtox)); ALLOCATE(data%mtox_lims(num_mtox))
     data%num_mtox = num_mtox
     DO i=1,data%num_mtox
       data%id_l_mtox(i) =  aed_locate_variable(mtox_vars(i))
       data%mtox_lims(i) =  mtox_lims(i)
       !print*,'Tox : ', TRIM(tfe_vars(i)), ' * ', data%tfe_varscale(i)
     ENDDO
   ENDIF

!   !-- GENERAL
!   IF( simGalaxiidSpawning .OR. simCharaHabitat .OR. simRuppiaHabitat) THEN
!     data%id_wettime = aed_define_sheet_diag_variable('wettime','d','time cell has been innundated')
!     data%id_drytime = aed_define_sheet_diag_variable('drytime','d','time cell has been exposed')
!   ENDIF

   ! Register environmental dependencies
   data%id_E_salt  = aed_locate_global('salinity')
   data%id_E_extc  = aed_locate_global('extc_coef')
   data%id_E_temp  = aed_locate_global('temperature')
   data%id_E_depth = aed_locate_global('layer_ht')
   data%id_E_bathy     = aed_locate_global_sheet('bathy')
   data%id_E_matz      = aed_locate_global_sheet('material')
   data%id_E_Io        = aed_locate_global_sheet('par_sf')
   data%id_E_airtemp   = aed_locate_global_sheet('air_temp')
   data%id_E_stress    = aed_locate_global_sheet('taub')
   data%id_E_nearlevel = aed_locate_global_sheet('nearest_depth')

END SUBROUTINE aed_define_habitat_water
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed_calculate_riparian_habitat_water(data,column,layer_idx,pc_wet)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED habitat.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_habitat_water_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(in) :: pc_wet
!
!LOCALS
   ! Environment
   AED_REAL :: temp, salt, wlevel, extc, bathy, matz, Io, Ig, vel, stress, tau0, stem05, stem25

   ! State
   AED_REAL :: depth, ph, hab, sdepth, uzaass, aass, conc, mtox, turb, grav, stem, svwc

   ! Temporary variables
   INTEGER  :: i
   AED_REAL :: euphotic, drytime, light

   ! Parameters
   AED_REAL, PARAMETER :: crit_soil_acidity = 100.000
   AED_REAL, PARAMETER :: crit_water_ph = 6.0
   AED_REAL, PARAMETER :: crit_soil_type = 5.5 ! (matz 6&7 is sand)
   AED_REAL, PARAMETER :: crit_rveg_depth = 0.6
   AED_REAL, PARAMETER :: sav_height = 0.1 !assume plant are 10cm to middle
   AED_REAL, PARAMETER :: crit_salinity = 50.0
   AED_REAL, PARAMETER :: crit_leg_depth = 0.12
   AED_REAL, PARAMETER :: crit_hab_conc = 500.

   AED_REAL :: fs_sdepth , fs_substr, fs_spntem, fs_stress, fs_dewatr, fs_mattem

   AED_REAL :: rhpl,rhfl,rhsd,rhtr,rhsp = 0.,falg
   AED_REAL :: pshpl, pshfl, pshsd, pass, height
   AED_REAL :: crns = 0.,creg = 0.,crht = 0.,crml = 0.
   AED_REAL :: limitation(5,6)

!-------------------------------------------------------------------------------
!BEGIN
   matz = 0.0 ; salt = 0.0 ; euphotic = 0.0 ; bathy = 0.0  !## CAB [-Wmaybe-uninitialized]

   !---------------------------------------------------------------------------+
   !-- HABITAT TEMPLATE X:
!   IF( data%simMosquitoRisk ) THEN
!
!   !  _DIAG_VAR_S_(data%id_fish) = mtox
!
!   ENDIF




END SUBROUTINE aed_calculate_riparian_habitat_water
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





END MODULE aed_habitat_water
