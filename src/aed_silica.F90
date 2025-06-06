!###############################################################################
!#                                                                             #
!# aed_silica.F90                                                              #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
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
!# Created 24 August 2011                                                      #
!#                                                                             #
!###############################################################################
!                                                                              !
!         .----------------.  .----------------.  .----------------.           !
!         | .--------------. || .--------------. || .--------------. |         !
!         | |    _______   | || |     _____    | || |   _____      | |         !
!         | |   /  ___  |  | || |    |_   _|   | || |  |_   _|     | |         !
!         | |  |  (__ \_|  | || |      | |     | || |    | |       | |         !
!         | |   '.___`-.   | || |      | |     | || |    | |   _   | |         !
!         | |  |`\____) |  | || |     _| |_    | || |   _| |__/ |  | |         !
!         | |  |_______.'  | || |    |_____|   | || |  |________|  | |         !
!         | |              | || |              | || |              | |         !
!         | '--------------' || '--------------' || '--------------' |         !
!         '----------------'  '----------------'  '----------------'           !
!                                                                              !
!###############################################################################

#include "aed.h"

MODULE aed_silica
!-------------------------------------------------------------------------------
! Silica module contains basic equations for uptake and release
!-------------------------------------------------------------------------------
   USE aed_core

   IMPLICIT NONE

   PRIVATE   ! By default make everything private.
!
   PUBLIC aed_silica_data_t
!
   TYPE,extends(aed_model_data_t) :: aed_silica_data_t
      !# Variable identifiers
      INTEGER  :: id_rsi,id_oxy
      INTEGER  :: id_Fsed_rsi
      INTEGER  :: id_temp
      INTEGER  :: id_sed_rsi

      !# Model parameters
      AED_REAL :: Fsed_rsi,Ksed_rsi,theta_sed_rsi
      LOGICAL  :: use_oxy,use_sed_model

     CONTAINS
         PROCEDURE :: define            => aed_define_silica
         PROCEDURE :: calculate         => aed_calculate_silica
         PROCEDURE :: calculate_benthic => aed_calculate_benthic_silica
!        PROCEDURE :: mobility          => aed_mobility_silica
!        PROCEDURE :: light_extinction  => aed_light_extinction_silica
!        PROCEDURE :: delete            => aed_delete_silica

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
SUBROUTINE aed_define_silica(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!  Here, the aed namelist is read and the variables exported
!  by the model are registered with AED
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)                      :: namlst
   CLASS(aed_silica_data_t),INTENT(inout) :: data

!
!LOCALS
   INTEGER  :: status

!  %% NAMELIST   %%  /aed_silica/
!  %% Last Checked 20/08/2021
   AED_REAL          :: rsi_initial=4.5
   AED_REAL          :: rsi_min=zero_
   AED_REAL          :: rsi_max=nan_
   AED_REAL          :: Fsed_rsi = 3.5
   AED_REAL          :: Ksed_rsi = 30.0
   AED_REAL          :: theta_sed_rsi = 1.0
   CHARACTER(len=64) :: silica_reactant_variable=''
   CHARACTER(len=64) :: Fsed_rsi_variable=''
! %% From Module Globals
!  INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
!                                             ! 1 = basic diagnostic outputs
!                                             ! 2 = flux rates, and supporitng
!                                             ! 3 = other metrics
!                                             !10 = all debug & checking outputs
!  %% END NAMELIST   %%  /aed_silica/

   NAMELIST /aed_silica/ rsi_initial,rsi_min,rsi_max,Fsed_rsi,Ksed_rsi,theta_sed_rsi, &
                         silica_reactant_variable,Fsed_rsi_variable, diag_level
!
!-------------------------------------------------------------------------------
!BEGIN
   print *,"        aed_silica configuration"

   ! Read the namelist
   read(namlst,nml=aed_silica,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist for &aed_silica'

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   data%Fsed_rsi  = Fsed_rsi/secs_per_day
   data%Ksed_rsi  = Ksed_rsi
   data%theta_sed_rsi = theta_sed_rsi
   data%use_oxy = silica_reactant_variable .NE. '' !This means oxygen module switched on

   ! Register state variables
   data%id_rsi = aed_define_variable('rsi','mmol Si/m3', 'silica',     &
                                    rsi_initial,minimum=rsi_min,maximum=rsi_max)

   ! Register external state variable dependencies
   IF (data%use_oxy) &
      data%id_oxy = aed_locate_variable(silica_reactant_variable)

   data%use_sed_model = Fsed_rsi_variable .NE. ''
   IF (data%use_sed_model) &
      data%id_Fsed_rsi = aed_locate_sheet_variable(Fsed_rsi_variable)

   ! Register diagnostic variables
   data%id_sed_rsi = aed_define_sheet_diag_variable('dsf_rsi','mmol/m2/d', &
                     'Si exchange across sed/water interface')

   ! Register environmental dependencies
   data%id_temp = aed_locate_global('temperature')
END SUBROUTINE aed_define_silica
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed_calculate_silica(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of aed_silica model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_silica_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!ARGUMENT
   !AED_REAL                   :: rsi,oxy,temp,tss !State variables
!
!-------------------------------------------------------------------------------
!BEGIN


END SUBROUTINE aed_calculate_silica
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_benthic_silica(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED silica.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_silica_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp

   ! State
   AED_REAL :: rsi,oxy

   ! Temporary variables
   AED_REAL :: rsi_flux, Fsed_rsi, fT, fDO
!
!-------------------------------------------------------------------------------
!BEGIN

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   temp = _STATE_VAR_(data%id_temp) ! local temperature

    ! Retrieve current (local) state variable values.
   rsi = _STATE_VAR_(data%id_rsi)   ! silica

   ! Compute the sediment flux dependent on overlying oxygen & temperature
   fT = data%theta_sed_rsi**(temp-20.0)
   IF (data%use_oxy) THEN
      ! Sediment flux dependent on oxygen and temperature
      oxy = _STATE_VAR_(data%id_oxy)
      fDO = data%Ksed_rsi/(data%Ksed_rsi+oxy)
   ELSE
      ! Sediment flux dependent on temperature only
      fDO = one_
   ENDIF

   ! Set the flux dependent on if and how it is linked
   IF (data%use_sed_model) THEN
     ! Linked to aed_sedflux, check if its constant or dynamically set
     IF ( aed_is_const_var(data%id_Fsed_rsi) ) THEN
        Fsed_rsi = _DIAG_VAR_S_(data%id_Fsed_rsi) * MIN(3.,fDO * fT) / secs_per_day
     ELSE
        Fsed_rsi = _DIAG_VAR_S_(data%id_Fsed_rsi) / secs_per_day
     ENDIF
   ELSE
     Fsed_rsi = data%Fsed_rsi * MIN(3.,fDO * fT)
   ENDIF
   rsi_flux = Fsed_rsi

   ! Set bottom fluxes for the pelagic (change per surface area per second)
   _FLUX_VAR_(data%id_rsi) = _FLUX_VAR_(data%id_rsi) + (rsi_flux)

   ! Set sink and source terms for the benthos (change per surface area per second)
   ! Note that this should include the fluxes to and from the pelagic.
   !_FLUX_VAR_B_(data%id_ben_rsi) = _FLUX_VAR_B_(data%id_ben_rsi) + (-rsi_flux)

   ! Also store sediment flux as diagnostic variable.
   _DIAG_VAR_S_(data%id_sed_rsi) = rsi_flux * secs_per_day


END SUBROUTINE aed_calculate_benthic_silica
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_silica
