!###############################################################################
!#                                                                             #
!# aed_light.F90                                                               #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
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
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Created April 2022                                                          #
!#                                                                             #
!###############################################################################

!#------------------------------##################------------------------------
!### This module is not used yet and is under development. Everything may change
!#------------------------------##################------------------------------

#include "aed.h"

!
MODULE aed_light

   USE aed_core
!  USE aed_common

   IMPLICIT NONE

   PRIVATE

   LOGICAL :: link_water_clarity = .FALSE.

!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE BioExtinction(column, count, extc)
!-------------------------------------------------------------------------------
!
! Calculate the specific light attenuation additions due to AED modules
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed_column_t), INTENT(inout) :: column(:)
   INTEGER,  INTENT(in)    :: count
   AED_REAL, INTENT(inout) :: extc(:)
!
#if 0
!LOCAL VARIABLES:
   INTEGER :: i
   AED_REAL :: localext
!
!-------------------------------------------------------------------------------
!BEGIN
   localext = zero_

   CALL aed_light_extinction(column, 1, localext)
   IF (link_water_clarity) THEN
      extc(1) = localext
   ELSE
      extc(1) = localext + Kw
   END IF

   IF (count <= 1) RETURN

   DO i = 2, count
      CALL aed_light_extinction(column, i, localext)
      IF (link_water_clarity) THEN
         extc(i) = localext
      ELSE
         extc(i) = localext + Kw
      END IF
   ENDDO
#endif
END SUBROUTINE BioExtinction
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE Light(column, count, Io, extc, par_, h_)
!-------------------------------------------------------------------------------
!
! Calculate photosynthetically active radiation over entire column
! based on surface radiation, and background and biotic extinction.
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed_column_t), INTENT(inout) :: column(:)
   INTEGER,  INTENT(in)    :: count
   AED_REAL, INTENT(in)    :: Io
   AED_REAL, INTENT(inout) :: extc(:)
   AED_REAL, INTENT(inout) :: par_(:)
   AED_REAL, INTENT(inout) :: h_(:)
!
!LOCAL VARIABLES:
   INTEGER :: i
   AED_REAL :: zz, localext, localshade
!
!-------------------------------------------------------------------------------
!BEGIN
   zz = zero_
   localext = zero_

   CALL BioExtinction(column,count,extc)

   localext = extc(1)
   zz = 0.001 !0.5*h_(1)    !MH: assume top of layer
   par_(1) = 0.45 * Io * EXP( -(localext) * zz )

   IF (count <= 1) RETURN

   DO i = 2, count
      localext = extc(i)

      !zz = zz + 0.5*h_(i)
      zz = h_(i)
      par_(i) = par_(i-1) * EXP( -(localext) * zz )
   ENDDO
END SUBROUTINE Light
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE update_light(column, nlev)
!-------------------------------------------------------------------------------
! Calculate photosynthetically active radiation over entire column based
! on surface radiation, attenuated based on background & biotic extinction
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed_column_t), INTENT(inout) :: column(:)
   INTEGER,INTENT(in)  :: nlev
!
!LOCALS
#if 0
   INTEGER :: i
   AED_REAL :: localext, localext_up
!
!-------------------------------------------------------------------------------
!BEGIN
   localext = zero_; localext_up = zero_

   ! Surface Kd
   CALL aed_light_extinction(column, nlev, localext)

   ! Surface PAR
   par(nlev) = par_fraction * rad(nlev) * EXP( -(lKw+localext)*1e-6*dz(nlev) )

   ! Now set the top of subsequent layers, down to the bottom
   DO i = (nlev-1),1,-1
      localext_up = localext
      CALL aed_light_extinction(column, i, localext)

      par(i) = par(i+1) * EXP( -(lKw + localext_up) * dz(i+1) )

      IF (bioshade_feedback) extc_coef(i) = lKw + localext
   ENDDO
#endif
END SUBROUTINE update_light
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_light
