!###############################################################################
!#                                                                             #
!# aed+.h                                                                      #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2013 - 2020 -  The University of Western Australia               #
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
!###############################################################################
#ifndef _AED_PLUS_H_
#define _AED_PLUS_H_

#include <aed.h>

#define AED_PLUS_VERSION  "2.0.0"

!# for aed_geochemistry
#define MAX_GC_COMPONENTS 20
#define MAX_GC_MINERALS   20

!# for aed_vegetation
#define MAX_VEG_TYPES   256

!# for aed_macrophytes
#define MAX_ZONES       256

!# for aed_ass
#define MAX_ASS_PARAMS  20

#define REACTION_START_CH  '['
#define REACTION_END_CH    ']'
#define PLUS               '+'
#define MINUS              '-'
#define EQUALS             '='
#define L_PAREN            '('
#define R_PAREN            ')'

#define _FLUX_VAR_R_(id)  column(id)%flux_rip

#ifdef __GFORTRAN__
#define PAUSE write(*,'(''paused, type [enter] to continue'')') ; read (*,*)
#endif

#endif
