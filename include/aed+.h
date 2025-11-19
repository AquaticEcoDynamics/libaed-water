/*#############################################################################*
 #                                                                             #
 # aed+.h                                                                      #
 #                                                                             #
 #  Developed by :                                                             #
 #      AquaticEcoDynamics (AED) Group                                         #
 #      School of Agriculture and Environment                                  #
 #      The University of Western Australia                                    #
 #                                                                             #
 #      http://aquatic.science.uwa.edu.au/                                     #
 #                                                                             #
 #  Copyright 2013-2025 - The University of Western Australia                  #
 #                                                                             #
 #   AED is free software: you can redistribute it and/or modify               #
 #   it under the terms of the GNU General Public License as published by      #
 #   the Free Software Foundation, either version 3 of the License, or         #
 #   (at your option) any later version.                                       #
 #                                                                             #
 #   AED is distributed in the hope that it will be useful,                    #
 #   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
 #   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
 #   GNU General Public License for more details.                              #
 #                                                                             #
 #   You should have received a copy of the GNU General Public License         #
 #   along with this program.  If not, see <http://www.gnu.org/licenses/>.     #
 #                                                                             #
 *#############################################################################*/
#ifndef _AED_PLUS_H_
#define _AED_PLUS_H_

#include <aed.h>

#define AED_PLUS_VERSION  "2.3.1"

#ifndef __STD_C__

#define _FLUX_VAR_R_(id)  column(id)%flux_rip

#ifdef __GFORTRAN__
#define PAUSE write(*,'(''paused, type [enter] to continue'')') ; read (*,*)
#endif

#endif

#endif
