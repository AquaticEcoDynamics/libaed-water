!###############################################################################
!#                                                                             #
!# aed_spectra.F90                                                             #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2022-2026 : The University of Western Australia                  #
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
!# Created April 2022 (as aed_maths), moved to libaed-water Sept 2026          #
!#                                                                             #
!###############################################################################

#include "aed.h"

!
MODULE aed_spectra
!-------------------------------------------------------------------------------
! Spectral utilities and optical property tables shared by AED light and
! water-quality modules :
!   - interp                     : linear interpolation between spectral grids
!   - calculate_integral_weights : band integration weights
!   - pure water absorption and scattering (a_w, b_w on lambda_w)
!   - phytoplankton, mineral sediment and seagrass leaf IOPs
!-------------------------------------------------------------------------------

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC interp, calculate_integral_weights
   PUBLIC nlambda_w, a_w, b_w, lambda_w, a_Chlorophytes, a_Coccolithophores, &
          a_Cyanobacteria, a_Diatoms, a_Dinoflagellates, a_he6chl, &
          b_Chlorophytes, b_Coccolithophores, b_Cyanobacteria, b_Diatoms, &
          b_Dinoflagellates, b_he6chl, lambda_Chlorophytes, &
          lambda_Coccolithophores, lambda_Cyanobacteria, lambda_Diatoms, &
          lambda_Dinoflagellates, lambda_he6chl, a_Averagesediment, &
          a_Brownearth, a_Calcareoussand, a_Redclay, a_Yellowclay, &
          b_Averagesediment, b_Brownearth, b_Calcareoussand, b_Redclay, &
          b_Yellowclay, lambda_Averagesediment, lambda_Brownearth, &
          lambda_Calcareoussand, lambda_Redclay, lambda_Yellowclay, &
          a_seagrassleaf, lambda_seagrassleaf
!

   INTERFACE interp
      MODULE PROCEDURE interp_0d
      MODULE PROCEDURE interp_1d
      MODULE PROCEDURE interp_0d_scalar
      MODULE PROCEDURE interp_1d_scalar
   END INTERFACE


!-- Pure water optical properties

integer, parameter :: nlambda_w = 33

AED_REAL, parameter :: a_w(33) = (/ &
      6.112000e-01, 7.620000e-02, 4.610000e-02, 1.820000e-02, 6.300000e-03, 5.100000e-03, 8.300000e-03, 1.190000e-02, 2.150000e-02, 4.070000e-02, &
      5.500000e-02, 8.490000e-02, 1.995000e-01, 2.850000e-01, 3.512000e-01, 4.559000e-01, 6.433000e-01, 1.444900e+00, 2.390000e+00, 3.738200e+00, &
      2.748050e+01, 1.934700e+01, 6.718000e+01, 9.499760e+01, 3.631256e+02, 1.118607e+03, 9.448757e+02, 5.195995e+02, 6.467179e+02, 3.768561e+03, &
      2.628083e+03, 4.376230e+05, 1.338404e+06 &
   /)
AED_REAL, parameter :: b_w(33) = (/ &
      5.670000e-02, 1.870000e-02, 1.350000e-02, 1.000000e-02, 7.600000e-03, 5.800000e-03, 4.500000e-03, 3.600000e-03, 2.900000e-03, 2.300000e-03, &
      1.900000e-03, 1.600000e-03, 1.400000e-03, 1.200000e-03, 9.000000e-04, 7.000000e-04, 7.000000e-04, 6.000000e-04, 4.000000e-04, 2.000000e-04, &
      0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, &
      0.000000e+00, 0.000000e+00, 0.000000e+00 &
   /)
AED_REAL, parameter :: lambda_w(33) = (/ &
      2.500000e+02, 3.250000e+02, 3.500000e+02, 3.750000e+02, 4.000000e+02, 4.250000e+02, 4.500000e+02, 4.750000e+02, 5.000000e+02, 5.250000e+02, &
      5.500000e+02, 5.750000e+02, 6.000000e+02, 6.250000e+02, 6.500000e+02, 6.750000e+02, 7.000000e+02, 7.250000e+02, 7.750000e+02, 8.500000e+02, &
      9.500000e+02, 1.050000e+03, 1.150000e+03, 1.250000e+03, 1.350000e+03, 1.450000e+03, 1.550000e+03, 1.650000e+03, 1.750000e+03, 1.900000e+03, &
      2.200000e+03, 2.900000e+03, 3.700000e+03 &
   /)




!-- Updated phytoplankton IOPs

AED_REAL, parameter :: a_Chlorophytes(19) = (/ &
      2.760000e-03, 2.760000e-02, 2.760000e-02, 2.860000e-02, 3.040000e-02, 3.410000e-02, 3.380000e-02, 3.040000e-02, 2.390000e-02, 1.230000e-02, &
      7.100000e-03, 6.300000e-03, 6.800000e-03, 6.900000e-03, 9.600000e-03, 1.510000e-02, 3.700000e-03, 3.000000e-04, 0.000000e+00 &
   /)
AED_REAL, parameter :: a_Coccolithophores(19) = (/ &
      5.630000e-03, 5.630000e-02, 5.630000e-02, 5.720000e-02, 6.450000e-02, 7.210000e-02, 7.090000e-02, 5.960000e-02, 4.120000e-02, 2.290000e-02, &
      1.230000e-02, 8.500000e-03, 8.200000e-03, 8.100000e-03, 8.300000e-03, 1.800000e-02, 4.900000e-03, 8.000000e-04, 0.000000e+00 &
   /)
AED_REAL, parameter :: a_Cyanobacteria(19) = (/ &
      5.780000e-03, 5.780000e-02, 5.780000e-02, 6.010000e-02, 6.510000e-02, 7.670000e-02, 7.080000e-02, 5.230000e-02, 4.780000e-02, 4.080000e-02, &
      4.400000e-02, 4.510000e-02, 2.940000e-02, 3.120000e-02, 2.980000e-02, 3.730000e-02, 2.200000e-02, 1.320000e-02, 0.000000e+00 &
   /)
AED_REAL, parameter :: a_Diatoms(19) = (/ &
      3.090000e-03, 3.090000e-02, 3.090000e-02, 3.200000e-02, 3.400000e-02, 3.840000e-02, 3.790000e-02, 3.240000e-02, 2.720000e-02, 2.060000e-02, &
      1.400000e-02, 8.700000e-03, 7.100000e-03, 9.100000e-03, 1.040000e-02, 2.230000e-02, 4.200000e-03, 3.000000e-04, 0.000000e+00 &
   /)
AED_REAL, parameter :: a_Dinoflagellates(19) = (/ &
      2.690000e-03, 2.690000e-02, 2.690000e-02, 2.800000e-02, 2.990000e-02, 3.370000e-02, 3.450000e-02, 3.160000e-02, 2.700000e-02, 2.110000e-02, &
      1.670000e-02, 1.020000e-02, 9.000000e-03, 8.900000e-03, 9.300000e-03, 1.850000e-02, 6.500000e-03, 3.000000e-04, 0.000000e+00 &
   /)
AED_REAL, parameter :: a_he6chl(19) = (/ &
      2.800000e-03, 2.800000e-02, 2.800000e-02, 3.000000e-02, 3.635000e-02, 4.695000e-02, 4.720000e-02, 4.195000e-02, 3.340000e-02, 2.520000e-02, &
      1.785000e-02, 1.340000e-02, 1.180000e-02, 1.495000e-02, 1.780000e-02, 2.720000e-02, 1.075000e-02, 1.500000e-03, 0.000000e+00 &
   /)
AED_REAL, parameter :: b_Chlorophytes(19) = (/ &
      1.199000e-01, 1.199000e-01, 9.230000e-02, 9.130000e-02, 8.950000e-02, 8.590000e-02, 8.610000e-02, 8.950000e-02, 9.600000e-02, 1.076000e-01, &
      1.128000e-01, 1.137000e-01, 1.131000e-01, 1.130000e-01, 1.103000e-01, 1.048000e-01, 1.163000e-01, 1.196000e-01, 1.196000e-01 &
   /)
AED_REAL, parameter :: b_Coccolithophores(19) = (/ &
      8.391000e-01, 7.732000e-01, 6.882000e-01, 6.712000e-01, 6.416000e-01, 6.115000e-01, 5.904000e-01, 5.793000e-01, 5.753000e-01, 5.713000e-01, &
      5.595000e-01, 5.409000e-01, 5.189000e-01, 4.965000e-01, 4.740000e-01, 4.420000e-01, 4.326000e-01, 4.144000e-01, 4.144000e-01 &
   /)
AED_REAL, parameter :: b_Cyanobacteria(19) = (/ &
      5.092000e-01, 4.502000e-01, 3.668000e-01, 3.501000e-01, 3.251000e-01, 2.935000e-01, 2.794000e-01, 2.779000e-01, 2.624000e-01, 2.494000e-01, &
      2.262000e-01, 2.051000e-01, 2.008000e-01, 1.790000e-01, 1.604000e-01, 1.328000e-01, 1.282000e-01, 1.170000e-01, 1.170000e-01 &
   /)
AED_REAL, parameter :: b_Diatoms(19) = (/ &
      3.972000e-01, 3.655000e-01, 3.208000e-01, 3.120000e-01, 2.992000e-01, 2.841000e-01, 2.738000e-01, 2.686000e-01, 2.631000e-01, 2.589000e-01, &
      2.547000e-01, 2.494000e-01, 2.402000e-01, 2.275000e-01, 2.154000e-01, 1.927000e-01, 2.001000e-01, 1.932000e-01, 1.932000e-01 &
   /)
AED_REAL, parameter :: b_Dinoflagellates(19) = (/ &
      4.620000e-02, 4.620000e-02, 1.930000e-02, 1.820000e-02, 1.640000e-02, 1.250000e-02, 1.170000e-02, 1.460000e-02, 1.920000e-02, 2.510000e-02, &
      2.950000e-02, 3.600000e-02, 3.720000e-02, 3.730000e-02, 3.690000e-02, 2.770000e-02, 3.970000e-02, 4.590000e-02, 4.590000e-02 &
   /)
AED_REAL, parameter :: b_he6chl(19) = (/ &
      6.600000e-01, 5.076923e-01, 4.714286e-01, 4.400000e-01, 4.125000e-01, 3.882353e-01, 3.666667e-01, 3.473684e-01, 3.300000e-01, 3.142857e-01, &
      3.000000e-01, 2.869565e-01, 2.750000e-01, 2.640000e-01, 2.538462e-01, 2.444444e-01, 2.357143e-01, 2.275862e-01, 2.200000e-01 &
   /)
AED_REAL, parameter :: lambda_Chlorophytes(19) = (/ &
      2.500000e+02, 3.250000e+02, 3.500000e+02, 3.750000e+02, 4.000000e+02, 4.250000e+02, 4.500000e+02, 4.750000e+02, 5.000000e+02, 5.250000e+02, &
      5.500000e+02, 5.750000e+02, 6.000000e+02, 6.250000e+02, 6.500000e+02, 6.750000e+02, 7.000000e+02, 7.250000e+02, 7.500000e+02 &
   /)
AED_REAL, parameter :: lambda_Coccolithophores(19) = (/ &
      2.500000e+02, 3.250000e+02, 3.500000e+02, 3.750000e+02, 4.000000e+02, 4.250000e+02, 4.500000e+02, 4.750000e+02, 5.000000e+02, 5.250000e+02, &
      5.500000e+02, 5.750000e+02, 6.000000e+02, 6.250000e+02, 6.500000e+02, 6.750000e+02, 7.000000e+02, 7.250000e+02, 7.500000e+02 &
   /)
AED_REAL, parameter :: lambda_Cyanobacteria(19) = (/ &
      2.500000e+02, 3.250000e+02, 3.500000e+02, 3.750000e+02, 4.000000e+02, 4.250000e+02, 4.500000e+02, 4.750000e+02, 5.000000e+02, 5.250000e+02, &
      5.500000e+02, 5.750000e+02, 6.000000e+02, 6.250000e+02, 6.500000e+02, 6.750000e+02, 7.000000e+02, 7.250000e+02, 7.500000e+02 &
   /)
AED_REAL, parameter :: lambda_Diatoms(19) = (/ &
      2.500000e+02, 3.250000e+02, 3.500000e+02, 3.750000e+02, 4.000000e+02, 4.250000e+02, 4.500000e+02, 4.750000e+02, 5.000000e+02, 5.250000e+02, &
      5.500000e+02, 5.750000e+02, 6.000000e+02, 6.250000e+02, 6.500000e+02, 6.750000e+02, 7.000000e+02, 7.250000e+02, 7.500000e+02 &
   /)
AED_REAL, parameter :: lambda_Dinoflagellates(19) = (/ &
      2.500000e+02, 3.250000e+02, 3.500000e+02, 3.750000e+02, 4.000000e+02, 4.250000e+02, 4.500000e+02, 4.750000e+02, 5.000000e+02, 5.250000e+02, &
      5.500000e+02, 5.750000e+02, 6.000000e+02, 6.250000e+02, 6.500000e+02, 6.750000e+02, 7.000000e+02, 7.250000e+02, 7.500000e+02 &
   /)
AED_REAL, parameter :: lambda_he6chl(19) = (/ &
      2.500000e+02, 3.250000e+02, 3.500000e+02, 3.750000e+02, 4.000000e+02, 4.250000e+02, 4.500000e+02, 4.750000e+02, 5.000000e+02, 5.250000e+02, &
      5.500000e+02, 5.750000e+02, 6.000000e+02, 6.250000e+02, 6.500000e+02, 6.750000e+02, 7.000000e+02, 7.250000e+02, 7.500000e+02 &
   /)

!-- New HE6 mineral sediment IOPs

AED_REAL, parameter :: a_Averagesediment(19) = (/ &
      1.511500e-01, 1.362900e-01, 1.187500e-01, 1.024100e-01, 8.750000e-02, 7.325000e-02, 5.775000e-02, 4.925000e-02, 4.150000e-02, 3.350000e-02, &
      2.625000e-02, 1.975000e-02, 1.650000e-02, 1.475000e-02, 1.375000e-02, 1.250000e-02, 1.162000e-02, 1.111000e-02, 1.078000e-02 &
   /)
AED_REAL, parameter :: a_Brownearth(19) = (/ &
      2.100000e-01, 1.969700e-01, 1.800000e-01, 1.579900e-01, 1.370000e-01, 1.200000e-01, 9.800000e-02, 8.500000e-02, 7.300000e-02, 6.000000e-02, &
      4.700000e-02, 3.400000e-02, 2.700000e-02, 2.300000e-02, 2.100000e-02, 1.900000e-02, 1.800000e-02, 1.725000e-02, 1.649000e-02 &
   /)
AED_REAL, parameter :: a_Calcareoussand(19) = (/ &
      7.000000e-02, 5.924000e-02, 5.000000e-02, 4.244000e-02, 3.600000e-02, 3.000000e-02, 2.400000e-02, 2.000000e-02, 1.800000e-02, 1.600000e-02, &
      1.400000e-02, 1.300000e-02, 1.200000e-02, 1.200000e-02, 1.200000e-02, 1.100000e-02, 1.000000e-02, 9.490000e-03, 9.280000e-03 &
   /)
AED_REAL, parameter :: a_Redclay(19) = (/ &
      2.200000e-01, 1.844900e-01, 1.550000e-01, 1.323900e-01, 1.130000e-01, 9.300000e-02, 7.300000e-02, 6.100000e-02, 5.000000e-02, 4.000000e-02, &
      3.000000e-02, 2.000000e-02, 1.600000e-02, 1.400000e-02, 1.300000e-02, 1.200000e-02, 1.100000e-02, 1.044000e-02, 1.020000e-02 &
   /)
AED_REAL, parameter :: a_Yellowclay(19) = (/ &
      1.200000e-01, 1.044500e-01, 9.000000e-02, 7.684000e-02, 6.400000e-02, 5.000000e-02, 3.600000e-02, 3.100000e-02, 2.500000e-02, 1.800000e-02, &
      1.400000e-02, 1.200000e-02, 1.100000e-02, 1.000000e-02, 9.000000e-03, 8.000000e-03, 7.500000e-03, 7.270000e-03, 7.150000e-03 &
   /)
AED_REAL, parameter :: b_Averagesediment(19) = (/ &
      7.888200e-01, 7.960000e-01, 8.075000e-01, 8.212800e-01, 8.435000e-01, 8.762500e-01, 8.945000e-01, 8.867500e-01, 8.655000e-01, 8.485000e-01, &
      8.297500e-01, 8.045000e-01, 7.742500e-01, 7.357500e-01, 6.962500e-01, 6.537500e-01, 6.290000e-01, 5.946400e-01, 5.695400e-01 &
   /)
AED_REAL, parameter :: b_Brownearth(19) = (/ &
      6.000000e-01, 6.113600e-01, 6.300000e-01, 6.545300e-01, 6.970000e-01, 7.580000e-01, 7.880000e-01, 8.080000e-01, 8.180000e-01, 8.230000e-01, &
      8.330000e-01, 8.280000e-01, 7.980000e-01, 7.580000e-01, 7.070000e-01, 6.570000e-01, 6.160000e-01, 5.825600e-01, 5.541400e-01 &
   /)
AED_REAL, parameter :: b_Calcareoussand(19) = (/ &
      1.050000e+00, 1.051390e+00, 1.050000e+00, 1.044930e+00, 1.041000e+00, 1.041000e+00, 1.033000e+00, 1.008000e+00, 9.750000e-01, 9.430000e-01, &
      9.020000e-01, 8.610000e-01, 8.440000e-01, 7.990000e-01, 7.620000e-01, 7.380000e-01, 7.130000e-01, 6.720000e-01, 6.480000e-01 &
   /)
AED_REAL, parameter :: b_Redclay(19) = (/ &
      8.000000e-01, 8.049700e-01, 8.100000e-01, 8.149100e-01, 8.210000e-01, 8.360000e-01, 8.660000e-01, 8.510000e-01, 8.210000e-01, 7.910000e-01, &
      7.690000e-01, 7.460000e-01, 7.160000e-01, 6.790000e-01, 6.420000e-01, 5.900000e-01, 5.670000e-01, 5.370000e-01, 5.220000e-01 &
   /)
AED_REAL, parameter :: b_Yellowclay(19) = (/ &
      7.000000e-01, 7.162600e-01, 7.400000e-01, 7.707600e-01, 8.150000e-01, 8.700000e-01, 8.910000e-01, 8.800000e-01, 8.480000e-01, 8.370000e-01, &
      8.150000e-01, 7.830000e-01, 7.390000e-01, 7.070000e-01, 6.740000e-01, 6.300000e-01, 6.200000e-01, 5.870000e-01, 5.540000e-01 &
   /)
AED_REAL, parameter :: lambda_Averagesediment(19) = (/ &
      2.500000e+02, 3.250000e+02, 3.500000e+02, 3.750000e+02, 4.000000e+02, 4.250000e+02, 4.500000e+02, 4.750000e+02, 5.000000e+02, 5.250000e+02, &
      5.500000e+02, 5.750000e+02, 6.000000e+02, 6.250000e+02, 6.500000e+02, 6.750000e+02, 7.000000e+02, 7.250000e+02, 7.500000e+02 &
   /)
AED_REAL, parameter :: lambda_Brownearth(19) = (/ &
      2.500000e+02, 3.250000e+02, 3.500000e+02, 3.750000e+02, 4.000000e+02, 4.250000e+02, 4.500000e+02, 4.750000e+02, 5.000000e+02, 5.250000e+02, &
      5.500000e+02, 5.750000e+02, 6.000000e+02, 6.250000e+02, 6.500000e+02, 6.750000e+02, 7.000000e+02, 7.250000e+02, 7.500000e+02 &
   /)
AED_REAL, parameter :: lambda_Calcareoussand(19) = (/ &
      2.500000e+02, 3.250000e+02, 3.500000e+02, 3.750000e+02, 4.000000e+02, 4.250000e+02, 4.500000e+02, 4.750000e+02, 5.000000e+02, 5.250000e+02, &
      5.500000e+02, 5.750000e+02, 6.000000e+02, 6.250000e+02, 6.500000e+02, 6.750000e+02, 7.000000e+02, 7.250000e+02, 7.500000e+02 &
   /)
AED_REAL, parameter :: lambda_Redclay(19) = (/ &
      2.500000e+02, 3.250000e+02, 3.500000e+02, 3.750000e+02, 4.000000e+02, 4.250000e+02, 4.500000e+02, 4.750000e+02, 5.000000e+02, 5.250000e+02, &
      5.500000e+02, 5.750000e+02, 6.000000e+02, 6.250000e+02, 6.500000e+02, 6.750000e+02, 7.000000e+02, 7.250000e+02, 7.500000e+02 &
   /)
AED_REAL, parameter :: lambda_Yellowclay(19) = (/ &
      2.500000e+02, 3.250000e+02, 3.500000e+02, 3.750000e+02, 4.000000e+02, 4.250000e+02, 4.500000e+02, 4.750000e+02, 5.000000e+02, 5.250000e+02, &
      5.500000e+02, 5.750000e+02, 6.000000e+02, 6.250000e+02, 6.500000e+02, 6.750000e+02, 7.000000e+02, 7.250000e+02, 7.500000e+02 &
   /)

! leaf light absorbance at selected wavelengths; from Baird et al (2016)
AED_REAL, parameter ::       a_seagrassleaf(8) = (/ 0.72, 0.72, 0.68, 0.38, 0.38, 0.60, 0.04, 0.0 /) 
AED_REAL, parameter :: lambda_seagrassleaf(8) = (/ 300., 390., 500., 530., 640., 680., 705., 800. /) 


CONTAINS

!###############################################################################
SUBROUTINE interp_0d(nsource,x,y,ntarget,targetx,targety)
!-------------------------------------------------------------------------------
! 0D interpolation, extrapolates beyond boundaries
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)                      :: nsource, ntarget
   AED_REAL,DIMENSION(nsource),INTENT(in)  :: x, y
   AED_REAL,DIMENSION(ntarget),INTENT(in)  :: targetx
   AED_REAL,DIMENSION(ntarget),INTENT(out) :: targety
!
!LOCALS
   INTEGER  :: i, j
   AED_REAL :: frac

!-------------------------------------------------------------------------------
!BEGIN
   i = 1
   DO j = 1,ntarget
      DO while (i+1<nsource)
         IF (x(i+1)>=targetx(j)) exit
         i = i+1
      ENDDO
      frac = (targetx(j)-x(i))/(x(i+1)-x(i))
      targety(j) = y(i) + frac*(y(i+1)-y(i))
   ENDDO
END SUBROUTINE interp_0d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE interp_0d_scalar(nsource,x,y,targetx,targety)
!-------------------------------------------------------------------------------
! 0D interpolation, extrapolates beyond boundaries
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)                     :: nsource
   AED_REAL,DIMENSION(nsource),INTENT(in) :: x, y
   AED_REAL,INTENT(in)                    :: targetx
   AED_REAL,INTENT(out)                   :: targety
!
!LOCALS
   INTEGER  :: i
   AED_REAL :: frac

!-------------------------------------------------------------------------------
!BEGIN
   i = 1
   DO while (i+1<nsource)
      IF (x(i+1)>=targetx) exit
      i = i+1
   ENDDO
   frac = (targetx-x(i))/(x(i+1)-x(i))
   targety = y(i) + frac*(y(i+1)-y(i))
END SUBROUTINE interp_0d_scalar
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE interp_1d(m,nsource,x,y,ntarget,targetx,targety)
!-------------------------------------------------------------------------------
! 1D interpolation, extrapolates beyond boundaries
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)                        :: nsource, ntarget, m
   AED_REAL,DIMENSION(nsource),INTENT(in)    :: x
   AED_REAL,DIMENSION(m,nsource),INTENT(in)  :: y
   AED_REAL,DIMENSION(ntarget),INTENT(in)    :: targetx
   AED_REAL,DIMENSION(m,ntarget),INTENT(out) :: targety
!
!LOCALS
   INTEGER  :: i, j
   AED_REAL :: frac

!-------------------------------------------------------------------------------
!BEGIN
   i = 1
   DO j = 1,ntarget
      DO while (i+1<nsource)
         IF (x(i+1)>=targetx(j)) exit
         i = i+1
      ENDDO
      frac = (targetx(j)-x(i))/(x(i+1)-x(i))
      targety(:,j) = y(:,i) + frac*(y(:,i+1)-y(:,i))
   ENDDO
END SUBROUTINE interp_1d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE interp_1d_scalar(m,nsource,x,y,targetx,targety)
!-------------------------------------------------------------------------------
! 1D interpolation, extrapolates beyond boundaries
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)                       :: nsource, m
   AED_REAL,DIMENSION(nsource),INTENT(in)   :: x
   AED_REAL,DIMENSION(m,nsource),INTENT(in) :: y
   AED_REAL,INTENT(in)                      :: targetx
   AED_REAL,DIMENSION(m),INTENT(out)        :: targety
!
!LOCALS
   INTEGER  :: i
   AED_REAL :: frac

!-------------------------------------------------------------------------------
!BEGIN
   i = 1
   DO while (i+1<nsource)
      IF (x(i+1)>=targetx) exit
      i = i+1
   ENDDO
   frac = (targetx-x(i))/(x(i+1)-x(i))
   targety(:) = y(:,i) + frac*(y(:,i+1)-y(:,i))
END SUBROUTINE interp_1d_scalar
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE calculate_integral_weights(xl, xr, n, x, w)
   !-------------------------------------------------------------------------------
   !-------------------------------------------------------------------------------
   !ARGUMENTS
      INTEGER, INTENT(in)  :: n
      AED_REAL,INTENT(in)  :: x(n), xl, xr
      AED_REAL,INTENT(out) :: w(n)
   !
   !LOCALS
      INTEGER  :: i
      AED_REAL :: deltax, f
   !-------------------------------------------------------------------------------
   !BEGIN
      w = 0.

      IF (x(1) > xr) THEN
         ! All points to right of desired range.
         w(1) = xr - xl
         return
      ELSEIF (x(n) < xl) THEN
         ! All points to left of desired range.
         w(n) =  xr - xl
         return
      ENDIF

      DO i = 2, n
         deltax = x(i) - x(i-1)
         IF (x(i) >= xl .and. x(i) <= xr) THEN
            ! Right-hand point in desired range.
            IF (x(i-1) >= xl) THEN
               ! Whole segment in range
               ! Integral: 0.5*(y1+y2)*(x2-x1)
               w(i-1:i) = w(i-1:i) + 0.5 * deltax
            ELSE
               ! Segment crosses left boundary
               ! y at boundary: yb = y1 + (y2-y1)/(x2-x1)*(xl-x1) = y1*(1-(xl-x1)/(x2-x1)) + y2*(xl-x1)/(x2-x1)
               ! integral = 0.5*(yb+y2)*(x2-xl)
               ! yb+y2 = y1*(1-(xl-x1)/(x2-x1)) + y2*(1+(xl-x1)/(x2-x1))
               f = (xl - x(i-1)) / deltax
               w(i-1) = w(i-1) + (1.0 - f) * (x(i) - xl) * 0.5
               w(i  ) = w(i  ) + (1.0 + f) * (x(i) - xl) * 0.5
            ENDIF
         ELSEIF (x(i) > xr .and. x(i-1) < xr) THEN
            ! Right-hand point beyond desired range, left-hand point before right range boundary.
            IF (x(i-1) >= xl) THEN
               ! Segment crosses right boundary
               ! y at boundary: yb = y1 + (y2-y1)/(x2-x1)*(xr-x1) = y1*(1-(xr-x1)/(x2-x1)) + y2*(xr-x1)/(x2-x1)
               ! integral = 0.5*(y1+yb)*(xr-x1)
               ! y1+yb = y1*(2-(xr-x1)/(x2-x1)) + y2*(xr-x1)/(x2-x1)
               f = (xr - x(i-1)) / deltax
               w(i-1) = w(i-1) + (2.0 - f) * (xr - x(i-1)) * 0.5
               w(i  ) = w(i  ) + (         f) * (xr - x(i-1)) * 0.5
            ELSE
               ! Segment crosses both boundaries
               ! y at centre: yc = y1 + (y2-y1)/(x2-x1)*(0.5*(xl+xr)-x1) = y1*(1-(0.5*(xl+xr)-x1)/(x2-x1)) + y2*(0.5*(xl+xr)-x1)/(x2-x1)
               ! integral: (xr-xl)*yc
               f = (0.5*(xl+xr)-x(i-1))/deltax
               w(i-1) = w(i-1) + (1.-f)*(xr-xl)
               w(i  ) = w(i  ) + (      f)*(xr-xl)
            ENDIF
         ENDIF
      ENDDO
      IF (x(1)>xl) w(1) = w(1) + (x(1)-xl)
      IF (x(n)<xr) w(n) = w(n) + (xr-x(n))
   END SUBROUTINE calculate_integral_weights
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_spectra
