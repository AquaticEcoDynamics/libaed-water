!###############################################################################
!#                                                                             #
!# aed_pibm_utils.F90                                                          #
!#                                                                             #
!#  Developed by :                                                             #
!#      Bingzhang Chen and Iria Sala                                           #
!#      University of Strathclyde, United Kingdom                              #
!#                                                                             #
!#      https://github.com/BingzhangChen/PIBM                                  #
!#                                                                             #
!#  Copyright (c) 2025 Bingzhang Chen and Iria Sala                            #
!#  University of Strathclyde, United Kingdom                                  #
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
!# Created Jul 2018                                                            #
!#                                                                             #
!###############################################################################

!###############################################################################
MODULE PARAMS
implicit NONE

! Parameters
real :: mu0     =  5.00  !Maximal growth rate normalized to 15 C (d-1)
real :: aI0     =  1.20  !Chl-specific Slope of P-I curve ((W m-2)-1 (gChl molC)-1 d-1)
real :: KN      =  0.50  !Half saturation constant of phytoplankton growth on N
real :: gmax    =  2.00  !Maximal zooplankton grazing rate
real :: Kp      =  0.15  !Half saturation constant of zooplankton grazing
real :: mz      =  0.002  !Zooplankton mortality grazing coefficient
real :: GGE     =  0.30  !Zooplankton Gross Growth Efficiency [dimensionless]
real :: unass   =  0.24  !Fraction of unassimilated food ingested by zooplankton [nd]
real :: RdN     =  0.1  !Rate of detritus converting to DIN
real :: wDET    =  0.5 ! Sinking rate of detritus

! Standard deviation of log zooplankton feeding preference
real :: SDZoo   =  0.5 

! Acclimation parameter
real, parameter :: nx = 1.d0

! Activation energy for phytoplankton growth
real, parameter :: Ep  = 0.32d0

! Activation energy for zooplankton grazing
real, parameter :: Ez  = 0.65d0

! Maixmal Chl:N ratio (g:mol) for phytoplankton
real, parameter :: thetaNmax = 3d0

!the value of rhochl before last sunset
real            :: rhochl_L  = 0d0

!QNmin and QNmax are allometric functions of Vol (Ward et al. 2012) [mol N: mol C]:
real, parameter   :: QNmin_a = 0.07d0     !Normalization constant for QNmin [molN:molC]
real, parameter   :: QNmin_b = -0.17d0    !Allometric exponent for QNmin
real, parameter   :: QNmax_a = 0.166d0    !Normalization constant for QNmax [molN:molC]
real, parameter   :: QNmax_b = 0.0d0      !Allometric exponent for QNmax which is independent from size

!Kn is an allometric function of Vol (Cdiv) (Edwards et al. 2012) [uM]:
real, parameter   :: KN_a   = 0.14d0      !Normalization constant for KN
real, parameter   :: KN_b   = 0.33d0      !Allometric exponent for KN

real, parameter   :: pi= 3.1415926535

integer, parameter :: NTrait = 3
real :: nu(NTrait) = [1d-12, 1d-12, 1d-12] !Probability per generation per cell
real :: sigma(NTrait) = [0.1, 0.1, 0.1]    !Standard deviation of mutation of the three traits

end module params

!###############################################################################
!
Module Trait_functions
!-------------------------------------------------------------------------------
! Trait_functions --- utility functions for the PIBM GMK98_ToptSizeLight model
!
!-------------------------------------------------------------------------------
use params

!This module provides several functions calculating phytoplankton physiological rates as a function of environmental conditions (e.g., temperature) and traits
implicit none

private

public :: TEMPBOL, temp_Topt, palatability, PHY_C2Vol, PHY_ESD2C, Ainf
public :: Pmax_size, respiration

CONTAINS

!Function converting phytoplankton ESD (micron) to Carbon (unit: pmol C per cell) with the parameters obtained from Maranon et al. (2013)
real function PHY_ESD2C(ESD) result(y)
implicit none
real, intent(in)  :: ESD !Phytoplankton ESD (micron)

real, parameter :: a = -0.69
real, parameter :: b = 0.88

real :: Vol = 0d0

!Calculate volume from ESD
Vol = pi/6d0 * ESD**3

!the parameters of a and b are derived from pg C cell-1 to volume; so we need to
!convert the carbon unit from pmol C cell-1 to pg C cell-1

!Calculate carbon (pmol C per cell) from volume
y = 10.d0**a * Vol**b / 12.d0

return
end function PHY_ESD2C

!Function converting phytoplankton carbon to volume (unit: micron^3) with the parameters obtained from Maranon et al. (2013)
pure real function PHY_C2Vol(p_C) result(y)
implicit none
real, intent(in)  :: p_C !Phytoplankton carbon (pmol C/cell)

real, parameter :: a = -0.69
real, parameter :: b = 0.88

!the parameters of a and b are derived from pg C cell-1 to volume; so we need to
!convert the carbon unit from pmol C cell-1 to pg C cell-1

y = (12d0 * p_C/10.d0**a)**(1d0/b)
return
end function PHY_C2Vol

!Function calculating the prey palatability based on Ward et al. L&O 2012 (Eq. A21) and Banas Ecol. Mod. (2011) Table 2
real function palatability(Vpred, Vprey, SDpref_Z) result(y)
implicit none

!Predator Volume
real, intent(in)  :: Vpred

!Prey Volume
real, intent(in)  :: Vprey

! Standard deviation of log zooplankton feeding preference
real, intent(in)  :: SDpref_Z

!The actual predator:prey volume ratio
real  :: R_real = 0d0

!Optimal prey Volume of the predator
real :: Vprey_opt = 1d3 

real :: Xpred, Xprey !ESD of predator and prey
real :: Xprey_opt = 0. !Optimal prey ESD

!Optimal predator:prey volume ratio
real :: R_opt = 1d3 

!Maximal predator:prey volume ratio following copepods (Hansen et al. 1994)
real, parameter :: R_opt_max = 18.**3

real :: cff = 0d0

!First calculate prey and predator ESD (micron) from volume
Xpred = (6d0*Vpred/pi)**(.33333333333333)
Xprey = (6d0*Vprey/pi)**(.33333333333333)

!Then calculate optimal prey ESD (micron)
Xprey_opt = 0.65 * Xpred**0.56

!Convert prey ESD to volume
Vprey_opt = pi/6d0 * Xprey_opt**3

R_opt  = Vpred/Vprey_opt

R_opt  = min(R_opt, R_opt_max)
R_real = Vpred/Vprey

cff = log(R_real/R_opt)
cff = cff**2/(2.d0 * SDpref_Z**2)

!To avoid underflow
if (cff > 5d2) then
   y = 0d0
else
   y = exp(-cff)
endif
return
end function palatability

pure real function Pmax_size(ESD, Pmax0) result(y)
implicit none
real, intent(in) :: ESD

!Pmax0: Maximal photosynthesis rate (d-1)
real, intent(in) :: Pmax0

!Constant in Eqn. 14 of Wirtz (2011) (a' = (rho*/rho)^.333*a) in which a = 0.34, rho* = 0.25, and rho = 0.5
real, parameter :: a_p = 0.27

!End of declaration

y = Pmax0/(1.d0 + a_p *ESD)

return
End function Pmax_size

pure real function respiration(ESD, r_s) result(y)
implicit none
real, intent(in):: ESD
real, intent(in):: r_s !respiration rate (d-1) at V_s

real, parameter :: b_rho = 0.d0 !Size scaling of C density

!Cell volume when rho_dia = rho_green
real, parameter :: V_s = 8.d0 !micron^3

real :: V, ESD_s

V = pi/6.d0*ESD**3 
ESD_s = (6.d0*V_s/pi)**0.333333

y = r_s*ESD_s/ESD * (V/V_s)**b_rho

return
end function 

pure real function TEMPBOL(Ea,tC)
implicit none
!DESCRIPTION:
!The temperature dependence of plankton rates are fomulated according to the Arrhenuis equation. 
! tC: in situ temperature
! Tr: reference temperature
!
!INPUT PARAMETERS:
real, intent (in) :: Ea, tC

! boltzman constant constant [ eV /K ]
real, parameter   :: kb = 8.62d-5, Tr = 15D0

TEMPBOL = exp(-(Ea/kb)*(1D0/(273.15 + tC)-1D0/(273.15 + Tr)))
return 
end function TEMPBOL

REAL function temp_Topt(tC, mumax0, Topt_) result(y)
!Function of a rate depending on Temperature and optimal temperature (Topt_) modified from Chen Ecol. Mod. (2022)
IMPLICIT NONE
real, intent(in) :: mumax0    !Maximal rate normalized to an optimal temperature of 15 ºC
real, intent(in) :: tC         !Environmental temperature in ºC
real, intent(in) :: Topt_   !Optimal temperature in ºC

real, parameter   :: Ea0   = 0.98  
real, parameter   :: Ed0   = 2.3
real, parameter   :: Ei    = 0.22  
real, parameter   :: beta  =-0.2  !Exponent for Ea0
real, parameter   :: phi   = 0.27  !Exponent for Ed

real :: Ed, Ea, mumax

mumax = alloscale(Topt_, mumax0,  Ei) 
Ea    = alloscale(Topt_, Ea0,  beta) 
Ed    = alloscale(Topt_, Ed0,  phi) 
y     = JOHNSON(tC, mumax, Ea, Ed, Topt_)
return
END function temp_Topt

REAL FUNCTION JOHNSON(tC, mumax, Ea, Ed, Topt_) RESULT(y)
!Temperature function following Dell et al. PNAS (2011) and Chen & Laws L&O (2017)
IMPLICIT NONE
!Both tC and Topt_ are in ºC
real,   intent(in)     :: tC, mumax, Ea, Ed, Topt_
real,   parameter   :: kb   = 8.62D-5
real,   parameter   :: T0   = 273.15D0
real,   parameter   :: Tref = 15D0
real                         :: Eh, x, theta, b

if (Ed .le. 0d0) stop "Ed must be greater than zero!"
Eh = Ed+Ea
x    = TK(TC)
theta = TK(Topt_)
b = x - theta
y = mumax*(Ea/Ed + 1.d0) * exp(Ea*b)/(1.D0+Ea/ED*exp(Eh*b))   
return
END FUNCTION JOHNSON

PURE REAL FUNCTION TK(TC)
IMPLICIT NONE
!DESCRIPTION:
!The temperature dependence of plankton rates are fomulated according to the Arrhenuis equation. 
! tC: in situ temperature
! Tr: reference temperature
!
!INPUT PARAMETERS:
REAL, INTENT (IN) :: TC
! boltzman constant constant [ eV /K ]
REAL, PARAMETER   :: kb = 8.62d-5, Tr = 15.0

TK = -(1./kb)*(1./(273.15 + tC) - 1./(273.15 + Tr))
return 
END FUNCTION TK

PURE REAL FUNCTION alloscale(Topt_, mu0p, alpha)
IMPLICIT NONE
real, intent(in) :: Topt_     !Topt in ºC
real, intent(in) :: mu0p  !Normalized growth rate
real, intent(in) :: alpha    !Exponent of thermal traits normalized to z
alloscale =  mu0p * exp(TK(Topt_) * alpha) 
END FUNCTION alloscale

!------------------------------------------------------------------------------------------------
!Function to estimate photoinhibition following Nikolaou et al. (2016) (J. Theor. Biol.), and
!Han (2001) (J. Theor. Biol.)
!Assuming that acclimation to photoinhibition is at the time-scale of ms.
!------------------------------------------------------------------------------------------------
PURE REAL FUNCTION Ainf(PAR_, alpha_, QN_, QNmin_, QNmax_, theta_)

implicit none

!Declaration of variables:
real, intent(in) :: PAR_           !Irradiance [W m-2]
real, intent(in) :: alpha_         !Slope of the P-I curve [Unit: molC/gChl m2/uE]
real, intent(in) :: QN_            !N:C ratio of the phyto. super-individual [mol N mol C-1]
real, intent(in) :: QNmin_         !Minimal N:C ratio
real, intent(in) :: QNmax_         !Maximal N:C ratio
real, intent(in) :: theta_         !Chl:C ratio [mg Chl mmol C]

real, parameter  :: Tau   = 5.5d-3 !Turnover time of the electron transport chain [s]
real, parameter  :: Beta  = 0.492  !Pre-exponential factor of effective cross-section eq [m2 uE-1 (g Chl)^(1/Kappa) (g C)^(-1/Kappa)]
real, parameter  :: Kappa = 0.469  !Exponent of effective cross-section equation [nd]
real, parameter  :: Kd    = 5d-6   !Damage constant of a photosynthetic unit [nd]

real, parameter  :: WtouE = 4.57   !Constant to convert PAR units from [Wm-2] to [uE m-2 s-1]
real, parameter  :: a_ = 2d-5      !The constant a in the equation relating Kr and alphaChl
real, parameter  :: b_ = 5d-7      !The constant b in the equation relating Kr and alphaChl
real, parameter  :: v_ = -6.64     !The constant v in the equation relating Kr and alphaChl
real             :: Kr0            !Repair constant of a photosynthetic unit [s-1] under nutrient saturated conditions which depends on alpha to impose a tradeoff
real             :: Kr             !Nutrient dependent Repair constant of a photosynthetic unit [s-1]
real             :: K              !Ratio of damage to repair constants [s]
real             :: Sigma          !Effective cross-section of the PSU [m2 uE-1]

real             :: thetaA         !Chl:C ratio (g Chl g C-1) to be consistent with Han (2001)
real             :: PARWm2
real             :: Lno3           !Nutrient limitation index
real             :: alpha_new      !alphaChl with the correct unit
!End of declaration

!PAR, unit conversion:
PARWm2 = PAR_ * WtouE    ![W m-2] to [uE m-2]

!Convert the unit of alpha from molC/gChl (W m-2)-1 d-1 to  molC/gChl m2/uE
alpha_new = alpha_ /WtouE/864d2

!Carbon-specific chlorophyll quota, uit conversion:
thetaA = theta_ / 12d0     ![mg Chl mmol C] to [g Chl g C-1]

!Effective cross-section of the PSU [m2 uE-1] Nikolaou et al. (2016)
Sigma = Beta * thetaA**Kappa

!Repair constant of a photosynthetic unit [s-1], following Han et al. (2001):

!Nutrient limitation index
Lno3 = (QN_ - QNmin_) / (QNmax_ - QNmin_)

!Kr0 depends on alpha_ using an empirical equation
Kr0 = a_ * (alpha_new / b_)**v_

Kr = Kr0 * Lno3

if (Kr < 1d-10) then
  Ainf = 0.d0
else
  !Ratio of damage to repair constants [s]:
  K  = Kd / Kr
  
  !Calculate photoinhibition [nd]:
  Ainf = 1d0 / (1d0 + Tau * Sigma * PARWm2 + K * Tau * Sigma**2 * PARWm2**2)
endif

return 
END FUNCTION Ainf
!------------------------------------------------------------------------------------------------

END MODULE
!###############################################################################

!###############################################################################
!
Module gmk
!-------------------------------------------------------------------------------
! Trait_functions --- utility functions for the PIBM GMK98_ToptSizeLight model
!
!-------------------------------------------------------------------------------
use params

!This module provides several functions calculating phytoplankton physiological rates as a function of environmental conditions (e.g., temperature) and traits
implicit none

private

public :: GMK98_Ind_TempSizeLight

CONTAINS

SUBROUTINE GMK98_Ind_TempSizeLight(Temp, PAR, NO3, Topt_, C, N, Chl, Cdiv, alphaChl_, dC, dN, dChl)
!-------------------------------------------------------------------------------
USE Trait_functions, only : temp_Topt, PHY_C2Vol, Ainf, Pmax_size, respiration
!USE params,          only : thetaNmax, mu0, rhoChl_L, QNmin_a, QNmin_b
!USE params,          only : QNmax_a, QNmax_b, KN_a, KN_b, nx, pi
!ML USE state_variables, only : NO3_min

implicit none

! ML took from state_variables.f90
real,    parameter :: NO3_min = 0.01  !Minimal NO3 concentration

!Declaration on variables:
real, intent(in)  :: Temp             !Associated temperarure [degree C]
real, intent(in)  :: PAR              !Associated PAR [W m-2]
real, intent(in)  :: NO3              !Associated NO3 concentration [mmol N m-3]
real, intent(in)  :: C                !Current cellular carbon [pmol C cell-1]
real, intent(in)  :: N                !Current cellular nitrogen [pmol N cell-1]
real, intent(in)  :: Chl              !Current cellular Chl [pg C cell-1]
real, intent(in)  :: Cdiv             !Cellular carbon content threshold for division [pmol cell-1]

real, intent(in)  :: Topt_            !Optimal temperature [degree C]
real, intent(in)  :: alphaChl_        !Slope of the P-I curve [Unit the same as aI0]
real, intent(out) :: dN               !Changes in the cellular nitrogen content [pmol N cell-1 d-1]
real, intent(out) :: dC               !Changes in the cellular carbon content [pmol C cell-1 d-1]
real, intent(out) :: dChl             !Changes in the cellular Chl content [pg Chl cell-1 d-1]
real              :: Vol    = 0d0     !Cell volume of phytoplankton [um3]
real              :: ESD_   = 0d0     !ESD of phytoplankton [um]
real              :: QNmin  = 0.05    !Minimal N:C ratio [mmol N mmol C]
real              :: QNmax  = 0.18    !Maximal N:C ratio [mmol N mmol C]
real              :: dQN    = 0.13    !(Qmax - Qmin) [mmol N mmol C]
real              :: QN     = 0.      !Current N:C ratio [mmol N mmol C]
real              :: theta  = 0.      !Current Chl:C ratio [mg Chl mmol C]
real              :: VCN    = 0.      !DIN uptake rate by phytoplankton [mol N mol C-1 d-1]
real              :: Lno3   = 0.      !Nutrient limitation [pmol N m-3]
real              :: SI     = 0.      !Light limitation index [fpar]
real              :: PCmax  = 0.      !Maximal photosynthesis rate (regulated by QN) [d-1]
real              :: PC     = 0.      !Carbon specific rate of photosynthesis [d-1]
real              :: rhoChl = 0.      !Phyto C production devoted to Chl synthesis [mg Chl mmol C-1]
real              :: Ik     = 0.      !Saturation parameter for the PI curve [W m-2 s-1]
real              :: A      = 0.      !Photoinhibition, following Nikolau et al. (2016)

real, parameter   :: RC     = 0.1d0   !Basic C respiration rate [d-1]
real, parameter   :: RN     = 0.1d0   !Basic N respiration rate [d-1]
real, parameter   :: RChl   = 0.1d0   !Basic Chl respiration rate [d-1]
real              :: RcT    = 0.d0    !Temperature dependent C-based respiration rate [d-1]
real              :: RNT    = 0.d0    !Temperature dependent N-based respiration rate [d-1]
real              :: RChlT  = 0.d0    !Temperature dependent Chl-based respiration rate [d-1]

real, parameter   :: zeta   = 3.0d0   !Cost of biosynthesis [mol C mol N-1]

! Maximal specific N uptake as a function of temp. under resource (nutrient and light) replete conditions [mol N mol C-1 d-1]:
real              :: Vcref  = 0.

!Maximal growth rate as a function of temperature under resource (nutrient and light) replete conditions [uM]:
real              :: muT    = 0.

!Kn is an allometric function of Vol (Cdiv) (Edwards et al. 2012) [uM]:
real              :: KN     = 0.      !Half-saturation constant [uM]
real, parameter   :: a1 = 0.d0        !Allometric exponent between mumax and alphaChl
real              :: CDiv1 = 0.       !CDiv value with ug C per cell
!End of declaration

if (C .le. 0d0) then
  print *, 'C was 0 and I returned!'
   dN   = 0.d0
   dC   = 0.d0
   dChl = 0.d0
   return
endif

!Current N:C ratio [mmol N mmol C]:
QN = N/C

!Current Chl:C ratio [mg Chl mmol C]
theta = Chl/C

!Convert phytoplankton CDiv to Volume:
Vol = PHY_C2Vol(CDiv)

!Convert Volume to ESD:
ESD_ = (6.d0*Vol/pi)**0.3333333 

!Nitrate half-saturation constant of phyto growth based on cell volume [uM]:
KN = KN_a * Vol**KN_b

!Minimal N:C ratio [mmol N mmol C] following Ward et al. (2012):
QNmin = QNmin_a * Vol**QNmin_b

!Maximal N:C ratio [mmol N mmol C] following Maranon et al. (2013)
QNmax = QNmax_a * Vol**QNmax_b

!Constrain QN between QNmin and QNmax due to numerical issues
QN = max(min(QN, QNmax), QNmin)

!(Qmax - Qmin) [mmol N mmol C]:
dQN = QNmax - QNmin

!Maximal growth rate as a function of temperature under resource (nutrient and light) replete conditions:
!mu0 should be a function of alphaChl
muT = mu0 * exp(a1 * (alphaChl_ - .1)) !0.1 is the average alphaChl value

!Temperature dependent maximal growth rate at 1 ug C cell-1
muT = temp_Topt(Temp, muT, Topt_)

!Apply the size-scaling relationship following Wirtz (2011)
muT = Pmax_size(ESD_, muT)

!Assuming the same temperature dependence of nutrient uptake rate as on photosynthesis rate.
!QNmax/QNmin may be also a function of temperature which needs being further investigated.
Vcref  = muT * QNmax

!Assume the same temperature dependence of respiration as on photosynthetic rate (long-term adaptation; Barton et al. 2020):
RcT   = temp_Topt(Temp, Rc,   Topt_)
RcT   = respiration(ESD_, RcT)

RNT   = temp_Topt(Temp, RN,   Topt_)
RNT   = respiration(ESD_, RNT)

RChlT = temp_Topt(Temp, RChl, Topt_)
RChlT = respiration(ESD_, RChlT)

!Nutrient limitation [nd]:
Lno3 = (QN - QNmin) / dQN

if (Lno3 .le. 0d0) then
   PC = 0d0
else
   !Maximal photosynthesis rate (regulated by QN) [d-1]:
   PCmax = muT * Lno3

   !Light saturation parameter [W m-2 d-1]:
   Ik = PCmax / alphachl_ / theta

   !Calculate the fraction of open PSU [nd]:
   if (PAR > 0.) then !Photoinhibition
      A = Ainf(PAR, alphachl_, QN, QNmin, QNmax, theta)
   else
      A = 1d0
   endif

   !Light limitation index [nd]:
   SI = -A*PAR/Ik

   if (abs(SI) < 1d-10) then
      SI = 0.d0
   else
      SI = 1.d0 - exp(SI)
   endif

   !Photosynthesis rate [d-1]:
   PC = PCmax * SI
Endif

!Define rhochl [g Chl mol C-1]: fraction of phytoplankton carbon production that is devoted to Chl synthesis.
!If dark, assume that rhochl equaled the value calculated for the end of the preceding light period.
if (PAR <= 0d0) then
   rhochl   = rhoChl_L
else
   rhochl   = thetaNmax * PC / alphachl_ / theta / PAR
   rhoChl_L = rhochl
endif

!DIN uptake rate by phytoplankton [mol N mol C-1 d-1]:
VCN = Vcref * (NO3 - NO3_min)/ (NO3 + KN) * ((QNmax - QN) / dQN)**nx  !Vcref already temperature dependent
VCN = max(VCN, 0d0)

!Changes of cellular carbon [d-1]:
dC = C * (PC - zeta*VCN - RcT)

!Changes of cellular nitrogen [pmol N cell-1 d-1]:
!RNT has to be zero to avoid continuous decline of N per cell
dN = N * (VCN/QN - RNT)

!Changes of cellular Chl [d-1]:
dChl = Chl * (rhochl*VCN / theta - RChlT)

return
END subroutine GMK98_Ind_TempSizeLight

END MODULE
!###############################################################################

module mGf90
implicit none

private
public srand_mtGaus
contains

real function srand_mtGaus(N, mean, cvm) result(y)
implicit none

!The dimension of the multivariate vector
integer, intent(in)   :: N

real,    intent(in)   :: mean(N)  !Mean values of the multivariate distribution

real,    intent(in)   :: cvm(N,N) !Covariance matrix of the multivariate distribution
dimension             :: y(N)
integer    :: nullty  ! i/o for cholesky subroutine
integer    :: error  
integer    :: i, j
real       :: cvm_(N*(N+1)/2)
real       :: chol(N*(N+1)/2)

!Convert cvm to a vector with length N*(N+1)/2
do i = 1, N
  do j = 1,i
    cvm_(i*(i-1)/2+j) = cvm(i,j)
  enddo
enddo

call cholesky(cvm_,N,N*(N+1)/2,chol,nullty,error)

y = multiGauss(chol,mean,N)
end function srand_mtGaus

real function gasdev() result(gval)
!This function generates an array of independent Gaussian variates
implicit none
real :: FAC, R, V1, V2, X
real, save :: GSET
integer, save :: ISET = 0

IF (ISET.EQ.0) THEN
   R = 99 
   do while( R .ge. 1.0d0 )
      call random_number(V1)
      V1= 2d0*V1 - 1d0

      call random_number(V2)
      V2= 2d0*V2 - 1d0

      R = V1**2 + V2**2
   end do
   FAC  = SQRT( -2.0d0*LOG(R)/R )
   GSET = V1*FAC
   gval = V2*FAC
   ISET = 1
ELSE
   gval=GSET
   ISET=0
ENDIF

RETURN
END function gasdev

subroutine cholesky ( a, n, nn, u, nullty, ifault )
!*****************************************************************************80
!
!! CHOLESKY computes the Cholesky factorization of a PDS matrix.
!
!  Discussion:
!
!    For a positive definite symmetric matrix A, the Cholesky factor U
!    is an upper triangular matrix such that A = U' * U.
!
!    This routine was originally named "CHOL", but that conflicted with
!    a built in MATLAB routine name.
!
!  Modified:
!
!    01 February 2008
!
!  Author:
!
!    Michael Healy
!    Modifications by AJ Miller.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Michael Healy,
!    Algorithm AS 6:
!    Triangular decomposition of a symmetric matrix,
!    Applied Statistics,
!    Volume 17, Number 2, 1968, pages 195-197.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A((N*(N+1))/2), a positive definite matrix 
!    stored by rows in lower triangular form as a one dimensional array, 
!    in the sequence
!    A(1,1),
!    A(2,1), A(2,2),
!    A(3,1), A(3,2), A(3,3), and so on.
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Input, integer NN, the dimension of the array used to store A, 
!    which should be at least (N*(N+1))/2.
!
!    Output, real ( kind = 8 ) U((N*(N+1))/2), an upper triangular matrix,
!    stored by columns, which is the Cholesky factor of A.  The program is
!    written in such a way that A and U can share storage.
!    It is also the same as lower triangular matrix stored by rows (by B. Chen on 12 Aug 2019)
!
!    Output, integer ( kind = 4 ) NULLTY, the rank deficiency of A.  If NULLTY is zero,
!    the matrix is judged to have full rank.
!
!    Output, integer ( kind = 4 ) IFAULT, an error indicator.
!    0, no error was detected;
!    1, if N < 1;
!    2, if A is not positive semi-definite.
!    3, NN < (N*(N+1))/2.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) ETA, should be set equal to the smallest positive
!    value such that 1.0 + ETA is calculated as being greater than 1.0 in the
!    accuracy being used.
!
  implicit none

  integer ( kind = 4 ), intent(in) :: n
  integer ( kind = 4 ), intent(in) :: nn

  real    ( kind = 8 ), intent(in) :: a(nn)
  real    ( kind = 8 ), parameter :: eta = 1.0D-19
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ), intent(out) :: ifault
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ), intent(out) :: nullty
  real    ( kind = 8 ), intent(out) :: u(nn)
  real    ( kind = 8 ) w
  real    ( kind = 8 ) x

  ifault = 0
  nullty = 0

  if ( n <= 0 ) then
    ifault = 1
    return
  end if

  if ( nn < ( n * ( n + 1 ) ) / 2 ) then
    ifault = 3
    return
  end if

  j = 1
  k = 0
  ii = 0
!
!  Factorize column by column, ICOL = column number.
!
  do icol = 1, n

    ii = ii + icol
    x = eta * eta * a(ii)
    l = 0
    kk = 0
!
!  IROW = row number within column ICOL.
!
    do irow = 1, icol

      kk = kk + irow
      k = k + 1
      w = a(k)
      m = j

      do i = 1, irow - 1
        l = l + 1
        w = w - u(l) * u(m)
        m = m + 1
      end do

      l = l + 1

      if ( irow == icol ) then
        exit
      end if

      if ( u(l) /= 0.0D+00 ) then

        u(k) = w / u(l)

      else

        u(k) = 0.0D+00

        if ( abs ( x * a(k) ) < w * w ) then
          ifault = 2
          return
        end if

      end if

    end do
!
!  End of row, estimate relative accuracy of diagonal element.
!
    if ( abs ( w ) <= abs ( eta * a(k) ) ) then

      u(k) = 0.0D+00
      nullty = nullty + 1

    else

      if ( w < 0.0D+00 ) then
        ifault = 2
        return
      end if

      u(k) = sqrt ( w )

    end if

    j = j + icol

  end do

  return
end subroutine cholesky

!-----------------------------------------------------------------------------
subroutine subchl ( a, b, n, u, nullty, ifault, ndim, det )

!*****************************************************************************80
! 
!! SUBCHL computes the Cholesky factorization of a (subset of a) PDS matrix.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    FORTRAN77 version by Michael Healy, PR Freeman
!    FORTRAN90 version by  John Burkardt
!
!  Reference:
!
!    PR Freeman,
!    Remark AS R44:
!    A Remark on AS 6 and AS7: Triangular decomposition of a symmetric matrix
!    and Inversion of a positive semi-definite symmetric matrix,
!    Applied Statistics,
!    Volume 31, Number 3, 1982, pages 336-339.
!
!    Michael Healy,
!    Algorithm AS 6:
!    Triangular decomposition of a symmetric matrix,
!    Applied Statistics,
!    Volume 17, Number 2, 1968, pages 195-197.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A((M*(M+1))/2), a positive definite matrix 
!    stored by rows in lower triangular form as a one dimensional array, 
!    in the sequence
!    A(1,1),
!    A(2,1), A(2,2),
!    A(3,1), A(3,2), A(3,3), and so on.  
!    In the simplest case, M, the order of A, is equal to N.
!
!    Input, integer ( kind = 4 ) B(N), indicates the order in which the
!    rows and columns of A are to be used.  In the simplest case, 
!    B = (1,2,3...,N).
!
!    Input, integer ( kind = 4 ) N, the order of the matrix, that is, 
!    the matrix formed by using B to select N rows and columns of A.
!
!    Output, real ( kind = 8 ) U((N*(N+1))/2), an upper triangular matrix,
!    stored by columns, which is the Cholesky factor of A.  The program is
!    written in such a way that A and U can share storage.
!
!    Output, integer ( kind = 4 ) NULLTY, the rank deficiency of A.  
!    If NULLTY is zero, the matrix is judged to have full rank.
!
!    Output, integer ( kind = 4 ) IFAULT, an error indicator.
!    0, no error was detected;
!    1, if N < 1;
!    2, if A is not positive semi-definite.
!
!    Input, integer ( kind = 4 ) NDIM, the dimension of A and U, which might 
!    be presumed to be (N*(N+1))/2.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndim

  real    ( kind = 8 ) a(ndim)
  integer ( kind = 4 ) b(n)
  real    ( kind = 8 ) det
  real    ( kind = 8 ), parameter :: eta = 1.0D-09
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nullty
  real    ( kind = 8 ) u(ndim)
  real    ( kind = 8 ) w
  real    ( kind = 8 ) x

  ifault = 0
  nullty = 0
  det = 1.0D+00

  if ( n <= 0 ) then
    ifault = 1
    return
  end if

  ifault = 2
  j = 1
  k = 0

  do icol = 1, n

    ij = ( b(icol) * ( b(icol) - 1 ) ) / 2
    ii = ij + b(icol)
    x = eta * eta * a(ii)
    l = 0

    do irow = 1, icol

      kk = ( b(irow) * ( b(irow) + 1 ) ) / 2
      k = k + 1
      jj = ij + b(irow)
      w = a(jj)
      m = j

      do i = 1, irow - 1
        l = l + 1
        w = w - u(l) * u(m)
        m = m + 1
      end do

      l = l + 1

      if ( irow == icol ) then
        exit
      end if

      if ( u(l) /= 0.0D+00 ) then

        u(k) = w / u(l)

      else

        if ( abs ( x * a(kk) ) < w * w ) then
          ifault = 2
          return
        end if

        u(k) = 0.0D+00

      end if

    end do

    if ( abs ( eta * a(kk) ) <= abs ( w ) ) then

      if ( w < 0.0D+00 ) then
        ifault = 2
        return
      end if

      u(k) = sqrt ( w )

    else

      u(k) = 0.0D+00
      nullty = nullty + 1

    end if

    j = j + icol
    det = det * u(k) * u(k)

  end do

  return
end subroutine subchl

real function multigauss(R,mu,n) result(G)
!!$ Generates a one dimensional array of length n containing multivariate Gaussian pseudo-random numbers
!!$ where R is the lower-triangular Cholesky factor of the covariance matrix of the desired distribution 
!!$ and mu is a one dimensional array of length n, containing the mean value for each random variable. 
!!$ R must be a one dimensional array, containing the lower-triangular matrix stored by row, starting from the 
!!$ uppermost, leftmost entry (first row, first column). 
   implicit none
   integer, intent(in) :: n
   real, intent(in) :: R(n*(n+1)/2)
   real, intent(in) :: mu(n)
   real :: Nu(n)
   dimension :: G(n)
   integer :: i, j
   
   if( n*(n+1)/2 .ne. size(R) ) then
      write(6,*) ' n*(n+1)/2 != size(R) in multiGauss '
      write(6,*) ' Cholesky factor matrix has size ', size(R) 
      write(6,*) ' but you are requesting a vector of size ',n
      stop
   else
      
!!$ generate an array of independent Gaussian variates
      do j = 1, n
         Nu(j) = gasdev()
      end do

!!$ start with the desired mean, and add the product of 
!!$ [the lower-triangular Cholesky factor of the covariance matrix] x [ the vector of iid Gaussian variates]
      do i = 1, n
         G(i) = mu(i)
         do j =  1,i
            G(i) = G(i) + R( i*(i-1)/2 + j ) * Nu(j)
         end do
      end do

   end if
   return
 end function multigauss
end module mGf90
