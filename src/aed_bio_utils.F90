!###############################################################################
!#                                                                             #
!# aed_bio_utils.F90                                                           #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2013-2026 : The University of Western Australia                  #
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
!# Created August 2011                                                         #
!#                                                                             #
!###############################################################################

#include "aed.h"

MODULE aed_bio_utils
!-------------------------------------------------------------------------------
!  aed_bio_utils --- utility functions for phytoplankton & macroalgae models
!-------------------------------------------------------------------------------
   USE aed_core

   USE aed_util,ONLY : exp_integral, &
                       aed_bio_temp_function, &
                       fTemp_function

   IMPLICIT NONE

   PRIVATE   ! By default make everything private
!
   PUBLIC phyto_data_t, phyto_param_t
   PUBLIC phyto_salinity, phyto_fN, phyto_fP, phyto_fSi
   PUBLIC phyto_internal_nitrogen, phyto_internal_phosphorus
   PUBLIC photosynthesis_irradiance, bio_respiration
   PUBLIC ino3, inh4, idon, in2, ifrp, idop
   PUBLIC findMin
!
   TYPE phyto_data_t
      ! General Attributes
      CHARACTER(64) :: p_name
      AED_REAL :: p0,Xcc
      ! Growth rate parameters
      INTEGER  :: fT_Method
      AED_REAL :: R_growth, theta_growth, T_std, T_opt, T_max, kTn, aTn, bTn
      ! Light configuration and parameters
      INTEGER  :: lightModel
      AED_REAL :: I_K, I_S, KePHY
      ! Respiration parameters
      AED_REAL :: f_pr, R_resp, k_fdom, k_fres, theta_resp
      ! Salinity parameters
      INTEGER  :: salTol
      AED_REAL :: S_bep, S_maxsp, S_opt
      ! Nitrogen parameters
      INTEGER  :: simDINUptake, simDONUptake, simNFixation, simINDynamics
      AED_REAL :: N_o, K_N, X_nmin, X_nmax, X_ncon, R_nuptake, k_nfix, R_nfix
      ! Phosphorus parameters
      INTEGER  :: simDIPUptake, simIPDynamics
      AED_REAL :: P_0, K_P, X_pmin, X_pmax, X_pcon, R_puptake
      ! Silica parameters
      INTEGER  :: simSiUptake
      AED_REAL :: Si_0, K_Si, X_sicon
      ! Carbon parameters
      INTEGER  :: simCUptake, dic_mode
      ! Sedimentation parameters
      INTEGER  :: settling
      AED_REAL :: w_p, d_phy, rho_phy, f1, f2, c1, c3
      ! Resuspension parameters
      AED_REAL  :: resuspension, tau_0
      ! Growth form (benthic, water, surface)
      INTEGER  :: growth_form, slough_model
      ! Physiology model selector and per-group sloughing parameters.
      !   NOTE(2026-09-23): used by aed_macroalgae only. physiology_model
      !   chooses the growth formulation for the group (0 = generic AED,
      !   1 = CGM, 2 = GLCMv3); slough_rate/_stress/_burial hold the
      !   per-group values resolved from the module-wide &aed_macroalgae
      !   scalars and the optional *_p override arrays.
      INTEGER  :: physiology_model
      AED_REAL :: slough_rate, slough_stress, slough_burial
      ! Particle parameters
      INTEGER  :: simSplit
      AED_REAL :: X_cinit, X_ninit, X_pinit, X_chlinit, Cdiv, n0, Lnalphachl, mort_prob
      AED_REAL :: RC, RN, RP, RChl, zeta_N, zeta_P, a1
      AED_REAL :: nx, thetaNmax, thetaPmax, QNmin_a, QNmin_b, QNmax_a, QNmax_b, QPmin_a
      AED_REAL :: QPmin_b, QPmax_a, QPmax_b, KN_a, KN_b, KPho_a, KPho_b
      AED_REAL :: a_c2vol, b_c2vol, a_pmax, rho, rho_star, b_rho, V_s, Ea0, Ed0, Ei, beta
      AED_REAL :: phi, Tau, Beta_Ainf, Kappa, Kd_Ainf, a_, b_, v_
      AED_REAL :: n_split
      !# 2026-07: see f_colony in phyto_param_t below. Effective hydrodynamic size
      !# multiplier, applied to ESD in the Stokes settling term only.
      AED_REAL :: f_colony
      !# 2026-07 (S5): dark survival / dormancy. f_dorm = carbon fraction of Cdiv
      !# below which a cell becomes dormant instead of continuing to starve;
      !# r_dorm = maintenance respiration multiplier while dormant; I_germ =
      !# PAR threshold for germination. See handoff S5.
      AED_REAL :: f_dorm, r_dorm, I_germ
      AED_REAL :: dorm_width   !# 2026-07 (S5b): width of the smooth dormancy ramp
      !# Non-photochemical quenching. NPQ_max is percent reduction of effective
      !# PAR at saturation; rates are per day; E_thresh is W/m2 after spectral
      !# weighting.
      AED_REAL :: NPQ_max, NPQ_k_on, NPQ_k_off, NPQ_E_thresh
      !# 2026-08: NPQ_mode selects the quenching formulation.
      !#   0 = legacy rate-based relaxation toward an equilibrium set by INSTANTANEOUS
      !#       irradiance (k_on/k_off timescales ~5/20 min, so almost no memory).
      !#   1 = Ranjbar et al. (2024) cumulative light DOSE, which is what makes this a
      !#       light-HISTORY model: CL_t = max(I_t - PRR + CL_{t-1}, 0), and NPQ is
      !#       integrated from CL rather than solved for instantaneous light.
      !# NPQ_PRR is the photo-recovery rate in W/m2 (Ranjbar quote 800 umol/m2/s; divide
      !# by WtouE = 4.57 because every PAR variable in this module is W/m2 - see the unit
      !# note in aed_phyto_abm.F90). NPQ_alpha is their fitted scaling on d(NPQ)/dt.
      INTEGER  :: NPQ_mode
      AED_REAL :: NPQ_PRR, NPQ_alpha, NPQ_dtref_d
      !# 2026-08: NPQ_apply separates WHERE quenching acts from WHICH equation computes it.
      !#   0 = fluorescence diagnostic only. This is what Ranjbar et al. (2024) actually
      !#       did: NPQ is an observation operator mapping biomass onto the phycocyanin
      !#       signal a probe would read, and biomass itself is untouched.
      !#   1 = additionally reduce the PAR handed to the growth engine, on the grounds
      !#       that quenching physically diverts absorbed excitation away from
      !#       photochemistry. More mechanistic, but beyond what the paper validated.
      INTEGER  :: NPQ_apply
      !# Platt et al. (1980) photoinhibition on the GMK98 P-I curve. 0 = OFF and is
      !# bit-identical to the un-inhibited curve, so it doubles as the control.
      AED_REAL :: beta_inh
      !# Water Research 2022 colony aggregation/disaggregation. simColony 0 = static
      !# f_colony (unchanged). colony_dD_dtref_h is the time unit of dD in their Eq 10,
      !# which the paper does not state - see the handoff.
      INTEGER  :: simColony
      AED_REAL :: colony_dD_dtref_h, colony_d_max
      !# 2026-08-19: de-ballasting (density DECREASE, shedding carbohydrate ballast to
      !# regain buoyancy) was gated on r_dorm - the SAME rate that suppresses general
      !# C/N/P/Chl maintenance loss while dormant. A particle that ballasts down and then
      !# goes dormant at the bed needs a genuinely independent recovery rate: r_dorm tuned
      !# low for general survival makes its one escape route correspondingly slow. Default
      !# equals the ORIGINAL r_dorm value (0.15), a neutral starting point distinct from
      !# whatever r_dorm is tuned to for maintenance loss. Appended at the END of the type
      !# (not inserted mid-structure) after a mid-structure insertion produced unrelated
      !# runtime corruption elsewhere (see handoff) - append-only is the safe pattern here.
      AED_REAL :: r_deballast
      !# 2026-08-19: I_Kb is the light scale of the BUOYANCY (carbohydrate ballasting)
      !# response, and rho_max is a PER-GROUP density ceiling. Both exist because the
      !# quantities they replace were global/shared and carried the wrong meaning for
      !# gas-vacuolate cyanobacteria - see aed_phyto_abm.F90 for the measurements.
      AED_REAL :: I_Kb, rho_max
      !# 2026-08-21: rho_min is the PER-GROUP density FLOOR, the mirror of rho_max above.
      !# Added because the global aed.nml min_rho = 985 was binding on both gas-vacuolate
      !# groups - they sat at exactly 985.00 - which made real buoyancy regulation
      !# impossible. Reynolds, Oliver & Walsby (1987) Table 2 gives Anabaena flos-aquae
      !# (now Dolichospermum) 920-1030 and Oscillatoria rubescens (now Planktothrix
      !# rubescens) 990-1065, so a single global floor cannot represent both.
      !# <= 0 falls back to the global min_rho, exactly as rho_max does, so every group is
      !# bit-identical until a value is set. Appended last to preserve type layout.
      AED_REAL :: rho_min
      !# 2026-08-22: Wallace & Hamilton lagged buoyancy regulation. Kromkamp & Walsby (1990)
      !# and Visser et al. (1997) assume cell density responds to irradiance INSTANTLY; Wallace
      !# & Hamilton, Limnol. Oceanogr. 44(2):273-281 (1999) measured a ~20 min physiological lag
      !# while carbohydrate storage ramps up, and Wallace & Hamilton, J. Plankton Res.
      !# 22(6):1127-1138 (2000) showed that a colony mixed through the light gradient faster
      !# than that lag never fully ballasts and ends the mixing event MORE buoyant - the
      !# mechanism behind persistent surface blooms.
      !#   buoy_model  0 = legacy Webb c1*(1-exp(-I/I_Kb)) - c3   (default; bit-identical)
      !#               1 = Wallace & Hamilton lagged model, integrated in glm_ptm.c
      !# Rates are kg/m3/day and irradiance is W/m2, matching ip_par. Converted from the paper's
      !# kg/m3/min and umol/m2/s with 1 W/m2 = 4.57 umol photons/m2/s; see
      !# analysis/test_buoyancy_wh.py, which reproduces their Eq 7 and Figs 2-3 and asserts the
      !# conversion. DO NOT take c1 from the 2000 paper: its 0.146 ug/g/s is exactly 100x low
      !# (its c3 converts correctly, which is how the typo was identified).
      !# Appended last to preserve type layout, as rho_min was.
      INTEGER  :: buoy_model
      AED_REAL :: c1_buoy, c2_buoy, c3_buoy, KI_buoy, tau_buoy, I_dark, tau_dose
      !# k_resp: Visser et al. (1997) Eq 3 first-order ballast decay, /day. Appended last.
      AED_REAL :: k_resp
   END TYPE phyto_data_t


   ! %% NAMELIST   %% phyto_param_t
   TYPE phyto_param_t
      CHARACTER(64) :: p_name
      AED_REAL :: p_initial
      AED_REAL :: p0, w_p, Xcc, R_growth
      INTEGER  :: fT_Method
      AED_REAL :: theta_growth, T_std, T_opt, T_max
      INTEGER  :: lightModel
      AED_REAL :: I_K, I_S, KePHY
      ! Respiration parameters
      AED_REAL :: f_pr, R_resp, theta_resp, k_fres, k_fdom
      ! Salinity parameters
      INTEGER  :: salTol
      AED_REAL :: S_bep, S_maxsp, S_opt
      ! Nitrogen parameters
      INTEGER  :: simDINUptake, simDONUptake, simNFixation, simINDynamics
      AED_REAL :: N_o, K_N, X_ncon, X_nmin, X_nmax, R_nuptake, k_nfix, R_nfix
      ! Phosphorus parameters
      INTEGER  :: simDIPUptake, simIPDynamics
      AED_REAL :: P_0, K_P, X_pcon, X_pmin, X_pmax, R_puptake
      ! Silica parameters
      INTEGER  :: simSiUptake
      AED_REAL :: Si_0, K_Si, X_sicon
      !
      AED_REAL :: c1, c3, f1, f2, d_phy
      ! Particle parameters
      INTEGER  :: simSplit
      AED_REAL :: X_cinit, X_ninit, X_pinit, X_chlinit, Cdiv, n0, Lnalphachl, mort_prob
      AED_REAL :: RC, RN, RP, RChl, zeta_N, zeta_P, a1
      AED_REAL :: nx, thetaNmax, thetaPmax, QNmin_a, QNmin_b, QNmax_a, QNmax_b, QPmin_a
      AED_REAL :: QPmin_b, QPmax_a, QPmax_b, KN_a, KN_b, KPho_a, KPho_b
      AED_REAL :: a_c2vol, b_c2vol, a_pmax, rho, rho_star, b_rho, V_s, Ea0, Ed0, Ei, beta
      AED_REAL :: phi, Tau, Beta_Ainf, Kappa, Kd_Ainf, a_, b_, v_
      AED_REAL :: n_split
      !# 2026-07: effective hydrodynamic size multiplier applied to ESD in the STOKES
      !# settling term ONLY. Cell ESD (from a_c2vol/b_c2vol) correctly drives allometric
      !# nutrient kinetics, but filamentous/colonial taxa sink and float as aggregates far
      !# larger than one cell. Stokes velocity scales as ESD^2, so at the ~2 um cell scale
      !# buoyancy is ~1200x weaker than turbulent mixing and no density regulation can
      !# hold station (handoff 33-Z). Default 1.0 = unchanged behaviour.
      AED_REAL :: f_colony
      !# 2026-07 (S5): dark survival / dormancy. f_dorm = carbon fraction of Cdiv
      !# below which a cell becomes dormant instead of continuing to starve;
      !# r_dorm = maintenance respiration multiplier while dormant; I_germ =
      !# PAR threshold for germination. See handoff S5.
      AED_REAL :: f_dorm, r_dorm, I_germ
      AED_REAL :: dorm_width   !# 2026-07 (S5b): width of the smooth dormancy ramp
      !# Non-photochemical quenching. NPQ_max is percent reduction of effective
      !# PAR at saturation; rates are per day; E_thresh is W/m2 after spectral
      !# weighting.
      AED_REAL :: NPQ_max, NPQ_k_on, NPQ_k_off, NPQ_E_thresh
      !# 2026-08: see phyto_data_t above. 0 = legacy rate-based, 1 = Ranjbar cumulative
      !# dose. NPQ_PRR in W/m2 (not umol/m2/s), NPQ_alpha dimensionless.
      INTEGER  :: NPQ_mode
      AED_REAL :: NPQ_PRR, NPQ_alpha, NPQ_dtref_d
      !# 2026-08: NPQ_apply separates WHERE quenching acts from WHICH equation computes it.
      !#   0 = fluorescence diagnostic only. This is what Ranjbar et al. (2024) actually
      !#       did: NPQ is an observation operator mapping biomass onto the phycocyanin
      !#       signal a probe would read, and biomass itself is untouched.
      !#   1 = additionally reduce the PAR handed to the growth engine, on the grounds
      !#       that quenching physically diverts absorbed excitation away from
      !#       photochemistry. More mechanistic, but beyond what the paper validated.
      INTEGER  :: NPQ_apply
      !# Platt et al. (1980) photoinhibition on the GMK98 P-I curve. 0 = OFF and is
      !# bit-identical to the un-inhibited curve, so it doubles as the control.
      AED_REAL :: beta_inh
      !# Water Research 2022 colony aggregation/disaggregation. simColony 0 = static
      !# f_colony (unchanged). colony_dD_dtref_h is the time unit of dD in their Eq 10,
      !# which the paper does not state - see the handoff.
      INTEGER  :: simColony
      AED_REAL :: colony_dD_dtref_h, colony_d_max
      !# 2026-08-19: de-ballasting (density DECREASE, shedding carbohydrate ballast to
      !# regain buoyancy) was gated on r_dorm - the SAME rate that suppresses general
      !# C/N/P/Chl maintenance loss while dormant. A particle that ballasts down and then
      !# goes dormant at the bed needs a genuinely independent recovery rate: r_dorm tuned
      !# low for general survival makes its one escape route correspondingly slow. Default
      !# equals the ORIGINAL r_dorm value (0.15), a neutral starting point distinct from
      !# whatever r_dorm is tuned to for maintenance loss. Appended at the END of the type
      !# (not inserted mid-structure) after a mid-structure insertion produced unrelated
      !# runtime corruption elsewhere (see handoff) - append-only is the safe pattern here.
      AED_REAL :: r_deballast
      !# 2026-08-19: I_Kb is the light scale of the BUOYANCY (carbohydrate ballasting)
      !# response, and rho_max is a PER-GROUP density ceiling. Both exist because the
      !# quantities they replace were global/shared and carried the wrong meaning for
      !# gas-vacuolate cyanobacteria - see aed_phyto_abm.F90 for the measurements.
      AED_REAL :: I_Kb, rho_max
      !# 2026-08-21: rho_min is the PER-GROUP density FLOOR, the mirror of rho_max above.
      !# Added because the global aed.nml min_rho = 985 was binding on both gas-vacuolate
      !# groups - they sat at exactly 985.00 - which made real buoyancy regulation
      !# impossible. Reynolds, Oliver & Walsby (1987) Table 2 gives Anabaena flos-aquae
      !# (now Dolichospermum) 920-1030 and Oscillatoria rubescens (now Planktothrix
      !# rubescens) 990-1065, so a single global floor cannot represent both.
      !# <= 0 falls back to the global min_rho, exactly as rho_max does, so every group is
      !# bit-identical until a value is set. Appended last to preserve type layout.
      AED_REAL :: rho_min
      !# 2026-08-22: Wallace & Hamilton lagged buoyancy regulation. Kromkamp & Walsby (1990)
      !# and Visser et al. (1997) assume cell density responds to irradiance INSTANTLY; Wallace
      !# & Hamilton, Limnol. Oceanogr. 44(2):273-281 (1999) measured a ~20 min physiological lag
      !# while carbohydrate storage ramps up, and Wallace & Hamilton, J. Plankton Res.
      !# 22(6):1127-1138 (2000) showed that a colony mixed through the light gradient faster
      !# than that lag never fully ballasts and ends the mixing event MORE buoyant - the
      !# mechanism behind persistent surface blooms.
      !#   buoy_model  0 = legacy Webb c1*(1-exp(-I/I_Kb)) - c3   (default; bit-identical)
      !#               1 = Wallace & Hamilton lagged model, integrated in glm_ptm.c
      !# Rates are kg/m3/day and irradiance is W/m2, matching ip_par. Converted from the paper's
      !# kg/m3/min and umol/m2/s with 1 W/m2 = 4.57 umol photons/m2/s; see
      !# analysis/test_buoyancy_wh.py, which reproduces their Eq 7 and Figs 2-3 and asserts the
      !# conversion. DO NOT take c1 from the 2000 paper: its 0.146 ug/g/s is exactly 100x low
      !# (its c3 converts correctly, which is how the typo was identified).
      !# Appended last to preserve type layout, as rho_min was.
      INTEGER  :: buoy_model
      AED_REAL :: c1_buoy, c2_buoy, c3_buoy, KI_buoy, tau_buoy, I_dark, tau_dose
      !# k_resp: Visser et al. (1997) Eq 3 first-order ballast decay, /day. Appended last.
      AED_REAL :: k_resp
   END TYPE phyto_param_t
   ! %% END NAMELIST   %% phyto_param_t

!Module Locals
   INTEGER,PARAMETER :: ino3 = 1, inh4 = 2, idon = 3, in2 = 4, ifrp = 1, idop = 2

!===============================================================================
CONTAINS


!###############################################################################
SUBROUTINE phyto_internal_phosphorus(phytos,group,npup,phy,IP,primprod,        &
                                                 fT,pup,respiration,exudation, &
                                                     uptake,excretion,mortality)
!-------------------------------------------------------------------------------
! Calculates the biotic group internal phosphorus stores and fluxes
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(phyto_data_t),DIMENSION(:),INTENT(in)  :: phytos
   INTEGER,INTENT(in)                          :: group
   INTEGER,INTENT(in)                          :: npup
   AED_REAL,INTENT(in)                         :: phy
   AED_REAL,INTENT(in)                         :: IP
   AED_REAL,INTENT(in)                         :: primprod
   AED_REAL,INTENT(in)                         :: fT,pup,respiration,exudation
   AED_REAL,INTENT(out)                        :: uptake(:),excretion,mortality
!CONSTANTS
   AED_REAL,PARAMETER :: one_e_neg5 = 1e-5
!LOCALS
   AED_REAL :: tmpary1,tmpary2,theX_pcon
   INTEGER  :: c
!
!-------------------------------------------------------------------------------
!BEGIN
   uptake     = zero_
   excretion  = zero_
   mortality  = zero_

   ! Uptake of phosphorus
   IF (phytos(group)%simIPDynamics == 0 .OR. phytos(group)%simIPDynamics == 1) THEN
      ! Static phosphorus uptake function
      ! uptake = X_pcon * mu * phy

      theX_pcon = phytos(group)%X_pcon * phy
      DO c = 1,npup
         ! uptake is spread over relevant sources (assumes evenly)
         uptake(c) = - (theX_pcon/npup) * primprod
      END DO
   ELSEIF (phytos(group)%simIPDynamics == 2) THEN

      ! Dynamic phosphorus uptake function
      ! R_puptake * fT * phy * (X_pmax-[IP/phy])/(X_pmax-X_pmin) * (PO4/K_P+PO4])

      theX_pcon = IP
      tmpary1   = phytos(group)%R_puptake * fT * phy
      tmpary2   = MAX(zero_, phytos(group)%X_pmax - (IP / phy))
      tmpary1   = tmpary1 * tmpary2 / (phytos(group)%X_pmax-phytos(group)%X_pmin)
      uptake(1) =-tmpary1 * phyto_fP(phytos,group,frp=pup)      ! FRP
      uptake(2) = zero_                                         ! DOP
   ELSE
      ! Unknown phosphorus uptake function
      print *,'STOP: unknown simIPDynamics (',phytos(group)%simIPDynamics,') for: ',phytos(group)%p_name
      ERROR STOP 1
   ENDIF

   ! Release of phosphorus due to excretion from phytoplankton and
   ! contribution of mortality and excretion to OM

   excretion = (respiration*phytos(group)%k_fdom + exudation)*theX_pcon
   mortality = respiration*(1.0-phytos(group)%k_fdom)*theX_pcon
END SUBROUTINE phyto_internal_phosphorus
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE phyto_internal_nitrogen(phytos,group,do_N2uptake,phy,IN,primprod,   &
                                   fT,no3up,nh4up,a_nfix,respiration,exudation,&
                                   PNf,uptake,excretion,mortality)
!-------------------------------------------------------------------------------
! Calculates the biotic group internal nitrogen stores and fluxes
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(phyto_data_t),DIMENSION(:),INTENT(in)  :: phytos
   INTEGER,INTENT(in)                          :: group
   LOGICAL,INTENT(in)                          :: do_N2uptake
   AED_REAL,INTENT(in)                         :: phy
   AED_REAL,INTENT(in)                         :: IN
   AED_REAL,INTENT(in)                         :: primprod
   AED_REAL,INTENT(in)                         :: fT,no3up,nh4up
   AED_REAL,INTENT(inout)                      :: a_nfix
   AED_REAL,INTENT(in)                         :: respiration,exudation
   AED_REAL,INTENT(out)                        :: PNf
   AED_REAL,INTENT(out)                        :: uptake(:),excretion,mortality
!
!CONSTANTS
   AED_REAL,PARAMETER :: one_e_neg5 = 1e-5
!
!LOCALS
   AED_REAL  :: tmpary1,tmpary2,theX_ncon
!
!-------------------------------------------------------------------------------
!BEGIN
   uptake     = zero_
   excretion  = zero_
   mortality  = zero_

   ! Uptake of nitrogen
   IF (phytos(group)%simINDynamics == 0 .OR. phytos(group)%simINDynamics == 1) THEN
      ! Static nitrogen uptake function (assuming fixed stoichiometry)
      ! uptake = X_ncon * mu * phy

      theX_ncon = phytos(group)%X_ncon * phy
      uptake(1)  = -theX_ncon * primprod
   ELSEIF (phytos(group)%simINDynamics == 2) THEN
      ! Dynamic nitrogen uptake function
      ! R_nuptake * fT * phy * (X_nmax-IN/phy)/(X_nmax-X_nmin) * (DIN/[K_N+DIN])

      theX_ncon = IN
      tmpary1   = phytos(group)%R_nuptake * fT * phy
      tmpary2   = MAX(phytos(group)%X_nmax - (IN / phy),zero_)
      tmpary1   = tmpary1 * tmpary2 / (phytos(group)%X_nmax-phytos(group)%X_nmin)
      uptake(1) = tmpary1 * phyto_fN(phytos,group,din=no3up+nh4up)
      uptake(1) = -uptake(1)
   ELSE
      ! Unknown nitrogen uptake function
      print *,'STOP: unknown simINDynamics (',phytos(group)%simINDynamics,') for: ',phytos(group)%p_name
      ERROR STOP 1
   ENDIF

   ! Allocate a portion of N uptake to N fixation, where relevant:
   IF (phytos(group)%simNFixation /= 0) THEN
      a_nfix = phytos(group)%R_nfix * a_nfix * phy
      IF (a_nfix > ABS(uptake(1))) THEN
         ! Extreme case:
         a_nfix = -uptake(1)
         uptake(1) = zero_
      ELSE
         ! Reduce n-uptake by the amount fixed:
         IF ( uptake(1) /= 0. ) &
            uptake(1) = uptake(1) * (ABS(uptake(1))-a_nfix) / ABS(uptake(1))
      ENDIF
   ENDIF

   ! Disaggregate N sources to NO3, NH4, DON and N2, based on configuraiton
   PNf = phyto_pN(phytos,group,nh4up,no3up)

   IF (phytos(group)%simDINUptake /= 0) THEN
     !uptake(inh4) = uptake(1) * (1.0-PNf) !inh4 == 2
     !uptake(ino3) = uptake(1) * PNf       !ino3 == 1
      uptake(inh4) = uptake(1) * (PNf)     !inh4 == 2
      uptake(ino3) = uptake(1) * (1.-PNf)  !ino3 == 1
   ENDIF
   IF (phytos(group)%simDONUptake /= 0) THEN
      uptake(idon) = 0.0                   !MH to fix  (idon == 3)
   ENDIF
!   IF (phytos(group)%simNFixation /= 0 .AND. do_N2uptake) THEN
   IF (phytos(group)%simNFixation /= 0) THEN
      uptake(iN2) = -a_nfix                 ! iN2 == 4
   ENDIF

   ! Release of nitrogen due to excretion from phytoplankton and
   ! contribution of mortality and excretion OM:
   ! (/day +/day)* mg N/ mg C * mgC

   excretion = (respiration*phytos(group)%k_fdom + exudation)*theX_ncon
   mortality = respiration*(1.0-phytos(group)%k_fdom)*theX_ncon

   ! should check here e or m is not exceeding X_nmin
END SUBROUTINE phyto_internal_nitrogen
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION phyto_fN(phytos, group, IN, din, don) RESULT(fN)
!-------------------------------------------------------------------------------
! Nitrogen limitation of phytoplankton.
! Michaelis-Menton type formulation or droop model for species with IN
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(phyto_data_t),DIMENSION(:),INTENT(in)  :: phytos
   INTEGER,INTENT(in)                          :: group
   AED_REAL,INTENT(in),OPTIONAL                :: IN
   AED_REAL,INTENT(in),OPTIONAL                :: din
   AED_REAL,INTENT(in),OPTIONAL                :: don
!
!LOCALS
   AED_REAL :: fN
   AED_REAL :: nup

!-------------------------------------------------------------------------------
!BEGIN
   fN=one_

   IF (PRESENT(din) .OR. PRESENT(don)) THEN
     ! Calculate external nutrient limitation factor
     nup = 0.0
     IF (PRESENT(din) .AND. phytos(group)%simDINUptake == 1) THEN
       nup = nup + din
     ENDIF
     IF (PRESENT(don) .AND. phytos(group)%simDONUptake == 1) THEN
       nup = nup + don
     ENDIF
     ! FIX 2026-08-29: the denominator (nup-N_o+K_N) was unguarded, so when
     ! N_o-nup > K_N it turned negative and fN came out positive; the existing
     ! upper clamp below then pinned a starved cell to fN=1.0 (max growth)
     ! instead of 0. Guard the denominator as phyto_fP (~L533) already does.
     fN = (nup-phytos(group)%N_o) / &
           (phytos(group)%K_N + (MAX(zero_, (nup-phytos(group)%N_o))))
   ELSE
     ! Calculate internal nutrient limitation factor
     fN =   phytos(group)%X_nmax*(1.0-phytos(group)%X_nmin/IN) / &
            (phytos(group)%X_nmax-phytos(group)%X_nmin)
   ENDIF

   IF ( fN < zero_ ) fN = zero_
   IF ( fN > 1.000 ) fN = 1.000
END FUNCTION phyto_fN
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION phyto_fP(phytos, group, IP, frp) RESULT(fP)
!-------------------------------------------------------------------------------
! Phosphorus limitation of phytoplankton
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(phyto_data_t),DIMENSION(:),INTENT(in)  :: phytos
   INTEGER,INTENT(in)                          :: group
   AED_REAL,INTENT(in), OPTIONAL               :: IP
   AED_REAL,INTENT(in), OPTIONAL               :: frp
!
!LOCALS
   AED_REAL :: fP
!
!-------------------------------------------------------------------------------
!BEGIN
   fP=one_

   IF(PRESENT(frp)) THEN
     fP = (frp-phytos(group)%P_0) / &
             (phytos(group)%K_P + (MAX(zero_, (frp-phytos(group)%P_0))))
   ELSE
     fP = phytos(group)%X_pmax * (1.0 - phytos(group)%X_pmin/IP) / &
                          (phytos(group)%X_pmax-phytos(group)%X_pmin)
   ENDIF

   IF( fP < zero_ ) fP = zero_
   IF( fP > 1.000 ) fP = 1.000
END FUNCTION phyto_fP
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION phyto_fSi(phytos, group, Si) RESULT(fSi)
!-------------------------------------------------------------------------------
! Silica limitation (eg. for diatoms)
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(phyto_data_t),DIMENSION(:),INTENT(in)  :: phytos
   INTEGER,INTENT(in)                          :: group
   AED_REAL,INTENT(in)                         :: Si
!
!LOCALS
   AED_REAL :: fSi
!
!-------------------------------------------------------------------------------
!BEGIN
   fSi = one_

   IF (phytos(group)%simSiUptake == 1) THEN
     ! FIX 2026-08-29: the denominator (Si-Si_0+K_Si) was unguarded, so when
     ! Si_0-Si > K_Si it went negative and fSi came out POSITIVE and >1 (or
     ! infinite at the zero crossing) for a silica-starved cell - i.e. maximum
     ! growth from no silica. Mirror the guard already used by phyto_fP (~L533)
     ! and add the missing upper clamp, matching phyto_fN / phyto_fP.
     fSi = (Si-phytos(group)%Si_0) / &
           (phytos(group)%K_Si + (MAX(zero_, (Si-phytos(group)%Si_0))))
     IF ( fSi < zero_ ) fSi=zero_
     IF ( fSi > 1.000 ) fSi=1.000
   ENDIF
END FUNCTION phyto_fSi
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION phyto_pN(phytos,group,NH4,NO3) RESULT(pN)
!-------------------------------------------------------------------------------
! Calculates the relative preference of uptake by phytoplankton of
! ammonia uptake over nitrate.
!-------------------------------------------------------------------------------
   TYPE(phyto_data_t),DIMENSION(:),INTENT(in)  :: phytos
   INTEGER,INTENT(in)                          :: group
   AED_REAL,INTENT(in)                         :: NH4
   AED_REAL,INTENT(in)                         :: NO3
!
!LOCALS
   AED_REAL :: pN
!
!-------------------------------------------------------------------------------
!BEGIN
   pN = zero_

   IF (NH4 > 0.0) THEN
      pN = NH4*NO3 / ((NH4+phytos(group)%K_N)*(NO3+phytos(group)%K_N)) &
         + NH4*phytos(group)%K_N / ((NH4+NO3)*(NO3+phytos(group)%K_N))
   ENDIF
END FUNCTION phyto_pN
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION findMin(a1,a2,a3,a4) RESULT(theMin)
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in) :: a1,a2,a3,a4
!LOCALS
   AED_REAL     :: theMin
!
!-------------------------------------------------------------------------------
!BEGIN
   theMin = a1
   IF(a2 < theMin)      theMin = a2
   IF(a3 < theMin)      theMin = a3
   IF(a4 < theMin)      theMin = a4

   IF( theMin<zero_ )  theMin=zero_
END FUNCTION findMin
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION phyto_salinity(phytos,group,salinity) RESULT(fSal)
!-------------------------------------------------------------------------------
! Salinity tolerance of phytoplankton
!
! CAEDYM Implmentation based on Griffin et al (2001) and Robson and Hamilton
! (2004), plus Lassiter option also added from Zhu et al (2017).
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(phyto_data_t),DIMENSION(:),INTENT(in)  :: phytos
   INTEGER,INTENT(in)                          :: group
   AED_REAL,INTENT(in)                         :: salinity
!
!LOCALS
   AED_REAL :: fSal ! Returns the salinity function
   AED_REAL :: tmp1,tmp2,tmp3,fSa,fSb,fSc,fSo
   INTEGER  :: sal_model
!
!-------------------------------------------------------------------------------
!BEGIN
   fSal = 0.
   ! User can set salTol as negative (growth supressison) or positive
   ! (respiration enhancement). The equations are the same but will range from
   ! 0-1 for supression and >1 for enhancment. Users must ensure S_bep is set
   ! accordingly - Generally, S_bep>1 for salTol>0 and vice versa.
   sal_model = ABS(phytos(group)%salTol)

   IF (sal_model == 0) THEN
      fSal = one_
   ELSEIF (sal_model == 1) THEN
      !# f(S) = 1 at S=S_opt, f(S) = S_bep at S=S_maxsp.

      IF (phytos(group)%salTol<0 .and. phytos(group)%S_bep>1) &
        PRINT *,'WARNING: salTol flag for phyto group: ',group, &
                   ' is set for growth supression, but S_bep is >1: ', phytos(group)%S_bep

      ! change fSal above Sopt for increasing stress on freshwater species
      IF (salinity>phytos(group)%S_opt) THEN
        tmp1 = (phytos(group)%S_bep-1.0) / ((phytos(group)%S_maxsp - phytos(group)%S_opt)**2.0)
        tmp2 = (phytos(group)%S_bep-1.0) * 2.0*phytos(group)%S_opt / &
              ((phytos(group)%S_maxsp-phytos(group)%S_opt)**2.0)
        tmp3 = (phytos(group)%S_bep-1.0) * phytos(group)%S_opt*phytos(group)%S_opt / &
              ((phytos(group)%S_maxsp-phytos(group)%S_opt)**2.0) + 1.0

        fSal = tmp1*(salinity**2.0) - tmp2*salinity + tmp3
      ELSE
        fSal = 1.0
      ENDIF
   ELSEIF (sal_model == 2) THEN
      !# f(S) = 1 at S>=S_opt, f(S) = S_bep at S=0.

      IF (phytos(group)%salTol<0 .and. phytos(group)%S_bep>1) &
        PRINT *,'WARNING: salTol flag for phyto group: ',group, &
                      ' is set for growth supression, but S_bep is >1: ', phytos(group)%S_bep

      IF (salinity<phytos(group)%S_opt) THEN
         fSal = (phytos(group)%S_bep-1.0) * (salinity**2.0)/(phytos(group)%S_opt**2.0)   &
                - 2.0*(phytos(group)%S_bep-1.0)*salinity/phytos(group)%S_opt &
                + phytos(group)%S_bep
      ELSE
        fSal = 1.0
      ENDIF
   ELSEIF (sal_model == 3) THEN
      ! f(S) = 1 at S=S_opt, f(S) = S_bep at S=0 and 2*S_opt.

      IF (phytos(group)%salTol<0 .and. phytos(group)%S_bep>1) &
        PRINT *,'WARNING: salTol flag for phyto group: ',group, &
                      ' is set for growth supression, but S_bep is >1: ', phytos(group)%S_bep

      IF (salinity < phytos(group)%S_opt) THEN
         fSal = (phytos(group)%S_bep-1.0)*(salinity**2.0)/(phytos(group)%S_opt**2.0)-  &
                      2.0*(phytos(group)%S_bep-1.0)*salinity/phytos(group)%S_opt+phytos(group)%S_bep
      ENDIF
      IF ((salinity>phytos(group)%S_maxsp) .AND. (salinity<(phytos(group)%S_maxsp + phytos(group)%S_opt))) THEN
         fSal = (phytos(group)%S_bep - one_)*(phytos(group)%S_maxsp + phytos(group)%S_opt - salinity)**2  &
             / (phytos(group)%S_opt**2) -                                                                   &
             2 * (phytos(group)%S_bep - one_) * (phytos(group)%S_maxsp + phytos(group)%S_opt - salinity)  &
             / phytos(group)%S_opt + phytos(group)%S_bep
      ENDIF
      IF ((salinity >=  phytos(group)%S_opt) .AND. (salinity <= phytos(group)%S_maxsp) ) fSal = 1
      IF ( salinity >= (phytos(group)%S_maxsp + phytos(group)%S_opt) ) fSal = phytos(group)%S_bep
    ELSEIF (sal_model == 4) THEN
       ! Lassiter.
       ! This is used to control growth on species that like brackish water
       fSa = phytos(group)%S_bep
       fSb = 1.
       fSc = phytos(group)%S_maxsp
       fSo = phytos(group)%S_opt

       IF(salinity>fSc)THEN
         fSal = zero_
       ELSE
         fSal = fSb*EXP(fSa*(salinity-fSo))*((fSc-salinity)/(fSc-fSo))**(fSa*(fSc-fSo))
       ENDIF

       ! check if its used as resipration enhancement or growth supression
       IF (phytos(group)%salTol>0 ) fSal = (one_-fSal) + one_
   ELSE
      fSal = one_
      PRINT *,'WARNING: Unsupported salTol flag for phyto group: ',group,'=', phytos(group)%salTol
   ENDIF

   IF( fSal < zero_ ) fSal = zero_
END FUNCTION phyto_salinity
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION photosynthesis_irradiance(lightModel, I_K, I_S, par, extc, Io, dz) RESULT(fI)
!-------------------------------------------------------------------------------
! Light limitation of pytoplankton via various model approaches. Refer to
! overview presented in Table 1 of:
!
! Baklouti, M., Diaz, F., Pinazo, C., Faure, V., Quéguiner, B., 2006.
!  Investigation of mechanistic formulations depicting phytoplankton dynamics
!    for models of marine pelagic ecosystems and description of a new model.
!  Progress in Oceanography 71 (1), 1-33.
!
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER, INTENT(in) :: lightModel
   AED_REAL,INTENT(in) :: I_K
   AED_REAL,INTENT(in) :: I_S
   AED_REAL,INTENT(in) :: par
   AED_REAL,INTENT(in) :: extc
   AED_REAL,INTENT(in) :: Io
   AED_REAL,INTENT(in) :: dz
   AED_REAL            :: fI !-- Returns the light limitation
!
!CONSTANTS
   AED_REAL,PARAMETER  :: one_e_neg3 = 1e-3
!
!LOCALS
   AED_REAL,PARAMETER  :: A = 5.0, eps = 0.5
   AED_REAL            :: par_t,par_b,par_c
   AED_REAL            :: z1,z2,x
   LOGICAL,SAVE        :: warned_lightModel = .FALSE.  ! 2026-08-29 : warn-once
!
!-------------------------------------------------------------------------------
!BEGIN
   fI    = 0.0
   IF (Io == zero_) RETURN

   ! Light at different points within the layer (MH: to be replaced by function)
   par_t = par
   par_b = par_t * EXP( -extc * dz )
   par_c = par_t * EXP( -extc * dz/2. )

   SELECT CASE (lightModel)
      CASE ( 0 )
         ! Light limitation without photoinhibition.
         !   This is the Webb et al (1974) model solved using the numerical
         !   integration approach as in CAEDYM (Hipsey and Hamilton, 2008)

         z1 = -par_t / I_K
         z2 = -par_b / I_K

         z1 = exp_integral(z1)
         z2 = exp_integral(z2)

         fI = 1.0 + (z2 - z1) / MAX(extc * dz,one_e_neg3)

         IF (par_t < 5e-5 .OR. fI < 5e-5) fI = 0.0        ! A simple check

      CASE ( 1 )
         ! Light limitation without photoinhibition.
         ! This is the Monod (1950) model.

         x = par_c/I_K
         fI = x / (one_ + x)

      CASE ( 2 )
         ! Light limitation with photoinhibition.
         ! This is the Steele (1962) model.

         x = par_c/I_S
         fI = x * EXP(one_ - x)
         IF (par_t < 5e-5 .OR. fI < 5e-5) fI = 0.0

      CASE ( 3 )
         ! Light limitation without photoinhibition.
         ! This is the Webb et al. (1974) model.

         x = par_c/I_K
         fI = one_ - EXP(-x)

      CASE ( 4 )
         ! Light limitation without photoinhibition.
         ! This is the Jassby and Platt (1976) model.

         x = par_c/I_K
         fI = TANH(x)

      CASE ( 5 )
         ! Light limitation without photoinhibition.
         ! This is the Chalker (1980) model.

         x = par_c/I_K
         fI = (EXP(x * (one_ + eps)) - one_) / &
              (EXP(x * (one_ + eps)) + eps)

      CASE ( 6 )
         ! Light limitation with photoinhibition.
         ! This is the Klepper et al. (1988) / Ebenhoh et al. (1997) model.
         x = par_c/I_S
         fI = ((2.0 + A) * x) / ( one_ + (A * x) + (x * x) )

      CASE ( 7 )
         ! Light limitation with photoinhibition.
         ! This is an integrated form of Steele model.

         fI = ( EXP(1-par_b/I_S) - &
                EXP(1-par_t/I_S)   ) / (extc * dz)

      CASE ( 10 )
         ! Light limitation without photoinhibition.
         ! This is the Webb et al. (1974) model.

         x = par_b/I_K          ! Uses BOTTOM light
         fI = one_ - EXP(-x)

      CASE ( 11 )
        ! Light limitation without photoinhibition.
        ! This is the Webb et al. (1974) model.

        x = par_t/I_K           ! Uses SURFACE light
        fI = one_ - EXP(-x)

      CASE DEFAULT
        ! FIX 2026-08-29: there was no CASE DEFAULT. An out-of-range lightModel
        ! (e.g. a typo, or 8/9 which are not implemented) fell through the
        ! SELECT leaving fI at its initialised 0.0, so the group had zero light
        ! limitation forever with no message. Warn once and fall back to the
        ! CASE(0) Webb et al. (1974) integrated form.
        IF ( .NOT. warned_lightModel ) THEN
          PRINT *,'WARNING: unsupported lightModel = ',lightModel, &
                  ' ; falling back to lightModel=0 (Webb et al. 1974, integrated)'
          warned_lightModel = .TRUE.
        ENDIF

        z1 = -par_t / I_K
        z2 = -par_b / I_K

        z1 = exp_integral(z1)
        z2 = exp_integral(z2)

        fI = 1.0 + (z2 - z1) / MAX(extc * dz,one_e_neg3)

        IF (par_t < 5e-5 .OR. fI < 5e-5) fI = 0.0        ! A simple check

  END SELECT

   IF ( fI < zero_ ) fI = zero_
END FUNCTION photosynthesis_irradiance
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bio_respiration(R_resp,theta_resp,temp) RESULT(respiration)
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in) :: R_resp
   AED_REAL,INTENT(in) :: theta_resp
   AED_REAL,INTENT(in) :: temp
!
!LOCALS
   AED_REAL :: respiration ! Returns the phytoplankton respiration.
!
!-------------------------------------------------------------------------------
!BEGIN
   respiration = R_resp * theta_resp**(temp-20.0)
END FUNCTION bio_respiration
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE aed_bio_utils
