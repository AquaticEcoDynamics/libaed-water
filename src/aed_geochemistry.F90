!###############################################################################
!#                                                                             #
!# aed_geochemistry.F90                                                        #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2012-2026 : The University of Western Australia                  #
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
!#                                                                             #
!###############################################################################

#include "aed.h"


MODULE aed_geochemistry
!-------------------------------------------------------------------------------
! aed_geochemistry --- geochemistry model
!
! The AED module geochemistry contains equations that describe ....
!-------------------------------------------------------------------------------
   USE aed_core

   USE aed_gclib,ONLY : AED_GC_Input,printTables
   USE aed_gcsolver,ONLY : ConfigEquilibriumSolver, &
                            GetListOfGeochemDiagnostics, &
                            InitialiseGCProperties, &
                            UpdateEquilibration, &
                            returnGCDerivedVector, &
                            simManganRedox, &
                            simArsenicRedox, &
                            simSulfurRedox, &
                            simCarbonRedox, &
                            simIronRedox


   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed_geochemistry_data_t
!
   TYPE,extends(aed_model_data_t) :: aed_geochemistry_data_t
      !# Variable identifiers
      INTEGER  :: id_comp(MAX_GC_COMPONENTS), id_mins(MAX_GC_MINERALS)
      INTEGER  :: id_cdep(MAX_GC_COMPONENTS), id_mdep(MAX_GC_MINERALS)
      INTEGER  :: id_ubalchg ,id_pH, id_c_pco2, id_o_oxy
      INTEGER  :: id_temp,id_sal
      INTEGER  :: id_sed_dic
      INTEGER  :: id_gcdiag(MAX_GC_COMPONENTS), id_noncon
      INTEGER  :: id_totC
      INTEGER  :: id_feii, id_feiii
      INTEGER  :: id_h2s, id_so4
      INTEGER  :: id_tss, id_compd(MAX_GC_COMPONENTS)

      !# Model parameters
      INTEGER  :: num_comp, num_mins
      INTEGER  :: inflow_pH_update
      LOGICAL  :: component_linked(MAX_GC_COMPONENTS),mineral_linked(MAX_GC_MINERALS)
      LOGICAL  :: simEq
      AED_REAL :: Riron_red, theta_iron_red, Kiron_red
      AED_REAL :: Riron_aox, Riron_box, theta_iron_ox, Kiron_ox
      AED_REAL :: Rsulf_red, theta_sulf_red, Ksulf_red
      AED_REAL :: Rsulf_ox, theta_sulf_ox, Ksulf_ox
      INTEGER  :: iron_ox_model
      AED_REAL :: speciation_dt
      AED_REAL :: Fsed_gch(MAX_GC_COMPONENTS)
      AED_REAL :: Ksed_gch_o2(MAX_GC_COMPONENTS)
      AED_REAL :: Ksed_gch_pH(MAX_GC_COMPONENTS)
      AED_REAL :: w_gch(MAX_GC_MINERALS)
      AED_REAL,          DIMENSION(:), ALLOCATABLE :: DissComp,PartComp
      CHARACTER(len=64), DIMENSION(:), ALLOCATABLE :: listDissTransVars,listPartTransVars

      INTEGER  :: MeAdsorptionModel(MAX_GC_COMPONENTS), resuspension
      LOGICAL  :: simMeAdsorption, ads_use_pH, ads_use_external_tss
      AED_REAL :: KMep(MAX_GC_COMPONENTS), Kadsratio(MAX_GC_COMPONENTS), Qmax(MAX_GC_COMPONENTS), theta_KMe, K_sal   ! Adsoprtion
      AED_REAL :: w_gch_ads(MAX_GC_COMPONENTS)                                 ! Sorbed pool settling
      LOGICAL  :: link_ads_vvel(MAX_GC_COMPONENTS)                             ! adopt sorbent particle vvel
      INTEGER  :: id_compd_vvel                                                ! sorbent _vvel diagnostic
      INTEGER  :: id_compd_frc(MAX_GC_COMPONENTS)                              ! sorbed fraction diagnostic
      INTEGER  :: id_compd_set(MAX_GC_COMPONENTS)                              ! sorbed settling flux diagnostic
      INTEGER  :: id_compd_srp(MAX_GC_COMPONENTS)                              ! kinetic sorption rate diagnostic
      INTEGER  :: id_dz
      AED_REAL :: Rsorp_gch(MAX_GC_COMPONENTS)                                 ! kinetic sorption rate (0 = equilibrium)


     CONTAINS
         PROCEDURE :: define            => aed_define_geochemistry
         PROCEDURE :: initialize        => aed_initialize_geochemistry
         PROCEDURE :: calculate         => aed_calculate_geochemistry
         PROCEDURE :: calculate_benthic => aed_calculate_benthic_geochemistry
         PROCEDURE :: equilibrate       => aed_equilibrate_geochemistry
         PROCEDURE :: mobility          => aed_mobility_geochemistry
!        PROCEDURE :: light_extinction  => aed_light_extinction_geochemistry
         PROCEDURE :: inflow_update     => aed_inflow_update_geochemistry
!        PROCEDURE :: delete            => aed_delete_geochemistry

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
SUBROUTINE aed_define_geochemistry(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and the variables exported
!  by the model are registered with AED.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_geochemistry_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst
!
!LOCALS
   INTEGER           :: status, i

!  %% NAMELIST   %%  /aed_geochemistry/
!  %% Last Checked 20/08/2021
   INTEGER           :: speciation_dt
   INTEGER           :: num_components,num_minerals
   INTEGER           :: nDissTransportables, nPartTransportables
   INTEGER           :: inflow_pH_update = 0
   LOGICAL           :: simEq = .TRUE.
   AED_REAL          :: min
   AED_REAL          :: dis_initial(MAX_GC_COMPONENTS) = 0.0
   AED_REAL          :: Fsed_gch(MAX_GC_COMPONENTS) = 0.0
   AED_REAL          :: Ksed_gch_o2(MAX_GC_COMPONENTS) = 0.0  ! 0 => no O2 effect on Fsed_gch
   AED_REAL          :: Ksed_gch_pH(MAX_GC_COMPONENTS) = 0.0  ! 0 => no pH effect on Fsed_gch
   AED_REAL          :: min_initial(MAX_GC_MINERALS) = 0.0
   AED_REAL          :: w_gch(MAX_GC_MINERALS)
   AED_REAL          :: pH_initial = 7.5
   AED_REAL          :: pe_initial = 8.0
   !-- Kinetic redox parameters. These are only referenced when the matching
   !-- couple is active (see setupWCRedoxConfiguration), but they must still be
   !-- initialised here: without a default an omitted namelist entry leaves the
   !-- variable undefined rather than zero.
   AED_REAL          :: Riron_red      = zero_  ! /day  FeIII -> FeII (biotic)
   AED_REAL          :: theta_iron_red = one_   ! 1.0 => no T correction
   AED_REAL          :: Kiron_red      = 100.0  ! mmol/m3 O2 inhibition of redn
   AED_REAL          :: Riron_aox      = zero_  ! /day  abiotic FeII oxidation
   AED_REAL          :: Riron_box      = zero_  ! /day  biotic  FeII oxidation
   AED_REAL          :: theta_iron_ox  = one_   ! 1.0 => no T correction
   AED_REAL          :: Kiron_ox       = 100.0  ! mmol/m3 O2 half-sat of oxidn
   AED_REAL          :: Rsulf_red      = zero_  ! /day  SO4 -> H2S
   AED_REAL          :: theta_sulf_red = one_   ! 1.0 => no T correction
   AED_REAL          :: Ksulf_red      = 100.0  ! mmol/m3 O2 inhibition of redn
   AED_REAL          :: Rsulf_ox       = zero_  ! /day  H2S -> SO4
   AED_REAL          :: theta_sulf_ox  = one_   ! 1.0 => no T correction
   AED_REAL          :: Ksulf_ox       = 100.0  ! mmol/m3 O2 half-sat of oxidn
   !-- 0 = biotic only, pH independent (default, legacy behaviour)
   !-- 1 = Singer-Stumm speciation based abiotic oxidation, plus the biotic term
   INTEGER           :: iron_ox_model  = 0
   CHARACTER(len=64) :: geochem_file = ''
   CHARACTER(len=64) :: dis_components(MAX_GC_COMPONENTS) = ''
   CHARACTER(len=64) :: component_link(MAX_GC_COMPONENTS) = ''
   CHARACTER(len=64) :: speciesOutput(10) = ''
   CHARACTER(len=64) :: the_minerals(MAX_GC_MINERALS) = ''
   CHARACTER(len=64) :: mineral_link(MAX_GC_MINERALS) = ''
   CHARACTER(len=64) :: ph_link = 'CAR_pH'
   CHARACTER(len=64) :: pco2_link = 'CAR_pCO2'
   CHARACTER(len=64) :: oxy_link = 'OXY_oxy'
   ! Adsorption
   LOGICAL           :: simMeAdsorption = .FALSE.
   LOGICAL           :: ads_use_external_tss = .FALSE.
   INTEGER           :: MeAdsorptionModel(MAX_GC_COMPONENTS) = 0
   LOGICAL           :: ads_use_pH   = .FALSE.
   AED_REAL          :: KMep(MAX_GC_COMPONENTS) = 1.05
   AED_REAL          :: theta_KMe   = 1.0    ! 1.0 => no T correction of KMep
   AED_REAL          :: K_sal       = zero_  ! 0   => no salinity correction of KMep
   AED_REAL          :: Kadsratio(MAX_GC_COMPONENTS)    = 1.05
   AED_REAL          :: Qmax(MAX_GC_COMPONENTS)         = 1.05
   AED_REAL          :: w_gch_ads(MAX_GC_COMPONENTS)    = zero_  ! m/day; -ve = settling
   AED_REAL          :: Rsorp_gch(MAX_GC_COMPONENTS)    = zero_  ! /day relaxation to sorption
                                                                 ! equilibrium; 0 => instantaneous
   CHARACTER(len=64) :: sorption_target_variable=''

! %% From Module Globals
!  INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
!                                             ! 1 = basic diagnostic outputs
!                                             ! 2 = flux rates, and supporitng
!                                             ! 3 = other metrics
!                                             !10 = all debug & checking outputs
!  %% END NAMELIST   %%  /aed_geochemistry/

   CHARACTER(len=64), DIMENSION(:), ALLOCATABLE :: diagnosticList

   NAMELIST /aed_geochemistry/ speciation_dt, geochem_file,                    &
                    num_components, dis_components, component_link, Fsed_gch,  &
                    Ksed_gch_o2, Ksed_gch_pH,                                  &
                    dis_initial, num_minerals, the_minerals, mineral_link,     &
                    w_gch, min_initial, pH_initial, speciesOutput, simEq,      &
                    Riron_red, theta_iron_red, Kiron_red,                      &
                    Riron_aox, Riron_box, theta_iron_ox, Kiron_ox,             &
                    iron_ox_model,                                             &
                    Rsulf_red, theta_sulf_red, Ksulf_red,                      &
                    Rsulf_ox, theta_sulf_ox, Ksulf_ox,                         &
                    ph_link, inflow_pH_update, pco2_link, diag_level, &
                    simMeAdsorption, &
                    ads_use_external_tss, &
                    sorption_target_variable, &
                    MeAdsorptionModel, &
                    KMep, &
                    theta_KMe, &
                    K_sal, &
                    ads_use_pH, &
                    Kadsratio, &
                    Qmax, &
                    w_gch_ads, &
                    Rsorp_gch
!-------------------------------------------------------------------------------
!BEGIN
   print *,"        aed_geochemistry configuration"
   print *,"         WARNING! aed_geochemistry model is under development"


   ! MH:JOBS
   ! remove pe
   ! species outputs

   !----------------------------------------------------------------------------
   ! Initialise variables
   component_link(:) = ''
   data%component_linked(:) = .FALSE.
   data%mineral_linked(:) = .FALSE.
   data%id_compd(:) = -1
   data%link_ads_vvel(:) = .FALSE.
   data%w_gch_ads(:) = zero_
   data%id_compd_frc(:) = -1
   data%id_compd_set(:) = -1
   data%id_compd_srp(:) = -1
   data%Rsorp_gch(:) = zero_

! Initialisation now done in declaration
!  dis_initial = 0.0  ! default, overwritten by namelist
!  min_initial = 0.0  ! default, overwritten by namelist
!  pH_initial  = 7.5  ! default, overwritten by namelist
!  pe_initial  = 8.0  ! default, not used

!  ph_link     = 'CAR_pH'
!  pco2_link   = 'CAR_pCO2'

   !----------------------------------------------------------------------------
   ! Read the namelist
   read(namlst,nml=aed_geochemistry,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_geochemistry'

   data%simEq = simEq
   data%speciation_dt = speciation_dt  ! Note this is now managed in FV_AED or GLM_AED

   data%Riron_red = Riron_red             ; data%Kiron_red= Kiron_red
   data%theta_iron_red= theta_iron_red    ; data%theta_iron_ox=theta_iron_ox
   data%Riron_aox=Riron_aox/ secs_per_day ; data%Riron_box=Riron_box/ secs_per_day
   data%Kiron_ox = Kiron_ox               ; data%iron_ox_model = iron_ox_model

   data%Rsulf_red = Rsulf_red             ; data%Ksulf_red= Ksulf_red
   data%theta_sulf_red= theta_sulf_red    ; data%theta_sulf_ox=theta_sulf_ox
   data%Rsulf_ox=Rsulf_ox / secs_per_day

   data%inflow_pH_update = inflow_pH_update

   speciesOutput = ''
   speciesOutput(1) = 'NONCON'
   !speciesOutput(1) = 'HCO3-'

   data%simMeAdsorption =  simMeAdsorption
   data%MeAdsorptionModel = MeAdsorptionModel
   data%KMep = KMep
   data%ads_use_pH = ads_use_pH
   data%ads_use_external_tss = ads_use_external_tss
   data%theta_KMe = theta_KMe
   data%K_sal = K_sal
   data%Kadsratio = Kadsratio
   data%Qmax = Qmax

   print *,"simMeAdsorption",simMeAdsorption
   print *,"MeAdsorptionModel",MeAdsorptionModel

   !----------------------------------------------------------------------------
   ! Now load the geochem database

   !CALL aed_geochem_load_params(data, geochem_file, modelinfo)
   CALL AED_GC_Input(geochem_file)
   !CALL printTables()

   !----------------------------------------------------------------------------
   ! Configure the geochemical solver
   CALL ConfigEquilibriumSolver( num_components,  num_minerals,                &
                                 dis_components(1:num_components),             &
                                 the_minerals(1:num_minerals),                 &
                                 nDissTransportables, nPartTransportables,     &
                                 data%listDissTransVars,data%listPartTransVars)

   data%num_comp = nDissTransportables
   data%num_mins = nPartTransportables

   !-- Report which kinetic redox couples the component set has activated. Each
   !-- requires BOTH of its redox states to be configured components, and is
   !-- otherwise skipped entirely - previously without saying so, which made a
   !-- missing member look like a rate-constant problem.
   print *,"        aed_geochemistry kinetic redox couples:"
   CALL reportRedoxCouple("iron     (FeII/FeIII)", simIronRedox)
   CALL reportRedoxCouple("sulfur   (H2S/SO4)   ", simSulfurRedox)
   CALL reportRedoxCouple("mangan   (MnII/MnIV) ", simManganRedox)
   CALL reportRedoxCouple("arsenic  (AsIII/AsV) ", simArsenicRedox)
   CALL reportRedoxCouple("carbon   (CH4/CO3)   ", simCarbonRedox)
   IF ( simIronRedox .AND. iron_ox_model == 1 ) THEN
      print *,"          FeII oxidation: Singer-Stumm speciation (pH dependent)"
      print *,"          WARNING! iron_ox_model=1 is not yet mass balance verified"
   ELSEIF ( simIronRedox ) THEN
      print *,"          FeII oxidation: biotic Monod only (pH INDEPENDENT)"
   ENDIF
   IF ( .NOT. (simIronRedox .OR. simSulfurRedox .OR. simManganRedox .OR.       &
               simArsenicRedox .OR. simCarbonRedox) ) THEN
      print *,"          none active: geochemistry is equilibrium-only"
   ENDIF

   !----------------------------------------------------------------------------
   ! Initialise the module level geochemical values ready for the registration

   ALLOCATE(data%DissComp(nDissTransportables))
   ALLOCATE(data%PartComp(nPartTransportables))

   data%DissComp = zero_
   DO i=1,nDissTransportables
     data%DissComp(i) = dis_initial(i)
     data%Fsed_gch(i) = Fsed_gch(i) / secs_per_day
     data%Ksed_gch_o2(i) = Ksed_gch_o2(i)
     data%Ksed_gch_pH(i) = Ksed_gch_pH(i)
     IF ( w_gch_ads(i) < -999. ) THEN
        ! sentinel: adopt the sorbent particle's settling velocity (_vvel link)
        data%link_ads_vvel(i) = .TRUE.
        data%w_gch_ads(i) = zero_
     ELSE
        data%link_ads_vvel(i) = .FALSE.
        data%w_gch_ads(i) = w_gch_ads(i) / secs_per_day
     ENDIF
     data%Rsorp_gch(i) = Rsorp_gch(i) / secs_per_day
   END DO
   component_link(num_components+1) = ph_link  ! Special pH var
   data%DissComp(num_components+1) = pH_initial
   data%DissComp(num_components+2) = pe_initial

   data%PartComp = zero_
   DO i=1,num_minerals
     data%PartComp(i) = min_initial(i)
     data%w_gch(i) = w_gch(i) / secs_per_day
   END DO

   CALL InitialiseGCProperties(data%DissComp, data%PartComp, 2)

   !print *,'data%DissComp',data%DissComp
   !print *,'data%PartComp',data%PartComp

   !----------------------------------------------------------------------------
   ! Process dis components adding as state vars or dependancies as appropriate
   DO i=1,nDissTransportables
      !print *,'i,',i,component_link(i)
      IF ( component_link(i) .EQ. '' ) THEN
         min = zero_
         IF ( TRIM(data%listDissTransVars(i)) == 'ubalchg' ) min=nan_
         ! Register state variables
         data%id_comp(i) = aed_define_variable(                                &
                                    TRIM(data%listDissTransVars(i)),           &
                                    'mmol/m**3','geochemistry',                &
                                    data%DissComp(i),                          &
                                    minimum=min)
         IF ( simMeAdsorption .and. MeAdsorptionModel(i)>0 ) THEN
            data%id_compd(i) = aed_define_variable(                            &
                                      TRIM(data%listDissTransVars(i))//'_ads', &
                                      'mmol/m**3','geochemistry',              &
                                      zero_,                                   &
                                      minimum=min,                             &
                                      mobility=data%w_gch_ads(i))
            IF ( diag_level>0 ) THEN
               data%id_compd_frc(i) = aed_define_diag_variable(                &
                                 TRIM(data%listDissTransVars(i))//'_ads_frc',  &
                                 '-','sorbed fraction of total')
               data%id_compd_set(i) = aed_define_diag_variable(                &
                                 TRIM(data%listDissTransVars(i))//'_ads_set',  &
                                 'mmol/m**3/d','sorbed pool settling flux')
               IF ( Rsorp_gch(i) > zero_ )                                     &
                  data%id_compd_srp(i) = aed_define_diag_variable(             &
                                 TRIM(data%listDissTransVars(i))//'_ads_srp',  &
                                 'mmol/m**3/d','kinetic sorption rate')
            ENDIF
         ENDIF

      ELSE
         ! Register external state variable dependencies
         data%id_cdep(i) = aed_locate_variable( TRIM(component_link(i)) )
         data%component_linked(i) = .true.
         IF ( simMeAdsorption .AND. MeAdsorptionModel(i)>0 )                   &
            print *,'        WARNING: MeAdsorptionModel is set for linked ',   &
                    'component ',TRIM(data%listDissTransVars(i)),              &
                    ' - sorption of linked components is not supported'
      ENDIF
      IF( TRIM(data%listDissTransVars(i)) .EQ. 'FeII' ) data%id_feii = data%id_comp(i)
      IF( TRIM(data%listDissTransVars(i)) .EQ. 'FeIII') data%id_feiii= data%id_comp(i)
      IF( TRIM(data%listDissTransVars(i)) .EQ. 'H2S' ) data%id_h2s   = data%id_comp(i)
      IF( TRIM(data%listDissTransVars(i)) .EQ. 'SO4') data%id_so4    = data%id_comp(i)
   ENDDO

   !----------------------------------------------------------------------------
   ! Process minerals adding as state vars or dependancies as appropriate
   DO i=1,num_minerals
      IF ( mineral_link(i) .EQ. '' ) THEN
         ! Register state variables
         data%id_mins(i) = aed_define_variable(                               &
                                    TRIM(the_minerals(i)),                     &
                                    'mmol/m**3','geochemistry',                &
                                    data%PartComp(i),                          &
                                    minimum=zero_,                             &
                                    mobility=data%w_gch(i))
      ELSE
         ! Register external state variable dependencies
         data%id_mdep(i) = aed_locate_variable( mineral_link(i))
         data%mineral_linked(i) = .true.
      ENDIF
   ENDDO


   !----------------------------------------------------------------------------
   ! Register links to other modules

   data%id_o_oxy = aed_locate_variable(oxy_link)

   IF ( ph_link .EQ. '' ) THEN
     ! Register as a diagnostic variable
     data%id_pH = aed_define_diag_variable( 'pH', 'pH', 'pH')
   ELSE
     ! Link to module
     data%id_pH = aed_locate_variable(ph_link)
   ENDIF

   data%id_compd_vvel = -1
   IF ( simMeAdsorption ) THEN
      IF ( ads_use_external_tss ) THEN
         print *,'        Me adsorption is configured to use the host TSS'
         data%id_tss = aed_locate_global('tss')
         IF ( ANY(data%link_ads_vvel) ) THEN
            print *,'  ERROR w_gch_ads vvel link requires an internal ',      &
                    'sorption_target_variable, not external TSS'
            STOP
         ENDIF
      ELSEIF ( sorption_target_variable .NE. '' ) THEN
         print *,'        Me is adsorbing to ',TRIM(sorption_target_variable)
         data%id_tss = aed_locate_variable(sorption_target_variable)
         IF ( ANY(data%link_ads_vvel) ) THEN
            print *,'        _ads pools adopt settling of ',                  &
                    TRIM(sorption_target_variable)//'_vvel'
            data%id_compd_vvel =                                              &
                aed_locate_variable(TRIM(sorption_target_variable)//'_vvel')
         ENDIF
      ELSE
         print *,'  ERROR Me adsorption is configured but no internal or ',   &
                 'external sorption target is set'
         STOP
      ENDIF
   ENDIF

   ! solution to get aed_carbon's pCO2 updated for atm exchange ...
   data%id_c_pco2 = aed_locate_variable(pco2_link)

   !----------------------------------------------------------------------------
   ! Register diagnostic variables

   CALL GetListOfGeochemDiagnostics(speciesOutput,diagnosticList)

   DO i=1,SIZE(diagnosticList)
     data%id_gcdiag(i) = aed_define_diag_variable( diagnosticList(i), &
                         '?mmol/m**3?', 'Geochemistry Diagnostic')
   END DO

   data%id_noncon = aed_define_diag_variable( 'noncon_mh', &
                         '-', 'non-convergence flag (0=converged, 1=failed)')

   !----------------------------------------------------------------------------

   ! Register environmental dependencies
   data%id_temp = aed_locate_global( 'temperature' )
   data%id_sal = aed_locate_global( 'salinity' )
   data%id_dz = aed_locate_global( 'layer_ht' )

   !----------------------------------------------------------------------------

END SUBROUTINE aed_define_geochemistry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




!###############################################################################
SUBROUTINE aed_initialize_geochemistry(data, column, layer_idx)
!-------------------------------------------------------------------------------
! Routine to update the dynamics of "Acid Sulfate Soils" (ASS) and determine   !
! the flux to the water column from exposed or re-wetted sediment              !
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_geochemistry_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL :: temp, tss

   ! State
   AED_REAL,   DIMENSION(SIZE(data%DissComp))  :: dissConcs
   AED_REAL,   DIMENSION(SIZE(data%PartComp))  :: partConcs
   ! Temporary variables
   INTEGER  :: i
!-------------------------------------------------------------------------------
!BEGIN

   !-- Retrieve current environmental conditions for the cell.
   temp = _STATE_VAR_(data%id_temp) ! local temperature

   !-- Retrieve current (local) state variable values into array for the gcsolver
   DO i=1,data%num_comp
      IF (.NOT.data%component_linked(i)) THEN
          dissConcs(i) = _STATE_VAR_(data%id_comp(i))
      ELSE
          dissConcs(i) = _STATE_VAR_(data%id_cdep(i))
      ENDIF
   ENDDO
   DO i=1,data%num_mins
      IF (.NOT.data%mineral_linked(i)) THEN
          partConcs(i) = _STATE_VAR_(data%id_mins(i))
      ELSE
          partConcs(i) = _STATE_VAR_(data%id_mdep(i))
      ENDIF
   ENDDO

   !-- Redo geochemical equilibration, now spatial initialisation is done
   CALL InitialiseGCProperties(dissConcs, partConcs, 2, inTemp=REAL(temp))

   !-- Copy back into main AED arrays
   DO i=1,data%num_comp
      IF (.NOT.data%component_linked(i)) THEN
         _STATE_VAR_(data%id_comp(i)) =  dissConcs(i)
      ELSE
         _STATE_VAR_(data%id_cdep(i)) =  dissConcs(i)
      ENDIF
   ENDDO
   DO i=1,data%num_mins
      IF (.NOT.data%mineral_linked(i)) THEN
         _STATE_VAR_(data%id_mins(i)) =  partConcs(i)
      ELSE
         _STATE_VAR_(data%id_mdep(i)) =  partConcs(i)
      ENDIF
   ENDDO


END SUBROUTINE aed_initialize_geochemistry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed_calculate_geochemistry(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of aed_geochemistry model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_geochemistry_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL           :: reduction,oxidation
   AED_REAL           :: feii,feiii,h2s,so4,oxy,temp
   AED_REAL           :: MeDis,MePar,sorption
   INTEGER            :: i
   AED_REAL,PARAMETER :: one_e_neg_ten = 1e-10
!-------------------------------------------------------------------------------
!BEGIN

   reduction = zero_
   oxidation = zero_

   !-- 1. Iron  ---------------------------------------------------------------!
   IF (simIronRedox) THEN

       oxy   = _STATE_VAR_(data%id_o_oxy) ! oxygen
       temp  = _STATE_VAR_(data%id_temp)
       feii  = _STATE_VAR_(data%id_feii)
       feiii = _STATE_VAR_(data%id_feiii)

       !-- Reduction
       reduction = calcIronReduction(data,feiii,oxy,temp)

       !-- Oxidation
       oxidation = calcIronOxidation(data,feii,oxy,temp)

       !-- Update FeII, O2, FeIII
       _FLUX_VAR_(data%id_feii) = _FLUX_VAR_(data%id_feii) + reduction - oxidation

       !-- 4Fe++ + O2 + 4H+ -> 4Fe+++ + 2H2O : 1/4 mol O2 per mol Fe(II) oxidised
       _FLUX_VAR_(data%id_o_oxy) = _FLUX_VAR_(data%id_o_oxy) - oxidation *0.25

       _FLUX_VAR_(data%id_feiii) = _FLUX_VAR_(data%id_feiii) + oxidation - reduction

   END IF

   reduction = zero_
   oxidation = zero_
   !-- 2. Sulfur --------------------------------------------------------------!
   IF (simSulfurRedox) THEN

      oxy   = _STATE_VAR_(data%id_o_oxy) ! oxygen
      temp  = _STATE_VAR_(data%id_temp)
      h2s   = _STATE_VAR_(data%id_h2s)
      so4   = _STATE_VAR_(data%id_so4)

      !-- Reduction
      reduction = calcSulfurReduction(data,so4,oxy,temp)

      !-- Oxidation
      oxidation = calcSulfurOxidation(data,h2s,oxy,temp)

      !-- Update
      _FLUX_VAR_(data%id_h2s) = _FLUX_VAR_(data%id_h2s) + reduction - oxidation

      !-- HS- + 2O2 -> SO4-- + H+ : 2 mol O2 consumed per mol S oxidised
      _FLUX_VAR_(data%id_o_oxy) = _FLUX_VAR_(data%id_o_oxy) - oxidation * 2.0

      _FLUX_VAR_(data%id_so4) = _FLUX_VAR_(data%id_so4) + oxidation - reduction

   END IF

   !-- 3. Sorption kinetics ----------------------------------------------------!
   !   Components with Rsorp_gch>0 relax first-order toward the equilibrium
   !   partition; components with Rsorp_gch=0 are re-partitioned instantly in
   !   aed_equilibrate_geochemistry instead. NB with the Euler host solution
   !   keep Rsorp_gch*dt <= 1 or the relaxation can overshoot.
   IF (data%simMeAdsorption) THEN
      DO i=1,data%num_comp
         IF ( data%MeAdsorptionModel(i) == 0 ) CYCLE
         IF ( data%component_linked(i) ) CYCLE
         IF ( data%id_compd(i) <= 0 ) CYCLE
         IF ( data%Rsorp_gch(i) <= zero_ ) CYCLE

         CALL calcSorptionPartition(data,column,layer_idx,i,MeDis,MePar)

         ! relax the dissolved pool toward its equilibrium value
         ! (+ve = adsorbing, -ve = desorbing)
         sorption = data%Rsorp_gch(i) *                                        &
                            ( _STATE_VAR_(data%id_comp(i)) - MeDis )

         _FLUX_VAR_(data%id_comp(i))  = _FLUX_VAR_(data%id_comp(i))  - sorption
         _FLUX_VAR_(data%id_compd(i)) = _FLUX_VAR_(data%id_compd(i)) + sorption

         IF ( data%id_compd_srp(i) > 0 ) &
            _DIAG_VAR_(data%id_compd_srp(i)) = sorption * secs_per_day
         IF ( data%id_compd_frc(i) > 0 ) &
            _DIAG_VAR_(data%id_compd_frc(i)) = _STATE_VAR_(data%id_compd(i)) / &
                   MAX( _STATE_VAR_(data%id_comp(i))                           &
                       +_STATE_VAR_(data%id_compd(i)), one_e_neg_ten )
      ENDDO
   ENDIF

END SUBROUTINE aed_calculate_geochemistry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_benthic_geochemistry(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED
! geochemistry. Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_geochemistry_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
     AED_REAL :: temp,ph
   ! State
     AED_REAL :: dic,oxy
   ! Parameters
     AED_REAL, PARAMETER :: KpHS = 2.0   ! about half of maximum flux at pH5
     AED_REAL, PARAMETER :: KDOs = 2.0*1e3/16.
   ! Temporary variables
     AED_REAL :: gch_flux,oxyEffect,pHEffect
     INTEGER :: i

!-------------------------------------------------------------------------------
!BEGIN

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   temp = _STATE_VAR_(data%id_temp)   ! local temperature

   ! Retrieve current (local) state variable values.
   ph = _STATE_VAR_(data%id_ph)       ! local pH
   oxy = _STATE_VAR_(data%id_o_oxy)

   DO i=1,data%num_comp

     oxyEffect = one_
     pHEffect = one_

     ! Sediment flux dependent on oxygen
      IF( data%Ksed_gch_o2(i)<-1e-8 ) THEN
        oxyEffect = ( data%Ksed_gch_o2(i) / (data%Ksed_gch_o2(i) + oxy) )
      ELSEIF( data%Ksed_gch_o2(i)>1e-8 ) THEN
        oxyEffect = ( oxy / (data%Ksed_gch_o2(i) + oxy) )
      ENDIF
     ! Sediment flux dependent on pH
      IF( data%Ksed_gch_pH(i)<-1e-8 ) THEN
        pHEffect = ( data%Ksed_gch_pH(i) / (data%Ksed_gch_pH(i) + abs(ph-7.0)) )
      ELSEIF( data%Ksed_gch_pH(i)>1e-8 ) THEN
        pHEffect = ( abs(ph-7.0)  / (data%Ksed_gch_pH(i) + abs(ph-7.0)) )
      ENDIF

      gch_flux = data%Fsed_gch(i) * 1.05**(temp-20.0) * oxyEffect * pHEffect

      IF (.NOT.data%component_linked(i)) THEN
        ! geochem module variables, fluxed here as mmol/m2/day
         _FLUX_VAR_(data%id_comp(i)) =  _FLUX_VAR_(data%id_comp(i)) + gch_flux
      ELSE
        ! other module variables, can be fluxed by them
        !_FLUX_VAR_(data%id_cdep(i)) =  _FLUX_VAR_(data%id_comp(i)) + gch_flux
      ENDIF
   ENDDO

END SUBROUTINE aed_calculate_benthic_geochemistry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_equilibrate_geochemistry(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Update partitioning of phosphate between dissolved and particulate pools
! after kinetic transformations are applied
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_geochemistry_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: temp, sal
   ! State
   AED_REAL,   DIMENSION(SIZE(data%DissComp))  :: dissConcs
   AED_REAL,   DIMENSION(SIZE(data%PartComp))  :: partConcs
   ! Temporary variables
   INTEGER  :: i
   AED_REAL :: pco2,nc, MeDis,MePar
   AED_REAL,PARAMETER :: one_e_neg_ten = 1e-10

!-------------------------------------------------------------------------------
!BEGIN

   !-- Retrieve current environmental conditions for the cell.
   temp = _STATE_VAR_(data%id_temp) ! local temperature
   sal  = _STATE_VAR_(data%id_sal)  ! local salinity

   !-- Retrieve current gch state variable values into work array for the gcsolver
   DO i=1,data%num_comp
      IF (.NOT.data%component_linked(i)) THEN
          dissConcs(i) = _STATE_VAR_(data%id_comp(i))
      ELSE
          dissConcs(i) = _STATE_VAR_(data%id_cdep(i))
      ENDIF
   ENDDO
   DO i=1,data%num_mins
      IF (.NOT.data%mineral_linked(i)) THEN
          partConcs(i) = _STATE_VAR_(data%id_mins(i))
      ELSE
          partConcs(i) = _STATE_VAR_(data%id_mdep(i))
      ENDIF
   ENDDO

   !-- Do geochemical equilibration
   IF (data%simEq) &
      CALL UpdateEquilibration(dissConcs, partConcs, concMode=2, &
                               inTemp=REAL(temp), inSalt=REAL(sal), &
                               stoEq=.true., upDerv=.true.)


   !-- Copy back into main AED arrays
   DO i=1,data%num_comp
      IF (.NOT.data%component_linked(i)) THEN
         _STATE_VAR_(data%id_comp(i)) =  dissConcs(i)
      ELSE
         _STATE_VAR_(data%id_cdep(i)) =  dissConcs(i)
      ENDIF
   ENDDO
   DO i=1,data%num_mins
      IF (.NOT.data%mineral_linked(i)) THEN
         _STATE_VAR_(data%id_mins(i)) =  partConcs(i)
      ELSE
         _STATE_VAR_(data%id_mdep(i)) =  partConcs(i)
      ENDIF
   ENDDO

   !-- Update diagnostic arrays
   IF( returnGCDerivedVector("pCO2",pco2) > 0 .AND. data%simEq ) THEN
     !print *,'pco2: ',pco2
     _DIAG_VAR_(data%id_c_pco2) = pco2
     !_DIAG_VAR_(data%id_gcdiag(6)) = pco2
   ENDIF
   !-- Per-cell solver failure flag: 0 = accepted, 1 = rejected. Written only
   !-- to the dedicated noncon_mh diagnostic; the id_gcdiag() index for NONCON
   !-- is configuration dependent so it must not be hardwired here.
   IF( returnGCDerivedVector("NONCON",nc) > 0 .AND. data%simEq ) THEN
     _DIAG_VAR_(data%id_noncon) = nc
   ENDIF



   !Adsorption - instantaneous equilibrium partitioning. Components with
   !Rsorp_gch>0 are excluded here; they relax kinetically toward the same
   !partition in aed_calculate_geochemistry.
   IF( data%simMeAdsorption ) THEN
     DO i=1,data%num_comp
        IF ( data%MeAdsorptionModel(i) ==0 ) CYCLE
        IF ( data%component_linked(i) ) CYCLE   ! linked components have no _ads pool
        IF ( data%id_compd(i) <= 0 ) CYCLE      ! no sorbed state was registered
        IF ( data%Rsorp_gch(i) > zero_ ) CYCLE  ! kinetic: handled in calculate

        CALL calcSorptionPartition(data,column,layer_idx,i,MeDis,MePar)

        ! Set back to core variables
        _STATE_VAR_(data%id_comp(i)) = MeDis
        _STATE_VAR_(data%id_compd(i)) = MePar

        IF ( data%id_compd_frc(i) > 0 ) &
           _DIAG_VAR_(data%id_compd_frc(i)) = MePar / MAX(MeDis+MePar, one_e_neg_ten)
     ENDDO
   ENDIF

END SUBROUTINE aed_equilibrate_geochemistry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_mobility_geochemistry(data,column,layer_idx,mobility)
!-------------------------------------------------------------------------------
! Set vertical movement of the sorbed (_ads) pools (+ve up; -ve down).
! Components with link_ads_vvel adopt the sorbent particle's _vvel diagnostic;
! others keep their constant w_gch_ads. Mineral mobilities are define-time
! constants and are not touched here.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_geochemistry_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: mobility(:)
!
!LOCALS
   AED_REAL :: vvel, dz
   INTEGER  :: i
!-------------------------------------------------------------------------------
!BEGIN
   IF (.NOT. data%simMeAdsorption) RETURN

   vvel = zero_
   IF ( data%id_compd_vvel > 0 ) &
      vvel = _DIAG_VAR_(data%id_compd_vvel) / secs_per_day

   dz = _STATE_VAR_(data%id_dz)

   DO i=1,data%num_comp
      IF ( data%MeAdsorptionModel(i) == 0 ) CYCLE
      IF ( data%component_linked(i) ) CYCLE
      IF ( data%id_compd(i) <= 0 ) CYCLE

      IF ( data%link_ads_vvel(i) ) THEN
         mobility(data%id_compd(i)) = vvel
      ELSE
         mobility(data%id_compd(i)) = data%w_gch_ads(i)
      ENDIF

      IF ( data%id_compd_set(i) > 0 ) &
         _DIAG_VAR_(data%id_compd_set(i)) = (mobility(data%id_compd(i))/dz)   &
                          * _STATE_VAR_(data%id_compd(i)) * secs_per_day
   ENDDO

END SUBROUTINE aed_mobility_geochemistry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_inflow_update_geochemistry(data, wqinf, temp, salt)
!-------------------------------------------------------------------------------
   CLASS (aed_geochemistry_data_t),INTENT(in) :: data
   AED_REAL,DIMENSION(:),INTENT(inout) :: wqinf
   AED_REAL,             INTENT(inout) :: temp, salt
!-------------------------------------------------------------------------------
!
!LOCALS
   ! State
   AED_REAL,   DIMENSION(SIZE(data%DissComp))  :: dissConcs
   AED_REAL,   DIMENSION(SIZE(data%PartComp))  :: partConcs
   ! Temporary variables
   INTEGER  :: i
!-------------------------------------------------------------------------------
!BEGIN
   IF ( data%inflow_pH_update == 0 ) RETURN
!print *,'wqinf',wqinf
   !-- Reset inflow state variable values into array for the gcsolver
   DO i=1,data%num_comp
      IF (.NOT.data%component_linked(i)) THEN
          dissConcs(i) = wqinf(data%id_comp(i))
      ELSE
          dissConcs(i) = wqinf(data%id_cdep(i))
      ENDIF
  !    print *,'IN i',i,dissConcs(i),data%id_comp(i),data%id_comp(i)
   ENDDO
   DO i=1,data%num_mins
      IF (.NOT.data%mineral_linked(i)) THEN
          partConcs(i) = wqinf(data%id_mins(i))
      ELSE
          partConcs(i) = wqinf(data%id_mdep(i))
      ENDIF
   ENDDO

   !-- Redo geochemical equilibration, now spatial initialisation is done
   CALL InitialiseGCProperties(dissConcs, partConcs, 2, inTemp=REAL(temp))

   !-- Copy back into main AED arrays
   DO i=1,data%num_comp
      IF (.NOT.data%component_linked(i)) THEN
         wqinf(data%id_comp(i)) =  dissConcs(i)
      ELSE
         wqinf(data%id_cdep(i)) =  dissConcs(i)
      ENDIF
!      print *,'OUT i',i,dissConcs(i)

   ENDDO
   DO i=1,data%num_mins
      IF (.NOT.data%mineral_linked(i)) THEN
         wqinf(data%id_mins(i)) =  partConcs(i)
      ELSE
         wqinf(data%id_mdep(i)) =  partConcs(i)
      ENDIF
   ENDDO
END SUBROUTINE aed_inflow_update_geochemistry
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
 FUNCTION calcIronOxidation(data, ironII, oxygen, temp) RESULT(OxidRate)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_geochemistry_data_t),INTENT(in) :: data
   AED_REAL, INTENT(IN)  :: ironII
   AED_REAL, INTENT(IN)  :: oxygen
   AED_REAL, INTENT(IN)  :: temp
   AED_REAL      :: OxidRate
!LOCALS
   AED_REAL, PARAMETER :: ko = 2.14e-5  !@25C
   AED_REAL, PARAMETER :: k1 = 6.78e1   !@25C
   AED_REAL, PARAMETER :: k2 = 2.14e7   !@25C
   AED_REAL, PARAMETER :: FEII_MolWeight = 55.79
   AED_REAL  :: Fe_2, FeOH, FeOH2
   AED_REAL  :: kb, kf
   INTEGER   :: status
!-------------------------------------------------------------------------------
!BEGIN

   kb = data%Riron_box
   kf = data%Riron_aox

   IF ( data%iron_ox_model == 1 ) THEN
      !-- Singer-Stumm: abiotic Fe(II) oxidation is first order in O2 and is
      !-- summed over the hydrolysis series, since FeOH+ and Fe(OH)2 oxidise
      !-- orders of magnitude faster than Fe+2. The pH dependence is implicit,
      !-- carried entirely by the speciation - which is why the rate falls by
      !-- ~10^6 from pH 7 to pH 3 and keeps Fe(II) mobile in acid drainage.
      !-- The biotic term kb is retained: it is what dominates at low pH.
      !-- CAUTION: restored verbatim from the long commented out original and
      !-- NOT yet mass balance verified. The speciation terms carry moles from
      !-- allSpecies scaled by FEII_MolWeight*1e3, while the biotic term uses
      !-- ironII straight from the state variable in mmol/m3 - confirm the two
      !-- bases agree before trusting a calibration against Riron_aox.
      status = returnGCDerivedVector("Fe+2      ", Fe_2)
      IF(status/=1)Fe_2 = zero_
      status = returnGCDerivedVector("FeOH+     ", FeOH)
      IF(status/=1)FeOH = zero_
      status = returnGCDerivedVector("Fe(OH)2   ", FeOH2)
      IF(status/=1)FeOH2 = zero_

      OxidRate = oxygen *  (                                                   &
                 kf*ko*Fe_2  *FEII_MolWeight*1e3  +                            &
                 kf*k1*FeOH  *FEII_MolWeight*1e3  +                            &
                 kf*k2*FeOH2 *FEII_MolWeight*1e3  +                            &
                 kb*ironII ) * (data%theta_iron_ox**(temp-20.0))
   ELSE
      !-- Biotic only, Monod in O2. NB pH independent: Fe(II) oxidises at the
      !-- same rate at pH 3 as at pH 7, which is not realistic for acid
      !-- sulfate systems. Use iron_ox_model = 1 where that matters.
      OxidRate = kb*ironII * (data%theta_iron_ox**(temp-20.0))                 &
                           * oxygen/(data%Kiron_ox + oxygen)
   ENDIF

 END FUNCTION calcIronOxidation
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
 FUNCTION calcIronReduction(data, ironIII, oxygen, temp) RESULT(RednRate)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_geochemistry_data_t),INTENT(in) :: data
   AED_REAL, INTENT(IN)  :: ironIII
   AED_REAL, INTENT(IN)  :: oxygen
   AED_REAL, INTENT(IN)  :: temp
   !-- Outgoing
   AED_REAL :: RednRate
   !-- Local
   AED_REAL :: light
   AED_REAL :: FeOH2
   AED_REAL :: ko
   INTEGER  :: status
!-------------------------------------------------------------------------------
!BEGIN
   RednRate = zero_

   !-- Biotic rate: iron reducers
   RednRate  = (data%Riron_red  * (data%theta_iron_red**(temp-20.0)) *         &
                data%Kiron_red / (data%Kiron_red + oxygen)) * ironIII

   !-- Photo-reduction
!   IF() THEN
!     light = zero_
!     light = _STATE_VAR_(data%id)
!
!     !-- Photo-reduction rate (/day) is a linear function of PAR
!     ko = data%kFeRpr * light/2.5e3
!
!     FeOH2 = zero_
!~status = returnGCDerivedVector("Fe(OH)2   ", FeOH2)
!~IF(status/=1)FeOH2 = zero_
!     !#MH: Causing ELCD crash 20090817. Temporailiy Disabled
!     !RednRate = RednRate +                                                   &
!     !           ko(vdo) * FeOH2(vdo) * VarDetails(DICHM(FEIII))%MolWeight*1e3
!     RednRate = zero_
!   END IF

 END FUNCTION calcIronReduction
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
 FUNCTION calcSulfurReduction(data,sulfate,oxygen,temp) RESULT(RednRate)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_geochemistry_data_t),INTENT(in) :: data
   AED_REAL, INTENT(IN)  :: sulfate
   AED_REAL, INTENT(IN)  :: oxygen
   AED_REAL, INTENT(IN)  :: temp
   !-- Outgoing
   AED_REAL              :: RednRate
!-------------------------------------------------------------------------------
!BEGIN

   RednRate  = zero_
   RednRate  = data%Rsulf_red * (data%theta_sulf_red**(temp-20.0))             &
                * sulfate * (data%Ksulf_red / (data%Ksulf_red + oxygen))

 END FUNCTION calcSulfurReduction
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
 FUNCTION calcSulfurOxidation(data,sulfide,oxygen,temp) RESULT(OxdnRate)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_geochemistry_data_t),INTENT(in) :: data
   AED_REAL, INTENT(IN)  :: sulfide
   AED_REAL, INTENT(IN)  :: oxygen
   AED_REAL, INTENT(IN)  :: temp
   !-- Outgoing
   AED_REAL              :: OxdnRate
!-------------------------------------------------------------------------------
!BEGIN

   OxdnRate  = data%Rsulf_ox * (data%theta_sulf_ox**(temp-20.0))               &
                  * sulfide * (oxygen / (data%Ksulf_ox + oxygen))

 END FUNCTION calcSulfurOxidation
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE calcSorptionPartition(data,column,layer_idx,i,MeDis,MePar)
!-------------------------------------------------------------------------------
! Compute the equilibrium partition (dissolved vs sorbed) of component i for
! the current cell conditions. This is the target state shared by both
! sorption pathways: applied directly by aed_equilibrate_geochemistry
! (instantaneous) or relaxed toward by aed_calculate_geochemistry (kinetic).
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_geochemistry_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx, i
   AED_REAL,INTENT(out) :: MeDis, MePar
!
!LOCALS
   AED_REAL :: temp, sal, tss, pH, KMep, inDis, inPar
!-------------------------------------------------------------------------------
!BEGIN
   temp = _STATE_VAR_(data%id_temp)
   sal  = _STATE_VAR_(data%id_sal)

   tss = zero_
   IF (data%id_tss>0) tss = _STATE_VAR_(data%id_tss)  ! sorbent (TSS or particle)

   inDis = _STATE_VAR_(data%id_comp(i))
   inPar = _STATE_VAR_(data%id_compd(i))

   ! Adjust local sorption coefficient for temperature and salinity
   ! (only used by MeAdsorptionModel = 1, the linear Kd model)
   KMep = data%KMep(i) * KMe_fT_fSal(data%theta_KMe, data%K_sal, sal, temp)

   IF (data%ads_use_pH) THEN
      pH = _STATE_VAR_(data%id_pH)
      CALL MetalAdsorptionFraction(data%MeAdsorptionModel(i),                  &
                                   inDis+inPar, tss,                           &
                                   KMep, data%Kadsratio(i), data%Qmax(i),      &
                                   MeDis, MePar,                               &
                                   thepH=pH)
   ELSE
      CALL MetalAdsorptionFraction(data%MeAdsorptionModel(i),                  &
                                   inDis+inPar, tss,                           &
                                   KMep, data%Kadsratio(i), data%Qmax(i),      &
                                   MeDis, MePar,                               &
                                   temp_=temp, salt_=sal)
   ENDIF

END SUBROUTINE calcSorptionPartition
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
 SUBROUTINE MetalAdsorptionFraction(MetalAdsorptionModel, &
                                    MeTot_,               &
                                    ParticleConc_,        &
                                    KMep,K,Qm,            &
                                    MeDis,MePar,          &
                                    thepH,                &
                                    temp_,salt_)
!-------------------------------------------------------------------------------
! Routine to compute fraction of Me adsorped to sediment/particulate concentration
!-------------------------------------------------------------------------------
INTEGER,  INTENT(IN)  :: MetalAdsorptionModel
AED_REAL, INTENT(IN)  :: MeTot_, ParticleConc_
AED_REAL, INTENT(IN)  :: KMep,K,Qm
AED_REAL, INTENT(OUT) :: MeDis, MePar
AED_REAL, INTENT(IN), OPTIONAL :: thepH
AED_REAL, INTENT(in), OPTIONAL :: temp_
AED_REAL, INTENT(in), OPTIONAL :: salt_

!-------------------------------------------------------------------------------
!LOCALS
AED_REAL :: MeTot, ParticleConc
AED_REAL :: buffer, f_pH, pH                                        ! for option 2
AED_REAL :: parA, parB, Pexch, EPC0, Ppar, temp, salt, KMep_fT_fSal ! for option 3
AED_REAL,PARAMETER :: one_e_neg_ten = 1e-10
!
!-------------------------------------------------------------------------------
!BEGIN
MeDis   = zero_
MePar   = zero_
buffer   = zero_
f_pH     = one_

! calculate the total possible Me for sorption, and solids
MeTot        = MAX(one_e_neg_ten, MeTot_ )        ! Co in Chao (mg)
ParticleConc  = MAX(one_e_neg_ten, ParticleConc_ ) ! s in Chao  (mg = mol/L * g/mol * mg/g)


IF(MetalAdsorptionModel == 1) THEN
!-----------------------------------------------------
! This is the model for PO4 sorption from Ji 2008:
!
! Ji, Z-G. 2008. Hydrodynamics and Water Quality. Wiley Press.

MePar = (KMep*ParticleConc) / (one_+KMep*ParticleConc) * MeTot
MeDis = one_ / (one_+KMep*ParticleConc) * MeTot


ELSEIF(MetalAdsorptionModel == 2) THEN
!-----------------------------------------------------
! This is the model for PO4 sorption from Chao et al. 2010:
!
! Chao, X. et al. 2010. Three-dimensional numerical simulation of
!   water quality and sediment associated processes with application
!   to a Mississippi delta lake. J. Environ. Manage. 91 p1456-1466.
IF(PRESENT(thepH)) THEN
pH = MIN(MAX(2.0,thepH),12.0)

! -0.0094x2 + 0.0428x + 0.9574
! (ursula.salmon@uwa.edu.au: fPH for PO4 sorption to Fe in Mine Lakes)
f_pH = -0.0094*pH*pH + 0.0428*pH + 0.9574
ELSE
f_pH = one_
END IF

! calculate particulate fraction based on quadratic solution

! Chao Eq 16
buffer = SQRT(((MeTot+(1./K)-(ParticleConc*Qm*f_pH)))**2. + (4.*f_pH*ParticleConc*Qm/K))
MePar  = 0.5 * ((MeTot+(1./K)+(ParticleConc*Qm*f_pH))  - buffer  )

! Check for stupid solutions
IF(MePar > MeTot) MePar = MeTot
IF(MePar < zero_) MePar = zero_

! Now set dissolved portion
MeDis = MeTot - MePar


ELSE
!-----------------------------------------------------
! No model is selected

MeDis = MeTot
MePar = zero_

END IF

END SUBROUTINE MetalAdsorptionFraction
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION KMe_fT_fSal(theta_KMe, K_sal, sal, temp)
!-------------------------------------------------------------------------------
! Temperature / salinity sensitivity of the linear sorption coefficient,
! following the P-adsorption form of Zhang and Huang 2011 (as in
! aed_phosphorus::Kpo4p_fT_fSal): sorption increases with temperature and
! decreases with salinity. Defaults (theta_KMe=1, K_sal=0) are neutral.
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in) :: theta_KMe ! theta for T dependence; 1.0 => no effect
   AED_REAL,INTENT(in) :: K_sal     ! half-saturation of salinity; 0 => no effect
   AED_REAL,INTENT(in) :: sal
   AED_REAL,INTENT(in) :: temp
   AED_REAL,PARAMETER  :: Topt = 45.
!
!LOCALS
   AED_REAL :: fT, fSal
!
!-------------------------------------------------------------------------------
!BEGIN
   fT = (theta_KMe**(temp-Topt))

   IF(K_sal==zero_)THEN
     fSal = one_
   ELSE
     fSal = K_sal/(K_sal+sal)
   ENDIF

   KMe_fT_fSal = fT * fSal

END FUNCTION KMe_fT_fSal
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE reportRedoxCouple(label, active)
!-------------------------------------------------------------------------------
! Print whether one kinetic redox couple was activated by the component set
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(len=*),INTENT(in) :: label
   LOGICAL,INTENT(in)          :: active
!
!-------------------------------------------------------------------------------
!BEGIN
   !-- label is space padded by the caller for alignment, so do not TRIM it
   IF ( active ) THEN
      print *,"          ",label," : ACTIVE"
   ELSE
      print *,"          ",label," : off (both redox states not configured)"
   ENDIF

END SUBROUTINE reportRedoxCouple
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_geochemistry
