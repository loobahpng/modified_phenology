!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
module mo_cbal_cpools
    USE mo_kind,             ONLY: dp
    USE mo_cbal_parameters,  ONLY: cn_green, cn_woods, cn_litter_green, cn_litter_wood, cn_slow
    USE mo_land_surface,     ONLY: fract_small
    USE mo_jsbach_constants, ONLY: Tmelt, molarMassC_kg, molarMassN_kg

    implicit none

    public :: update_Cpools
    public :: N_process
    public :: relocate_CarbonAndNitrogen
    public :: C_relocation_from_LUtransitions
    public :: C_relocation_from_harvest
    PUBLIC :: C_relocation_from_forest_management_clear_cuts, C_relocation_from_product_pool_decay
    public :: relocate_carbon_desert
    public :: relocate_carbon_fire
    public :: relocate_carbon_damage
    public :: printCpoolParameters
    public :: printNpoolParameters
    private :: yasso
    private :: N_loss_lcc

    REAL(dp), parameter,public  :: frac_wood_aboveGround  = 0.7_dp  ! Fraction of Carbon above ground in wood pool (for separation
                                                                    !    of woody litter into above and below ground litter pools)
    REAL(dp), parameter,public  :: frac_green_aboveGround = 0.5_dp  ! Fraction of Carbon above ground in green pool (for separation
                                                                    !    of green litter into above and below ground litter pools)
   

    ! === END OF PUBLIC PART === BEGIN OF PRIVATE PART =============================================================================

    private

    REAL(dp), parameter :: Q10                    = 1.8_dp       ! Empirical parameter for temperature dependence of heterotrophic
                                                                 !    respiration rate
    REAL(dp), parameter :: kappa                  = 1.0_dp       ! Empirical parameter for soil moisture dependence of heterotrophic
                                                                 !    respiration rate
    REAL(dp), parameter :: referenceTemp_Q10      = Tmelt            ! Reference temperature in the Q10 formula [Kelvin]
    REAL(dp), parameter :: tau_Cpool_reserve      = 1.0_dp *365._dp  ! time constant by which reserve pool is depreciated [days]

    REAL(dp), parameter :: tau_Cpool_slow         = 100._dp*365._dp  ! time constant by which slow pool is depreciated  [days]
    REAL(dp), parameter :: tau_Cpool_crop_harvest = 1._dp*365._dp    ! time constant by which the crop harvest pool decays [days]
    REAL(dp), parameter :: greenC2leafC           = 4.0_dp           ! ratio of carbon green pool (leaves, fine roots, starches,
                                                                     !    sugars) to leaf carbon
    REAL(dp), parameter :: frac_C_litter_wood2atmos = 0.2_dp     ! fraction of heterotrophic loss of woody litter emitted to the
                                                                 !    atmosphere
    REAL(dp), parameter :: frac_C_faeces2_LG      = 0.3_dp       ! fraction of C from faces put into the litter green pool
    REAL(dp), parameter :: frac_C_crop_harvest    = 0.5_dp       ! fraction of C from ag litter flux put to the crop harvest pool 
    REAL(dp), parameter :: alpha_critical         = 0.35_dp      ! critical value of relative soil moisture below which
                                                                 !    heterotrophic respiration almost stops
    REAL(dp), parameter :: alpha_min              = 0.1_dp       ! effective lowest value for soil moisture entering the
                                                                 !    heterotrophic respiration
    REAL(dp), parameter :: N2O_rate_nitrification = 0.001_dp     ! JSBACH assumes lower bound of estimation in the literture as: 
                                                                 !   (0.1%) due to nitrification 
    REAL(dp), parameter :: N2O_rate_denitrification = 0.00125_dp ! 0.0025_dp  ! & (0.2%) due to denitrification 
                                                                 ! actually 0.2% of denitrified N, not of N inputs
                                                                 ! assumption: 50% loss of NO3 by leaching
    REAL(dp), parameter :: N2O_rate               = 0.005_dp     ! JSBACH assumes 0.5-1%
                                                                 !    N2O emissions rate per day (<0.1%-0.2% nitrification ); 
                                                                 !    0.2-4.7% due to denit in Xu Ri, 2008 * therein
    REAL(dp), parameter :: sminn_NH4_fraction     = 0.4_dp       ! soil mineral NH4 is 40% of SMINN pool and rest in the form of NO3
                                                                 ! Ref: Xu Ri, 2008

    REAL(dp), parameter :: sec_per_day            = 86400._dp    ! seconds per day
    REAL(dp), parameter :: days_per_year          = 365.25_dp
    REAL(dp), parameter :: sec_per_year           = days_per_year*sec_per_day

contains

  ! --- printCpoolParameters() -----------------------------------------------------------------------------------------------
  !
  ! Prints out the parameters used for the carbon pool model
  !
  SUBROUTINE printCpoolParameters
    USE mo_exception,  ONLY: message, message_text

    WRITE (message_text,*) "=== CpoolParameters ========================================="
    CALL message("",message_text)
    WRITE (message_text,*) "                     Q10=", Q10
    CALL message("",message_text)
    WRITE (message_text,*) "                   kappa=", kappa
    CALL message("",message_text)
    WRITE (message_text,*) "       referenceTemp_Q10=", referenceTemp_Q10, " Kelvin"
    CALL message("",message_text)
    WRITE (message_text,*) "       tau_Cpool_reserve=", tau_Cpool_reserve/365._dp, " years"
    CALL message("",message_text)
    WRITE (message_text,*) "          tau_Cpool_slow=", tau_Cpool_slow/365._dp, " years"
    CALL message("",message_text)
    WRITE (message_text,*) "  tau_Cpool_crop_harvest=", tau_Cpool_crop_harvest/365._dp, " years"
    CALL message("",message_text)
    WRITE (message_text,*) "            greenC2leafC=", greenC2leafC
    CALL message("",message_text)
    WRITE (message_text,*) "frac_C_litter_wood2atmos=", frac_C_litter_wood2atmos
    CALL message("",message_text)
    WRITE (message_text,*) "  frac_green_aboveGround=", frac_green_aboveGround
    CALL message("",message_text)
    WRITE (message_text,*) "   frac_wood_aboveGround=", frac_wood_aboveGround
    CALL message("",message_text)
    WRITE (message_text,*) "       frac_C_faeces2_LG=", frac_C_faeces2_LG
    CALL message("",message_text)
    WRITE (message_text,*) "     frac_C_crop_harvest=", frac_C_crop_harvest
    CALL message("",message_text)
    WRITE (message_text,*) "          alpha_critical=", alpha_critical
    CALL message("",message_text)
    WRITE (message_text,*) "               alpha_min=", alpha_min
    CALL message("",message_text)
    WRITE (message_text,*) "============================================================="
    CALL message("",message_text)

  END SUBROUTINE printCpoolParameters

 ! --- printNpoolParameters() ------------------------------------------------------------------------------------------------------
  !
  ! Prints out the parameters used for the nitrogen pool model
  !
  SUBROUTINE printNpoolParameters
    USE mo_exception,       ONLY: message, message_text

    WRITE (message_text,*) "=== NpoolParameters ========================================="
    CALL message("",message_text)
    WRITE (message_text,*) "        cn_green=", cn_green
    CALL message("",message_text)
    WRITE (message_text,*) "        cn_woods=", cn_woods
    CALL message("",message_text)
    WRITE (message_text,*) " cn_litter_green=", cn_litter_green
    CALL message("",message_text)
    WRITE (message_text,*) "  cn_litter_wood=", cn_litter_wood
    CALL message("",message_text)
    WRITE (message_text,*) "         cn_slow=", cn_slow
    CALL message("",message_text)
    WRITE (message_text,*) "============================================================="

  END SUBROUTINE printNpoolParameters 


  ! --- update_Cpools() -----------------------------------------------------------------------------------------------------------
  !
  ! Updates carbon and nitrogen pools and computes soil respiration rate and Net Ecosystem Product [mol(C)/m^2 s].
  ! Has to be called once each day for each covertype, where carbon and nitrogen pools exists. 
  !

  elemental subroutine update_Cpools(LAI, LAI_previous, NPPrate, topSoilTemp, alpha,                          &
                                          frac_npp_2_woodPool,frac_npp_2_reservePool, frac_npp_2_exudates,    &
                                          frac_green_2_herbivory,                                             &
                                          tau_Cpool_litter_green,tau_Cpool_litter_wood,tau_Cpool_woods,       &
                                          LAI_shed_constant,                                                  &
                                          frac_C_litter_green2atmos,Max_C_content_woods,                      &
                                          specific_leaf_area_C, reserveC2leafC,                               &
                                          max_LAI, is_vegetation, is_crop,                                    &
                                          Cpool_green, Cpool_woods, Cpool_reserve,                            &
                                          Cpool_crop_harvest,                                                 &
                                          Cpool_litter_green_ag,Cpool_litter_green_bg,                        &
                                          Cpool_litter_wood_ag,Cpool_litter_wood_bg,                          & 
                                          Cpool_slow,                                                         &
                                          soilResp_rate, soilResp_rate_pot,                                   &
                                          NPP_flux_correction, excess_NPP,root_exudates,                      &
                                          Cflx_litterTotal,                                                   &
                                          Cflx_herbivory,Cflx_herbivory_LG, Cflx_herbivory_2_atm,             &
                                          Cflx_2_crop_harvest, Cflx_crop_harvest_2_atm,                       &
                                          NPP_act,&
                                          leaf_shedding_debug,   & ! HW   
                                          Cpool_green_minus_shed_leaves,     & ! HW
                                          excess_carbon_debug, & !HW  
                                          litter_leaf,   & ! HW   
                                          leaf_shedding_rate, & ! HW  
                                          frac_litter_wood_new,                                      &
                                          ! Yasso variables
                                          temp2_30d, precip_30d,                                              &
                                          YCpool_acid_ag1, YCpool_water_ag1, YCpool_ethanol_ag1,              &
                                          YCpool_nonsoluble_ag1,                                              & 
                                          YCpool_acid_bg1, YCpool_water_bg1, YCpool_ethanol_bg1,              &
                                          YCpool_nonsoluble_bg1, YCpool_humus_1,                              & 
                                          YCpool_acid_ag2, YCpool_water_ag2, YCpool_ethanol_ag2,              &
                                          YCpool_nonsoluble_ag2,                                              & 
                                          YCpool_acid_bg2, YCpool_water_bg2, YCpool_ethanol_bg2,              &
                                          YCpool_nonsoluble_bg2, YCpool_humus_2,                              & 
                                          LeafLit_coef1, LeafLit_coef2, LeafLit_coef3, LeafLit_coef4,         &
                                          LeafLit_coef5, WoodLit_coef1, WoodLit_coef2,                        &
                                          WoodLit_coef3, WoodLit_coef4, WoodLit_coef5,                        &
                                          WoodLitterSize,                                                     &
                                          ! Nitrogen variables
                                          redFact_Nlimit,                                                     &
                                          Npool_green, Npool_woods, Npool_mobile,                             &
                                          Npool_crop_harvest,                                                 &
                                          Npool_litter_green_ag,Npool_litter_green_bg,                        &
                                          Npool_litter_wood_ag,Npool_litter_wood_bg,                          &
                                          Npool_slow, SMINN_pool,                                             &
                                          minNflux_litter_green_ag,minNflux_litter_green_bg,                  &
                                          minNflux_litter_wood_ag,minNflux_litter_wood_bg,                    &
                                          minNflux_slow,Nplant_demand,Nsoil_demand,                           &
                                          Ntotal_demand,SMINN_herbivory,N2O_emissions_mineraliz,              &      
                                          N2O_emissions_slow, N2O_emissions_grazing, Nflx_2_crop_harvest      &
                                           ) 

    real(dp),intent(in)    :: LAI                    !! Yesterdays mean LAI
    real(dp),intent(in)    :: LAI_previous           !! The day before yesterdays mean LAI
    real(dp),intent(in)    :: NPPrate                !! Yesterdays mean NPPrate [mol(C)/m^2 s]
    real(dp),intent(in)    :: topSoilTemp            !! Yesterdays mean temperature of upper soil layer [degree Kelvin]
    real(dp),intent(in)    :: alpha                  !! Yesterdays mean water stress factor (between 0 and 1)
    real(dp),intent(in)    :: frac_npp_2_woodPool    !! Fraction of NPP to be put maximally into the green pool
    real(dp),intent(in)    :: frac_npp_2_reservePool !! Fraction of NPP to be put into the optimally into the reserve pool
    real(dp),intent(in)    :: frac_npp_2_exudates    !! Fraction of NPP to be put into the optimally into the root exudates
    real(dp),intent(in)    :: frac_green_2_herbivory
    real(dp),intent(in)    :: tau_Cpool_litter_green !! Time constant by which green litter pool is depreciated  [days]
    real(dp),intent(in)    :: tau_Cpool_litter_wood  !! Time constant by which woody litter pool is depreciated  [days]
    real(dp),intent(in)    :: tau_Cpool_woods        !! Time constant by which woods Pool is depreciated  [days]
    real(dp),intent(in)    :: LAI_shed_constant      !! Leaf shedding at a constant rate for evergreens etc.
    real(dp),intent(in)    :: frac_C_litter_green2atmos !! Fraction of heterotrophic loss of green litter pools emitted to the 
                                                        !!    atmosphere
    real(dp),intent(in)    :: Max_C_content_woods    !! Maximum carbon content of wood pool (from lctLib)
    real(dp),intent(in)    :: specific_leaf_area_C   !! Specific leaf area (from lctLib) [m^2(leaf)/mol(Carbon)]
    real(dp),intent(in)    :: reserveC2leafC         !! Ratio of max. carbon in reserve pool to max. carbon of leaves
    real(dp),intent(in)    :: max_LAI                !! Maximum value of LAI
    logical, intent(in)    :: is_vegetation          !! logical vegetation mask
    logical, intent(in)    :: is_crop                !! logical crop mask
    real(dp),intent(inout) :: Cpool_green            !! Green carbon pool: on input last value; updated on output 
                                                     !!    [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_woods            !! Wood carbon pool: on input last value; updated on output
                                                     !!    [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_reserve          !! Reserve carbon pool: on input last value; updated on output
                                                     !!    [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_crop_harvest     !! Crop harvest carbon pool: on input last value; updated on output
                                                     !!    [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_litter_green_ag  !! Above ground green litter C-pool: updated on output [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_litter_green_bg  !! Below ground green litter C-pool: updated on output [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_litter_wood_ag   !! Above ground woody litter C-pool: updated on output [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_litter_wood_bg   !! Below ground woody litter C-pool: updated on output [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_slow             !! Slow soil carbon pool: on input last value; updated on output
                                                     !!    [mol(C)/m^2(canopy)]
    real(dp),intent(out)   :: soilResp_rate          !! Soil (=heterotrophic) respiration rate  [mol(C)/m^2 s] Note: this is a loss,
                                                     !!    therefore negative!
    real(dp),intent(out)   :: soilResp_rate_pot      !! soil respiration without N limitation     
    real(dp),intent(out)   :: NPP_flux_correction    !! Amount by which the NPP rate entering the routine has to be corrected. This
                                                     !!    correction arises either because otherwise the reserve pool would get
                                                     !!    negative (positive correction), or the wood pool would exceed its
                                                     !!    maximum value (negative correction). [mol(C)/m^2 s]
    real(dp),intent(out)   :: excess_NPP             !! Part of NPP that could not be stored in one of the plant carbon pools 
                                                     !!    (green, wood, reserve), but had to be  thrown away (a posteriori 
                                                     !!    reduction of NPP) [mol(C)/m^2 s]
    real(dp),intent(out)   :: root_exudates          !! Value of NPP_2_rootExudates had to be dropped into the litter_green_bg pool.
                                                     !!    [mol(C)/m^2 s]
    real(dp),intent(out)   :: Cflx_litterTotal       !! Total carbon flux from vegetation to litter (green+woody, ag+bg)

    real(dp),intent(out)   :: leaf_shedding_debug       !! HW
    real(dp),intent(out)   :: Cpool_green_minus_shed_leaves     !! HW
    real(dp),intent(out)   :: excess_carbon_debug       !! HW
    real(dp),intent(out)   :: litter_leaf       !! HW
    real(dp),intent(in)   :: leaf_shedding_rate       !! HW

    real(dp),intent(out)   :: Cflx_herbivory         !!
    real(dp),intent(out)   :: Cflx_herbivory_LG      !!
    real(dp),intent(out)   :: Cflx_herbivory_2_atm   !!
    real(dp),intent(out)   :: Cflx_2_crop_harvest    !! C flux from crops to the crop harvest pool [mol(C)m-2(canopy)s-1]
    real(dp),intent(out)   :: Cflx_crop_harvest_2_atm!! C flux from the crop harvest pool to the atmosphere [mol(C)m-2(canopy)s-1]
    real(dp),intent(out)   :: NPP_act                !! Actual NPP after N-limitation and excess carbon drop [mol(C)/m^2 s]

    ! added for thonicke fire algorithm (spitfire) needed to correctly split the carbon pools into the different fuel classes
    real(dp),optional,intent(inout) :: frac_litter_wood_new   !! new fraction in above ground wood litter pool  

    ! Meteorology for Yasso
    real(dp),optional,intent(in)    :: temp2_30d              !! 30 day mean temperature
    real(dp),optional,intent(in)    :: precip_30d             !! 30 day mean precipitation

    ! Yasso pools
    ! Size class 1; green litter
    !   Aboveground  
    real(dp),optional,intent(inout) :: YCpool_acid_ag1        !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_water_ag1       !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_ethanol_ag1     !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_nonsoluble_ag1  !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   Belowground  
    real(dp),optional,intent(inout) :: YCpool_acid_bg1        !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_water_bg1       !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_ethanol_bg1     !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_nonsoluble_bg1  !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_humus_1         !! Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    ! Size class 2; woody litter
    !   Aboveground  
    real(dp),optional,intent(inout) :: YCpool_acid_ag2        !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_water_ag2       !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_ethanol_ag2     !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_nonsoluble_ag2  !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   Belowground  
    real(dp),optional,intent(inout) :: YCpool_acid_bg2        !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_water_bg2       !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_ethanol_bg2     !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_nonsoluble_bg2  !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_humus_2         !! Humus litter pool for Yasso [mol(C)/m^2(canopy)]

    ! Parameters for Yasso
       ! Vegetation dependent coefficients to separate leaf litter into classes of chemical composition
    real(dp),optional,intent(in)    :: LeafLit_coef1         !! fraction going to the acid soluble pools 
    real(dp),optional,intent(in)    :: LeafLit_coef2         !! fraction going to the water soluble pools
    real(dp),optional,intent(in)    :: LeafLit_coef3         !! fraction going to the ethanol soluble pools
    real(dp),optional,intent(in)    :: LeafLit_coef4         !! fraction going to the non soluble pools
    real(dp),optional,intent(in)    :: LeafLit_coef5         !! fraction going to the humus pool
       ! Vegetation dependent coefficients to separate wood litter into classes of chemical composition
    real(dp),optional,intent(in)    :: WoodLit_coef1         !! fraction going to the acid soluble pools 
    real(dp),optional,intent(in)    :: WoodLit_coef2         !! fraction going to the water soluble pools
    real(dp),optional,intent(in)    :: WoodLit_coef3         !! fraction going to the ethanol soluble pools
    real(dp),optional,intent(in)    :: WoodLit_coef4         !! fraction going to the non soluble pools
    real(dp),optional,intent(in)    :: WoodLit_coef5         !! fraction going to the humus pool
    real(dp),optional,intent(in)    :: WoodLitterSize        !! size of coarse debris


    !! NOTE: all optional arguments are needed with runs with nitrogen
    real(dp),optional,intent(out)   :: redFact_Nlimit         !! Reduction factor for NPP due to N-limitation [1]
    real(dp),optional,intent(inout) :: Npool_green            !! Green N-pool: updated on output [mol(N)/m^2(canopy)]
    real(dp),optional,intent(inout) :: Npool_woods            !! Woods N-pool: updated on output [mol(N)/m^2(canopy)]
    real(dp),optional,intent(inout) :: Npool_mobile           !! Mobile plant N-pool: updated on output [mol(N)/m^2(canopy)]
    real(dp),optional,intent(inout) :: Npool_crop_harvest     !! Crop harvest N-pool: updated on output [mol(N)/m^2(canopy)]
    real(dp),optional,intent(inout) :: Npool_litter_green_ag  !! Above ground green litter N-pool: updated on output
                                                              !!    [mol(N)/m^2(canopy)]
    real(dp),optional,intent(inout) :: Npool_litter_green_bg  !! Below ground green litter N-pool: updated on output
                                                              !!    [mol(N)/m^2(canopy)]
    real(dp),optional,intent(inout) :: Npool_litter_wood_ag   !! Above ground woody litter N-pool: updated on output
                                                              !!    [mol(N)/m^2(canopy)]
    real(dp),optional,intent(inout) :: Npool_litter_wood_bg   !! Below ground woody litter N-pool: updated on output
                                                              !!    [mol(N)/m^2(canopy)]
    real(dp),optional,intent(inout) :: Npool_slow             !! Slow soil organic N-pool: updated on output [mol(N)/m^2(canopy)]
    real(dp),optional,intent(inout) :: SMINN_pool             !! Soil mineral nitrogen pool: updated on output [mol(N)/m^2(canopy)]
    real(dp),optional,intent(out)   :: minNflux_litter_green_ag !! Uptake of mineral N from above ground green litter pool
                                                                !!    [mol(N)/m^2 s]
    real(dp),optional,intent(out)   :: minNflux_litter_green_bg !! Uptake of mineral N from below ground green litter pool
                                                                !!    [mol(N)/m^2 s]
    real(dp),optional,intent(out)   :: minNflux_litter_wood_ag  !! Uptake of mineral N from above ground wood litter pool
                                                                !!    [mol(N)/m^2 s]
    real(dp),optional,intent(out)   :: minNflux_litter_wood_bg  !! Uptake of mineral N from below ground wood litter pool
                                                                !!    [mol(N)/m^2 s]
    real(dp),optional,intent(out)   :: minNflux_slow          !! Uptake of mineral N from slow soil pool  [mol(N)/m^2 s]
    real(dp),optional,intent(out)   :: Nplant_demand          !! Uptake of N for plant from soil [mol(N)/m^2 s] 
    real(dp),optional,intent(out)   :: Nsoil_demand           !! Uptake of N for soil microbes [mol(N)/m^2 s] 
    real(dp),optional,intent(out)   :: Ntotal_demand          !! Total N demand for both plant and soil [mol(N)/m^2 s]
    real(dp),optional,intent(out)   :: SMINN_herbivory        !! SMINN gained from herbivory  [mol(N)/m^2 s]
    real(dp),optional,intent(out)   :: N2O_emissions_mineraliz  !! Nitros oxide emission [mol(N)/m^2(canopy)s] due to mineralization
    real(dp),optional,intent(out)   :: N2O_emissions_slow     !! N2O emission [mol(N)/m^2(canopy)s] due to mineralization from
                                                              !!      slow pool
    real(dp),optional,intent(out)   :: N2O_emissions_grazing  !! N2O emission [mol(N)/m^2(canopy)s] due to N input by herbivores
    real(dp),optional,intent(out)   :: Nflx_2_crop_harvest    !! N flux from crops to the crop harvest pool [mol(N)m-2(canopy)s-1]


    ! locals

    real(dp) :: d_Cpool_slow
    real(dp) :: d_Cpool_litter_green_ag, d_Cpool_litter_green_bg
    real(dp) :: d_Cpool_litter_wood_ag, d_Cpool_litter_wood_bg
    real(dp) :: Cpool_green_max,Cpool_reserve_optimal
    real(dp) :: NPP_2_greenPool,NPP_2_woodPool,NPP_2_reservePool,NPP_2_rootExudates,Green_2_herbivory
    real(dp) :: C_2_litter_greenPool,C_2_litter_woodPool
    real(dp) :: excess_carbon
    real(dp) :: leaf_shedding
    real(dp) :: alpha_mod

    real(dp) :: Cpool_green_pot, Cpool_woods_pot         !! values of potential carbon pools i.e. before accounting for N-limitation
    real(dp) :: Cpool_slow_pot, Cpool_reserve_pot
    real(dp) :: Cpool_litter_green_ag_pot, Cpool_litter_green_bg_pot
    real(dp) :: Cpool_litter_wood_ag_pot, Cpool_litter_wood_bg_pot

    real(dp) :: Cflx_NPP_2_green_pot,Cflx_NPP_2_wood_pot,Cflx_NPP_2_reserve_pot 
    real(dp) :: Cflx_litterGreen_ag_2_atmos_pot, Cflx_litterGreen_bg_2_atmos_pot
    real(dp) :: Cflx_litterGreen_ag_2_slow_pot, Cflx_litterGreen_bg_2_slow_pot
    real(dp) :: Cflx_litter_wood_ag_2_atmos_pot,Cflx_litter_wood_bg_2_atmos_pot
    real(dp) :: Cflx_litter_wood_ag_2_slow_pot,Cflx_litter_wood_bg_2_slow_pot
    real(dp) :: Cflx_slow_2_atmos
    real(dp) :: Cflx_faeces_2_LG, Cflx_faeces_2_atm
    real(dp) :: litter_wood_new, litter_wood_total

    ! local variables for Yasso
    real(dp), dimension(18)  :: Yasso_io_pools
    real(dp), dimension(2)  :: Weather
    real(dp), dimension(16) :: Yasso_out 
    real(dp), dimension(5)  :: LeafLit_coefV
    real(dp), dimension(5)  :: WoodLit_coefV   
    real(dp)                :: leafLitter      ! amount of leafLitter 
    real(dp)                :: sizeZero

    real(dp) :: Nflux_NPP_2_green_pot,Nflux_NPP_2_wood_pot 
    real(dp) :: Nflux_litterGreen_ag_2_slow_pot, Nflux_litterGreen_bg_2_slow_pot
    real(dp) :: Nflux_litter_wood_ag_2_slow_pot,Nflux_litter_wood_bg_2_slow_pot
    real(dp) :: minNflux_litter_green_ag_pot,minNflux_litter_green_bg_pot
    real(dp) :: minNflux_litter_wood_ag_pot,minNflux_litter_wood_bg_pot
    real(dp) :: N2O_emissions_litter_pot 

    real(dp) :: minN_plant_demand,minN_soil_demand,total_minN_demand
    real(dp) :: decompRate_litter_green

    real(dp) :: NPP_2_green_via_mobileN, NPP_2_wood_via_mobileN
    real(dp) :: Nmobile_2_green, Nmobile_2_wood, Nmobile_2_green_pot, Nmobile_2_wood_pot
    real(dp) :: redFact_Nmobile, redFact_Nlimitation
    real(dp) :: Nflx_faeces_2_SMINN, minNflux_litter_active
    real(dp) :: Nflx_crop_harvest_2_SMINN

    LOGICAL  :: with_Nitrogen        !! .FALSE.: only Carbon pools are updated; .TRUE.: in addition also Nitrogen pools are updated 
    LOGICAL  :: with_yasso           !! .FALSE.: cbalance model for litter and soil carbon is used; .TRUE.: yasso model is used 
    REAL(dp), PARAMETER :: N2O_ef_grazing = 0._dp  !! N2O grazing switched off; reasonable value for the emission factor: 0.0125_dp 


    !! Determine whether call is with/without nitrogen
    IF (PRESENT(Npool_green)) THEN
       with_Nitrogen = .true.
    ELSE
       with_Nitrogen = .false.
    END IF

    !! Determine whether call is with/without yasso model
    IF (PRESENT(YCpool_acid_ag1)) THEN
       with_yasso = .true.
    ELSE
       with_yasso = .false.
    END IF

    ! Note: All fields with intent(out) are non-defined on return, even if the calling routine has set them.

    soilResp_rate            = 0.0_dp
    soilResp_rate_pot        = 0.0_dp
    NPP_act                  = 0.0_dp
    redFact_Nlimitation      = 1.0_dp

    IF (with_Nitrogen) THEN
       redFact_Nlimit           = 1.0_dp
       minNflux_litter_green_ag = 0.0_dp
       minNflux_litter_green_bg = 0.0_dp
       minNflux_litter_wood_ag  = 0.0_dp
       minNflux_litter_wood_bg  = 0.0_dp
       minNflux_slow            = 0.0_dp
       Nplant_demand            = 0.0_dp
       Nsoil_demand             = 0.0_dp
       Ntotal_demand            = 0.0_dp
       SMINN_herbivory          = 0.0_dp
       N2O_emissions_mineraliz  = 0.0_dp
       N2O_emissions_slow       = 0.0_dp
       N2O_emissions_grazing    = 0.0_dp
       Nflx_2_crop_harvest      = 0.0_dp
    END IF
    IF (PRESENT(frac_litter_wood_new)) THEN
       litter_wood_new = 0.0_dp
    END IF

    ! Preparations
    !-------------------------------------------

    alpha_mod = max(alpha_min,(alpha-alpha_critical)/(1._dp-alpha_critical)) !! modified soil moist. such that heterotrophic
                                                                             !!     respiration is zero below alpha_critical

    ! Initializations
    NPP_flux_correction   = 0.0_dp
    Cpool_green_max       = 0.0_dp
    Cpool_reserve_optimal = 0.0_dp
    excess_NPP            = 0.0_dp
    root_exudates         = 0.0_dp
    Cflx_litterTotal      = 0.0_dp

    leaf_shedding_debug      = 0.0_dp ! HW
    Cpool_green_minus_shed_leaves     = 0.0_dp ! HW
    excess_carbon_debug      = 0.0_dp ! HW
    litter_leaf      = 0.0_dp ! HW

    Cflx_herbivory        = 0.0_dp
    Cflx_herbivory_LG     = 0.0_dp
    Cflx_herbivory_2_atm  = 0.0_dp
    Cflx_2_crop_harvest   = 0.0_dp
    Cflx_crop_harvest_2_atm = 0.0_dp

    IF (with_yasso) THEN
       Yasso_io_pools(1)=YCpool_acid_ag1
       Yasso_io_pools(2)=YCpool_water_ag1
       Yasso_io_pools(3)=YCpool_ethanol_ag1
       Yasso_io_pools(4)=YCpool_nonsoluble_ag1
       Yasso_io_pools(5)=YCpool_acid_bg1
       Yasso_io_pools(6)=YCpool_water_bg1
       Yasso_io_pools(7)=YCpool_ethanol_bg1
       Yasso_io_pools(8)=YCpool_nonsoluble_bg1
       Yasso_io_pools(9)=YCpool_humus_1
       Yasso_io_pools(10)=YCpool_acid_ag2
       Yasso_io_pools(11)=YCpool_water_ag2
       Yasso_io_pools(12)=YCpool_ethanol_ag2
       Yasso_io_pools(13)=YCpool_nonsoluble_ag2
       Yasso_io_pools(14)=YCpool_acid_bg2
       Yasso_io_pools(15)=YCpool_water_bg2
       Yasso_io_pools(16)=YCpool_ethanol_bg2
       Yasso_io_pools(17)=YCpool_nonsoluble_bg2
       Yasso_io_pools(18)=YCpool_humus_2

       Weather(1) = temp2_30d
       Weather(2) = precip_30d

       LeafLit_coefV(1) = LeafLit_coef1
       LeafLit_coefV(2) = LeafLit_coef2
       LeafLit_coefV(3) = LeafLit_coef3
       LeafLit_coefV(4) = LeafLit_coef4
       LeafLit_coefV(5) = LeafLit_coef5
       WoodLit_coefV(1) = WoodLit_coef1
       WoodLit_coefV(2) = WoodLit_coef2
       WoodLit_coefV(3) = WoodLit_coef3
       WoodLit_coefV(4) = WoodLit_coef4
       WoodLit_coefV(5) = WoodLit_coef5

       Yasso_out(1:16) = 0.0_dp           ! array to store the Yasso pools/fluxes 
       leafLitter = 0.0_dp                ! amount of leafLitter 
       sizeZero = 0.0_dp
    END IF

    if (is_vegetation) then      !! only for covertypes with positive specific leaf area computations are meaningful

       Cpool_green_max = greenC2leafC * LAI / specific_leaf_area_C
       Cpool_reserve_optimal = reserveC2leafC * max_LAI / specific_leaf_area_C

       !! ============== NON-N-LIMITED CARBON AND NITROGEN ALLOCATION =========================================================
       !!
       !! Perform those C- and N-allocation steps that are not restricted by N-availabilty
       !! 1. Leaf shedding
       !! 2. Wood shedding
       !! 3. Depletion of the reserve pool
       !! 4. Decomposition of slow soil pool (whether this is really independent of N-availability is debatable!!)
       !!
       !! 
       !! ----------------------------------------------------------------------------
       !! 1. Leaf shedding: Transfer of C from green pool to leaf litter litter pool

!!$ tr
!!$       IF (LAI >= LAI_previous) THEN                              !! If LAI is increasing or constant the leaf shedding is given
!!$          leaf_shedding = LAI * LAI_shed_constant                 !!    by a constant loss of leaves for evergreens, raingreens
!!$                                                                  !!    and grasses.
!!$       ELSE                                                       !! Otherwise
!!$          leaf_shedding = LAI_previous - LAI                      !!    the leaf shedding is given by the decrease in LAI.
!!$       END IF
!!$ tr

! ------HW
       leaf_shedding=leaf_shedding_rate*LAI_previous*1.0_dp
       leaf_shedding_debug=leaf_shedding
! HW---- 
!       leaf_shedding = MAX(LAI * LAI_shed_constant,LAI_previous - LAI)  !! Leaf shedding is assured at a minimum constant rate 
!                                                                        !! which maybe exceeded by the decrease in LAI.       

       C_2_litter_greenPool =  &                                     !! Maximally the whole carbon from the green pool is transfered
          min(2.0_dp * leaf_shedding / specific_leaf_area_C, & !!    by leaf and root shedding to the green litter pool
          Cpool_green)                                               


       Cpool_green_minus_shed_leaves = Cpool_green- 2.0_dp * leaf_shedding /specific_leaf_area_C      ! HW

!       Cpool_green = Cpool_green - C_2_litter_greenPool              !! This removes carbon from the green pool

!!       Cpool_green = max(Cpool_green - greenC2leafC * (LAI_previous - LAI) / specific_leaf_area_C,0.0_dp)     !! HW 
         !! This removes net carbon change from the green pool HW

       Cpool_green = max(Cpool_green - 2.0_dp * leaf_shedding / specific_leaf_area_C,0.0_dp)     !! HW
        !! This removes litter carbon from the green pool HW
       excess_carbon = max(0.0_dp,Cpool_green - Cpool_green_max)     !! If by the reduction of LAI the maximum value of the green
                                                                     !!    pool is still smaller than the current value, there is 

       excess_carbon_debug=excess_carbon     ! HW

       C_2_litter_greenPool = C_2_litter_greenPool + excess_carbon   !!    excess C that has also to be shedded to the green litter
       Cpool_green = Cpool_green - excess_carbon                     !!    pool and also substracted from the green pool

       IF (is_crop) THEN
          !!  Move shedded Carbon pro rata to crop harvest flux and to ag green litter pools
          Cflx_2_crop_harvest = frac_green_aboveGround * frac_C_crop_harvest * C_2_litter_greenPool
          Cpool_litter_green_ag = Cpool_litter_green_ag + &
             frac_green_aboveGround * (1._dp - frac_C_crop_harvest) * C_2_litter_greenPool
       ELSE
          !! Move shedded Carbon to ag green litter pools
          Cpool_litter_green_ag = Cpool_litter_green_ag + frac_green_aboveGround * C_2_litter_greenPool
       END IF
       !! Move shedded Carbon to bg green litter pools
       Cpool_litter_green_bg = Cpool_litter_green_bg + (1._dp - frac_green_aboveGround) * C_2_litter_greenPool
       
       IF (with_Nitrogen) THEN

          Npool_green = Npool_green - C_2_litter_greenPool / cn_green   !! Shed N according to C/N-ratio of green litter
          IF (is_crop) THEN
             Nflx_2_crop_harvest = frac_green_aboveGround * frac_C_crop_harvest * C_2_litter_greenPool / cn_litter_green
             Npool_crop_harvest = Npool_crop_harvest + Nflx_2_crop_harvest
             Nflx_2_crop_harvest = Nflx_2_crop_harvest / sec_per_day ! flux from day-1 to s-1
             Npool_litter_green_ag = Npool_litter_green_ag + &
                 frac_green_aboveGround * (1._dp - frac_C_crop_harvest) * C_2_litter_greenPool / cn_litter_green
          ELSE
             Npool_litter_green_ag = &
                 Npool_litter_green_ag + frac_green_aboveGround * C_2_litter_greenPool / cn_litter_green
          END IF
          Npool_litter_green_bg = &
              Npool_litter_green_bg + (1._dp - frac_green_aboveGround) * C_2_litter_greenPool / cn_litter_green

          !! retranslocation the surplus N from shedded green parts to mobile N   ("leaf-N retranslocation")             
          Npool_mobile = Npool_mobile + (1._dp/cn_green - 1._dp/cn_litter_green) * C_2_litter_greenPool  !! makes CN constant
                                                                                                         !!    b/w green-litterpool
       END IF

       Cflx_litterTotal = C_2_litter_greenPool
       litter_leaf = C_2_litter_greenPool ! HW

       !! ---------------------------------------------------------------------------------
       !! Herbivory loss from Cpool_green & Npool_green w.r.t PFTs in lctlib

       Green_2_herbivory  = Cpool_green * frac_green_2_herbivory             !! Grazing flux indirect comput.from NPP losses 
       Cpool_green        = Cpool_green - Green_2_herbivory                  !! Loss due to grazing

       Cflx_faeces_2_LG  = Green_2_herbivory * frac_C_faeces2_LG             !! Cflux_faeces_2_LG to be stored in litter green pool
                                                                             !!      (not considered in Cflx_litterTotal) 
       Cflx_faeces_2_atm = Green_2_herbivory * (1.0_dp - frac_C_faeces2_LG)  !! Cflux_faeces_2_atm to be lost to atmosphere 
                                                                             !!      (dk: not added to CO2?)
       Cpool_litter_green_ag = Cpool_litter_green_ag + Cflx_faeces_2_LG      !! Cflux_faeces_2_LG put into LG pool  

       IF (with_Nitrogen) THEN
          Npool_green = Npool_green - Green_2_herbivory /cn_green

          !! It is assumed that there is no N loss to the atmosphere due to grazing. All N goes to the animals faeces according
          !! to the C/N-ration of green litter.
          !! The N of faeces corresponding to the Cflx_faeces_2_LG is assigned to the above ground green litter pool, whereas N
          !! corresponding to respired faeces is assigned to the mineral N pool.
          Npool_litter_green_ag = Npool_litter_green_ag + Cflx_faeces_2_LG /cn_litter_green

          !! retranslocation the surplus N from faeces of grazing animals (urine) to SMINN pool
          !! note: C/N-ratio of herbivory pool is equal to C/N-ratio of green pool as there is no N loss assumed
          Nflx_faeces_2_SMINN =   1._dp/cn_green                          * Cflx_faeces_2_atm &
                               + (1._dp/cn_green - 1._dp/cn_litter_green) * Cflx_faeces_2_LG      
          SMINN_pool = SMINN_pool + (1._dp - N2O_ef_grazing) * Nflx_faeces_2_SMINN

          !! N2O emissions from mineral N released due to grazing; first assumption: EF similar to fertilizer application
          N2O_emissions_grazing =  N2O_ef_grazing * Nflx_faeces_2_SMINN

          !! diagnostic output: sminN gain from herbivory faeces and dungs
          SMINN_herbivory = SMINN_herbivory + (1._dp - N2O_ef_grazing) * Nflx_faeces_2_SMINN

          !! conversion from flux per day to flux per second
          SMINN_herbivory = SMINN_herbivory / sec_per_day
       END IF

       !! ---------------------------------------------------------------------------------
       !! 2. Wood shedding: Transfer of C from wood pool to wood litter pools
       !! 
                                                            !! Assuming that forests continously die at the inverse lifetime of trees
       C_2_litter_woodPool = Cpool_woods / tau_Cpool_woods              !! .. the shedded wood (MAX FUNCTION FOR NONWOODY PFTS)
       Cpool_litter_wood_ag = Cpool_litter_wood_ag &                    !! .. is put partly into the above ground woody litter pool
            + frac_wood_aboveGround*C_2_litter_woodPool
       if (present(frac_litter_wood_new)) &
          litter_wood_new = frac_wood_aboveGround * C_2_litter_woodPool  
       Cpool_litter_wood_bg = Cpool_litter_wood_bg &                    !! .. and the rest into the below ground woody litter pool
            + (1.0_dp-frac_wood_aboveGround)*C_2_litter_woodPool
       Cpool_woods = Cpool_woods -  C_2_litter_woodPool                 !! .. and then all is substracted from the wood pool.

       Cflx_litterTotal = Cflx_litterTotal + C_2_litter_woodPool

       if(with_Nitrogen) then 
          !! Handle the associated organic N-fluxes
          Npool_woods = Npool_woods - C_2_litter_woodPool / cn_woods            !! N removed from wood pool..
          Npool_litter_wood_ag = Npool_litter_wood_ag &                         !! .. is partly put to the ag woody N litter pool..
               + frac_wood_aboveGround*C_2_litter_woodPool /cn_litter_wood
          Npool_litter_wood_bg = Npool_litter_wood_bg &                         !! .. and the rest to the bg woody N litter pool..
               + (1.0_dp-frac_wood_aboveGround)*C_2_litter_woodPool /cn_litter_wood
          SMINN_pool = SMINN_pool &                                             !! .. while the soil mineral N-pool gains the 
               + (1._dp/cn_woods - 1._dp/cn_litter_wood) * C_2_litter_woodPool  !!    N-difference
          !! ASSUMPTION: cn_woods <= cn_litter_wood
       end if

       !! ----------------------------------------------------------------------------
       !! 3. Depletion of reserve pool: Transfer of C from reserve pool to green litter pool
       !!    Some organic carbon of the reserve pool is always lost with mortality of plants. 
       !!    Note that the reserve pool contains no Nitrogen (starches and sugar are free of N).
       !!    Note also that the green litter pool has no fixed C/N-ratio (therefore no compensation
       !!    flux is needed)

       C_2_litter_greenPool = Cpool_reserve / tau_Cpool_reserve
       IF (is_crop) THEN
          Cflx_2_crop_harvest = Cflx_2_crop_harvest + &
              frac_green_aboveGround * frac_C_crop_harvest * C_2_litter_greenPool
          Cpool_litter_green_ag = Cpool_litter_green_ag + &
              frac_green_aboveGround * (1._dp - frac_C_crop_harvest) * C_2_litter_greenPool
       ELSE
          Cpool_litter_green_ag = Cpool_litter_green_ag + frac_green_aboveGround * C_2_litter_greenPool
       END IF
       Cpool_litter_green_bg = Cpool_litter_green_bg + (1._dp - frac_green_aboveGround) * C_2_litter_greenPool
       Cpool_reserve = Cpool_reserve - C_2_litter_greenPool
       Cpool_crop_harvest = Cpool_crop_harvest + Cflx_2_crop_harvest

       Cflx_litterTotal = Cflx_litterTotal + C_2_litter_greenPool

       !! ----------------------------------------------------------------------------
       !! 4. Decomposition of slow soil pool (i) C-flux (ii) Mineral N-flux
       !!
       d_Cpool_slow = &
            -MIN(alpha_mod**kappa * q10**((topSoilTemp - referenceTemp_Q10)/10._dp)  / & !! Slow soil respiration
            tau_Cpool_slow,1._dp)  * Cpool_slow                                          !! .. according to Q10-model
       Cpool_slow = Cpool_slow + d_Cpool_slow                 !! Remove respired Carbon (note the negative sign!)

       Cflx_slow_2_atmos            = -d_Cpool_slow          !! Remember the decomposition flux to the atmosphere (CO2-release)


       !!
       !! ============== START OF POTENTIAL PLANT-CARBON ALLOCATION TO ESTIMATE PLANT-N-DEMAND =================================
       !!                       (exception: case NPP<0 handles actual allocation)

       Cpool_green_pot   = Cpool_green
       Cpool_woods_pot   = Cpool_woods
       Cpool_reserve_pot = Cpool_reserve

       !! Preliminary determination of distribution of NPP to the various pools (may be corrected afterwards according to 
       !!    avaialable Nitrogen)

       if(NPPrate <= 0.0_dp) then                                      !! In case of negative or zero NPP
          NPP_2_reservePool =  NPPrate  * sec_per_day                  !! .. the plant looses C and  
          Cpool_reserve     = Cpool_reserve + NPP_2_reservePool        !! .. we try to take it from the reserve pool 
          if(Cpool_reserve < 0.0_dp) then                              !! But if thereby the reserve pool gets negative 
             NPP_2_reservePool = NPP_2_reservePool - Cpool_reserve     !! .. we absorb negative NPP only up to zero reseve pool and 
             NPP_flux_correction = -Cpool_reserve / sec_per_day        !! .. give the deficit as a flux correction back to the
             Cpool_reserve = 0.0_dp                                    !! .. calling routine. Hence all Carbon from the reserve pool
          end if                                                       !! .. is used up.

          !! For correct handling of nitrogen limitation below (which handles only the case of positive NPP) set NPP to zero:
          NPP_2_reservePool = 0.0_dp
          NPP_2_greenPool   = 0.0_dp   !! All other transfer rates are zero
          NPP_2_woodPool    = 0.0_dp
          NPP_2_rootExudates= 0.0_dp
          Cpool_reserve_pot = Cpool_reserve      !! remember modified reserve pool


       else  !! NPP is positive                                                            

          NPP_2_woodPool    = frac_npp_2_woodPool * NPPrate * sec_per_day    !! NPP it is distributed according to predefined 
                                                                             !!    relative fraction
          NPP_2_reservePool = frac_npp_2_reservePool * NPPrate * sec_per_day
          NPP_2_rootExudates= frac_npp_2_exudates * NPPrate * sec_per_day
          NPP_2_greenPool   = (1.0_dp - frac_npp_2_woodPool - frac_npp_2_reservePool                 & 
               - frac_npp_2_exudates) * NPPrate * sec_per_day

          !! Growth of Wood Pool (structural carbon of living plants)

          Cpool_woods_pot = Cpool_woods_pot + NPP_2_woodPool                 !! Then it is attempted to put all NPP share into it
          excess_carbon = max(0.0_dp,Cpool_woods_pot-Max_C_content_woods)    !! .. but thereby the pool may get too large
          Cpool_woods_pot = Cpool_woods_pot - excess_carbon                  !! .. so that the excess carbon has once more to be 
                                                                             !! .. subtracted
          NPP_2_greenPool = NPP_2_greenPool + excess_carbon                  !! .. and instead made available to the green pool
          NPP_2_woodPool = NPP_2_woodPool - excess_carbon                    !! .. and the actual amount of NPP put to the wood pool
                                                                             !! .. is less


          !! Growth of Reserve Pool (sugars and starches) and determination how much will enter the Green Pool
          if(Cpool_reserve_pot < Cpool_reserve_optimal) then                 !! .. If the reserve pool is smaller than optimal,
             Cpool_reserve_pot = Cpool_reserve_pot +  NPP_2_reservePool      !! .... it is filled by the available NPP
             excess_carbon =                                         &       !! .... Thereby it may happen that it gets larger than
                  max(0.0_dp,Cpool_reserve_pot - Cpool_reserve_optimal)      !! .... optimal so that there is excess carbon
             Cpool_reserve_pot = Cpool_reserve_pot - excess_carbon           !! .... that needs not be taken up
             NPP_2_greenPool = NPP_2_greenPool + excess_carbon               !! .... but can better be used to increase the green
             NPP_2_reservePool = NPP_2_reservePool - excess_carbon           !! .... pool and the actual amount of NPP put to the
                                                                             !! .... reserve pool is less.
          else                                                               !! .... Otherwise (reserve pool is larger than optimal)
             NPP_2_greenPool = NPP_2_greenPool + NPP_2_reservePool           !! .... all NPP is left for the green pool
             NPP_2_reservePool = 0.0_dp                                      !! .... so that nothing is stored in the reserve pool.
          end if

          !! Growth of Green Pool (leaves and fine roots): (in case of too much NPP, try to put it into the reserve pool).

          Cpool_green_pot = Cpool_green_pot + NPP_2_greenPool                !! Green pool is filled by the available NPP.
          excess_carbon = max(0.0_dp,Cpool_green_pot - Cpool_green_max)      !! .. Thereby it may get larger than appropriate for 
                                                                             !!    current LAI.
          NPP_2_greenPool = NPP_2_greenPool - excess_carbon                  !! .. Hence the actual amount of NPP put to the green
                                                                             !!    pool is less.
          Cpool_green_pot = Cpool_green_pot - excess_carbon                  !! .. This excess carbon needs not be taken up but
          if(Cpool_reserve_pot < Cpool_reserve_optimal) then                 !! .... if the reserve pool is smaller than optimal
             Cpool_reserve_pot = Cpool_reserve_pot + excess_carbon           !! .... it is tried to put the carbon there,
             NPP_2_reservePool = NPP_2_reservePool + excess_carbon           !! .... which means that additional NPP is put to the 
                                                                             !!      reserve pool.
             excess_carbon =                                          &      !! .... Thereby it may happen that the reserve pool
                  max(0.0_dp,Cpool_reserve_pot - Cpool_reserve_optimal)      !!      increases beyond the optimal value.
             Cpool_reserve_pot = Cpool_reserve_pot - excess_carbon           !! .... In that case the excess carbon is once more
                                                                             !!      removed from the reserve pool
             NPP_2_reservePool = NPP_2_reservePool - excess_carbon           !! .... so that the actual amount of NPP put into the
                                                                             !!      reserve pool is less.
          end if                                                             !!

          Cpool_litter_green_bg = Cpool_litter_green_bg + NPP_2_rootExudates
          ! THT 30.10.2012: This taken into account in Yasso too.. Put NPP_2_rootExudates to Yasso-water-pool. Added 20.11.2012
          ! Yasso_io_pools(2) = Yasso_io_pools(2) + NPP_2_rootExudates ! DSG: exudates goes into yasso

          excess_NPP = excess_carbon / sec_per_day

          root_exudates = NPP_2_rootExudates / sec_per_day

       end if !! NPP > 0 end


       !! ============== START OF POTENTIAL LITTER-CARBON ALLOCATION ==============================================

       Cpool_litter_green_ag_pot   = Cpool_litter_green_ag
       Cpool_litter_green_bg_pot   = Cpool_litter_green_bg
       Cpool_litter_wood_ag_pot    = Cpool_litter_wood_ag
       Cpool_litter_wood_bg_pot    = Cpool_litter_wood_bg

       ! Update Green Litter Pools (litter from dead leaves, fruits, debris from bark and fine roots)

       decompRate_litter_green =  &                                     !! decomposition rate of green litter according to Q10 model
            -MIN(alpha_mod**kappa * q10**((topSoilTemp - referenceTemp_Q10)/10._dp)/ tau_Cpool_litter_green,1._dp)

       d_Cpool_litter_green_ag      = decompRate_litter_green * Cpool_litter_green_ag_pot     !! Decomposition of green plant litter
       Cpool_litter_green_ag_pot    = Cpool_litter_green_ag_pot + d_Cpool_litter_green_ag
       d_Cpool_litter_green_bg      = decompRate_litter_green * Cpool_litter_green_bg_pot     !! Decomposition of green plant litter
       Cpool_litter_green_bg_pot    = Cpool_litter_green_bg_pot + d_Cpool_litter_green_bg

       ! Update Wood Litter Pools (litter from all dead woody material, above and below ground)

       d_Cpool_litter_wood_ag = -MIN(q10**((topSoilTemp - referenceTemp_Q10)/10._dp)  / &   !! Decomposition of woody litter
            tau_Cpool_litter_wood,1._dp) * Cpool_litter_wood_ag_pot                         !! .. according to Q10-model
       d_Cpool_litter_wood_bg = -MIN(q10**((topSoilTemp - referenceTemp_Q10)/10._dp)  / &   !! Decomposition of woody litter
            tau_Cpool_litter_wood,1._dp) * Cpool_litter_wood_bg_pot                         !! .. according to Q10-model
       !! .. without soil moisture dependence
       Cpool_litter_wood_ag_pot = Cpool_litter_wood_ag_pot + d_Cpool_litter_wood_ag
       Cpool_litter_wood_bg_pot = Cpool_litter_wood_bg_pot + d_Cpool_litter_wood_bg


       !! ============== START OF POTENTIAL ALLOCATION OF SLOW SOIL POOL ==============================================


       ! Update Slow Soil Pool (uptake of C from litter pools, decomposition has already been handled above)
       Cpool_slow_pot = Cpool_slow  
       Cpool_slow_pot = Cpool_slow_pot                                                           &
            - (1._dp - frac_C_litter_wood2atmos) * (d_Cpool_litter_wood_ag + d_Cpool_litter_wood_bg)  &
            - (1._dp - frac_C_litter_green2atmos) * (d_Cpool_litter_green_ag + d_Cpool_litter_green_bg)           



       !! ============== use nitrogen from N mobile pool to satisfy (at least partly) for potential carbon fluxes ===
       ! Remarks: Nmobile_2_green + Nmobile_2_wood <= Npool_mobile
       ! Npool_mobile accounts minN_plant_demand by trasfering N into green and wood pool 

       Cflx_NPP_2_green_pot          = NPP_2_greenPool                                       !! potential growth of green carbon 
       Cflx_NPP_2_wood_pot           = NPP_2_woodPool                                        !! potential growth of wood carbon 

       NPP_2_green_via_mobileN = 0.0_dp        
       NPP_2_wood_via_mobileN = 0.0_dp 

       IF (with_Nitrogen) THEN

          Nmobile_2_green_pot = Cflx_NPP_2_green_pot/cn_green
          Nmobile_2_wood_pot  = Cflx_NPP_2_wood_pot /cn_woods 

          !... CN mode 
          if (NPPrate > 0.0_dp) then                                   !! (Npool_mobile> 0) & avoid Npool_mobile negative 
             if(Nmobile_2_green_pot + Nmobile_2_wood_pot <= Npool_mobile) then
                Nmobile_2_green = Nmobile_2_green_pot
                Nmobile_2_wood  = Nmobile_2_wood_pot
                NPP_2_green_via_mobileN= Cflx_NPP_2_green_pot         
                NPP_2_wood_via_mobileN = Cflx_NPP_2_wood_pot    
             else
                redFact_Nmobile = Npool_mobile/(Nmobile_2_green_pot + Nmobile_2_wood_pot)  !! avoid Npool_mobile < 0; NPPrate > 0
                Nmobile_2_green = redFact_Nmobile * Nmobile_2_green_pot
                Nmobile_2_wood  = redFact_Nmobile * Nmobile_2_wood_pot
                NPP_2_green_via_mobileN= Nmobile_2_green * cn_green
                NPP_2_wood_via_mobileN = Nmobile_2_wood  * cn_woods
             end if

             !! update C-N green and wood pools 
             Cpool_green = Cpool_green + NPP_2_green_via_mobileN
             Cpool_woods = Cpool_woods + NPP_2_wood_via_mobileN
             Npool_green = Npool_green + Nmobile_2_green
             Npool_woods = Npool_woods + Nmobile_2_wood

             Npool_mobile = max(0.0_dp, Npool_mobile - Nmobile_2_green - Nmobile_2_wood)           !!  update Npool_mobile
          end if
       END IF

       Cflx_NPP_2_green_pot = Cflx_NPP_2_green_pot - NPP_2_green_via_mobileN
       Cflx_NPP_2_wood_pot  = Cflx_NPP_2_wood_pot  - NPP_2_wood_via_mobileN

       Nflux_NPP_2_green_pot = Cflx_NPP_2_green_pot/cn_green
       Nflux_NPP_2_wood_pot  = Cflx_NPP_2_wood_pot/cn_woods


       !! ========= COMPUTE MINERAL-PLANT N DEMAND =================

       minN_plant_demand = Nflux_NPP_2_green_pot + Nflux_NPP_2_wood_pot 


       !! ============== COLLECT POTENTIAL C-FLUXES ===============================================================

       Cflx_NPP_2_reserve_pot        = NPP_2_reservePool                                     !! potential growth of reserve carbon
       Cflx_litterGreen_ag_2_atmos_pot = -frac_C_litter_green2atmos * d_Cpool_litter_green_ag !! ag green litter decomp. -> atm.
       Cflx_litterGreen_bg_2_atmos_pot = -frac_C_litter_green2atmos * d_Cpool_litter_green_bg !! bg green litter decomp. -> atm.
       Cflx_litterGreen_ag_2_slow_pot  = &
           - (1._dp -frac_C_litter_green2atmos) * d_Cpool_litter_green_ag                    !! ag green litter decomp. -> slow pool
       Cflx_litterGreen_bg_2_slow_pot  = &
            -(1._dp -frac_C_litter_green2atmos) * d_Cpool_litter_green_bg                    !! bg green litter decomp. -> slow pool
       Cflx_litter_wood_ag_2_atmos_pot  = -frac_C_litter_wood2atmos * d_Cpool_litter_wood_ag !! ag-wood litter decomp.  -> atm.
       Cflx_litter_wood_bg_2_atmos_pot  = -frac_C_litter_wood2atmos * d_Cpool_litter_wood_bg !! bg-wood litter decomp.  -> atm.
       Cflx_litter_wood_ag_2_slow_pot   = &
            -(1._dp - frac_C_litter_wood2atmos) * d_Cpool_litter_wood_ag                     !! ag-wood litter decomp.  -> slow pool
       Cflx_litter_wood_bg_2_slow_pot   = &
            -(1._dp - frac_C_litter_wood2atmos) * d_Cpool_litter_wood_bg                     !! bg_wood litter decomp.  -> slow pool

       if(with_Nitrogen) then
          IF (with_yasso) THEN

            ! 0. reset all variables which were set by CTL already
            Cflx_slow_2_atmos               = 0.0_dp

            ! 1. compute potential soil decomposition fluxes 
             !! 1.1 call yasso for leaf litter 
             CALL yasso(Yasso_io_pools(1:9), Weather, leafLitter, LeafLit_coefV,               & 
                              sizeZero, Yasso_out, frac_green_aboveGround, NPP_2_rootExudates)

            ! 2.1 Update only the varibles which are needed
              Cflx_litterGreen_ag_2_slow_pot  = Yasso_out(11) 
              Cflx_litterGreen_bg_2_slow_pot  = Yasso_out(12) 
              Cflx_slow_2_atmos               = Yasso_out(13)
              Cflx_litterGreen_ag_2_atmos_pot = Yasso_out(14)
              Cflx_litterGreen_bg_2_atmos_pot = Yasso_out(15)

              ! As Litter green Pool has a flexible CN ratio we need the
              ! decomposition rate to derive the mineralization flux;
              decompRate_litter_green = Yasso_out(16) ! should be negative

              ! update potential soil respiration rate
              soilResp_rate_pot =  Yasso_out(10) / sec_per_day 
 
             !! 1.2 call yasso for woody litter
             CALL yasso(Yasso_io_pools(10:18), Weather, C_2_litter_woodPool, WoodLit_coefV,       & 
                           WoodLitterSize, Yasso_out, frac_wood_aboveGround, 0.0_dp) ! no exudates

            ! 2.2 Update only varibles which are needed
              Cflx_litter_wood_ag_2_slow_pot  = Yasso_out(11)
              Cflx_litter_wood_bg_2_slow_pot  = Yasso_out(12)
              Cflx_slow_2_atmos               = Cflx_slow_2_atmos + Yasso_out(13) ! don't forget the leaf humus from call before
              Cflx_litter_wood_ag_2_atmos_pot = Yasso_out(14)
              Cflx_litter_wood_bg_2_atmos_pot = Yasso_out(15)
              ! update potential soil respiration rate; add respiration from
              ! woody litter:
              soilResp_rate_pot =  soilResp_rate_pot  + Yasso_out(10) / sec_per_day 
          ENDIF


          minNflux_slow = Cflx_slow_2_atmos/cn_slow             !! Mineral N-flux from slow pool to soil mineral N-pool

          !! Handle the associated organic N-flux
          !DSG: we need to take the minNflux_slow here: Npool_slow = max(0.0_dp,Npool_slow + d_Cpool_slow / cn_slow)  !! check
          ! in case YASSO is active d_Cpool_slow/cn_slow would be wrong 
          Npool_slow = max(0.0_dp, Npool_slow - minNflux_slow)

          !! part of the minNflux_slow flux emits as N2O gas  due to nitrification & denitrification
          !! by considering implicit NH4 & NO3 availability in soil.
          !! with    Temperature optimum at 38. degC Xu-Ri (2006)
          !!          alpha (soil moisture) for soil moisture regulation of denitrification
          N2O_emissions_slow = (minNflux_slow * sminn_NH4_fraction * N2O_rate_nitrification                      &
                              + minNflux_slow * (1._dp - sminn_NH4_fraction) * N2O_rate_denitrification * alpha) &
                              * MAX(0._dp, MIN(topSoilTemp/38._dp, 1._dp))

          minNflux_slow = minNflux_slow - N2O_emissions_slow        !! Remained N-flux after N2O loss to atm.
          SMINN_pool    = SMINN_pool + minNflux_slow                !! Put released mineral N to SMINN_pool

       end if

       IF (.NOT. with_Nitrogen) THEN

          redFact_Nlimitation = 1.0_dp   !! This assures absence of nitrogen limitation in this mixed C-N code

       ELSE

          !! ============== DERIVE ASSOCIATED POTENTIAL FLUXES OF ORGANIC N ==============================================
          
          Nflux_litterGreen_ag_2_slow_pot = Cflx_litterGreen_ag_2_slow_pot / cn_slow !! N-content determined by N-sink 
          Nflux_litterGreen_bg_2_slow_pot = Cflx_litterGreen_bg_2_slow_pot / cn_slow !! N-content determined by N-sink 
          Nflux_litter_wood_ag_2_slow_pot  = Cflx_litter_wood_ag_2_slow_pot  / cn_slow !! N-content determined by N-sink 
          Nflux_litter_wood_bg_2_slow_pot  = Cflx_litter_wood_bg_2_slow_pot  / cn_slow !! N-content determined by N-sink 

          !! ============== DETERMINE POTENTIAL FLUXES OF MINERAL N ==============================================
          !
          ! Remark: Mineral-N-fluxes are considerd positive when they increase the SMINN-pool
          !
          ! How to derive the mineral N-fluxes:
          ! -----------------------------------
          !
          ! Mineral fluxes are determined as a kind of residues: assuming different but fixed C/N-ratios for the organic pools these
          ! could not be maintained without compensating Nitrogen fluxes: Therefore the strategy is to compute the mineral N-fluxes
          ! such that they assure N-mass conservation. 
          !
          ! Derivation of mineral N-flux (also called "immobilization flux"):
          !
          ! Consider a pool of organic carbon with N-content N and C-content C. Assume an influx of (organic) C called c1 associated
          ! with a N-influx called n1. Assume further two Carbon outfluxes, an organic one called c2 associated with a N-outflux n2,
          ! and a mineral Carbon outflux called c0.
          ! In addition there is a mineral N outflux m. This mineral flux m is coming from the mineral N-pool that must compensate
          ! for potential C/N-deviations. The aim is to determine m:
          !
          ! The change in C and N is then given by
          !
          !      dC     
          ! (1)  -- = c1 - c2 -c0
          !      dt 
          !
          !      dN   
          ! (2)  -- = n1 - n2 - m
          !      dt   
          !
          ! It is further assumed that the N-to-C ratios nc11 and nc2 for in- and outfluxes, respectively, is fixed:
          !
          ! (3)  n1 = nc1 * c1,   n2 = nc2 * c2
          !
          ! Combining the last 2 equations and solving for the compensating mineral flux gives:
          !
          !                                dN
          ! (4)  m = nc1 * c1 - nc2 * c2 - --
          !                                dt
          !
          ! If the N/C ratio of the considered pool shall be kept constant under the prescribed fluxes c1,c2,n1,n2, i.e. 
          !
          !           N   dN      
          ! (5) nc := - = -- 
          !           C   dC
          !
          ! so that
          !
          !      dN        dC
          ! (6)  -- = nc * --
          !      dt        dt
          !
          ! equations (4), (6) and (1) give
          !
          ! (7)  m = nc1 * c1 - nc2 * c2 - nc*(c1 - c2 -c0) = (nc1 - nc)*c1 + (nc - nc2)*c2 + nc*c0.
          !
          ! If in addition it is known that the process by which the organic C-pool is depleted is determined by the rate equation
          !
          !      dC        
          ! (8)  -- = - r * C
          !      dt        
          !
          ! an alternative equation to (7) can be derived: First, using (5) in (8) shows that N is depleted by the same rate:
          !
          !      dN        
          ! (9)  -- = - r * N
          !      dt        
          !
          ! and entering this into (4) gives the alternative result
          !
          ! (10) m = nc1 * c1 - nc2 * c2 + r*N.
          !
          ! This latter form is especially useful when C and N vary independently so that the N/C-ratio changes in time. (Note that
          ! this is no contradiction to (5), because other processes than the depletion (8) may change N/C-ratio, e.g. influxes of
          ! C and N with different N/C-ratio.)
          !
          ! ---------------------------------------------------------------
          ! Now Compute the potential mineral N fluxes (Note: fluxes into SMINN-pool are positive) 

          minNflux_litter_green_ag_pot =    &      !! ag green litter -> SMINN (use eq. (10), because green litter C/N is not fixed)
               - Nflux_litterGreen_ag_2_slow_pot - decompRate_litter_green * Npool_litter_green_ag

          minNflux_litter_green_bg_pot =    &      !! ag green litter -> SMINN (use eq. (10), because green litter C/N is not fixed)
               - Nflux_litterGreen_bg_2_slow_pot - decompRate_litter_green * Npool_litter_green_bg

          minNflux_litter_wood_ag_pot =    &  !! above ground woody litter -> SMINN (use eq. (7), because woody litter C/N is fixed)
               (1._dp/cn_litter_wood - 1._dp/cn_slow) * Cflx_litter_wood_ag_2_slow_pot &
               + Cflx_litter_wood_ag_2_atmos_pot/cn_litter_wood 

          minNflux_litter_wood_bg_pot =    &  !! below ground woody litter -> SMINN (use eq. (7), because woody litter C/N is fixed)
               (1._dp/cn_litter_wood - 1._dp/cn_slow) * Cflx_litter_wood_bg_2_slow_pot &
               + Cflx_litter_wood_bg_2_atmos_pot/cn_litter_wood

          !! ============== COMPUTE N2O emissions==============================================
          ! --- Nitros oxide emissions(N2O)--N2O_Emissions() due to N mineralization ----
          !! N2O emissions rate per day (<0.1%-0.2% nitrification ); 0.2-4.7% due to denit in Xu Ri, 2008 * therein

          !! minNflux_litter can be negativ, as cn_slow much smaller as cn_litter_wood and release of N in decomp not always 
          !!       sufficient

          IF ( (minNflux_litter_green_ag_pot + minNflux_litter_green_bg_pot + minNflux_litter_wood_ag_pot &
                   + minNflux_litter_wood_bg_pot) .gt. 0._dp ) THEN
            !!if N is released in litter decomposition

            ! dk: new calculation of N2O from litter:
            !  above ground: treated as fertilizer: same emission factor
            !  below ground: decreasing SOC with soil depth represented by assuming that 50% of SOC are within top soil layers
            !                relevant for N2O production (usually 30 cm); 
            !  plus: moisture dependency of denitrification: alpha

            minNflux_litter_active = (  minNflux_litter_green_ag_pot + 0.5_dp * minNflux_litter_green_bg_pot  &
                                      + minNflux_litter_wood_ag_pot  + 0.5_dp * minNflux_litter_wood_bg_pot ) &
                                     * MAX(0._dp, MIN(topSoilTemp/38._dp, 1._dp))
            N2O_emissions_litter_pot =          sminn_NH4_fraction  * minNflux_litter_active * N2O_rate_nitrification           &
                                     + (1._dp - sminn_NH4_fraction) * minNflux_litter_active * N2O_rate_denitrification * alpha
          ELSE
            N2O_emissions_litter_pot = 0._dp
          END IF

          !!
          !! ============== COMPUTE MINERAL-N DEMAND ==============================================
          !!

          minN_soil_demand = -( &   !! "-" because of sign convention for fluxes to SMINN-pool. Note: Slow soil pool only releases N
                                 minNflux_litter_green_ag_pot + minNflux_litter_green_bg_pot  &  
                               + minNflux_litter_wood_ag_pot  + minNflux_litter_wood_bg_pot   &
                              )
          !!
          !! ============== DETERMINE REDUCTION FACTOR FROM NITROGEN LIMITATION ========================================
          !!

          ! If minN_soil_demand is negative, N is released to the sminN-pool. It can be considered in total_minN_demand as available N,
          ! BUT litter decay will be reduced by redFAct, and LESS N would be availabe as assumed here,
          ! BUT if negative minN_soil_demand would be excluded here, it would not be added to the sminn_pool lateron

          total_minN_demand = minN_plant_demand + minN_soil_demand + N2O_emissions_litter_pot
          
          ! N2O emissions do not determine soil N demand, however, N2O 
          ! will be subtracted from sminn, hence considered in total demand

          if( total_minN_demand  <= SMINN_pool + epsilon(1.0_dp) ) then     !! enough mineral N available
             redFact_Nlimitation = 1.0_dp
          else
             redFact_Nlimitation = SMINN_pool / total_minN_demand
          end if

          Nplant_demand = max(0.0_dp,minN_plant_demand * redFact_Nlimitation )
          Nsoil_demand  = max(0.0_dp,minN_soil_demand  * redFact_Nlimitation )
          Ntotal_demand = max(0.0_dp,total_minN_demand * redFact_Nlimitation )

       END IF

       !! ============== UPDATE ALL C-POOLS =======================================================================

       Cpool_green       = Cpool_green       + redFact_Nlimitation * Cflx_NPP_2_green_pot
       Cpool_woods       = Cpool_woods       + redFact_Nlimitation * Cflx_NPP_2_wood_pot
       Cpool_reserve     = Cpool_reserve     +  Cflx_NPP_2_reserve_pot        !!dk removed Nlimit, as no dependent on N

       Cpool_litter_green_ag = Cpool_litter_green_ag       &
            - redFact_Nlimitation * (Cflx_litterGreen_ag_2_atmos_pot + Cflx_litterGreen_ag_2_slow_pot)

       Cpool_litter_green_bg = Cpool_litter_green_bg       &
            - redFact_Nlimitation * (Cflx_litterGreen_bg_2_atmos_pot + Cflx_litterGreen_bg_2_slow_pot)

       Cpool_litter_wood_ag = Cpool_litter_wood_ag       &
            - redFact_Nlimitation * (Cflx_litter_wood_ag_2_atmos_pot + Cflx_litter_wood_ag_2_slow_pot)

       Cpool_litter_wood_bg = Cpool_litter_wood_bg       &
            - redFact_Nlimitation * (Cflx_litter_wood_bg_2_atmos_pot + Cflx_litter_wood_bg_2_slow_pot)

       Cpool_slow        = Cpool_slow     &
            + redFact_Nlimitation * (  Cflx_litterGreen_ag_2_slow_pot + Cflx_litterGreen_bg_2_slow_pot &
                                     + Cflx_litter_wood_ag_2_slow_pot  + Cflx_litter_wood_bg_2_slow_pot  )

       IF (with_Nitrogen) THEN

          !! ============== UPDATE ORGANIC-N-POOLS =======================================================================
          !!
          !! N-Pools should be related to respective C-pools simply by division with approriate C/N-ratio (except for green litter
          !! pool). But use here more complicated update following the logic of the above allocation scheme to allow for checking
          !! of N-conservation!
          !!

          Npool_green          = Npool_green       + redFact_Nlimitation * Nflux_NPP_2_green_pot

          Npool_woods          = Npool_woods       + redFact_Nlimitation * Nflux_NPP_2_wood_pot

          Npool_litter_green_ag   = Npool_litter_green_ag                              &
               - redFact_Nlimitation * (Nflux_litterGreen_ag_2_slow_pot + minNflux_litter_green_ag_pot)

          Npool_litter_green_bg   = Npool_litter_green_bg                              &
               - redFact_Nlimitation * (Nflux_litterGreen_bg_2_slow_pot + minNflux_litter_green_bg_pot)

          Npool_litter_wood_ag = Npool_litter_wood_ag                              &
               - redFact_Nlimitation * (Nflux_litter_wood_ag_2_slow_pot + minNflux_litter_wood_ag_pot)

          Npool_litter_wood_bg = Npool_litter_wood_bg                              &
               - redFact_Nlimitation * (Nflux_litter_wood_bg_2_slow_pot + minNflux_litter_wood_bg_pot)

          Npool_slow           = Npool_slow       &
              + redFact_Nlimitation *  ( Nflux_litterGreen_ag_2_slow_pot + Nflux_litterGreen_bg_2_slow_pot &
                                       + Nflux_litter_wood_ag_2_slow_pot + Nflux_litter_wood_bg_2_slow_pot ) 

          !! ============== UPDATE MINERAL-N-POOL =======================================================================

          SMINN_pool        = max(SMINN_pool - redFact_Nlimitation * total_minN_demand,0.0_dp)

          !!... Actual N2O_emissions_mineraliz() based on the N-availability
          N2O_emissions_mineraliz =  N2O_emissions_slow     &                           !! Net N2O emissions  
                                   + redFact_Nlimitation * N2O_emissions_litter_pot
       end IF

       !! ============== COMPUTE ACTUAL NPP (is needed outside the routine e.g. to recompute transpiration through stomata) =======

       if(NPPrate > 0.0_dp) then 
          NPP_act = (NPP_2_green_via_mobileN + NPP_2_wood_via_mobileN + NPP_2_rootExudates +               &
               redFact_Nlimitation * (Cflx_NPP_2_green_pot + Cflx_NPP_2_wood_pot) + Cflx_NPP_2_reserve_pot &
               )/sec_per_day
       else
          NPP_act = NPPrate + NPP_flux_correction  !dk: NPPrate negatave, NPP_flux_correction could not be removed from reserve pool
       endif


       !! ============== COMPUTE SOIL RESPIRATION ===================================================================

       IF (.NOT. with_yasso) THEN
          soilResp_rate = -(                                                                              &
               redFact_Nlimitation * (  Cflx_litterGreen_ag_2_atmos_pot + Cflx_litterGreen_bg_2_atmos_pot &
                                      + Cflx_litter_wood_ag_2_atmos_pot + Cflx_litter_wood_bg_2_atmos_pot &
                                    ) + Cflx_slow_2_atmos                                                 &
                           ) / sec_per_day

          soilResp_rate_pot = -( (  Cflx_litterGreen_ag_2_atmos_pot + Cflx_litterGreen_bg_2_atmos_pot    &
                                   + Cflx_litter_wood_ag_2_atmos_pot +  Cflx_litter_wood_bg_2_atmos_pot  &
                                  ) + Cflx_slow_2_atmos                                                  &
                        ) / sec_per_day
  
       ELSE 
       
          leafLitter    = Cflx_litterTotal - C_2_litter_woodPool + Cflx_faeces_2_LG - Cflx_2_crop_harvest

          !! call yasso for leaf litter; this time with nutrient effects
          CALL yasso(Yasso_io_pools(1:9), Weather, leafLitter, LeafLit_coefV,               & 
                           sizeZero, Yasso_out, frac_green_aboveGround, NPP_2_rootExudates,redFact_Nlimitation)

          YCpool_acid_ag1         =Yasso_out(1)
          YCpool_water_ag1        =Yasso_out(2)
          YCpool_ethanol_ag1      =Yasso_out(3)
          YCpool_nonsoluble_ag1   =Yasso_out(4)
          YCpool_acid_bg1         =Yasso_out(5)
          YCpool_water_bg1        =Yasso_out(6)
          YCpool_ethanol_bg1      =Yasso_out(7)
          YCpool_nonsoluble_bg1   =Yasso_out(8)
          YCpool_humus_1          =Yasso_out(9)


          !! update respiration with respiration from leaf litter decomposition
          soilResp_rate = Yasso_out(10) / sec_per_day 

          !! call yasso for woody litter; this time with nutrient effects
          CALL yasso(Yasso_io_pools(10:18), Weather, C_2_litter_woodPool, WoodLit_coefV,       & 
                           WoodLitterSize, Yasso_out, frac_wood_aboveGround, 0.0_dp, redFact_Nlimitation) ! no exudates
  
          YCpool_acid_ag2         =Yasso_out(1)
          YCpool_water_ag2        =Yasso_out(2)
          YCpool_ethanol_ag2      =Yasso_out(3)
          YCpool_nonsoluble_ag2   =Yasso_out(4)
          YCpool_acid_bg2         =Yasso_out(5)
          YCpool_water_bg2        =Yasso_out(6)
          YCpool_ethanol_bg2      =Yasso_out(7)
          YCpool_nonsoluble_bg2   =Yasso_out(8)
          YCpool_humus_2          =Yasso_out(9)
       
          !! Yasso total respiration
          ! add respiration from wood litter decomposition
           soilResp_rate = soilResp_rate  + Yasso_out(10) / sec_per_day 

       END IF

       !! ============== Decay of crop harvest =============

       Cflx_crop_harvest_2_atm = Cpool_crop_harvest / tau_Cpool_crop_harvest
       Cpool_crop_harvest = Cpool_crop_harvest - Cflx_crop_harvest_2_atm
       soilResp_rate = soilResp_rate - Cflx_crop_harvest_2_atm / sec_per_day

       IF (with_Nitrogen) THEN
          soilResp_rate_pot = soilResp_rate_pot - Cflx_crop_harvest_2_atm / sec_per_day
          Nflx_crop_harvest_2_SMINN = Npool_crop_harvest / tau_Cpool_crop_harvest
          Npool_crop_harvest = Npool_crop_harvest - Nflx_crop_harvest_2_SMINN
          SMINN_pool = SMINN_pool + Nflx_crop_harvest_2_SMINN
       END IF

       Cflx_2_crop_harvest = Cflx_2_crop_harvest / sec_per_day ! flux from day-1 to s-1
       Cflx_crop_harvest_2_atm = Cflx_crop_harvest_2_atm / sec_per_day ! flux from day-1 to s-1

       !! ============== COMPUTE Herbivory RESPIRATION =====
        
       Cflx_herbivory_2_atm = -(Cflx_faeces_2_atm) /sec_per_day

       !! ============== Litter and herbivory diagnostics (flux from day-1 to s-1)
       Cflx_litterTotal = Cflx_litterTotal /sec_per_day
       litter_leaf = litter_leaf /sec_per_day ! HW
       Cflx_herbivory    = Green_2_herbivory / sec_per_day
       Cflx_herbivory_LG = Cflx_faeces_2_LG / sec_per_day   ! DSG: goes into yasso


       IF (with_Nitrogen) THEN

          redFact_Nlimit = redFact_Nlimitation

          !! ============== OTHER DIAGNOSTIC OUTPUTS ===================================================================

         !! redFact_Nlimit was missing first, but minNflux is reduced with N limitation
          minNflux_litter_green_ag   = redFact_Nlimitation * minNflux_litter_green_ag_pot / sec_per_day
          minNflux_litter_green_bg   = redFact_Nlimitation * minNflux_litter_green_bg_pot / sec_per_day
          minNflux_litter_wood_ag    = redFact_Nlimitation * minNflux_litter_wood_ag_pot  / sec_per_day
          minNflux_litter_wood_bg    = redFact_Nlimitation * minNflux_litter_wood_bg_pot  / sec_per_day
          minNflux_slow              = minNflux_slow / sec_per_day                                      ! no N limited
          Nplant_demand = Nplant_demand / sec_per_day                                                   ! demands already N limited
          Nsoil_demand  = Nsoil_demand  / sec_per_day
          Ntotal_demand = Ntotal_demand / sec_per_day
          N2O_emissions_mineraliz =  N2O_emissions_mineraliz / sec_per_day                              ! already N limited 
                                                                                              ! only N2O_emissions_litter_pot	   
          N2O_emissions_slow =  N2O_emissions_slow / sec_per_day                  ! N2O lost from slow to atmosphere; not N limited
          N2O_emissions_grazing = N2O_emissions_grazing / sec_per_day

       end IF

       ! fraction of new wood litter, assume old and new litter are respired in the same way, respiration does not change the fraction
       if (present(frac_litter_wood_new)) then
          if (with_yasso) then
            litter_wood_total = YCpool_water_ag2 + YCpool_acid_ag2 + YCpool_ethanol_ag2 + YCpool_nonsoluble_ag2
          else
            litter_wood_total = Cpool_litter_wood_ag
          endif
          if ((litter_wood_total >  EPSILON(1._dp)) .AND. (litter_wood_new < litter_wood_total)) then
             frac_litter_wood_new = litter_wood_new / litter_wood_total
          else
             frac_litter_wood_new = 1._dp    ! frac_litter_wood_new=1 means fuel fractions will be reset to original values
          endif
       endif
 
    end if

  end subroutine update_Cpools


  ! --- N_process() ---------------------------------------------------------------------------------------------------  
  !!
  !! ====Biological N Fixation, Deposition, Denitrification(N-forms gases loss) & Leaching===========
  !DESCRIPTION: On the radiation time step, update the nitrogen fixation rate as a function of annual total NPP (60PgC/Y) & ET. 
  !This rate gets updated once per year.! All N fixation goes to the soil mineral N pool. 
  !Forcing Nfix to global constant = 0.4 gN/m2/yr; BNF= 105 TgN/y; Range by Cleveland, 1999 is 100-290 TgN/y. 
  !BNF_NPP = 0.001666*NPP_annual & BNF_ET = 0.008 * (ET -40); unit: 40 cm/y or 400/365*86400 mm/s

  !BNF Calibration: For 100 Tg N/Y; -----------------------------------
  !                                     100 Tg   100 * e12
  !                             factor  = ------ = ---------- = 0.00166
  !                                     60 Pg    60  * e15 
  !converts to BNF in mole ->12(14); calibrated as 12 (to gN)->100Tg N/y =0.00166 -> more in NCAR CLM model codes 
  !->BNF_NPP = 1.8 *(1.0- exp(-0.003 * NPP_annual))/12 

  !                   OR
  ! BNF_NPP = 0.005*NPP (CLM-if NPP is 1PgC); BNF_NPP = 0.0005*NPP (TRIPLEX)
  ! BNF_NPP = 0.005*NPP/365   per day

  !!..Fixation based on ET -> Century Model, Schimel, 1996; Xu-Ri, 2008 
  !->BNF_ET=0.008 * (ET - 400)/(365 * 86400) 

  !....N Deposition (wet + Dry)----------
  ! Global const = 0.5 g/m2/y = 1.5*10-8 g/m2/s (Ref -> Class model, 2007 & Dickison, 2002 
  !->Ndepo= (0.4/(365* 86400))/14  !! g - mole Nitrogen 
  !Reference -> Frank Dentener, 2006 & Galloway et al., 2004 // Inputs for 1860, 1990, 2000, 2030, 2050 
  !Data downloaded from ftp://daac.ornl.gov/data/global_climate/global_N_deposition_maps
  !Global Maps of Atmospheric Nitrogen Deposition, 1860, 1993, and 2050 -> TM3 model: 3-D Chemistry trasport 
  !->Data used: Dentener, 2006 GBC V20 --1860, 2000(S1), 2030(S4-SRES A2) (by Personal communicattion)

  !....Denitrification (gaseous loss of Nitrogen: NO, N2O, N2)--- 
  !!-> Denit = 1.15*SMINN**0.57   -> Ref: Century model, Thomas,2002 per day
  !!-> Denit = dnp*SMINN          -> Ref: CLM model dnp = 0.01 per day
  
  !---Leaching descriptions-----------------
  !! L = (sf*SMINN)*runoff/tot_liquid_water -> Ref:CLM model and also in Century (Thomas, 2002); sf = 0.1
  !! L = f_leach * (SMINN_pool * 0.1_dp)    -> Ref:Xu Ri, 2008 
  !! unit conversion -> 1m water = 1000 kg water & 1mm water = 1 kg water 
  !!->CLM model
  !! L_online = (SMINN_pool * sf) * runoff_acc/moisture  !!.. multiplied with 0.001 m->kg H2O--(63 code)
  !! L_offline = SMINN * sf * R/M                        !!-> maxmoist in m --runoff from Echam output 160 code
  !! L = SMINN * sf * R/(alpha*max_moisture)
  !! The above formulation results in an extremely strong n limitation in grid boxes in high northern lattitudes that 
  !! feature shallow soils.To mitigate this problem leaching is formulated as a function of the drainage exclusivly. 
  
  ! --- Nitros oxide emissions(N2O)--N2O_Emissions() ----
  !! ====N2O Emissions from Natural Ecosystems Soils & Agricultural Ecosystems Soils
  !! N2O_emissions = 1% of N Deposition per day -->> IPCC methodology; Skiba et al, 1998; Bloemerts, M., 2009.
  !! & from BNF, the emission factor is 1.25% in IPCC methodology
  !! BNF converts N2 to NH4+ which undergone nitrification & denitrification process

  elemental pure subroutine N_process (NPP_act_yDayMean, nfix_to_sminn, ndep_to_sminn,                &
                                       ndep_forc, nfert_forc, nfert_to_sminn, sminn_to_denit, sminn_leach, &
                                       SMINN_pool, NetEcosyst_N_flux, drainage, alpha,   &
                                       max_moisture, NPP_run_mean, N2O_emissions_depfix, is_vegetation, &
                                       N2O_emissions_nfert) 

    real(dp),intent(in)    :: max_moisture              !! Depth of soil bucket in [m]
    real(dp),intent(out)   :: nfix_to_sminn             !! SMINN increase from biological nitrogen fixation [mol(N)/m^2(canopy)s]
                                                        !!    based on NPP
    real(dp),intent(out)   :: ndep_to_sminn             !! SMINN increase from N-deposition [mol(N)/m^2(canopy) s]
    real(dp),intent(in)    :: ndep_forc                 !! Deposition of N [Mol(N)/m^2(canopy) s] based on atm forcing.
    real(dp),intent(in)    :: nfert_forc                !! Fertilizer of N [g(N)/m^2(grid box) year] based on external forcing.
    real(dp),intent(out)   :: nfert_to_sminn            !! SMINN increase from N-fertilizer [mol(N)/m^2(canopy) s] only in crops   
    real(dp),intent(out)   :: sminn_to_denit            !! SMINN loss by dentrification [mol(N)/m^2(canopy)s]
    real(dp),intent(out)   :: sminn_leach               !! SMINN loss by leaching [mol(N)/m^2(canopy)s] 
    real(dp),intent(in)    :: NPP_act_yDayMean          !! Mean NPP at the current day (this routine is called once a day) 
                                                        !!    [mol(C)/m^2 s]
    real(dp),intent(inout) :: SMINN_pool                !! mineral N in soils [mol(N)/m^2(canopy)]
    real(dp),intent(in)    :: drainage                  !! Yesterdays mean drainage [m/s] 
    real(dp),intent(in)    :: alpha                     !! Yesterdays mean water stress factor (between 0 and 1)
    real(dp),intent(out)   :: NetEcosyst_N_flux         !! Total balance of N gains and losses (positive for ecosystem gain) 
                                                        !!    [mol(N)/m^2(canopy)s]
    real(dp),intent(inout) :: NPP_run_mean              !! exponential running mean of NPP used for nitrogen fixation estimates
    real(dp),intent(out)   :: N2O_emissions_depfix      !! Nitros oxide emission [mol(N)/m^2(canopy)s]
    logical, intent(in)    :: is_vegetation             !!
    real(dp),intent(out)   :: N2O_emissions_nfert       !! Nitros oxide emission [mol(N)/m^2(canopy)s] from N fertilizer use    

    !! parameters

    REAL(dp),parameter  :: delta_time = 1.0_dp            !! time step (=calling interval) of this routine [days]
    REAL(dp),parameter  :: decay_time = 365.0_dp          !! Decay time of memory loss for exponential running mean of NPP [days]

    ! awi: this is the standard value with cbalance
    !REAL(dp),parameter  :: sf  = 0.0005_dp                !! fraction of soil mineral N in soluble form < 0.05% !awi: original
    !awi: the following values for the sf parameter are used with yasso:
    !REAL(dp),parameter  :: sf  = 0.05_dp                 !! fraction of soil mineral N in soluble form < 5%
    !REAL(dp),parameter  :: sf  = 0.08_dp                  !! fraction of soil mineral N in soluble form < 8%
    !REAL(dp),parameter  :: sf  = 0.12_dp                  !! fraction of soil mineral N in soluble form < 12%
    !REAL(dp),parameter  :: sf  = 0.16_dp                  !! fraction of soil mineral N in soluble form < 16%
    REAL(dp),parameter  :: sf  = 0.2_dp                  !! fraction of soil mineral N in soluble form < 20%
    !awi: the following values for the dnp parameter are used with yasso:
    !REAL(dp),parameter  :: dnp = 0.0001_dp                !! denitrification rate per day (< .01%)
    !REAL(dp),parameter  :: dnp = 0.00015_dp               !! denitrification rate per day (< .015%) !
    !REAL(dp),parameter  :: dnp = 0.0002_dp                !! denitrification rate per day (< .02%) ! awi: this is also the standard
                                                          !!   value with cbalance
    REAL(dp),parameter  :: dnp = 0.000275_dp                !! denitrification rate per day (< .0275%)
    REAL(dp),parameter  :: N2O_ef_Ndepo = 0.01_dp         !! N2O emissions Factor (.01); 1% in IPCC methodology
    REAL(dp),parameter  :: N2O_ef_Nfix  = 0.0125_dp       !! N2O emissions Factor (.0125); 1.25% in IPCC methodology
    REAL(dp),parameter  :: N2O_ef_Nfert = 0.0125_dp       !! N2O emissions Factor (.025); 2.5% in Davidson E.A (2009)

    !! locals

    REAL(dp) :: F_NPP_act               !! Exponential weighting factor for computing NPP_run_mean [1]
    REAL(dp) :: N_NPP_act               !! Normalization for computing NPP_run_mean [1] 
    REAL(dp) :: f_leach                 !! Fraction of leach

    nfix_to_sminn           = 0.0_dp
    ndep_to_sminn           = 0.0_dp
    sminn_to_denit          = 0.0_dp
    sminn_leach             = 0.0_dp
    NetEcosyst_N_flux       = 0.0_dp
    N2O_emissions_depfix    = 0.0_dp
    nfert_to_sminn          = 0.0_dp
    N2O_emissions_nfert     = 0.0_dp    

    !! -- Computation of an exponentially weighted running mean of NPP
    !!
    !!                        NPP(t+1)         delta_T
    !!    NPP_runmean(t+1) =  -------- + exp(- ------- ) * NPP_runmean(t)
    !!                          Norm             tau
    !!
    !! where the normalization Norm is given by
    !!
    !!                                delta_T
    !!    Norm = SUM(k=0,infty) exp(- ------- k) = (1 - exp(-delta_T/tau))^(-1)
    !!                                  tau

    F_NPP_act = exp(-delta_time/decay_time)
    N_NPP_act = 1._dp / (1._dp - F_NPP_act) 
    NPP_run_mean = NPP_act_yDayMean/N_NPP_act + F_NPP_act * NPP_run_mean

    !! --- Biological nitrogen fixation  based on NPP 

    !! awi: 0.7 reduced to 0.53 for YASSO-CN version 
    nfix_to_sminn = &                               !! 0.7 is calibrated to match global n2fix (present day):e.g  see Cleveland 1999;
         (0.53_dp * (1.0_dp - exp(-0.003_dp * NPP_run_mean * sec_per_year)))/sec_per_year   &
        *(molarMassN_kg/molarMassC_kg)              !! Nitrogen fixation in [mol(N)/m^2 s]
    nfix_to_sminn = max(0.0_dp,nfix_to_sminn)       !! To prevent n-fixation for negative NPP 

    !! -- N Deposition (wet + Dry) ---------- 

    IF (is_vegetation) ndep_to_sminn = ndep_forc

    !! -- Denitrification--------
                      
    sminn_to_denit = SMINN_pool * dnp/sec_per_day * alpha  !!dk: dnp depends on T and H2O; future plan: scale with N2O

    !! -- Nitros oxide emissions(N2O)-------

    N2O_emissions_depfix = ndep_to_sminn * N2O_ef_Ndepo           !! N2O emissions due to N deposition

    !! -- Leaching ----------------

    f_leach = min(1.0_dp, drainage*sec_per_day/max(1e-13_dp, alpha*max_moisture)) !! fraction of water removed during one day 
                                                                                !! (leaching limited to total bucket content)
    sminn_leach = f_leach * SMINN_pool * sf/sec_per_day !! Mol(N) leached per second (only soluble fraction  can be leached)
    !!-- Distribute N fertiliser into all crops types; note: N-fertilizer read in [g (N)/m^2(grid box) year]
    nfert_to_sminn = (nfert_forc * 0.5_dp/sec_per_year)*1000._dp/molarMassN_kg     !!.. gm/14 = conversion fact. g(N) to mol(N)
                                                                       !!..& half lost as harvested N
    N2O_emissions_nfert = (nfert_forc /sec_per_year)*1000._dp/molarMassN_kg * N2O_ef_Nfert   !! 2.5% of total Nfert use in crops  


    !--------------update SMINN_pool----

    SMINN_pool = SMINN_pool + nfix_to_sminn * sec_per_day           !! N-Fixation        
    SMINN_pool = SMINN_pool + ndep_to_sminn * sec_per_day           !! N-Deposition
    SMINN_pool = SMINN_pool + nfert_to_sminn * sec_per_day          !! N-Fertilizer application to crops
    SMINN_pool = SMINN_pool - sminn_to_denit* sec_per_day           !! N-Denitrification
    SMINN_pool = SMINN_pool - sminn_leach   * sec_per_day           !! N-Leaching
    SMINN_pool = SMINN_pool - N2O_emissions_depfix * sec_per_day    !! update SMINN_pool
    SMINN_pool = SMINN_pool - N2O_emissions_nfert * sec_per_day     !! update SMINN_pool


 !--------------Assessment of N inputs & lossess--- NetEcosyst_N_flux
 
    NetEcosyst_N_flux = nfix_to_sminn + ndep_to_sminn + nfert_to_sminn - (sminn_to_denit + sminn_leach) &
                       - N2O_emissions_depfix  - N2O_emissions_nfert 

  end subroutine N_process


  ! --- relocate_carbonAndNitrogen() -----------------------------------------------------------------------------------------------
  !
  ! This routine models the consequences of shifts in the cover fractions of the tiles for the carbon (C) and nitrogen (N) pools.
  ! This change in vegetation composition is either due to human land-cover change or due to vegetation dynamics.
  ! Concerning the C and N stored on those partitions of shrinking tiles that are lost to extending tiles during the time interval
  ! considered: In the case of landcover change the C and N from living plant pools (Cpool_green, Cpool_reserve and Cpool_woods plus
  ! associated N-pools) is relocated, partly by a release to the atmosphere (e.g. by fire stubbing) and the remaining C and N into
  ! soil pools.
  ! The subroutine is structured in 4 parts:
  ! 1. determine amount of C and N released from living plant pools to atmosphere, litter, and soil (in case of land cover change)
  !    or rescale C and N pools of living plants on shrinking tiles (in case of dynamic vegetation)
  ! 2. determine amount of C and N to be relocated from slow pool and litter pools of shrinking to extending tiles (as a result of
  !    the change in area)
  ! 3. lower down the C and N density for extending tiles (to account for the change in area)
  ! 4. distribute freed C and N (from part 1 and part 2 of the routine) to extending tiles
  ! ATTENTION:
  ! In the case of vegetation dynamics it is assumed that the C and N from living plant pools is already lost to other
  ! reservoirs (atmosphere, litter, soil) by those processes that removed the vegetation (general carbon flow reflecting minimum
  ! mortality, fire, wind break) -- therefore in the case of vegetation dynamics this routine has to be called in combination with
  ! other routines that handle those destructive processes in order to assure consistency.
  ! ATTENTION: 
  ! 1. If this routine is called with N-pools, then not only carbon is relocated, but also nitrogen.
  ! 2. To handle C and N relocation for landcover change, the necessary arrays must be present in the call.
  !    In absence of N-pools ("carbon only mode") only carbon is reshuffled.  
  !

  subroutine relocate_CarbonAndNitrogen(lctlib, surface, cf_old, cf_new, veg_fract_correction,      &
                                        ! Carbon
                                        Cpool_green, Cpool_woods, Cpool_reserve,                    &
                                        Cpool_crop_harvest,                                         &
                                        ! soil carbon pools of the old cbalance scheme
                                        Cpool_litter_green_ag, Cpool_litter_green_bg,               &
                                        Cpool_litter_wood_ag, Cpool_litter_wood_bg,                 &
                                        Cpool_slow,                                                 &
                                        ! Yasso litter and soil carbon pools
                                        YCpool_acid_ag1,                                             &
                                        YCpool_water_ag1,                                            &
                                        YCpool_ethanol_ag1,                                          &
                                        YCpool_nonsoluble_ag1,                                       &
                                        YCpool_acid_bg1,                                             &
                                        YCpool_water_bg1,                                            &
                                        YCpool_ethanol_bg1,                                          &
                                        YCpool_nonsoluble_bg1,                                       &
                                        YCpool_humus_1,                                              &
                                        YCpool_acid_ag2,                                             &
                                        YCpool_water_ag2,                                            &
                                        YCpool_ethanol_ag2,                                          &
                                        YCpool_nonsoluble_ag2,                                       &
                                        YCpool_acid_bg2,                                             &
                                        YCpool_water_bg2,                                            &
                                        YCpool_ethanol_bg2,                                          &
                                        YCpool_nonsoluble_bg2,                                       &
                                        YCpool_humus_2,                                              &
                                        ! YASSO coefficients, only needed with LCC
                                        LeafLit_coef,                                               &
                                        WoodLit_coef,                                               &
                                        ! LCC: lcc_scheme=1 (standard JSBACH) or lcc_scheme=2 (grand slam protocol)
                                        C_2_atmos,                                                  &
                                        C_2_litterGreenPools,                                       &
                                        C_2_litterWoodPool_ag,C_2_litterWoodPool_bg,                &
                                        ! Nitrogen
                                        Npool_green, Npool_woods, Npool_mobile,                     &
                                        Npool_crop_harvest,                                         &
                                        Npool_litter_green_ag, Npool_litter_green_bg,               &
                                        Npool_litter_wood_ag, Npool_litter_wood_bg,                 &
                                        Npool_slow,                                                 &
                                        SMINN_pool,                                                 &
                                        ! nitrogen variables for LCC
                                        Nitrogen_2_atmos, Nitrogen_2_litterGreenPools,              &
                                        Nitrogen_2_litterWoodPool_ag, Nitrogen_2_litterWoodPool_bg, &
                                        Nitrogen_2_SMINNpool,                                       &
                                        ! LCC scheme 2: fluxes after Houghton
                                        Cpool_onSite,Cpool_paper,Cpool_construction,                &
                                        C_onSite_2_atmos, C_paper_2_atmos, C_construction_2_atmos,  &
                                        C_2_onSite, C_2_paper, C_2_construction,                    &
                                        lcc_scheme)
    
    USE mo_exception,        ONLY: finish
    USE mo_cbal_parameters,  ONLY: frac_wood_2_atmos, frac_green_2_atmos, frac_reserve_2_atmos, frac_mobile_2_atmos
    USE mo_jsbach_lctlib,    ONLY: lctlib_type 
    USE mo_land_surface,     ONLY: land_surface_type

    type(lctlib_type),intent(in) :: lctlib               !! PFT-specific constants
    type(land_surface_type),intent(in) :: surface        !! Access to cover_type
    real(dp),intent(in)    :: cf_old(:,:)                !! Cover fractions before landcover change or vegetation dynamics [] 
    real(dp),intent(in)    :: cf_new(:,:)                !! Cover fraction after landcover change or vegetation dynamics []
    real(dp),intent(in)    :: veg_fract_correction(:,:)  !! Correction factor for cover fractions 1-exp(-LAI_max/2) (accounts for
                                                         !!    sparseness of vegetation)
    real(dp),intent(inout) :: Cpool_green(:,:)           !! Value of green carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_woods(:,:)           !! Value of wood carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_reserve(:,:)         !! Value of reserve carbon pool [mol(C)/m^(canopy)2]
    real(dp),intent(inout) :: Cpool_crop_harvest(:,:)    !! Value of crop harvest carbon pool [mol(C)/m^(canopy)2]
    
    ! Cbalance litter & soil pools:
    real(dp), optional, intent(inout) :: Cpool_litter_green_ag(:,:) !! above ground green litter carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: Cpool_litter_green_bg(:,:) !! below ground green litter carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: Cpool_litter_wood_ag(:,:)  !! above ground wood litter carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: Cpool_litter_wood_bg(:,:)  !! below ground wood litter carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: Cpool_slow(:,:)            !! slow soil carbon pool [mol(C)/m^2(canopy)]
    
    ! Yasso litter & soil pools:
    ! Size class 1; green litter
    real(dp), optional, intent(inout) :: YCpool_acid_ag1(:,:)       !! above ground acid soluble carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_water_ag1(:,:)      !! above ground water soluble carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_ethanol_ag1(:,:)    !! above ground ethanol soluble carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_nonsoluble_ag1(:,:) !! above ground non carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_acid_bg1(:,:)       !! below ground acid soluble carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_water_bg1(:,:)      !! below ground water soluble carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_ethanol_bg1(:,:)    !! below ground ethanol soluble carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_nonsoluble_bg1(:,:) !! below ground non soluble carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_humus_1(:,:)         !! humus carbon pool [mol(C)/m^(canopy)2]
    ! Size class 2; woody litter
    real(dp), optional, intent(inout) :: YCpool_acid_ag2(:,:)       !! above ground acid soluble carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_water_ag2(:,:)      !! above ground water soluble carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_ethanol_ag2(:,:)    !! above ground ethanol soluble carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_nonsoluble_ag2(:,:) !! above ground non carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_acid_bg2(:,:)       !! below ground acid soluble carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_water_bg2(:,:)      !! below ground water soluble carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_ethanol_bg2(:,:)    !! below ground ethanol soluble carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_nonsoluble_bg2(:,:) !! below ground non soluble carbon pool [mol(C)/m^(canopy)2]
    real(dp), optional, intent(inout) :: YCpool_humus_2(:,:)         !! humus carbon pool [mol(C)/m^(canopy)2]

    ! Variables needed with standard jsbach LCC scheme
    real(dp),optional,intent(out) :: C_2_atmos(:)                 !! Amount of carbon directly emitted to atmosphere in this
                                                                  !!    timestep [mol(C)/m^2(vegetated area)]
    real(dp),optional,intent(out) :: C_2_litterGreenPools(:)      !! Amount of carbon relocated by land cover change from green
                                                                  !!    and reserve pool to below and above ground green litter
    real(dp),optional,intent(out) :: C_2_litterWoodPool_ag(:)     !! Amount of carbon relocated by land cover change from wood
                                                                  !!    pool to above ground woody litter pool
                                                                  !!    [mol(C)/m^2(vegetated area)]
    real(dp),optional,intent(out) :: C_2_litterWoodPool_bg(:)     !! Amount of carbon relocated by land cover change from wood
                                                                  !!    pool to below ground woody litter pool
                                                                  !!    [mol(C)/m^2(vegetated area)]
    real(dp),optional,intent(inout) ::  Cpool_onSite(:)           !! Amount of carbon remains in onSite anthro annual pool 
                                                                  !!  [mol(C)/m^2(vegetated area)] 
    real(dp),optional,intent(inout) ::  Cpool_paper(:)            !! Amount of carbon remains in anthro decadal pool
    real(dp),optional,intent(inout) ::  Cpool_construction(:)     !! Amount of carbon remains in anthro centinnial pool
    real(dp),optional,intent(out)   ::  C_onSite_2_atmos(:)       !! Carbon flux from anthro annual pool to atmos 
                                                                  !!  [mol(C)/m^2(vegetated area)/day] 
    real(dp),optional,intent(out)   ::  C_paper_2_atmos(:)              !! Carbon flux from green anthro annual pool to atmos 
                                                                        !!   [mol(C)/m^2(vegetated area)/day] 
    real(dp),optional,intent(out)   ::  C_construction_2_atmos(:)       !! Carbon flux from woody anthro annual pool to atmos
                                                                        !!   [mol(C)/m^2(vegetated area)/day]
    real(dp),optional,intent(out)   ::  C_2_onSite(:)                   !! Carbon flux to onsite [mol(C)/m^2(vegetated area)/day]
    real(dp),optional,intent(out)   ::  C_2_construction(:)             !! Carbon flux from woody anthro annual pool to paper
                                                                        !!   [mol(C)/m^2(vegetated area)/day] 
    real(dp),optional,intent(out)   ::  C_2_paper(:)                    !! Carbon flux from woody anthro annual pool to
                                                                        !!   construction [mol(C)/m^2(vegetated area)/day]

    real(dp),intent(inout),optional :: Npool_green(:,:)           !! Green nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_woods(:,:)           !! Wood nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_mobile(:,:)          !! Plant mobile nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_crop_harvest(:,:)    !! Crop harvest nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_litter_green_ag(:,:) !! Above ground green litter nitrogen pool [mol(N)/m^(canopy)2]
    real(dp),intent(inout),optional :: Npool_litter_green_bg(:,:) !! Below ground green litter nitrogen pool [mol(N)/m^(canopy)2]
    real(dp),intent(inout),optional :: Npool_litter_wood_ag(:,:)  !! Wood litter nitrogen pool [mol(N)/m^(canopy)2]
    real(dp),intent(inout),optional :: Npool_litter_wood_bg(:,:)  !! Wood litter nitrogen pool [mol(N)/m^(canopy)2]
    real(dp),intent(inout),optional :: Npool_slow(:,:)            !! Slow soil nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: SMINN_pool(:,:)            !! Soil mineral nitrogen pool [mol(N)/m^2(canopy)]

    real(dp),optional,intent(out) :: Nitrogen_2_atmos(:)             !! Amount of nitrogen directly emitted to atmosphere by lcc 
                                                                     !! .. in this timestep [mol(N)/m^2(vegetated area)]
    real(dp),optional,intent(out) :: Nitrogen_2_litterGreenPools(:)  !! Amount of nitrogen relocated by land cover change from green
                                                                     !!    and reserve pool to above and below ground green litter
                                                                     !!    pools [mol(N)/m^2(vegetated area)]
    real(dp),optional,intent(out) :: Nitrogen_2_litterWoodPool_ag(:) !! Amount of nitrogen relocated by land cover change from wood
                                                                     !!    pool to above ground woody litter pool 
                                                                     !!    [mol(N)/m^2(vegetated area)]
    real(dp),optional,intent(out) :: Nitrogen_2_litterWoodPool_bg(:) !! Amount of nitrogen relocated by land cover change from wood
                                                                     !!    pool to below ground woody litter pool
                                                                     !!    [mol(N)/m^2(vegetated area)]
    real(dp),optional,intent(out) :: Nitrogen_2_SMINNpool(:)         !! Amount of nitrogen relocated by land cover change from wood
                                                                     !!    pool to soil mineral N pool (this is the surplus N that 
                                                                     !!    the woody litter N-pool cannot take up)   

    REAL(dp),intent(in),optional   :: LeafLit_coef(:,:,:)            !! fractions to seperate fresh litter to non woody yasso pools
    REAL(dp),intent(in),optional   :: WoodLit_coef(:,:,:)            !! fractions to seperate fresh litter to woody yasso pools
    
    integer,optional, intent(in) :: lcc_scheme
    
    !! locals
    !! yasso variables
    real(dp) :: Chumus_1_2_humus_1(1:size(cf_old,DIM=1))                 !! Total carbon relocated from humus pools of shrinking 
                                                                     !! ..  to humus pools of extending tiles
    real(dp) :: Cacid_ag1_2_acid_ag1(1:size(cf_old,DIM=1))             !! Total carbon relocated from above ground acid pools of
                                                                     !! ..  shrinking to above ground acid pools of extending tiles
    real(dp) :: Cacid_bg1_2_acid_bg1(1:size(cf_old,DIM=1))             !! Total carbon relocated from below ground acid pools of
                                                                     !! ..  shrinking to below ground acid pools of extending tiles
    real(dp) :: Cwater_ag1_2_water_ag1(1:size(cf_old,DIM=1))           !! Total carbon relocated from above ground water pools of
                                                                     !! ..  shrinking to above ground water pools of extending tiles
    real(dp) :: Cwater_bg1_2_water_bg1(1:size(cf_old,DIM=1))           !! Total carbon relocated from below ground water pools of
                                                                     !! ..  shrinking to below ground water pools of extending tiles
    real(dp) :: Cethanol_ag1_2_ethanol_ag1(1:size(cf_old,DIM=1))       !! Total carbon relocated from above ground ethanol pools of
                                                                     !! ..  shrinking to above ground ethanol pools of ext. tiles
    real(dp) :: Cethanol_bg1_2_ethanol_bg1(1:size(cf_old,DIM=1))       !! Total carbon relocated from below ground ethanol pools of
                                                                     !! ..  shrinking to below ground ethanol pools of ext. tiles
    real(dp) :: Cnonsoluble_ag1_2_nonsoluble_ag1(1:size(cf_old,DIM=1)) !! Total carbon relocated from above ground nonsoluble pools
                                                                       !! of shrinking to above ground nonsoluble pools of ext. tiles
    real(dp) :: Cnonsoluble_bg1_2_nonsoluble_bg1(1:size(cf_old,DIM=1)) !! Total carbon relocated from below ground nonsoluble pools
                                                                       !! of shrinking to below ground nonsoluble pools of ext. tiles
    real(dp) :: Chumus_2_2_humus_2(1:size(cf_old,DIM=1))                 !! Total carbon relocated from humus pools of shrinking 
                                                                     !! ..  to humus pools of extending tiles
    real(dp) :: Cacid_ag2_2_acid_ag2(1:size(cf_old,DIM=1))             !! Total carbon relocated from above ground acid pools of
                                                                     !! ..  shrinking to above ground acid pools of extending tiles
    real(dp) :: Cacid_bg2_2_acid_bg2(1:size(cf_old,DIM=1))             !! Total carbon relocated from below ground acid pools of
                                                                     !! ..  shrinking to below ground acid pools of extending tiles
    real(dp) :: Cwater_ag2_2_water_ag2(1:size(cf_old,DIM=1))           !! Total carbon relocated from above ground water pools of
                                                                     !! ..  shrinking to above ground water pools of extending tiles
    real(dp) :: Cwater_bg2_2_water_bg2(1:size(cf_old,DIM=1))           !! Total carbon relocated from below ground water pools of
                                                                     !! ..  shrinking to below ground water pools of extending tiles
    real(dp) :: Cethanol_ag2_2_ethanol_ag2(1:size(cf_old,DIM=1))       !! Total carbon relocated from above ground ethanol pools of
                                                                     !! ..  shrinking to above ground ethanol pools of ext. tiles
    real(dp) :: Cethanol_bg2_2_ethanol_bg2(1:size(cf_old,DIM=1))       !! Total carbon relocated from below ground ethanol pools of
                                                                     !! ..  shrinking to below ground ethanol pools of ext. tiles
    real(dp) :: Cnonsoluble_ag2_2_nonsoluble_ag2(1:size(cf_old,DIM=1)) !! Total carbon relocated from above ground nonsoluble pools
                                                                       !! of shrinking to above ground nonsoluble pools of ext. tiles
    real(dp) :: Cnonsoluble_bg2_2_nonsoluble_bg2(1:size(cf_old,DIM=1)) !! Total carbon relocated from below ground nonsoluble pools
                                                                       !! of shrinking to below ground nonsoluble pools of ext. tiles
    real(dp) :: C_2_litterGreenPool_ag(1:size(cf_old,DIM=1))    !! Carbon relocated by land cover change from green and reserve
                                                                !!    pool to above ground green litter
    real(dp) :: C_2_litterGreenPool_bg(1:size(cf_old,DIM=1))    !! Carbon relocated by land cover change from green and reserve
                                                                !!    pool to below ground green litter
    real(dp) :: Nitrogen_2_litterGreenPool_ag(1:size(cf_old,DIM=1))  !! Nitrogen relocated by land cover change from green and 
                                                                     !!    reserve pool to above ground green litter
    real(dp) :: Nitrogen_2_litterGreenPool_bg(1:size(cf_old,DIM=1))  !! Nitrogen relocated by land cover change from green and
                                                                     !!    reserve pool to below ground green litter
    real(dp) :: Cslow_2_slow(1:size(cf_old,DIM=1))                   !! Total carbon relocated from slow pools of shrinking 
                                                                     !! ..  to slow pools of extending tiles
    real(dp) :: Nslow_2_slow(1:size(cf_old,DIM=1))                   !! Total nitrogen relocated from slow pools of shrinking 
                                                                     !! ..  to slow pools of extending tiles
    real(dp) :: Ccrop_harvest_2_crop_harvest(1:size(cf_old,DIM=1))   !! Total carbon relocated from crop harvest pools of shrinking 
                                                                     !! ..  to crop harvest pools of extending tiles
    real(dp) :: Ncrop_harvest_2_crop_harvest(1:size(cf_old,DIM=1))   !! Total nitrogen relocated from crop harvest pools of shrinking 
                                                                     !! ..  to crop harvest pools of extending tiles
    real(dp) :: ClitGreen_ag_2_litGreen_ag(1:size(cf_old,DIM=1))  !! Total carbon relocated from above ground green litter pools of 
                                                                  !! shrinking to above ground green litter pool of extending tiles
    real(dp) :: ClitGreen_bg_2_litGreen_bg(1:size(cf_old,DIM=1))  !! Total carbon relocated from below ground green litter pools of 
                                                                  !! shrinking to below ground green litter pool of extending tiles
    
    real(dp) :: NlitGreen_ag_2_litGreen_ag(1:size(cf_old,DIM=1))  !! Total N relocated from above ground green litter pools of
                                                                  !! shrinking to above ground green litter pools of extending tiles
    real(dp) :: NlitGreen_bg_2_litGreen_bg(1:size(cf_old,DIM=1))  !! Total N relocated from below ground green litter pools of
                                                                  !! shrinking to below ground green litter pools of extending tiles
    real(dp) :: ClitterWood_ag_2_litterWood_ag(1:size(cf_old,DIM=1))    !! Total C relocated from ag wood litter pools of shrinking
                                                                        !! ..  to ag wood litter pools of extending tiles
    real(dp) :: ClitterWood_bg_2_litterWood_bg(1:size(cf_old,DIM=1))    !! Total C relocated from bg wood litter pools of shrinking
                                                                        !! ..  to bg wood litter pools of extending tiles
    real(dp) :: NlitterWood_ag_2_litterWood_ag(1:size(cf_old,DIM=1))    !! Total N relocated from ag wood litter pools of shrinking
                                                                        !! ..  to ag wood litter pools of extending tiles
    real(dp) :: NlitterWood_bg_2_litterWood_bg(1:size(cf_old,DIM=1))    !! Total N relocated from bg wood litter pools of shrinking
                                                                        !! ..  to bg wood litter pools of extending tiles
    real(dp) :: SMINN_2_SMINN(1:size(cf_old,DIM=1))                     !! Total nitrogen relocated from soil mineral N pools of
                                                                        !! ..  shrinking to SMINN pools of extending tile

    logical  :: tiles_extend(1:size(cf_old,DIM=1),1:size(cf_old,DIM=2)) !! Logical mask indicating those tiles that are extending
                                                                        !! by landcover change
    real(dp) :: cf_rescale(1:size(cf_old,DIM=1),1:size(cf_old,DIM=2))   !! Scale factor depending on changing landcover
    real(dp) :: shrinkRatio (1:size(cf_old,DIM=1),1:size(cf_old,DIM=2)) !! Ratio old and new cover fractions of extending tiles
    real(dp) :: flx_2_litter(1:size(cf_old,DIM=1))                      !! flux to litter pools (auxiliary variable, may be C or N)
    real(dp) :: sum_cf_new_old(1:size(cf_old,DIM=1))                    !! Sum of cover fractions added to old cover fractions of
                                                                        !! extending tiles
    logical  :: land_cover_change                                       !! true:  C relocation due to land cover change
                                                                        !! false: C relocation due to (natural) vegetation dynamics 
    integer  :: ntiles, i, nidx
    logical  :: nitrogenMode !! Flag: whether nitrogen shall be handled: only if N-pools are present in the call
    logical  :: with_yasso   !! Flag: whether yasso litter and soil pools shall be handled

    !! PREPARATIONS !!
    land_cover_change = .false.
    NitrogenMode      = .false.
    with_yasso        = .false.

    !! Avoid wrong use of the routine!
    !! This routine handles either land cover change or vegetation dynamics:
    !! For land cover change pass frac_wood_2_atmos, frac_green_2_atmos, frac_reserve_2_atmos, carbon_2_atmos,
    !! carbon_2_litter_green_bgSoilPool and carbon_2_slowSoilPool to this routine.

    !! Check for landcover change in carbon mode
    IF (PRESENT(C_2_atmos)             .OR. PRESENT(C_2_litterGreenPools) .OR. &
        PRESENT(C_2_litterWoodPool_ag) .OR. PRESENT(C_2_litterWoodPool_bg)) THEN
        IF (.NOT. (PRESENT(C_2_atmos) .AND. PRESENT(C_2_litterGreenPools) .AND. &
             PRESENT(C_2_litterWoodPool_ag) .AND. PRESENT(C_2_litterWoodPool_bg))) &
           CALL finish('relocate_CarbonAndNitrogen()','at least one variable missing to handle land cover change')
        land_cover_change=.TRUE.
    END IF
 
    !! Check for nitrogen mode
    if(present(Npool_green)           .or. present(Npool_woods)           .or. &
       present(Npool_mobile)          .or. present(Npool_litter_green_ag) .or. &
       present(Npool_litter_wood_ag)  .or. present(Npool_litter_wood_bg)  .or. &
       present(Npool_litter_green_bg) .or. present(Npool_slow)            .or. &
       present(SMINN_pool)            .or. present(Npool_crop_harvest) ) then
       nitrogenMode=.true.
       if (.not. (present(Npool_green)           .and. present(Npool_woods)           .and. &
                  present(Npool_mobile)          .and. present(Npool_litter_green_ag) .and. &
                  present(Npool_litter_wood_ag)  .and. present(Npool_litter_wood_bg)  .and. &
                  present(Npool_litter_green_bg) .and. present(Npool_slow)            .and. &
                  present(SMINN_pool)            .and. present(Npool_crop_harvest))) &
          call finish('relocate_CarbonAndNitrogen()','at least one variable missing to handle nitrogen')
    end if

    !!check for .not. yasso mode
    if (present(Cpool_litter_green_ag) .or. present(Cpool_litter_green_bg) .or. &
        present(Cpool_litter_wood_ag)  .or. present(Cpool_litter_wood_bg)  .or. &
        present(Cpool_slow)) then
       if (.not.(present(Cpool_litter_green_ag) .and. present(Cpool_litter_green_bg) .and. &
          present(Cpool_litter_wood_ag)  .and. present(Cpool_litter_wood_bg)         .and. &
          present(Cpool_slow))) &
          call finish('relocate_CarbonAndNitrogen()','at least one variable missing to handle cbalance')
     end if

    !! Check if yasso pools shall be handled
    if(present(YCpool_acid_ag1)        .or. present(YCpool_acid_bg1)            .or. &
       present(YCpool_water_ag1)       .or. present(YCpool_water_bg1)           .or. &
       present(YCpool_ethanol_ag1)     .or. present(YCpool_ethanol_bg1)         .or. &
       present(YCpool_nonsoluble_ag1)  .or. present(YCpool_nonsoluble_bg1)      .or. &
       present(YCpool_humus_1)         .or.                                          &
       present(YCpool_acid_ag2)        .or. present(YCpool_acid_bg2)            .or. &
       present(YCpool_water_ag2)       .or. present(YCpool_water_bg2)           .or. &
       present(YCpool_ethanol_ag2)     .or. present(YCpool_ethanol_bg2)         .or. &
       present(YCpool_nonsoluble_ag2)  .or. present(YCpool_nonsoluble_bg2)      .or. &
       present(YCpool_humus_2)) then
       with_yasso=.true.
       if (.not. (present(YCpool_acid_ag1)        .and. present(YCpool_acid_bg1)        .and. &
                  present(YCpool_water_ag1)       .and. present(YCpool_water_bg1)       .and. &
                  present(YCpool_ethanol_ag1)     .and. present(YCpool_ethanol_bg1)     .and. &
                  present(YCpool_nonsoluble_ag1)  .and. present(YCpool_nonsoluble_bg1)  .and. &
                  present(YCpool_humus_1)                                               .and. &
                  present(YCpool_acid_ag2)        .and. present(YCpool_acid_bg2)        .and. &
                  present(YCpool_water_ag2)       .and. present(YCpool_water_bg2)       .and. &
                  present(YCpool_ethanol_ag2)     .and. present(YCpool_ethanol_bg2)     .and. &
                  present(YCpool_nonsoluble_ag2)  .and. present(YCpool_nonsoluble_bg2)  .and. &
                  present(YCpool_humus_2))) &
        call finish('relocate_CarbonAndNitrogen()','at least one variable missing to handle yasso pools')
    end if

    !! preparations

    ntiles = size(cf_old,DIM=2)
    nidx = size(cf_old,DIM=1)

    !! determine extending and shrinking tiles

    WHERE (cf_new(:,:) > cf_old(:,:))  
       tiles_extend(:,:) = .TRUE.
       cf_rescale(:,:) = 0._dp
    ELSEWHERE 
       tiles_extend(:,:) = .FALSE.
       cf_rescale(:,:) = (cf_old(:,:) - cf_new(:,:)) * veg_fract_correction(:,:)
    END WHERE

    !! 1. PART OF THE SUBROUTINE relocate_CarbonAndNitrogen !!
    !! determine amount of C and N released from living plant pools on shrinking tiles
    !! to atmosphere, litter (or anthropogenic pools), and soil (land_cover_change = true)
    !! or rescale living plant pools of shrinking tiles (land_cover_change = false ==> dynamic vegetation)

    IF (land_cover_change) THEN
       
       IF (lcc_scheme==1) THEN ! Standard JSBACH scheme

          !! C to be relocate from living vegetation pools of shrinking tiles to litter pools of expanding tiles
          !! ... notice cf_rescale is now 0 for expanding tiles -> no need for usage of mask
          flx_2_litter(:) = SUM(cf_rescale(:,:) * Cpool_woods(:,:),DIM=2) * (1._dp - frac_wood_2_atmos )
          C_2_litterWoodPool_ag(:) = flx_2_litter(:) * frac_wood_aboveGround
          C_2_litterWoodPool_bg(:) = flx_2_litter(:) * (1._dp - frac_wood_aboveGround)

          flx_2_litter(:) = SUM(cf_rescale(:,:) * Cpool_green(:,:), DIM=2) * (1._dp - frac_green_2_atmos) &
                          + SUM(cf_rescale(:,:) * Cpool_reserve(:,:), DIM=2) * (1._dp - frac_reserve_2_atmos)
          C_2_litterGreenPool_ag(:) = flx_2_litter(:) * frac_green_aboveGround
          C_2_litterGreenPool_bg(:) = flx_2_litter(:) * (1._dp - frac_green_aboveGround)
          C_2_litterGreenPools(:) = C_2_litterGreenPool_ag(:) + C_2_litterGreenPool_bg(:)         

          !! Living vegetation C from shrinking tiles to atmosphere
          C_2_atmos(:) = SUM((  frac_green_2_atmos   * Cpool_green  (:,:) &
                              + frac_wood_2_atmos    * Cpool_woods  (:,:) &
                              + frac_reserve_2_atmos * Cpool_reserve(:,:) &
                            ) * cf_rescale(:,:),DIM=2)
          
          IF (nitrogenMode) THEN
             
             !! N to be relocated from living vegetation pools of shrinking tiles to litter pools of expanding tiles

             !! ... Note: assuming cn_litter_wood >= cn_woods the amount of N from the wood pool to the litter pool is reduced here 
             !! ... by the factor cn_woods / cn_litter_wood (see below)
             flx_2_litter(:) =  SUM(cf_rescale(:,:) * Npool_woods(:,:),DIM=2) * (1.0_dp - frac_wood_2_atmos ) &
                                * cn_woods / cn_litter_wood
             Nitrogen_2_litterWoodPool_ag(:) = flx_2_litter(:) * frac_wood_aboveGround
             Nitrogen_2_litterWoodPool_bg(:) = flx_2_litter(:) * (1._dp - frac_wood_aboveGround)

             !! ... Note: reserve pool contains no N,
             !! ...       N in mobile pool goes to SMINN pool and atmos
             flx_2_litter(:) = SUM(cf_rescale(:,:) * Npool_green(:,:),DIM=2) * (1.0_dp - frac_green_2_atmos)
             Nitrogen_2_litterGreenPool_ag(:) = flx_2_litter(:) * frac_green_aboveGround
             Nitrogen_2_litterGreenPool_bg(:) = flx_2_litter(:) * ( 1._dp - frac_green_aboveGround)
             Nitrogen_2_litterGreenPools(:) = Nitrogen_2_litterGreenPool_ag(:) + Nitrogen_2_litterGreenPool_bg(:)

             !! ... Note: assuming cn_litter_wood >= cn_woods the surplus N from the transfer of wood to the woody litter pool (see above)
             !! ... is put into the soil mineral N pool
             Nitrogen_2_SMINNpool(:) = SUM(((1._dp - frac_wood_2_atmos) * Npool_woods(:,:)     & 
                                          * (1._dp - cn_woods / cn_litter_wood)                &
                                          + (1.0_dp - frac_mobile_2_atmos) * Npool_mobile(:,:) &
                                           ) * cf_rescale(:,:),DIM=2)
            
             !! Living N from shrinking tiles to atmosphere
             Nitrogen_2_atmos(:) = SUM(  (  frac_green_2_atmos  * Npool_green(:,:)    &
                                          + frac_wood_2_atmos   * Npool_woods(:,:)    &
                                          + frac_mobile_2_atmos * Npool_mobile(:,:) ) &
                                       * cf_rescale(:,:),DIM=2)
             
          END IF ! nitrogenMode
         
       ELSE ! lcc_scheme ne 1
       !! Anthropogenic C-reloc (lcc_scheme==2). details: T Kato et al. 2009 based on grand slam protocol by Houghton et al. (1983)

          !! by the call of C_loss_and_update_anthro_pools: put all the above ground litter ...
          !! ... and the above ground part of green pool, reserve pool, and wood pool in anthro pools and update anthro pools
          IF (.NOT. with_yasso) THEN
             CALL C_loss_and_update_anthro_pools(nidx,ntiles,lcc_scheme,lctlib,surface,       &
                                                 Cpool_green,Cpool_woods,Cpool_reserve,       &
                 Cpool_litter_green_ag        =  Cpool_litter_green_ag,                       &
                 Cpool_litter_wood_ag         =  Cpool_litter_wood_ag,                        &
                 Cpool_onSite                 =  Cpool_onSite     ,                           &
                 Cpool_paper                  =  Cpool_paper,                                 &
                 Cpool_construction           =  Cpool_construction,                          &
                 scale_fac                    =  cf_rescale,                                  &
                 C_2_onSite                   =  C_2_onSite,                                  &
                 C_2_paper                    =  C_2_paper,                                   &
                 C_2_construction             =  C_2_construction,                            &
                 C_onSite_2_atmos             =  C_onSite_2_atmos,                            &
                 C_paper_2_atmos              =  C_paper_2_atmos,                             &
                 C_construction_2_atmos       =  C_construction_2_atmos,                      &
                 C_2_atmos                    =  C_2_atmos)
          ELSE ! yasso
             CALL C_loss_and_update_anthro_pools(nidx,ntiles,lcc_scheme,lctlib,surface,       &
                                                 Cpool_green,Cpool_woods,Cpool_reserve,       &
                 YCpool_acid_ag1              =  YCpool_acid_ag1,                             &
                 YCpool_water_ag1             =  YCpool_water_ag1,                            &
                 YCpool_ethanol_ag1           =  YCpool_ethanol_ag1,                          &
                 YCpool_nonsoluble_ag1        =  YCpool_nonsoluble_ag1,                       &
                 YCpool_acid_ag2              =  YCpool_acid_ag2,                             &
                 YCpool_water_ag2             =  YCpool_water_ag2,                            &
                 YCpool_ethanol_ag2           =  YCpool_ethanol_ag2,                          &
                 YCpool_nonsoluble_ag2        =  YCpool_nonsoluble_ag2,                       &
                 Cpool_onSite                 =  Cpool_onSite,                                &
                 Cpool_paper                  =  Cpool_paper,                                 &
                 Cpool_construction           =  Cpool_construction,                          &
                 scale_fac                    =  cf_rescale,                                  &
                 C_2_onSite                   =  C_2_onSite,                                  &
                 C_2_paper                    =  C_2_paper,                                   &
                 C_2_construction             =  C_2_construction,                            &
                 C_onSite_2_atmos             =  C_onSite_2_atmos,                            &
                 C_paper_2_atmos              =  C_paper_2_atmos,                             &
                 C_construction_2_atmos       =  C_construction_2_atmos,                      &
                 C_2_atmos                    =  C_2_atmos)
          END IF ! yasso

          C_2_litterWoodPool_ag (:) = 0._dp
          C_2_litterWoodPool_bg (:) = SUM(cf_rescale(:,:) * Cpool_woods(:,:), DIM=2) * (1._dp - frac_wood_aboveGround)
          C_2_litterGreenPool_ag(:) = 0._dp
          C_2_litterGreenPool_bg(:) = SUM(cf_rescale(:,:) * (Cpool_green(:,:) + Cpool_reserve(:,:)), DIM=2) * &
                                      (1._dp - frac_green_aboveGround)
          C_2_litterGreenPools(:) = C_2_litterGreenPool_ag(:) + C_2_litterGreenPool_bg(:)

          IF (nitrogenMode) THEN
             !! Note: N_loss_lcc is called here only to determine the amount of N associated with the C put in the anthro pools.
             !! This N is assumed to be put immediately in the atmosphere (just to keep mass conservation)
             CALL N_loss_lcc(nidx, ntiles, lcc_scheme,                    &
                             Npool_green, Npool_woods, Npool_mobile,      &
                             scale_fac = cf_rescale(:,:),                 &
!!$ TR: now above ground litter is only redistributed by lcc and not put in onSite pool anymore
!!$                             Npool_litter_green_ag = Npool_litter_green_ag(:,:), &
!!$                             Npool_litter_wood_ag = Npool_litter_wood_ag(:,:), &
                             N_2_atmos = Nitrogen_2_atmos(:))

             !! ... Note: assuming cn_litter_wood >= cn_woods the amount of N from the wood pool to the litter pool is reduced here 
             !! ... by the factor cn_woods / cn_litter_wood
             Nitrogen_2_litterWoodPool_ag(:) = 0._dp
             Nitrogen_2_litterWoodPool_bg(:) = SUM(cf_rescale(:,:) * Npool_woods(:,:), DIM=2) &
                                               * (1._dp - frac_wood_aboveGround) * cn_woods / cn_litter_wood
             Nitrogen_2_litterGreenPool_ag(:) = 0._dp 
             Nitrogen_2_litterGreenPool_bg(:) = SUM(cf_rescale(:,:) * Npool_green(:,:), DIM=2) * (1._dp - frac_green_aboveGround)
             Nitrogen_2_litterGreenPools(:) = Nitrogen_2_litterGreenPool_ag(:) + Nitrogen_2_litterGreenPool_bg(:)
             !! ... Note: assuming cn_litter_wood >= cn_woods the surplus N from the transfer of wood to the woody litter pool
             !! ... is put into the soil mineral N pool
             Nitrogen_2_SMINNpool(:) = SUM(cf_rescale(:,:) * Npool_mobile(:,:), DIM=2) * (1._dp - frac_green_aboveGround) &
                                     + SUM(cf_rescale(:,:) * Npool_woods(:,:), DIM=2) * (1._dp - frac_wood_aboveGround) &
                                     * (1._dp - cn_woods / cn_litter_wood)
          END IF ! nitrogenMode

       ENDIF ! update anthro pools
    
    ELSE ! land_cover_change = false ==> cover frac changes from dynveg
    !! C-relocation from cover fraction change from the dynamic vegetation

       !! rescale living plant pools of shrinking tiles 
       
       WHERE(.NOT. tiles_extend(:,:) .AND. cf_new(:,:) >= fract_small)
          shrinkRatio(:,:) = cf_old(:,:) / cf_new(:,:)
       ELSEWHERE
          shrinkRatio(:,:) = 1._dp
       END WHERE
       
       Cpool_green  (:,:) = Cpool_green  (:,:) * shrinkRatio(:,:)
       Cpool_reserve(:,:) = Cpool_reserve(:,:) * shrinkRatio(:,:)
       Cpool_woods  (:,:) = Cpool_woods  (:,:) * shrinkRatio(:,:)

       if (nitrogenMode) then
         Npool_green  (:,:) = Npool_green  (:,:) * shrinkRatio(:,:)
         Npool_mobile (:,:) = Npool_mobile (:,:) * shrinkRatio(:,:)
         Npool_woods  (:,:) = Npool_woods  (:,:) * shrinkRatio(:,:)
       end if ! nitrogenMode
       
    END IF ! land cover change/dynveg 

    !! 2. PART OF THE SUBROUTINE relocate_CarbonAndNitrogen !!
    !! determine C and N to be relocated from slow pool and litter pools of shrinking to extending tiles (all schemes and reasons
    !! for cover frac changes)
    !! relocation from slow pool

    ! fluxes of soil carbon
    IF (.NOT. with_yasso) THEN
       Cslow_2_slow(:) = SUM(cf_rescale(:,:) * Cpool_slow(:,:),DIM=2)
    ELSE
       Chumus_1_2_humus_1(:) = SUM(cf_rescale(:,:) * YCpool_humus_1(:,:),DIM=2)
       Chumus_2_2_humus_2(:) = SUM(cf_rescale(:,:) * YCpool_humus_2(:,:),DIM=2)
    END IF
    Ccrop_harvest_2_crop_harvest(:) = SUM(cf_rescale(:,:) * Cpool_crop_harvest(:,:),DIM=2)

    IF (nitrogenMode) THEN
      Nslow_2_slow(:) = SUM(cf_rescale(:,:) * Npool_slow(:,:),DIM=2)
      Ncrop_harvest_2_crop_harvest(:) = SUM(cf_rescale(:,:) * Npool_crop_harvest(:,:),DIM=2)
      SMINN_2_SMINN(:) = SUM(cf_rescale(:,:) * SMINN_pool(:,:),DIM=2)
    ENDIF

    ! fluxes above ground litter
!!$ TR: In case of lcc_scheme=2 ag litter remains in ag litter pools (is only redistributed from shrinking to extending tiles)
!!$ TR: and is not put in onSite pool anymore.
!!$    IF (lcc_scheme==1 .OR. .NOT. land_cover_change) THEN !! dynveg or no anthro pools
       
       !! relocation from litter pools
       IF (.NOT. with_yasso) THEN
          ClitGreen_ag_2_litGreen_ag    (:) = SUM(cf_rescale(:,:) * Cpool_litter_green_ag(:,:),DIM=2)
          ClitterWood_ag_2_litterWood_ag(:) = SUM(cf_rescale(:,:) * Cpool_litter_wood_ag (:,:),DIM=2)
       ELSE ! yasso
          Cacid_ag1_2_acid_ag1            (:) = SUM(cf_rescale(:,:) * YCpool_acid_ag1      (:,:),DIM=2)
          Cwater_ag1_2_water_ag1          (:) = SUM(cf_rescale(:,:) * YCpool_water_ag1     (:,:),DIM=2)
          Cethanol_ag1_2_ethanol_ag1      (:) = SUM(cf_rescale(:,:) * YCpool_ethanol_ag1   (:,:),DIM=2)
          Cnonsoluble_ag1_2_nonsoluble_ag1(:) = SUM(cf_rescale(:,:) * YCpool_nonsoluble_ag1(:,:),DIM=2)
          Cacid_ag2_2_acid_ag2            (:) = SUM(cf_rescale(:,:) * YCpool_acid_ag2      (:,:),DIM=2)
          Cwater_ag2_2_water_ag2          (:) = SUM(cf_rescale(:,:) * YCpool_water_ag2     (:,:),DIM=2)
          Cethanol_ag2_2_ethanol_ag2      (:) = SUM(cf_rescale(:,:) * YCpool_ethanol_ag2   (:,:),DIM=2)
          Cnonsoluble_ag2_2_nonsoluble_ag2(:) = SUM(cf_rescale(:,:) * YCpool_nonsoluble_ag2(:,:),DIM=2)
       END IF ! yasso

       IF (nitrogenMode) THEN
          NlitGreen_ag_2_litGreen_ag    (:) = SUM(cf_rescale(:,:) * Npool_litter_green_ag(:,:),DIM=2)
          NlitterWood_ag_2_litterWood_ag(:) = SUM(cf_rescale(:,:) * Npool_litter_wood_ag (:,:),DIM=2)  
       END IF ! nitrogenMode
       
!!$    END IF ! dynveg or no anthro pools

    ! fluxes below ground litter
    IF (.NOT. with_yasso) THEN
       ClitGreen_bg_2_litGreen_bg    (:) = SUM(cf_rescale(:,:) * Cpool_litter_green_bg(:,:),DIM=2)
       ClitterWood_bg_2_litterWood_bg(:) = SUM(cf_rescale(:,:) * Cpool_litter_wood_bg (:,:),DIM=2)
    ELSE ! with yasso
       Cacid_bg1_2_acid_bg1            (:) = SUM(cf_rescale(:,:) * YCpool_acid_bg1      (:,:),DIM=2)
       Cwater_bg1_2_water_bg1          (:) = SUM(cf_rescale(:,:) * YCpool_water_bg1     (:,:),DIM=2)
       Cethanol_bg1_2_ethanol_bg1      (:) = SUM(cf_rescale(:,:) * YCpool_ethanol_bg1   (:,:),DIM=2)
       Cnonsoluble_bg1_2_nonsoluble_bg1(:) = SUM(cf_rescale(:,:) * YCpool_nonsoluble_bg1(:,:),DIM=2)
       Cacid_bg2_2_acid_bg2            (:) = SUM(cf_rescale(:,:) * YCpool_acid_bg2      (:,:),DIM=2)
       Cwater_bg2_2_water_bg2          (:) = SUM(cf_rescale(:,:) * YCpool_water_bg2     (:,:),DIM=2)
       Cethanol_bg2_2_ethanol_bg2      (:) = SUM(cf_rescale(:,:) * YCpool_ethanol_bg2   (:,:),DIM=2)
       Cnonsoluble_bg2_2_nonsoluble_bg2(:) = SUM(cf_rescale(:,:) * YCpool_nonsoluble_bg2(:,:),DIM=2)
    ENDIF ! yasso

    IF (nitrogenMode) THEN
       NlitGreen_bg_2_litGreen_bg    (:) = SUM(cf_rescale(:,:) * Npool_litter_green_bg(:,:),DIM=2)
       NlitterWood_bg_2_litterWood_bg(:) = SUM(cf_rescale(:,:) * Npool_litter_wood_bg (:,:),DIM=2)
    END IF ! nitrogenMode

    !! 3. PART OF THE SUBROUTINE relocate_CarbonAndNitrogen !!
    !! lower down the C and N density [mol/m^2(canopy)] for extending tiles (all schemes and reasons for cover frac changes)
    WHERE(tiles_extend(:,:) .AND. cf_new(:,:) >= fract_small)
       shrinkRatio(:,:) = cf_old(:,:) / cf_new(:,:)
    ELSEWHERE
       shrinkRatio(:,:) = 1._dp
    END WHERE

    Cpool_green(:,:)           = Cpool_green(:,:)           * shrinkRatio(:,:) 
    Cpool_reserve(:,:)         = Cpool_reserve(:,:)         * shrinkRatio(:,:)
    Cpool_woods(:,:)           = Cpool_woods(:,:)           * shrinkRatio(:,:)
    Cpool_crop_harvest(:,:)    = Cpool_crop_harvest(:,:)    * shrinkRatio(:,:)

    IF (.NOT. with_yasso) THEN
      Cpool_litter_green_ag(:,:) = Cpool_litter_green_ag(:,:) * shrinkRatio(:,:)
      Cpool_litter_green_bg(:,:) = Cpool_litter_green_bg(:,:) * shrinkRatio(:,:)
      Cpool_litter_wood_ag(:,:)  = Cpool_litter_wood_ag(:,:)  * shrinkRatio(:,:)
      Cpool_litter_wood_bg(:,:)  = Cpool_litter_wood_bg(:,:)  * shrinkRatio(:,:)
      Cpool_slow(:,:)            = Cpool_slow(:,:)            * shrinkRatio(:,:)
    ELSE ! yasso
      YCpool_acid_ag1      (:,:) = YCpool_acid_ag1      (:,:)   * shrinkRatio(:,:) 
      YCpool_water_ag1     (:,:) = YCpool_water_ag1     (:,:)   * shrinkRatio(:,:) 
      YCpool_ethanol_ag1   (:,:) = YCpool_ethanol_ag1   (:,:)   * shrinkRatio(:,:) 
      YCpool_nonsoluble_ag1(:,:) = YCpool_nonsoluble_ag1(:,:)   * shrinkRatio(:,:) 
      YCpool_acid_bg1      (:,:) = YCpool_acid_bg1      (:,:)   * shrinkRatio(:,:) 
      YCpool_water_bg1     (:,:) = YCpool_water_bg1     (:,:)   * shrinkRatio(:,:) 
      YCpool_ethanol_bg1   (:,:) = YCpool_ethanol_bg1   (:,:)   * shrinkRatio(:,:) 
      YCpool_nonsoluble_bg1(:,:) = YCpool_nonsoluble_bg1(:,:)   * shrinkRatio(:,:) 
      YCpool_humus_1       (:,:) = YCpool_humus_1       (:,:)   * shrinkRatio(:,:) 
      YCpool_acid_ag2      (:,:) = YCpool_acid_ag2      (:,:)   * shrinkRatio(:,:) 
      YCpool_water_ag2     (:,:) = YCpool_water_ag2     (:,:)   * shrinkRatio(:,:) 
      YCpool_ethanol_ag2   (:,:) = YCpool_ethanol_ag2   (:,:)   * shrinkRatio(:,:) 
      YCpool_nonsoluble_ag2(:,:) = YCpool_nonsoluble_ag2(:,:)   * shrinkRatio(:,:) 
      YCpool_acid_bg2      (:,:) = YCpool_acid_bg2      (:,:)   * shrinkRatio(:,:) 
      YCpool_water_bg2     (:,:) = YCpool_water_bg2     (:,:)   * shrinkRatio(:,:) 
      YCpool_ethanol_bg2   (:,:) = YCpool_ethanol_bg2   (:,:)   * shrinkRatio(:,:) 
      YCpool_nonsoluble_bg2(:,:) = YCpool_nonsoluble_bg2(:,:)   * shrinkRatio(:,:) 
      YCpool_humus_2       (:,:) = YCpool_humus_2       (:,:)   * shrinkRatio(:,:) 
    END IF ! yasso
     
    IF (nitrogenMode) THEN 
       Npool_green(:,:)           = Npool_green(:,:)           * shrinkRatio(:,:) 
       Npool_woods(:,:)           = Npool_woods(:,:)           * shrinkRatio(:,:)
       Npool_mobile(:,:)          = Npool_mobile(:,:)          * shrinkRatio(:,:)
       Npool_crop_harvest(:,:)    = Npool_crop_harvest(:,:)    * shrinkRatio(:,:)
       Npool_litter_green_ag(:,:) = Npool_litter_green_ag(:,:) * shrinkRatio(:,:)
       Npool_litter_green_bg(:,:) = Npool_litter_green_bg(:,:) * shrinkRatio(:,:)
       Npool_litter_wood_ag(:,:)  = Npool_litter_wood_ag(:,:)  * shrinkRatio(:,:)
       Npool_litter_wood_bg(:,:)  = Npool_litter_wood_bg(:,:)  * shrinkRatio(:,:)
       Npool_slow(:,:)            = Npool_slow(:,:)            * shrinkRatio(:,:)
       SMINN_pool(:,:)            = SMINN_pool(:,:)            * shrinkRatio(:,:)
    END IF

    !! 4. PART OF THE SUBROUTINE relocate_CarbonAndNitrogen !!
    !! distribute freed C and N to extending tiles

    !! sum of cover fractions added to cover fractions of extending tiles
    sum_cf_new_old(:) = 0._dp
    DO i = 1,ntiles
       WHERE (cf_new(:,i) - cf_old(:,i) > EPSILON(1._dp))
          sum_cf_new_old(:) = sum_cf_new_old(:) + (cf_new(:,i) - cf_old(:,i)) * veg_fract_correction(:,i)
       END WHERE
    END DO

    !! new temporary usage of cf_rescale for fewer computations and thus faster execution
    DO i = 1,ntiles
       WHERE (cf_new(:,i) - cf_old(:,i) > EPSILON(1._dp))
          cf_rescale(:,i) = (cf_new(:,i) - cf_old(:,i)) / (sum_cf_new_old(:) * cf_new(:,i))
       ELSEWHERE
          cf_rescale(:,i) = 0._dp
       END WHERE
    END DO

    !! distribute freed C and N to soil/humus/SMINN pools of extending tiles
    DO i = 1,ntiles

       !! soil pools are treated equally, no matter the source for the cf-change
       Cpool_crop_harvest(:,i) = Cpool_crop_harvest(:,i) + Ccrop_harvest_2_crop_harvest(:) * cf_rescale(:,i)
       IF (.NOT. with_yasso) THEN
          Cpool_slow(:,i) = Cpool_slow(:,i) + Cslow_2_slow(:) * cf_rescale(:,i)
       ELSE                   
          YCpool_humus_1(:,i) = YCpool_humus_1(:,i) + Chumus_1_2_humus_1(:) * cf_rescale(:,i)
          YCpool_humus_2(:,i) = YCpool_humus_2(:,i) + Chumus_2_2_humus_2(:) * cf_rescale(:,i)
       END IF
       IF (nitrogenMode) THEN
          Npool_slow(:,i) = Npool_slow(:,i) + Nslow_2_slow(:) * cf_rescale(:,i)
          Npool_crop_harvest(:,i) = Npool_crop_harvest(:,i) + Ncrop_harvest_2_crop_harvest(:) * cf_rescale(:,i)
          SMINN_pool(:,i) = SMINN_pool(:,i) + SMINN_2_SMINN(:) * cf_rescale(:,i)
       ENDIF

       IF (land_cover_change) THEN

          IF (nitrogenMode) THEN
             SMINN_pool(:,i) = SMINN_pool(:,i) + Nitrogen_2_SMINNpool(:) * cf_rescale(:,i)
          ENDIF

!!$ TR: In case of lcc_scheme=2 ag litter remains in ag litter pools (is only redistributed from shrinking to extending tiles)
!!$ TR: and is not put in onSite pool anymore.
!!$          IF (lcc_scheme==1) THEN

             !! distribute freed C and N to litter pools of extending tiles
             IF (.NOT. with_yasso) THEN
                WHERE (cf_new(:,i) - cf_old(:,i) > EPSILON(1._dp))
                   Cpool_litter_green_ag(:,i) =  &
                   Cpool_litter_green_ag(:,i) + (C_2_litterGreenPool_ag(:) + ClitGreen_ag_2_litGreen_ag    (:)) * cf_rescale(:,i)
                   Cpool_litter_wood_ag (:,i) =  &
                   Cpool_litter_wood_ag (:,i) + (C_2_litterWoodPool_ag (:) + ClitterWood_ag_2_litterWood_ag(:)) * cf_rescale(:,i)
                END WHERE
             ELSE ! yasso
                WHERE (cf_new(:,i) - cf_old(:,i) > EPSILON(1._dp))

                   YCpool_acid_ag1(:,i) =                                                        &
                        YCpool_acid_ag1(:,i) + (C_2_litterGreenPool_ag(:) * LeafLit_coef(:,i,1)  & ! Input leaf litter
                        + Cacid_ag1_2_acid_ag1(:)) * cf_rescale(:,i)                               ! Input from other tile

                   YCpool_acid_ag2(:,i) =                                                        &
                        YCpool_acid_ag2(:,i) + (C_2_litterWoodPool_ag(:)  * WoodLit_coef(:,i,1)  & ! Input woody litter
                        + Cacid_ag2_2_acid_ag2(:)) * cf_rescale(:,i)                                ! Input from other tile

                   YCpool_water_ag1(:,i) =                                                         &
                        YCpool_water_ag1(:,i)  + (C_2_litterGreenPool_ag(:) * LeafLit_coef(:,i,2)  &
                        + Cwater_ag1_2_water_ag1(:)) * cf_rescale(:,i)

                   YCpool_water_ag2(:,i) =                                                         &
                        YCpool_water_ag2(:,i)  + (C_2_litterWoodPool_ag(:)  * WoodLit_coef(:,i,2)  &
                        + Cwater_ag2_2_water_ag2(:)) * cf_rescale(:,i)

                   YCpool_ethanol_ag1(:,i) =                                                        &
                        YCpool_ethanol_ag1(:,i) + (C_2_litterGreenPool_ag(:) * LeafLit_coef(:,i,3)  &
                        + Cethanol_ag1_2_ethanol_ag1(:)) * cf_rescale(:,i)  

                   YCpool_ethanol_ag2(:,i) =                                                        &
                        YCpool_ethanol_ag2(:,i) + (C_2_litterWoodPool_ag(:)  * WoodLit_coef(:,i,3)  &
                        + Cethanol_ag2_2_ethanol_ag2(:)) * cf_rescale(:,i)

                   YCpool_nonsoluble_ag1(:,i) =                                                        &
                        YCpool_nonsoluble_ag1(:,i) + (C_2_litterGreenPool_ag(:) * LeafLit_coef(:,i,4)  &
                        + Cnonsoluble_ag1_2_nonsoluble_ag1(:)) * cf_rescale(:,i)

                   YCpool_nonsoluble_ag2(:,i) =                                                        &
                        YCpool_nonsoluble_ag2(:,i) + (C_2_litterWoodPool_ag(:)  * WoodLit_coef(:,i,4)  &
                        + Cnonsoluble_ag2_2_nonsoluble_ag2(:)) * cf_rescale(:,i)

                   YCpool_humus_1(:,i) =                                                               &
                        YCpool_humus_1(:,i) + C_2_litterGreenPool_ag(:) * LeafLit_coef(:,i,5) * cf_rescale(:,i)

                   YCpool_humus_2(:,i) =                                                               &
                        YCpool_humus_2(:,i) + C_2_litterWoodPool_ag(:)  * WoodLit_coef(:,i,5) * cf_rescale(:,i)

                END WHERE
             END IF ! yasso

             IF (nitrogenMode) THEN 
                WHERE (cf_new(:,i) - cf_old(:,i) > EPSILON(1._dp))
                   Npool_litter_green_ag(:,i)= &
                   Npool_litter_green_ag(:,i)+(Nitrogen_2_litterGreenPool_ag(:)+NlitGreen_ag_2_litGreen_ag    (:))*cf_rescale(:,i)
                   Npool_litter_wood_ag (:,i)= &
                   Npool_litter_wood_ag (:,i)+(Nitrogen_2_litterWoodPool_ag (:)+NlitterWood_ag_2_litterWood_ag(:))*cf_rescale(:,i)
                END WHERE
             END IF ! nitrogen
!!$          ELSE ! lcc_scheme ne 1
!!$             !! In the case lcc and anthropogenic pools are used, the freed carbon goes into the anthropogenic pools
!!$             !! instead of above ground litter pools (see above)
!!$          END IF ! lcc_scheme
          IF (.NOT. with_yasso) THEN
             WHERE (cf_new(:,i) - cf_old(:,i) > EPSILON(1._dp)) ! StW: Will C_2_litterXPool_bg always be 0 for shrinking tiles?
                Cpool_litter_green_bg(:,i) =  &
                Cpool_litter_green_bg(:,i) + (C_2_litterGreenPool_bg(:) + ClitGreen_bg_2_litGreen_bg    (:)) * cf_rescale(:,i)
                Cpool_litter_wood_bg (:,i) =  &
                Cpool_litter_wood_bg (:,i) + (C_2_litterWoodPool_bg (:) + ClitterWood_bg_2_litterWood_bg(:)) * cf_rescale(:,i)
             END WHERE
          ELSE ! with yasso
             WHERE (cf_new(:,i) - cf_old(:,i) > EPSILON(1._dp)) ! StW: Will C_2_litterXPool_bg always be 0 for shrinking tiles?
                YCpool_acid_bg1(:,i) =                                                        &
                     YCpool_acid_bg1(:,i) + (C_2_litterGreenPool_bg(:) * LeafLit_coef(:,i,1)  & ! Input leaf litter
                     + Cacid_bg1_2_acid_bg1(:)) * cf_rescale(:,i)                               ! Input from other tile

                YCpool_acid_bg2(:,i) =                                                        &
                     YCpool_acid_bg2(:,i) + (C_2_litterWoodPool_bg(:)  * WoodLit_coef(:,i,1)  & ! Input woody litter
                     + Cacid_bg2_2_acid_bg2(:)) * cf_rescale(:,i)                               ! Input from other tile

                YCpool_water_bg1(:,i) =                                                        &
                     YCpool_water_bg1(:,i) + (C_2_litterGreenPool_bg(:) * LeafLit_coef(:,i,2)  &
                     + Cwater_bg1_2_water_bg1(:)) * cf_rescale(:,i)

                YCpool_water_bg2(:,i) =                                                        &
                     YCpool_water_bg2(:,i) + (C_2_litterWoodPool_bg(:)  * WoodLit_coef(:,i,2)  &
                     + Cwater_bg2_2_water_bg2(:)) * cf_rescale(:,i)

                YCpool_ethanol_bg1(:,i) =                                                        &
                     YCpool_ethanol_bg1(:,i) + (C_2_litterGreenPool_bg(:) * LeafLit_coef(:,i,3)  &
                     + Cethanol_bg1_2_ethanol_bg1(:)) * cf_rescale(:,i)

                YCpool_ethanol_bg2(:,i) =                                                        &
                     YCpool_ethanol_bg2(:,i) + (C_2_litterWoodPool_bg(:)  * WoodLit_coef(:,i,3)  &
                     + Cethanol_bg2_2_ethanol_bg2(:)) * cf_rescale(:,i)

                YCpool_nonsoluble_bg1(:,i) =                                                        &
                     YCpool_nonsoluble_bg1(:,i) + (C_2_litterGreenPool_bg(:) * LeafLit_coef(:,i,4)  &
                     + Cnonsoluble_bg1_2_nonsoluble_bg1(:)) * cf_rescale(:,i)

                YCpool_nonsoluble_bg2(:,i) =                                                        &
                     YCpool_nonsoluble_bg2(:,i) + (C_2_litterWoodPool_bg(:)  * WoodLit_coef(:,i,4)  &
                     + Cnonsoluble_bg2_2_nonsoluble_bg2(:)) * cf_rescale(:,i)

                YCpool_humus_1(:,i) =                                                               &
                     YCpool_humus_1(:,i) + C_2_litterGreenPool_bg(:) * LeafLit_coef(:,i,5) * cf_rescale(:,i)

                YCpool_humus_2(:,i) =                                                               &
                     YCpool_humus_2(:,i) + C_2_litterWoodPool_bg(:)  * WoodLit_coef(:,i,5) * cf_rescale(:,i)

             END WHERE
          ENDIF ! yasso

          IF (nitrogenMode) THEN 
             WHERE (cf_new(:,i) - cf_old(:,i) > EPSILON(1._dp))
                Npool_litter_green_bg(:,i)= &
                Npool_litter_green_bg(:,i)+(Nitrogen_2_litterGreenPool_bg(:)+NlitGreen_bg_2_litGreen_bg    (:))*cf_rescale(:,i)
                Npool_litter_wood_bg (:,i)= &
                Npool_litter_wood_bg (:,i)+(Nitrogen_2_litterWoodPool_bg (:)+NlitterWood_bg_2_litterWood_bg(:))*cf_rescale(:,i)
             END WHERE
          END IF ! nitrogen

       ELSE ! land_cover_change = false ==> cover frac changes from dynveg

          IF (.NOT. with_yasso) THEN
             WHERE((cf_new(:,i) - cf_old(:,i)) > EPSILON(1._dp) .AND. sum_cf_new_old(:) > EPSILON(1._dp))
                Cpool_litter_green_ag(:,i) = Cpool_litter_green_ag(:,i) + ClitGreen_ag_2_litGreen_ag    (:) * cf_rescale(:,i)
                Cpool_litter_green_bg(:,i) = Cpool_litter_green_bg(:,i) + ClitGreen_bg_2_litGreen_bg    (:) * cf_rescale(:,i)
                Cpool_litter_wood_ag (:,i) = Cpool_litter_wood_ag (:,i) + ClitterWood_ag_2_litterWood_ag(:) * cf_rescale(:,i)
                Cpool_litter_wood_bg (:,i) = Cpool_litter_wood_bg (:,i) + ClitterWood_bg_2_litterWood_bg(:) * cf_rescale(:,i)
             END WHERE
          ELSE ! yasso
             WHERE((cf_new(:,i) - cf_old(:,i)) > EPSILON(1._dp) .AND. sum_cf_new_old(:) > EPSILON(1._dp))
                YCpool_acid_ag1(:,i)       = YCpool_acid_ag1(:,i)       + Cacid_ag1_2_acid_ag1                (:) * cf_rescale(:,i)
                YCpool_acid_bg1(:,i)       = YCpool_acid_bg1(:,i)       + Cacid_bg1_2_acid_bg1                (:) * cf_rescale(:,i)
                YCpool_water_ag1(:,i)      = YCpool_water_ag1(:,i)      + Cwater_ag1_2_water_ag1              (:) * cf_rescale(:,i)
                YCpool_water_bg1(:,i)      = YCpool_water_bg1(:,i)      + Cwater_bg1_2_water_bg1              (:) * cf_rescale(:,i)
                YCpool_ethanol_ag1(:,i)    = YCpool_ethanol_ag1(:,i)    + Cethanol_ag1_2_ethanol_ag1          (:) * cf_rescale(:,i)
                YCpool_ethanol_bg1(:,i)    = YCpool_ethanol_bg1(:,i)    + Cethanol_bg1_2_ethanol_bg1          (:) * cf_rescale(:,i)
                YCpool_nonsoluble_ag1(:,i) = YCpool_nonsoluble_ag1(:,i) + Cnonsoluble_ag1_2_nonsoluble_ag1    (:) * cf_rescale(:,i)
                YCpool_nonsoluble_bg1(:,i) = YCpool_nonsoluble_bg1(:,i) + Cnonsoluble_bg1_2_nonsoluble_bg1    (:) * cf_rescale(:,i)
                YCpool_acid_ag2(:,i)       = YCpool_acid_ag2(:,i)       + Cacid_ag2_2_acid_ag2                (:) * cf_rescale(:,i)
                YCpool_acid_bg2(:,i)       = YCpool_acid_bg2(:,i)       + Cacid_bg2_2_acid_bg2                (:) * cf_rescale(:,i)
                YCpool_water_ag2(:,i)      = YCpool_water_ag2(:,i)      + Cwater_ag2_2_water_ag2              (:) * cf_rescale(:,i)
                YCpool_water_bg2(:,i)      = YCpool_water_bg2(:,i)      + Cwater_bg2_2_water_bg2              (:) * cf_rescale(:,i)
                YCpool_ethanol_ag2(:,i)    = YCpool_ethanol_ag2(:,i)    + Cethanol_ag2_2_ethanol_ag2          (:) * cf_rescale(:,i)
                YCpool_ethanol_bg2(:,i)    = YCpool_ethanol_bg2(:,i)    + Cethanol_bg2_2_ethanol_bg2          (:) * cf_rescale(:,i)
                YCpool_nonsoluble_ag2(:,i) = YCpool_nonsoluble_ag2(:,i) + Cnonsoluble_ag2_2_nonsoluble_ag2    (:) * cf_rescale(:,i)
                YCpool_nonsoluble_bg2(:,i) = YCpool_nonsoluble_bg2(:,i) + Cnonsoluble_bg2_2_nonsoluble_bg2    (:) * cf_rescale(:,i)
             END WHERE
          END IF ! yasso
          IF (nitrogenMode) THEN 
             WHERE((cf_new(:,i) - cf_old(:,i)) > EPSILON(1._dp) .AND. sum_cf_new_old(:) > EPSILON(1._dp))
                Npool_litter_green_ag(:,i) = Npool_litter_green_ag(:,i) + NlitGreen_ag_2_litGreen_ag    (:) * cf_rescale(:,i)
                Npool_litter_green_bg(:,i) = Npool_litter_green_bg(:,i) + NlitGreen_bg_2_litGreen_bg    (:) * cf_rescale(:,i)
                Npool_litter_wood_ag (:,i) = Npool_litter_wood_ag (:,i) + NlitterWood_ag_2_litterWood_ag(:) * cf_rescale(:,i)
                Npool_litter_wood_bg (:,i) = Npool_litter_wood_bg (:,i) + NlitterWood_bg_2_litterWood_bg(:) * cf_rescale(:,i)
             END WHERE
          END IF ! nitrogen

       END IF ! land_cover_change

    END DO ! tile loop

  END SUBROUTINE relocate_CarbonAndNitrogen

  ! --- relocate_carbon_desert() ---------------------------------------------------------------------------------------------------
  !
  ! This routine models the consequences of desert dynamics for the carbon pools. More precisely:
  ! It is assumed that no carbon fluxes have to be considered as these are already taken into account by the other routines in this
  ! module. So, the amount of carbon per unit area is only rescaled to conserve the carbon mass.

  subroutine relocate_carbon_desert(nidx, ntiles, is_present, is_glacier, veg_ratio_max, veg_ratio_max_old, &
                                         Cpool_green, Cpool_woods, Cpool_reserve, Cpool_crop_harvest,   &
                                         Cpool_litter_green_ag, Cpool_litter_green_bg,                  &
                                         Cpool_litter_wood_ag,Cpool_litter_wood_bg, Cpool_slow,         &
                                         ! Nitrogen pools
                                         Npool_green, Npool_woods, Npool_mobile, Npool_crop_harvest,    &
                                         Npool_litter_green_ag, Npool_litter_green_bg,                  &
                                         Npool_litter_wood_ag, Npool_litter_wood_bg,                    &
                                         Npool_slow, SMINN_pool,                                        &
                                         ! Yasso litter and soil carbon pools
                                         YCpool_acid_ag1,                                                &
                                         YCpool_water_ag1,                                               &
                                         YCpool_ethanol_ag1,                                             &
                                         YCpool_nonsoluble_ag1,                                          &
                                         YCpool_acid_bg1,                                                &
                                         YCpool_water_bg1,                                               &
                                         YCpool_ethanol_bg1,                                             &
                                         YCpool_nonsoluble_bg1,                                          &
                                         YCpool_humus_1,                                                 &
                                         YCpool_acid_ag2,                                                &
                                         YCpool_water_ag2,                                               &
                                         YCpool_ethanol_ag2,                                             &
                                         YCpool_nonsoluble_ag2,                                          &
                                         YCpool_acid_bg2,                                                &
                                         YCpool_water_bg2,                                               &
                                         YCpool_ethanol_bg2,                                             &
                                         YCpool_nonsoluble_bg2,                                          &
                                         YCpool_humus_2,                                                 &
                                         ! Anthropogenic carbon pools
                                         Cpool_onSite,                                                   &
                                         Cpool_paper,                                                    &
                                         Cpool_construction,                                             &
                                         Cpool_onSite_harvest,                                           &
                                         Cpool_paper_harvest,                                            &
                                         Cpool_construction_harvest)
                                       
    USE mo_exception,        ONLY: finish

    integer,intent(in)     :: nidx                        !! Vector length
    integer,intent(in)     :: ntiles                      !! Number of tiles
    logical,intent(in)     :: is_present(:,:)             !! Logical mask for grid points treated by jsbach
    logical,intent(in)     :: is_glacier(:,:)             !! Glacier mask
    real(dp),intent(in)    :: veg_ratio_max(:)            !! Fractional cover of vegetated area
    real(dp),intent(in)    :: veg_ratio_max_old(:)        !! Fractional cover of vegetated area of last year
    real(dp),intent(inout) :: Cpool_green(:,:)            !! Value of green carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_woods(:,:)            !! Value of wood carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_reserve(:,:)          !! Value of reserve carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_crop_harvest(:,:)     !! Value of crop harvest carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: Cpool_litter_green_ag(:,:)    !! Above ground green litter carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: Cpool_litter_green_bg(:,:)    !! Below ground green litter carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: Cpool_litter_wood_ag(:,:)     !! Above ground woody litter carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: Cpool_litter_wood_bg(:,:)     !! Below ground woody litter carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: Cpool_slow(:,:)               !! Value of slow soil carbon pool [mol(C)/m^2(canopy)]
    ! Nitrogen
    real(dp),intent(inout),optional :: Npool_green(:,:)              !! Green nitrogen pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_woods(:,:)              !! Wood nitrogen pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_mobile(:,:)             !! Mobile nitrogen pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_crop_harvest(:,:)       !! Crop harvest nitrogen pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_litter_green_ag(:,:)    !! Above ground green litter nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_litter_green_bg(:,:)    !! Below ground green litter nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_litter_wood_ag(:,:)     !! Above ground woody litter nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_litter_wood_bg(:,:)     !! Below ground woody litter nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_slow(:,:)               !! Slow soil nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: SMINN_pool(:,:)               !! Mineral soil N [mol(N)/m^2(canopy)]
    !   YASSO 
    ! SIZE CLASS 1 : green litter
    !   above ground C pools
    real(dp),intent(inout),optional  :: YCpool_acid_ag1(:,:)          !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional  :: YCpool_water_ag1(:,:)         !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional  :: YCpool_ethanol_ag1(:,:)       !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional  :: YCpool_nonsoluble_ag1(:,:)    !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools  
    real(dp),intent(inout),optional  :: YCpool_acid_bg1(:,:)          !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional  :: YCpool_water_bg1(:,:)         !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional  :: YCpool_ethanol_bg1(:,:)       !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional  :: YCpool_nonsoluble_bg1(:,:)    !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional  :: YCpool_humus_1(:,:)           !! Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    ! SIZE CLASS 2 : woody litter
    !   above ground C pools
    real(dp),intent(inout),optional  :: YCpool_acid_ag2(:,:)          !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional  :: YCpool_water_ag2(:,:)         !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional  :: YCpool_ethanol_ag2(:,:)       !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional  :: YCpool_nonsoluble_ag2(:,:)    !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools  
    real(dp),intent(inout),optional  :: YCpool_acid_bg2(:,:)          !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional  :: YCpool_water_bg2(:,:)         !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional  :: YCpool_ethanol_bg2(:,:)       !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional  :: YCpool_nonsoluble_bg2(:,:)    !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional  :: YCpool_humus_2(:,:)           !! Humus litter pool for Yasso [mol(C)/m^2(canopy)]

    ! antropogenic C pools [mol m-2(veg)]
    REAL(dp), INTENT(inout), OPTIONAL :: Cpool_onSite(:)         !! C left on the ground by landcover changes [mol(C)/m^2(veg)]
    REAL(dp), INTENT(inout), OPTIONAL :: Cpool_paper(:)          !! C stored in short/interme. term ant. pool from landcover changes
    REAL(dp), INTENT(inout), OPTIONAL :: Cpool_construction(:)   !! C stored in long term anthro pool from landcover changes
    REAL(dp), INTENT(inout), OPTIONAL :: Cpool_onSite_harvest(:) !! C stored in short term anthro pool from harvest
    REAL(dp), INTENT(inout), OPTIONAL :: Cpool_paper_harvest(:)  !! C stored in intermediate term anthro pool from harvest
    REAL(dp), INTENT(inout), OPTIONAL :: Cpool_construction_harvest(:) !! C stored in long term anthro pool from harvest

    !! locals

    real(dp)  :: veg_ratio_old_new(nidx)          ! ratio of old fraction of vegetated land to new fraction of vegetated land
    LOGICAL   :: nitrogenMode = .FALSE.
    LOGICAL   :: with_yasso = .FALSE.
    LOGICAL   :: with_anthro_pools = .FALSE.

    !! check for nitrogen pools 
    IF(present(Npool_green)           .or. present(Npool_woods)           .or. &
       present(Npool_mobile)          .or. present(Npool_litter_green_ag) .or. &
       present(Npool_litter_wood_ag)  .or. present(Npool_litter_wood_bg)  .or. &
       present(Npool_litter_green_bg) .or. present(Npool_slow)            .or. &
       present(SMINN_pool)            .or. present(Npool_crop_harvest)) THEN
       nitrogenMode=.true.
       IF (.NOT. (PRESENT(Npool_green)           .AND. PRESENT(Npool_woods)           .AND. &
                  PRESENT(Npool_mobile)          .AND. PRESENT(Npool_litter_green_ag) .AND. &
                  PRESENT(Npool_litter_wood_ag)  .AND. PRESENT(Npool_litter_wood_bg)  .AND. &
                  PRESENT(Npool_litter_green_bg) .AND. PRESENT(Npool_slow)            .AND. &
                  PRESENT(SMINN_pool)            .AND. PRESENT(Npool_crop_harvest))) &
          CALL finish('relocate_carbon_desert()','at least one variable missing to handle nitrogen pools')       
    END IF

    !! check if yasso or cbalance litter and soil pools should be handled
    IF (PRESENT(YCpool_acid_ag1)        .OR. PRESENT(YCpool_acid_bg1)            .OR. &
        PRESENT(YCpool_water_ag1)       .OR. PRESENT(YCpool_water_bg1)           .OR. &
        PRESENT(YCpool_ethanol_ag1)     .OR. PRESENT(YCpool_ethanol_bg1)         .OR. &
        PRESENT(YCpool_nonsoluble_ag1)  .OR. PRESENT(YCpool_nonsoluble_bg1)      .OR. &
        PRESENT(YCpool_humus_1)         .OR.                                          &
        PRESENT(YCpool_acid_ag2)        .OR. PRESENT(YCpool_acid_bg2)            .OR. &
        PRESENT(YCpool_water_ag2)       .OR. PRESENT(YCpool_water_bg2)           .OR. &
        PRESENT(YCpool_ethanol_ag2)     .OR. PRESENT(YCpool_ethanol_bg2)         .OR. &
        PRESENT(YCpool_nonsoluble_ag2)  .OR. PRESENT(YCpool_nonsoluble_bg2)      .OR. &
        PRESENT(YCpool_humus_2)) THEN
        with_yasso=.true.
        IF (.NOT. (PRESENT(YCpool_acid_ag1)        .AND. PRESENT(YCpool_acid_bg1)        .AND. &
                   PRESENT(YCpool_water_ag1)       .AND. PRESENT(YCpool_water_bg1)       .AND. &
                   PRESENT(YCpool_ethanol_ag1)     .AND. PRESENT(YCpool_ethanol_bg1)     .AND. &
                   PRESENT(YCpool_nonsoluble_ag1)  .AND. PRESENT(YCpool_nonsoluble_bg1)  .AND. &
                   PRESENT(YCpool_humus_1)         .AND.                                       &
                   PRESENT(YCpool_acid_ag2)        .AND. PRESENT(YCpool_acid_bg2)        .AND. &
                   PRESENT(YCpool_water_ag2)       .AND. PRESENT(YCpool_water_bg2)       .AND. &
                   PRESENT(YCpool_ethanol_ag2)     .AND. PRESENT(YCpool_ethanol_bg2)     .AND. &
                   PRESENT(YCpool_nonsoluble_ag2)  .AND. PRESENT(YCpool_nonsoluble_bg2)  .AND. &
                   PRESENT(YCpool_humus_2))) THEN 
            CALL finish('relocate_carbon_desert()','at least one variable missing to handle yasso pools')
       END IF
    END IF

    !! check for anthropogenic pools 
    IF (PRESENT(Cpool_onSite)        .OR.  PRESENT(Cpool_paper)          .OR.  &
       PRESENT(Cpool_construction)   .OR.  PRESENT(Cpool_onSite_harvest) .OR.  &
       PRESENT(Cpool_paper_harvest)  .OR.  PRESENT(Cpool_construction_harvest)) THEN
       with_anthro_pools = .TRUE.
       IF (.NOT. (PRESENT(Cpool_onSite)        .AND.  PRESENT(Cpool_paper)          .AND.  &
                  PRESENT(Cpool_construction)  .AND.  PRESENT(Cpool_onSite_harvest) .AND.  &
                  PRESENT(Cpool_paper_harvest) .AND.  PRESENT(Cpool_construction_harvest))) THEN
          CALL finish('relocate_carbon_desert','at least one variable missing to handle anthropogenic pools')       
        END IF
    END IF

    WHERE (veg_ratio_max(:) > 0.5_dp*fract_small)
       veg_ratio_old_new(:) = veg_ratio_max_old(:) / veg_ratio_max(:)
    ELSEWHERE
       veg_ratio_old_new(:) = 1._dp                ! never happens on non-glacier land points
    ENDWHERE

    WHERE (.NOT. is_glacier(:,:) .AND. is_present(:,:))
       cpool_green(:,:) = cpool_green(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
       cpool_reserve(:,:) = cpool_reserve(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
       cpool_crop_harvest(:,:) = cpool_crop_harvest(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
       cpool_woods(:,:) = cpool_woods(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
    END WHERE

    IF (.NOT. with_yasso) THEN
       WHERE (.NOT. is_glacier(:,:) .AND. is_present(:,:))
          cpool_litter_green_ag(:,:) = cpool_litter_green_ag(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          cpool_litter_green_bg(:,:) = cpool_litter_green_bg(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          cpool_litter_wood_ag(:,:) = cpool_litter_wood_ag(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          cpool_litter_wood_bg(:,:) = cpool_litter_wood_bg(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          cpool_slow(:,:) = cpool_slow(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
       END WHERE
    ELSE ! yasso
       WHERE (.NOT. is_glacier(:,:) .AND. is_present(:,:))
          YCpool_acid_ag1       (:,:) = YCpool_acid_ag1       (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          YCpool_water_ag1      (:,:) = YCpool_water_ag1      (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          YCpool_ethanol_ag1    (:,:) = YCpool_ethanol_ag1    (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          YCpool_nonsoluble_ag1 (:,:) = YCpool_nonsoluble_ag1 (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          YCpool_humus_1        (:,:) = YCpool_humus_1        (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          YCpool_acid_bg1       (:,:) = YCpool_acid_bg1       (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          YCpool_water_bg1      (:,:) = YCpool_water_bg1      (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          YCpool_ethanol_bg1    (:,:) = YCpool_ethanol_bg1    (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          YCpool_nonsoluble_bg1 (:,:) = YCpool_nonsoluble_bg1 (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          YCpool_acid_ag2       (:,:) = YCpool_acid_ag2       (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          YCpool_water_ag2      (:,:) = YCpool_water_ag2      (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          YCpool_ethanol_ag2    (:,:) = YCpool_ethanol_ag2    (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          YCpool_nonsoluble_ag2 (:,:) = YCpool_nonsoluble_ag2 (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          YCpool_humus_2        (:,:) = YCpool_humus_2        (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          YCpool_acid_bg2       (:,:) = YCpool_acid_bg2       (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          YCpool_water_bg2      (:,:) = YCpool_water_bg2      (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          YCpool_ethanol_bg2    (:,:) = YCpool_ethanol_bg2    (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          YCpool_nonsoluble_bg2 (:,:) = YCpool_nonsoluble_bg2 (:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
       END WHERE
    END IF ! yasso

    IF (with_anthro_pools) THEN
       Cpool_onSite(:)               = Cpool_onSite(:)               * veg_ratio_old_new(:)
       Cpool_paper(:)                = Cpool_paper(:)                * veg_ratio_old_new(:)
       Cpool_construction(:)         = Cpool_construction(:)         * veg_ratio_old_new(:)
       Cpool_onSite_harvest(:)       = Cpool_onSite_harvest(:)       * veg_ratio_old_new(:)
       Cpool_paper_harvest(:)        = Cpool_paper_harvest(:)        * veg_ratio_old_new(:)
       Cpool_construction_harvest(:) = Cpool_construction_harvest(:) * veg_ratio_old_new(:)
    END IF

    IF (nitrogenMode) THEN
       WHERE (.NOT. is_glacier(:,:) .AND. is_present(:,:))
          npool_green(:,:) = npool_green(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          npool_mobile(:,:) = npool_mobile(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          npool_crop_harvest(:,:) = npool_crop_harvest(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          npool_woods(:,:) = npool_woods(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          npool_litter_green_ag(:,:) = npool_litter_green_ag(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          npool_litter_green_bg(:,:) = npool_litter_green_bg(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          npool_litter_wood_ag(:,:) = npool_litter_wood_ag(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          npool_litter_wood_bg(:,:) = npool_litter_wood_bg(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          npool_slow(:,:) = npool_slow(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2)
          sminn_pool(:,:) = sminn_pool(:,:) * SPREAD(veg_ratio_old_new(:),NCOPIES=ntiles,DIM=2) 
       END WHERE
    END IF ! nitrogenMode

  end subroutine relocate_carbon_desert


  ! --- relocate_carbon_fire() ---------------------------------------------------------------------------------------------------
  !
  ! This routine models the consequences of fire for the carbon pools. More precisely:
  ! It is assumed that for the burned area the carbon from the above ground litter pools (Cpool_litter_green_ag,
  ! Cpool_litter_wood_ag) is completely released to the atmosphere and the carbon from the living plant pools
  ! (Cpool_green, Cpool_reserve, Cpool_woods) is partly released to the atmosphere and partly relocated
  ! into the litter pools.
  !
  subroutine relocate_carbon_fire(nidx, ntiles, with_yasso, ldiag, frac_burned, cf,             &
                                       veg_fract_correction, frac_wood_2_atmos_fire,            &
                                       Cpool_green, Cpool_woods, Cpool_reserve,                 &
                                       carbon_2_GreenLitterPools,                               &
                                       carbon_2_WoodLitterPools,                                &
                                       carbon_2_atmos_tiled_ext, carbon_2_atmos,                &
                                       Cpool_litter_green_ag,Cpool_litter_green_bg,             &
                                       Cpool_litter_wood_ag, Cpool_litter_wood_bg,              &
                                       ! Yasso litter and soil carbon pools
                                       YCpool_acid_ag1,      YCpool_acid_ag2,                   &
                                       YCpool_water_ag1,     YCpool_water_ag2,                  &
                                       YCpool_ethanol_ag1,   YCpool_ethanol_ag2,                &
                                       YCpool_nonsoluble_ag1,YCpool_nonsoluble_ag2,             &
                                       YCpool_acid_bg1,      YCpool_acid_bg2,                   &
                                       YCpool_water_bg1,     YCpool_water_bg2,                  &
                                       YCpool_ethanol_bg1,   YCpool_ethanol_bg2,                &
                                       YCpool_nonsoluble_bg1,YCpool_nonsoluble_bg2,             &
                                       YCpool_humus_1,       YCpool_humus_2,                    &
                                       LeafLitcoef,          WoodLitcoef,                       &
                                       !Nitrogen Pools and fluxes
                                       Npool_green, Npool_woods, Npool_mobile, SMINN_pool,      &
                                       Npool_litter_green_ag,Npool_litter_green_bg,             &
                                       Npool_litter_wood_ag,Npool_litter_wood_bg,               &
                                       nitrogen_2_GreenLitterPools, nitrogen_2_WoodLitterPools, &
                                       nitrogen_2_atmos, nitrogen_2_sminn)
!!
!! CHR 09-11-03: DIESE ROUTINE MUSS AUS DYNVEG HERAUSGEZOGEN WERDEN UM AUCH AGRARFLAECHEN IN DAS FEUER MIT EINZUBEZIEHEN
!!  

    USE mo_exception,        ONLY: finish

    integer,intent(in)     :: nidx                        !! Vector length
    integer,intent(in)     :: ntiles                      !! Number of tiles
    logical,intent(in)     :: with_yasso
    logical,intent(in)     :: ldiag
    real(dp),intent(in)    :: frac_burned(:,:)            !! Fraction of the vegetated area of each tile burned since the last call
                                                          !!    of this routine
    real(dp),intent(in)    :: cf(:,:)                     !! Cover fractions
    real(dp),intent(in)    :: veg_fract_correction(:,:)   !! Correction factor for cover fractions 1-exp(-LAI_max/2) (accounts for
                                                          !!    sparseness of vegetation)
    real(dp),intent(in)    :: frac_wood_2_atmos_fire      !! Fraction of above ground wood immediately emitted to atmosphere by fire
    real(dp),intent(inout) :: Cpool_green(:,:)            !! Value of green carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_woods(:,:)            !! Value of wood carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_reserve(:,:)          !! Value of reserve carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: carbon_2_GreenLitterPools(:)!! Amount of carbon relocated by wind damage (in) and by wind and fire
                                                          !! .. (out) to the green litter pools [mol(C)/m^2(vegetated area)]
    real(dp),intent(inout) :: carbon_2_WoodLitterPools(:) !! Amount of carbon relocated by wind damage (in) and by wind and fire
                                                          !! .. (out) to the wood litter pools [mol(C)/m^2(vegetated area)]
    real(dp),intent(inout),target :: carbon_2_atmos_tiled_ext(:,:) !! Carbon immedieatly released by fire per tile [mol(C)m-2(veg)]
    real(dp),intent(out)   :: carbon_2_atmos(:)           !! Amount of carbon immediately released by fire
    ! Cbalance litter pools
    real(dp),intent(inout),optional :: Cpool_litter_green_ag(:,:)    !! Above ground green litter carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: Cpool_litter_green_bg(:,:)    !! Below ground green litter carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: Cpool_litter_wood_ag(:,:)     !! Wood litter carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: Cpool_litter_wood_bg(:,:)     !! Wood litter carbon pool [mol(C)/m^2(canopy)]

    !   YASSO 
    ! SIZE CLASS 1 ; GREEN
    !   above ground C pools
    real(dp),intent(inout),optional :: YCpool_acid_ag1(:,:)          !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: YCpool_water_ag1(:,:)         !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: YCpool_ethanol_ag1(:,:)       !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: YCpool_nonsoluble_ag1(:,:)    !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools
    real(dp),intent(inout),optional :: YCpool_acid_bg1(:,:)          !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: YCpool_water_bg1(:,:)         !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: YCpool_ethanol_bg1(:,:)       !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: YCpool_nonsoluble_bg1(:,:)    !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: YCpool_humus_1(:,:)           !! Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    ! SIZE CLASS 2 ; WOOD
    !   above ground C pools
    real(dp),intent(inout),optional :: YCpool_acid_ag2(:,:)          !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: YCpool_water_ag2(:,:)         !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: YCpool_ethanol_ag2(:,:)       !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: YCpool_nonsoluble_ag2(:,:)    !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools  
    real(dp),intent(inout),optional :: YCpool_acid_bg2(:,:)          !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: YCpool_water_bg2(:,:)         !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: YCpool_ethanol_bg2(:,:)       !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: YCpool_nonsoluble_bg2(:,:)    !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: YCpool_humus_2(:,:)           !! Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    !   parameters 
    real(dp),intent(in),optional    :: LeafLitcoef(:,:,:)            !! Factor to spread non woody litterfall into yasso pools [ ]
    real(dp),intent(in),optional    :: WoodLitcoef(:,:,:)            !! Factor to spread woody litterfall into yasso pools [ ]

    ! Nitrogen
    real(dp),intent(inout),optional :: Npool_green(:,:)              !! Green nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_woods(:,:)              !! Wood nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_mobile(:,:)             !! Mobile nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: SMINN_pool(:,:)               !! Mineral nitrogen pool [mol(N)/m^2(canopy)]    
    real(dp),intent(inout),optional :: Npool_litter_green_ag(:,:)    !! Above ground green litter nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_litter_green_bg(:,:)    !! Below ground green litter nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_litter_wood_ag(:,:)     !! Wood litter nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: Npool_litter_wood_bg(:,:)     !! Wood litter nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional :: nitrogen_2_GreenLitterPools(:)!! Amount of nitrogen relocated by fire to the green litter
                                                                     !!    pools [mol(N)/m^2(vegetated area)]
    real(dp),intent(inout),optional :: nitrogen_2_WoodLitterPools(:) !! Amount of nitrogen relocated by fire and wind damage 
                                                                     !!    to the wood litter pools [mol(N)/m^2(vegetated area)]
    real(dp),intent(out),  optional :: nitrogen_2_atmos(:)           !! Amount of nitrogen immediately released by fire and wind 
                                                                     !!    damage to the atmosphere [mol(N)/m^2(vegetated area)] 
    real(dp),intent(inout),optional :: nitrogen_2_sminn(:)           !! Amount of nitrogen relocated by fire and wind damage to the
                                                                     !!    soil mineral N pool [mol(N)/m^2(vegetated area)]

    !! locals

    REAL(dp), POINTER :: carbon_2_atmos_tiled_local(:,:), carbon_2_atmos_tiled(:,:)
    real(dp)  ::  carbon_2_GreenLitterPools_tiled(nidx,ntiles)
    real(dp)  ::  carbon_2_WoodLitterPools_tiled(nidx,ntiles)
    real(dp)  ::  N_2_GreenLitterPools_tiled(nidx,ntiles)
    real(dp)  ::  N_2_WoodLitterPools_tiled(nidx,ntiles)
    real(dp)  ::  N_2_atmos_tiled(nidx,ntiles)
    real(dp)  ::  N_2_sminn_tiled(nidx,ntiles)
    logical   ::  nitrogenMode

    nitrogenMode=.false.
    if(present(Npool_green)           .or. present(Npool_woods)           .or. &
       present(Npool_mobile)          .or. present(Npool_litter_green_ag) .or. &
       present(Npool_litter_wood_ag)  .or. present(Npool_litter_wood_bg)  .or. &
       present(Npool_litter_green_bg) .or. present(SMINN_pool) ) then
       nitrogenMode=.true.
       if (.not. (present(Npool_green)           .and. present(Npool_woods)           .and. &
                  present(Npool_mobile)          .and. present(Npool_litter_green_ag) .and. &
                  present(Npool_litter_wood_ag)  .and. present(Npool_litter_wood_bg)  .and. &
                  present(Npool_litter_green_bg) .and. present(SMINN_pool) ) ) &
                  CALL finish('relocate_carbon_fire()','at least one variable missing to handle nitrogen dynamics')
    end if    

    !! preparations
    IF (ldiag) THEN
      carbon_2_atmos_tiled => carbon_2_atmos_tiled_ext
    ELSE
      ALLOCATE(carbon_2_atmos_tiled_local(nidx,ntiles))
      carbon_2_atmos_tiled => carbon_2_atmos_tiled_local
    ENDIF
    carbon_2_GreenLitterPools_tiled(:,:) = 0._dp
    carbon_2_WoodLitterPools_tiled(:,:) = 0._dp
    N_2_GreenLitterPools_tiled(:,:) = 0._dp
    N_2_WoodLitterPools_tiled(:,:) = 0._dp
    N_2_atmos_tiled(:,:) = 0._dp   
    N_2_sminn_tiled(:,:)   = 0._dp

    !! diagnose amount of carbon released from the green and reserve pool to the green litter pools
    !! diagnose amount of carbon released from the wood pool to the woody litter pools
    !! determine amount of carbon released from the wood pool, living tissue pools, and litter pools to the atmosphere

    carbon_2_GreenLitterPools_tiled(:,:) = (Cpool_green(:,:) + Cpool_reserve(:,:)) * (1._dp - frac_green_aboveGround) * &
                                           cf(:,:) * veg_fract_correction(:,:) * frac_burned(:,:)
    carbon_2_WoodLitterPools_tiled(:,:) = Cpool_woods(:,:) * (1._dp - frac_wood_2_atmos_fire * frac_wood_aboveGround) * &
                                          cf(:,:) * veg_fract_correction(:,:) * frac_burned(:,:)


    IF (.NOT. with_yasso) THEN
       carbon_2_atmos_tiled(:,:) = ((Cpool_green(:,:) + Cpool_reserve(:,:)) * frac_green_aboveGround + &
                                   Cpool_litter_green_ag(:,:) + Cpool_litter_wood_ag(:,:) + &
                                   Cpool_woods(:,:) * frac_wood_2_atmos_fire * frac_wood_aboveGround) * &
                                   cf(:,:) * veg_fract_correction(:,:) * frac_burned(:,:)
    ELSE
       carbon_2_atmos_tiled(:,:) = ((Cpool_green(:,:) + Cpool_reserve(:,:)) * frac_green_aboveGround      + &
                                   YCpool_water_ag1(:,:) + YCpool_acid_ag1(:,:) + YCPool_ethanol_ag1(:,:) + &
                                   YCpool_nonsoluble_ag1(:,:) +                                             &
                                   YCpool_water_ag2(:,:) + YCpool_acid_ag2(:,:) + YCPool_ethanol_ag2(:,:) + &
                                   YCpool_nonsoluble_ag2(:,:) +                                             &
                                   Cpool_woods(:,:) * frac_wood_2_atmos_fire * frac_wood_aboveGround) *     &
                                   cf(:,:) * veg_fract_correction(:,:) * frac_burned(:,:)
    END IF
    
    if (nitrogenMode) then
      !Npool_mobile above ground to atmosphere
      !Npool_mobile below ground to sminn_pool
      !Npool_litter_ag to atmosphere
      !SMINN pool not affected by fire, only through N released by other pools
        N_2_atmos_tiled(:,:) = ((Npool_green(:,:) + Npool_mobile(:,:)) * frac_green_aboveGround + &
                                   Npool_litter_green_ag(:,:) + Npool_litter_wood_ag(:,:) + &                
                                   Npool_woods(:,:) * frac_wood_2_atmos_fire * frac_wood_aboveGround) * &    
                                   cf(:,:) * veg_fract_correction(:,:) * frac_burned(:,:) 
        N_2_GreenLitterPools_tiled(:,:) = Npool_green(:,:) * (1._dp - frac_green_aboveGround) * &
                                              cf(:,:) * veg_fract_correction(:,:) * frac_burned(:,:)
        N_2_WoodLitterPools_tiled(:,:) = Cpool_woods(:,:) * (1._dp - frac_wood_2_atmos_fire * frac_wood_aboveGround) * &
                                             1._dp/cn_litter_wood * &
                                             cf(:,:) * veg_fract_correction(:,:) * frac_burned(:,:)
        N_2_sminn_tiled(:,:) = ( Npool_mobile(:,:)* (1._dp - frac_green_aboveGround) + &
                                  Cpool_woods(:,:)  * (1._dp - frac_wood_2_atmos_fire * frac_wood_aboveGround) * &
                                  (1._dp/cn_woods - 1._dp/cn_litter_wood) )                    * &
                                    cf(:,:) * veg_fract_correction(:,:) * frac_burned(:,:)
    end if

    !! sum carbon fluxes for all tiles

    carbon_2_GreenLitterPools(:) = carbon_2_GreenLitterPools(:) + SUM(carbon_2_GreenLitterPools_tiled(:,:),DIM=2)
    carbon_2_WoodLitterPools(:) = carbon_2_WoodLitterPools(:) + SUM(carbon_2_WoodLitterPools_tiled(:,:),DIM=2)
    carbon_2_atmos(:) = SUM(carbon_2_atmos_tiled(:,:),DIM=2)
    
    if (nitrogenMode) then
      nitrogen_2_GreenLitterPools(:) = nitrogen_2_GreenLitterPools(:) + SUM(N_2_GreenLitterPools_tiled(:,:),DIM=2)
      nitrogen_2_WoodLitterPools(:) = nitrogen_2_WoodLitterPools(:) + SUM(N_2_WoodLitterPools_tiled(:,:),DIM=2)
      nitrogen_2_atmos(:) = SUM(N_2_atmos_tiled(:,:),DIM=2)
      nitrogen_2_sminn(:) = SUM(N_2_sminn_tiled(:,:),DIM=2)
    end if 


    IF (.NOT. with_yasso) THEN
       !! lower down the carbon density of above ground litter pools
       Cpool_litter_green_ag(:,:) = Cpool_litter_green_ag(:,:) * (1._dp - frac_burned(:,:))
       Cpool_litter_wood_ag (:,:) = Cpool_litter_wood_ag (:,:) * (1._dp - frac_burned(:,:))

       !! transfer carbon from the living tissue pools (green and reserve) to the below ground green litter pool
       Cpool_litter_green_bg(:,:) = Cpool_litter_green_bg(:,:) + &
            (Cpool_green(:,:) + Cpool_reserve(:,:)) * (1._dp - frac_green_aboveGround) * frac_burned(:,:)

       !! transfer carbon from the wood pool to the above ground woody litter pool
       Cpool_litter_wood_ag(:,:) = Cpool_litter_wood_ag(:,:) + &
            Cpool_woods(:,:) * ((1._dp - frac_wood_2_atmos_fire) * frac_wood_aboveGround) * frac_burned(:,:)

       !! transfer carbon from the wood pool to the below ground woody litter pool
       Cpool_litter_wood_bg(:,:) = Cpool_litter_wood_bg(:,:) + &
            Cpool_woods(:,:) * (1._dp - frac_wood_aboveGround) * frac_burned(:,:)
    ELSE

       ! Lower down density of aboveground litter pools
       YCpool_acid_ag1      (:,:) = YCpool_acid_ag1       (:,:) * (1._dp - frac_burned(:,:))   ! lower down density
       YCpool_water_ag1     (:,:) = YCpool_water_ag1      (:,:) * (1._dp - frac_burned(:,:))   ! lower down density
       YCpool_ethanol_ag1   (:,:) = YCpool_ethanol_ag1    (:,:) * (1._dp - frac_burned(:,:))   ! lower down density
       YCpool_nonsoluble_ag1(:,:) = YCpool_nonsoluble_ag1 (:,:) * (1._dp - frac_burned(:,:))   ! lower down density

       YCpool_acid_ag2      (:,:) = YCpool_acid_ag2       (:,:) * (1._dp - frac_burned(:,:))   ! lower down density
       YCpool_water_ag2     (:,:) = YCpool_water_ag2      (:,:) * (1._dp - frac_burned(:,:))   ! lower down density
       YCpool_ethanol_ag2   (:,:) = YCpool_ethanol_ag2    (:,:) * (1._dp - frac_burned(:,:))   ! lower down density
       YCpool_nonsoluble_ag2(:,:) = YCpool_nonsoluble_ag2 (:,:) * (1._dp - frac_burned(:,:))   ! lower down density

       ! Add green and reserve vegetation carbon to belowground green litter pools
       YCpool_acid_bg1      (:,:) = YCpool_acid_bg1       (:,:) +                                   &
                                    ( 1._dp-frac_green_aboveGround) * frac_burned(:,:)              &
                                     * (Cpool_green(:,:) + Cpool_reserve(:,:)) * LeafLitcoef(:,:,1)
       YCpool_water_bg1     (:,:) = YCpool_water_bg1      (:,:) +                                   &
                                    ( 1._dp-frac_green_aboveGround) * frac_burned(:,:)              &
                                     * (Cpool_green(:,:) + Cpool_reserve(:,:)) * LeafLitcoef(:,:,2)
       YCpool_ethanol_bg1   (:,:) = YCpool_ethanol_bg1    (:,:) +                                   &
                                    ( 1._dp-frac_green_aboveGround) * frac_burned(:,:)              &
                                     * (Cpool_green(:,:) + Cpool_reserve(:,:)) * LeafLitcoef(:,:,3)
       YCpool_nonsoluble_bg1(:,:) = YCpool_nonsoluble_bg1 (:,:) +                                   &
                                    ( 1._dp-frac_green_aboveGround) * frac_burned(:,:)              &
                                     * (Cpool_green(:,:) + Cpool_reserve(:,:)) * LeafLitcoef(:,:,4)
       YCpool_humus_1       (:,:) = YCpool_humus_1        (:,:) +                                   &
                                    ( 1._dp-frac_green_aboveGround) * frac_burned(:,:)              &
                                     * (Cpool_green(:,:) + Cpool_reserve(:,:)) * LeafLitcoef(:,:,5)

       ! Add the remaings of burned wood to the aboveground wood pools:
       YCpool_acid_ag2      (:,:) = YCpool_acid_ag2      (:,:)                                    & 
                                   +  ((1._dp - frac_wood_2_atmos_fire) * frac_wood_aboveGround)  &
                                    * frac_burned(:,:) * Cpool_woods(:,:) * WoodLitcoef(:,:,1)
       YCpool_water_ag2     (:,:) = YCpool_water_ag2                                              & 
                                   +  ((1._dp - frac_wood_2_atmos_fire) * frac_wood_aboveGround)  &
                                    * frac_burned(:,:) * Cpool_woods(:,:) * WoodLitcoef(:,:,2)
       YCpool_ethanol_ag2   (:,:) = YCpool_ethanol_ag2   (:,:)                                    & 
                                   +  ((1._dp - frac_wood_2_atmos_fire) * frac_wood_aboveGround)  &
                                    * frac_burned(:,:) * Cpool_woods(:,:) * WoodLitcoef(:,:,3)
       YCpool_nonsoluble_ag2(:,:) = YCpool_nonsoluble_ag2(:,:)                                    & 
                                    +  ((1._dp - frac_wood_2_atmos_fire) * frac_wood_aboveGround) &
                                     * frac_burned(:,:) * Cpool_woods(:,:) * WoodLitcoef(:,:,4)
       YCpool_humus_2       (:,:) = YCpool_humus_2       (:,:)                                    &
                                    +  ((1._dp - frac_wood_2_atmos_fire) * frac_wood_aboveGround) &
                                     * frac_burned(:,:) * Cpool_woods(:,:) * WoodLitcoef(:,:,5)

       ! Add the remaing of burned wood to the belowground wood pools:
       YCpool_acid_bg2      (:,:) = &
       YCpool_acid_bg2      (:,:) + (1._dp-frac_wood_aboveGround) * frac_burned(:,:) * Cpool_woods(:,:) * WoodLitcoef(:,:,1)
       YCpool_water_bg2     (:,:) = &
       YCpool_water_bg2     (:,:) + (1._dp-frac_wood_aboveGround) * frac_burned(:,:) * Cpool_woods(:,:) * WoodLitcoef(:,:,2)
       YCpool_ethanol_bg2   (:,:) = &
       YCpool_ethanol_bg2   (:,:) + (1._dp-frac_wood_aboveGround) * frac_burned(:,:) * Cpool_woods(:,:) * WoodLitcoef(:,:,3)
       YCpool_nonsoluble_bg2(:,:) = &
       YCpool_nonsoluble_bg2(:,:) + (1._dp-frac_wood_aboveGround) * frac_burned(:,:) * Cpool_woods(:,:) * WoodLitcoef(:,:,4)
       YCpool_humus_2       (:,:) = &
       YCpool_humus_2       (:,:) + (1._dp-frac_wood_aboveGround) * frac_burned(:,:) * Cpool_woods(:,:) * WoodLitcoef(:,:,5)
    END IF

    IF (nitrogenMode) THEN
       !! lower down the nitrogen density of above ground litter pools
       Npool_litter_green_ag(:,:) = Npool_litter_green_ag(:,:) * (1._dp - frac_burned(:,:))
       Npool_litter_wood_ag (:,:) = Npool_litter_wood_ag (:,:) * (1._dp - frac_burned(:,:))

       !! transfer nitrogen from the green pool to the below ground green litter pool
       Npool_litter_green_bg(:,:) = Npool_litter_green_bg(:,:)     &
                                    + Npool_green(:,:) * (1._dp - frac_green_aboveGround) * frac_burned(:,:)

       !! transfer nitrogen from the wood pool to the above ground woody litter pool
       Npool_litter_wood_ag(:,:) =  Npool_litter_wood_ag(:,:)      &
                                    + Cpool_woods(:,:) * (1._dp - frac_wood_2_atmos_fire) * frac_wood_aboveGround  &
                                                       *  1._dp/cn_litter_wood * frac_burned(:,:)

       !! transfer nitrogen from the wood pool to the below ground woody litter pool
       Npool_litter_wood_bg(:,:) =  Npool_litter_wood_bg(:,:)      &
                                    + Cpool_woods(:,:) * (1._dp - frac_wood_aboveGround) &
                                                       *  1._dp/cn_litter_wood * frac_burned(:,:)

       !! transfer nitrogen to the mineral nitrogen pool  
       sminN_pool(:,:)  =  sminN_pool(:,:) + (  Npool_mobile(:,:) * (1._dp - frac_green_aboveGround)                         &
                                              + Cpool_woods(:,:)  * (1._dp - frac_wood_2_atmos_fire * frac_wood_aboveGround) &
                                                                  * (1._dp/cn_woods - 1._dp/cn_litter_wood)                  &
                                             ) * frac_burned(:,:)
    END IF

    !! lower down the carbon density of living plant pools
    Cpool_green(:,:) = Cpool_green(:,:) * (1._dp - frac_burned(:,:))
    Cpool_woods(:,:) = Cpool_woods(:,:) * (1._dp - frac_burned(:,:))
    Cpool_reserve(:,:) = Cpool_reserve(:,:) * (1._dp - frac_burned(:,:))
     
    IF (nitrogenMode) THEN
       !! lower down the nitrogen density of living plant pools
       Npool_green(:,:)  = Npool_green(:,:)  * (1._dp - frac_burned(:,:))
       Npool_woods(:,:)  = Npool_woods(:,:)  * (1._dp - frac_burned(:,:))
       Npool_mobile(:,:) = Npool_mobile(:,:) * (1._dp - frac_burned(:,:))
    END IF

    if (.not. ldiag) deallocate(carbon_2_atmos_tiled_local)

  end subroutine relocate_carbon_fire

  ! --- relocate_carbon_damage() ---------------------------------------------------------------------------------------------------
  !
  ! This routine models the consequences of damages to the vegetation (e.g. wind break) for the carbon pools. More precisely:
  ! It is assumed that for the damaged area the carbon from the living plant pools (Cpool_green, Cpool_reserve and Cpool_woods)
  ! is partly relocated into the litter pools.
  !
  subroutine relocate_carbon_damage(nidx, ntiles, with_yasso, frac_damaged, cf,             &
                                         veg_fract_correction,                              &
                                         Cpool_green, Cpool_woods, Cpool_reserve,           &
                                         carbon_2_GreenLitterPools,                         &
                                         carbon_2_WoodLitterPools,                          &            
                                         Cpool_litter_green_ag,Cpool_litter_green_bg,       &
                                         Cpool_litter_wood_ag, Cpool_litter_wood_bg,        &
                                         ! Yasso pools
                                         YCpool_acid_ag1,      YCpool_acid_ag2,             &
                                         YCpool_water_ag1,     YCpool_water_ag2,            &
                                         YCpool_ethanol_ag1,   YCpool_ethanol_ag2,          &
                                         YCpool_nonsoluble_ag1,YCpool_nonsoluble_ag2,       &
                                         YCpool_acid_bg1,      YCpool_acid_bg2,             &
                                         YCpool_water_bg1,     YCpool_water_bg2,            &
                                         YCpool_ethanol_bg1,   YCpool_ethanol_bg2,          &
                                         YCpool_nonsoluble_bg1,YCpool_nonsoluble_bg2,       &
                                         YCpool_humus_1,       YCpool_humus_2,              &
                                         LeafLitcoef,          WoodLitcoef,                 &
                                         !Nitrogen pools and fluxes
                                         Npool_green, Npool_woods, Npool_mobile, SMINN_pool,&
                                         Npool_litter_green_ag,Npool_litter_green_bg,       &
                                         Npool_litter_wood_ag, Npool_litter_wood_bg,        &
                                         nitrogen_2_GreenLitterPools,                       &
                                         nitrogen_2_WoodLitterPools,                        &
                                         nitrogen_2_sminn)

    USE mo_exception,        ONLY: finish

    integer,intent(in)     :: nidx                        !! Vector length
    integer,intent(in)     :: ntiles                      !! Number of tiles
    logical,intent(in)     :: with_yasso
    real(dp),intent(in)    :: frac_damaged(:,:)           !! Fraction of the vegetated area of each tile damaged since the last call
                                                          !!    of this routine
    real(dp),intent(in)    :: cf(:,:)                     !! Cover fractions
    real(dp),intent(in)    :: veg_fract_correction(:,:)   !! Correction factor for cover fractions 1-exp(-LAI_max/2) (accounts for
                                                          !!    sparseness of vegetation)
    real(dp),intent(inout) :: Cpool_green(:,:)            !! Value of green carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_woods(:,:)            !! Value of wood carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_reserve(:,:)          !! Value of reserve carbon pool [mol(C)/m^(canopy)2]
    real(dp),intent(out)   :: carbon_2_GreenLitterPools(:)!! Amount of carbon relocated by wind damage
                                                          !! .. to the green litter pools [mol(C)/m^2(vegetated area)]
    real(dp),intent(out)   :: carbon_2_WoodLitterPools(:) !! Amount of carbon relocated by wind damage and fire
    
    real(dp),intent(inout),optional :: Cpool_litter_green_ag(:,:)  !! Value of above ground green litter carbon pool 
                                                          !!    [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: Cpool_litter_green_bg(:,:)  !! Value of below ground green litter carbon pool 
                                                          !!    [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: Cpool_litter_wood_ag(:,:)   !! Value of wood litter carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: Cpool_litter_wood_bg(:,:)   !! Value of wood litter carbon pool [mol(C)/m^2(canopy)]
    !   YASSO 
    ! Size class 1: green
    !   above ground C pools
    real(dp),intent(inout),optional :: YCpool_acid_ag1(:,:)        !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: YCpool_water_ag1(:,:)       !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: YCpool_ethanol_ag1(:,:)     !! Ethanol soluble litter pool for Yasso [mol(C)/m2(canopy)]
    real(dp),intent(inout),optional :: YCpool_nonsoluble_ag1(:,:)  !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools  
    real(dp),intent(inout),optional :: YCpool_acid_bg1(:,:)        !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: YCpool_water_bg1(:,:)       !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: YCpool_ethanol_bg1(:,:)     !! Ethanol soluble litter pool for Yasso [mol(C)/m2(canopy)]
    real(dp),intent(inout),optional :: YCpool_nonsoluble_bg1(:,:)  !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: YCpool_humus_1(:,:)         !! Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    ! Size class 2: wood 
    !   above ground C pools
    real(dp),intent(inout),optional :: YCpool_acid_ag2(:,:)        !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: YCpool_water_ag2(:,:)       !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: YCpool_ethanol_ag2(:,:)     !! Ethanol soluble litter pool for Yasso [mol(C)/m2(canopy)]
    real(dp),intent(inout),optional :: YCpool_nonsoluble_ag2(:,:)  !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools  
    real(dp),intent(inout),optional :: YCpool_acid_bg2(:,:)        !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: YCpool_water_bg2(:,:)       !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: YCpool_ethanol_bg2(:,:)     !! Ethanol soluble litter pool for Yasso [mol(C)/m2(canopy)]
    real(dp),intent(inout),optional :: YCpool_nonsoluble_bg2(:,:)  !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),intent(inout),optional :: YCpool_humus_2(:,:)         !! Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    !   parameters
    real(dp),intent(in),optional    :: LeafLitcoef(:,:,:)          !! Factor to spread non woody litter into yasso pools [ ]
    real(dp),intent(in),optional    :: WoodLitcoef(:,:,:)          !! Factor to spread woody litterfall into yasso pools [ ]
    ! nitrogen pools
    real(dp),intent(inout),optional  :: Npool_green(:,:)                !! Value of green nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional  :: Npool_woods(:,:)                !! Value of wood nitrogen pool [mol(N)/m^2(canopy)]
    real(dp),intent(inout),optional  :: Npool_mobile(:,:)               !! Value of mobile nitrogen pool [mol(N)/m^(canopy)2]
    real(dp),intent(inout),optional  :: SMINN_pool(:,:)                 !! Value of mineral nitrogen pool [mol(N)/m^(canopy)2]   
    real(dp),intent(inout),optional  :: Npool_litter_green_ag(:,:)      !! Value of above ground green litter nitrogen pool
                                                                        !!    [mol(N)/m^(canopy)2]
    real(dp),intent(inout),optional  :: Npool_litter_green_bg(:,:)      !! Value of below ground green litter nitrogen pool
                                                                        !!    [mol(N)/m^(canopy)2]
    real(dp),intent(inout),optional  :: Npool_litter_wood_ag(:,:)       !! Value of above ground wood litter nitrogen pool
                                                                        !!    [mol(N)/m^(canopy)2]
    real(dp),intent(inout),optional  :: Npool_litter_wood_bg(:,:)       !! Value of below wood litter nitrogen pool
                                                                        !!    [mol(N)/m^(canopy)2]
    ! nitrogen fluxes
    real(dp),intent(out),optional    :: nitrogen_2_GreenLitterPools(:)  !! Amount of nitrogen relocated by wind damage
                                                                        !! .. to the green litter pools [mol(N)/m^2(vegetated area)]
    real(dp),intent(out),optional    :: nitrogen_2_WoodLitterPools(:)   !! Amount of nitrogen relocated by wind damage
                                                                        !! .. to the wood litter pools [mol(N)/m^2(vegetated area)]
    real(dp),intent(out),optional    :: nitrogen_2_sminn(:)             !! Amount of nitrogen relocated by wind damage
                                                                        !! ... to the mineral N pool [mol(N)/m2(vegetated area) ]
    !! locals

    real(dp)  ::  carbon_2_WoodLitterPools_tiled(nidx,ntiles)
    real(dp)  ::  carbon_2_GreenLitterPools_tiled(nidx,ntiles)
    real(dp)  ::  N_2_WoodLitterPools_tiled(nidx,ntiles)
    real(dp)  ::  N_2_GreenLitterPools_tiled(nidx,ntiles) 
    real(dp)  ::  N_2_SMINN_tiled(nidx,ntiles)
    real(dp)  ::  Cpool_temp(nidx,ntiles)
    logical   ::  nitrogenMode
    
    nitrogenMode=.false.
    if (present(Npool_green)          .or. present(Npool_woods)           .or. &
       present(Npool_mobile)          .or. present(Npool_litter_green_ag) .or. &
       present(Npool_litter_wood_ag)  .or. present(Npool_litter_wood_bg)  .or. &
       present(Npool_litter_green_bg) .or. present(SMINN_pool) ) then
       nitrogenMode = .true.
       if (.not. (present(Npool_green)           .and. present(Npool_woods)           .and. &
                  present(Npool_mobile)          .and. present(Npool_litter_green_ag) .and. &
                  present(Npool_litter_wood_ag)  .and. present(Npool_litter_wood_bg)  .and. &
                  present(Npool_litter_green_bg) .and. present(SMINN_pool) ) )&
            CALL finish('relocate_carbon_damage()','at least one variable missing to handle nitrogen dynamics')
    end if    

    !! preparations

    carbon_2_WoodLitterPools_tiled(:,:) = 0._dp
    carbon_2_greenLitterPools_tiled(:,:) = 0._dp
    N_2_WoodLitterPools_tiled(:,:) = 0._dp
    N_2_greenLitterPools_tiled(:,:) = 0._dp  
    N_2_SMINN_tiled(:,:) = 0._dp  

    !! transfer carbon from the wood pool to the wood litter pools nd
    !! transfer carbon from the green and reserve pool to the green litter pool
    !! determine amount of carbon relocated

    carbon_2_WoodLitterPools_tiled(:,:) = Cpool_woods(:,:) &
                                          * cf(:,:) * veg_fract_correction(:,:) * frac_damaged(:,:)
    carbon_2_GreenLitterPools_tiled(:,:) = (Cpool_green(:,:) + Cpool_reserve(:,:)) &
                                          * cf(:,:) * veg_fract_correction(:,:) * frac_damaged(:,:)
    
    if (nitrogenMode) then
!dk Npool_mobile assigned to litter green pool    
       Npool_litter_wood_bg(:,:) = Npool_litter_wood_bg(:,:) &
                               + Cpool_woods(:,:) * (1._dp - frac_wood_aboveGround)* 1._dp/cn_litter_wood  &
                               * frac_damaged(:,:)
       Npool_litter_wood_ag(:,:) = Npool_litter_wood_ag(:,:) &
                               + Cpool_woods(:,:) * frac_wood_aboveGround * 1._dp/cn_litter_wood  &
                               * frac_damaged(:,:)
       Npool_litter_green_bg(:,:) = Npool_litter_green_bg(:,:) &
                                + (Npool_green(:,:) + Npool_mobile(:,:)) * (1._dp - frac_green_aboveGround) &
                                * frac_damaged(:,:)
       Npool_litter_green_ag(:,:) = Npool_litter_green_ag(:,:) &
                                + (Npool_green(:,:) + Npool_mobile(:,:)) * frac_green_aboveGround &
                                * frac_damaged(:,:)
       SMINN_pool(:,:) = SMINN_pool(:,:) &
                               + Cpool_woods(:,:) * (1._dp/cn_woods - 1._dp/cn_litter_wood)     &
                               * frac_damaged(:,:)
       N_2_WoodLitterPools_tiled(:,:) = Cpool_woods(:,:) * 1._dp/cn_litter_wood &
                                           * cf(:,:) * veg_fract_correction(:,:) * frac_damaged(:,:)
       N_2_GreenLitterPools_tiled(:,:) = (Npool_green(:,:) + Npool_mobile(:,:)) &
                                           * cf(:,:) * veg_fract_correction(:,:) * frac_damaged(:,:)
       N_2_sminn_tiled(:,:) = Cpool_woods(:,:) * (1._dp/cn_woods - 1._dp/cn_litter_wood) * frac_damaged(:,:)
    end if

    IF (.NOT. with_yasso) THEN
       Cpool_litter_wood_bg(:,:) = Cpool_litter_wood_bg(:,:) &
                                 + Cpool_woods(:,:) * (1._dp - frac_wood_aboveGround)*frac_damaged(:,:)
       Cpool_litter_wood_ag(:,:) = Cpool_litter_wood_ag(:,:) &
                                 + Cpool_woods(:,:) * frac_wood_aboveGround * frac_damaged(:,:)
       Cpool_litter_green_bg(:,:) = Cpool_litter_green_bg(:,:) &
                                  + (Cpool_green(:,:) + Cpool_reserve(:,:)) * (1._dp - frac_green_aboveGround) &
                                  * frac_damaged(:,:)
       Cpool_litter_green_ag(:,:) = Cpool_litter_green_ag(:,:) &
                                  + (Cpool_green(:,:) + Cpool_reserve(:,:)) * frac_green_aboveGround &
                                  * frac_damaged(:,:)
    ELSE

       ! Add aboveground damaged green and reserve
       Cpool_temp(:,:) = Cpool_green(:,:) + Cpool_reserve(:,:)
       YCpool_acid_ag1      (:,:) = &
       YCpool_acid_ag1      (:,:) + frac_green_aboveGround * frac_damaged(:,:) * LeafLitcoef(:,:,1) * Cpool_temp(:,:)
       YCpool_water_ag1     (:,:) = &
       YCpool_water_ag1     (:,:) + frac_green_aboveGround * frac_damaged(:,:) * LeafLitcoef(:,:,2) * Cpool_temp(:,:)
       YCpool_ethanol_ag1   (:,:) = &
       YCpool_ethanol_ag1   (:,:) + frac_green_aboveGround * frac_damaged(:,:) * LeafLitcoef(:,:,3) * Cpool_temp(:,:)
       YCpool_nonsoluble_ag1(:,:) = &
       YCpool_nonsoluble_ag1(:,:) + frac_green_aboveGround * frac_damaged(:,:) * LeafLitcoef(:,:,4) * Cpool_temp(:,:)

       ! Add belowground damaged green and reserve
       YCpool_acid_bg1      (:,:) = &
       YCpool_acid_bg1      (:,:) + (1._dp - frac_green_aboveGround) * frac_damaged(:,:) * LeafLitcoef(:,:,1) * Cpool_temp(:,:)
       YCpool_water_bg1     (:,:) = &
       YCpool_water_bg1     (:,:) + (1._dp - frac_green_aboveGround) * frac_damaged(:,:) * LeafLitcoef(:,:,2) * Cpool_temp(:,:)
       YCpool_ethanol_bg1   (:,:) = &
       YCpool_ethanol_bg1   (:,:) + (1._dp - frac_green_aboveGround) * frac_damaged(:,:) * LeafLitcoef(:,:,3) * Cpool_temp(:,:)
       YCpool_nonsoluble_bg1(:,:) = &
       YCpool_nonsoluble_bg1(:,:) + (1._dp - frac_green_aboveGround) * frac_damaged(:,:) * LeafLitcoef(:,:,4) * Cpool_temp(:,:)
       YCpool_humus_1       (:,:) = &
       YCpool_humus_1       (:,:) +                                    frac_damaged(:,:) * LeafLitcoef(:,:,5) * Cpool_temp(:,:)

       ! Add aboveground damaged wood
       YCpool_acid_ag2      (:,:) = &
       YCpool_acid_ag2      (:,:) + frac_wood_aboveGround  * frac_damaged(:,:) * WoodLitcoef(:,:,1) * Cpool_woods(:,:)
       YCpool_water_ag2     (:,:) = &
       YCpool_water_ag2     (:,:) + frac_wood_aboveGround  * frac_damaged(:,:) * WoodLitcoef(:,:,2) * Cpool_woods(:,:)
       YCpool_ethanol_ag2   (:,:) = &
       YCpool_ethanol_ag2   (:,:) + frac_wood_aboveGround  * frac_damaged(:,:) * WoodLitcoef(:,:,3) * Cpool_woods(:,:)
       YCpool_nonsoluble_ag2(:,:) = &
       YCpool_nonsoluble_ag2(:,:) + frac_wood_aboveGround  * frac_damaged(:,:) * WoodLitcoef(:,:,4) * Cpool_woods(:,:)

       ! Add belowground damaged wood
       YCpool_acid_bg2      (:,:) = &
       YCpool_acid_bg2      (:,:) + (1._dp -frac_wood_aboveGround)  * frac_damaged(:,:) * WoodLitcoef(:,:,1) * Cpool_woods(:,:)
       YCpool_water_bg2     (:,:) = &
       YCpool_water_bg2     (:,:) + (1._dp -frac_wood_aboveGround)  * frac_damaged(:,:) * WoodLitcoef(:,:,2) * Cpool_woods(:,:)
       YCpool_ethanol_bg2   (:,:) = &
       YCpool_ethanol_bg2   (:,:) + (1._dp -frac_wood_aboveGround)  * frac_damaged(:,:) * WoodLitcoef(:,:,3) * Cpool_woods(:,:)
       YCpool_nonsoluble_bg2(:,:) = &
       YCpool_nonsoluble_bg2(:,:) + (1._dp -frac_wood_aboveGround)  * frac_damaged(:,:) * WoodLitcoef(:,:,4) * Cpool_woods(:,:)
       YCpool_humus_2       (:,:) = &
       YCpool_humus_2       (:,:) +                                   frac_damaged(:,:) * WoodLitcoef(:,:,5) * Cpool_woods(:,:)

    END IF

    !! sum carbon fluxes for all tiles

    carbon_2_WoodLitterPools(:) = SUM(carbon_2_WoodLitterPools_tiled(:,:),DIM=2)
    carbon_2_GreenLitterPools(:) = SUM(carbon_2_GreenLitterPools_tiled(:,:),DIM=2)
    if (nitrogenMode) then
      nitrogen_2_WoodLitterPools(:) = SUM(N_2_WoodLitterPools_tiled(:,:),DIM=2)
      nitrogen_2_GreenLitterPools(:) = SUM(N_2_GreenLitterPools_tiled(:,:),DIM=2) 
      nitrogen_2_sminn(:) = SUM(N_2_sminn_tiled(:,:),DIM=2)   
    end if

    !! lower down the carbon density of living plant pools
    Cpool_green(:,:) = Cpool_green(:,:) * (1._dp - frac_damaged(:,:))
    Cpool_woods(:,:) = Cpool_woods(:,:) * (1._dp - frac_damaged(:,:))
    Cpool_reserve(:,:) = Cpool_reserve(:,:) * (1._dp - frac_damaged(:,:))
    
    if (nitrogenMode) then
      Npool_green(:,:) = Npool_green(:,:) * (1._dp - frac_damaged(:,:))
      Npool_woods(:,:) = Npool_woods(:,:) * (1._dp - frac_damaged(:,:))
      Npool_mobile(:,:) = Npool_mobile(:,:) * (1._dp - frac_damaged(:,:))        
    end if

  end subroutine relocate_carbon_damage

  ! --- C_relocation_from_LUtransitions() -----------------------------------------------------------------------------------------
  !
  ! This routine models the consequences of landcover transitions, i.e. here the landcover changes is given in form of a transition
  ! matrix. Otherwise the logic of this routine is quite similar to the routine relocate_carbonAndNitrogen() from above.
  !

  SUBROUTINE C_relocation_from_LUtransitions(lctlib,surface,nidx,ntiles,is_vegetation,                &
                                             cf_old, cf_new, Tile_TransMtrx, veg_fract_correction,    &
                                             frac_wood_2_atmos, frac_green_2_atmos,                   &
                                             frac_reserve_2_atmos,                                    &
                                             C_2_atmos,                                               &
                                             Cpool_green, Cpool_woods, Cpool_reserve,                 &
                                             Cpool_crop_harvest,                                      &
                                             Cpool_litter_green_ag,Cpool_litter_green_bg,             &
                                             Cpool_litter_wood_ag,Cpool_litter_wood_bg,               &
                                             Cpool_slow,                                              &
                                             !for N option
                                             Npool_green, Npool_woods, Npool_mobile, Npool_crop_harvest, &
                                             Npool_litter_green_ag, Npool_litter_green_bg,               &
                                             Npool_litter_wood_ag, Npool_litter_wood_bg,                 &
                                             Npool_slow,SMINN_pool,Nitrogen_2_atmos,frac_mobile_2_atmos, &    
                                             !yasso
                                             YCpool_acid_ag1, YCpool_acid_bg1,                          &
                                             YCpool_water_ag1, YCpool_water_bg1,                        &
                                             YCpool_ethanol_ag1, YCpool_ethanol_bg1,                    &
                                             YCpool_nonsoluble_ag1, YCpool_nonsoluble_bg1,              &
                                             YCpool_humus_1,                                            &
                                             YCpool_acid_ag2, YCpool_acid_bg2,                          &
                                             YCpool_water_ag2, YCpool_water_bg2,                        &
                                             YCpool_ethanol_ag2, YCpool_ethanol_bg2,                    &
                                             YCpool_nonsoluble_ag2, YCpool_nonsoluble_bg2,              &
                                             YCpool_humus_2,                                            &
                                             LeafLitcoef, WoodLitcoef,                                &
                                             ! LCC fluxes after Houghton
                                             C_2_litterGreenPools,                                    &
                                             C_2_litterWoodPool_ag,C_2_litterWoodPool_bg,             &         
                                             Cpool_onSite,Cpool_paper,Cpool_construction,             &
                                             C_onSite_2_atmos,C_paper_2_atmos,C_construction_2_atmos, &
                                             C_2_onSite,C_2_paper,C_2_construction,                   &                             
                                             lcc_scheme)   

    USE mo_jsbach_lctlib,   ONLY: lctlib_type
    USE mo_land_surface,    ONLY: land_surface_type
    USE mo_exception,       ONLY: finish

    type(lctlib_type),       intent(in) :: lctlib           !! PFT-specific constants
    type(land_surface_type), intent(in) :: surface
    integer,  intent(in)   :: lcc_scheme                    !! scheme for land cover change 
    integer, intent(in)    :: nidx                          !! Number of gridpoints in vectors
    integer, intent(in)    :: ntiles                        !! Number of tiles
    logical, intent(in)    :: is_vegetation(:,:)            !! Logical mask indcating vegetated tiles
    real(dp),intent(in)    :: cf_old(:,:)                   !! Cover fractions before landcover change or vegetation dynamics [] 
    real(dp),intent(in)    :: cf_new(:,:)                   !! Cover fraction after landcover change or vegetation dynamics []
    real(dp),intent(in)    :: Tile_TransMtrx(:,:,:)         !! Transition matrix between all tiles describing the landcover change
                                                            !!    (noOfgridPoints x noOfTiles x noOfTiles)
    real(dp),intent(in)    :: veg_fract_correction(:,:)     !! Correction factor for cover fractions 1-exp(-LAI_max/2) (accounts for
                                                            !!    sparseness of vegetation)
    real(dp),intent(in)    :: frac_wood_2_atmos             !! Fraction of wood pool immediately emitted to atmosphere (e.g. by 
                                                            !!    burning vegetation down)
    real(dp),intent(in)    :: frac_green_2_atmos            !! Fraction of green pool immediately emitted to atmosphere (e.g. by 
                                                            !!    burning vegetation down)
    real(dp),intent(in)    :: frac_reserve_2_atmos          !! Fraction of reserve pool immediately emitted to atmosphere (e.g. by
                                                            !!    burning vegetation down)
    real(dp),intent(out)   :: C_2_atmos(:)                  !! Amount of carbon directly emitted to atmosphere in this timestep 
                                                            !!    [mol(C)/m^2(vegetated area)]
    real(dp),intent(inout) :: Cpool_green(:,:)              !! Value of green carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_woods(:,:)              !! Value of wood carbon pool [mol(C)/m^2(canopy)]
    real(dp),intent(inout) :: Cpool_reserve(:,:)            !! Value of reserve carbon pool [mol(C)/m^(canopy)2]
    real(dp),intent(inout) :: Cpool_crop_harvest(:,:)       !! Value of crop harvest carbon pool [mol(C)/m^(canopy)2]

    ! litter and soil C pools of old cbalance scheme
    real(dp),optional,intent(inout) :: Cpool_litter_green_ag(:,:)      !! above ground green litter carbon pool [mol(C)/m^(canopy)2]
    real(dp),optional,intent(inout) :: Cpool_litter_green_bg(:,:)      !! below ground green litter carbon pool [mol(C)/m^(canopy)2]
    real(dp),optional,intent(inout) :: Cpool_litter_wood_ag(:,:)       !! above ground wood litter carbon pool [mol(C)/m^(canopy)2]
    real(dp),optional,intent(inout) :: Cpool_litter_wood_bg(:,:)       !! below ground wood litter carbon pool [mol(C)/m^(canopy)2]
    real(dp),optional,intent(inout) :: Cpool_slow(:,:)                 !! slow soil carbon pool [mol(C)/m^2(canopy)]

    ! Variables only needed for nitrogen option
    real(dp),intent(inout),optional :: Npool_green(:,:)
    real(dp),intent(inout),optional :: Npool_woods(:,:)
    real(dp),intent(inout),optional :: Npool_mobile(:,:)
    real(dp),intent(inout),optional :: Npool_crop_harvest(:,:)
    real(dp),intent(inout),optional :: Npool_litter_green_ag(:,:)
    real(dp),intent(inout),optional :: Npool_litter_green_bg(:,:)
    real(dp),intent(inout),optional :: Npool_litter_wood_ag(:,:)
    real(dp),intent(inout),optional :: Npool_litter_wood_bg(:,:)
    real(dp),intent(inout),optional :: Npool_slow(:,:)
    real(dp),intent(inout),optional :: SMINN_pool(:,:)
    real(dp),intent(out),optional :: Nitrogen_2_atmos(:)
    real(dp),intent(in),optional  :: frac_mobile_2_atmos

    ! variables only needed with yasso
    ! size class 1: green
    !   above ground C pools  
    real(dp),optional,intent(inout) :: YCpool_acid_ag1(:,:)            !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_water_ag1(:,:)           !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_ethanol_ag1(:,:)         !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_nonsoluble_ag1(:,:)      !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools  
    real(dp),optional,intent(inout) :: YCpool_acid_bg1(:,:)            !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_water_bg1(:,:)           !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_ethanol_bg1(:,:)         !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_nonsoluble_bg1(:,:)      !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_humus_1(:,:)             !! Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    ! size class 2: wood 
    !   above ground C pools  
    real(dp),optional,intent(inout) :: YCpool_acid_ag2(:,:)            !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_water_ag2(:,:)           !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_ethanol_ag2(:,:)         !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_nonsoluble_ag2(:,:)      !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools  
    real(dp),optional,intent(inout) :: YCpool_acid_bg2(:,:)            !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_water_bg2(:,:)           !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_ethanol_bg2(:,:)         !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_nonsoluble_bg2(:,:)      !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp),optional,intent(inout) :: YCpool_humus_2(:,:)             !! Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    !   parameters 
    real(dp),optional,intent(in)    :: LeafLitcoef(:,:,:)              !! Factor to spread non woody litterfall into yasso pools [ ]
    real(dp),optional,intent(in)    :: WoodLitcoef(:,:,:)              !! Factor to spread woody litterfall into yasso pools [ ]

    ! variables only needed with the original landuse transitions scheme (lcc_scheme==1)
    real(dp),optional,intent(out)   :: C_2_litterGreenPools(:)         !! Amount of carbon relocated by land cover change from green
                                                                       !!    and reserve pool to below and above ground green litter
    real(dp),optional,intent(out)   :: C_2_litterWoodPool_ag(:)        !! Amount of carbon relocated by land cover change from wood
                                                                       !!    pool to above ground woody litter pool
                                                                       !!    [mol(C)/m^2(vegetated area)]
    real(dp),optional,intent(out)   :: C_2_litterWoodPool_bg(:)        !! Amount of carbon relocated by land cover change from wood
                                                                       !!    pool to below ground woody litter pool
                                                                       !!    [mol(C)/m^2(vegetated area)]

    ! variables only needed with the landuse transitions scheme with antropogenic pools (lcc_scheme 2 or 3)
    real(dp),optional,intent(inout) ::  Cpool_onSite(:)                !! Amount of carbon remains in anthro annual pool
                                                                        !!   [mol(C)/m^2(vegetated area)] 
    real(dp),optional,intent(inout) ::  Cpool_paper(:)                  !! Amount of carbon remains in anthro decadal pool
    real(dp),optional,intent(inout) ::  Cpool_construction(:)           !! Amount of carbon remains in anthro centinnial pool
    real(dp),optional,intent(out)   ::  C_onSite_2_atmos(:)    !! C flux from anthro annual pool to atmos [mol(C)/m^2(veg)/day] 
    real(dp),optional,intent(out)   ::  C_paper_2_atmos(:)              !! Carbon flux from green anthro annual pool to atmos
                                                                        !!   [mol(C)/m^2(vegetated area)/day] 
    real(dp),optional,intent(out)   ::  C_construction_2_atmos(:)       !! Carbon flux from woody anthro annual pool to atmos
                                                                        !!   [mol(C)/m^2(vegetated area)/day] 
    real(dp),optional,intent(out)   ::  C_2_onSite(:)                   !! C flux to anthro annual pool [mol(C)/m^2(veg area)/day]
    real(dp),optional,intent(out)   ::  C_2_construction(:)             !! Carbon flux from woody anthro annual pool to paper
                                                                        !!   [mol(C)/m^2(vegetated area)/day] 
    real(dp),optional,intent(out)   ::  C_2_paper(:)                    !! Carbon flux from woody anthro annual pool to
                                                                        !!   construction [mol(C)/m^2(vegetated area)/day]
    
    !! locals

    real(dp) :: Cpool_green_loss(1:nidx,1:ntiles)    !! Carbon loss of green pool for each tile (per corrected vegetated area)
    real(dp) :: Cpool_reserve_loss(1:nidx,1:ntiles)  !! Carbon loss of reserve pool for each tile (per corrected vegetated area)
    real(dp) :: Cpool_woods_loss(1:nidx,1:ntiles)    !! Carbon loss of wood pool for each tile (per corrected vegetated area)
!!$    real(dp) :: Cpool_litter_green_ag_loss(1:nidx,1:ntiles)
!!$    real(dp) :: Cpool_litter_wood_ag_loss(1:nidx,1:ntiles)      
!!$    real(dp) :: YCpool_acid_ag1_loss(1:nidx,1:ntiles)      
!!$    real(dp) :: YCpool_water_ag1_loss(1:nidx,1:ntiles)      
!!$    real(dp) :: YCpool_ethanol_ag1_loss(1:nidx,1:ntiles)      
!!$    real(dp) :: YCpool_nonsoluble_ag1_loss(1:nidx,1:ntiles)      
!!$    real(dp) :: YCpool_acid_ag2_loss(1:nidx,1:ntiles)      
!!$    real(dp) :: YCpool_water_ag2_loss(1:nidx,1:ntiles)      
!!$    real(dp) :: YCpool_ethanol_ag2_loss(1:nidx,1:ntiles)      
!!$    real(dp) :: YCpool_nonsoluble_ag2_loss(1:nidx,1:ntiles)      
    real(dp) :: scale_fac(1:nidx,1:ntiles)

    integer  :: n, i
    real(dp) :: LitterGreen(1:nidx,1:ntiles)    !! Carbon gain of green litter pools (below and above ground)    
    real(dp) :: LitterWood(1:nidx,1:ntiles)     !! Carbon gain of wood litter pools (below and above ground)
    real(dp) :: Npool_green_loss(1:nidx,1:ntiles) 
    real(dp) :: Npool_wood_loss(1:nidx,1:ntiles) 
    real(dp) :: Npool_mobile_loss(1:nidx,1:ntiles)
!!$    real(dp) :: N_litter_green_ag_loss(1:nidx,1:ntiles)
!!$    real(dp) :: N_litter_wood_ag_loss(1:nidx,1:ntiles)
    real(dp) :: LitterGreenN(1:nidx,1:ntiles)
    real(dp) :: LitterWoodN(1:nidx,1:ntiles)
    real(dp) :: SMINNgain(1:nidx,1:ntiles)
   
    logical  :: with_yasso
    logical :: nitrogenMode

    !! --- GO! ---------------------------------------

    !! Check for nitrogen mode
    nitrogenMode = .false.
    if(present(Npool_green)           .or. present(Npool_woods)           .or. &
       present(Npool_mobile)          .or. present(Npool_litter_green_ag) .or. &
       present(Npool_litter_wood_ag)  .or. present(Npool_litter_wood_bg)  .or. &
       present(Npool_litter_green_bg) .or. present(Npool_slow)            .or. &
       present(SMINN_pool)            .or. present(Npool_crop_harvest) ) then
       nitrogenMode = .true.
       if (.not. (present(Npool_green)           .and. present(Npool_woods)           .and. &
                  present(Npool_mobile)          .and. present(Npool_litter_green_ag) .and. &
                  present(Npool_litter_wood_ag)  .and. present(Npool_litter_wood_bg)  .and. &
                  present(Npool_litter_green_bg) .and. present(Npool_slow)            .and. &
                  present(SMINN_pool)            .and. present(Npool_crop_harvest) ) ) &
              CALL finish('C_relocation_from_LUtransitions()','at least one variable missing to handle nitrogen dynamics') 
    end if

    !! Check if yasso or cbalance litter and soil pools should be handled
    with_yasso  = .FALSE.  
    IF (PRESENT(YCpool_acid_ag1)        .OR. PRESENT(YCpool_acid_bg1)            .OR. &
        PRESENT(YCpool_water_ag1)       .OR. PRESENT(YCpool_water_bg1)           .OR. &
        PRESENT(YCpool_ethanol_ag1)     .OR. PRESENT(YCpool_ethanol_bg1)         .OR. &
        PRESENT(YCpool_nonsoluble_ag1)  .OR. PRESENT(YCpool_nonsoluble_bg1)      .OR. &
        PRESENT(YCpool_humus_1)         .OR.                                          &
        PRESENT(YCpool_acid_ag2)        .OR. PRESENT(YCpool_acid_bg2)            .OR. &
        PRESENT(YCpool_water_ag2)       .OR. PRESENT(YCpool_water_bg2)           .OR. &
        PRESENT(YCpool_ethanol_ag2)     .OR. PRESENT(YCpool_ethanol_bg2)         .OR. &
        PRESENT(YCpool_nonsoluble_ag2)  .OR. PRESENT(YCpool_nonsoluble_bg2)      .OR. &
        PRESENT(YCpool_humus_2)         .OR. PRESENT(LeafLitcoef)                .OR. &
         PRESENT(WoodLitcoef)) THEN
        with_yasso=.TRUE.
        IF (.NOT. (PRESENT(YCpool_acid_ag1)        .AND. PRESENT(YCpool_acid_bg1)        .AND. &
                   PRESENT(YCpool_water_ag1)       .AND. PRESENT(YCpool_water_bg1)       .AND. &
                   PRESENT(YCpool_ethanol_ag1)     .AND. PRESENT(YCpool_ethanol_bg1)     .AND. &
                   PRESENT(YCpool_nonsoluble_ag1)  .AND. PRESENT(YCpool_nonsoluble_bg1)  .AND. &
                   PRESENT(YCpool_humus_1)         .AND.                                       &
                   PRESENT(YCpool_acid_ag2)        .AND. PRESENT(YCpool_acid_bg2)        .AND. &
                   PRESENT(YCpool_water_ag2)       .AND. PRESENT(YCpool_water_bg2)       .AND. &
                   PRESENT(YCpool_ethanol_ag2)     .AND. PRESENT(YCpool_ethanol_bg2)     .AND. &
                   PRESENT(YCpool_nonsoluble_ag2)  .AND. PRESENT(YCpool_nonsoluble_bg2)  .AND. &
                   PRESENT(YCpool_humus_2)          .AND. PRESENT(LeafLitcoef)           .AND. &
                   PRESENT(WoodLitcoef) ) &
            ) THEN
            CALL finish('C_relocation_from_LUtransitions()','at least one variable missing to handle yasso pools')
        END IF
    END IF
    
    ! Scaling factor for new cover fractions where tiles are shrinking
    FORALL(i=1:ntiles)
       scale_fac(:,i) = veg_fract_correction(:,i) * cf_old(:,i) * (1._dp-Tile_TransMtrx(:,i,i))
    END FORALL

    ! New sizes of plant pools
    ! Calculate carbon loss of all pools and update anthro pools if requested
    SELECT CASE (lcc_scheme)
    CASE(1) ! CN affected by land cover change is put in litter pools
    IF (.NOT. with_yasso) THEN
      CALL C_loss_and_update_anthro_pools(nidx, ntiles, lcc_scheme, lctlib, surface,     &
                                            Cpool_green, Cpool_woods, Cpool_reserve,     &
            Cpool_litter_green_ag        =  Cpool_litter_green_ag,                       &
            Cpool_litter_wood_ag         =  Cpool_litter_wood_ag,                        &
            scale_fac                    =  scale_fac,                                   &
            C_green_loss                 =  Cpool_green_loss,                            & 
            C_wood_loss                  =  Cpool_woods_loss,                            &
            C_reserve_loss               =  Cpool_reserve_loss                           &
!!$ TR: currently above ground litter pools are not affected by lcc
!!$            C_litter_green_ag_loss       =  Cpool_litter_green_ag_loss,                  &
!!$            C_litter_wood_ag_loss        =  Cpool_litter_wood_ag_loss                    &
            )
    ELSE  !! with_yasso
      CALL C_loss_and_update_anthro_pools(nidx, ntiles, lcc_scheme, lctlib, surface,     &
                                            Cpool_green, Cpool_woods, Cpool_reserve,     &
            YCpool_acid_ag1              =  YCpool_acid_ag1,                             &
            YCpool_water_ag1             =  YCpool_water_ag1,                            &
            YCpool_ethanol_ag1           =  YCpool_ethanol_ag1,                          &
            YCpool_nonsoluble_ag1        =  YCpool_nonsoluble_ag1,                       &
            YCpool_acid_ag2              =  YCpool_acid_ag2,                             &
            YCpool_water_ag2             =  YCpool_water_ag2,                            &
            YCpool_ethanol_ag2           =  YCpool_ethanol_ag2,                          &
            YCpool_nonsoluble_ag2        =  YCpool_nonsoluble_ag2,                       &
            scale_fac                    =  scale_fac,                                   &
            C_green_loss                 =  Cpool_green_loss,                            &
            C_wood_loss                  =  Cpool_woods_loss,                            &
            C_reserve_loss               =  Cpool_reserve_loss                           &
!!$ TR: currently above ground litter pools are not affected by lcc
!!$            C_acid_ag1_loss               =  YCpool_acid_ag1_loss,                       &
!!$            C_water_ag1_loss              =  YCpool_water_ag1_loss,                      &
!!$            C_ethanol_ag1_loss            =  YCpool_ethanol_ag1_loss,                    &
!!$            C_nonsoluble_ag1_loss         =  YCpool_nonsoluble_ag1_loss,                 &
!!$            C_acid_ag2_loss               =  YCpool_acid_ag2_loss,                       &
!!$            C_water_ag2_loss              =  YCpool_water_ag2_loss,                      &
!!$            C_ethanol_ag2_loss            =  YCpool_ethanol_ag2_loss,                    &
!!$            C_nonsoluble_ag2_loss         =  YCpool_nonsoluble_ag2_loss                  &
            )
    END IF  ! with_yasso
    if (nitrogenMode) then
      CALL N_loss_lcc(nidx,ntiles, lcc_scheme,                &
                      Npool_green, Npool_woods, Npool_mobile, &
                      scale_fac,                              &
                      Npool_green_loss  = Npool_green_loss,   &
                      Npool_wood_loss   = Npool_wood_loss,    &
                      Npool_mobile_loss = Npool_mobile_loss)
    end if ! nitrogenMode

    CASE(2) ! LCC scheme Grand Slam after Houghton et al. (lcc_scheme==2)
            ! CN affected by land cover change is put in anthro pools, N is lost and is not tracked anymore
    IF (.NOT. with_yasso) THEN
      CALL C_loss_and_update_anthro_pools(nidx, ntiles, lcc_scheme, lctlib, surface,     &
                                            Cpool_green, Cpool_woods, Cpool_reserve,     &
            Cpool_litter_green_ag        =  Cpool_litter_green_ag,                       &
            Cpool_litter_wood_ag         =  Cpool_litter_wood_ag,                        &
            Cpool_onSite                 =  Cpool_onSite,                                &
            Cpool_paper                  =  Cpool_paper,                                 &
            Cpool_construction           =  Cpool_construction,                          &
            scale_fac                    =  scale_fac,                                   &
            C_2_onSite                   =  C_2_onSite,                                  &
            C_2_paper                    =  C_2_paper,                                   &
            C_2_construction             =  C_2_construction,                            &
            C_onSite_2_atmos             =  C_onSite_2_atmos,                            &
            C_paper_2_atmos              =  C_paper_2_atmos,                             &
            C_construction_2_atmos       =  C_construction_2_atmos,                      &
            C_2_atmos                    =  C_2_atmos,                                   &
            C_green_loss                 =  Cpool_green_loss,                            & 
            C_wood_loss                  =  Cpool_woods_loss,                            &
            C_reserve_loss               =  Cpool_reserve_loss                           &
!!$ TR: currently above ground litter pools are not affected by lcc
!!$            C_litter_green_ag_loss       =  Cpool_litter_green_ag_loss,                  &
!!$            C_litter_wood_ag_loss        =  Cpool_litter_wood_ag_loss                    &
            )
    ELSE  !! with_yasso
      CALL C_loss_and_update_anthro_pools(nidx, ntiles, lcc_scheme, lctlib, surface,     &
                                            Cpool_green, Cpool_woods, Cpool_reserve,     &
            YCpool_acid_ag1              =  YCpool_acid_ag1,                             &
            YCpool_water_ag1             =  YCpool_water_ag1,                            &
            YCpool_ethanol_ag1           =  YCpool_ethanol_ag1,                          &
            YCpool_nonsoluble_ag1        =  YCpool_nonsoluble_ag1,                       &
            YCpool_acid_ag2              =  YCpool_acid_ag2,                             &
            YCpool_water_ag2             =  YCpool_water_ag2,                            &
            YCpool_ethanol_ag2           =  YCpool_ethanol_ag2,                          &
            YCpool_nonsoluble_ag2        =  YCpool_nonsoluble_ag2,                       &
            Cpool_onSite                 =  Cpool_onSite,                                &
            Cpool_paper                  =  Cpool_paper,                                 &
            Cpool_construction           =  Cpool_construction,                          &
            scale_fac                    =  scale_fac,                                   &
            C_2_onSite                   =  C_2_onSite,                                  &
            C_2_paper                    =  C_2_paper,                                   &
            C_2_construction             =  C_2_construction,                            &
            C_onSite_2_atmos             =  C_onSite_2_atmos,                            &
            C_paper_2_atmos              =  C_paper_2_atmos,                             &
            C_construction_2_atmos       =  C_construction_2_atmos,                      &
            C_2_atmos                    =  C_2_atmos,                                   &
            C_green_loss                 =  Cpool_green_loss,                            &
            C_wood_loss                  =  Cpool_woods_loss,                            &
            C_reserve_loss               =  Cpool_reserve_loss                           &
!!$ TR: currently above ground litter pools are not affected by lcc
!!$            C_acid_ag1_loss               =  YCpool_acid_ag1_loss,                       &
!!$            C_water_ag1_loss              =  YCpool_water_ag1_loss,                      &
!!$            C_ethanol_ag1_loss            =  YCpool_ethanol_ag1_loss,                    &
!!$            C_nonsoluble_ag1_loss         =  YCpool_nonsoluble_ag1_loss,                 &
!!$            C_acid_ag2_loss               =  YCpool_acid_ag2_loss,                       &
!!$            C_water_ag2_loss              =  YCpool_water_ag2_loss,                      &
!!$            C_ethanol_ag2_loss            =  YCpool_ethanol_ag2_loss,                    &
!!$            C_nonsoluble_ag2_loss         =  YCpool_nonsoluble_ag2_loss                  &
            )
    END IF  ! with_yasso
    IF (nitrogenMode) THEN
      CALL N_loss_lcc(nidx,ntiles, lcc_scheme,                &
                      Npool_green, Npool_woods, Npool_mobile, &
                      scale_fac,                              &
!!$ TR: currently above ground litter is not affected by lcc
!!$ TR                      Npool_litter_green_ag = Npool_litter_green_ag(:,:), &
!!$ TR                      Npool_litter_wood_ag = Npool_litter_wood_ag(:,:), & 
                      Npool_green_loss  = Npool_green_loss,   &
                      Npool_wood_loss   = Npool_wood_loss,    &
                      Npool_mobile_loss = Npool_mobile_loss,  &
!!$ TR                      N_litter_green_ag_loss = N_litter_green_ag_loss, &
!!$ TR                      N_litter_wood_ag_loss = N_litter_wood_ag_loss, &
                      N_2_atmos = Nitrogen_2_atmos(:))
    END IF ! nitrogenMode

    END SELECT !! lcc_scheme    

    !! Update vegetation C-pools and N-pools 
    WHERE(is_vegetation(:,:))
       Cpool_green            =   (cf_old * veg_fract_correction * Cpool_green           - Cpool_green_loss           ) &
                                / (cf_new * veg_fract_correction)
       Cpool_woods            =   (cf_old * veg_fract_correction * Cpool_woods           - Cpool_woods_loss           ) &
                                / (cf_new * veg_fract_correction)
       Cpool_reserve          =   (cf_old * veg_fract_correction * Cpool_reserve         - Cpool_reserve_loss         ) &
                                / (cf_new * veg_fract_correction)
    END WHERE

    if (nitrogenMode) then
       where(is_vegetation)
           Npool_green = ( cf_old * veg_fract_correction*Npool_green   - Npool_green_loss )   &
                                          / ( cf_new*veg_fract_correction )
           Npool_woods = ( cf_old * veg_fract_correction*Npool_woods   - Npool_wood_loss )   &
                                          / ( cf_new*veg_fract_correction )
           Npool_mobile = ( cf_old * veg_fract_correction*Npool_mobile   - Npool_mobile_loss )   &
                                          / ( cf_new*veg_fract_correction )
       end where
    end if

    IF (lcc_scheme == 2) THEN

       IF (.NOT. with_yasso) THEN
!!$          WHERE(is_vegetation(:,:))
!!$             Cpool_litter_green_ag  =   (cf_old * veg_fract_correction * Cpool_litter_green_ag - Cpool_litter_green_ag_loss ) &
!!$                                      / (cf_new * veg_fract_correction)
!!$             Cpool_litter_wood_ag   =   (cf_old * veg_fract_correction * Cpool_litter_wood_ag  - Cpool_litter_wood_ag_loss  ) &
!!$                                      / (cf_new * veg_fract_correction)
!!$          END WHERE
          FORALL(n=1:nidx, i=1:ntiles, is_vegetation(n,i))
             Cpool_litter_green_ag(n,i) = (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                &
                                               * Cpool_litter_green_ag(n,:))                                                  &
                                          ) / (cf_new(n,i) * veg_fract_correction(n,i))
             Cpool_litter_wood_ag (n,i) = (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                &
                                               * Cpool_litter_wood_ag (n,:))                                                  &
                                          ) / (cf_new(n,i) * veg_fract_correction(n,i))
             Cpool_litter_green_bg(n,i) = (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                &
                                               * Cpool_litter_green_bg(n,:))                                                  &
                                           + (Cpool_green_loss(n,i)+Cpool_reserve_loss(n,i)) * (1._dp-frac_green_aboveGround) &
                                          ) / (cf_new(n,i) * veg_fract_correction(n,i))
             Cpool_litter_wood_bg (n,i) = (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                &
                                               * Cpool_litter_wood_bg (n,:))                                                  &
                                           +  Cpool_woods_loss(n,i)                          * (1._dp-frac_wood_aboveGround)  &
                                          ) / (cf_new(n,i) * veg_fract_correction(n,i))
          END FORALL
       ELSE  !! with yasso 
!!$          WHERE(is_vegetation(:,:))
!!$             YCpool_acid_ag1       =   (cf_old * veg_fract_correction * YCpool_acid_ag1       - YCpool_acid_ag1_loss )         &
!!$                                    / (cf_new * veg_fract_correction) 
!!$             YCpool_water_ag1      =   (cf_old * veg_fract_correction * YCpool_water_ag1      - YCpool_water_ag1_loss )        &
!!$                                    / (cf_new * veg_fract_correction) 
!!$             YCpool_ethanol_ag1    =   (cf_old * veg_fract_correction * YCpool_ethanol_ag1    - YCpool_ethanol_ag1_loss )      &
!!$                                    / (cf_new * veg_fract_correction) 
!!$             YCpool_nonsoluble_ag1 =   (cf_old * veg_fract_correction * YCpool_nonsoluble_ag1 - YCpool_nonsoluble_ag1_loss )   &
!!$                                    / (cf_new * veg_fract_correction) 
!!$             YCpool_acid_ag2       =   (cf_old * veg_fract_correction * YCpool_acid_ag2       - YCpool_acid_ag2_loss )         &
!!$                                    / (cf_new * veg_fract_correction) 
!!$             YCpool_water_ag2      =   (cf_old * veg_fract_correction * YCpool_water_ag2      - YCpool_water_ag2_loss )        &
!!$                                    / (cf_new * veg_fract_correction) 
!!$             YCpool_ethanol_ag2      =   (cf_old * veg_fract_correction * YCpool_ethanol_ag2  - YCpool_ethanol_ag2_loss )      &
!!$                                    / (cf_new * veg_fract_correction) 
!!$             YCpool_nonsoluble_ag2 =   (cf_old * veg_fract_correction * YCpool_nonsoluble_ag2 - YCpool_nonsoluble_ag2_loss )   &
!!$                                    / (cf_new * veg_fract_correction) 
!!$          END WHERE
          FORALL(n=1:nidx, i=1:ntiles, is_vegetation(n,i))
           ! new pool content =  inputs due to cover fraction change scaled with new fraction
            YCpool_acid_ag1(n,i)    = (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                    &
                                               * YCpool_acid_ag1(n,:))                                                        &
                                      ) / (cf_new(n,i) * veg_fract_correction(n,i))
            YCpool_water_ag1(n,i)   = (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                    &
                                               * YCpool_water_ag1(n,:))                                                       &
                                      ) / (cf_new(n,i) * veg_fract_correction(n,i))
            YCpool_ethanol_ag1(n,i) = (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                    &
                                               * YCpool_ethanol_ag1(n,:))                                                     &
                                      ) / (cf_new(n,i) * veg_fract_correction(n,i))
            YCpool_nonsoluble_ag1(n,i) = (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                 &
                                               * YCpool_nonsoluble_ag1(n,:))                                                  &
                                      ) / (cf_new(n,i) * veg_fract_correction(n,i))
            YCpool_acid_ag2(n,i)  = (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                      &
                                               * YCpool_acid_ag2(n,:))                                                        &
                                    ) / (cf_new(n,i) * veg_fract_correction(n,i))
            YCpool_water_ag2(n,i) = (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                      &
                                               * YCpool_water_ag2(n,:))                                                       &
                                    ) / (cf_new(n,i) * veg_fract_correction(n,i))
            YCpool_ethanol_ag2(n,i)=(SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                      &
                                               * YCpool_ethanol_ag2(n,:))                                                     &
                                    ) / (cf_new(n,i) * veg_fract_correction(n,i))
            YCpool_nonsoluble_ag2(n,i)= (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                  &
                                               * YCpool_nonsoluble_ag2(n,:))                                                  &
                                    ) / (cf_new(n,i) * veg_fract_correction(n,i))
           ! new pool content =  inputs due to cover fraction change 
           !                   + inputs from vegetation pools (only belowground & acid soluble part)
           !                     scaled with new fraction
            YCpool_acid_bg1(n,i)    = (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                    &
                                               * YCpool_acid_bg1(n,:))                                                        &
                                           + (Cpool_green_loss(n,i)+Cpool_reserve_loss(n,i)) * (1._dp-frac_green_aboveGround) &
                                           * LeafLitcoef(n,i,1))/ (cf_new(n,i) * veg_fract_correction(n,i))
            YCpool_water_bg1(n,i)   = (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                    &
                                               * YCpool_water_bg1(n,:))                                                       &
                                           + (Cpool_green_loss(n,i)+Cpool_reserve_loss(n,i)) * (1._dp-frac_green_aboveGround) &
                                           * LeafLitcoef(n,i,2))/ (cf_new(n,i) * veg_fract_correction(n,i))
            YCpool_ethanol_bg1(n,i) = (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                    &
                                               * YCpool_ethanol_bg1(n,:))                                                     &
                                           + (Cpool_green_loss(n,i)+Cpool_reserve_loss(n,i)) * (1._dp-frac_green_aboveGround) &
                                           * LeafLitcoef(n,i,3))/ (cf_new(n,i) * veg_fract_correction(n,i))
            YCpool_nonsoluble_bg1(n,i) = (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                 &
                                               * YCpool_nonsoluble_bg1(n,:))                                                  &
                                           + (Cpool_green_loss(n,i)+Cpool_reserve_loss(n,i)) * (1._dp-frac_green_aboveGround) &
                                           * LeafLitcoef(n,i,4))/ (cf_new(n,i) * veg_fract_correction(n,i))
            YCpool_humus_1(n,i)   =  (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                     &
                                               * YCpool_humus_1(n,:))                                                         &
                                           + (Cpool_green_loss(n,i)+Cpool_reserve_loss(n,i)) * (1._dp-frac_green_aboveGround) &
                                           * LeafLitcoef(n,i,5))/ (cf_new(n,i) * veg_fract_correction(n,i))
            YCpool_acid_bg2(n,i)  = (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                      &
                                               * YCpool_acid_bg2(n,:))                                                        &
                                           +  Cpool_woods_loss(n,i)                          * (1._dp-frac_wood_aboveGround)  &
                                           * WoodLitcoef(n,i,1))/ (cf_new(n,i) * veg_fract_correction(n,i))
            YCpool_water_bg2(n,i) = (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                      &
                                               * YCpool_water_bg2(n,:))                                                       &
                                           +  Cpool_woods_loss(n,i)                          * (1._dp-frac_wood_aboveGround)  &
                                           * WoodLitcoef(n,i,2))/ (cf_new(n,i) * veg_fract_correction(n,i))
            YCpool_ethanol_bg2(n,i)=(SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                      &
                                               * YCpool_ethanol_bg2(n,:))                                                     &
                                           +  Cpool_woods_loss(n,i)                          * (1._dp-frac_wood_aboveGround)  &
                                           * WoodLitcoef(n,i,3))/ (cf_new(n,i) * veg_fract_correction(n,i))
            YCpool_nonsoluble_bg2(n,i)= (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                  &
                                               * YCpool_nonsoluble_bg2(n,:))                                                  &
                                           +  Cpool_woods_loss(n,i)                          * (1._dp-frac_wood_aboveGround)  &
                                           * WoodLitcoef(n,i,4))/ (cf_new(n,i) * veg_fract_correction(n,i))
            YCpool_humus_2(n,i)   =  (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)                     &
                                               * YCpool_humus_2(n,:))                                                         &
                                           +  Cpool_woods_loss(n,i)                          * (1._dp-frac_wood_aboveGround)  &
                                           * WoodLitcoef(n,i,5))/ (cf_new(n,i) * veg_fract_correction(n,i))
          END FORALL
       END IF

       IF (nitrogenMode) THEN
!!$          WHERE(is_vegetation(:,:))
!!$             Npool_litter_green_ag  =   (cf_old * veg_fract_correction * Npool_litter_green_ag - N_litter_green_ag_loss ) &
!!$                                      / (cf_new * veg_fract_correction)
!!$             Npool_litter_wood_ag   =   (cf_old * veg_fract_correction * Npool_litter_wood_ag  - N_litter_wood_ag_loss  ) &
!!$                                      / (cf_new * veg_fract_correction)
!!$          END WHERE

          FORALL(n=1:nidx, i=1:ntiles, is_vegetation(n,i))

             Npool_litter_green_ag(n,i) = (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)  &
                                               * Npool_litter_green_ag(n,:))                                    &
                                          ) / (cf_new(n,i) * veg_fract_correction(n,i))
             Npool_litter_wood_ag (n,i) = (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)  &
                                               * Npool_litter_wood_ag(n,:))                                     &
                                          ) / (cf_new(n,i) * veg_fract_correction(n,i))

             Npool_litter_green_bg(n,i) = (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)  &
                                               * Npool_litter_green_bg(n,:))                                    &
                                           +  Npool_green_loss(n,i) * (1._dp-frac_green_aboveGround)            &
                                          ) / (cf_new(n,i) * veg_fract_correction(n,i))
             !! ... Note: assuming cn_litter_wood >= cn_woods the surplus N from the transfer of wood to the woody litter pool
             !! ... is put into the soil mineral N pool
             Npool_litter_wood_bg (n,i) = (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)  &
                                               * Npool_litter_wood_bg(n,:))                                     &
                                           +  Npool_wood_loss(n,i) * cn_woods / cn_litter_wood                  &
                                           * (1._dp-frac_wood_aboveGround)                                      &
                                          ) / (cf_new(n,i) * veg_fract_correction(n,i))
             SMINN_pool(n,i) = (SUM(Tile_TransMtrx(n,i,:) * cf_old(n,:) * veg_fract_correction(n,:)             &
                                    * SMINN_pool(n,:))                                                          &
                                + (Npool_wood_loss(n,i) * (1._dp - frac_wood_aboveGround)                       &
                                * (1._dp - cn_woods / cn_litter_wood)                                           &
                                + Npool_mobile_loss(n,i) * (1._dp-frac_green_aboveGround))                      &
                               ) / (cf_new(n,i) * veg_fract_correction(n,i))
          END FORALL
       END IF ! nitrogenMode

    ELSE ! lcc_scheme=1

       !! Carbon loss to atmosphere
       C_2_atmos(:) =   frac_green_2_atmos   * SUM(Cpool_green_loss  (:,:),DIM=2) &
                      + frac_wood_2_atmos    * SUM(Cpool_woods_loss  (:,:),DIM=2) &
                      + frac_reserve_2_atmos * SUM(Cpool_reserve_loss(:,:),DIM=2)

       !!nitrogen loss to atmosphere
       IF (nitrogenMode) THEN
           Nitrogen_2_atmos(:) =  frac_green_2_atmos  * SUM(Npool_green_loss(:,:), DIM=2) &
                                + frac_wood_2_atmos   * SUM(Npool_wood_loss(:,:),  DIM=2) &
                                + frac_mobile_2_atmos * SUM(Npool_mobile_loss(:,:),DIM=2)
       END IF

       !! Carbon gains of litter pools
       LitterGreen(:,:) =  (1._dp-frac_green_2_atmos  ) * Cpool_green_loss  (:,:) &
                         + (1._dp-frac_reserve_2_atmos) * Cpool_reserve_loss(:,:)
       LitterWood(:,:)  =  (1._dp-frac_wood_2_atmos  )  * Cpool_woods_loss  (:,:)

       IF (nitrogenMode) THEN
          !! ... Note: assuming cn_litter_wood >= cn_woods the surplus N from the transfer of wood to the woody litter pool
          !! ... is put into the soil mineral N pool
          LitterGreenN(:,:) = (1.0_dp - frac_green_2_atmos) * Npool_green_loss(:,:)
          SMINNgain(:,:) = (1.0_dp - frac_mobile_2_atmos) * Npool_mobile_loss(:,:)  &
                           + (1.0_dp - frac_wood_2_atmos) * Npool_wood_loss(:,:)    &
                           * (1.0_dp - cn_woods / cn_litter_wood)
          LitterWoodN(:,:) = (1.0_dp - frac_wood_2_atmos) * Npool_wood_loss(:,:)    &
                           * cn_woods / cn_litter_wood
       END IF
    
       !! Compute output arrays for litter
       C_2_litterGreenPools (:) = SUM(LitterGreen(:,:),DIM=2)
       C_2_litterWoodPool_ag(:) = SUM(LitterWood (:,:),DIM=2) *        frac_wood_aboveGround
       C_2_litterWoodPool_bg(:) = SUM(LitterWood (:,:),DIM=2) * (1._dp-frac_wood_aboveGround)

       !! Update litter pools
       IF (.NOT. with_yasso) THEN
          FORALL(n=1:nidx, i=1:ntiles, is_vegetation(n,i))
             Cpool_litter_green_ag(n,i) = &
                  ( frac_green_aboveGround * LitterGreen(n,i) &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*Cpool_litter_green_ag(n,:) ) &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
             Cpool_litter_green_bg(n,i) = &
                  ( (1._dp-frac_green_aboveGround) * LitterGreen(n,i) &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*Cpool_litter_green_bg(n,:) ) &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
             Cpool_litter_wood_ag(n,i) = &
                  ( frac_wood_aboveGround * LitterWood(n,i) &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*Cpool_litter_wood_ag(n,:) ) &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
             Cpool_litter_wood_bg(n,i) = &
                  ( (1._dp-frac_wood_aboveGround) * LitterWood(n,i) &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*Cpool_litter_wood_bg(n,:) ) &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
          END FORALL
       ELSE ! with_yasso
          FORALL (n=1:nidx, i=1:ntiles, is_vegetation(n,i))
             YCpool_acid_ag1(n,i)       =                                                                     &  ! Aboveground
                  (  frac_green_aboveGround * LitterGreen(n,i) * LeafLitcoef(n,i,1)                          &  ! Input leaf litter
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:) *veg_fract_correction(n,:)*YCpool_acid_ag1(n,:) )      &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_water_ag1(n,i)      =                                                                                 &
                  (  frac_green_aboveGround * LitterGreen(n,i) * LeafLitcoef(n,i,2)                                      &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_water_ag1(n,:) )      &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_ethanol_ag1(n,i)    =                                                                                 &
                  (  frac_green_aboveGround * LitterGreen(n,i) * LeafLitcoef(n,i,3)                                      &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_ethanol_ag1(n,:) )    &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_nonsoluble_ag1(n,i) =                                                                                 &
                  (  frac_green_aboveGround * LitterGreen(n,i) * LeafLitcoef(n,i,4)                                      &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_nonsoluble_ag1(n,:) ) &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_acid_bg1(n,i)       =                                                                                 &
                  (  (1._dp - frac_green_aboveGround) * LitterGreen(n,i) * LeafLitcoef(n,i,1)                            &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_acid_bg1(n,:) )       &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_water_bg1(n,i)      =                                                                                 &
                  (  (1._dp - frac_green_aboveGround) * LitterGreen(n,i) * LeafLitcoef(n,i,2)                            &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_water_bg1(n,:) )      &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_ethanol_bg1(n,i)    =                                                                                 &
                  (  (1._dp - frac_green_aboveGround) * LitterGreen(n,i) * LeafLitcoef(n,i,3)                            &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_ethanol_bg1(n,:) )    &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_nonsoluble_bg1(n,i) =                                                                                 &
                  (  (1._dp - frac_green_aboveGround) * LitterGreen(n,i) * LeafLitcoef(n,i,4)                            &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_nonsoluble_bg1(n,:) ) &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_humus_1(n,i) =                                                                                        &
                  (  (1._dp - frac_green_aboveGround) * LitterGreen(n,i) * LeafLitcoef(n,i,5)                            &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_humus_1(n,:) )        &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))

            YCpool_acid_ag2(n,i)       =                                                                     &  ! Aboveground
                   ( frac_wood_aboveGround  * LitterWood(n,i)  * WoodLitcoef(n,i,1)                          &  ! Input wood litter
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:) *veg_fract_correction(n,:)*YCpool_acid_ag2(n,:) )      &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_water_ag2(n,i)      =                                                                                 &
                   ( frac_wood_aboveGround  * LitterWood(n,i)  * WoodLitcoef(n,i,2)                                      &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_water_ag2(n,:) )      &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_ethanol_ag2(n,i)    =                                                                                 &
                   ( frac_wood_aboveGround  * LitterWood(n,i)  * WoodLitcoef(n,i,3)                                      &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_ethanol_ag2(n,:) )    &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_nonsoluble_ag2(n,i) =                                                                                 &
                   ( frac_wood_aboveGround  * LitterWood(n,i)  * WoodLitcoef(n,i,4)                                      &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_nonsoluble_ag2(n,:) ) &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_acid_bg2(n,i)       =                                                                                 &
                   ( (1._dp - frac_wood_aboveGround) * LitterWood(n,i)  * WoodLitcoef(n,i,1)                             &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_acid_bg2(n,:) )       &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_water_bg2(n,i)      =                                                                                 &
                   ( (1._dp - frac_wood_aboveGround) * LitterWood(n,i)  * WoodLitcoef(n,i,2)                             &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_water_bg2(n,:) )      &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_ethanol_bg2(n,i)    =                                                                                 &
                   ( (1._dp - frac_wood_aboveGround) * LitterWood(n,i)  * WoodLitcoef(n,i,3)                             &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_ethanol_bg2(n,:) )    &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_nonsoluble_bg2(n,i) =                                                                                 &
                   ( (1._dp - frac_wood_aboveGround) * LitterWood(n,i)  * WoodLitcoef(n,i,4)                             &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_nonsoluble_bg2(n,:) ) &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
            YCpool_humus_2(n,i) =                                                                                        &
                  (  (1._dp - frac_green_aboveGround) * LitterWood(n,i) * WoodLitcoef(n,i,5)                             &
                       + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*YCpool_humus_2(n,:) )        &
                  ) / (cf_new(n,i)*veg_fract_correction(n,i))
          END FORALL
       END IF ! with_yasso

       IF (nitrogenMode) THEN
           FORALL(n=1:nidx, i=1:ntiles, is_vegetation(n,i))
              Npool_litter_green_ag(n,i) = &
               ( frac_green_aboveGround * LitterGreenN(n,i) &
                    + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*Npool_litter_green_ag(n,:) ) &
               ) / (cf_new(n,i)*veg_fract_correction(n,i))
              Npool_litter_green_bg(n,i) = &
               ( (1._dp-frac_green_aboveGround) * LitterGreenN(n,i) &
                    + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*Npool_litter_green_bg(n,:) ) &
               ) / (cf_new(n,i)*veg_fract_correction(n,i))
              Npool_litter_wood_ag(n,i) = &
               ( frac_wood_aboveGround * LitterWoodN(n,i) &
                    + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*Npool_litter_wood_ag(n,:) ) &
               ) / (cf_new(n,i)*veg_fract_correction(n,i))	    
              Npool_litter_wood_bg(n,i) = &
               ( (1._dp-frac_wood_aboveGround) * LitterWoodN(n,i) &
                    + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*Npool_litter_wood_bg(n,:) ) &
               ) / (cf_new(n,i)*veg_fract_correction(n,i))
              SMINN_pool(n,i) = ( SMINNgain(n,i)                   &
                           + SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*SMINN_pool(n,:) ) &
                           ) / (cf_new(n,i)*veg_fract_correction(n,i))	          
           END FORALL
       END IF ! nitrogenMode

    ENDIF !lcc_scheme

    !! New size of slow soil pool - same treatment for all lcc_schemes
    !! DSG: YASSO - treatment of humus is not the same for llc_schemes
    IF (.NOT. with_yasso) THEN
       FORALL(n=1:nidx, i=1:ntiles, is_vegetation(n,i))
          Cpool_slow(n,i) = SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*Cpool_slow(n,:) ) &
                             / (cf_new(n,i)*veg_fract_correction(n,i))
       END FORALL
    END IF
    FORALL(n=1:nidx, i=1:ntiles, is_vegetation(n,i))
       Cpool_crop_harvest(n,i) = SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*Cpool_crop_harvest(n,:) ) &
                          / (cf_new(n,i)*veg_fract_correction(n,i))
    END FORALL    

    IF (nitrogenMode) THEN
       FORALL(n=1:nidx, i=1:ntiles, is_vegetation(n,i))
          Npool_slow(n,i) = SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*Npool_slow(n,:) ) &
                             / (cf_new(n,i)*veg_fract_correction(n,i))
          Npool_crop_harvest(n,i) = SUM( Tile_TransMtrx(n,i,:) * cf_old(n,:)*veg_fract_correction(n,:)*Npool_crop_harvest(n,:) ) &
                             / (cf_new(n,i)*veg_fract_correction(n,i))
       END FORALL
    END IF

  END SUBROUTINE C_relocation_from_LUtransitions

  ! --- C_relocation_from_harvest() -----------------------------------------------------------------------------------------
  !
  ! This routine models the consequences of harvest for the carbon pools and nitrogen pools.
  !
  SUBROUTINE C_relocation_from_harvest(nidx,ntiles,lcc_scheme,lctlib,surface, &
                         cf,veg_fract_correction, is_woody,                   &
                         Cpool_green,Cpool_woods,Cpool_reserve,               &
                         harvest,                                             &
                         C2litter, C2atmos,                                   &
                         Cpool_litter_green_bg, Cpool_litter_wood_bg,         &
                         YCpool_acid_bg1, YCpool_water_bg1, YCpool_ethanol_bg1,  &
                         YCpool_nonsoluble_bg1, YCpool_humus_1,                  &
                         YCpool_acid_bg2, YCpool_water_bg2, YCpool_ethanol_bg2,  &
                         YCpool_nonsoluble_bg2, YCpool_humus_2,                  &
                         LeafLit_coef, WoodLit_coef,                          &
                         Cpool_onSite, Cpool_paper, Cpool_construction,       &
                         Carbon_2_onSite, Carbon_2_paper, Carbon_2_construction, &
                         Carbon_onSite_2_atmos,                               &
                         Carbon_paper_2_atmos, Carbon_construction_2_atmos,   &
                         frac_harvest_2_slash, frac_harvest_2_atmos,          &
                         Npool_green, Npool_woods, Npool_mobile,              &
                         Npool_litter_green_bg, Npool_litter_wood_bg,         &
                         SMINN_pool,                                          &
                         N2atmos, frac_mo_harv_2_atmos)

    USE mo_jsbach_lctlib,   ONLY: lctlib_type
    USE mo_land_surface,    ONLY: land_surface_type
    USE mo_cbal_parameters, ONLY: tau_onSite, tau_paper, tau_construction
    USE mo_exception,       ONLY: finish

    integer,  intent(in)                :: nidx                       !! Number of gridpoints in vectors
    integer,  intent(in)                :: ntiles                     !! Number of tiles
    integer,  intent(in)                :: lcc_scheme
    type(lctlib_type), intent(in)       :: lctlib
    type(land_surface_type), intent(in) :: surface
    real(dp), intent(in)                :: cf(:,:)                    !! Cover fractions before landcover change or dynveg [] 
    real(dp), intent(in)                :: veg_fract_correction(:,:)  !! Correction factor for cover fractions 1-exp(-LAI_max/2)
                                                                      !!   (accounts for sparseness of vegetation)
    logical,  intent(in)                :: is_woody(:,:)              !! Indicates tiles with woody (non-agricultural) vegetation
    real(dp), intent(inout)             :: Cpool_green(:,:)           !! Value of green carbon pool [mol(C)/m^2(canopy)]
    real(dp), intent(inout)             :: Cpool_woods(:,:)           !! Value of wood carbon pool [mol(C)/m^2(canopy)]
    real(dp), intent(inout)             :: Cpool_reserve(:,:)         !! Value of reserve carbon pool [mol(C)/m^2(canopy)]
    real(dp), intent(in)                :: harvest(:)                 !! Carbon harvested [mol(C)/m^2(corr. veget. area)/timestep]
    real(dp), intent(out)               :: C2litter(:)                !! Amount of harvested carbon relocated to the below ground 
                                                                      !!     green litter pools in this timestep 
                                                                      !!     [mol(C)/m^2(vegetated area)/timestep]
    real(dp), intent(out)               :: C2atmos(:)                 !! Amount of harvested C directly emitted to atmosphere in 
                                                                      !!     this timestep [mol(C)/m^2(vegetated area)/timestep]
    real(dp), intent(inout), optional   :: Cpool_litter_green_bg(:,:) !! Value of below ground green litter C pool [mol/m^2(canopy)]
    real(dp), intent(inout), optional   :: Cpool_litter_wood_bg(:,:)  !! Value of below ground wood litter C pool [mol/m^2(canopy)]
    real(dp), intent(inout), optional   :: YCpool_acid_bg1(:,:)       !! Value of below ground acid soluble C pool [mol/m^2(canopy)]
    real(dp), intent(inout), optional   :: YCpool_water_bg1(:,:)      !! Value of below ground acid soluble C pool [mol/m^2(canopy)]
    real(dp), intent(inout), optional   :: YCpool_ethanol_bg1(:,:)    !! Value of below ground acid soluble C pool [mol/m^2(canopy)]
    real(dp), intent(inout), optional   :: YCpool_nonsoluble_bg1(:,:) !! Value of below ground acid soluble C pool [mol/m^2(canopy)]
    real(dp), intent(inout), optional   :: YCpool_humus_1(:,:)        !! Value of below ground acid soluble C pool [mol/m^2(canopy)]
    real(dp), intent(inout), optional   :: YCpool_acid_bg2(:,:)       !! Value of below ground acid soluble C pool [mol/m^2(canopy)]
    real(dp), intent(inout), optional   :: YCpool_water_bg2(:,:)      !! Value of below ground acid soluble C pool [mol/m^2(canopy)]
    real(dp), intent(inout), optional   :: YCpool_ethanol_bg2(:,:)    !! Value of below ground acid soluble C pool [mol/m^2(canopy)]
    real(dp), intent(inout), optional   :: YCpool_nonsoluble_bg2(:,:) !! Value of below ground acid soluble C pool [mol/m^2(canopy)]
    real(dp), intent(inout), optional   :: YCpool_humus_2(:,:)        !! Value of below ground acid soluble C pool [mol/m^2(canopy)]
    real(dp), intent(in),    optional   :: LeafLit_coef(:,:,:)        !! Factor to seperate non woody litte into EWAN pools [ ]
    real(dp), intent(in),    optional   :: WoodLit_coef(:,:,:)        !! Factor to seperate woody litte into EWAN pools     [ ]
    real(dp), intent(inout), optional   :: Cpool_onSite(:)            !! Anthropogenic pools & fluxes : ...
    real(dp), intent(inout), optional   :: Cpool_paper(:)
    real(dp), intent(inout), optional   :: Cpool_construction(:)
    real(dp), intent(out),   optional   :: Carbon_2_onSite(:)
    real(dp), intent(out),   optional   :: Carbon_2_paper(:)
    real(dp), intent(out),   optional   :: Carbon_2_construction(:)
    real(dp), intent(out),   optional   :: Carbon_onSite_2_atmos(:)
    real(dp), intent(out),   optional   :: Carbon_paper_2_atmos(:)
    real(dp), intent(out),   optional   :: Carbon_construction_2_atmos(:)
    real(dp), intent(in),    optional   :: frac_harvest_2_slash       !! Fraction of harvested carbon put to slash 
    real(dp), intent(in),    optional   :: frac_harvest_2_atmos       !! Fraction of harvested carbon immediately emitted to atmos. 
    real(dp), intent(inout), optional   :: Npool_green(:,:)           !! Npools ...
    real(dp), intent(inout), optional   :: Npool_woods(:,:)
    real(dp), intent(inout), optional   :: Npool_mobile(:,:)
    real(dp), intent(inout), optional   :: Npool_litter_green_bg(:,:)
    real(dp), intent(inout), optional   :: Npool_litter_wood_bg(:,:)
    real(dp), intent(inout), optional   :: SMINN_pool(:,:)
    real(dp), intent(out),   optional   :: N2atmos(:)    !! N lost to the atmosphere
    real(dp), intent(in),    optional   :: frac_mo_harv_2_atmos

    !! Local variables
    real(dp) :: Cpool_G_ag_corr(1:nidx,1:ntiles)     !! Above ground carbon in green pool (per corrected vegetated area)
    real(dp) :: Cpool_R_ag_corr(1:nidx,1:ntiles)     !! Above ground carbon in reserve pool (per corrected vegetated area)
    real(dp) :: Cpool_W_ag_corr(1:nidx,1:ntiles)     !! Above ground carbon in wood pool (per corrected vegetated area)

    real(dp) :: Cpool_nat_ag_corr(1:nidx)     !! Total above ground carbon from woody vegetation (per corrected vegetated area)
    real(dp) :: Harvest_G(1:nidx,1:ntiles)  !! Harvest from green pool
    real(dp) :: Harvest_R(1:nidx,1:ntiles)  !! Harvest from reserve pool
    real(dp) :: Harvest_W(1:nidx,1:ntiles)  !! Harvest from wood pools

    real(dp) :: Harvest_G_rel(1:nidx,1:ntiles)  !! fraction harvested from green pool (Harvest_G/Cpool_green)
    real(dp) :: Harvest_W_rel(1:nidx,1:ntiles)  !! fraction harvested from wood pool (Harvest_W/Cpool_woods)
    real(dp) :: Npool_green_loss(1:nidx,1:ntiles)  !! amount of green N harvested
    real(dp) :: Npool_mobile_loss(1:nidx,1:ntiles) !! amount of mobile N harvested
    real(dp) :: Npool_wood_loss(1:nidx,1:ntiles)   !! amount of wood N harvested

    real(dp) :: scale_factor(1:nidx,1:ntiles)

    integer  :: n, k, ct
    real(dp) :: Tmp
    
    logical  :: with_yasso 
    logical :: nitrogenMode

    !! --- GO! ---------------------------------------

    !! Check for nitrogen mode
    nitrogenMode=.false.
    IF (PRESENT(Npool_green)           .OR. &
        PRESENT(Npool_woods)           .OR. &
        PRESENT(Npool_mobile)          .OR. &
        PRESENT(Npool_litter_wood_bg)  .OR. &
        PRESENT(Npool_litter_green_bg) .OR. &
        PRESENT(SMINN_pool) ) THEN
        nitrogenMode=.TRUE.
        IF (.NOT. (PRESENT(Npool_green)           .AND. &
                   PRESENT(Npool_woods)           .AND. &
                   PRESENT(Npool_mobile)          .AND. &
                   PRESENT(Npool_litter_wood_bg)  .AND. &
                   PRESENT(Npool_litter_green_bg) .AND. &
                   PRESENT(SMINN_pool) ) ) THEN
            CALL finish('C_relocation_from_harvest()','at least one variable missing to handle nitrogen pools')
       END IF
    END IF

    !! Check if yasso or cbalance litter and soil pools should be handled
    with_yasso  = .FALSE.  
    IF (PRESENT(YCpool_acid_bg1)            .OR. &
        PRESENT(YCpool_water_bg1)           .OR. &
        PRESENT(YCpool_ethanol_bg1)         .OR. &
        PRESENT(YCpool_nonsoluble_bg1)      .OR. &
        PRESENT(YCpool_humus_2)             .OR. &
        PRESENT(YCpool_acid_bg2)            .OR. &
        PRESENT(YCpool_water_bg2)           .OR. &
        PRESENT(YCpool_ethanol_bg2)         .OR. &
        PRESENT(YCpool_nonsoluble_bg2)      .OR. &
        PRESENT(YCpool_humus_2)             .OR. &
        PRESENT(LeafLit_coef)               .OR. &
        PRESENT(WoodLit_coef) ) THEN
        with_yasso=.TRUE.
        IF (.NOT. (PRESENT(YCpool_acid_bg1)        .AND. &
                   PRESENT(YCpool_water_bg1)       .AND. &
                   PRESENT(YCpool_ethanol_bg1)     .AND. &
                   PRESENT(YCpool_nonsoluble_bg1)  .AND. &
                   PRESENT(YCpool_humus_1)         .AND. &
                   PRESENT(YCpool_acid_bg2)        .AND. &
                   PRESENT(YCpool_water_bg2)       .AND. &
                   PRESENT(YCpool_ethanol_bg2)     .AND. &
                   PRESENT(YCpool_nonsoluble_bg2)  .AND. &
                   PRESENT(YCpool_humus_2)         .AND. &
                   PRESENT(LeafLit_coef)           .AND. &
                   PRESENT(WoodLit_coef) ) ) THEN
            CALL finish('C_relocation_from_harvest()','at least one variable missing to handle yasso pools')
        END IF
    END IF

    !! Above ground carbon in vegetation pools
    Cpool_G_ag_corr(1:nidx,1:ntiles) = frac_green_aboveGround * veg_fract_correction(:,:) * cf(:,:) * Cpool_green  (:,:)
    Cpool_R_ag_corr(1:nidx,1:ntiles) = frac_green_aboveGround * veg_fract_correction(:,:) * cf(:,:) * Cpool_reserve(:,:)
    Cpool_W_ag_corr(1:nidx,1:ntiles) = frac_wood_aboveGround  * veg_fract_correction(:,:) * cf(:,:) * Cpool_woods  (:,:)

    !! Total above carbon in pools of woody vegetation of a grid cell
    Cpool_nat_ag_corr(1:nidx) = 0.0_dp
    Cpool_nat_ag_corr(1:nidx) = SUM( Cpool_G_ag_corr(1:nidx,1:ntiles) &
                                   + Cpool_R_ag_corr(1:nidx,1:ntiles) &
                                   + Cpool_W_ag_corr(1:nidx,1:ntiles) &
                                   , mask = is_woody(1:nidx,1:ntiles), DIM=2)

    !! Harvest from vegetation pools assuming harvest proportional to above ground carbon
    Harvest_G(1:nidx,1:ntiles) = 0.0_dp
    Harvest_R(1:nidx,1:ntiles) = 0.0_dp
    Harvest_W(1:nidx,1:ntiles) = 0.0_dp
    Harvest_G_rel(1:nidx,1:ntiles) = 0.0_dp
    Harvest_W_rel(1:nidx,1:ntiles) = 0.0_dp

    FORALL(n=1:nidx,k=1:ntiles, is_woody(n,k) .AND. Cpool_nat_ag_corr(n) > EPSILON(1.0_dp) &
                                                   .AND. Cpool_green(n,k) > EPSILON(1.0_dp))
       Harvest_G(n,k) = MIN(harvest(n),Cpool_nat_ag_corr(n))*Cpool_G_ag_corr(n,k)/Cpool_nat_ag_corr(n)
    END FORALL
    FORALL(n=1:nidx,k=1:ntiles, is_woody(n,k) .AND. Cpool_nat_ag_corr(n) > EPSILON(1.0_dp))
       Harvest_R(n,k) = MIN(harvest(n),Cpool_nat_ag_corr(n))*Cpool_R_ag_corr(n,k)/Cpool_nat_ag_corr(n)
    END FORALL
    FORALL(n=1:nidx,k=1:ntiles, is_woody(n,k) .AND. Cpool_nat_ag_corr(n) > EPSILON(1.0_dp) &
                                                   .AND. Cpool_woods(n,k) > EPSILON(1.0_dp))
       Harvest_W(n,k) = MIN(harvest(n),Cpool_nat_ag_corr(n))*Cpool_W_ag_corr(n,k)/Cpool_nat_ag_corr(n) 
    END FORALL

    WHERE (veg_fract_correction(:,:) > 0._dp)
       scale_factor(:,:) = 1.0_dp / (veg_fract_correction(:,:) * MAX(cf(:,:),fract_small))
    ELSEWHERE
       scale_factor(:,:) = 0.0_dp
    ENDWHERE

    !! relative loss of C by harvest, to be used to calculate loss of N by harvest
    FORALL(n=1:nidx,k=1:ntiles, is_woody(n,k) .AND. Cpool_nat_ag_corr(n) > EPSILON(1.0_dp) &
                                                   .AND. Cpool_green(n,k) > EPSILON(1.0_dp))
       Harvest_G_rel(n,k) = Harvest_G(n,k) * scale_factor(n,k) / Cpool_green(n,k)
    END FORALL
    FORALL(n=1:nidx,k=1:ntiles, is_woody(n,k) .AND. Cpool_nat_ag_corr(n) > EPSILON(1.0_dp) &
                                                   .AND. Cpool_woods(n,k) > EPSILON(1.0_dp))
       Harvest_W_rel(n,k) = Harvest_W(n,k) * scale_factor(n,k) / Cpool_woods(n,k)
    END FORALL

    !! Loss of nitrogen from living plant pools
    Npool_green_loss(:,:) = 0.0_dp
    Npool_mobile_loss(:,:) = 0.0_dp
    Npool_wood_loss = 0.0_dp

    IF (nitrogenMode) THEN
       FORALL(n=1:nidx,k=1:ntiles, is_woody(n,k) .AND. Cpool_nat_ag_corr(n) > EPSILON(1.0_dp) &
                                                      .AND. Cpool_green(n,k) > EPSILON(1.0_dp))
          Npool_green_loss(n,k)  = Npool_green(n,k)  * Harvest_G_rel(n,k)
          Npool_mobile_loss(n,k) = Npool_mobile(n,k) * Harvest_G_rel(n,k)
       END FORALL
       FORALL(n=1:nidx,k=1:ntiles, is_woody(n,k) .AND. Cpool_nat_ag_corr(n) > EPSILON(1.0_dp) &
                                                      .AND. Cpool_woods(n,k) > EPSILON(1.0_dp))
          Npool_wood_loss(n,k)   = Npool_woods(n,k)  * Harvest_W_rel(n,k)
       END FORALL
    END IF

    !! Depletion of living plant C-pools (note: veg_fract_correction(:,:)=1 in glacier boxes) 
    FORALL(n=1:nidx,k=1:ntiles, is_woody(n,k) .AND. Cpool_nat_ag_corr(n) > EPSILON(1.0_dp))
       Cpool_green  (n,k) = Cpool_green  (n,k) - Harvest_G(n,k) * scale_factor(n,k)
       Cpool_reserve(n,k) = Cpool_reserve(n,k) - Harvest_R(n,k) * scale_factor(n,k)
       Cpool_woods  (n,k) = Cpool_woods  (n,k) - Harvest_W(n,k) * scale_factor(n,k)
    END FORALL

    !! depletion of living plant N-pools
    IF (nitrogenMode) THEN
       FORALL(n=1:nidx,k=1:ntiles, is_woody(n,k) .AND. Cpool_nat_ag_corr(n) > EPSILON(1.0_dp))
          Npool_green(n,k)  = Npool_green(n,k)  - Npool_green_loss(n,k)
          Npool_mobile(n,k) = Npool_mobile(n,k) - Npool_mobile_loss(n,k)
          Npool_woods(n,k)  = Npool_woods(n,k)  - Npool_wood_loss(n,k)
       END FORALL
    END IF

    IF (lcc_scheme==1) THEN

       !! Increase of carbon in below ground litter pools
       IF (.NOT. with_yasso) THEN
          FORALL(n=1:nidx,k=1:ntiles, is_woody(n,k) .AND. Cpool_nat_ag_corr(n) > EPSILON(1.0_dp))
             Cpool_litter_green_bg(n,k) = &
             Cpool_litter_green_bg(n,k) + (1._dp-frac_harvest_2_atmos) * scale_factor(n,k) * (Harvest_G(n,k)+Harvest_R(n,k)) 
             Cpool_litter_wood_bg (n,k) = &
             Cpool_litter_wood_bg (n,k) + (1._dp-frac_harvest_2_atmos) * scale_factor(n,k) *  Harvest_W(n,k)
          END FORALL

       ELSE ! yasso

          FORALL(n=1:nidx,k=1:ntiles, is_woody(n,k) .AND. Cpool_nat_ag_corr(n) > EPSILON(1.0_dp))
             YCpool_acid_bg1      (n,k) =                                                                  &
             YCpool_acid_bg1      (n,k) +  (1._dp-frac_harvest_2_atmos) * scale_factor(n,k)                &
                                       * ( LeafLit_coef(n,k,1)    * (Harvest_G(n,k)+Harvest_R(n,k))  )      ! Non woody litter input
             YCpool_water_bg1     (n,k) =                                                                  &
             YCpool_water_bg1     (n,k) +  (1._dp-frac_harvest_2_atmos) * scale_factor(n,k)                &
                                       * ( LeafLit_coef(n,k,2)     * (Harvest_G(n,k)+Harvest_R(n,k)) )      ! Non woody litter input
             YCpool_ethanol_bg1   (n,k) =                                                                  &
             YCpool_ethanol_bg1   (n,k) +  (1._dp-frac_harvest_2_atmos) * scale_factor(n,k)                &
                                       * ( LeafLit_coef(n,k,3)     * (Harvest_G(n,k)+Harvest_R(n,k)) )      ! Non woody litter input
             YCpool_nonsoluble_bg1(n,k) =                                                                  &
             YCpool_nonsoluble_bg1(n,k) +  (1._dp-frac_harvest_2_atmos) * scale_factor(n,k)                &
                                       * ( LeafLit_coef(n,k,4)     * (Harvest_G(n,k)+Harvest_R(n,k)) )      ! Non woody litter input
             YCpool_humus_1        (n,k) =                                                                 &
             YCpool_humus_1        (n,k) +  (1._dp-frac_harvest_2_atmos) * scale_factor(n,k)               &
                                       * ( LeafLit_coef(n,k,5)     * (Harvest_G(n,k)+Harvest_R(n,k)) )      ! Non woody litter input


             YCpool_acid_bg2      (n,k) =                                                                  &
             YCpool_acid_bg2      (n,k) +  (1._dp-frac_harvest_2_atmos) * scale_factor(n,k)                &
                                           * WoodLit_coef(n,k,1)  * Harvest_W(n,k)                          ! Woody litter input
             YCpool_water_bg2     (n,k) =                                                                  &
             YCpool_water_bg2     (n,k) +  (1._dp-frac_harvest_2_atmos) * scale_factor(n,k)                &
                                           * WoodLit_coef(n,k,2)   * Harvest_W(n,k)                         ! Woody litter input
             YCpool_ethanol_bg2   (n,k) =                                                                  &
             YCpool_ethanol_bg2   (n,k) +  (1._dp-frac_harvest_2_atmos) * scale_factor(n,k)                &
                                           * WoodLit_coef(n,k,3)   * Harvest_W(n,k)                         ! Woody litter input
             YCpool_nonsoluble_bg2(n,k) =                                                                  &
             YCpool_nonsoluble_bg2(n,k) +  (1._dp-frac_harvest_2_atmos) * scale_factor(n,k)                &
                                           * WoodLit_coef(n,k,4)   * Harvest_W(n,k)                         ! Woody litter input
             YCpool_humus_2        (n,k) =                                                                 &
             YCpool_humus_2        (n,k) +  (1._dp-frac_harvest_2_atmos) * scale_factor(n,k)               &
                                           * WoodLit_coef(n,k,5)   * Harvest_W(n,k)                         ! Woody litter input
          END FORALL

       END IF ! yasso

       IF (nitrogenMode) THEN
          N2atmos(:) = 0.0_dp
          DO n = 1,nidx
             DO k = 1,ntiles
                IF (is_woody(n,k) .AND. Cpool_nat_ag_corr(n) > EPSILON(1.0_dp)) THEN
                   !! increase of nitrogen in below ground litter pools and SMINN pool
                   !! valid for Cbalance as well as for YASSO, as N-pool structure is the same for both soil carbon models
                   Npool_litter_green_bg(n,k) = Npool_litter_green_bg(n,k) &
                       + (1._dp - frac_harvest_2_atmos) * Npool_green_loss(n,k)
                   Npool_litter_wood_bg(n,k)  = Npool_litter_wood_bg(n,k)  &
                       + (1._dp - frac_harvest_2_atmos) * Npool_wood_loss(n,k) * (cn_woods / cn_litter_wood)
                   SMINN_pool(n,k) = SMINN_pool(n,k)                       &
                       + (1._dp - frac_harvest_2_atmos) * Npool_wood_loss(n,k) * (1.0_dp - cn_woods / cn_litter_wood) &
                       + (1._dp - frac_mo_harv_2_atmos) * Npool_mobile_loss(n,k)
                   !! sum direct nitrogen emissions to atmosphere
                   N2atmos(n) = N2atmos(n) + (veg_fract_correction(n,k) * MAX(cf(n,k),fract_small)) &
                       * (frac_harvest_2_atmos * Npool_green_loss(n,k)     &
                       +  frac_harvest_2_atmos * Npool_wood_loss(n,k)      &
                       +  frac_mo_harv_2_atmos * Npool_mobile_loss(n,k))
                END IF
             END DO
          END DO
       END IF

       !! Direct carbon emissions to atmosphere
       C2atmos(1:nidx) = frac_harvest_2_atmos * &
            SUM(Harvest_G(1:nidx,1:ntiles) + Harvest_R(1:nidx,1:ntiles) + Harvest_W(1:nidx,1:ntiles), DIM=2)

       !! Harvest into litter pools
       C2litter(1:nidx) = (1._dp-frac_harvest_2_atmos) * &
            SUM(Harvest_G(1:nidx,1:ntiles) + Harvest_R(1:nidx,1:ntiles) + Harvest_W(1:nidx,1:ntiles), DIM=2)
  
    ELSE ! lcc_scheme==2
       !! Relocate Carbon into anthropogenic pools (paper and construction)
       Carbon_2_onSite      (:) = 0.0_dp
       Carbon_2_paper       (:) = 0.0_dp
       Carbon_2_construction(:) = 0.0_dp
       C2litter             (:) = 0.0_dp
       IF (nitrogenMode) THEN
          N2atmos(:) = 0.0_dp
       END IF
       DO k = 1, ntiles
          DO n=1,nidx
             IF ( is_woody(n,k) .AND. Cpool_nat_ag_corr(n) > EPSILON(1.0_dp) ) THEN
                ct  = surface%cover_type(n,k)
                Tmp = (Harvest_G(n,k) + Harvest_R(n,k) + Harvest_W(n,k)) * (1.0_dp - frac_harvest_2_slash)
                Carbon_2_onSite(n) = Carbon_2_onSite(n) + Tmp &
                  * (1._dp - lctlib%frac_harv_C_2_paper(ct) - lctlib%frac_harv_C_2_construction(ct))
                Carbon_2_paper(n) = Carbon_2_paper(n) + Tmp * lctlib%frac_harv_C_2_paper(ct)
                Carbon_2_construction(n) = Carbon_2_construction(n) + Tmp * lctlib%frac_harv_C_2_construction(ct)
                IF (nitrogenMode) THEN
                   N2atmos(n) = N2atmos(n) + (veg_fract_correction(n,k) * MAX(cf(n,k),fract_small)) &
                      * (Npool_green_loss(n,k) + Npool_wood_loss(n,k) + Npool_mobile_loss(n,k)) &
                      * (1._dp - frac_harvest_2_slash)
                END IF
                IF (.NOT. with_yasso) THEN
                   Cpool_litter_green_bg(n,k) = Cpool_litter_green_bg(n,k) + scale_factor(n,k) &
                                                * frac_harvest_2_slash * (Harvest_G(n,k) + Harvest_R(n,k))
                   Cpool_litter_wood_bg (n,k) = Cpool_litter_wood_bg (n,k) + scale_factor(n,k) &
                                                * frac_harvest_2_slash * Harvest_W(n,k)
                ELSE ! with yasso
                   YCpool_acid_bg1      (n,k) = YCpool_acid_bg1      (n,k) + LeafLit_coef(n,k,1) * frac_harvest_2_slash &
                                               * scale_factor(n,k) * (Harvest_G(n,k) + Harvest_R(n,k))
                   YCpool_water_bg1     (n,k) = YCpool_water_bg1     (n,k) + LeafLit_coef(n,k,2) * frac_harvest_2_slash &
                                               * scale_factor(n,k) * (Harvest_G(n,k) + Harvest_R(n,k))
                   YCpool_ethanol_bg1   (n,k) = YCpool_ethanol_bg1   (n,k) + LeafLit_coef(n,k,3) * frac_harvest_2_slash &
                                               * scale_factor(n,k) * (Harvest_G(n,k) + Harvest_R(n,k))
                   YCpool_nonsoluble_bg1(n,k) = YCpool_nonsoluble_bg1(n,k) + LeafLit_coef(n,k,4) * frac_harvest_2_slash &
                                               * scale_factor(n,k) * (Harvest_G(n,k) + Harvest_R(n,k))
                   YCpool_humus_1       (n,k) = YCpool_humus_1       (n,k) + LeafLit_coef(n,k,5) * frac_harvest_2_slash &
                                               * scale_factor(n,k) * (Harvest_G(n,k) + Harvest_R(n,k))
                   YCpool_acid_bg2      (n,k) = YCpool_acid_bg2      (n,k) + WoodLit_coef(n,k,1) * frac_harvest_2_slash &
                                               * scale_factor(n,k) * Harvest_W(n,k)
                   YCpool_water_bg2     (n,k) = YCpool_water_bg2     (n,k) + WoodLit_coef(n,k,2) * frac_harvest_2_slash &
                                               * scale_factor(n,k) * Harvest_W(n,k)
                   YCpool_ethanol_bg2   (n,k) = YCpool_ethanol_bg2   (n,k) + WoodLit_coef(n,k,3) * frac_harvest_2_slash &
                                               * scale_factor(n,k) * Harvest_W(n,k)
                   YCpool_nonsoluble_bg2(n,k) = YCpool_nonsoluble_bg2(n,k) + WoodLit_coef(n,k,4) * frac_harvest_2_slash &
                                               * scale_factor(n,k) * Harvest_W(n,k)
                   YCpool_humus_2       (n,k) = YCpool_humus_2      (n,k)  + WoodLit_coef(n,k,5) * frac_harvest_2_slash &
                                               * scale_factor(n,k) * Harvest_W(n,k)
                END IF ! yasso
                C2litter                  (n) = C2litter(n) + frac_harvest_2_slash * scale_factor(n,k) &
                                               * (Harvest_G(n,k) + Harvest_R(n,k) + Harvest_W(n,k))
                IF (nitrogenMode) THEN
                   !! increase of nitrogen in below ground litter pools and SMINN pool
                   !! valid for Cbalance as well as for YASSO, as N-pool structure is the same for both soil carbon models
                   Npool_litter_green_bg(n,k) = Npool_litter_green_bg(n,k) + Npool_green_loss(n,k) * frac_harvest_2_slash
                   Npool_litter_wood_bg(n,k)  = Npool_litter_wood_bg(n,k)           &
                      + Npool_wood_loss(n,k) * (cn_woods / cn_litter_wood) * frac_harvest_2_slash
                   SMINN_pool(n,k) = SMINN_pool(n,k)                                                       &
                      + Npool_wood_loss(n,k) * (1.0_dp - cn_woods / cn_litter_wood) * frac_harvest_2_slash &
                      + Npool_mobile_loss(n,k) * frac_harvest_2_slash
                END IF ! nitrogen
             ENDIF ! woody
          ENDDO
       ENDDO

       !! Carbon from anthropogenic pools to atmosphere
       Carbon_onSite_2_atmos      (:) = Cpool_onSite      (:) / tau_onSite
       Carbon_paper_2_atmos       (:) = Cpool_paper       (:) / tau_paper
       Carbon_construction_2_atmos(:) = Cpool_construction(:) / tau_construction

       !! New sizes of anthropogenic pools
       Cpool_onSite      (:) = Cpool_onSite      (:) + Carbon_2_onSite      (:) - Carbon_onSite_2_atmos      (:)
       Cpool_paper       (:) = Cpool_paper       (:) + Carbon_2_paper       (:) - Carbon_paper_2_atmos       (:)
       Cpool_construction(:) = Cpool_construction(:) + Carbon_2_construction(:) - Carbon_construction_2_atmos(:)

       !! Total CO2 emission from harvest
       C2atmos(:) = Carbon_onSite_2_atmos(:) + Carbon_paper_2_atmos(:) + Carbon_construction_2_atmos(:)

    ENDIF ! lcc_scheme

  END SUBROUTINE C_relocation_from_harvest


  ! Distribute carbon from clear cuts according to the given fractions
  ! (implementation following the example of C_relocation_from_harvest -- where only Yasso bg pools are filled)
  SUBROUTINE C_relocation_from_forest_management_clear_cuts(nidx, ntiles,                                                        &
      & veg_ratio_max, cover_fract, veg_fract_correction, clearcutThisCellAndTile, LeafLit_coef, WoodLit_coef,                   &
      & FOM_fm_frac_woodC_to_burning_this_year, FOM_fm_frac_woodC_to_paper_this_year,                                            &
      & FOM_fm_frac_woodC_to_construction_this_year, FOM_fm_frac_woodC_to_eternity_this_year,                                    &
      & Cpool_woods, Cpool_green, Cpool_reserve, YCpool_humus_1, YCpool_acid_bg1, YCpool_water_bg1, YCpool_ethanol_bg1,          &
      & YCpool_nonsoluble_bg1, YCpool_humus_2, YCpool_acid_bg2, YCpool_water_bg2, YCpool_ethanol_bg2, YCpool_nonsoluble_bg2,     &
      & FOM_fm_Box_burning_pool, FOM_fm_Box_paper_pool, FOM_fm_Box_construction_pool, FOM_fm_Box_eternity_pool,                  &
      & FOM_fm_Box_wood_c_to_product_pools_this_year, FOM_fm_Box_green_c_harvested_this_year,                                    &
      & FOM_fm_Box_reserve_c_harvested_this_year, FOM_fm_Box_wood_c_to_litter_this_year,FOM_fm_Box_wood_c_to_atmos_this_year, & 
      & YCpool_acid_ag2, YCpool_water_ag2,  &
      & YCpool_ethanol_ag2, YCpool_nonsoluble_ag2, surface)


    USE mo_land_surface,    ONLY: land_surface_type
    USE mo_cbal_parameters, ONLY: frac_harvest_2_atmos
    USE mo_exception,       ONLY: message, message_text

    INTEGER,  INTENT(IN)    :: nidx                          !! Number of gridpoints in vectors
    INTEGER,  INTENT(IN)    :: ntiles                        !! Number of tiles
    REAL(dp), INTENT(IN)    :: veg_ratio_max(:)
    REAL(dp), INTENT(IN)    :: cover_fract(:,:)
    REAL(dp), INTENT(IN)    :: veg_fract_correction(:,:)     !! Correction factor for cover fractions 1-exp(-LAI_max/2)
                                                                      !!   (accounts for sparseness of vegetation)
    LOGICAL(dp), INTENT(IN) :: clearcutThisCellAndTile(:,:)  !! Boolean indicating if pft in cell is subject to clear cut
    REAL(dp), INTENT(IN)    :: LeafLit_coef(:,:,:)           !! Factor to seperate non woody litte into EWAN pools
    REAL(dp), INTENT(IN)    :: WoodLit_coef(:,:,:)           !! Factor to seperate woody litter into EWAN pools
    !REAL(dp), INTENT(IN)    :: frac_harvest_2_atmos          !! Fraction of harvested carbon immediately emitted to atmos. 

    REAL(dp), INTENT(IN)    :: FOM_fm_frac_woodC_to_burning_this_year(:,:)      !! Fraction of ag wood C into burning pool
    REAL(dp), INTENT(IN)    :: FOM_fm_frac_woodC_to_paper_this_year(:,:)        !! Fraction of ag wood C into paper pool
    REAL(dp), INTENT(IN)    :: FOM_fm_frac_woodC_to_construction_this_year(:,:) !! Fraction of ag wood C into construction pool
    REAL(dp), INTENT(IN)    :: FOM_fm_frac_woodC_to_eternity_this_year(:,:)     !! Fraction of ag wood C into eternity pool

    REAL(dp), INTENT(INOUT) :: Cpool_woods(:,:)           !! Value of wood carbon pool [mol(C)/m^2(canopy)]
    REAL(dp), INTENT(INOUT) :: Cpool_green(:,:)           !! Value of green carbon pool [mol(C)/m^2(canopy)]
    REAL(dp), INTENT(INOUT) :: Cpool_reserve(:,:)         !! Value of reserve carbon pool [mol(C)/m^2(canopy)]
    REAL(dp), INTENT(INOUT) :: YCpool_humus_1(:,:)
    REAL(dp), INTENT(INOUT) :: YCpool_acid_bg1(:,:)
    REAL(dp), INTENT(INOUT) :: YCpool_water_bg1(:,:)
    REAL(dp), INTENT(INOUT) :: YCpool_ethanol_bg1(:,:)
    REAL(dp), INTENT(INOUT) :: YCpool_nonsoluble_bg1(:,:)
    REAL(dp), INTENT(INOUT) :: YCpool_humus_2(:,:)
    REAL(dp), INTENT(INOUT) :: YCpool_acid_bg2(:,:)
    REAL(dp), INTENT(INOUT) :: YCpool_water_bg2(:,:)
    REAL(dp), INTENT(INOUT) :: YCpool_ethanol_bg2(:,:)
    REAL(dp), INTENT(INOUT) :: YCpool_nonsoluble_bg2(:,:)
    real(dp), intent(inout) :: YCpool_acid_ag2(:,:)             !! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp), intent(inout) :: YCpool_water_ag2(:,:)            !! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp), intent(inout) :: YCpool_ethanol_ag2(:,:)          !! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(dp), intent(inout) :: YCpool_nonsoluble_ag2(:,:)       !! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]

    REAL(dp), INTENT(INOUT) :: FOM_fm_Box_burning_pool(:)
    REAL(dp), INTENT(INOUT) :: FOM_fm_Box_paper_pool(:)
    REAL(dp), INTENT(INOUT) :: FOM_fm_Box_construction_pool(:)
    REAL(dp), INTENT(INOUT) :: FOM_fm_Box_eternity_pool(:)

    REAL(dp), INTENT(INOUT) :: FOM_fm_Box_wood_c_to_product_pools_this_year(:,:)
    REAL(dp), INTENT(INOUT) :: FOM_fm_Box_green_c_harvested_this_year(:,:)
    REAL(dp), INTENT(INOUT) :: FOM_fm_Box_reserve_c_harvested_this_year(:,:)
    REAL(dp), INTENT(INOUT) :: FOM_fm_Box_wood_c_to_litter_this_year(:,:)
    REAL(dp), INTENT(INOUT) :: FOM_fm_Box_wood_c_to_atmos_this_year(:)  !KN: this was introduced to put part of the AG wood fraction into 
                                                                        !atmosphere to mimic fire for unmanaged forest (set-up JV)
    TYPE(land_surface_type), INTENT(IN) :: surface   ! Surface parameters
    
    ! local
    INTEGER  :: n_gridcell, i_tile


    DO i_tile = 1,ntiles
      !harvest of managed forest
       
      WHERE(clearcutThisCellAndTile(1:nidx,i_tile) .AND. surface%is_managed(1:nidx,i_tile))
        !WRITE(message_text,*) 'mo_cbal_cpools: cutting this pft', i_tile
        ! All green and reserve carbon (below and above ground) go into the litter
        YCpool_acid_bg1(:,i_tile) = YCpool_acid_bg1(:,i_tile) +  LeafLit_coef(:,i_tile,1) &
          & * (Cpool_green(:,i_tile) + Cpool_reserve(:,i_tile))
        YCpool_water_bg1(:,i_tile) = YCpool_water_bg1(:,i_tile) +  LeafLit_coef(:,i_tile,2) &
          & * (Cpool_green(:,i_tile) + Cpool_reserve(:,i_tile))
        YCpool_ethanol_bg1(:,i_tile) = YCpool_ethanol_bg1(:,i_tile) +  LeafLit_coef(:,i_tile,3) &
          & * (Cpool_green(:,i_tile) + Cpool_reserve(:,i_tile))
        YCpool_nonsoluble_bg1(:,i_tile) = YCpool_nonsoluble_bg1(:,i_tile) +  LeafLit_coef(:,i_tile,4) &
          & * (Cpool_green(:,i_tile) + Cpool_reserve(:,i_tile))
        YCpool_humus_1(:,i_tile) = YCpool_humus_1(:,i_tile) +  LeafLit_coef(:,i_tile,5) &
          & * (Cpool_green(:,i_tile) + Cpool_reserve(:,i_tile))

        FOM_fm_Box_green_c_harvested_this_year(:,i_tile) = veg_ratio_max(:) &
          & * cover_fract(:,i_tile) * veg_fract_correction(:,i_tile) * Cpool_green(:,i_tile)

        FOM_fm_Box_reserve_c_harvested_this_year(:,i_tile) = veg_ratio_max(:) &
          & * cover_fract(:,i_tile) * veg_fract_correction(:,i_tile) * Cpool_reserve(:,i_tile)

        Cpool_green(:,i_tile) = 0._dp
        Cpool_reserve(:,i_tile) = 0._dp

        ! Below ground wood carbon goes into the litter
        YCpool_acid_bg2(:,i_tile) = YCpool_acid_bg2(:,i_tile) +  WoodLit_coef(:,i_tile,1) &
          & * (1._dp -frac_wood_aboveGround)* Cpool_woods(:,i_tile)
        YCpool_water_bg2(:,i_tile) = YCpool_water_bg2(:,i_tile) +  WoodLit_coef(:,i_tile,2) &
          & * (1._dp -frac_wood_aboveGround)* Cpool_woods(:,i_tile)
        YCpool_ethanol_bg2(:,i_tile) = YCpool_ethanol_bg2(:,i_tile) +  WoodLit_coef(:,i_tile,3) &
          & * (1._dp -frac_wood_aboveGround)* Cpool_woods(:,i_tile)
        YCpool_nonsoluble_bg2(:,i_tile) = YCpool_nonsoluble_bg2(:,i_tile) +  WoodLit_coef(:,i_tile,4) &
          & * (1._dp -frac_wood_aboveGround)* Cpool_woods(:,i_tile)
        YCpool_humus_2(:,i_tile) = YCpool_humus_2(:,i_tile) +  WoodLit_coef(:,i_tile,5) &
          & * (1._dp -frac_wood_aboveGround)* Cpool_woods(:,i_tile)

          
        ! A fraction of AG wood C is going to the litter to account for damage during harvest (slash fraction)
        ! Fraction is according to CMIP6 implementation (Thomas Raddatz), 30% of 130%=23%
        YCpool_acid_bg2(:,i_tile) = YCpool_acid_bg2(:,i_tile) +  WoodLit_coef(:,i_tile,1) &
          & * frac_wood_aboveGround* Cpool_woods(:,i_tile)* 0.23_dp
        YCpool_water_bg2(:,i_tile) = YCpool_water_bg2(:,i_tile) +  WoodLit_coef(:,i_tile,2) &
          & * frac_wood_aboveGround* Cpool_woods(:,i_tile)* 0.23_dp
        YCpool_ethanol_bg2(:,i_tile) = YCpool_ethanol_bg2(:,i_tile) +  WoodLit_coef(:,i_tile,3) &
          & * frac_wood_aboveGround* Cpool_woods(:,i_tile)* 0.23_dp
        YCpool_nonsoluble_bg2(:,i_tile) = YCpool_nonsoluble_bg2(:,i_tile) +  WoodLit_coef(:,i_tile,4) &
          & * frac_wood_aboveGround* Cpool_woods(:,i_tile)* 0.23_dp
        YCpool_humus_2(:,i_tile) = YCpool_humus_2(:,i_tile) +  WoodLit_coef(:,i_tile,5) &
          & * frac_wood_aboveGround* Cpool_woods(:,i_tile)* 0.23_dp

        ! AG and BG
        FOM_fm_Box_wood_c_to_litter_this_year(:,i_tile) = veg_ratio_max(:) &
          & * cover_fract(:,i_tile) * veg_fract_correction(:,i_tile) &
          & * (1._dp -frac_wood_aboveGround)* Cpool_woods(:,i_tile) + veg_ratio_max(:) &
          & * cover_fract(:,i_tile) * veg_fract_correction(:,i_tile) &
          & * frac_wood_aboveGround* Cpool_woods(:,i_tile)* 0.23_dp 

        ! The rest is distributed according to the given fractions
        FOM_fm_Box_burning_pool(:) = FOM_fm_Box_burning_pool(:) + veg_ratio_max(:) &
          & * cover_fract(:,i_tile) * veg_fract_correction(:,i_tile) &
          & * FOM_fm_frac_woodC_to_burning_this_year(:,i_tile) * frac_wood_aboveGround* Cpool_woods(:,i_tile)* 0.77_dp
        FOM_fm_Box_paper_pool(:) = FOM_fm_Box_paper_pool(:) + veg_ratio_max(:) &
          & * cover_fract(:,i_tile) * veg_fract_correction(:,i_tile) &
          & * FOM_fm_frac_woodC_to_paper_this_year(:,i_tile) * frac_wood_aboveGround* Cpool_woods(:,i_tile)* 0.77_dp
        FOM_fm_Box_construction_pool(:) = FOM_fm_Box_construction_pool(:) + veg_ratio_max(:) &
          & * cover_fract(:,i_tile) * veg_fract_correction(:,i_tile) &
          & * FOM_fm_frac_woodC_to_construction_this_year(:,i_tile) * frac_wood_aboveGround * Cpool_woods(:,i_tile)* 0.77_dp
        FOM_fm_Box_eternity_pool(:) = FOM_fm_Box_eternity_pool(:) + veg_ratio_max(:) &
          & * cover_fract(:,i_tile) * veg_fract_correction(:,i_tile) &
          & * FOM_fm_frac_woodC_to_eternity_this_year(:,i_tile) * frac_wood_aboveGround* Cpool_woods(:,i_tile)* 0.77_dp

        FOM_fm_Box_wood_c_to_product_pools_this_year(:,i_tile) = veg_ratio_max(:) &
          & * cover_fract(:,i_tile) * veg_fract_correction(:,i_tile) &
          & * frac_wood_aboveGround * Cpool_woods(:,i_tile)* 0.77_dp

        !Nothing going to atmosphere in case of managed forest
        FOM_fm_Box_wood_c_to_atmos_this_year(:)=  FOM_fm_Box_wood_c_to_atmos_this_year(:)+  0._dp

        Cpool_woods(:,i_tile) = 0._dp
      
      !stand replacing disturbance for unmanaged forest as defined in the lctlib 
      ! (used for application of JV; we have a rotation cycle for unmanaged forest to mimic disturbances)
      ELSEWHERE((.NOT. surface%is_managed(1:nidx,i_tile)) .AND. clearcutThisCellAndTile(1:nidx,i_tile))
                        

         ! All green and reserve carbon (below and above ground) go into the litter (both for fire and wind)
        YCpool_acid_bg1(:,i_tile) = YCpool_acid_bg1(:,i_tile) +  LeafLit_coef(:,i_tile,1) &
          & * (Cpool_green(:,i_tile) + Cpool_reserve(:,i_tile))
        YCpool_water_bg1(:,i_tile) = YCpool_water_bg1(:,i_tile) +  LeafLit_coef(:,i_tile,2) &
          & * (Cpool_green(:,i_tile) + Cpool_reserve(:,i_tile))
        YCpool_ethanol_bg1(:,i_tile) = YCpool_ethanol_bg1(:,i_tile) +  LeafLit_coef(:,i_tile,3) &
          & * (Cpool_green(:,i_tile) + Cpool_reserve(:,i_tile))
        YCpool_nonsoluble_bg1(:,i_tile) = YCpool_nonsoluble_bg1(:,i_tile) +  LeafLit_coef(:,i_tile,4) &
          & * (Cpool_green(:,i_tile) + Cpool_reserve(:,i_tile))
        YCpool_humus_1(:,i_tile) = YCpool_humus_1(:,i_tile) +  LeafLit_coef(:,i_tile,5) &
          & * (Cpool_green(:,i_tile) + Cpool_reserve(:,i_tile))

        FOM_fm_Box_green_c_harvested_this_year(:,i_tile) = veg_ratio_max(:) &
          & * cover_fract(:,i_tile) * veg_fract_correction(:,i_tile) * Cpool_green(:,i_tile)

        FOM_fm_Box_reserve_c_harvested_this_year(:,i_tile) = veg_ratio_max(:) &
          & * cover_fract(:,i_tile) * veg_fract_correction(:,i_tile) * Cpool_reserve(:,i_tile)

        Cpool_green(:,i_tile) = 0._dp
        Cpool_reserve(:,i_tile) = 0._dp

        ! Below ground wood carbon goes into the litter (both for fire and wind)
        YCpool_acid_bg2(:,i_tile) = YCpool_acid_bg2(:,i_tile) +  WoodLit_coef(:,i_tile,1) &
          & * (1._dp -frac_wood_aboveGround)* Cpool_woods(:,i_tile)
        YCpool_water_bg2(:,i_tile) = YCpool_water_bg2(:,i_tile) +  WoodLit_coef(:,i_tile,2) &
          & * (1._dp -frac_wood_aboveGround)* Cpool_woods(:,i_tile)
        YCpool_ethanol_bg2(:,i_tile) = YCpool_ethanol_bg2(:,i_tile) +  WoodLit_coef(:,i_tile,3) &
          & * (1._dp -frac_wood_aboveGround)* Cpool_woods(:,i_tile)
        YCpool_nonsoluble_bg2(:,i_tile) = YCpool_nonsoluble_bg2(:,i_tile) +  WoodLit_coef(:,i_tile,4) &
          & * (1._dp -frac_wood_aboveGround)* Cpool_woods(:,i_tile)
        YCpool_humus_2(:,i_tile) = YCpool_humus_2(:,i_tile) +  WoodLit_coef(:,i_tile,5) &
          & * (1._dp -frac_wood_aboveGround)* Cpool_woods(:,i_tile)

       
        !KN: Above ground wood disturbed by wind is all going to litter,
        ! a fraction of the AG wood disturbed by fire is going directly to the atmosphere
        ! the fractions to calculate the amount of AG wood that comes from wind and fire is 
        ! determined by the ratio of the global area burned to global area damaged by wind
        ! in the DECK picontrol, weighted by occurrence of forest types (=3.6). 
        ! This ratio is very specific to our application, so I hardcoded it for now...  
   
        YCpool_acid_ag2(:,i_tile) = YCpool_acid_ag2(:,i_tile) +  WoodLit_coef(:,i_tile,1) &
          & * frac_wood_aboveGround* Cpool_woods(:,i_tile)*( 1._dp/4.6_dp + (3.6_dp/4.6_dp)* (1._dp -frac_harvest_2_atmos)) 
        YCpool_water_ag2(:,i_tile) =  YCpool_water_ag2(:,i_tile) +  WoodLit_coef(:,i_tile,2) &
          & * frac_wood_aboveGround* Cpool_woods(:,i_tile)*( 1._dp/4.6_dp + (3.6_dp/4.6_dp)* (1._dp -frac_harvest_2_atmos))          
        YCpool_ethanol_ag2(:,i_tile) =   YCpool_ethanol_ag2(:,i_tile) +  WoodLit_coef(:,i_tile,3) &
          & * frac_wood_aboveGround* Cpool_woods(:,i_tile)*( 1._dp/4.6_dp + (3.6_dp/4.6_dp)* (1._dp -frac_harvest_2_atmos))       
        YCpool_nonsoluble_ag2(:,i_tile) = YCpool_nonsoluble_bg2(:,i_tile) +  WoodLit_coef(:,i_tile,4) &
          & * frac_wood_aboveGround* Cpool_woods(:,i_tile)*( 1._dp/4.6_dp + (3.6_dp/4.6_dp)* (1._dp -frac_harvest_2_atmos)) 

        !Total (ag + bg) wood carbon going to the litter 
        FOM_fm_Box_wood_c_to_litter_this_year(:,i_tile) = veg_ratio_max(:) &
          & * cover_fract(:,i_tile) * veg_fract_correction(:,i_tile) &
          & * ( frac_wood_aboveGround*Cpool_woods(:,i_tile)*( 1._dp/4.6_dp + (3.6_dp/4.6_dp)* (1._dp -frac_harvest_2_atmos))&
          & + (1._dp -frac_wood_aboveGround)* Cpool_woods(:,i_tile))
                        
        FOM_fm_Box_wood_c_to_atmos_this_year(:) = FOM_fm_Box_wood_c_to_atmos_this_year(:)  + veg_ratio_max(:) &
          & * cover_fract(:,i_tile) * veg_fract_correction(:,i_tile) &
          & * ( frac_wood_aboveGround*Cpool_woods(:,i_tile)*((3.6_dp/4.6_dp)* frac_harvest_2_atmos))
             
        Cpool_woods(:,i_tile) = 0._dp

      END WHERE
    ENDDO

  END SUBROUTINE C_relocation_from_forest_management_clear_cuts


  ! Decay of wood C in a decaying product pool (forest_management)
  ELEMENTAL SUBROUTINE C_relocation_from_product_pool_decay(this_tau, this_pool, this_pool_to_atmos_flux_per_sec)

    REAL(dp), INTENT(IN)    :: this_tau
    REAL(dp), INTENT(INOUT) :: this_pool
    REAL(dp), INTENT(INOUT) :: this_pool_to_atmos_flux_per_sec

    IF(this_tau .GT. 0._dp) THEN
      this_pool_to_atmos_flux_per_sec = this_pool / this_tau
      this_pool = this_pool - this_pool_to_atmos_flux_per_sec

      ! flux output variables are per second
      this_pool_to_atmos_flux_per_sec = this_pool_to_atmos_flux_per_sec / 86400._dp
    ENDIF

  END SUBROUTINE C_relocation_from_product_pool_decay

! Determine loss of carbon from living vegetation pools (green, reserve, wood) and above ground litter pools
! on shrinking tiles due to landcover change (maps or transitions).
! Update of anthropogenic pools: add carbon from shrinking tiles (living vegetation and above ground litter),
! substract carbon from anthro pools due to decay (this is put into the atmosphere)
! ATTENTION ! Only the anthro pools are changed here. All other carbon pools are unchanged (intent(in)),
! ATTENTION ! so that the losses determined here have to be substracted from these pools later.
! ATTENTION ! Only the above ground carbon lost from living vegetation pools (green, reserve, wood) is put in the anthro pools
! AATENTION ! the below ground carbon lost from living vegetation pools has still to be handled.
  SUBROUTINE C_loss_and_update_anthro_pools(nidx,ntiles,lcc_scheme,lctlib,surface,        &
                                            ! living biomass pools
                                            Cpool_green,Cpool_wood,Cpool_reserve,         &
                                            ! Cbalance above ground litter pools
                                            Cpool_litter_green_ag,Cpool_litter_wood_ag,   &
                                            ! YASSO above ground litter pools
                                            YCpool_acid_ag1, YCpool_water_ag1,            &
                                            YCpool_ethanol_ag1, YCpool_nonsoluble_ag1,    &
                                            YCpool_acid_ag2, YCpool_water_ag2,            &
                                            YCpool_ethanol_ag2, YCpool_nonsoluble_ag2,    &
                                            ! anthro pools
                                            Cpool_onSite,Cpool_paper,Cpool_construction,  &
                                            ! scale factor
                                            scale_fac,                                    &
                                            ! fluxes to anthro pools
                                            C_2_onSite,C_2_paper,C_2_construction,        &
                                            ! fluxes from anthro pools to atmosphere
                                            C_onSite_2_atmos,                             &
                                            C_paper_2_atmos,C_construction_2_atmos,       &
                                            C_2_atmos,                                    &
                                            ! flux of carbon lost from living biomass pools
                                            C_green_loss,                                 &
                                            C_wood_loss,                                  &
                                            C_reserve_loss                                &
!!$ TR: currently above ground litter is not affected by lcc
!!$                                            ! flux of carbon lost from Cbalance above ground litter pools 
!!$                                            C_litter_green_ag_loss,                       &
!!$                                            C_litter_wood_ag_loss,                        &
!!$                                            ! flux of carbon lost from YASSO above ground litter pools
!!$                                            C_acid_ag1_loss,                              &
!!$                                            C_water_ag1_loss,                             &
!!$                                            C_ethanol_ag1_loss,                           &
!!$                                            C_nonsoluble_ag1_loss,                        &
!!$                                            C_acid_ag2_loss,                              &
!!$                                            C_water_ag2_loss,                             &
!!$                                            C_ethanol_ag2_loss,                           &
!!$                                            C_nonsoluble_ag2_loss                         &
                                            )

  USE mo_cbal_parameters, ONLY : tau_onSite, tau_paper, tau_construction
  USE mo_jsbach_lctlib,   ONLY : lctlib_type
  USE mo_land_surface,    ONLY : land_surface_type
  USE mo_exception,       ONLY : finish

  INTEGER,                  INTENT(in)    :: nidx, ntiles     ! Grid size
  INTEGER,                  INTENT(in)    :: lcc_scheme       ! landcover change scheme number
  TYPE (lctlib_type),       INTENT(in)    :: lctlib           ! PFT specific constants
  TYPE (land_surface_type), INTENT(in)    :: surface          ! Surface cover types
  REAL(dp),                 INTENT(in)    ::               &
    Cpool_green(:,:), Cpool_reserve(:,:), Cpool_wood(:,:)     ! Carbon pools of living material
  REAL(dp), OPTIONAL,       INTENT(in)    ::               &  ! CBALANCE
    Cpool_litter_green_ag(:,:), Cpool_litter_wood_ag(:,:)     !   Above ground carbon litter pools
  REAL(dp), OPTIONAL,       INTENT(in)    ::               &  ! YASSO
    YCpool_acid_ag1(:,:), YCpool_water_ag1(:,:),           &  !   Above ground carbon litter pools
    YCpool_ethanol_ag1(:,:), YCpool_nonsoluble_ag1(:,:),   &  !   Above ground carbon litter pools
    YCpool_acid_ag2(:,:), YCpool_water_ag2(:,:),           &  !   Above ground carbon litter pools
    YCpool_ethanol_ag2(:,:), YCpool_nonsoluble_ag2(:,:)       !   Above ground carbon litter pools
  REAL(dp), OPTIONAL,       INTENT(inout) ::               &  ! Optional since optional to C_relocation_from_LUtransitions
    Cpool_onSite(:), Cpool_paper(:), Cpool_construction(:)    ! Anthropogenic carbon pools
  REAL(dp),                 INTENT(in)    ::               &
    scale_fac(:,:)                                            ! Scaling factor for C fluxes relating to lcc
  REAL(dp), OPTIONAL,       INTENT(out)   ::               &  ! Optional since optional to C_relocation_from_LUtransitions
    C_2_onSite(:), C_2_paper(:), C_2_construction(:),      &  ! Carbon floating into anthropogenic pools
    C_onSite_2_atmos(:),C_paper_2_atmos(:),C_Construction_2_atmos(:), & ! Carbon from anthropogenic pools to atmosphere
    C_2_atmos(:)                                              ! Total carbon floating from anthropogenic pools to atmosphere
! These variables are needed only for lcc specified by transition matrices
  REAL(dp), OPTIONAL,       INTENT(out)   ::               &
    C_green_loss(:,:),                                     &  ! Carbon lost from green pool
    C_wood_loss(:,:),                                      &  ! Carbon lost from wood pool
    C_reserve_loss(:,:)                                       ! Carbon lost from reserve pool
!!$ TR: currently above ground litter is not affected by lcc
!!$    C_litter_green_ag_loss(:,:),                           &  ! Carbon lost from green litter pool above ground
!!$    C_litter_wood_ag_loss(:,:),                            &  ! Carbon lost from green litter pool above ground
!!$    ! YASSO   
!!$    C_acid_ag1_loss(:,:),                                  &  ! Carbon lost from acid soluble pool above ground
!!$    C_water_ag1_loss(:,:),                                 &  ! Carbon lost from water soluble pool above ground
!!$    C_ethanol_ag1_loss(:,:),                               &  ! Carbon lost from ethanol soluble pool above ground
!!$    C_nonsoluble_ag1_loss(:,:),                            &  ! Carbon lost from non soluble pool above ground
!!$    C_acid_ag2_loss(:,:),                                  &  ! Carbon lost from acid soluble pool above ground
!!$    C_water_ag2_loss(:,:),                                 &  ! Carbon lost from water soluble pool above ground
!!$    C_ethanol_ag2_loss(:,:),                               &  ! Carbon lost from ethanol soluble pool above ground
!!$    C_nonsoluble_ag2_loss(:,:)                                ! Carbon lost from non soluble pool above ground
! locals
  INTEGER  :: i, itile, ct
  LOGICAL  :: with_yasso  
  REAL(dp) :: CGreenLoss(nidx,ntiles), CWoodLoss(nidx,ntiles), CReserveLoss(nidx,ntiles)
!!$ TR: currently above ground litter is not affected by lcc
!!$  REAL(dp) :: CLitterGreenAGLoss(nidx,ntiles), CLitterWoodAGLoss(nidx,ntiles)
!!$  REAL(dp) :: CAcidAG1Loss(nidx,ntiles),                                  &
!!$              CWaterAG1Loss(nidx,ntiles),                                 &
!!$              CEthanolAG1Loss(nidx,ntiles),                               &
!!$              CNonsolubleAG1Loss(nidx,ntiles),                            &
!!$              CAcidAG2Loss(nidx,ntiles),                                  &
!!$              CWaterAG2Loss(nidx,ntiles),                                 &
!!$              CEthanolAG2Loss(nidx,ntiles),                               &
!!$              CNonsolubleAG2Loss(nidx,ntiles)

    !! Check if yasso or cbalance litter and soil pools should be handled
    with_yasso  = .FALSE.  
    IF (PRESENT(YCpool_acid_ag1)        .OR.  &
        PRESENT(YCpool_water_ag1)       .OR.  &
        PRESENT(YCpool_ethanol_ag1)     .OR.  &
        PRESENT(YCpool_nonsoluble_ag1)  .OR.  &
        PRESENT(YCpool_acid_ag2)        .OR.  &
        PRESENT(YCpool_water_ag2)       .OR.  &
        PRESENT(YCpool_ethanol_ag2)     .OR.  &
        PRESENT(YCpool_nonsoluble_ag2)        &
        ) THEN
        with_yasso=.TRUE.
        IF (.NOT. (PRESENT(YCpool_acid_ag1)        .AND. &
                   PRESENT(YCpool_water_ag1)       .AND. &
                   PRESENT(YCpool_ethanol_ag1)     .AND. &
                   PRESENT(YCpool_nonsoluble_ag1)  .AND. &
                   PRESENT(YCpool_acid_ag2)        .AND. &
                   PRESENT(YCpool_water_ag2)       .AND. &
                   PRESENT(YCpool_ethanol_ag2)     .AND. &
                   PRESENT(YCpool_nonsoluble_ag2)        &
                   ) &
            ) CALL finish('C_loss_and_update_anthro_pools()','at least one variable missing to handle yasso pools')
    END IF

    ! calculate loss of carbon from living biomass pools due to land cover change
    CGreenLoss   (:,:) = scale_fac(:,:) * Cpool_green  (:,:) 
    CWoodLoss    (:,:) = scale_fac(:,:) * Cpool_wood   (:,:)
    CReserveLoss (:,:) = scale_fac(:,:) * Cpool_reserve(:,:)
!!$ TR: currently ag litter is not affected by lcc
!!$    ! calculate loss of carbon from above ground litter pools
!!$    IF (.NOT. with_yasso) THEN
!!$       CLitterGreenAGLoss(:,:) = scale_fac(:,:) * Cpool_litter_green_ag(:,:)  
!!$       CLitterWoodAGLoss (:,:) = scale_fac(:,:) * Cpool_litter_wood_ag (:,:)  
!!$    ELSE
!!$       CAcidAG1Loss       (:,:) = scale_fac(:,:) * YCpool_acid_ag1       (:,:) 
!!$       CWaterAG1Loss      (:,:) = scale_fac(:,:) * YCpool_water_ag1      (:,:)
!!$       CEthanolAG1Loss    (:,:) = scale_fac(:,:) * YCpool_ethanol_ag1    (:,:)
!!$       CNonsolubleAG1Loss (:,:) = scale_fac(:,:) * YCpool_nonsoluble_ag1 (:,:)
!!$       CAcidAG2Loss       (:,:) = scale_fac(:,:) * YCpool_acid_ag2       (:,:) 
!!$       CWaterAG2Loss      (:,:) = scale_fac(:,:) * YCpool_water_ag2      (:,:)
!!$       CEthanolAG2Loss    (:,:) = scale_fac(:,:) * YCpool_ethanol_ag2    (:,:)
!!$       CNonsolubleAG2Loss (:,:) = scale_fac(:,:) * YCpool_nonsoluble_ag2 (:,:)
!!$    END IF

    ! update anthro pools
    IF (lcc_scheme == 2) THEN

      C_2_onSite      (:) = 0._dp
      C_2_paper       (:) = 0._dp
      C_2_construction(:) = 0._dp

      ! Add aboveground wood loss to anthropogenic pools
      DO itile=1,ntiles
        DO i=1,nidx
          ct=surface%cover_type(i,itile)
          C_2_construction(i) = C_2_construction(i) &
                                + CWoodLoss(i,itile) * frac_wood_aboveGround * lctlib%frac_lcc_C_2_construction(ct)
          C_2_paper(i)        = C_2_paper(i) &
                                + CWoodLoss(i,itile) * frac_wood_aboveGround * lctlib%frac_lcc_C_2_paper(ct)
          C_2_onSite(i)       = C_2_onSite(i) &
                                + CWoodLoss(i,itile) * frac_wood_aboveGround * &
                                (1._dp - lctlib%frac_lcc_C_2_paper(ct) - lctlib%frac_lcc_C_2_construction(ct)) &
                                + (CGreenLoss(i,itile) + CReserveLoss(i,itile)) * frac_green_aboveGround
        ENDDO
      ENDDO

!!$ TR: currently ag litter is not put in anthro pools
!!$      ! Add litter inputs
!!$      IF (.NOT. with_yasso) THEN
!!$        DO itile=1,ntiles
!!$          DO i=1,nidx
!!$            ct=surface%cover_type(i,itile)
!!$            C_2_onSite(i) = C_2_onSite(i) &
!!$                            + CLitterGreenAGLoss(i,itile) + CLitterWoodAGLoss (i,itile)
!!$          ENDDO
!!$        ENDDO
!!$      ELSE  ! with_yasso
!!$        DO itile=1,ntiles
!!$          DO i=1,nidx
!!$            ct=surface%cover_type(i,itile)
!!$            C_2_onSite(i) = C_2_onSite(i)                                               &  
!!$                            + CAcidAG1Loss(i,itile) + CAcidAG2Loss(i,itile)             &
!!$                            + CWaterAG1Loss(i,itile) + CWaterAG2Loss(i,itile)           & 
!!$                            + CEthanolAG1Loss(i,itile) + CEthanolAG2Loss(i,itile)       & 
!!$                            + CNonsolubleAG1Loss(i,itile) + CNonsolubleAG2Loss(i,itile)
!!$          ENDDO
!!$        ENDDO
!!$      END IF

      ! Fluxes out of anthropogenic pools
      C_construction_2_atmos(:) = Cpool_construction(:) / tau_construction ! ... construction
      C_paper_2_atmos       (:) = Cpool_paper       (:) / tau_paper        ! ... paper
      C_onSite_2_atmos      (:) = Cpool_onSite      (:) / tau_onSite

      ! Carbon loss to atmosphere
      C_2_atmos(:) = C_onSite_2_atmos(:) + C_paper_2_atmos(:) + C_construction_2_atmos(:)
          
      ! New sizes of anthropogenic pools
      Cpool_onSite      (:) = Cpool_onSite      (:) + C_2_onSite      (:) - C_onSite_2_atmos      (:)
      Cpool_paper       (:) = Cpool_paper       (:) + C_2_paper       (:) - C_paper_2_atmos       (:)
      Cpool_construction(:) = Cpool_construction(:) + C_2_construction(:) - C_construction_2_atmos(:)

    ENDIF ! lcc_scheme = 2

    ! pass back the variables specifying the loss of carbon if requested by the calling routine
    IF (PRESENT(C_green_loss))           C_green_loss(:,:) = CGreenLoss(:,:)
    IF (PRESENT(C_wood_loss))            C_wood_loss(:,:) = CWoodLoss(:,:)
    IF (PRESENT(C_reserve_loss))         C_reserve_loss(:,:) = CReserveLoss(:,:)
!!$ TR: currently above ground litter is not affected by lcc
!!$    IF (PRESENT(C_litter_green_ag_loss)) C_litter_green_ag_loss(:,:) = CLitterGreenAGLoss(:,:)
!!$    IF (PRESENT(C_litter_wood_ag_loss))  C_litter_wood_ag_loss(:,:) = CLitterWoodAGLoss(:,:)
!!$    IF (PRESENT(C_acid_ag1_loss))        C_acid_ag1_loss(:,:) = CAcidAG1Loss(:,:)
!!$    IF (PRESENT(C_water_ag1_loss))       C_water_ag1_loss(:,:) = CWaterAG1Loss(:,:)
!!$    IF (PRESENT(C_ethanol_ag1_loss))     C_ethanol_ag1_loss(:,:) = CEthanolAG1Loss(:,:)
!!$    IF (PRESENT(C_nonsoluble_ag1_loss))  C_nonsoluble_ag1_loss(:,:) = CNonsolubleAG1Loss(:,:)
!!$    IF (PRESENT(C_acid_ag2_loss))        C_acid_ag2_loss(:,:) = CAcidAG2Loss(:,:)
!!$    IF (PRESENT(C_water_ag2_loss))       C_water_ag2_loss(:,:) = CWaterAG2Loss(:,:)
!!$    IF (PRESENT(C_ethanol_ag2_loss))     C_ethanol_ag2_loss(:,:) = CEthanolAG2Loss(:,:)
!!$    IF (PRESENT(C_nonsoluble_ag2_loss))  C_nonsoluble_ag2_loss(:,:) = CNonsolubleAG2Loss(:,:)

  END SUBROUTINE C_loss_and_update_anthro_pools

! Calculate loss of nitrogen from living plant pools and above ground litter pools due to landcover change on shrinking tiles.
! This routine also determines the flux of N to the atmosphere by landcover change, if the grand slam protocol is used 
! (lcc_scheme=2).
! It is analogue to routine C_loss_and_update_anthro_pools.
! ATTENTION ! All nitrogen pools are not changed in this routine (all intent(in)).
! ATTENTION ! The losses determined here must be substracted from the nitrogen pools later. 
  SUBROUTINE N_loss_lcc(nidx,ntiles, lcc_scheme,                    &
                       ! nitrogen pools
                       Npool_green, Npool_woods, Npool_mobile,      &
                       ! scale factor
                       scale_fac,                                   &
                       Npool_litter_green_ag, Npool_litter_wood_ag, &
                       ! flux of nitrogen lost from living biomass pools
                       Npool_green_loss,                            &
                       Npool_wood_loss,                             &
                       Npool_mobile_loss,                           &
                       ! flux of nitrogen lost from above ground litter pools
                       N_litter_green_ag_loss,                      &
                       N_litter_wood_ag_loss,                       &
                       ! flux of nitrogen lost to atmosphere
                       N_2_atmos)

  USE mo_exception,        ONLY: finish

  INTEGER,            INTENT(in)   ::  nidx, ntiles        ! number of grid boxes, number of tiles
  INTEGER,            INTENT(in)   ::  lcc_scheme          ! 
  REAL(dp),           INTENT(in)   ::  &
    Npool_green(:,:), Npool_woods(:,:), Npool_mobile(:,:)  ! Nitrogen pools of living material
  REAL(dp),           INTENT(in)   ::  scale_fac(:,:)      ! Scaling factor for N fluxes relating to lcc
  REAL(dp), OPTIONAL, INTENT(in)   ::  &
    Npool_litter_green_ag(:,:), Npool_litter_wood_ag(:,:)  ! Nitrogen pools of above ground litter
  REAL(dp), OPTIONAL, INTENT(out)  ::  &
    Npool_green_loss(:,:),             &                   ! Nitrogen lost from green pool
    Npool_wood_loss(:,:),              &                   ! Nitrogen lost from wood pool
    Npool_mobile_loss(:,:),            &                   ! Nitrogen lost from mobile plant pool
    N_litter_green_ag_loss(:,:),       &                   ! Nitrogen lost from green above bround litter
    N_litter_wood_ag_loss(:,:)                             ! Nitrogen lost from woody above ground litter
  REAL(dp), OPTIONAL, INTENT(out)  ::  N_2_atmos(:)        ! Nitrogen lost to atmosphere
  ! locals
  REAL(dp)  ::  NGreenLoss(nidx,ntiles), NWoodLoss(nidx,ntiles), NMobileLoss(nidx,ntiles)
  REAL(dp)  ::  NLitterGreenAGLoss(nidx,ntiles), NLitterWoodAGLoss(nidx,ntiles)

  NGreenLoss (:,:) = scale_fac(:,:) * Npool_green(:,:)
  NWoodLoss  (:,:) = scale_fac(:,:) * Npool_woods(:,:)
  NMobileLoss(:,:) = scale_fac(:,:) * Npool_mobile(:,:)

  IF (PRESENT(Npool_litter_green_ag)) NLitterGreenAGLoss(:,:) = scale_fac(:,:) * Npool_litter_green_ag(:,:)
  IF (PRESENT(Npool_litter_wood_ag)) NLitterWoodAGLoss(:,:) = scale_fac(:,:) * Npool_litter_wood_ag(:,:)

  IF (PRESENT(Npool_green_loss))  Npool_green_loss (:,:) = NGreenLoss (:,:)
  IF (PRESENT(Npool_wood_loss))   Npool_wood_loss  (:,:) = NWoodLoss  (:,:)
  IF (PRESENT(Npool_mobile_loss)) Npool_mobile_loss(:,:) = NMobileLoss(:,:)
  IF (PRESENT(N_litter_green_ag_loss)) N_litter_green_ag_loss(:,:) = NLitterGreenAGLoss(:,:)
  IF (PRESENT(N_litter_wood_ag_loss))  N_litter_wood_ag_loss(:,:) = NLitterWoodAGLoss(:,:)
  IF (lcc_scheme == 2) THEN
     IF (.NOT. PRESENT(N_2_atmos)) &
        CALL finish('N_loss_lcc()','flux to atmos must be determined, if anthro pools are used')
     IF (PRESENT(Npool_litter_green_ag) .AND. PRESENT(Npool_litter_wood_ag)) THEN
        N_2_atmos(:) = SUM((NGreenLoss (:,:) + NMobileLoss(:,:)) * frac_green_aboveGround + &
                       NWoodLoss(:,:) * frac_wood_aboveGround +                         &
                       NLitterGreenAGLoss(:,:) + NLitterWoodAGLoss(:,:), DIM=2)
     ELSE
        N_2_atmos(:) = SUM((NGreenLoss (:,:) + NMobileLoss(:,:)) * frac_green_aboveGround + &
                       NWoodLoss(:,:) * frac_wood_aboveGround, DIM=2)
     END IF
  END IF

  END SUBROUTINE N_loss_lcc


  PURE SUBROUTINE yasso (Yasso_io_pools, Weather, litter, Lit_coefV,WoodLitterSize, Yasso_out, frac_aboveground, &
                         root_exudates, redFact_Nlimit)
    !! DSG: 11.03.2013
    !! modified Yasso soil C model using an Euler scheme and readable coding style
    !! JSBACH specific modification of model structure: AWEN pools are separated into aboveground and belowground pools
    !! 
    !! Herbivory flux is treated like normal litter, because the herbivory fluxes in JSBACH are not evaluated, yet.
    !! If the represenation of herbivory in JSBACH is evaluated I recommend to treat herbivory flux in YASSO different from general
    !! litter.
    !!
    !! The routine needs as input 1. daily litter flux
    !!                            2. Yasso pools
    !!                            3. PFT specific parameters
    !!
    !! output is on a daily timestep: respiration
    !! output pools: AWEN + humus 
    !! The routine must be called twice, first for non-woody litter
    !! then for woody litter
    !!-------------------------------------------------------------------------------------------------------------------------
    !! Yasso - Soil carbon model (Jari Liski)
    !!
    !! The code for the yasso model has been developed and is provided by the Finnish Environment Institute SYKE. It is 
    !! distributed as part of the MPI Earth System Model under the MPI-M license conditions as stated in the header of this module. 
    !!
    !! model to calculate the amount of soil organic carbon, changes in the amount of soil organic carbon and
    !! heterotrophic soil respiration.
    !! For documention of yasso see Liski et al 2005, Tuomi et al 2008,2009,2011
    !!
    !! implementation: Tea Thum (FMI), Petri Räisänen (FMI), Daniel Goll (MPI-M)
    !!-------------------------------------------------------------------------------------------------------------------------

    IMPLICIT NONE
    
    REAL(dp), DIMENSION(9),  INTENT(IN)  :: Yasso_io_pools    ! Yasso pools IN                                        [mol(c)/m2]
    REAL(dp),                INTENT(IN)  :: WoodLitterSize    ! Size of woody litter; 
                                                              !     is zero when the subroutine is called for non-woody litter
    REAL(dp), DIMENSION(2),  INTENT(IN)  :: Weather           ! climatic drivers: air temperature and precipitation, 15-d means
    REAL(dp),                INTENT(IN)  :: litter            ! fresh litter inputs  (above and belowground)         [mol(C)/m2/day]
    REAL(dp), DIMENSION(5),  INTENT(IN)  :: Lit_coefV         ! fractions to seperate incoming litter to the yasso pool    [ ]
    REAL(dp),                INTENT(IN)  :: frac_aboveground  ! parameter to seperate litter influx into above- and belowground 
                                                              !     part                                                   [ ]
    REAL(dp),                INTENT(IN)  :: root_exudates     ! root exudates                                        [mol(C)/m2/day]

    REAL(dp), DIMENSION(16), INTENT(OUT) :: Yasso_out         ! updated Yasso pools(1:9) [mol(c)/m2] & respiration(10) [mol(c)/m2/d]
                                                              ! & Cflx_2_humus (11) [mol(c)/m2/d] ... see update_Cpools for more 
                                                              ! information

    REAL(dp), OPTIONAL, INTENT(IN)       ::  redFact_Nlimit

    !local
    ! Yasso pools   [mol(c)/m2]
    REAL(dp)     :: YCpool_acid_ag    ! Aboveground  
    REAL(dp)     :: YCpool_water_ag  
    REAL(dp)     :: YCpool_ethanol_ag 
    REAL(dp)     :: YCpool_nonsoluble_ag
    REAL(dp)     :: YCpool_acid_bg    ! Belowground  
    REAL(dp)     :: YCpool_water_bg  
    REAL(dp)     :: YCpool_ethanol_bg 
    REAL(dp)     :: YCpool_nonsoluble_bg
    
    REAL(dp)     :: YCpool_humus   

    ! Yasso fluxes (internally used, therefore annual fluxes)  
    REAL(dp)     :: Cflx_litter_2_acid             ! Litter influx to acid soluble pools    [mol(c)/m2/a]   
    REAL(dp)     :: Cflx_litter_2_water            ! Litter influx to water soluble pools   [mol(c)/m2/a]   
    REAL(dp)     :: Cflx_litter_2_ethanol          ! Litter influx to ethanol soluble pools [mol(c)/m2/a]   
    REAL(dp)     :: Cflx_litter_2_nonsoluble       ! Litter influx to non soluble pools     [mol(c)/m2/a]   
    REAL(dp)     :: Cflx_litter_2_humus            ! Litter influx to humus soluble pool    [mol(c)/m2/a]   
  
    ! Aboveground
    REAL(dp)     :: Cflx_from_acid_ag              ! Loss flux of acid soluble pool        [mol(c)/m2/a]
    REAL(dp)     :: Cflx_from_water_ag             ! Loss flux of water soluble pool       [mol(c)/m2/a]
    REAL(dp)     :: Cflx_from_ethanol_ag           ! Loss flux of ethanol soluble pool     [mol(c)/m2/a]
    REAL(dp)     :: Cflx_from_nonsoluble_ag        ! Loss flux of non soluble pool         [mol(c)/m2/a]
    
    REAL(dp)     :: Cflx_2_acid_ag                 ! Gain flux of acid soluble pool        [mol(c)/m2/a]
    REAL(dp)     :: Cflx_2_water_ag                ! Gain flux of water soluble pool       [mol(c)/m2/a]              
    REAL(dp)     :: Cflx_2_ethanol_ag              ! Gain flux of ethanol soluble pool     [mol(c)/m2/a]
    REAL(dp)     :: Cflx_2_nonsoluble_ag           ! Gain flux of non soluble pool         [mol(c)/m2/a]
     
    ! Belowground
    REAL(dp)     :: Cflx_from_acid_bg              ! Loss flux of acid soluble pool        [mol(c)/m2/a]
    REAL(dp)     :: Cflx_from_water_bg             ! Loss flux of water soluble pool       [mol(c)/m2/a]
    REAL(dp)     :: Cflx_from_ethanol_bg           ! Loss flux of ethanol soluble pool     [mol(c)/m2/a]
    REAL(dp)     :: Cflx_from_nonsoluble_bg        ! Loss flux of non soluble pool         [mol(c)/m2/a]
    REAL(dp)     :: Cflx_from_humus                ! Loss flux of humus pool               [mol(c)/m2/a]
    
    REAL(dp)     :: Cflx_2_acid_bg                 ! Gain flux of acid soluble pool        [mol(c)/m2/a]
    REAL(dp)     :: Cflx_2_water_bg                ! Gain flux of water soluble pool       [mol(c)/m2/a]              
    REAL(dp)     :: Cflx_2_ethanol_bg              ! Gain flux of ethanol soluble pool     [mol(c)/m2/a]
    REAL(dp)     :: Cflx_2_nonsoluble_bg           ! Gain flux of non soluble pool         [mol(c)/m2/a]
    REAL(dp)     :: Cflx_2_humus                   ! Gain flux of humus pool               [mol(c)/m2/a]
    REAL(dp)     :: Cflx_2_humusAG                 ! Gain flux of humus pool; from ag only [mol(c)/m2/a]
    REAL(dp)     :: Cflx_2_humusBG                 ! Gain flux of humus pool; from bg only [mol(c)/m2/a]
    
    ! Respiration is output of Yasso therefore a daily flux 
    REAL(dp)     :: soilResp_rateYasso             ! Flux to the atmosphere                [mol(c)/m2/d]  
    REAL(dp)     :: soilResp_rateLitterAG          ! Flux to the atmosphere; from AWEN ag only  [mol(c)/m2/d]
    REAL(dp)     :: soilResp_rateLitterBG          ! Flux to the atmosphere; from AWEN bg only  [mol(c)/m2/d]

    ! Decomposition rates [1/a]
    REAL(dp)     :: d_acid
    REAL(dp)     :: d_water
    REAL(dp)     :: d_ethanol
    REAL(dp)     :: d_nonsoluble
    REAL(dp)     :: d_humus

    ! Dcomposition rate for N litter green (needed for N cycle)
    REAL(dp)     :: d_litter_green              ! the rate with which the sum of AWEN pools decays
    REAL(dp)     :: Pseudo_litter_green         ! Sum of all AWEN pools

    ! PFT-specific parameters
    REAL(dp)     :: frac_litter_2_acid        ! Fraction of incoming material which enters the acid-soluble fraction of litter
    REAL(dp)     :: frac_litter_2_water       ! Fraction of incoming material which enters the water-soluble fraction of litter
    REAL(dp)     :: frac_litter_2_ethanol     ! Fraction of incoming material which enters the ethanol-soluble fraction of litter
    REAL(dp)     :: frac_litter_2_nonsoluble  ! Fraction of incoming material which enters the nonsoluble fraction of litter
    REAL(dp)     :: frac_litter_2_humus       ! Fraction of incoming material which enters the humus pool

    ! Yasso parameters (Tuomi '11)
    ! reference decompositation rates for T = 0C and unlimited water availability  in [1/yr]... 
    REAL(dp), PARAMETER   :: ref_decomp_rate_acid       = -0.72_dp    !  ... acid-soluble fraction 
    REAL(dp), PARAMETER   :: ref_decomp_rate_water      = -5.9_dp     !  ... water-soluble fraction 
    REAL(dp), PARAMETER   :: ref_decomp_rate_ethanol    = -0.28_dp    !  ... ethanol-soluble fraction 
    REAL(dp), PARAMETER   :: ref_decomp_rate_nonsoluble = -0.031_dp   !  ... nonsoluble-soluble fraction 
    REAL(dp), PARAMETER   :: ref_decomp_rate_humus      = -0.0016_dp  !  ... humus 
    ! parameters for temperature dependence
    REAL(dp), PARAMETER   :: temp_p1                    = 0.095_dp    ! first  order temperature dependence [1/C]
    REAL(dp), PARAMETER   :: temp_p2                    = -0.0014_dp  ! second order temperature dependence [1/C]
    ! parameter for precipitation dependence
    REAL(dp), PARAMETER   :: precip_p1                  = -1.21_dp    ! first  order precipitation dependence [a/m]
    ! parameters for size dependence 
    REAL(dp), PARAMETER   :: size_p1                    = -1.71_dp    ! first  order size dependence [1/cm]
    REAL(dp), PARAMETER   :: size_p2                    = 0.86_dp     ! second order size dependence [1/cm]
    REAL(dp), PARAMETER   :: size_p3                    = -0.306_dp   ! size dependence power [ ]
    
    ! Fractions to split the decomposition flux of each pool into fluxes to the other AWEN pools and atmosphere 
    ! REMARK: there is no a priori (to the calibration procedure of YASSO) assumption about possible fluxes between pools, 
    ! therefore all fluxes are theoretically possible, however some of the fraction are zero which means the corresponding flux 
    ! does not exists. However, after a recalibration of the model this fraction may change.
    !
    ! fractions of decomposition flux of acid soluble pool going to the other AWEN pools [ ]
    ! (8,11,14)
    REAL(dp), PARAMETER   :: A_2_W                      =  0.99_dp
    REAL(dp), PARAMETER   :: A_2_E                      =  0.00_dp
    REAL(dp), PARAMETER   :: A_2_N                      =  0.00_dp
    ! fractions of decomposition flux of water soluble pool going to the other AWEN pools [ ]
    ! (5,12,15)
    REAL(dp), PARAMETER   :: W_2_A                      = 0.48_dp
    REAL(dp), PARAMETER   :: W_2_E                      = 0.00_dp
    REAL(dp), PARAMETER   :: W_2_N                      = 0.015_dp
    ! fractions of decomposition flux of ethanol soluble pool going to the other AWEN pools [ ]
    ! (6,9,16)
    REAL(dp), PARAMETER   :: E_2_A                      = 0.01_dp
    REAL(dp), PARAMETER   :: E_2_W                      = 0.00_dp
    REAL(dp), PARAMETER   :: E_2_N                      = 0.95_dp
    ! fractions of decomposition flux of non soluble pool going to the other AWEN pools [ ]
    ! (7,10,13)
    REAL(dp), PARAMETER   :: N_2_A                      = 0.83_dp
    REAL(dp), PARAMETER   :: N_2_W                      = 0.01_dp
    REAL(dp), PARAMETER   :: N_2_E                      = 0.02_dp
    ! fraction of decomposition fluxes of AWEN pools which enters the humus pool [ ]
    REAL(dp), PARAMETER   :: AWEN_2_H                   = 0.0045_dp

    ! climatic drivers
    REAL(dp)     :: precip                  ! 15 runnning mean of precipitation      [m/a]
    REAL(dp)     :: temp                    ! 15 runnning mean of 2m air temperature [C]

    REAL(dp)     :: d_temp                  ! term which accounts for the temperature influence on litter decomposition
    REAL(dp)     :: d_precip                ! term which accounts for the precipitation influence on litter decomposition
    REAL(dp)     :: d_size                  ! term which accounts for the litter size influence on litter decomposition
     
    REAL(dp)     :: redFactor  ! dummy for redFact_Nlimit

    ! time stepping; the parameterization of yasso is done on a annual time
    ! step, thus the fluxes have to be scaled up to annual rates. 
    REAL(dp)     :: dt                      ! time step of mo_update_CPools [year]             
    dt = 1.0_dp/days_per_year   

    IF(PRESENT(redFact_Nlimit)) THEN ! nitrogen cycle active & we know already the Nlimitation factor
        redFactor  = redFact_Nlimit
    ELSE
        redFactor  = 1.0_dp
    ENDIF

    ! 0. Preparations
     ! initialize parameters
     frac_litter_2_acid       = Lit_coefV(1)
     frac_litter_2_water      = Lit_coefV(2)
     frac_litter_2_ethanol    = Lit_coefV(3)
     frac_litter_2_nonsoluble = Lit_coefV(4)
     frac_litter_2_humus      = Lit_coefV(5)

     ! initialize yasso pools (AWEN + humus)
     YCpool_acid_ag        = Yasso_io_pools(1)
     YCpool_water_ag       = Yasso_io_pools(2)
     YCpool_ethanol_ag     = Yasso_io_pools(3)
     YCpool_nonsoluble_ag  = Yasso_io_pools(4)
     
     YCpool_acid_bg        = Yasso_io_pools(5)
     YCpool_water_bg       = Yasso_io_pools(6)
     YCpool_ethanol_bg     = Yasso_io_pools(7)
     YCpool_nonsoluble_bg  = Yasso_io_pools(8)
     
     YCpool_humus          = Yasso_io_pools(9)
     
    ! get the sum of the AWEN pools before decomposition
    Pseudo_litter_green = YCpool_acid_ag          &
                        + YCpool_water_ag         &
                        + YCpool_ethanol_ag       &
                        + YCpool_nonsoluble_ag    &
                        + YCpool_acid_bg          &
                        + YCpool_water_bg         &
                        + YCpool_ethanol_bg       &
                        + YCpool_nonsoluble_bg

     ! Change units of climatic forcing variables
     precip = Weather(2)*sec_per_year/1000._dp  ! mm/s -> m/a 
     temp   = Weather(1) - Tmelt                ! K -> C
    
     ! Calculate annual litter influxes [mol(c)/m2/a]
     Cflx_litter_2_acid       = frac_litter_2_acid        *litter * days_per_year
     Cflx_litter_2_water      = frac_litter_2_water       *litter * days_per_year
     Cflx_litter_2_ethanol    = frac_litter_2_ethanol     *litter * days_per_year
     Cflx_litter_2_nonsoluble = frac_litter_2_nonsoluble  *litter * days_per_year
     Cflx_litter_2_humus      = frac_litter_2_humus       *litter * days_per_year


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!! START CALCULATION OF SOIL CARBON !!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! 1. Calculate decomposition rates 
     
     ! Temperature dependence of decomposition
     d_temp   = EXP(temp_p1*temp + temp_p2*temp**2.0_dp)  
     ! Precipitation dependence of decomposition
     d_precip = 1.0_dp - EXP(precip_p1*precip)
     ! Litter size dependence of decomposition -- no effect if WoodlitterSize = 0.0 
     d_size   = MIN(1.0_dp,(1.0_dp + size_p1 * WoodLitterSize + size_p2 * WoodLitterSize**2.0_dp)**size_p3)

     ! decomposition rates accounting for temperature, precipitation, litter size, and nutrient limitation 
     d_acid       =  redFactor * ref_decomp_rate_acid       *d_temp *d_precip *d_size
     d_water      =  redFactor * ref_decomp_rate_water      *d_temp *d_precip *d_size
     d_ethanol    =  redFactor * ref_decomp_rate_ethanol    *d_temp *d_precip *d_size
     d_nonsoluble =  redFactor * ref_decomp_rate_nonsoluble *d_temp *d_precip *d_size

     ! the decomposition of humus is not limited by nutrients; thus no redFactor
     d_humus      =  ref_decomp_rate_humus      *d_temp *d_precip           ! no size effect on humus 
     
    ! 2. Calculate fluxes 
     
     ! loss fluxes (negative values): 
     Cflx_from_acid_ag         = YCpool_acid_ag       * d_acid        
     Cflx_from_water_ag        = YCpool_water_ag      * d_water                  
     Cflx_from_ethanol_ag      = YCpool_ethanol_ag    * d_ethanol             
     Cflx_from_nonsoluble_ag   = YCpool_nonsoluble_ag * d_nonsoluble       
     
     Cflx_from_acid_bg         = YCpool_acid_bg       * d_acid        
     Cflx_from_water_bg        = YCpool_water_bg      * d_water                  
     Cflx_from_ethanol_bg      = YCpool_ethanol_bg    * d_ethanol             
     Cflx_from_nonsoluble_bg   = YCpool_nonsoluble_bg * d_nonsoluble       
     
     Cflx_from_humus           = YCpool_humus         * d_humus         

     ! gain fluxes (positive values): 
     ! fixed fractions of each loss flux enters another pool; REMARK: the fraction can be zero (see above why) 
     Cflx_2_acid_ag           = ABS(Cflx_from_water_ag      * W_2_A  &    ! returns positive fluxes
                                  + Cflx_from_ethanol_ag    * E_2_A  &
                                  + Cflx_from_nonsoluble_ag * N_2_A) 

     Cflx_2_water_ag          = ABS(Cflx_from_acid_ag       * A_2_W  &
                                  + Cflx_from_ethanol_ag    * E_2_W  &
                                  + Cflx_from_nonsoluble_ag * N_2_W) 

     Cflx_2_ethanol_ag        = ABS(Cflx_from_acid_ag       * A_2_E  &
                                  + Cflx_from_water_ag      * W_2_E  &
                                  + Cflx_from_nonsoluble_ag * N_2_E) 

     Cflx_2_nonsoluble_ag     = ABS(Cflx_from_acid_ag       * A_2_N  &
                                  + Cflx_from_water_ag      * W_2_N  &
                                  + Cflx_from_ethanol_ag    * E_2_N) 

     Cflx_2_acid_bg           = ABS(Cflx_from_water_bg      * W_2_A  &    
                                  + Cflx_from_ethanol_bg    * E_2_A  &
                                  + Cflx_from_nonsoluble_bg * N_2_A) 

     Cflx_2_water_bg          = ABS(Cflx_from_acid_bg       * A_2_W  &
                                  + Cflx_from_ethanol_bg    * E_2_W  &
                                  + Cflx_from_nonsoluble_bg * N_2_W) 

     Cflx_2_ethanol_bg        = ABS(Cflx_from_acid_bg       * A_2_E  &
                                  + Cflx_from_water_bg      * W_2_E  &
                                  + Cflx_from_nonsoluble_bg * N_2_E) 

     Cflx_2_nonsoluble_bg     = ABS(Cflx_from_acid_bg       * A_2_N  &
                                  + Cflx_from_water_bg      * W_2_N  &
                                  + Cflx_from_ethanol_bg    * E_2_N) 

     Cflx_2_humus          = ABS(Cflx_from_acid_ag                &
                               + Cflx_from_water_ag               &
                               + Cflx_from_ethanol_ag             &
                               + Cflx_from_nonsoluble_ag          &
                               + Cflx_from_acid_bg                &
                               + Cflx_from_water_bg               &
                               + Cflx_from_ethanol_bg             &
                               + Cflx_from_nonsoluble_bg          &
                               ) * AWEN_2_H  

     Cflx_2_humusAG        = ABS(Cflx_from_acid_ag                &
                               + Cflx_from_water_ag               &
                               + Cflx_from_ethanol_ag             &
                               + Cflx_from_nonsoluble_ag          &
                               ) * AWEN_2_H

     Cflx_2_humusBG        = ABS(Cflx_from_acid_bg                &
                               + Cflx_from_water_bg               &
                               + Cflx_from_ethanol_bg             &
                               + Cflx_from_nonsoluble_bg          &
                               ) * AWEN_2_H

     ! the remaining fractions of the loss fluxes enter the atmosphere as respiration
     soilResp_rateYasso    =  (Cflx_from_acid_ag + Cflx_from_acid_bg)               &
                              * (1.0_dp - A_2_W - A_2_E - A_2_N - AWEN_2_H)         &
                              + (Cflx_from_water_ag + Cflx_from_water_bg)           &      
                              * (1.0_dp - W_2_A - W_2_E - W_2_N - AWEN_2_H)         &
                              + (Cflx_from_ethanol_ag + Cflx_from_ethanol_bg)       &     
                              * (1.0_dp - E_2_A - E_2_W - E_2_N - AWEN_2_H)         &
                              + (Cflx_from_nonsoluble_ag + Cflx_from_nonsoluble_bg) &  
                              * (1.0_dp - N_2_A - N_2_W - N_2_E - AWEN_2_H)         &
                              + Cflx_from_humus 

     ! litter ag & bg respiration (needed for N cycle)                       
     soilResp_rateLitterAG   =  (Cflx_from_acid_ag )                                &
                              * (1.0_dp - A_2_W - A_2_E - A_2_N - AWEN_2_H)         &
                              + (Cflx_from_water_ag )                               &
                              * (1.0_dp - W_2_A - W_2_E - W_2_N - AWEN_2_H)         &
                              + (Cflx_from_ethanol_ag )                             &
                              * (1.0_dp - E_2_A - E_2_W - E_2_N - AWEN_2_H)         &
                              + (Cflx_from_nonsoluble_ag )                          &
                              * (1.0_dp - N_2_A - N_2_W - N_2_E - AWEN_2_H)         

     soilResp_rateLitterBG   =  (Cflx_from_acid_bg )                                &
                              * (1.0_dp - A_2_W - A_2_E - A_2_N - AWEN_2_H)         &
                              + (Cflx_from_water_bg )                               &
                              * (1.0_dp - W_2_A - W_2_E - W_2_N - AWEN_2_H)         &
                              + (Cflx_from_ethanol_bg )                             &
                              * (1.0_dp - E_2_A - E_2_W - E_2_N - AWEN_2_H)         &
                              + (Cflx_from_nonsoluble_bg )                          &
                              * (1.0_dp - N_2_A - N_2_W - N_2_E - AWEN_2_H)

   ! 3. update Yasso pools
     YCpool_acid_ag           = MAX(0.0_dp,                                         &
                              YCpool_acid_ag                                        & ! old pool   
                              + (Cflx_from_acid_ag                                  & ! incoming flux from AWEN pools
                                + Cflx_2_acid_ag                                    & ! outgoing flux
                                + frac_aboveground * Cflx_litter_2_acid) *dt)         ! incoming flux from litter
     YCpool_water_ag          = MAX(0.0_dp,                                         & ! and so on .....
                              YCpool_water_ag                                       &
                              + (Cflx_from_water_ag                                 &
                                + Cflx_2_water_ag                                   &
                                + frac_aboveground * Cflx_litter_2_water) *dt)
     YCpool_ethanol_ag        = MAX(0.0_dp,                                         &
                              YCpool_ethanol_ag                                     &
                              + (Cflx_from_ethanol_ag                               &
                                + Cflx_2_ethanol_ag                                 &
                                + frac_aboveground * Cflx_litter_2_ethanol) *dt)
     YCpool_nonsoluble_ag     = MAX(0.0_dp,                                         &
                              YCpool_nonsoluble_ag                                  &
                              + (Cflx_from_nonsoluble_ag                            &
                                + Cflx_2_nonsoluble_ag                              &
                                + frac_aboveground * Cflx_litter_2_nonsoluble) *dt)
    

     YCpool_acid_bg           = MAX(0.0_dp,                                         &
                              YCpool_acid_bg                                        & ! old pool   
                              + (Cflx_from_acid_bg                                  & ! incoming flux from AWEN pools
                                + Cflx_2_acid_bg                                    & ! outgoing flux
                                + (1._dp - frac_aboveground) *                      & ! 
                                  Cflx_litter_2_acid) *dt)                            ! incoming flux from litter
     YCpool_water_bg          = MAX(0.0_dp,                                         & ! ...
                              YCpool_water_bg                                       &
                              + (Cflx_from_water_bg                                 &
                                + Cflx_2_water_bg                                   &
                                + (1._dp - frac_aboveground) *                      &
                                  Cflx_litter_2_water) *dt                          &
                              + root_exudates)                                    ! exudates are carbohydrates only and belowground 
     YCpool_ethanol_bg        = MAX(0.0_dp,                                         &
                              YCpool_ethanol_bg                                     &
                              + (Cflx_from_ethanol_bg                               &
                                + Cflx_2_ethanol_bg                                 &
                                + (1._dp - frac_aboveground) *                      &
                                     Cflx_litter_2_ethanol) *dt)
     YCpool_nonsoluble_bg     = MAX(0.0_dp,                                         &
                              YCpool_nonsoluble_bg                                  &
                              + (Cflx_from_nonsoluble_bg                            &
                                + Cflx_2_nonsoluble_bg                              &
                                + (1._dp - frac_aboveground) *                      &
                                  Cflx_litter_2_nonsoluble) *dt)

     YCpool_humus          = MAX(0.0_dp,                            &
                              YCpool_humus                          &
                              + (Cflx_from_humus                    &  
                                + Cflx_2_humus                      &
                                + Cflx_litter_2_humus) *dt)

     ! compute d_litter_green; this is given by the change in AWEN pools
     !          Flx                 Pool(t)                   Pool(t+1)
     !    d_P = ----     &   Flx = ----------   --->   d_P = -----------  - 1
     !          Pool(t)             Pool(t+1)                 Pool(t)

     IF (Pseudo_litter_green .GT. 0.0_dp) THEN
        d_litter_green = (YCpool_acid_ag            &
                          + YCpool_water_ag         &
                          + YCpool_ethanol_ag       &
                          + YCpool_nonsoluble_ag    &
                          + YCpool_acid_bg          &
                          + YCpool_water_bg         &
                          + YCpool_ethanol_bg       &
                          + YCpool_nonsoluble_bg)   &
                          / Pseudo_litter_green - 1._dp
      ELSE
         d_litter_green = 0.0_dp
      ENDIF

   ! 4.1 Write pools into the output array
    Yasso_out(1)           = YCpool_acid_ag 
    Yasso_out(2)           = YCpool_water_ag
    Yasso_out(3)           = YCpool_ethanol_ag
    Yasso_out(4)           = YCpool_nonsoluble_ag
    Yasso_out(5)           = YCpool_acid_bg 
    Yasso_out(6)           = YCpool_water_bg
    Yasso_out(7)           = YCpool_ethanol_bg
    Yasso_out(8)           = YCpool_nonsoluble_bg
    Yasso_out(9)           = YCPool_humus
   ! 4.2 Write respiration into the output array
    Yasso_out(10)           = soilResp_rateYasso * dt ! convert back to daily

    IF (.NOT.PRESENT(redFact_Nlimit)) THEN
       ! the routine might have been called to diagnose the nutrient demand, so we give out 
       ! 1. flux for immobilisation: this is the sum of all carbon coming from the non-humus pools to humus
       Yasso_out(11)           = Cflx_2_humusAG * dt ! daily ! the flux is positive
       Yasso_out(12)           = Cflx_2_humusBG * dt ! daily ! the flux is positive
       ! 2. flux for mineralisation:  this is the flux from humus to atmosphere (Cflx_slow_2_atmos)
       Yasso_out(13)           =  - Cflx_from_humus * dt ! daily ! sign CHECKED!
       ! ... as well as all fluxes from AWEN pools to atmosphere
       Yasso_out(14)           = -soilResp_rateLitterAG * dt !daily ! the flux is negative
       Yasso_out(15)           = -soilResp_rateLitterBG * dt !daily ! the flux is negative
       ! 3. decomposition rate for Npool litter green 
       Yasso_out(16)           = d_litter_green
    ENDIF

  END SUBROUTINE yasso

end module mo_cbal_cpools
