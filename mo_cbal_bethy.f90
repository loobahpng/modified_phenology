!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
module mo_cbal_bethy
!
! Computation of the net primary productivity (NPP) as in BETHY. 
!

  USE mo_jsbach_grid,      ONLY: grid_type, domain_type, kstart, kend, nidx
  USE mo_kind,             ONLY: dp
  USE mo_exception,        ONLY: finish, message
  USE mo_linked_list,      ONLY: t_stream
  USE mo_netCDF,           ONLY: FILE_INFO, NF_MAX_NAME
  USE mo_io,               ONLY: IO_open, IO_read, IO_close
  USE mo_mpi,              ONLY: p_parallel_io, p_io, p_parallel, p_bcast
  USE mo_jsbach_constants, ONLY: molarMassCO2_kg, molarMassN_kg, molarMassN2_kg, molarMassN2O_kg
  USE mo_time_control,     ONLY: delta_time    !! time step in seconds
  USE mo_jsbach,           ONLY: debug, debug_Cconservation, test_Cconservation, debug_Nconservation, test_Nconservation

  IMPLICIT NONE

  ! === BEGIN OF PUBLIC PART =======================================================================================================


  TYPE cbalance_type  ! contains the state variables of the cbalance module (dims: domain%nland x vegetation%ntiles)
     INTEGER          :: ntiles
     REAL(dp),POINTER :: Cpool_green(:,:)    !! C-pool for leaves, fine roots, vegetative organs and other green (living) parts 
                                             !! .. of vegetation [mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: Cpool_reserve(:,:)  !! C-pool for carbohydrate reserve (sugars, starches) that allows plants to survive 
                                             !! .. bad times[mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: Cpool_woods(:,:)    !! C-pool for stems, thick roots and other (dead) structural material of living
                                             !! .. plants [mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: Cpool_crop_harvest(:,:)    !! C-pool for biomass harvested from crops [mol(C)/m^2(grid box)]     
     REAL(dp),POINTER :: Cpool_litter_green_ag(:,:) !! C-pool for above ground non-woody litter especially from leaves
                                                    !!     [mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: Cpool_litter_green_bg(:,:) !! C-pool for below ground non-woody litter especially from fine roots
                                                    !!     [mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: Cpool_litter_wood_ag(:,:)  !! C-pool for litter originating from above ground woody parts of the plants
                                                    !!     [mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: Cpool_litter_wood_bg(:,:)  !! C-pool for litter originating from below ground woody parts of the plants
                                                    !!     [mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: Cpool_slow(:,:)            !! C-pool for below ground organic material (coming from the litter pools)
                                                    !!     [mol(C)/m^2(canopy)]

     REAL(dp),POINTER :: pseudo_temp(:)             !! a proxy for 15-day mean air temperature, [K] needed with yasso
     REAL(dp),POINTER :: pseudo_precip(:)           !! a proxy for 15-day mean preciptation rate [Kg/m^2 s] needed with yasso
     !SIZE CLASS 1
     REAL(dp),POINTER :: YCpool_acid_ag1(:,:)       !! Yasso above ground litter-pool for acid soluble litter [mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: YCpool_water_ag1(:,:)      !! Yasso above ground litter-pool for water soluble litter [mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: YCpool_ethanol_ag1(:,:)    !! Yasso above ground litter-pool for ethanol soluble litter [mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: YCpool_nonsoluble_ag1(:,:) !! Yasso above ground litter-pool for non-soluble litter [mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: YCpool_acid_bg1(:,:)       !! Yasso below ground litter-pool for acid soluble litter [mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: YCpool_water_bg1(:,:)      !! Yasso below ground litter-pool for water soluble litter [mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: YCpool_ethanol_bg1(:,:)    !! Yasso below ground litter-pool for ethanol soluble litter [mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: YCpool_nonsoluble_bg1(:,:) !! Yasso below ground litter-pool for non-soluble litter [mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: YCpool_humus_1(:,:)        !! Yasso below ground litter-pool for slow C compartment [mol(C)^2/m(canopy)]
     !SIZE CLASS 2
     REAL(dp),POINTER :: YCpool_acid_ag2(:,:)       !! Yasso above ground litter-pool for acid soluble litter [mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: YCpool_water_ag2(:,:)      !! Yasso above ground litter-pool for water soluble litter [mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: YCpool_ethanol_ag2(:,:)    !! Yasso above ground litter-pool for ethanol soluble litter [mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: YCpool_nonsoluble_ag2(:,:) !! Yasso above ground litter-pool for non-soluble litter [mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: YCpool_acid_bg2(:,:)       !! Yasso below ground litter-pool for acid soluble litter [mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: YCpool_water_bg2(:,:)      !! Yasso below ground litter-pool for water soluble litter [mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: YCpool_ethanol_bg2(:,:)    !! Yasso below ground litter-pool for ethanol soluble litter [mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: YCpool_nonsoluble_bg2(:,:) !! Yasso below ground litter-pool for non-soluble litter [mol(C)/m^2(canopy)]
     REAL(dp),POINTER :: YCpool_humus_2(:,:)        !! Yasso below ground litter-pool for slow C compartment [mol(C)^2/m(canopy)]

     REAL(dp),POINTER :: NPP_Rate(:,:)       !! The instantaneous NPP rate [mol(C)/m^2(canopy) s]
     REAL(dp),POINTER :: LAI_sum(:,:)        !! used to accumulate LAI over a day
     REAL(dp),POINTER :: NPP_sum(:,:)        !! used to accumulated NPP-Rate over a day 
     REAL(dp),POINTER :: GPP_sum(:,:)        !! used to accumulated GPP-Rate over a day     
     REAL(dp),POINTER :: topSoilTemp_sum(:)  !! used to accumulated upper layer soil temperature over a day
     REAL(dp),POINTER :: alpha_sum(:,:)      !! used to accumulated water stress coefficient alpha over a day
     REAL(dp),POINTER :: drainage_sum(:)     !! used to accumulate drainage over a day
     
     REAL(dp),POINTER :: soil_respiration(:,:)    !! mean daily rate of heterotrophic (soil) respiration [mol(CO2)/m^2(ground)] 
     REAL(dp),POINTER :: soil_respiration_pot(:,:) !! mean daily hetrespiration without N limitation  [mol(CO2)/m^2(ground)]
     REAL(dp),POINTER :: NPP_flux_correction(:,:) !! Daily updated flux correction from yesterdays carbon balance 
                                                  !!    [mol(CO2)/m^2(canopy) s]
     REAL(dp),POINTER :: excess_NPP(:,:)          !! That part of NPP that because of structural limits could not be allocated in
                                                  !!    carbon pools [mol(CO2)/m^2(canopy) s]
     REAL(dp),POINTER :: NPP_act_yDayMean(:,:)    !! mean value of actual NPP-Rate yesterday, i.e. after N-limitation
                                                  !!    [mol(CO2)/(m^2(canopy) s)]

     REAL(dp),POINTER :: root_exudates(:,:)       !! Total root exudates entering to the litter green pools [mol(C)/m^2(canopy) s]
     REAL(dp),POINTER :: Cflux_herbivory(:,:)        
     REAL(dp),POINTER :: Cflux_herbivory_LG(:,:)    
     REAL(dp),POINTER :: Cflux_herbivory_2_atm(:,:)
     REAL(dp),POINTER :: Cflx_2_crop_harvest(:,:)
     REAL(dp),POINTER :: Cflx_crop_harvest_2_atm(:,:)

     ! variables needed with Spitfire to compute fuel classes from wood and green pools and litter
     REAL(dp), POINTER :: frac_litter_wood_new(:,:)

     ! Anthropogenic pools [mol m-2(veg)]
     REAL(dp),POINTER :: Cpool_onSite(:)          !! Carbon left on the ground by landcover changes
     REAL(dp),POINTER :: Cpool_paper(:)             !! Carbon stored in short/intermediate term anthro pool from landcover changes
     REAL(dp),POINTER :: Cpool_construction(:)      !! Carbon stored in long term anthro pool from landcover changes
     REAL(dp),POINTER :: Cpool_onSite_harvest(:)    !! Carbon left on the ground by harvest
     REAL(dp),POINTER :: Cpool_paper_harvest(:)     !! Carbon stored in short/intermediate term anthro pool from harvest
     REAL(dp),POINTER :: Cpool_construction_harvest(:) !! Carbon stored in long term anthro pool from harvest

     ! FM production pools - grid-cell mass variables in m^2 grid box [mol(C)/m^2(grid box)]
     REAL(dp),POINTER :: FOM_fm_Box_burning_pool(:)            !! Carbon in the burning pool [mol(C)/m^2(grid box)]
     REAL(dp),POINTER :: FOM_fm_Box_paper_pool(:)              !! Carbon in the paper pool [mol(C)/m^2(grid box)]
     REAL(dp),POINTER :: FOM_fm_Box_construction_pool(:)       !! Carbon in the construction pool [mol(C)/m^2(grid box)]
     REAL(dp),POINTER :: FOM_fm_Box_eternity_pool(:)           !! Carbon in the eternity pool [mol(C)/m^2(grid box)]

  END TYPE cbalance_type

  TYPE cbalance_diag_type  !! contains diagnostic variables  (1.dim: domain%nland, 2.dim: vegetation%ntiles)
     REAL(dp),POINTER :: boxC_green(:,:)              !! Cpool_green() relative to grid box area [mol(C)/m^2(grid box)]
     REAL(dp),POINTER :: boxC_reserve(:,:)            !! Cpool_reserve() relative to grid box area [mol(C)/m^2(grid box)]
     REAL(dp),POINTER :: boxC_woods(:,:)              !! Cpool_woods() relative to grid box area [mol(C)/m^2(grid box)]
     REAL(dp),POINTER :: boxC_crop_harvest(:,:)       !! Cpool_crop_harvest() relative to grid box area [mol(C)/m^2(grid box)]
     REAL(dp),POINTER :: boxC_litter_green_ag(:,:)    !! Cpool_litter_green_ag() relative to grid box area [mol(C)/m^2(grid box)]
     REAL(dp),POINTER :: boxC_litter_green_bg(:,:)    !! Cpool_litter_green_bg() relative to grid box area [mol(C)/m^2(grid box)] 
     REAL(dp),POINTER :: boxC_litter_wood(:,:)        !! sum of Cpool_litter_wood_ag and Cpool_litter_wood_bg relative to grid box
                                                      !!     area [mol(C)/m^2(grid box)]
     REAL(dp),POINTER :: boxC_slow(:,:)               !! Cpool_slow() relative to grid box area [mol(C)/m^2(grid box)]
                                                      !! yasso above ground carbon pools in kg and relative to grid box area
     REAL(dp),POINTER :: boxYC_acid_ag1(:,:)           !!    as YCpool_acid_ag() but in [kg(C)/m^2(grid box)] 
     REAL(dp),POINTER :: boxYC_water_ag1(:,:)          !!    as YCpool_water_ag() but in [kg(C)/m^2(grid box)] 
     REAL(dp),POINTER :: boxYC_ethanol_ag1(:,:)        !!    as YCpool_ethanol_ag() but in [kg(C)/m^2(grid box)] 
     REAL(dp),POINTER :: boxYC_nonsoluble_ag1(:,:)     !!    as YCpool_nonsoluble_ag() but in [kg(C)/m^2(grid box)] 
                                                      !! yasso below ground carbon pools in kg and relative to grid box area
     REAL(dp),POINTER :: boxYC_acid_bg1(:,:)           !!    as YCpool_acid_bg() but in [kg(C)/m^2(grid box)] 
     REAL(dp),POINTER :: boxYC_water_bg1(:,:)          !!    as YCpool_water_bg() but in [kg(C)/m^2(grid box)] 
     REAL(dp),POINTER :: boxYC_ethanol_bg1(:,:)        !!    as YCpool_ethanol_bg() but in [kg(C)/m^2(grid box)] 
     REAL(dp),POINTER :: boxYC_nonsoluble_bg1(:,:)     !!    as YCpool_nonsoluble_bg() but in [kg(C)/m^2(grid box)] 
     REAL(dp),POINTER :: boxYC_humus_1(:,:)            !!    as YCpool_humus() but in [kg(C)/m^2(grid box)] 
     REAL(dp),POINTER :: boxYC_acid_ag2(:,:)           !!    as YCpool_acid_ag() but in [kg(C)/m^2(grid box)] 
     REAL(dp),POINTER :: boxYC_water_ag2(:,:)          !!    as YCpool_water_ag() but in [kg(C)/m^2(grid box)] 
     REAL(dp),POINTER :: boxYC_ethanol_ag2(:,:)        !!    as YCpool_ethanol_ag() but in [kg(C)/m^2(grid box)] 
     REAL(dp),POINTER :: boxYC_nonsoluble_ag2(:,:)     !!    as YCpool_nonsoluble_ag() but in [kg(C)/m^2(grid box)] 
                                                      !! yasso below ground carbon pools in kg and relative to grid box area
     REAL(dp),POINTER :: boxYC_acid_bg2(:,:)           !!    as YCpool_acid_bg() but in [kg(C)/m^2(grid box)] 
     REAL(dp),POINTER :: boxYC_water_bg2(:,:)          !!    as YCpool_water_bg() but in [kg(C)/m^2(grid box)] 
     REAL(dp),POINTER :: boxYC_ethanol_bg2(:,:)        !!    as YCpool_ethanol_bg() but in [kg(C)/m^2(grid box)] 
     REAL(dp),POINTER :: boxYC_nonsoluble_bg2(:,:)     !!    as YCpool_nonsoluble_bg() but in [kg(C)/m^2(grid box)] 
     REAL(dp),POINTER :: boxYC_humus_2(:,:)            !!    as YCpool_humus() but in [kg(C)/m^2(grid box)] 
     
     REAL(dp),POINTER :: pseudo_temp_yDay(:)          !! a proxy for 15-day mean air temperature at new_day, [K] needed with yasso
     REAL(dp),POINTER :: pseudo_precip_yDay(:)        !! a proxy for 15-day mean precipitation rate at new_day [Kg/m^2 s] with yasso

     REAL(dp),POINTER :: LAI_previousDayMean(:,:)     !! mean value of LAI the day before yesterday
     REAL(dp),POINTER :: LAI_yDayMean(:,:)            !! mean value of LAI yesterday (from  LAI_sum())
     REAL(dp),POINTER :: NPP_yDayMean(:,:)            !! mean value of NPP-Rate yesterday (from NPP_sum()) [mol(CO2)/(m^2(canopy)s)]

     REAL(dp),POINTER :: GPP_yDayMean(:,:)            !! mean value of GPP-Rate yesterday (from GPP_sum()) [mol(CO2)/(m^2(canopy)s)]
     REAL(dp),POINTER :: topSoilTemp_yDayMean(:)      !! mean value of upper layer soil temperature yesterday (from 
                                                      !!    topSoilTemp_sum()) [K]
     REAL(dp),POINTER :: alpha_yDayMean(:,:)          !! mean value of water stress coefficient alpha yesterday (from alpha_sum())
     REAL(dp),POINTER :: drainage_yDayMean(:)         !! mean value of drainage yesterday (from drainage_sum())

     REAL(dp),POINTER :: box_soil_respiration(:,:)    !! Daily updated soil respiration rate [mol(CO2)/m^2(grid box) s]
     REAL(dp),POINTER :: box_soil_respiration_pot(:,:) !! daily updated soil resp without N limitation
     REAL(dp),POINTER :: box_NPP_yDayMean(:,:)        !! Daily updated mean rate of potential NPP [mol(CO2)/m^2(grid box) s]
     REAL(dp),POINTER :: box_NPP_act_yDayMean(:,:)    !! Daily updated mean rate of actual NPP [mol(CO2)/m^2(grid box) s]
     REAL(dp),POINTER :: box_NPP_flux_correction(:,:) !! Daily updated flux correction for NPP [mol(CO2)/m^2(grid box) s]
     REAL(dp),POINTER :: box_GPP_yDayMean(:,:)        !! Daily updated mean GPP rate [mol(CO2)/m^2(grid box) s] 

     REAL(dp),POINTER :: box_root_exudates(:,:)       !! Same as root_exudates, but relative to grid box area 
                                                      !!    [mol(C)/m^2(grid box)s]
     REAL(dp),POINTER :: box_Cflux_herbivory(:,:)     !! C flux taken from vegetation pools by herbivory [mol(C)/m^2(grid box) s]
     REAL(dp),POINTER :: box_Cflux_herbivory_2_atm(:,:) !! C flux from herbivory going to the atmosphere [mol(C)/m^2(grid box) s]
     REAL(dp),POINTER :: box_Cflx_2_crop_harvest(:)   !! C flux from vegetation to crop harvest pool [mol(C)/m^2(grid box) s]
     REAL(dp),POINTER :: box_Cflx_crop_harvest_2_atm(:) !! C flux from crop harvest pool to atmosphere [mol(C)/m^2(grid box) s]
     REAL(dp),POINTER :: litter_flux(:,:)             !! Carbon flux from the vegetation to the litter pools [mol(C)/m^2(canopy) s]
     REAL(dp),POINTER :: box_litter_flux(:,:)         !! Same as litter_flux, but relative to grid box area [mol(C)/m^2(grid box) s]
     ! -----HW
     REAL(dp),POINTER :: leaf_shedding_debug(:,:)             !!
     REAL(dp),POINTER :: Cpool_green_minus_shed_leaves(:,:)             !!
     REAL(dp),POINTER :: excess_carbon_debug(:,:)             !!
     REAL(dp),POINTER :: litter_leaf(:,:)             !! Carbon flux from the leaf to the litter pools [mol(C)/m^2(canopy) s]
     REAL(dp),POINTER :: box_litter_leaf(:,:)         !! Same as litter_leaf, but relative to grid box area [mol(C)/m^2(grid box) s]
     real(dp),pointer :: leaf_shedding_rate(:,:) ! HW
     real(dp),pointer :: relative_extractable_water(:,:) ! HW
     real(dp),pointer :: tau_Cpool_woods(:,:) ! HW
     ! HW-----
     REAL(dp),POINTER :: box_Cpools_total(:)          !! Sum of all carbon pools [mol(C)/m^2(grid box)]
     REAL(dp),POINTER :: box_Cpools_total_old(:)      !! Array needed for carbon conservation test
     REAL(dp),POINTER :: box_Cpools_total_old2(:)     !! Array needed for carbon conservation test
     REAL(dp),POINTER :: testCconserv(:,:)            !! Array for testing carbon mass conservation (output should be zero)
     REAL(dp),POINTER :: jsbachCconserv(:)            !! Array to test carbon conservation of jsbach (should be zero)
     REAL(dp),POINTER :: CO2_fluxes_accu(:)           !! Array needed for carbon conservation test
     
     ! Anthropogenic pools summed for temporal averages [mol m-2(grid box)]
     REAL(dp),POINTER :: boxC_onSite_avg(:)               !! Carbon left on the ground by lcc
     REAL(dp),POINTER :: boxC_paper_avg(:)                !! Carbon stored in short/intermediate term anthro pool by lcc
     REAL(dp),POINTER :: boxC_construction_avg(:)         !! Carbon stored in long term anthro pool by lcc
     REAL(dp),POINTER :: boxC_onSite_harvest_avg(:)       !! Carbon left on the ground from harvest
     REAL(dp),POINTER :: boxC_paper_harvest_avg(:)        !! Carbon stored in short/intermediate term anthro pool from harvest
     REAL(dp),POINTER :: boxC_construction_harvest_avg(:) !! Carbon stored in long term anthro pool from harvest

  END TYPE cbalance_diag_type

  TYPE nbalance_type  ! contains the state variables for nitrogen of the cbalance module (dims: domain%nland x vegetation%ntiles)
     INTEGER          :: ntiles
     REAL(dp),POINTER :: Npool_green(:,:)             !! Organic N-pool for leaves, fine roots, vegetative organs and other green
                                                      !!    (living) parts of vegetation [mol(N)/m^2(canopy)]
     REAL(dp),POINTER :: Npool_mobile(:,:)            !! plant N-pool mobile: buffer for plant internal retranslocation of N 
                                                      !!    [mol(N)/m^2(canopy)]
     REAL(dp),POINTER :: Npool_woods(:,:)             !! Organic N-pool for stems, thick roots and other (dead) structural material
                                                      !!    of living plants [mol(N)/m^2(canopy)]
     REAL(dp),POINTER :: Npool_crop_harvest(:,:)      !! Organic N-pool for biomass harvested from crops [mol(N)/m^2(grid box)]
     REAL(dp),POINTER :: Npool_litter_green_ag(:,:)   !! Organic N-pool for above ground non-woody litter especially from leaves
                                                      !!    [mol(N)/m^2(canopy)]
     REAL(dp),POINTER :: Npool_litter_green_bg(:,:)   !! Organic N-pool for below ground non-woody litter especially from fine roots
                                                      !!    [mol(N)/m^2(canopy)]
     REAL(dp),POINTER :: Npool_litter_wood_ag(:,:)    !! Organic N-pool for litter from above ground woody parts of the plants 
                                                      !!    [mol(N)/m^2(canopy)]
     REAL(dp),POINTER :: Npool_litter_wood_bg(:,:)    !! Organic N-pool for litter from below ground woody parts of the plants
                                                      !!    [mol(N)/m^2(canopy)]
     REAL(dp),POINTER :: Npool_slow(:,:)              !! Organic N-pool for below ground organic material (coming from the litter
                                                      !!    pools) [mol(N)/m^2(canopy)]
     REAL(dp),POINTER :: SMINN_pool(:,:)              !! Mineral N-pool for all types of mineral N in soils [mol(N)/m^2(canopy)]
     REAL(dp),POINTER :: redFact_Nlimitation(:,:)     !! redFact_Nlimitation(no dimension)
     REAL(dp),POINTER :: minNflux_litter_green_ag(:,:)!! Mineral N-flux associated with above ground green litter
                                                      !!    [mol(N)/m^2(canopy)s]
     REAL(dp),POINTER :: minNflux_litter_green_bg(:,:)!! Mineral N-flux associated with below ground green litter
                                                      !!    [mol(N)/m^2(canopy)s]  
     REAL(dp),POINTER :: minNflux_litter_wood_ag(:,:) !! Mineral N-flux associated with above ground woody litter
                                                      !!    [mol(N)/m^2(canopy)s]
     REAL(dp),POINTER :: minNflux_litter_wood_bg(:,:) !! Mineral N-flux associated with below ground woody litter
                                                      !!    [mol(N)/m^2(canopy)s]
     REAL(dp),POINTER :: minNflux_slow(:,:)           !! Mineral N-flux from slow pool [mol(N)/m^2(canopy)s]
     REAL(dp),POINTER :: Nplant_demand(:,:)           !! N demand of plants [mol(N)/m^2(canopy)s]
     REAL(dp),POINTER :: Nsoil_demand(:,:)            !! N demand of soil for litter decomposition [mol(N)/m^2(canopy)s]
     REAL(dp),POINTER :: Ntotal_demand(:,:)           !! total N demand of plants and soil [mol(N)/m^2(canopy)s]
     REAL(dp),POINTER :: SMINN_herbivory(:,:)         !! SMIIN from herbivory via animal faeces [mol(N)/m^2(canopy)s]

     REAL(dp),POINTER :: nfix_to_sminn(:,:)           !! BNF [mol(N)/m^2(canopy)s]
     REAL(dp),POINTER :: ndep_to_sminn(:,:)           !! Deposition [mol (N)/m^2(canopy)s]
     REAL(dp),POINTER :: ndep_forc(:)                 !! Atmospheric N deposition forcing [kg(N) m-2(grid box) s-1]
     REAL(dp),POINTER :: nfert_forc(:)                !! N fertilizer forcing [g(N)/m^2(grid box)y]
     REAL(dp),POINTER :: nfert_forc_2d(:,:)           !! N fertiliser use in agri [g(N)/m^2(grid box) y] based on external forcing.
     REAL(dp),POINTER :: nfert_to_sminn(:,:)          !! N fertilizer [mol (N)/m^2(canopy)s]     
     REAL(dp),POINTER :: sminn_to_denit(:,:)          !! Dentrification [mol(N)/m^2(canopy)s]
     REAL(dp),POINTER :: sminn_leach(:,:)             !! Leach [mol(N)/m^2(canopy)s] 
     REAL(dp),POINTER :: N2O_emissions_depfix(:,:)    !! N2O emissions  [mol(N)/m^2(canopy)s] due to deposition and fixation
     REAL(dp),POINTER :: N2O_emissions_mineraliz(:,:) !! N2O emissions  [mol(N)/m^2(canopy)s] due to N-mineralization
     REAL(dp),POINTER :: N2O_emissions_nfert(:,:)     !! N2O emissions  [mol(N)/m^2(canopy)s] due to N-fertilizer use     
     REAL(dp),POINTER :: N2O_emissions_slow(:,:)      !!
     REAL(dp),POINTER :: N2O_emissions_grazing(:,:)   !! N2O emissions from N input by herbivores
     REAL(dp),POINTER :: Nflx_2_crop_harvest(:,:)     !! Flux of N from vegetation to crop harvest pool [mol(N)/m^2(canopy)s]
     REAL(dp),POINTER :: NetEcosyst_N_flux(:,:)       !! Total balance of N gains and losses (positive for ecosystem gain)
                                                      !!    [mol(N)/m^2(canopy)s]
     REAL(dp),POINTER :: NPP_run_mean(:,:)            !! Exponential running mean of NPP used for nitrogen fixation model
                                                      !!    [mol(C)/m^2(canopy)s]
  END TYPE nbalance_type

  TYPE nbalance_diag_type  !! contains diagnostic variables  (1.dim: domain%nland, 2.dim: vegetation%ntiles)
     REAL(dp),POINTER :: boxN_green(:,:)                 !! Npool_green() relative to grid box area [mol(N)/m^2(grid box)] 
     REAL(dp),POINTER :: boxN_mobile(:,:)                !! Npool_mobile() relative to grid box area [mol(N)/m^2(grid box)]
     REAL(dp),POINTER :: boxN_woods(:,:)                 !! Npool_woods() relative to grid box area [mol(N)/m^2(grid box)]
     REAL(dp),POINTER :: boxN_crop_harvest(:,:)          !! Npool_crop_harvest() relative to grid box area [mol(N)/m^2(grid box)]
     REAL(dp),POINTER :: boxN_litter_green_ag(:,:)       !! Npool_litter_green_ag() relative to grid box area [mol(N)/m^2(grid box)]
     REAL(dp),POINTER :: boxN_litter_green_bg(:,:)       !! Npool_litter_green_bg() relative to grid box area [mol(N)/m^2(grid box)]
     REAL(dp),POINTER :: boxN_litter_wood_ag(:,:)        !! Npool_litter_wood_ag() relative to grid box area [mol(N)/m^2(grid box)]
     REAL(dp),POINTER :: boxN_litter_wood_bg(:,:)        !! Npool_litter_wood_bg() relative to grid box area [mol(N)/m^2(grid box)]
     REAL(dp),POINTER :: boxN_slow(:,:)                  !! Npool_slow() relative to grid box area [mol(N)/m^2(grid box)
     REAL(dp),POINTER :: boxN_SMINN(:,:)                 !! SMINN_pool() relative to grid box area [mol(N)/m^2(grid box)
     REAL(dp),POINTER :: box_N2O_total(:,:)              !! Sum of all N2O fluxes [mikro-mol(N)/m^2(grid box)s]
     REAL(dp),POINTER :: box_N2O_emissions_depfix(:,:)   !! N2O emissions due to depos. and fixation [mikro-mol(N) m-2(grid box)s-1]
     REAL(dp),POINTER :: box_N2O_emissions_mineraliz(:,:)!! N2O emissions due to N-mineralization [mikro-mol(N) m-2(grid box) s-1]
     REAL(dp),POINTER :: box_minNflux_litter_green(:,:)  !! Mineral N-flux associated with green litter [mol(N) m-2(grid box) s-1]
     REAL(dp),POINTER :: box_minNflux_litter_wood(:,:)   !! Mineral N-flux associated with woody litter [mol(N) m-2(grid box) s-1]
     REAL(dp),POINTER :: box_minNflux_slow(:,:)          !! Mineral N-flux from slow pool [mol(N) m-2(grid box) s-1]
     REAL(dp),POINTER :: box_Nplant_demand(:,:)          !! N demand of plants  [mol(N) m-2(grid box) s-1]
     REAL(dp),POINTER :: box_Nsoil_demand(:,:)           !! N demand of soil for litter decomposition [mol(N) m-2(grid box) s-1]
     REAL(dp),POINTER :: box_Ntotal_demand(:,:)          !! total N demand of plants and soil [mol(N) m-2(grid box) s-1]
     REAL(dp),POINTER :: box_SMINN_herbivory(:,:)        !! SMINN from herbivory via animal faeces [mol(N) m-2(grid box) s-1]
     REAL(dp),POINTER :: box_nfix_to_sminn(:,:)          !! BNF [mol(N) m-2(grid box) s-1]
     REAL(dp),POINTER :: box_ndep_to_sminn(:,:)          !! Deposition [mol(N) m-2(grid box) s-1]
     REAL(dp),POINTER :: box_nfert_to_sminn(:,:)         !! N fertilizer [mol(N) m-2(grid box) s-1]
     REAL(dp),POINTER :: box_sminn_to_denit(:,:)         !! Dentrification [mol(N) m-2(grid box) s-1]
     REAL(dp),POINTER :: box_sminn_leach(:,:)            !! Leach [mol(N) m-2(grid box) s-1]
!! TR     REAL(dp),POINTER :: box_tot_sminn_gain(:,:)         !! Total N influx to soil mineral N pool [mol(N) m-2(grid box) s-1]
     REAL(dp),POINTER :: box_tot_sminn_loss(:,:)         !! Total N eflux of soil mineral N pool [mol(N) m-2(grid box) s-1]
     REAL(dp),POINTER :: box_Nflx_2_crop_harvest(:)      !! Flux of N from vegetation to crop harvest pool [mol(N) m-2(grid box) s-1]
     REAL(dp),POINTER :: testNconserv(:,:)               !! Array for testing nitrogen mass conservation (output should be zero)
     REAL(dp),POINTER :: box_Npools_total(:)             !! Sum of all nitrogen pools [mol(N)/m^2(grid box)]
     REAL(dp),POINTER :: box_Npools_total_old(:)         !! Array needed for nitrogen conservation test
     REAL(dp),POINTER :: box_Npools_total_old2(:)        !! Array needed for nitrogen conservation test
     REAL(dp),POINTER :: jsbachNconserv(:)               !! Array to test nitrogen conservation of jsbach (should be zero)
     REAL(dp),POINTER :: N2_fluxes_accu(:)               !! Array needed for nitrogen conservation test

  END TYPE nbalance_diag_type

  PUBLIC :: cbalance_type
  PUBLIC :: cbalance_diag_type
  PUBLIC :: cbalance_diag
  PUBLIC :: update_cbalance_bethy      ! This subroutine computes NPP
  PUBLIC :: init_cbalance_bethy        ! This subroutine initializes this module
  PUBLIC :: NPP_rate_bethy
  PUBLIC :: nbalance_type
  PUBLIC :: nbalance_diag_type
  PUBLIC :: check_C_conservation, check_N_conservation
  PUBLIC :: sumCnatural, sumCanthro, sumNpools

  ! === END OF PUBLIC PART === BEGIN OF PRIVATE PART ===============================================================================

  PRIVATE 

  PUBLIC :: read_nitrogen_deposition_fluxes, read_ndepo
  PUBLIC :: read_nitrogen_fertilizer_fluxes, read_nfert
  PUBLIC :: read_cpools

  !! parameters

  INTEGER, PARAMETER              :: sec_per_day         = 86400  !! seconds per day

  !! module variables

  LOGICAL, SAVE                   :: module_initialized = .FALSE. !! Signifies whether module has been initialized
  INTEGER, SAVE                   :: ntiles
  INTEGER, SAVE                   :: time_steps_per_day           !! Number of time steps per day

  TYPE(cbalance_diag_type), SAVE  :: cbalance_diag
  TYPE(nbalance_diag_type), SAVE  :: nbalance_diag

  TYPE(t_stream), POINTER, SAVE   :: IO_cbalance   ! Memory stream for cbalance model state
  TYPE(t_stream), POINTER, SAVE   :: IO_ybalance   ! Memory stream for yasso model state
  TYPE(t_stream), POINTER, SAVE   :: IO_nbalance   ! Memory stream for nitrogen balance model state

  TYPE(FILE_INFO),   SAVE         :: cpool_file
  REAL(dp), POINTER, SAVE         :: init_Cpool_green(:,:),init_Cpool_woods(:,:),init_Cpool_reserve(:,:)
  REAL(dp), POINTER, SAVE         :: init_Cpool_crop_harvest(:,:)
  REAL(dp), POINTER, SAVE         :: init_Cpool_litter_green_ag(:,:),init_Cpool_litter_green_bg(:,:)
  REAL(dp), POINTER, SAVE         :: init_Cpool_litter_wood_ag(:,:),init_Cpool_litter_wood_bg(:,:)
  REAL(dp), POINTER, SAVE         :: init_Cpool_slow(:,:)
  REAL(dp), POINTER, SAVE         :: init_YCpool_water_ag1(:,:),init_YCpool_ethanol_ag1(:,:),init_YCpool_acid_ag1(:,:)
  REAL(dp), POINTER, SAVE         :: init_YCpool_nonsoluble_ag1(:,:)
  REAL(dp), POINTER, SAVE         :: init_YCpool_water_bg1(:,:),init_YCpool_ethanol_bg1(:,:),init_YCpool_acid_bg1(:,:)
  REAL(dp), POINTER, SAVE         :: init_YCpool_nonsoluble_bg1(:,:),init_YCpool_humus_1(:,:) 
  REAL(dp), POINTER, SAVE         :: init_YCpool_water_ag2(:,:),init_YCpool_ethanol_ag2(:,:),init_YCpool_acid_ag2(:,:)
  REAL(dp), POINTER, SAVE         :: init_YCpool_nonsoluble_ag2(:,:)
  REAL(dp), POINTER, SAVE         :: init_YCpool_water_bg2(:,:),init_YCpool_ethanol_bg2(:,:),init_YCpool_acid_bg2(:,:)
  REAL(dp), POINTER, SAVE         :: init_YCpool_nonsoluble_bg2(:,:),init_YCpool_humus_2(:,:) 
  REAL(dp), POINTER, SAVE         :: init_frac_litter_wood_new(:,:)

 ! --- constants needed in "update_cbalance_bethy" (initialized in "init_cbalance_bethy") 
 !                                                             (P. Räisänen, 16 June 2008)    
  REAL(dp), SAVE          :: N_pseudo_temp    !! Normalization for computing pseudo-15-day mean air temperature
  REAL(dp), SAVE          :: F_pseudo_temp    !! exponential factor used for updating the pseudo-15-day mean temperature
  REAL(dp), SAVE          :: N_pseudo_precip  !! Normalization for computing pseudo-15-day mean precipitation rate
  REAL(dp), SAVE          :: F_pseudo_precip  !! exponential factor used for updating the pseudo-15-day mean precip

  TYPE(FILE_INFO), SAVE  :: npool_file              !! initial file for N pools (with read_npools=.TRUE.)
  REAL(dp), POINTER, SAVE         :: init_Npool_green(:,:),init_Npool_woods(:,:),init_Npool_mobile(:,:)
  REAL(dp), POINTER, SAVE         :: init_Npool_crop_harvest(:,:)
  REAL(dp), POINTER, SAVE         :: init_Npool_litter_green_ag(:,:),init_Npool_litter_green_bg(:,:)
  REAL(dp), POINTER, SAVE         :: init_Npool_litter_wood_ag(:,:),init_Npool_litter_wood_bg(:,:)
  REAL(dp), POINTER, SAVE         :: init_Npool_slow(:,:)
  REAL(dp), POINTER, SAVE         :: init_SMINN_pool(:,:), init_NPP_run_mean(:,:)
  REAL(dp), POINTER               :: init_Cpool_onSite              (:), &
                                     init_Cpool_paper               (:), &
                                     init_Cpool_construction        (:), &
                                     init_Cpool_onSite_harvest      (:), &
                                     init_Cpool_paper_harvest       (:), &
                                     init_Cpool_construction_harvest(:)

  ! namelist parameters
  LOGICAL                :: read_cpools  !! initialize C pools with values from cpool_file 
  LOGICAL                :: init_npools  !! initialize N pools consistent with C pools using C/N ratios
  LOGICAL                :: read_npools  !! initialize N pools with values from npool_file 
  LOGICAL                :: read_ndepo   !! read nitrogen deposition from external file
  LOGICAL                :: read_nfert   !! read nitrogen fertilization from external file
  CHARACTER(NF_MAX_NAME) :: cpool_file_name, npool_file_name, ndepo_file_name, nfert_file_name

CONTAINS

  SUBROUTINE init_cbalance_bethy(grid, domain, cbalance, nbalance, with_nitrogen, lcc_scheme, with_forest_management, &
                                 UseLandcoverTransitions, fileformat, fileztype, stream, Nstream, Ystream, with_yasso)

    USE mo_linked_list,            ONLY: LAND, TILES, NETCDF
    USE mo_output,                 ONLY: veg_table, nitrogen_table, yasso_table
    USE mo_memory_base,            ONLY: new_stream,default_stream_setting, &
                                         add =>add_stream_element
    USE mo_tr_scatter,             ONLY: scatter_field
    USE mo_netcdf,                 ONLY: MAX_DIM_NAME, io_inq_varid, io_get_var_double
    USE mo_cbal_cpools,            ONLY: printCpoolParameters, printNpoolParameters
    USE mo_cbal_parameters,        ONLY: config_cbal_parameters
    USE mo_jsbach,                 ONLY: missing_value, nml_unit
    USE mo_temp,                   ONLY: zreal2d, zreal3d, zreal2d_ptr
    USE mo_namelist,               ONLY: position_nml, POSITIONED
    USE mo_io_units,               ONLY: nout
    USE mo_time_control,           ONLY: lresume

    TYPE(grid_type),    INTENT(in)     :: grid
    TYPE(domain_type),  INTENT(in)     :: domain
    TYPE(cbalance_type),intent(inout)  :: cbalance     !! Pointer to cbalance structure (allocated in the calling routine)
                                                       !! .. that shall be initialized
    TYPE(nbalance_type),intent(inout)  :: nbalance
    LOGICAL,            INTENT(in)     :: with_nitrogen!! If .TRUE. initialise also for cycling of Nitrogen in addition to Carbon 
    LOGICAL,            INTENT(in)     :: with_yasso   !! If .TRUE. use Yasso model 
    LOGICAL,            INTENT(in)     :: with_forest_management
    LOGICAL,            INTENT(in)     :: UseLandcoverTransitions
    INTEGER,            INTENT(in)     :: lcc_scheme
    INTEGER,            INTENT(in)     :: fileformat   !! output file format
    INTEGER,            INTENT(in)     :: fileztype    !! output file compression
    TYPE(t_stream), POINTER, OPTIONAL  :: stream, NStream, Ystream
    
    ! local variables

    INTEGER :: IO_file_id, IO_var_id

    INTEGER :: g_nland, l_nland
    INTEGER                     :: dim1p(1), dim1(1)
    INTEGER                     :: dim3p(2), dim3(2)
    CHARACTER(LEN=max_dim_name) :: dim1n(1), dim3n(2)

    INTEGER :: i, read_status, f_unit
    LOGICAL :: with_N2O = .TRUE. 

    INCLUDE 'cbalance_ctl.inc'
    
    !! Check initialization

    IF (module_initialized) CALL finish("init_cbalance_bethy()","Attempt to initialize twice.")

    !! Read cbal parameters from namelist/set to defaults
    call config_cbal_parameters(nml_unit)

    !! Printout parameters 

    CALL printCpoolParameters
    IF (with_nitrogen) CALL printNpoolParameters

    !! set parameters

    ntiles = cbalance%ntiles

    IF (PRESENT(stream)) THEN 
       IF (.NOT. ASSOCIATED(stream)) THEN
          ! Add new stream
          CALL new_stream(stream, 'cbalance', filetype=fileformat, ztype=fileztype)
          ! Set default stream options
          CALL default_stream_setting(stream, table=veg_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
       ENDIF
       IO_cbalance => stream
    ELSE
       ! Add new stream
       CALL new_stream(IO_cbalance, 'cbalance', filetype=fileformat, ztype=fileztype)
       ! Set default stream options
       CALL default_stream_setting(IO_cbalance, table=veg_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
    ENDIF

    IF (PRESENT(Nstream)) THEN 
       IF (.NOT. ASSOCIATED(Nstream)) THEN
          ! Add new stream
          CALL new_stream(Nstream, 'nbalance', filetype=fileformat, ztype=fileztype)
          ! Set default stream options
          CALL default_stream_setting(Nstream, table=nitrogen_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
       ENDIF
       IO_nbalance => Nstream
    ELSE
       ! Add new stream
       CALL new_stream(IO_nbalance, 'nbalance', filetype=fileformat, ztype=fileztype)
       ! Set default stream options
       CALL default_stream_setting(IO_nbalance, table=nitrogen_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
    ENDIF

    IF (PRESENT(Ystream)) THEN 
       IF (.NOT. ASSOCIATED(Ystream)) THEN
          ! Add new stream
          CALL new_stream(Ystream, 'ybalance', filetype=fileformat, ztype=fileztype)
          ! Set default stream options
          CALL default_stream_setting(Ystream, table=yasso_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
       ENDIF
       IO_ybalance => Ystream
    ELSE
       ! Add new stream
       CALL new_stream(IO_ybalance, 'ybalance', filetype=fileformat, ztype=fileztype)
       ! Set default stream options
       CALL default_stream_setting(IO_ybalance, table=yasso_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
    ENDIF



    !! Read namelist cbalance_ctl

    IF (p_parallel_io) THEN

       ! define default values
       read_cpools = .FALSE.
       cpool_file_name = 'Cpools.nc'

       read_npools = .FALSE.
       npool_file_name = 'Npools.nc'
       init_npools = .FALSE.
       
       read_ndepo = .FALSE.
       ndepo_file_name = 'Ndepo.nc'

       read_nfert = .FALSE.
       nfert_file_name = 'Nfert.nc'

       f_unit = position_nml ('CBALANCE_CTL', nml_unit, status=read_status)
       SELECT CASE (read_status)
       CASE (POSITIONED)
          READ (f_unit, cbalance_ctl)
          CALL message('init_cbalance_bethy', 'Namelist CBALANCE_CTL: ')
          WRITE(nout, cbalance_ctl)
       END SELECT
    ENDIF

    IF (p_parallel_io) THEN
       IF (read_cpools) THEN
          CALL message('init_cbalance_bethy',&
                       ' Initial values for the carbon pools will be taken from file '//TRIM(cpool_file_name))
       ELSE
          CALL message('init_cbalance_bethy',' Carbon pools will be initialized with default values.')
       END IF
    END IF
    IF (p_parallel) THEN
       CALL p_bcast(read_cpools, p_io)
    END IF

    IF (p_parallel_io) THEN
       IF (with_nitrogen) then
          IF (read_npools) THEN
             CALL message('init_cbalance_bethy',&
                          ' Initial values for the nitrogen pools are read taken from file '//TRIM(npool_file_name))
          ELSE IF (init_npools) THEN
             CALL message('init_cbalance_bethy',' Nitrogen pools are initialized consistent with carbon pools.')
          END IF
          IF (read_ndepo) THEN
             CALL message('init_cbalance_bethy',&
                  ' Nitrogen deposition fluxes will be taken from file '//TRIM(ndepo_file_name))
          ELSE
             CALL message('init_cbalance_bethy',' No nitrogen deposition fluxes considered.')
          END IF
          IF (read_nfert) THEN
             CALL message('init_cbalance_bethy',&
                  ' Nitrogen fertilization will be taken from file '//TRIM(nfert_file_name))
          ELSE
             CALL message('init_cbalance_bethy',' No nitrogen fertilization inputs considered.')
          END IF
       ELSE
          CALL message('init_cbalance_bethy',' This is a run without nitrogen cycling.')
          read_npools = .FALSE.
          init_npools = .FALSE.
          read_ndepo  = .FALSE.
          read_nfert  = .FALSE.
       END IF
    END IF
    IF (p_parallel) THEN
       CALL p_bcast(read_npools, p_io)
       CALL p_bcast(init_npools, p_io)
       CALL p_bcast(read_ndepo, p_io)
       CALL p_bcast(read_nfert, p_io)
    END IF

    ! --- compute the number of time steps per day

    if( mod(86400,int(delta_time)) .ne. 0) then ! For computing the mean day temperature (see below) it is assumed that each day ..
                                                ! is computed with a fixed number of time steps. Therefore the program is  ..
                                                ! stopped wheen day length is not an integer multiple of the time step.
       call finish("init_cbalance_bethy","ERROR: Day length is not an integer multiple of the time step!")
    else 
       time_steps_per_day = 86400/int(delta_time)
    end if

    NULLIFY(cbalance%Cpool_onSite,               &
            cbalance%Cpool_paper,                &
            cbalance%Cpool_construction,         &
            cbalance%Cpool_onSite_harvest,       &
            cbalance%Cpool_paper_harvest,        &
            cbalance%Cpool_construction_harvest, &
            cbalance%frac_litter_wood_new)

    g_nland = grid%nland
    l_nland = domain%nland

    dim1p = (/ l_nland /)
    dim1  = (/ g_nland /)
    dim1n = (/ 'landpoint' /)

    dim3p = (/ l_nland, ntiles /)
    dim3  = (/ g_nland, ntiles /)
    dim3n(1) = 'landpoint'
    dim3n(2) = 'tiles'

    ! --- Define state variables as stream elements
    ! code numbers of this routine are 74-75, 130-179, 210-223 and 253-254 using the GRIB veg_table
    ! code numbers of this routine are 128-130 and 205-252 using the GRIB nitrogen_table;
    ! code numbers of this routine are 1-11 using the GRIB yasso_table
    CALL add(IO_cbalance,'Cpool_green'         ,cbalance%Cpool_green, lpost=debug .AND. fileformat==NETCDF,   &
             longname='C-Pool for Green Parts of Vegetation',     units='mol(C) m-2(canopy)',                 &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=150,lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'Cpool_woods'         ,cbalance%Cpool_woods, lpost=debug .AND. fileformat==NETCDF,   &
             longname='C-Pool for Structural Material of Plants', units='mol(C) m-2(canopy)',                 &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=151,lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'Cpool_reserve'     ,cbalance%Cpool_reserve, lpost=debug .AND. fileformat==NETCDF,   &
             longname='C-Pool for reserve carbohydrates (starches, sugars)', units='mol(C) m-2(canopy)',      &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=152,lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'Cpool_crop_harvest',cbalance%Cpool_crop_harvest, lpost=debug .AND. fileformat==NETCDF, &
             longname='C-Pool for Harvested Material from Crops', units='mol(C) m-2(canopy)',                 &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=213, lrerun=.TRUE., contnorest=.TRUE.,             &
             lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'Cpool_litter_green_ag', cbalance%Cpool_litter_green_ag,                             &
             lpost=debug .AND. fileformat==NETCDF .AND. .NOT. with_yasso,                                     &
             longname='C-Pool for above ground non-woody litter (leaves)',     units='mol(C) m-2(canopy)',    &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=177,lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'Cpool_litter_green_bg', cbalance%Cpool_litter_green_bg,                             &                     
             lpost=debug .AND. fileformat==NETCDF .AND. .NOT. with_yasso,                                     &
             longname='C-Pool for below ground non-woody litter (fine roots)', units='mol(C) m-2(canopy)',    &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=153,lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'Cpool_litter_wood_ag', cbalance%Cpool_litter_wood_ag,                               &
             lpost=debug .AND. fileformat==NETCDF .AND. .NOT. with_yasso,                                     &
             longname='C-Pool for above ground woody litter',     units='mol(C) m-2(canopy)',                 &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=148,lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'Cpool_litter_wood_bg', cbalance%Cpool_litter_wood_bg,                               &
             lpost=debug .AND. fileformat==NETCDF .AND. .NOT. with_yasso,                                     &
             longname='C-Pool for below ground woody litter',     units='mol(C) m-2(canopy)',                 &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=149,lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'Cpool_slow'          ,cbalance%Cpool_slow,                                    &
             lpost=debug .AND. fileformat==NETCDF .AND. .NOT. with_yasso,                               &
             longname='C-Pool for slowly respirated soil organic material', units='mol(C) m-2(canopy)', &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=154,lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'NPP_Rate'            ,cbalance%NPP_Rate,     lpost=.FALSE.,                   &
             longname='Net Primary Production Rate', units='mol(CO2) m-2(canopy) s-1',                  &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=155,lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'LAI_sum'             ,cbalance%LAI_sum,             ldims=dim3p, gdims=dim3,  &
             dimnames=dim3n, code=145, lrerun=.TRUE., lpost=.FALSE.)
    CALL add(IO_cbalance,'NPP_sum'             ,cbalance%NPP_sum,             ldims=dim3p, gdims=dim3,  &
             dimnames=dim3n, code=146, lrerun=.TRUE., lpost=.FALSE.)
    CALL add(IO_cbalance,'GPP_sum'             ,cbalance%GPP_sum,             ldims=dim3p, gdims=dim3,  &
             dimnames=dim3n, code=147, lrerun=.TRUE., lpost=.FALSE.)
    CALL add(IO_cbalance,'topSoilTemp_sum'     ,cbalance%topSoilTemp_sum,     ldims=dim1p, gdims=dim1,  &
             dimnames=dim1n, code=143, lrerun=.TRUE., lpost=.FALSE.)
    CALL add(IO_cbalance,'alpha_sum'           ,cbalance%alpha_sum,           ldims=dim3p, gdims=dim3,  &
             dimnames=dim3n, code=144, lrerun=.TRUE., lpost=.FALSE.)
    CALL add(IO_cbalance,'drainage_sum'          ,cbalance%drainage_sum,          ldims=dim1p, gdims=dim1,  &
             dimnames=dim1n, code=235, lrerun=.TRUE., contnorest=.TRUE., lpost=.FALSE.)

    CALL add(IO_cbalance,'soil_respiration', cbalance%soil_respiration, lpost=debug .AND. fileformat==NETCDF,   &
             longname='Soil respiration', units='mol(C) m-2(canopy)',                                           &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=156, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'soil_respiration_pot', cbalance%soil_respiration_pot, lpost=debug .AND. fileformat==NETCDF, &
             longname='Soil respiration without N limitation', units='mol(C) m-2(canopy)',                            &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=253, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'NPP_flux_correction' ,cbalance%NPP_flux_correction,                                   &
             longname='Flux correction for NPP',   units='mol(CO2) m-2(canopy) s-1',                            &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=157, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'excess_NPP' ,cbalance%excess_NPP,                                                     &
             longname='NPP that could not be not stored in vegetation carbon pools',   units='mol(CO2) m-2(canopy) s-1',&
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=130, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'root_exudates'    ,cbalance%root_exudates,                                            &
             longname='root exudates', units='mol(C) m-2(canopy) s-1',                                          &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=131, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'Cflux_herbivory'    ,cbalance%Cflux_herbivory,                                          &
             longname='Cflux_herbivory', units='mol(C) m-2(canopy) s-1',                                          &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=132, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'Cflux_herbivory_LG'    ,cbalance%Cflux_herbivory_LG,                                     &
             longname='Cflux_herbivory_LG', units='mol(C) m-2(canopy) s-1',                                        &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=133, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'Cflux_herbivory_2_atm'    ,cbalance%Cflux_herbivory_2_atm,                              &
             longname='Cflux_herbivory_2_atm', units='mol(C) m-2(canopy) s-1',                                    &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=134, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'Cflx_2_crop_harvest'   ,cbalance%Cflx_2_crop_harvest, contnorest=.TRUE.,                &
             longname='Cflx_2_crop_harvest', units='mol(C) m-2(canopy) s-1', ldims=dim3p, gdims=dim3,             &
             dimnames=dim3n, code=158, lrerun=.TRUE., lmiss=.TRUE., lpost=.FALSE., missval=missing_value)
    CALL add(IO_cbalance,'Cflx_crop_harvest_2_atm',cbalance%Cflx_crop_harvest_2_atm, contnorest=.TRUE.,           &
             longname='Cflx_crop_harvest_2_atm', units='mol(C) m-2(canopy) s-1', ldims=dim3p, gdims=dim3,         &
             dimnames=dim3n, code=211, lrerun=.TRUE., lmiss=.TRUE., lpost=.FALSE., missval=missing_value)

    IF (with_forest_management) THEN
       ! codes intersect with mo_landcover_change if using lcc_scheme=2;
       ! Grid-cell mass variables in m^2 grid box [mol(C)/m^2(grid box)]
       CALL add(IO_cbalance,'FOM_fm_Box_burning_pool', cbalance%FOM_fm_Box_burning_pool, lpost=.TRUE., &
         & longname='Grid-cell carbon in the burning pool', units='[mol(C)/m^2(grid box)]', lmiss=.TRUE., &
         & ldims=dim1p, gdims=dim1, dimnames=dim1n, code=225, contnorest=.TRUE., lrerun=.TRUE., missval=missing_value)
       CALL add(IO_cbalance,'FOM_fm_Box_paper_pool', cbalance%FOM_fm_Box_paper_pool, lpost=.TRUE., &
         & longname='Grid-cell carbon in the paper pool', units='[mol(C)/m^2(grid box)]', lmiss=.TRUE., &
         & ldims=dim1p, gdims=dim1, dimnames=dim1n, code=226, contnorest=.TRUE., lrerun=.TRUE., missval=missing_value)
       CALL add(IO_cbalance,'FOM_fm_Box_construction_pool', cbalance%FOM_fm_Box_construction_pool, lpost=.TRUE., &
         & longname='Grid-cell carbon in the construction pool', units='[mol(C)/m^2(grid box)]', lmiss=.TRUE., &
         & ldims=dim1p, gdims=dim1, dimnames=dim1n, code=227, contnorest=.TRUE., lrerun=.TRUE., missval=missing_value)
       CALL add(IO_cbalance,'FOM_fm_Box_eternity_pool', cbalance%FOM_fm_Box_eternity_pool, lpost=.TRUE., &
         & longname='Grid-cell carbon in the not decaying eternity pool', units='[mol(C)/m^2(grid box)]', lmiss=.TRUE., &
         & ldims=dim1p, gdims=dim1, dimnames=dim1n, code=228, contnorest=.TRUE., lrerun=.TRUE., missval=missing_value)
    ENDIF

    IF (lcc_scheme==2) THEN
       ! Anthropogenic lcc cpools on vegetated area
       CALL add(IO_cbalance,'Cpool_onSite_LCC', cbalance%Cpool_onSite,units='mol(C) m-2(veg)',                &
                longname='Carbon left on ground from land use change',code=212,                               &
                contnorest=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n, &
                lpost=.FALSE.,laccu=.FALSE.,lrerun=.TRUE.)
       CALL add(IO_cbalance,'Cpool_paper_LCC', cbalance%Cpool_paper,units='mol(C) m-2(veg)',code=214,         &
                longname='Wood carbon in short/intermediate term anthropogenic pool from land use change',    &
                contnorest=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n, &
                lpost=.FALSE.,laccu=.FALSE.,lrerun=.TRUE.)
       CALL add(IO_cbalance,'Cpool_construction_LCC', cbalance%Cpool_construction,units='mol(C) m-2(veg)',    &
                longname='Wood carbon in long term anthropogenic pool from land use change',code=215,         &
                contnorest=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n, &
                lpost=.FALSE.,laccu=.FALSE.,lrerun=.TRUE.)
       IF (UseLandcoverTransitions) THEN
          ! Anthropogenic harvest cpools on vegetated area
          CALL add(IO_cbalance,'Cpool_onSite_harvest', cbalance%Cpool_onSite_harvest,units='mol(C) m-2(veg)',    &
                   longname='Carbon left on ground from harvest',code=219,                                       &
                   contnorest=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n, &
                   lpost=.FALSE.,laccu=.FALSE.,lrerun=.TRUE.)
          CALL add(IO_cbalance,'Cpool_paper_harvest', cbalance%Cpool_paper_harvest,units='mol(C) m-2(veg)',      &
                   longname='Wood carbon in short/intermediate term anthropogenic pool from harvest',code=216,   &
                   contnorest=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n, &
                   lpost=.FALSE.,laccu=.FALSE.,lrerun=.TRUE.)
          CALL add(IO_cbalance,'Cpool_construction_harvest', cbalance%Cpool_construction_harvest, code=217,      &
                   units='mol(C) m-2(veg)', longname='Wood carbon in long term anthropogenic pool from harvest', &
                   contnorest=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n, &
                   lpost=.FALSE.,laccu=.FALSE.,lrerun=.TRUE.)
       ENDIF
    ENDIF
    CALL add(IO_cbalance,'frac_litter_wood_new', cbalance%frac_litter_wood_new, code=75,                         &
             units='', longname='fraction of new litter in woody litter pool',                                   &
             contnorest=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim3p,gdims=dim3,dimnames=dim3n,       &
             lpost=.TRUE.,laccu=.FALSE.,lrerun=.TRUE.)

    ! --- Define diagnostic variables as stream elements

    CALL add(IO_cbalance,'boxC_green',cbalance_diag%boxC_green, contnorest=.TRUE.,                              &
             longname='C-Pool for Green Parts of Vegetation',     units='mol(C) m-2(grid box)',                 &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=160, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'boxC_woods',cbalance_diag%boxC_woods, contnorest=.TRUE.,                              &
             longname='C-Pool for Structural Material of Plants', units='mol(C) m-2(grid box)',                 &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=161, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'boxC_reserve',cbalance_diag%boxC_reserve, contnorest=.TRUE.,                          &
             longname='C-Pool for reserve carbohydrates (starches, sugars)', units='mol(C) m-2(grid box)',      &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=162, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'boxC_crop_harvest',cbalance_diag%boxC_crop_harvest, contnorest=.TRUE.,                &
             longname='C-Pool for material harvested from crops', units='mol(C) m-2(grid box)',                 &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=225, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'boxC_litter_green_ag',cbalance_diag%boxC_litter_green_ag, contnorest=.TRUE.,          &
             lpost=.NOT.with_yasso,                                                                             &
             longname='C-Pool for above ground non-woody litter (leaves)', units='mol(C) m-2(grid box)',        &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=179, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'boxC_litter_green_bg',cbalance_diag%boxC_litter_green_bg, contnorest=.TRUE.,          &
             lpost=.NOT.with_yasso,                                                                             &
             longname='C-Pool for below ground non-woody litter (fine roots)', units='mol(C) m-2(grid box)',    &
             ldims=dim3p,gdims= dim3, dimnames=dim3n, code=163, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'boxC_litter_wood', cbalance_diag%boxC_litter_wood, contnorest=.TRUE.,                 &
             lpost=.NOT.with_yasso,                                                                             &
             longname='C-Pool for woody litter', units='mol(C) m-2(grid box)',                                  &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=159, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'boxC_slow',cbalance_diag%boxC_slow, contnorest=.TRUE.,                                &
             lpost=.NOT.with_yasso,                                                                             &
             longname='C-Pool for slowly respirated soil organic material',  units='mol(C) m-2(grid box)',      &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=164, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)

    CALL add(IO_cbalance,'LAI_previousDayMean',cbalance_diag%LAI_previousDayMean, contnorest=.TRUE.,  &
             longname='Mean Leaf Area Index of the Day Before Previous Day', units='', lpost=.FALSE., &
             ldims=dim3p,gdims= dim3, dimnames=dim3n, code=142, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'LAI_yDayMean',cbalance_diag%LAI_yDayMean,                            &
             longname='Mean Leaf Area Index of the Previous Day', units='', contnorest=.TRUE., &
             ldims=dim3p,gdims= dim3, dimnames=dim3n, code=165, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'NPP_yDayMean', cbalance_diag%NPP_yDayMean, contnorest=.TRUE.,        &
             longname='Mean NPP Rate of the Previous Day',   units='mol(CO2) m-2(canopy) s-1', &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=166, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'GPP_yDayMean', cbalance_diag%GPP_yDayMean, contnorest=.TRUE., lpost=.FALSE., &
             longname='Mean GPP Rate of the Previous Day',   units='mol(CO2) m-2(canopy) s-1', &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=167, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'topSoilTemp_yDayMean', cbalance_diag%topSoilTemp_yDayMean, contnorest=.TRUE., & 
             longname='Previous Day Mean Temperature of the Uppermost Soil Layer', units='K',  &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=168, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'alpha_yDayMean',cbalance_diag%alpha_yDayMean, contnorest=.TRUE.,     &
             longname='Previous Day Mean Value of the Water Stress Coefficient', units='',     &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=169, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance, 'drainage_yDayMean', cbalance_diag%drainage_yDayMean, contnorest=.TRUE.,               &
            longname='Previous Day Mean drainage', units='Kg m-2 s-1',                             &
            ldims=dim1p, gdims=dim1, dimnames=dim1n, code=248, lrerun=.TRUE.,  lmiss=.TRUE., missval=missing_value)
          
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
             !DSG: yasso forcing should be provided irrespectively if yasso is  !!
             !     on or off : but there is no free spot in the veg table       !!
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    IF (with_yasso) then
       CALL add(IO_ybalance, 'YCpool_acid_ag1', cbalance%YCpool_acid_ag1, lpost=debug .AND. fileformat==NETCDF,  &
              longname='C-Pool for acid-soluble litter in Yasso (aboveground)', units='mol(C) m-2(canopy)',      &
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=1, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'YCpool_acid_bg1', cbalance%YCpool_acid_bg1, lpost=debug .AND. fileformat==NETCDF,  &
              longname='C-Pool for acid-soluble litter in Yasso (belowground)', units='mol(C) m-2(canopy)',      &
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=2, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'YCpool_water_ag1', cbalance%YCpool_water_ag1, lpost=debug .AND. fileformat==NETCDF, &
              longname='C-Pool for water-soluble litter in Yasso (aboveground)', units='mol(C) m-2(canopy)',      &
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=3, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'YCpool_water_bg1', cbalance%YCpool_water_bg1, lpost=debug .AND. fileformat==NETCDF, &
              longname='C-Pool for water-soluble litter in Yasso (belowground)', units='mol(C) m-2(canopy)',      &
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=4, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'YCpool_ethanol_ag1', cbalance%YCpool_ethanol_ag1, lpost=debug .AND. fileformat==NETCDF, &
              longname='C-Pool for ethanol-soluble litter in Yasso (aboveground)', units='mol(C) m-2(canopy)',        &
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=5, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'YCpool_ethanol_bg1', cbalance%YCpool_ethanol_bg1, lpost=debug .AND. fileformat==NETCDF, &
              longname='C-Pool for ethanol-soluble litter in Yasso (belowground)', units='mol(C) m-2(canopy)',        &
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=6, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'YCpool_nonsoluble_ag1', cbalance%YCpool_nonsoluble_ag1, lpost=debug .AND. fileformat==NETCDF, &
              longname='C-Pool for nonsoluble litter in Yasso (aboveground)', units='mol(C) m-2(canopy)',                   &
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=7, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'YCpool_nonsoluble_bg1', cbalance%YCpool_nonsoluble_bg1, lpost=debug .AND. fileformat==NETCDF, &
              longname='C-Pool for nonsoluble litter in Yasso (belowground)', units='mol(C) m-2(canopy)',                   &
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=8, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'YCpool_humus_1', cbalance%YCpool_humus_1, lpost=debug .AND. fileformat==NETCDF,    &
              longname='C-Pool for humus in Yasso', units='mol(C) m-2(canopy)',                                  &
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=9, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)

       CALL add(IO_ybalance, 'YCpool_acid_ag2', cbalance%YCpool_acid_ag2, lpost=debug .AND. fileformat==NETCDF,  &
              longname='C-Pool for acid-soluble litter in Yasso (aboveground)', units='mol(C) m-2(canopy)',      &
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=21, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'YCpool_acid_bg2', cbalance%YCpool_acid_bg2, lpost=debug .AND. fileformat==NETCDF,  &
              longname='C-Pool for acid-soluble litter in Yasso (belowground)', units='mol(C) m-2(canopy)',      &
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=22, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'YCpool_water_ag2', cbalance%YCpool_water_ag2, lpost=debug .AND. fileformat==NETCDF, &
              longname='C-Pool for water-soluble litter in Yasso (aboveground)', units='mol(C) m-2(canopy)',      &
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=23, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'YCpool_water_bg2', cbalance%YCpool_water_bg2, lpost=debug .AND. fileformat==NETCDF, &
              longname='C-Pool for water-soluble litter in Yasso (belowground)', units='mol(C) m-2(canopy)',      &
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=24, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'YCpool_ethanol_ag2', cbalance%YCpool_ethanol_ag2, lpost=debug .AND. fileformat==NETCDF, &
              longname='C-Pool for ethanol-soluble litter in Yasso (aboveground)', units='mol(C) m-2(canopy)',        &
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=25, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'YCpool_ethanol_bg2', cbalance%YCpool_ethanol_bg2, lpost=debug .AND. fileformat==NETCDF, &
              longname='C-Pool for ethanol-soluble litter in Yasso (belowground)', units='mol(C) m-2(canopy)',        &
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=26, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'YCpool_nonsoluble_ag2', cbalance%YCpool_nonsoluble_ag2, lpost=debug .AND. fileformat==NETCDF, &
              longname='C-Pool for nonsoluble litter in Yasso (aboveground)', units='mol(C) m-2(canopy)',                   &
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=27, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'YCpool_nonsoluble_bg2', cbalance%YCpool_nonsoluble_bg2, lpost=debug .AND. fileformat==NETCDF, &
              longname='C-Pool for nonsoluble litter in Yasso (belowground)', units='mol(C) m-2(canopy)',                   &
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=28, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'YCpool_humus_2', cbalance%YCpool_humus_2, lpost=debug .AND. fileformat==NETCDF,    &
              longname='C-Pool for humus in Yasso', units='mol(C) m-2(canopy)',                                  &
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=29, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)

       CALL add(IO_ybalance, 'boxYC_acid_ag1', cbalance_diag%boxYC_acid_ag1,                                     &
              longname='C-Pool for acid-soluble litte in Yasso (aboveground)', units='mol(C) m-2(grid box)',     &  
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=31, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'boxYC_acid_bg1', cbalance_diag%boxYC_acid_bg1,                                     &
              longname='C-Pool for acid-soluble litte in Yasso (belowground)', units='mol(C) m-2(grid box)',     &  
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=32, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'boxYC_water_ag1', cbalance_diag%boxYC_water_ag1,                                   &
              longname='C-Pool for water-soluble litter in Yasso (aboveground)', units='mol(C) m-2(grid box)',   &  
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=33, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'boxYC_water_bg1', cbalance_diag%boxYC_water_bg1,                                   &
              longname='C-Pool for water-soluble litter in Yasso (belowground)', units='mol(C) m-2(grid box)',   &  
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=34, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'boxYC_ethanol_ag1', cbalance_diag%boxYC_ethanol_ag1,                               &
              longname='C-Pool for ethanol-soluble litter in Yasso (aboveground)', units='mol(C) m-2(grid box)', & 
              ldims=dim3p,gdims= dim3, dimnames=dim3n, code=35, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)  
       CALL add(IO_ybalance, 'boxYC_ethanol_bg1', cbalance_diag%boxYC_ethanol_bg1,                               &
              longname='C-Pool for ethanol-soluble litter in Yasso (belowground)', units='mol(C) m-2(grid box)', & 
              ldims=dim3p,gdims= dim3, dimnames=dim3n, code=36, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)  
       CALL add(IO_ybalance, 'boxYC_nonsoluble_ag1', cbalance_diag%boxYC_nonsoluble_ag1,                         &
              longname='C-Pool for nonsoluble litter in Yasso (aboveground)',  units='mol(C) m-2(grid box)',     & 
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=37, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'boxYC_nonsoluble_bg1', cbalance_diag%boxYC_nonsoluble_bg1,                         &
              longname='C-Pool for nonsoluble litter in Yasso (belowground)',  units='mol(C) m-2(grid box)',     & 
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=38, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'boxYC_humus_1', cbalance_diag%boxYC_humus_1,                                       &
              longname='C-Pool for humus in Yasso',                            units='mol(C) m-2(grid box)',     & 
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=39, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)

     
       CALL add(IO_ybalance, 'boxYC_acid_ag2', cbalance_diag%boxYC_acid_ag2,                                     &
              longname='C-Pool for acid-soluble litte in Yasso (aboveground)', units='mol(C) m-2(grid box)',     &  
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=41, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'boxYC_acid_bg2', cbalance_diag%boxYC_acid_bg2,                                     &
              longname='C-Pool for acid-soluble litte in Yasso (belowground)', units='mol(C) m-2(grid box)',     &  
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=42, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'boxYC_water_ag2', cbalance_diag%boxYC_water_ag2,                                   &
              longname='C-Pool for water-soluble litter in Yasso (aboveground)', units='mol(C) m-2(grid box)',   &  
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=43, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'boxYC_water_bg2', cbalance_diag%boxYC_water_bg2,                                   &
              longname='C-Pool for water-soluble litter in Yasso (belowground)', units='mol(C) m-2(grid box)',   &  
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=44, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'boxYC_ethanol_ag2', cbalance_diag%boxYC_ethanol_ag2,                               &
              longname='C-Pool for ethanol-soluble litter in Yasso (aboveground)', units='mol(C) m-2(grid box)', & 
              ldims=dim3p,gdims= dim3, dimnames=dim3n, code=45, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)  
       CALL add(IO_ybalance, 'boxYC_ethanol_bg2', cbalance_diag%boxYC_ethanol_bg2,                               &
              longname='C-Pool for ethanol-soluble litter in Yasso (belowground)', units='mol(C) m-2(grid box)', & 
              ldims=dim3p,gdims= dim3, dimnames=dim3n, code=46, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)  
       CALL add(IO_ybalance, 'boxYC_nonsoluble_ag2', cbalance_diag%boxYC_nonsoluble_ag2,                         &
              longname='C-Pool for nonsoluble litter in Yasso (aboveground)',  units='mol(C) m-2(grid box)',     & 
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=47, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'boxYC_nonsoluble_bg2', cbalance_diag%boxYC_nonsoluble_bg2,                         &
              longname='C-Pool for nonsoluble litter in Yasso (belowground)',  units='mol(C) m-2(grid box)',     & 
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=48, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'boxYC_humus_2', cbalance_diag%boxYC_humus_2,                                       &
              longname='C-Pool for humus in Yasso',                            units='mol(C) m-2(grid box)',     & 
              ldims=dim3p, gdims=dim3, dimnames=dim3n, code=49, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'pseudo_temp', cbalance%pseudo_temp, lpost=.FALSE.,                            & 
              longname='Pseudo-mean air temperature over previous 15 days', units='K', contnorest=.TRUE.,   &
              ldims=dim1p, gdims=dim1, dimnames=dim1n, code=19, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'pseudo_precip', cbalance%pseudo_precip, lpost=.FALSE., contnorest=.TRUE.,     & 
              longname='Pseudo-mean precipitation rate over previous 15 days', units='Kg m-2 s-1',          &  
              ldims=dim1p, gdims=dim1, dimnames=dim1n, code=20, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'pseudo_temp_yDay', cbalance_diag%pseudo_temp_yDay, contnorest=.TRUE.,         & 
              longname='Pseudo-mean air temperature over previous 15 days', units='K',                      &
              ldims=dim1p, gdims=dim1, dimnames=dim1n, code=50, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_ybalance, 'pseudo_precip_yDay', cbalance_diag%pseudo_precip_yDay, contnorest=.TRUE.,     & 
              longname='Pseudo-mean precipitation rate over previous 15 days', units='Kg m-2 s-1',          &  
              ldims=dim1p, gdims=dim1, dimnames=dim1n, code=51, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    END IF

    CALL add(IO_cbalance,'box_soil_respiration', cbalance_diag%box_soil_respiration, contnorest=.TRUE.,     &
             longname='Soil respiration', units='mol(C) m-2(grid box) s-1',                                 &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=170, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'box_soil_respiration_pot', cbalance_diag%box_soil_respiration_pot, lpost=with_nitrogen, &
             longname='Soil respiration without N limitation', units='mol(C) m-2(grid box) s-1',                  &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=254, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'box_NPP_yDayMean', cbalance_diag%box_NPP_yDayMean, contnorest=.TRUE.,             &
             longname='Mean potential NPP of the Previous Day',   units='mol(CO2) m-2(grid box) s-1',       &
             ldims=dim3p,gdims= dim3, dimnames=dim3n, code=171, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'box_NPP_act_yDayMean', cbalance_diag%box_NPP_act_yDayMean, contnorest=.TRUE.,     &
             longname='Mean actual NPP the Previous Day',   units='mol(CO2) m-2(grid box) s-1',             &
             ldims=dim3p,gdims= dim3, dimnames=dim3n, code=178, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'box_NPP_flux_correction', cbalance_diag%box_NPP_flux_correction,                  &
             longname='Flux correction for NPP',   units='mol(CO2) m-2(grid box) s-1', contnorest=.TRUE.,   &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=172, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'box_GPP_yDayMean', cbalance_diag%box_GPP_yDayMean, contnorest=.TRUE.,             &
             longname='Mean GPP Rate of the Previous Day',   units='mol(CO2) m-2(grid box) s-1',            &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=173, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'box_root_exudates', cbalance_diag%box_root_exudates, contnorest=.TRUE.,           &
             longname='root exudates', units='mol(C) m-2(grid box) s-1',                                    &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=135, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'box_Cflux_herbivory', cbalance_diag%box_Cflux_herbivory, contnorest=.TRUE.,       &
             longname='Cflux_herbivory', units='mol(C) m-2(grid box) s-1',                                  &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=136, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'box_Cflux_herbivory_2_atm', cbalance_diag%box_Cflux_herbivory_2_atm, contnorest=.TRUE., &
             longname='Carbon flux from herbivory to the atmsphere', units='mol(C) m-2(grid box) s-1',      &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=143, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'box_Cflx_2_crop_harvest', cbalance_diag%box_Cflx_2_crop_harvest, contnorest=.TRUE., &
             longname='Carbon flux from crops to crop harvest pool', units='mol(C) m-2(grid box) s-1',      &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=158, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'box_Cflx_crop_harvest_2_atm', cbalance_diag%box_Cflx_crop_harvest_2_atm, contnorest=.TRUE., &
             longname='Carbon flux from crop harvest pool to the atmsphere', units='mol(C) m-2(grid box) s-1', &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=211, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'litter_flux',cbalance_diag%litter_flux, lpost=debug .AND. fileformat==NETCDF,     &
             longname='Total litter flux entering the soil pools', units='mol(CO2) m-2(canopy) s-1', contnorest=.TRUE., &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=174, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'box_litter_flux', cbalance_diag%box_litter_flux, contnorest=.TRUE.,          &
             longname='Total litter flux entering the soil pools', units='mol(CO2) m-2(grid box) s-1', &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=175, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    ! -----HW
    CALL add(IO_cbalance,'leaf_shedding_debug',cbalance_diag%leaf_shedding_debug, contnorest=.true. ,     &
             longname='leaf_shedding_debug', units='',  &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=212, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'Cpool_green_minus_green_litter',cbalance_diag%Cpool_green_minus_shed_leaves, contnorest=.true. ,  &
             longname='Cpool_green_minus_green_litter', units='',  &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=213, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'excess_carbon_debug',cbalance_diag%excess_carbon_debug, contnorest=.true. ,     &
             longname='excess_carbon_debug', units='',  &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=214, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'litter_green',cbalance_diag%litter_leaf, contnorest=.true. ,     &
             longname='green litter flux entering the soil pools', units='mol(CO2) m-2(canopy) s-1',  &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=215, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'box_litter_green', cbalance_diag%box_litter_leaf, contnorest=.TRUE.,          &
             longname='green litter flux entering the soil pools', units='mol(CO2) m-2(grid box) s-1', &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=179, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'leaf_shedding_rate', cbalance_diag%leaf_shedding_rate, contnorest=.TRUE.,          &
             longname='leaf_shedding_rate', units='', &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=180, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'relative_extractable_water', cbalance_diag%relative_extractable_water, contnorest=.TRUE.,          &
             longname='relative_extractable_water', units='', &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=177, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance,'tau_Cpool_woods', cbalance_diag%tau_Cpool_woods, contnorest=.TRUE.,          &
             longname='tau_Cpool_woods', units='', &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=174, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    ! HW-----
    CALL add(IO_cbalance,'box_Cpools_total',cbalance_diag%box_Cpools_total, contnorest=.TRUE.,              &
             longname='Sum of carbon from all carbon pools',   units='mol(CO2) m-2(grid box)',              &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=176, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)

    IF (debug_Cconservation) THEN
       CALL add(IO_cbalance,'testCconserv',cbalance_diag%testCconserv,                                 &
             longname='For testing conservation of carbon mass: should be zero.',   units='mol(CO2) m-2(canopy) day-1',&
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=137, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
    END IF
    IF (test_Cconservation) THEN
       CALL add(IO_cbalance,'jsbachCconserv',cbalance_diag%jsbachCconserv, laccu=.FALSE., contnorest=.TRUE., &
             longname='jsbach carbon conservation test: should be zero.',   units='mol(CO2) m-2(grid box)', &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=138, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_cbalance,'box_Cpools_total_old',cbalance_diag%box_Cpools_total_old, contnorest=.TRUE., &
             longname='array needed for carbon conservation test', units='mol(CO2) m-2(grid box)', lpost=.FALSE., &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=140, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_cbalance,'box_Cpools_total_old2',cbalance_diag%box_Cpools_total_old2, contnorest=.TRUE., &
             longname='array needed for carbon conservation test', units='mol(CO2) m-2(grid box)', lpost=.FALSE., &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=141, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_cbalance,'CO2_fluxes_accu',cbalance_diag%CO2_fluxes_accu, laccu=.FALSE., lpost=debug_Cconservation, &
             longname='CO2 fluxes, accumulated over a day',   units='mol(CO2) m-2(grid box)', contnorest=.TRUE., &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=139, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    END IF

    ! --- Define nitrogen state variables as stream elements

    CALL add(IO_cbalance,'NPP_act_yDayMean'          ,cbalance%NPP_act_yDayMean,                        &
         longname='NPP_act_yDayMean', units='mol(CO2) m-2(canopy) s-1',                                 &
         ldims=dim3p, gdims=dim3, dimnames=dim3n, code=210,lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
    CALL add(IO_cbalance, 'redFact_Nlimitation', nbalance%redFact_Nlimitation, contnorest=init_npools.OR.read_npools, &
         longname='redFact_Nlimitation', units='dimensionless)', lpost=with_nitrogen,                                 &
         ldims=dim3p, gdims=dim3, dimnames=dim3n, code=212, lrerun=with_nitrogen, lmiss=.TRUE., missval=missing_value)

    IF (with_nitrogen) then
       CALL add(IO_nbalance, 'Npool_green',         nbalance%Npool_green,     contnorest=init_npools.OR.read_npools,         &
            longname='N-Pool for Green Parts of Vegetation',     units='mol(N) m-2(canopy)',                                 &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=213, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value,           &
            lpost=debug .AND. fileformat==NETCDF)
       CALL add(IO_nbalance, 'Npool_woods',         nbalance%Npool_woods,     contnorest=init_npools.OR.read_npools,         &
            longname='N-Pool for Structural Material of Plants', units='mol(N) m-2(canopy)',                                 &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=214, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value,           &
            lpost=debug .AND. fileformat==NETCDF)
       CALL add(IO_nbalance, 'Npool_mobile',        nbalance%Npool_mobile,    contnorest=init_npools.OR.read_npools,         &
            longname='N-Pool for mobile N',                      units='mol(N) m-2(canopy)',                                 &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=215, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value,           &
            lpost=debug .AND. fileformat==NETCDF)
       CALL add(IO_nbalance, 'Npool_crop_harvest',  nbalance%Npool_crop_harvest, contnorest=init_npools.OR.read_npools,      &
            longname='N-Pool for Harvested Material from Crops', units='mol(N) m-2(canopy)',                                 &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=235, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value,           &
            lpost=debug .AND. fileformat==NETCDF)
       CALL add(IO_nbalance, 'Npool_litter_green_ag', nbalance%Npool_litter_green_ag, contnorest=init_npools.OR.read_npools, &
            longname='N-Pool for above ground non-woody litter (leaves)', units='mol(N) m-2(canopy)',                        &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=216, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value,           &
            lpost=debug .AND. fileformat==NETCDF)
       CALL add(IO_nbalance, 'Npool_litter_green_bg', nbalance%Npool_litter_green_bg, contnorest=init_npools.OR.read_npools, &
            longname='N-Pool for below ground non-woody litter (fine roots)', units='mol(N) m-2(canopy)',                    &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=217, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value,           &
            lpost=debug .AND. fileformat==NETCDF)
       CALL add(IO_nbalance, 'Npool_litter_wood_ag', nbalance%Npool_litter_wood_ag, contnorest=init_npools.OR.read_npools,   &
            longname='N-Pool for above ground woody litter',     units='mol(N) m-2(canopy)',                                 &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=218, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value,           &
            lpost=debug .AND. fileformat==NETCDF)
       CALL add(IO_nbalance, 'Npool_litter_wood_bg', nbalance%Npool_litter_wood_bg, contnorest=init_npools.OR.read_npools,   &
            longname='N-Pool for below ground woody litter',     units='mol(N) m-2(canopy)',                                 &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=219, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value,           &
            lpost=debug .AND. fileformat==NETCDF)
       CALL add(IO_nbalance, 'Npool_slow',          nbalance%Npool_slow,      contnorest=init_npools.OR.read_npools,         &
            longname='N-Pool for slowly respirated soil organic material', units='mol(N) m-2(canopy)',                       &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=220, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value,           &
            lpost=debug .AND. fileformat==NETCDF)
       CALL add(IO_nbalance, 'sminN_pool',          nbalance%SMINN_pool,      contnorest=init_npools.OR.read_npools,         &
            longname='N-Pool for soil mineral nitrogen',         units='mol(N) m-2(canopy)',                                 &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=221,lrerun=.TRUE., lmiss=.TRUE., missval=missing_value,            &
            lpost=debug .AND. fileformat==NETCDF)

       CALL add(IO_nbalance, 'minNflux_litter_green_ag', nbalance%minNflux_litter_green_ag, contnorest=init_npools.OR.read_npools, &
            longname='Mineral N-flux associated with above ground green litter', units='mol(N) m-2(canopy) s-1',                   &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=223, lpost=.FALSE., lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_nbalance, 'minNflux_litter_green_bg', nbalance%minNflux_litter_green_bg, contnorest=init_npools.OR.read_npools, &
            longname='Mineral N-flux associated with above ground woody litter', units='mol(N) m-2(canopy) s-1 ',                  &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=224, lpost=.FALSE., lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_nbalance, 'minNflux_litter_wood_ag', nbalance%minNflux_litter_wood_ag, contnorest=init_npools.OR.read_npools,   &
            longname='Soil mineral Nitrogen flux associated with above ground woody litter', units='mol(N) m-2(canopy)  s-1',      &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=225, lpost=.FALSE., lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_nbalance, 'minNflux_litter_wood_bg', nbalance%minNflux_litter_wood_bg,  contnorest=init_npools.OR.read_npools,  &
            longname='Soil mineral Nitrogen flux associated with below ground woody litter', units='mol(N) m-2(canopy)  s-1',      &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=226, lpost=.FALSE., lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_nbalance, 'minNflux_slow',       nbalance%minNflux_slow,   contnorest=init_npools.OR.read_npools,               &
            longname='Mineral N-flux from slow N-pool',          units='mol(N) m-2(canopy)  s-1',                                  &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=227, lpost=.FALSE., lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       ! Note: Nplant_demand, Nsoil_demand and Ntotal_demand are multiplied with redFact_Nlimitation in update_Cpools (=> _uptake)
       CALL add(IO_nbalance, 'Nplant_uptake',       nbalance%Nplant_demand,   contnorest=.TRUE.,                                   &
            longname='N uptake by plants',                       units='mol(N) m-2(canopy) s-1 ',                                  &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=228, lpost=.FALSE., lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_nbalance, 'Nsoil_uptake',        nbalance%Nsoil_demand,    contnorest=.TRUE.,                                   &
            longname='N uptake for litter decomposition ',       units='mol(N) m-2(canopy) s-1 ',                                  &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=229, lpost=.FALSE., lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_nbalance, 'Ntotal_uptake',       nbalance%Ntotal_demand,   contnorest=.TRUE.,                                   &
            longname='N uptake of plants and for litter decomposition', units='mol(N) m-2(canopy) s-1 ',                           &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=230, lpost=.FALSE., lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)

       CALL add(IO_nbalance, 'Nfix_to_sminN',       nbalance%nfix_to_sminn,   contnorest=init_npools.OR.read_npools,               &
            longname='N flux to mineral N pool due to N fixation (BNF)', units='mol(N) m-2(canopy) s-1',                           &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=231, lpost=.FALSE., lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_nbalance, 'Ndep_to_sminN',       nbalance%ndep_to_sminn,   contnorest=init_npools.OR.read_npools,               &
            longname='N flux to mineral N pool due to N deposition', units='mol (N) m-2(canopy) s-1',                              &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=232, lpost=.FALSE., lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_nbalance, 'sminN_to_denit',      nbalance%sminn_to_denit,  contnorest=init_npools.OR.read_npools,               &
            longname='N flux from mineral N pool due to denitrification', units='mol(N) m-2(canopy) s-1',                          &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=233, lpost=.FALSE., lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_nbalance, 'sminN_leach',         nbalance%sminn_leach,     contnorest=init_npools.OR.read_npools,               &
            longname='N flux to mineral N pool due to leaching', units='mol(N) m-2(canopy)s-1',                                    &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=234, lpost=.FALSE., lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)

       CALL add(IO_nbalance, 'Ndep_forc',           nbalance%ndep_forc,       contnorest=init_npools.OR.read_npools,               &
            longname='atmospheric N deposition',                 units='kg (N) m-2(grid box) s-1',                                 &
            ldims=dim1p, gdims=dim1, dimnames=dim1n, code=236, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_nbalance, 'Nfert_forc',          nbalance%nfert_forc,      contnorest=init_npools.OR.read_npools,               &
            longname='N fertilizer forcing',                     units='mg (N) m-2(cropland) y-1',                                 &
            ldims=dim1p, gdims=dim1, dimnames=dim1n, code=209, lrerun=.TRUE., lmiss=.TRUE., lpost=.FALSE., missval=missing_value)
       CALL add(IO_nbalance, 'Nfert_forc_2d',       nbalance%nfert_forc_2d,   contnorest=init_npools.OR.read_npools,               &
            longname='N fertilizer forcing',                     units='mg (N) m-2(cropland) y-1 tile-1',                          &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=206, lrerun=.TRUE., lmiss=.TRUE., lpost=.FALSE., missval=missing_value) 
       CALL add(IO_nbalance, 'Nfert_to_sminN',      nbalance%nfert_to_sminn,  contnorest=init_npools.OR.read_npools,               &
            longname='N flux to mineral N pool due to fertilizers', units='mol (N) m-2(canopy) s-1',                               &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=208, lpost=.FALSE., lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)

       CALL add(IO_nbalance, 'NetEcosyst_N_flux',   nbalance%NetEcosyst_N_flux, contnorest=init_npools.OR.read_npools,             &
            longname='Total balance of N gains and losses (positive for ecosystem gain)', units='mol(N) m-2(canopy) s-1',          &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=237, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_nbalance, 'NPP_run_mean',        nbalance%NPP_run_mean,    contnorest=init_npools.OR.read_npools,               &
            longname='Exponential running mean of NPP used for N-fixation model', units='mol(N) m-2(canopy) s-1',                  &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=238, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)

       CALL add(IO_nbalance, 'sminN_herbivory',     nbalance%SMINN_herbivory, contnorest=init_npools.OR.read_npools,               &
            longname='N flux from herbivory via animal faeces to mineral N pool', units='mol(N) m-2(canopy) s-1',                  &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=239, lpost=.FALSE., lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)

       CALL add(IO_nbalance, 'Nflx_2_crop_harvest', nbalance%Nflx_2_crop_harvest, contnorest=.TRUE.,                               &
            longname='N flux from vegetation to crop harvest pool', units='mol(N) m-2(canopy) s-1',                                &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=158, lpost=.FALSE., lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)

       CALL add(IO_nbalance, 'N2O_emissions_depfix', nbalance%N2O_emissions_depfix, contnorest=init_npools.OR.read_npools,         &
            longname='N2O emissions due to deposition and fixation', units='mol(N) m-2(canopy)s-1', lpost=.FALSE.,                 &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=210, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_nbalance, 'N2O_emissions_mineraliz', nbalance%N2O_emissions_mineraliz, contnorest=init_npools.OR.read_npools,   &
            longname='N2O emissions due to N-mineralization', units='mol(N) m-2(canopy)s-1', lpost=.FALSE.,                        &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=211, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_nbalance, 'N2O_emissions_slow',  nbalance%N2O_emissions_slow, contnorest=init_npools.OR.read_npools,            &
            longname='N2O emission due to mineralization from slow pool', units='mol(N) m-2(canopy)s-1', lpost=.FALSE.,           &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=128,lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_nbalance, 'N2O_emissions_nfert', nbalance%N2O_emissions_nfert, contnorest=init_npools.OR.read_npools,           &
            longname='N2O emissions due to N fertilizer use', units='mol(N) m-2(canopy)s-1', lpost=.FALSE.,                        &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=207, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_nbalance, 'N2O_emissions_grazing', nbalance%N2O_emissions_grazing, contnorest=init_npools.OR.read_npools,       &
            longname='N2O emissions from N input by herbivores', units='mol(N) m-2(canopy)s-1', lpost=.FALSE.,                    &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=205,lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)

       ! --- Define nitrogen diagnostic variables as stream elements
       CALL add(IO_nbalance, 'box_N2O_emissions_depfix', nbalance_diag%box_N2O_emissions_depfix,             &
            longname='N2O emissions due to deposition and fixation', units='mikro-mol(N) m-2(grid box)', lpost=with_N2O, &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=251, lrerun=.TRUE., lmiss=.TRUE.,                  &
            contnorest=init_npools.OR.read_npools, missval=missing_value)
       CALL add(IO_nbalance, 'box_N2O_emissions_mineraliz', nbalance_diag%box_N2O_emissions_mineraliz,       &
            longname='N2O emissions due to mineralization', units='mikro-mol(N) m-2(grid box)', lpost=with_N2O,    &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=252, lrerun=.TRUE., lmiss=.TRUE.,                  &
            contnorest=init_npools.OR.read_npools, missval=missing_value)
       CALL add(IO_nbalance, 'box_N2O_total', nbalance_diag%box_N2O_total,                                   &
            longname='Sum of all N2O fluxes', units='mikro-mol(N) m-2(grid box) s-1', lpost=with_N2O,        &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=212, lrerun=.FALSE., lmiss=.TRUE.,                 &
            contnorest=init_npools.OR.read_npools, missval=missing_value)
       CALL add(IO_nbalance, 'boxN_green', nbalance_diag%boxN_green,                                         &
            longname='N-Pool for Green Parts of Vegetation',     units='mol(N) m-2(grid box)',               &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=240, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_nbalance, 'boxN_woods', nbalance_diag%boxN_woods,                                         &
            longname='N-Pool for Structural Material of Plants', units='mol(N) m-2(grid box)',               &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=241, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_nbalance, 'boxN_mobile', nbalance_diag%boxN_mobile,                                       &
            longname='N-Pool for mobile N', units='mol(N) m-2(grid box)',                                    &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=242, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_nbalance, 'boxN_crop_harvest', nbalance_diag%boxN_crop_harvest,                           &
            longname='N-Pool for material harvested from crops', units='mol(N) m-2(grid box)',               &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=222, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_nbalance, 'boxN_litter_green_ag', nbalance_diag%boxN_litter_green_ag,                     &
            longname='N-Pool for above ground non-woody litter (leaves)', units='mol(N) m-2(grid box)',      &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=243, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_nbalance, 'boxN_litter_green_bg', nbalance_diag%boxN_litter_green_bg,                     &
            longname='N-Pool for below ground non-woody litter (fine roots)', units='mol(N) m-2(grid box)',  &
            ldims=dim3p,gdims= dim3, dimnames=dim3n, code=244, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_nbalance, 'boxN_litter_wood_ag', nbalance_diag%boxN_litter_wood_ag,                       &
            longname='N-Pool for above ground woody litter', units='mol(N) m-2(grid box)',                   &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=245, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_nbalance, 'boxN_litter_wood_bg', nbalance_diag%boxN_litter_wood_bg,                       &
            longname='N-Pool for below ground woody litter', units='mol(N) m-2(grid box)',                   &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=248, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_nbalance, 'boxN_slow', nbalance_diag%boxN_slow,                                           &
            longname='N-Pool for slowly respirated soil organic material',  units='mol(N) m-2(grid box)',    &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=246, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_nbalance, 'boxN_sminN', nbalance_diag%boxN_SMINN,                                         &
            longname='N-Pool for soil mineral nitrogen',  units='mol(N) m-2(grid box)',                      &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=250, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_nbalance, 'box_Npools_total', nbalance_diag%box_Npools_total,                             &
            longname='Sum of nitrogen from all nitrogen pools', units='mol(N) m-2(grid box)',                &
            ldims=dim1p, gdims=dim1, dimnames=dim1n, code=247, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
       CALL add(IO_nbalance, 'box_minNflux_litter_green', nbalance_diag%box_minNflux_litter_green,           &
            longname='Mineral N flux assoiated with green litter', units='mol(N) m-2(grid box) s-1 ',        &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=70,  lrerun=.TRUE.,  lmiss=.TRUE.,                 &
            contnorest=init_npools.OR.read_npools, missval=missing_value)
       CALL add(IO_nbalance, 'box_minNflux_litter_wood', nbalance_diag%box_minNflux_litter_wood,             &
            longname='Mineral N flux associated with woody litter', units='mol(N) m-2(grid box) s-1',        &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=71,  lrerun=.TRUE.,  lmiss=.TRUE.,                 &
            contnorest=init_npools.OR.read_npools, missval=missing_value)
       CALL add(IO_nbalance, 'box_minNflux_slow', nbalance_diag%box_minNflux_slow,                           &
            longname='Soil mineral N flux from slow pool', units='mol(N) m-2(grid box)  s-1',                &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=72,  lrerun=.TRUE.,  lmiss=.TRUE.,                 &
            contnorest=init_npools.OR.read_npools, missval=missing_value)
       CALL add(IO_nbalance, 'box_Nplant_uptake', nbalance_diag%box_Nplant_demand,                           &
            longname='N uptake by plants', units='mol(N) m-2(grid box) s-1 ',                                &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=73,  lrerun=.TRUE.,  lmiss=.TRUE.,                 &
            contnorest=.TRUE., missval=missing_value)
       CALL add(IO_nbalance, 'box_Nsoil_uptake', nbalance_diag%box_Nsoil_demand,                             &
            longname='N uptake by soil', units='mol(N) m-2(grid box) s-1 ',                                  &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=74,  lrerun=.TRUE.,  lmiss=.TRUE.,                 &
            contnorest=.TRUE., missval=missing_value)
       CALL add(IO_nbalance, 'box_Ntotal_uptake', nbalance_diag%box_Ntotal_demand,                           &
            longname='total N uptake by plants and soil', units='mol(N) m-2(grid box) s-1 ',                 &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=76,  lrerun=.TRUE.,  lmiss=.TRUE.,                 &
            contnorest=.TRUE., missval=missing_value)
       CALL add(IO_nbalance, 'box_Nfix_to_sminN', nbalance_diag%box_nfix_to_sminn,                           &
            longname='N flux to mineral N pool due to N fixation (BNF)', units='mol(N) m-2(grid box) s-1',   &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=77,  lrerun=.TRUE.,  lmiss=.TRUE.,                 &
            contnorest=init_npools.OR.read_npools, missval=missing_value)
       CALL add(IO_nbalance, 'box_Nfert_to_sminN', nbalance_diag%box_nfert_to_sminn,                         &
            longname='N flux to mineral N pool due to fertilizers', units='mol(N) m-2(grid box) s-1',        &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=78,  lrerun=.TRUE.,  lmiss=.TRUE., lpost=.FALSE.,  &
            contnorest=init_npools.OR.read_npools, missval=missing_value)
       CALL add(IO_nbalance, 'box_Ndep_to_sminN', nbalance_diag%box_ndep_to_sminn,                           &
            longname='N flux to mineral N pool due to N deposition', units='mol (N) m-2(grid box) s-1',      &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=79,  lrerun=.TRUE.,  lmiss=.TRUE.,                 &
            contnorest=init_npools.OR.read_npools, missval=missing_value)
       CALL add(IO_nbalance, 'box_sminN_to_denit', nbalance_diag%box_sminn_to_denit,                         &
            longname='N flux from mineral N pool due to denitrification', units='mol(N) m-2(grid box)s-1',   &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=80,  lrerun=.TRUE.,  lmiss=.TRUE.,                 &
            contnorest=init_npools.OR.read_npools, missval=missing_value)
       CALL add(IO_nbalance, 'box_sminN_leach', nbalance_diag%box_sminn_leach,                               &
            longname='N flux to mineral N pool due to leaching', units='mol(N) m-2(grid box)s-1',            &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=81,  lrerun=.TRUE.,  lmiss=.TRUE.,                 &
            contnorest=init_npools.OR.read_npools, missval=missing_value)
       CALL add(IO_nbalance, 'box_sminN_herbivory', nbalance_diag%box_SMINN_herbivory,                       &
            longname='N flux from herbivory via animal faeces to mineral N pool', units='mol(N) m-2(grid box) s-1 ', &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=82,  lrerun=.TRUE.,  lmiss=.TRUE.,                 &
            contnorest=init_npools.OR.read_npools, missval=missing_value)
       CALL add(IO_nbalance, 'box_Nflx_2_crop_harvest', nbalance_diag%box_Nflx_2_crop_harvest,               &
            longname='N flux from vegetation to crop harvest pool', units='mol(N) m-2(grid box) s-1 ',       &
            ldims=dim1p, gdims=dim1, dimnames=dim1n, code=158, lrerun=.TRUE.,  lmiss=.TRUE.,                 &
            contnorest=.TRUE., missval=missing_value)
!! TR       CALL add(IO_nbalance, 'box_sminN_gain', nbalance_diag%box_tot_sminn_gain,                             &
!! TR            longname='Total N influx to soil mineral N pool', units='mol(N) m-2(grid box)s-1',               &
!! TR            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=83,  lrerun=.TRUE.,  lmiss=.TRUE.,                 &
!! TR            contnorest=init_npools.OR.read_npools, missval=missing_value)
       CALL add(IO_nbalance, 'box_sminN_loss', nbalance_diag%box_tot_sminn_loss,                             &
            longname='Total N outflux of soil mineral N pool', units='mol(N) m-2(grid box)s-1',              &
            ldims=dim3p, gdims=dim3, dimnames=dim3n, code=84,  lrerun=.TRUE.,  lmiss=.TRUE.,                 &
            contnorest=init_npools.OR.read_npools, missval=missing_value)
       IF (debug_Nconservation) THEN
          CALL add(IO_nbalance, 'testNconserv', nbalance_diag%testNconserv, contnorest=.TRUE.,                          &
               longname='For testing conservation of nitrogen mass: should be zero.', units='mol(N) m-2(canopy) day-1', &
               ldims=dim3p, gdims=dim3, dimnames=dim3n, code=249, lrerun=.FALSE., lmiss=.TRUE., missval=missing_value)
       END IF
       IF (test_Nconservation) THEN
          CALL add(IO_nbalance, 'jsbachNconserv', nbalance_diag%jsbachNconserv, laccu=.FALSE., contnorest=.TRUE., &
               longname='jsbach nitrogen conservation test: should be zero.',   units='mol(N2) m-2(grid box)', &
               ldims=dim1p, gdims=dim1, dimnames=dim1n, code=138, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
          CALL add(IO_nbalance, 'box_npools_total_old', nbalance_diag%box_Npools_total_old, contnorest=.TRUE., &
               longname='array needed for nitrogen conservation test', units='mol(N2) m-2(grid box)', lpost=.FALSE., &
               ldims=dim1p, gdims=dim1, dimnames=dim1n, code=140, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
          CALL add(IO_nbalance, 'box_npools_total_old2', nbalance_diag%box_Npools_total_old2, contnorest=.TRUE., &
               longname='array needed for nitrogen conservation test', units='mol(N2) m-2(grid box)', lpost=.FALSE., &
               ldims=dim1p, gdims=dim1, dimnames=dim1n, code=141, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
          CALL add(IO_nbalance, 'N2_fluxes_accu', nbalance_diag%N2_fluxes_accu, laccu=.FALSE., lpost=debug_Nconservation, &
               longname='N2 fluxes, accumulated over a day', units='mol(N2) m-2(grid box)', contnorest=.TRUE., &
               ldims=dim1p, gdims=dim1, dimnames=dim1n, code=139, lrerun=.TRUE., lmiss=.TRUE., missval=missing_value)
       END IF
    END IF

    IF (lcc_scheme==2) THEN
       ! Anthropogenic lcc cpools on grid box area averaged over output period
       CALL add(IO_cbalance,'box_Cpool_onSite_avg_LCC', cbalance_diag%boxC_onSite_avg,                              &
                units='mol(C) m-2(grid box)', longname='Carbon left on ground from land use change',code=218,       &
                contnorest=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n,       &
                lpost=.TRUE.,laccu=.TRUE.,lrerun=.FALSE.)
       CALL add(IO_cbalance,'box_Cpool_paper_avg_LCC', cbalance_diag%boxC_paper_avg,units='mol(C) m-2(grid_box)',   &
                code=220,longname='Wood carbon in short/intermediate term anthropogenic pool from land use change', &
                contnorest=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n,       &
                lpost=.TRUE.,laccu=.TRUE.,lrerun=.FALSE.)
       CALL add(IO_cbalance,'box_Cpool_construction_avg_LCC', cbalance_diag%boxC_construction_avg,                  &
                longname='Wood carbon in long term anthropogenic pool from land use change',code=221,               &
                contnorest=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n,       &
                lpost=.TRUE.,laccu=.TRUE.,lrerun=.FALSE.,units='mol(C) m-2(grid box)')
       IF (UseLandcoverTransitions) THEN
          ! Anthropogenic harvest cpools on grid box area averaged over output period
          CALL add(IO_cbalance,'box_Cpool_onSite_harvest_avg',cbalance_diag%boxC_onSite_harvest_avg,units='mol(C)m-2(grid box)', &
                   longname='Carbon left on ground from harvest',code=224,                                                       &
                   contnorest=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n,                 &
                   lpost=.TRUE.,laccu=.TRUE.,lrerun=.FALSE.)
          CALL add(IO_cbalance,'box_Cpool_paper_harvest_avg', cbalance_diag%boxC_paper_harvest_avg,units='mol(C) m-2(grid box)', &
                   longname='Wood carbon in short/intermediate term anthropogenic pool from harvest',code=222,                   &
                   contnorest=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n,                 &
                   lpost=.TRUE.,laccu=.TRUE.,lrerun=.FALSE.)
          CALL add(IO_cbalance,'box_Cpool_construction_harvest_avg', cbalance_diag%boxC_construction_harvest_avg, code=223,      &
                   units='mol(C) m-2(grid box)', longname='Wood carbon in long term anthropogenic pool from harvest',            &
                   contnorest=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n,                 &
                   lpost=.TRUE.,laccu=.TRUE.,lrerun=.FALSE.)
       ENDIF
    ENDIF

    ! --- initializations -------------------------------

    IF (read_cpools) THEN

       CALL message('init_cbalance_bethy','Reading initial carbon pools from file')

       ALLOCATE(init_Cpool_green(l_nland,ntiles))
       ALLOCATE(init_Cpool_woods(l_nland,ntiles))
       ALLOCATE(init_Cpool_reserve(l_nland,ntiles))
       ALLOCATE(init_Cpool_crop_harvest(l_nland,ntiles))
       ALLOCATE(init_frac_litter_wood_new(l_nland,ntiles))
       IF (with_yasso) THEN
          ALLOCATE(init_YCpool_acid_ag1(l_nland,ntiles))
          ALLOCATE(init_YCpool_water_ag1(l_nland,ntiles))
          ALLOCATE(init_YCpool_ethanol_ag1(l_nland,ntiles))
          ALLOCATE(init_YCpool_nonsoluble_ag1(l_nland,ntiles))
          ALLOCATE(init_YCpool_acid_bg1(l_nland,ntiles))
          ALLOCATE(init_YCpool_water_bg1(l_nland,ntiles))
          ALLOCATE(init_YCpool_ethanol_bg1(l_nland,ntiles))
          ALLOCATE(init_YCpool_nonsoluble_bg1(l_nland,ntiles))
          ALLOCATE(init_YCpool_humus_1(l_nland,ntiles))
          ALLOCATE(init_YCpool_acid_ag2(l_nland,ntiles))
          ALLOCATE(init_YCpool_water_ag2(l_nland,ntiles))
          ALLOCATE(init_YCpool_ethanol_ag2(l_nland,ntiles))
          ALLOCATE(init_YCpool_nonsoluble_ag2(l_nland,ntiles))
          ALLOCATE(init_YCpool_acid_bg2(l_nland,ntiles))
          ALLOCATE(init_YCpool_water_bg2(l_nland,ntiles))
          ALLOCATE(init_YCpool_ethanol_bg2(l_nland,ntiles))
          ALLOCATE(init_YCpool_nonsoluble_bg2(l_nland,ntiles))
          ALLOCATE(init_YCpool_humus_2(l_nland,ntiles))
       ELSE
          ALLOCATE(init_Cpool_litter_green_ag(l_nland,ntiles))
          ALLOCATE(init_Cpool_litter_green_bg(l_nland,ntiles))
          ALLOCATE(init_Cpool_litter_wood_ag(l_nland,ntiles))
          ALLOCATE(init_Cpool_litter_wood_bg(l_nland,ntiles))
          ALLOCATE(init_Cpool_slow(l_nland,ntiles))
       END IF
       IF (lcc_scheme==2) THEN
          ALLOCATE(init_Cpool_onSite      (l_nland), &
                   init_Cpool_paper       (l_nland), &
                   init_Cpool_construction(l_nland))
          IF (UseLandcoverTransitions) THEN
            ALLOCATE(init_Cpool_onSite_harvest      (l_nland), &
                     init_Cpool_paper_harvest       (l_nland), &
                     init_Cpool_construction_harvest(l_nland))
          ENDIF
       ENDIF
       
       ALLOCATE(zreal2d(domain%ndim,domain%nblocks))

       !! --- Get green pool (and allocate further memory)
       IF (p_parallel_io) THEN
          cpool_file%opened = .FALSE.
          CALL IO_open(TRIM(cpool_file_name), cpool_file, IO_READ)
          IO_file_id = cpool_file%file_id

          ALLOCATE(zreal3d(grid%nlon,grid%nlat,ntiles))

          CALL IO_inq_varid(IO_file_id, 'Cpool_green', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
       END IF
       NULLIFY(zreal2d_ptr)
       DO i=1,ntiles
          IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
          CALL scatter_field(zreal2d_ptr, zreal2d)
          init_Cpool_green(:,i) = PACK(zreal2d, MASK=domain%mask)
       END DO

       !! --- Get wood pool
       IF (p_parallel_io) THEN
          CALL IO_inq_varid(IO_file_id, 'Cpool_woods', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
       END IF
       NULLIFY(zreal2d_ptr)
       DO i=1,ntiles
          IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
          CALL scatter_field(zreal2d_ptr, zreal2d)
          init_Cpool_woods(:,i) = PACK(zreal2d, MASK=domain%mask)
       END DO

       !! --- Get reserve pool
       IF (p_parallel_io) THEN
          CALL IO_inq_varid(IO_file_id, 'Cpool_reserve', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
       END IF
       NULLIFY(zreal2d_ptr)
       DO i=1,ntiles
          IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
          CALL scatter_field(zreal2d_ptr, zreal2d)
          init_Cpool_reserve(:,i) = PACK(zreal2d, MASK=domain%mask)
       END DO

       !! --- Get crop harvest pool
       IF (p_parallel_io) THEN
          CALL IO_inq_varid(IO_file_id, 'Cpool_crop_harvest', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
       END IF
       NULLIFY(zreal2d_ptr)
       DO i=1,ntiles
          IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
          CALL scatter_field(zreal2d_ptr, zreal2d)
          init_Cpool_crop_harvest(:,i) = PACK(zreal2d, MASK=domain%mask)
       END DO

       !! --- Get fraction of new litter in woody litter pool
       IF (p_parallel_io) THEN
          CALL IO_inq_varid(IO_file_id, 'frac_litter_wood_new', IO_var_id)
          CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
       END IF
       NULLIFY(zreal2d_ptr)
       DO i=1,ntiles
          IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
          CALL scatter_field(zreal2d_ptr, zreal2d)
          init_frac_litter_wood_new(:,i) = PACK(zreal2d, MASK=domain%mask)
       END DO

       IF (with_yasso) THEN

          !! --- Get green litter pools
          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'YCpool_acid_ag1', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)    
             init_YCpool_acid_ag1(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO

          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'YCpool_acid_bg1', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)    
             init_YCpool_acid_bg1(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO

          !! --- Get Yasso water-soluble carbon pools
          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'YCpool_water_ag1', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_YCpool_water_ag1(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO

          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'YCpool_water_bg1', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_YCpool_water_bg1(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO
       
          !! --- Get Yasso ethanol-soluble carbon pools
          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'YCpool_ethanol_ag1', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d(:,:,:))
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_YCpool_ethanol_ag1(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO
       
          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'YCpool_ethanol_bg1', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d(:,:,:))
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_YCpool_ethanol_bg1(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO

          !! --- Get Yasso nonsoluble carbon pools
          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'YCpool_nonsoluble_ag1', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d(:,:,:))
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_YCpool_nonsoluble_ag1(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO

          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'YCpool_nonsoluble_bg1', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d(:,:,:))
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_YCpool_nonsoluble_bg1(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO

          !! --- Get Yasso humus carbon pool
          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'YCpool_humus_1', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d(:,:,:))
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_YCpool_humus_1(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO

           ! SAME FOR SECOND SIZE CLASS:        
          !! --- Get Yasso acid-soluble carbon pools

          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'YCpool_acid_ag2', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)    
             init_YCpool_acid_ag2(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO

          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'YCpool_acid_bg2', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)    
             init_YCpool_acid_bg2(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO

          !! --- Get Yasso water-soluble carbon pools
          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'YCpool_water_ag2', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_YCpool_water_ag2(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO

          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'YCpool_water_bg2', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_YCpool_water_bg2(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO
       
          !! --- Get Yasso ethanol-soluble carbon pools
          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'YCpool_ethanol_ag2', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d(:,:,:))
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_YCpool_ethanol_ag2(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO
       
          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'YCpool_ethanol_bg2', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d(:,:,:))
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_YCpool_ethanol_bg2(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO

          !! --- Get Yasso nonsoluble carbon pools
          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'YCpool_nonsoluble_ag2', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d(:,:,:))
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_YCpool_nonsoluble_ag2(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO

          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'YCpool_nonsoluble_bg2', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d(:,:,:))
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_YCpool_nonsoluble_bg2(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO

          !! --- Get Yasso humus carbon pool
          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'YCpool_humus_2', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d(:,:,:))
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_YCpool_humus_2(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO

       ELSE ! without yasso

          !! --- Get green litter pools
          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'Cpool_litter_green_ag', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_Cpool_litter_green_ag(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO

          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'Cpool_litter_green_bg', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d(:,:,:))
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_Cpool_litter_green_bg(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO

          !! --- Get woody litter pools
          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'Cpool_litter_wood_ag', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_Cpool_litter_wood_ag(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO

          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'Cpool_litter_wood_bg', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_Cpool_litter_wood_bg(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO

          !! --- Get slow carbon pool
          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'Cpool_slow', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d(:,:,:))
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_Cpool_slow(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO
       END IF

       IF (lcc_scheme==2) THEN
          ALLOCATE(zreal2d_ptr(grid%nlon,grid%nlat))

          !! --- Get onSite anthropogenic pool
          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'Cpool_onSite_LCC', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal2d_ptr)
          END IF
          CALL scatter_field(zreal2d_ptr, zreal2d)
          init_Cpool_onSite(:) = PACK(zreal2d, MASK=domain%mask)

          !! --- Get anthropogenic paper pool
          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'Cpool_paper_LCC', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal2d_ptr)
          END IF
          CALL scatter_field(zreal2d_ptr, zreal2d)
          init_Cpool_paper(:) = PACK(zreal2d, MASK=domain%mask)

          !! --- Get anthropogenic construction pool
          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'Cpool_construction_LCC', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal2d_ptr)
          END IF
          CALL scatter_field(zreal2d_ptr, zreal2d)
          init_Cpool_construction(:) = PACK(zreal2d, MASK=domain%mask)

          IF (UseLandcoverTransitions) THEN

             !! --- Get anthropogenic onSite pool from harvest
             IF (p_parallel_io) THEN
                CALL IO_inq_varid(IO_file_id, 'Cpool_onSite_harvest', IO_var_id)
                CALL IO_get_var_double(IO_file_id, IO_var_id, zreal2d_ptr)
             END IF
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_Cpool_onSite_harvest(:) = PACK(zreal2d, MASK=domain%mask)

             !! --- Get anthropogenic paper pool from harvest
             IF (p_parallel_io) THEN
                CALL IO_inq_varid(IO_file_id, 'Cpool_paper_harvest', IO_var_id)
                CALL IO_get_var_double(IO_file_id, IO_var_id, zreal2d_ptr)
             END IF
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_Cpool_paper_harvest(:) = PACK(zreal2d, MASK=domain%mask)

             !! --- Get anthropogenic construction pool from harvest
             IF (p_parallel_io) THEN
                CALL IO_inq_varid(IO_file_id, 'Cpool_construction_harvest', IO_var_id)
                CALL IO_get_var_double(IO_file_id, IO_var_id, zreal2d_ptr)
             END IF
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_Cpool_construction_harvest(:) = PACK(zreal2d, MASK=domain%mask)
          ENDIF
          DEALLOCATE(zreal2d_ptr)
       ENDIF

       !! --- Finish

       IF (p_parallel_io) THEN
          CALL IO_close(cpool_file)
          DEALLOCATE(zreal3d)
       END IF
       DEALLOCATE(zreal2d)

       IF (ANY(ABS(init_Cpool_green) > 10000._dp)) THEN
          CALL message('init_cbalance_bethy','Found cpool values in ini file with absolute value > 10000')
          CALL message('                   ','  This is likely due to using a cpool file that was')
          CALL message('                   ','  generated with a different land-sea mask.')
          CALL message('                   ','  Check your setup!')
          CALL finish('init_cbalance_bethy','')
       END IF

    ELSE !! Default initialization of cbalance carbon pools

       cbalance%Cpool_green(:,:)   = 0.0_dp !![mol[C]/m2]
       cbalance%Cpool_woods(:,:)   = 0.0_dp
       cbalance%Cpool_reserve(:,:) = 0.0_dp
       cbalance%Cpool_crop_harvest(:,:) = 0.0_dp
       cbalance%Cpool_litter_green_ag(:,:) = 60.0_dp
       cbalance%Cpool_litter_green_bg(:,:) = 60.0_dp
       cbalance%Cpool_litter_wood_bg(:,:) = 180.0_dp
       cbalance%Cpool_litter_wood_ag(:,:) = 180.0_dp
       cbalance%Cpool_slow(:,:)        = 2400.0_dp
       IF (lcc_scheme==2) THEN
          cbalance%Cpool_onSite      (:) = 0._dp
          cbalance%Cpool_paper       (:) = 0._dp
          cbalance%Cpool_construction(:) = 0._dp
          IF (UseLandcoverTransitions) THEN
             cbalance%Cpool_onSite_harvest      (:) = 0._dp
             cbalance%Cpool_paper_harvest       (:) = 0._dp
             cbalance%Cpool_construction_harvest(:) = 0._dp
          ENDIF
       ENDIF
       CALL message('init_cbalance_bethy','Carbon pools initialized with default values')

       IF (with_yasso) THEN
          !! if not read from file, initialize yasso carbon pools with zero
          cbalance%YCpool_acid_ag1(:,:)        = 0.0_dp ![mol[C]/m2]  
          cbalance%YCpool_water_ag1(:,:)       = 0.0_dp
          cbalance%YCpool_ethanol_ag1(:,:)     = 0.0_dp
          cbalance%YCpool_nonsoluble_ag1(:,:)  = 0.0_dp
          cbalance%YCpool_acid_bg1(:,:)        = 0.0_dp 
          cbalance%YCpool_water_bg1(:,:)       = 0.0_dp
          cbalance%YCpool_ethanol_bg1(:,:)     = 0.0_dp
          cbalance%YCpool_nonsoluble_bg1(:,:)  = 0.0_dp
          cbalance%YCpool_humus_1(:,:)         = 0.0_dp
          cbalance%YCpool_acid_ag2(:,:)        = 0.0_dp
          cbalance%YCpool_water_ag2(:,:)       = 0.0_dp
          cbalance%YCpool_ethanol_ag2(:,:)     = 0.0_dp
          cbalance%YCpool_nonsoluble_ag2(:,:)  = 0.0_dp
          cbalance%YCpool_acid_bg2(:,:)        = 0.0_dp 
          cbalance%YCpool_water_bg2(:,:)       = 0.0_dp
          cbalance%YCpool_ethanol_bg2(:,:)     = 0.0_dp
          cbalance%YCpool_nonsoluble_bg2(:,:)  = 0.0_dp
          cbalance%YCpool_humus_2(:,:)         = 0.0_dp

          CALL message('init_cbalance_bethy','All Yasso Carbon Pools set to 0 mol(C)/m^2 !!!')
       END IF
  
    END IF !! (read_cpools)
    
    IF (with_forest_management) THEN
        cbalance%FOM_fm_Box_burning_pool = 0.0_dp
        cbalance%FOM_fm_Box_paper_pool = 0.0_dp
        cbalance%FOM_fm_Box_construction_pool = 0.0_dp
        cbalance%FOM_fm_Box_eternity_pool = 0.0_dp
    ENDIF

    cbalance%NPP_Rate = 0.0_dp
    cbalance%soil_respiration = 0.0_dp
    cbalance%soil_respiration_pot = 0.0_dp
    cbalance%Cflux_herbivory_2_atm = 0.0_dp
    cbalance%Cflx_2_crop_harvest = 0.0_dp
    cbalance%Cflx_crop_harvest_2_atm = 0.0_dp
    cbalance%NPP_flux_correction = 0.0_dp
    cbalance%excess_NPP = 0.0_dp
    cbalance%root_exudates = 0.0_dp
    cbalance%Cflux_herbivory = 0.0_dp
    cbalance%Cflux_herbivory_LG = 0.0_dp

    cbalance_diag%box_root_exudates = 0.0_dp
    cbalance_diag%box_Cflux_herbivory = 0.0_dp
    cbalance_diag%box_Cflux_herbivory_2_atm = 0.0_dp
    cbalance_diag%box_Cflx_2_crop_harvest = 0.0_dp
    cbalance_diag%box_Cflx_crop_harvest_2_atm = 0.0_dp
    cbalance_diag%boxC_green = 0.0_dp
    cbalance_diag%boxC_woods = 0.0_dp
    cbalance_diag%boxC_reserve = 0.0_dp
    cbalance_diag%boxC_crop_harvest = 0.0_dp
    cbalance_diag%boxC_litter_green_ag = 0.0_dp
    cbalance_diag%boxC_litter_green_bg = 0.0_dp
    cbalance_diag%boxC_litter_wood = 0.0_dp
    cbalance_diag%boxC_slow = 0.0_dp
    cbalance_diag%LAI_yDayMean = 0.0_dp
    cbalance_diag%LAI_previousDayMean = 0.0_dp
    cbalance_diag%NPP_yDayMean = 0.0_dp
    cbalance_diag%GPP_yDayMean = 0.0_dp
    cbalance_diag%topSoilTemp_yDayMean = 0.0_dp
    cbalance_diag%alpha_yDayMean = 0.0_dp
    cbalance_diag%drainage_yDayMean = 0.0_dp
    cbalance_diag%box_soil_respiration = 0.0_dp
    cbalance_diag%box_soil_respiration_pot = 0.0_dp
    cbalance_diag%box_NPP_yDayMean = 0.0_dp
    cbalance_diag%box_NPP_act_yDayMean = 0.0_dp
    cbalance_diag%box_NPP_flux_correction = 0.0_dp
    cbalance_diag%box_GPP_yDayMean = 0.0_dp
    cbalance_diag%litter_flux = 0.0_dp
    cbalance_diag%box_litter_flux = 0.0_dp
    ! ---HW
    cbalance_diag%leaf_shedding_debug = 0.0_dp
    cbalance_diag%Cpool_green_minus_shed_leaves = 0.0_dp
    cbalance_diag%excess_carbon_debug = 0.0_dp
    cbalance_diag%litter_leaf = 0.0_dp
    cbalance_diag%box_litter_leaf = 0.0_dp
    cbalance_diag%relative_extractable_water = 0.0_dp
    ! HW---
    cbalance_diag%box_Cpools_total = 0.0_dp
    IF (lcc_scheme==2) THEN
       cbalance_diag%boxC_onSite_avg      (:) = 0._dp
       cbalance_diag%boxC_paper_avg       (:) = 0._dp
       cbalance_diag%boxC_construction_avg(:) = 0._dp
       IF (UseLandcoverTransitions) THEN
          cbalance_diag%boxC_onsite_harvest_avg      (:) = 0._dp
          cbalance_diag%boxC_paper_harvest_avg       (:) = 0._dp
          cbalance_diag%boxC_construction_harvest_avg(:) = 0._dp
       ENDIF
    ENDIF
    IF (debug_Cconservation) THEN
       cbalance_diag%testCconserv = 0.0_dp
    END IF
    IF (test_Cconservation) THEN
       cbalance_diag%jsbachCconserv = 0.0_dp 
       cbalance_diag%box_Cpools_total_old = 0.0_dp
       cbalance_diag%box_Cpools_total_old2 = 0.0_dp
       cbalance_diag%CO2_fluxes_accu = 0.0_dp
    END IF

    ! Computation of pseudo-15-day mean values of air temperature and precipitation           
    IF (with_yasso) THEN

       F_pseudo_temp = EXP(-delta_time / (15._dp * 86400._dp))
       N_pseudo_temp = 1._dp / (1._dp - F_pseudo_temp)

       F_pseudo_precip = EXP(-delta_time / (15._dp * 86400._dp))
       N_pseudo_precip = 1._dp / (1._dp - F_pseudo_precip)

       cbalance%pseudo_temp   = 285.0_dp     !Arbitrary initial values
       cbalance_diag%pseudo_temp_yDay = 285.0_dp !Arbitrary initial values
       cbalance%pseudo_precip = 3.0E-5_dp 
       cbalance_diag%pseudo_precip_yDay = 3.0E-5_dp

    END IF


    !! -- initialization of Nitrogen dynamics ----------------------------------------
    
    cbalance_diag%NPP_yDayMean(:,:) = 0.0_dp
    cbalance%NPP_act_yDayMean  = 1.0E-13_dp        ! jsbach interface [photosyn (called by update_bethy)] before they are computed
                                                   ! in [update_cbalance_bethy]. Initializing to zero would result in a 
                                                   ! zero divide in the first time step of an initialized run.
    IF (with_nitrogen) THEN
       
       nbalance%SMINN_pool              = 0.5_dp    !! set to 0.5 at begining 

       nbalance%redFact_Nlimitation     = 1.0_dp        
       nbalance%minNflux_litter_green_ag= 0.0_dp 
       nbalance%minNflux_litter_green_bg= 0.0_dp 
       nbalance%minNflux_litter_wood_ag = 0.0_dp 
       nbalance%minNflux_litter_wood_bg = 0.0_dp 
       nbalance%minNflux_slow           = 0.0_dp 
       nbalance%Nplant_demand           = 0.0_dp 
       nbalance%Nsoil_demand            = 0.0_dp 
       nbalance%Ntotal_demand           = 0.0_dp 
       nbalance%SMINN_herbivory         = 0.0_dp 

       nbalance%nfix_to_sminn        = 0.0_dp 
       nbalance%ndep_to_sminn        = 0.0_dp 
       nbalance%ndep_forc            = 0.0_dp 
       nbalance%nfert_forc           = 0.0_dp
       nbalance%nfert_forc_2d        = 0.0_dp
       nbalance%nfert_to_sminn       = 0.0_dp
       nbalance%sminn_to_denit       = 0.0_dp   
       nbalance%sminn_leach          = 0.0_dp
       nbalance%Nflx_2_crop_harvest  = 0.0_dp
       nbalance%N2O_emissions_depfix = 0.0_dp
       nbalance%N2O_emissions_nfert  = 0.0_dp     
       nbalance%N2O_emissions_mineraliz= 0.0_dp
       nbalance%N2O_emissions_slow   = 0.0_dp
       nbalance%N2O_emissions_grazing = 0.0_dp
       nbalance%NetEcosyst_N_flux(:,:) = 0.0_dp

       IF (.NOT. lresume) THEN
          nbalance%NPP_run_mean(:,:) = 0.0_dp
       END IF
       IF (debug_Nconservation) THEN
          nbalance_diag%testNconserv = 0.0_dp 
       END IF
       IF (test_Nconservation) THEN
          nbalance_diag%jsbachNconserv = 0.0_dp 
          nbalance_diag%box_Npools_total_old = 0.0_dp
          nbalance_diag%box_Npools_total_old2 = 0.0_dp
          nbalance_diag%N2_fluxes_accu = 0.0_dp
       END IF

       
       IF (read_npools) THEN  ! all Nitrogen pools are read from file

          CALL message('init_cbalance_bethy','Reading initial nitrogen pools from file: only pools without const C/N')

          ALLOCATE(init_Npool_green(l_nland,ntiles))
          ALLOCATE(init_Npool_woods(l_nland,ntiles))
          ALLOCATE(init_Npool_mobile(l_nland,ntiles))
          ALLOCATE(init_Npool_crop_harvest(l_nland,ntiles))
          ALLOCATE(init_Npool_litter_green_ag(l_nland,ntiles))
          ALLOCATE(init_Npool_litter_green_bg(l_nland,ntiles))
          ALLOCATE(init_Npool_litter_wood_ag(l_nland,ntiles))
          ALLOCATE(init_Npool_litter_wood_bg(l_nland,ntiles))
          ALLOCATE(init_Npool_slow(l_nland,ntiles))
          ALLOCATE(init_SMINN_pool(l_nland,ntiles))
          ALLOCATE(init_NPP_run_mean(l_nland,ntiles))

          ALLOCATE(zreal2d(domain%ndim,domain%nblocks))

          !! --- Get green N pool

          IF (p_parallel_io) THEN
             npool_file%opened = .FALSE.
             CALL IO_open(TRIM(npool_file_name), npool_file, IO_READ)
             IO_file_id = npool_file%file_id

             ALLOCATE(zreal3d(grid%nlon,grid%nlat,ntiles))

             CALL IO_inq_varid(IO_file_id, 'Npool_green', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_Npool_green(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO

          !! --- Get woods N pool

          IF (p_parallel_io) THEN
             npool_file%opened = .FALSE.
             CALL IO_open(TRIM(npool_file_name), npool_file, IO_READ)
             IO_file_id = npool_file%file_id

             CALL IO_inq_varid(IO_file_id, 'Npool_woods', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_Npool_woods(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO

          !! --- Get mobile N (reserve) pool

          IF (p_parallel_io) THEN
             npool_file%opened = .FALSE.
             CALL IO_open(TRIM(npool_file_name), npool_file, IO_READ)
             IO_file_id = npool_file%file_id

             CALL IO_inq_varid(IO_file_id, 'Npool_mobile', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_Npool_mobile(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO

          !! --- Get crop harvest
          IF (p_parallel_io) THEN
             npool_file%opened = .FALSE.
             CALL IO_open(TRIM(npool_file_name), npool_file, IO_READ)
             IO_file_id = npool_file%file_id

             CALL IO_inq_varid(IO_file_id, 'Npool_crop_harvest', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_Npool_crop_harvest(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO

          !! --- Get green litter pools

          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'Npool_litter_green_ag', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_Npool_litter_green_ag(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO

          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'Npool_litter_green_bg', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d(:,:,:))
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_Npool_litter_green_bg(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO

          !! --- Get woody litter pools

          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'Npool_litter_wood_ag', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_Npool_litter_wood_ag(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO

          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'Npool_litter_wood_bg', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d)
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_Npool_litter_wood_bg(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO

          !! --- Get slow pool

          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'Npool_slow', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d(:,:,:))
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_Npool_slow(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO

          !! --- Get SMINN pool  check if it needed 

          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'sminN_pool', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d(:,:,:))
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_SMINN_pool(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO

          !! --- Get NPP running mean 

          IF (p_parallel_io) THEN
             CALL IO_inq_varid(IO_file_id, 'NPP_run_mean', IO_var_id)
             CALL IO_get_var_double(IO_file_id, IO_var_id, zreal3d(:,:,:))
          END IF
          NULLIFY(zreal2d_ptr)
          DO i=1,ntiles
             IF (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
             CALL scatter_field(zreal2d_ptr, zreal2d)
             init_NPP_run_mean(:,i) = PACK(zreal2d, MASK=domain%mask)
          END DO

          !! -- close file and deallocate temporary memory
          IF (p_parallel_io) THEN
             CALL IO_close(npool_file)
             DEALLOCATE(zreal3d)
          END IF
          DEALLOCATE(zreal2d)

       END IF  ! read_npools
    END IF  ! with nitrogen

  END SUBROUTINE init_cbalance_bethy

  ! --- update_cbalance_bethy()  --------------------------------------------------------------------------------------------------

  SUBROUTINE update_cbalance_bethy(cbalance, &
                                   nbalance, &
                                   with_nitrogen,&
                                   with_yasso, &
                                   with_landcover_maps, & 
                                   with_landcover_transitions, & 
                                   with_forest_management, &
                                   lctlib, &
                                   surface, &
                                   lcc_scheme, &
                                   grossAssimilation, &
                                   darkRespiration, &
                                   topSoilTemp, &
                                   alpha, &
                                   LAI, &
                                   LAI_max, &
                                   veg_fract_correction, &
                                   temp_air, precip, &
                                   drainage,  & 
                                   maxmoisture, &
                                   CO2_flx2atm_npp, &
                                   CO2_flx2atm_soilresp, &
                                   CO2_flx2atm_herbivory, &
                                   CO2_flx2atm_npp_test, &
                                   N2O_flx2atm_mineraliz, &
                                   N2O_flx2atm_depfix, &
                                   N2O_flx2atm_nfert, &
                                   N2O_flx2atm_grazing, &
                                   N2_flx2atm_ecosystem &
                                   , leaf_shedding_rate & ! HW
                                   , relative_extractable_water & ! HW
                                   )

    
    USE mo_jsbach,                 ONLY: new_day, new_year
    USE mo_time_control,           ONLY: lstart,lresume
    USE mo_land_surface,           ONLY: land_surface_type
    USE mo_jsbach_lctlib,          ONLY: lctlib_type
    USE mo_cbal_cpools,            ONLY: update_Cpools, N_process
    USE mo_cbal_parameters,        ONLY: cn_green, cn_woods, cn_litter_green, cn_litter_wood, cn_slow
#if defined (__SX__) && defined (_OPENMP)
    USE omp_lib,                   ONLY: omp_get_thread_num
#endif

    TYPE(cbalance_type),    INTENT(inout) :: cbalance
    TYPE(nbalance_type),    INTENT(inout) :: nbalance
    LOGICAL,                INTENT(in)    :: with_nitrogen !! If .true. Nitrogen will be cycled in addition to carbon
    LOGICAL,                INTENT(in)    :: with_yasso
    LOGICAL,                INTENT(in)    :: with_landcover_maps, with_landcover_transitions
    LOGICAL,                INTENT(in)    :: with_forest_management
    TYPE(lctlib_type),      INTENT(in)    :: lctlib
    TYPE(land_surface_type),INTENT(in)    :: surface
    INTEGER,                INTENT(in)    :: lcc_scheme

    REAL(dp),INTENT(in)               :: grossAssimilation(:,:)    !! Gross primary product (where photorespiration has already 
                                                                   !!    been accounted for) [mol(CO2)/(m^2 s)]
    REAL(dp),INTENT(in)               :: darkRespiration(:,:)      !! dark respiration of leaves [mol(CO2)/(m^2 s)]
    REAL(dp),INTENT(in)               :: topSoilTemp(:)            !! As temperature of the soil pool the temperature in the top 
                                                                   !!    soil layer (on all tiles identical) is taken [Kelvin]
    REAL(dp),INTENT(in)               :: alpha(:,:)                !! relative soil moisture (used for heterotrophic respiration)
    REAL(dp),INTENT(in)               :: LAI(:,:)                  !! leaf area index of (current time step)
    REAL(dp),INTENT(in)               :: LAI_max(:,:)              !! maximum leaf area index of (current time step)
    REAL(dp),INTENT(in)               :: veg_fract_correction(:,:) !! Correction factor for cover fractions 1-exp(-LAI_max/2)
    REAL(dp),INTENT(in)               :: temp_air(:)               !! Air temperature (lowest level of the atmospere) [Kelvin]
    REAL(dp),INTENT(in)               :: precip(:)                 !! Precipitation rate [kg/(m^2 s)]
    real(dp),intent(in)               :: drainage(:)               !! drainage in [m/s]
    real(dp),intent(in)               :: maxmoisture(:)            !! depth of soil water buckets [m]

    real(dp),intent(in)               :: leaf_shedding_rate(:,:)   !! HW

    real(dp),intent(in)               :: relative_extractable_water(:,:)   !! HW REW

    REAL(dp),INTENT(out)              :: CO2_flx2atm_npp(:)        !! grid cell averages of net CO2 fluxes between
    REAL(dp),INTENT(out)              :: CO2_flx2atm_soilresp(:)   !! .. biosphere (due to NPP, soil respiration and
    REAL(dp),INTENT(out)              :: CO2_flx2atm_herbivory(:)  !! .. grazing) and atmosphere [kg(CO2)/(m^2(ground) s)]
                                                                   !! .. (positive for emission to atmosphere,
                                                                   !! ..  negative for absorption by biosphere)
    REAL(dp),INTENT(out)              :: CO2_flx2atm_npp_test(:)   !! CO2 flux due to NPP, only needed for conservation test. 
    REAL(dp),INTENT(out)              :: N2O_flx2atm_mineraliz(:)
    REAL(dp),INTENT(out)              :: N2O_flx2atm_depfix(:)
    REAL(dp),INTENT(out)              :: N2O_flx2atm_nfert(:)
    REAL(dp),INTENT(out)              :: N2O_flx2atm_grazing(:)
    REAL(dp),INTENT(out)              :: N2_flx2atm_ecosystem(:)

    ! local variables

    integer :: i, itile,k ! HW add k index
    integer :: kidx0,kidx1

    real(dp) :: frac_npp_2_woodPool(1:nidx,1:ntiles)
    real(dp) :: frac_npp_2_reservePool(1:nidx,1:ntiles)
    real(dp) :: frac_npp_2_exudates(1:nidx,1:ntiles)
    real(dp) :: frac_green_2_herbivory(1:nidx,1:ntiles)
    real(dp) :: tau_Cpool_litter_green(1:nidx,1:ntiles)
    real(dp) :: tau_Cpool_litter_wood(1:nidx,1:ntiles)
    real(dp) :: tau_Cpool_woods(1:nidx,1:ntiles)
    real(dp) :: LAI_shed_constant(1:nidx,1:ntiles)
    real(dp) :: frac_C_litter_green2atmos(1:nidx,1:ntiles)
    real(dp) :: Max_C_content_woods(1:nidx,1:ntiles)
    real(dp) :: specificLeafArea_C(nidx,1:ntiles)
    real(dp) :: reserveC2leafC(nidx,1:ntiles)
    real(dp) :: areaWeightingFactor(1:nidx,1:ntiles)
    real(dp) :: ndep_forc_tiles(1:nidx,1:ntiles)
    real(dp) ::  LeafLit_coef(1:nidx,1:ntiles,5)     ! Yasso vectors
    real(dp) ::  WoodLit_coef(1:nidx,1:ntiles,5)     !
    real(dp) ::  WoodLitterSize(1:nidx,1:ntiles)     !


#if defined (__SX__) && defined (_OPENMP)
    INTEGER :: tid
#endif

#if defined (__SX__) && defined (_OPENMP)
    tid = omp_get_thread_num()
#endif

    ! --- set block range indices

    kidx0 = kstart    ! first value of kindex() (the index of the first element of a block in the processor domain)
    kidx1 = kend      ! last value of kindex()  (the index of the last element of a block in the processor domain)

    IF (debug .AND. new_day) CALL message('update_cbalance_bethy','First time step of new day')
    IF (debug .AND. new_year) CALL message('update_cbalance_bethy','First time step of new year')

    ! initializations

    CO2_flx2atm_npp(1:nidx)       = 0._dp
    CO2_flx2atm_soilresp(1:nidx)  = 0._dp
    CO2_flx2atm_herbivory(1:nidx) = 0._dp
    N2O_flx2atm_mineraliz(1:nidx) = 0._dp
    N2O_flx2atm_depfix(1:nidx)    = 0._dp
    N2O_flx2atm_nfert(1:nidx)     = 0._dp
    N2O_flx2atm_grazing(1:nidx)   = 0._dp
    N2_flx2atm_ecosystem(1:nidx)  = 0._dp

    cbalance_diag%leaf_shedding_rate=leaf_shedding_rate ! HW
    cbalance_diag%relative_extractable_water=relative_extractable_water ! HW

    IF (read_cpools .AND. (lstart .OR. lresume)) THEN
       cbalance%Cpool_green(kidx0:kidx1,:)   = init_Cpool_green(kidx0:kidx1,:)
       cbalance%Cpool_woods(kidx0:kidx1,:)   = init_Cpool_woods(kidx0:kidx1,:)
       cbalance%Cpool_reserve(kidx0:kidx1,:) = init_Cpool_reserve(kidx0:kidx1,:)
       cbalance%Cpool_crop_harvest(kidx0:kidx1,:) = init_Cpool_crop_harvest(kidx0:kidx1,:)
       cbalance%frac_litter_wood_new(kidx0:kidx1,:) = init_frac_litter_wood_new(kidx0:kidx1,:)
       IF ((with_landcover_maps .OR. with_landcover_transitions) .AND. lcc_scheme==2) THEN
          cbalance%Cpool_onSite      (:) = 0._dp
          cbalance%Cpool_paper       (:) = 0._dp
          cbalance%Cpool_construction(:) = 0._dp
          IF (with_landcover_transitions) THEN
             cbalance%Cpool_onSite_harvest      (:) = 0._dp
             cbalance%Cpool_paper_harvest       (:) = 0._dp
             cbalance%Cpool_construction_harvest(:) = 0._dp
          ENDIF
       ENDIF
       IF (with_yasso) THEN
          cbalance%YCpool_acid_ag1(kidx0:kidx1,:)       = init_YCpool_acid_ag1(kidx0:kidx1,:)
          cbalance%YCpool_water_ag1(kidx0:kidx1,:)      = init_YCpool_water_ag1(kidx0:kidx1,:)
          cbalance%YCpool_ethanol_ag1(kidx0:kidx1,:)    = init_YCpool_ethanol_ag1(kidx0:kidx1,:)
          cbalance%YCpool_nonsoluble_ag1(kidx0:kidx1,:) = init_YCpool_nonsoluble_ag1(kidx0:kidx1,:)
          cbalance%YCpool_acid_bg1(kidx0:kidx1,:)       = init_YCpool_acid_bg1(kidx0:kidx1,:)
          cbalance%YCpool_water_bg1(kidx0:kidx1,:)      = init_YCpool_water_bg1(kidx0:kidx1,:)
          cbalance%YCpool_ethanol_bg1(kidx0:kidx1,:)    = init_YCpool_ethanol_bg1(kidx0:kidx1,:)
          cbalance%YCpool_nonsoluble_bg1(kidx0:kidx1,:) = init_YCpool_nonsoluble_bg1(kidx0:kidx1,:)
          cbalance%YCpool_humus_1(kidx0:kidx1,:)        = init_YCpool_humus_1(kidx0:kidx1,:)
          cbalance%YCpool_acid_ag2(kidx0:kidx1,:)       = init_YCpool_acid_ag2(kidx0:kidx1,:)
          cbalance%YCpool_water_ag2(kidx0:kidx1,:)      = init_YCpool_water_ag2(kidx0:kidx1,:)
          cbalance%YCpool_ethanol_ag2(kidx0:kidx1,:)    = init_YCpool_ethanol_ag2(kidx0:kidx1,:)
          cbalance%YCpool_nonsoluble_ag2(kidx0:kidx1,:) = init_YCpool_nonsoluble_ag2(kidx0:kidx1,:)
          cbalance%YCpool_acid_bg2(kidx0:kidx1,:)       = init_YCpool_acid_bg2(kidx0:kidx1,:)
          cbalance%YCpool_water_bg2(kidx0:kidx1,:)      = init_YCpool_water_bg2(kidx0:kidx1,:)
          cbalance%YCpool_ethanol_bg2(kidx0:kidx1,:)    = init_YCpool_ethanol_bg2(kidx0:kidx1,:)
          cbalance%YCpool_nonsoluble_bg2(kidx0:kidx1,:) = init_YCpool_nonsoluble_bg2(kidx0:kidx1,:)
          cbalance%YCpool_humus_2(kidx0:kidx1,:)        = init_YCpool_humus_2(kidx0:kidx1,:)
       ELSE
          cbalance%Cpool_litter_green_ag(kidx0:kidx1,:) = init_Cpool_litter_green_ag(kidx0:kidx1,:)
          cbalance%Cpool_litter_green_bg(kidx0:kidx1,:) = init_Cpool_litter_green_bg(kidx0:kidx1,:)
          cbalance%Cpool_litter_wood_ag(kidx0:kidx1,:)  = init_Cpool_litter_wood_ag(kidx0:kidx1,:)
          cbalance%Cpool_litter_wood_bg(kidx0:kidx1,:)  = init_Cpool_litter_wood_bg(kidx0:kidx1,:)
          cbalance%Cpool_slow(kidx0:kidx1,:)            = init_Cpool_slow(kidx0:kidx1,:)
       END IF
       CALL message ('update_cbalance_bethy','Cpools overwritten with initial values')
    END IF

    IF (init_npools .AND. (lstart .OR. lresume)) THEN

       ! initialize N pools from current C pools

       nbalance%Npool_green(kidx0:kidx1,:)           = cbalance%Cpool_green(kidx0:kidx1,:) / cn_green
       nbalance%Npool_woods(kidx0:kidx1,:)           = cbalance%Cpool_woods(kidx0:kidx1,:) / cn_woods
       nbalance%Npool_crop_harvest(kidx0:kidx1,:)    = cbalance%Cpool_crop_harvest(kidx0:kidx1,:) / cn_green
       IF (with_yasso) THEN
          nbalance%Npool_litter_green_ag(kidx0:kidx1,:) = (  cbalance%YCpool_acid_ag1(kidx0:kidx1,:)        &
                                                           + cbalance%YCpool_water_ag1(kidx0:kidx1,:)       & 
                                                           + cbalance%YCpool_ethanol_ag1(kidx0:kidx1,:)     &
                                                           + cbalance%YCpool_nonsoluble_ag1(kidx0:kidx1,:)) &
                                                         / cn_litter_green
          nbalance%Npool_litter_green_bg(kidx0:kidx1,:) = (  cbalance%YCpool_acid_bg1(kidx0:kidx1,:)        &
                                                           + cbalance%YCpool_water_bg1(kidx0:kidx1,:)       &
                                                           + cbalance%YCpool_ethanol_bg1(kidx0:kidx1,:)     &
                                                           + cbalance%YCpool_nonsoluble_bg1(kidx0:kidx1,:)) &
                                                         / cn_litter_green
          nbalance%Npool_litter_wood_ag(kidx0:kidx1,:)  = (  cbalance%YCpool_acid_ag2(kidx0:kidx1,:)        &
                                                           + cbalance%YCpool_water_ag2(kidx0:kidx1,:)       &
                                                           + cbalance%YCpool_ethanol_ag2(kidx0:kidx1,:)     &
                                                           + cbalance%YCpool_nonsoluble_ag2(kidx0:kidx1,:)) &
                                                         / cn_litter_wood
          nbalance%Npool_litter_wood_bg(kidx0:kidx1,:)  = (  cbalance%YCpool_acid_bg2(kidx0:kidx1,:)        &
                                                           + cbalance%YCpool_water_bg2(kidx0:kidx1,:)       &
                                                           + cbalance%YCpool_ethanol_bg2(kidx0:kidx1,:)     &
                                                           + cbalance%YCpool_nonsoluble_bg2(kidx0:kidx1,:)) &
                                                         / cn_litter_wood
          nbalance%Npool_slow(kidx0:kidx1,:)            = (  cbalance%YCpool_humus_1(kidx0:kidx1,:)         &
                                                           + cbalance%YCpool_humus_2(kidx0:kidx1,:)) / cn_slow
       ELSE
          nbalance%Npool_litter_green_ag(kidx0:kidx1,:) = cbalance%Cpool_litter_green_ag(kidx0:kidx1,:) / cn_litter_green
          nbalance%Npool_litter_green_bg(kidx0:kidx1,:) = cbalance%Cpool_litter_green_bg(kidx0:kidx1,:) / cn_litter_green
          nbalance%Npool_litter_wood_ag(kidx0:kidx1,:)  = cbalance%Cpool_litter_wood_ag(kidx0:kidx1,:)  / cn_litter_wood
          nbalance%Npool_litter_wood_bg(kidx0:kidx1,:)  = cbalance%Cpool_litter_wood_bg(kidx0:kidx1,:)  / cn_litter_wood
          nbalance%Npool_slow(kidx0:kidx1,:)            = cbalance%Cpool_slow(kidx0:kidx1,:)            / cn_slow
       END IF 
       nbalance%Npool_mobile(kidx0:kidx1,:) = nbalance%Npool_green(kidx0:kidx1,:)             !! use N-pool green for initialization
       nbalance%SMINN_pool(kidx0:kidx1,:)   = nbalance%Npool_litter_green_bg(kidx0:kidx1,:) + 1._dp

       CALL message ('update_cbalance_bethy','Npools re-calculated via C/N-ratios from Cpools')

    END IF

    IF (read_npools .AND. (lstart .OR. lresume)) THEN

       ! initialize N pools with values from file
       nbalance%Npool_green(kidx0:kidx1,:)           = init_Npool_green(kidx0:kidx1,:)
       nbalance%Npool_woods(kidx0:kidx1,:)           = init_Npool_woods(kidx0:kidx1,:)
       nbalance%Npool_mobile(kidx0:kidx1,:)          = init_Npool_mobile(kidx0:kidx1,:)
       nbalance%Npool_crop_harvest(kidx0:kidx1,:)    = init_Npool_crop_harvest(kidx0:kidx1,:)
       nbalance%Npool_litter_green_ag(kidx0:kidx1,:) = init_Npool_litter_green_ag(kidx0:kidx1,:)
       nbalance%Npool_litter_green_bg(kidx0:kidx1,:) = init_Npool_litter_green_bg(kidx0:kidx1,:)
       nbalance%Npool_litter_wood_ag(kidx0:kidx1,:)  = init_Npool_litter_wood_ag(kidx0:kidx1,:)
       nbalance%Npool_litter_wood_bg(kidx0:kidx1,:)  = init_Npool_litter_wood_bg(kidx0:kidx1,:)
       nbalance%Npool_slow(kidx0:kidx1,:)            = init_Npool_slow(kidx0:kidx1,:)
       nbalance%SMINN_pool(kidx0:kidx1,:)            = init_SMINN_pool(kidx0:kidx1,:)
       nbalance%NPP_run_mean(kidx0:kidx1,:)          = init_NPP_run_mean(kidx0:kidx1,:)

       CALL message ('update_cbalance_bethy','Npools overwritten with initial values')

    END IF

    ! --- go ...

    ! compute net primary production rate
    cbalance%NPP_rate(kidx0:kidx1,:) = NPP_rate_bethy(grossAssimilation(1:nidx,:),darkRespiration(1:nidx,:))

    ! Prepare area weighting factor to rescale from 1/[m^2(canopy)] to 1/[m^2(grid box)]
    areaWeightingFactor(:,:) = veg_fract_correction(:,:) * surface%cover_fract(kidx0:kidx1,:) &
                               * SPREAD(surface%veg_ratio_max(kidx0:kidx1),DIM=2,NCOPIES=ntiles)

    ! Update pseudo-15day-mean air temperature and precipitation 
    IF (with_yasso) THEN
       cbalance%pseudo_temp(kidx0:kidx1) = &
            temp_air(1:nidx)/N_pseudo_temp + F_pseudo_temp * cbalance%pseudo_temp(kidx0:kidx1)
       cbalance%pseudo_precip(kidx0:kidx1) = &
            precip(1:nidx)/N_pseudo_precip + F_pseudo_precip * cbalance%pseudo_precip(kidx0:kidx1)
    END IF

    IF( .NOT. new_day .OR. lstart) THEN ! perform daily sums ---------------------------------------------------

       cbalance%LAI_sum(kidx0:kidx1,:)          = cbalance%LAI_sum(kidx0:kidx1,:)         + lai(1:nidx,:)
       cbalance%NPP_sum(kidx0:kidx1,:)          = cbalance%NPP_sum(kidx0:kidx1,:)         + cbalance%NPP_rate(kidx0:kidx1,:)
       cbalance%GPP_sum(kidx0:kidx1,:)          = cbalance%GPP_sum(kidx0:kidx1,:)         + grossAssimilation(1:nidx,:)
       cbalance%topSoilTemp_sum(kidx0:kidx1)    = cbalance%topSoilTemp_sum(kidx0:kidx1)   + topSoilTemp(1:nidx)
       cbalance%alpha_sum(kidx0:kidx1,:)        = cbalance%alpha_sum(kidx0:kidx1,:)       + alpha(1:nidx,:)
       cbalance%drainage_sum(kidx0:kidx1)       = cbalance%drainage_sum(kidx0:kidx1)    + drainage(1:nidx)

    ELSE ! A new day begins ==> perform carbon balance -------------------------------------------------------------------

       ! Compute previous days means

       cbalance_diag%LAI_previousDayMean(kidx0:kidx1,:)  = cbalance_diag%LAI_yDayMean(kidx0:kidx1,:)
       cbalance_diag%LAI_yDayMean(kidx0:kidx1,:)         = cbalance%LAI_sum(kidx0:kidx1,:)/time_steps_per_day
       cbalance_diag%NPP_yDayMean(kidx0:kidx1,:)         = cbalance%NPP_sum(kidx0:kidx1,:)/time_steps_per_day
       cbalance_diag%GPP_yDayMean(kidx0:kidx1,:)         = cbalance%GPP_sum(kidx0:kidx1,:)/time_steps_per_day
       cbalance_diag%topSoilTemp_yDayMean(kidx0:kidx1)   = cbalance%topSoilTemp_sum(kidx0:kidx1)/time_steps_per_day
       cbalance_diag%alpha_yDayMean(kidx0:kidx1,:)       = cbalance%alpha_sum(kidx0:kidx1,:)/time_steps_per_day
       cbalance_diag%drainage_yDayMean(kidx0:kidx1)      = cbalance%drainage_sum(kidx0:kidx1)/time_steps_per_day

       ! Restart summing of previous days values

      cbalance%LAI_sum(kidx0:kidx1,:)             = lai(1:nidx,:)
      cbalance%NPP_sum(kidx0:kidx1,:)             = cbalance%NPP_Rate(kidx0:kidx1,:)
      cbalance%GPP_sum(kidx0:kidx1,:)             = grossAssimilation(1:nidx,:)
      cbalance%topSoilTemp_sum(kidx0:kidx1)       = topSoilTemp(1:nidx)
      cbalance%alpha_sum(kidx0:kidx1,:)           = alpha(1:nidx,:)
      cbalance%drainage_sum(kidx0:kidx1)          = drainage(1:nidx)

      IF (with_yasso) THEN
         cbalance_diag%pseudo_temp_yDay(kidx0:kidx1)   = cbalance%pseudo_temp(kidx0:kidx1)
         cbalance_diag%pseudo_precip_yDay(kidx0:kidx1) = cbalance%pseudo_precip(kidx0:kidx1)
      END IF

      ! Construct arrays with land cover type (PFT) properties (SLA, max. wood mass, etc.)

      DO itile=1,ntiles
         frac_npp_2_woodPool(1:nidx,itile)       = lctlib%frac_npp_2_woodPool(surface%cover_type(kidx0:kidx1,itile))
         frac_npp_2_reservePool(1:nidx,itile)    = lctlib%frac_npp_2_reservePool(surface%cover_type(kidx0:kidx1,itile))
         frac_npp_2_exudates(1:nidx,itile)       = lctlib%frac_npp_2_exudates(surface%cover_type(kidx0:kidx1,itile))
         frac_green_2_herbivory(1:nidx,itile)    = lctlib%frac_green_2_herbivory(surface%cover_type(kidx0:kidx1,itile))
         tau_Cpool_litter_green(1:nidx,itile)    = lctlib%tau_Cpool_litter_leaf(surface%cover_type(kidx0:kidx1,itile))
         tau_Cpool_litter_wood(1:nidx,itile)     = lctlib%tau_Cpool_litter_wood(surface%cover_type(kidx0:kidx1,itile))
         tau_Cpool_woods(1:nidx,itile)           = lctlib%tau_Cpool_woods(surface%cover_type(kidx0:kidx1,itile))
         LAI_shed_constant(1:nidx,itile)         = lctlib%LAI_shed_constant(surface%cover_type(kidx0:kidx1,itile))
         frac_C_litter_green2atmos(1:nidx,itile) = lctlib%frac_C_litter_green2atmos(surface%cover_type(kidx0:kidx1,itile))
         Max_C_content_woods(1:nidx,itile)       = lctlib%Max_C_content_woods(surface%cover_type(kidx0:kidx1,itile))
         specificLeafArea_C(1:nidx,itile)        = lctlib%specificLeafArea_C(surface%cover_type(kidx0:kidx1,itile))
         reserveC2leafC(1:nidx,itile)            = lctlib%reserveC2leafC(surface%cover_type(kidx0:kidx1,itile))
         WoodLitterSize(1:nidx,itile)            = lctlib%WoodLitterSize(surface%cover_type(kidx0:kidx1,itile))
         DO i=1,5                ! Yasso vectors
            LeafLit_coef(1:nidx,itile,i)        = lctlib%LeafLit_coef(surface%cover_type(kidx0:kidx1,itile),i)
            WoodLit_coef(1:nidx,itile,i)        = lctlib%WoodLit_coef(surface%cover_type(kidx0:kidx1,itile),i)
         END DO
      END DO
! HW to scale tau_Cpool_wood
  DO i = 1,nidx
    DO k = 1,ntiles
        IF (relative_extractable_water(i,k) .lt. 0.5) THEN
          tau_Cpool_woods(i,k)=tau_Cpool_woods(i,k)/(-13.33*relative_extractable_water(i,k)+7.67)
        END IF
    END DO
  END DO
  cbalance_diag%tau_Cpool_woods=tau_Cpool_woods! HW
      !! Prepare for test of Carbon conservation across update Cpools by remembering the sum of the Carbon pools
      !!     at the beginning of the day

      IF (debug_Cconservation) THEN
         cbalance_diag%testCconserv(kidx0:kidx1,:) = sumCnatural(cbalance, with_yasso)
      END IF

      ! ---call update N_process() ------

      IF (with_nitrogen) THEN

         IF (debug_Nconservation) nbalance_diag%testNconserv(kidx0:kidx1,:) = sumNpools(nbalance)

         !! distribute N depostion equally to all vegetation types: no scaling with cover_fract
         !!   conversion from [kg(N) m-2(grid box) s-1] to [Mol(N) m-2(canopy) s-1]

         WHERE (surface%is_vegetation(kidx0:kidx1,:))
            ndep_forc_tiles(1:nidx,1:ntiles) =  SPREAD(nbalance%ndep_forc(kidx0:kidx1),DIM=2,NCOPIES=ntiles)    &
                                              / SPREAD(surface%veg_ratio_max(kidx0:kidx1),DIM=2,NCOPIES=ntiles) &
                                              / veg_fract_correction(:,:) / molarMassN_kg
         END WHERE

         call N_process (cbalance%NPP_act_yDayMean(kidx0:kidx1,:),                     &
                         nbalance%nfix_to_sminn(kidx0:kidx1,:),                        &
                         nbalance%ndep_to_sminn(kidx0:kidx1,:),                        &
                         ndep_forc_tiles(1:nidx,1:ntiles),                             &
                         nbalance%nfert_forc_2d(kidx0:kidx1,:),                        &          ! _2d is per tile
                         nbalance%nfert_to_sminn(kidx0:kidx1,:),                       &
                         nbalance%sminn_to_denit(kidx0:kidx1,:),                       &
                         nbalance%sminn_leach(kidx0:kidx1,:),                          &
                         nbalance%SMINN_pool(kidx0:kidx1,:),                           &
                         nbalance%NetEcosyst_N_flux(kidx0:kidx1,:),                    &
                         SPREAD(cbalance_diag%drainage_yDayMean(kidx0:kidx1),DIM=2,NCOPIES=ntiles),  &
                         cbalance_diag%alpha_yDayMean(kidx0:kidx1,:),                  &
                         SPREAD(maxmoisture(:),DIM=2,NCOPIES=ntiles),                  &
                         nbalance%NPP_run_mean(kidx0:kidx1,:),                         &
                         nbalance%N2O_emissions_depfix(kidx0:kidx1,:),                 &
                         surface%is_vegetation(kidx0:kidx1,:) ,                        &
                         nbalance%N2O_emissions_nfert(kidx0:kidx1,:)                   &
                         )
  
         ! update of carbon and nitrogen pools and compute soil respiration rate
                      
         IF (.NOT.with_yasso) THEN
           CALL update_Cpools(cbalance_diag%LAI_yDayMean(kidx0:kidx1,:),          &
                              cbalance_diag%LAI_previousDayMean(kidx0:kidx1,:),   &
                              cbalance_diag%NPP_yDayMean(kidx0:kidx1,:),          &
                              SPREAD(cbalance_diag%topSoilTemp_yDayMean(kidx0:kidx1),DIM=2,NCOPIES=ntiles),  &
                              cbalance_diag%alpha_yDayMean(kidx0:kidx1,:),        &
                              frac_npp_2_woodPool(1:nidx,:),                      &
                              frac_npp_2_reservePool(1:nidx,:),                   &
                              frac_npp_2_exudates(1:nidx,:),                      &
                              frac_green_2_herbivory(1:nidx,:),                   &
                              tau_Cpool_litter_green(1:nidx,:),                   &
                              tau_Cpool_litter_wood(1:nidx,:),                    &
                              tau_Cpool_woods(1:nidx,:),                          &
                              LAI_shed_constant(1:nidx,:),                        &
                              frac_C_litter_green2atmos(1:nidx,:),                &
                              Max_C_content_woods(1:nidx,:),                      &
                              specificLeafArea_C(1:nidx,:),                       &
                              reserveC2leafC(1:nidx,:),                           &
                              LAI_max(1:nidx,:),                                  &
                              surface%is_vegetation(kidx0:kidx1,:),               &
                              surface%is_crop(kidx0:kidx1,:),                     &
                              cbalance%Cpool_green(kidx0:kidx1,:),                &
                              cbalance%Cpool_woods(kidx0:kidx1,:),                &
                              cbalance%Cpool_reserve(kidx0:kidx1,:),              &
                              cbalance%Cpool_crop_harvest(kidx0:kidx1,:),         &
                              cbalance%Cpool_litter_green_ag(kidx0:kidx1,:),      &
                              cbalance%Cpool_litter_green_bg(kidx0:kidx1,:),      &
                              cbalance%Cpool_litter_wood_ag(kidx0:kidx1,:),       &
                              cbalance%Cpool_litter_wood_bg(kidx0:kidx1,:),       &
                              cbalance%Cpool_slow(kidx0:kidx1,:),                 &
                              cbalance%soil_respiration(kidx0:kidx1,:),           &
                              cbalance%soil_respiration_pot(kidx0:kidx1,:),       &
                              cbalance%NPP_flux_correction(kidx0:kidx1,:),        &
                              cbalance%excess_NPP(kidx0:kidx1,:),                 &
                              cbalance%root_exudates(kidx0:kidx1,:),              &
                              cbalance_diag%litter_flux(kidx0:kidx1,:),           &
                              cbalance%Cflux_herbivory(kidx0:kidx1,:),            &
                              cbalance%Cflux_herbivory_LG(kidx0:kidx1,:),         &
                              cbalance%Cflux_herbivory_2_atm(kidx0:kidx1,:),      &
                              cbalance%Cflx_2_crop_harvest(kidx0:kidx1,:),        &
                              cbalance%Cflx_crop_harvest_2_atm(kidx0:kidx1,:),    &
                              cbalance%NPP_act_yDayMean(kidx0:kidx1,:),           &
                              cbalance_diag%leaf_shedding_debug(kidx0:kidx1,:),           & !HW
                              cbalance_diag%Cpool_green_minus_shed_leaves(kidx0:kidx1,:),           & !HW
                              cbalance_diag%excess_carbon_debug(kidx0:kidx1,:),           & !HW
                              cbalance_diag%litter_leaf(kidx0:kidx1,:),           & !HW
                              cbalance_diag%leaf_shedding_rate(kidx0:kidx1,:),           & !HW
    frac_litter_wood_new      =cbalance%frac_litter_wood_new(kidx0:kidx1,:),      &
    redFact_Nlimit            =nbalance%redFact_Nlimitation(kidx0:kidx1,:),        &      
    Npool_green               =nbalance%Npool_green(kidx0:kidx1,:),                &     
    Npool_woods               =nbalance%Npool_woods(kidx0:kidx1,:),                &    
    Npool_mobile              =nbalance%Npool_mobile(kidx0:kidx1,:),               &
    Npool_crop_harvest        =nbalance%Npool_crop_harvest(kidx0:kidx1,:),         &
    Npool_litter_green_ag     =nbalance%Npool_litter_green_ag(kidx0:kidx1,:),      &  
    Npool_litter_green_bg     =nbalance%Npool_litter_green_bg(kidx0:kidx1,:),      & 
    Npool_litter_wood_ag      =nbalance%Npool_litter_wood_ag(kidx0:kidx1,:),       &    
    Npool_litter_wood_bg      =nbalance%Npool_litter_wood_bg(kidx0:kidx1,:),       &      
    Npool_slow                =nbalance%Npool_slow(kidx0:kidx1,:),                 &     
    SMINN_pool                =nbalance%SMINN_pool(kidx0:kidx1,:),                 &    
    minNflux_litter_green_ag  =nbalance%minNflux_litter_green_ag(kidx0:kidx1,:),   &  
    minNflux_litter_green_bg  =nbalance%minNflux_litter_green_bg(kidx0:kidx1,:),   & 
    minNflux_litter_wood_ag   =nbalance%minNflux_litter_wood_ag(kidx0:kidx1,:),    &     
    minNflux_litter_wood_bg   =nbalance%minNflux_litter_wood_bg(kidx0:kidx1,:),    &   
    minNflux_slow             =nbalance%minNflux_slow(kidx0:kidx1,:),              &  
    Nplant_demand             =nbalance%Nplant_demand(kidx0:kidx1,:),              & 
    Nsoil_demand              =nbalance%Nsoil_demand(kidx0:kidx1,:),               &     
    Ntotal_demand             =nbalance%Ntotal_demand(kidx0:kidx1,:),              &    
    SMINN_herbivory           =nbalance%SMINN_herbivory(kidx0:kidx1,:),            &  
    N2O_emissions_mineraliz   =nbalance%N2O_emissions_mineraliz(kidx0:kidx1,:) ,   &
    N2O_emissions_slow        =nbalance%N2O_emissions_slow(kidx0:kidx1,:),         &
    N2O_emissions_grazing     =nbalance%N2O_emissions_grazing(kidx0:kidx1,:),      &
    Nflx_2_crop_harvest       =nbalance%Nflx_2_crop_harvest(kidx0:kidx1,:)         &
                              )
         ELSE
            CALL update_Cpools(cbalance_diag%LAI_yDayMean(kidx0:kidx1,:),          &
                               cbalance_diag%LAI_previousDayMean(kidx0:kidx1,:),   &
                               cbalance_diag%NPP_yDayMean(kidx0:kidx1,:),          &
                               SPREAD(cbalance_diag%topSoilTemp_yDayMean(kidx0:kidx1),DIM=2,NCOPIES=ntiles),  &
                               cbalance_diag%alpha_yDayMean(kidx0:kidx1,:),        &
                               frac_npp_2_woodPool(1:nidx,:),                      &
                               frac_npp_2_reservePool(1:nidx,:),                   &
                               frac_npp_2_exudates(1:nidx,:),                      &
                               frac_green_2_herbivory(1:nidx,:),                   &
                               tau_Cpool_litter_green(1:nidx,:),                   &
                               tau_Cpool_litter_wood(1:nidx,:),                    &
                               tau_Cpool_woods(1:nidx,:),                          &
                               LAI_shed_constant(1:nidx,:),                        &
                               frac_C_litter_green2atmos(1:nidx,:),                &
                               Max_C_content_woods(1:nidx,:),                      &
                               specificLeafArea_C(1:nidx,:),                       &
                               reserveC2leafC(1:nidx,:),                           &
                               LAI_max(1:nidx,:),                                  &
                               surface%is_vegetation(kidx0:kidx1,:),               &
                               surface%is_crop(kidx0:kidx1,:),                     &
                               cbalance%Cpool_green(kidx0:kidx1,:),                &
                               cbalance%Cpool_woods(kidx0:kidx1,:),                &
                               cbalance%Cpool_reserve(kidx0:kidx1,:),              &
                               cbalance%Cpool_crop_harvest(kidx0:kidx1,:),         &
                               cbalance%Cpool_litter_green_ag(kidx0:kidx1,:),      &
                               cbalance%Cpool_litter_green_bg(kidx0:kidx1,:),      &
                               cbalance%Cpool_litter_wood_ag(kidx0:kidx1,:),       &
                               cbalance%Cpool_litter_wood_bg(kidx0:kidx1,:),       &
                               cbalance%Cpool_slow(kidx0:kidx1,:),                 &
                               cbalance%soil_respiration(kidx0:kidx1,:),           &
                               cbalance%soil_respiration_pot(kidx0:kidx1,:),       &
                               cbalance%NPP_flux_correction(kidx0:kidx1,:),        &
                               cbalance%excess_NPP(kidx0:kidx1,:),                 &
                               cbalance%root_exudates(kidx0:kidx1,:),              &
                               cbalance_diag%litter_flux(kidx0:kidx1,:),           &
                               cbalance%Cflux_herbivory(kidx0:kidx1,:),            &
                               cbalance%Cflux_herbivory_LG(kidx0:kidx1,:),         &
                               cbalance%Cflux_herbivory_2_atm(kidx0:kidx1,:),      &
                               cbalance%Cflx_2_crop_harvest(kidx0:kidx1,:),        &
                               cbalance%Cflx_crop_harvest_2_atm(kidx0:kidx1,:),    &
                               cbalance%NPP_act_yDayMean(kidx0:kidx1,:),           &
                              cbalance_diag%leaf_shedding_debug(kidx0:kidx1,:),           & !HW
                              cbalance_diag%Cpool_green_minus_shed_leaves(kidx0:kidx1,:),           & !HW
                              cbalance_diag%excess_carbon_debug(kidx0:kidx1,:),           & !HW
                               cbalance_diag%litter_leaf(kidx0:kidx1,:),           & !HW
                              cbalance_diag%leaf_shedding_rate(kidx0:kidx1,:),           & !HW
    frac_litter_wood_new    =  cbalance%frac_litter_wood_new(kidx0:kidx1,:),       &
                     ! variables only needed with yasso   
    temp2_30d               =  SPREAD(cbalance_diag%pseudo_temp_yDay(kidx0:kidx1),DIM=2,NCOPIES=ntiles),    &
    precip_30d              =  SPREAD(cbalance_diag%pseudo_precip_yDay(kidx0:kidx1),DIM=2,NCOPIES=ntiles),  &
    YCpool_acid_ag1         =  cbalance%YCpool_acid_ag1(kidx0:kidx1,:),             &
    YCpool_water_ag1        =  cbalance%YCpool_water_ag1(kidx0:kidx1,:),            &
    YCpool_ethanol_ag1      =  cbalance%YCpool_ethanol_ag1(kidx0:kidx1,:),          &
    YCpool_nonsoluble_ag1   =  cbalance%YCpool_nonsoluble_ag1(kidx0:kidx1,:),       &
    YCpool_acid_bg1         =  cbalance%YCpool_acid_bg1(kidx0:kidx1,:),             &
    YCpool_water_bg1        =  cbalance%YCpool_water_bg1(kidx0:kidx1,:),            &
    YCpool_ethanol_bg1      =  cbalance%YCpool_ethanol_bg1(kidx0:kidx1,:),          &
    YCpool_nonsoluble_bg1   =  cbalance%YCpool_nonsoluble_bg1(kidx0:kidx1,:),       &
    YCpool_humus_1          =  cbalance%YCpool_humus_1(kidx0:kidx1,:),              &
    YCpool_acid_ag2         =  cbalance%YCpool_acid_ag2(kidx0:kidx1,:),             &
    YCpool_water_ag2        =  cbalance%YCpool_water_ag2(kidx0:kidx1,:),            &
    YCpool_ethanol_ag2      =  cbalance%YCpool_ethanol_ag2(kidx0:kidx1,:),          &
    YCpool_nonsoluble_ag2   =  cbalance%YCpool_nonsoluble_ag2(kidx0:kidx1,:),       &
    YCpool_acid_bg2         =  cbalance%YCpool_acid_bg2(kidx0:kidx1,:),             &
    YCpool_water_bg2        =  cbalance%YCpool_water_bg2(kidx0:kidx1,:),            &
    YCpool_ethanol_bg2      =  cbalance%YCpool_ethanol_bg2(kidx0:kidx1,:),          &
    YCpool_nonsoluble_bg2   =  cbalance%YCpool_nonsoluble_bg2(kidx0:kidx1,:),       &
    YCpool_humus_2          =  cbalance%YCpool_humus_2(kidx0:kidx1,:),              &
    LeafLit_coef1           =  LeafLit_coef(1:nidx,:,1),                           &         ! Yasso vectors
    LeafLit_coef2           =  LeafLit_coef(1:nidx,:,2),                           &
    LeafLit_coef3           =  LeafLit_coef(1:nidx,:,3),                           &
    LeafLit_coef4           =  LeafLit_coef(1:nidx,:,4),                           &
    LeafLit_coef5           =  LeafLit_coef(1:nidx,:,5),                           &
    WoodLit_coef1           =  WoodLit_coef(1:nidx,:,1),                           &
    WoodLit_coef2           =  WoodLit_coef(1:nidx,:,2),                           &
    WoodLit_coef3           =  WoodLit_coef(1:nidx,:,3),                           &
    WoodLit_coef4           =  WoodLit_coef(1:nidx,:,4),                           &
    WoodLit_coef5           =  WoodLit_coef(1:nidx,:,5),                           &     
    WoodLitterSize          =  WoodLitterSize(1:nidx,:),                           &
    ! variables needed for nitrogen cycle
    redFact_Nlimit            =nbalance%redFact_Nlimitation(kidx0:kidx1,:),        &      
    Npool_green               =nbalance%Npool_green(kidx0:kidx1,:),                &     
    Npool_woods               =nbalance%Npool_woods(kidx0:kidx1,:),                &    
    Npool_mobile              =nbalance%Npool_mobile(kidx0:kidx1,:),               &
    Npool_crop_harvest        =nbalance%Npool_crop_harvest(kidx0:kidx1,:),         &
    Npool_litter_green_ag     =nbalance%Npool_litter_green_ag(kidx0:kidx1,:),      &  
    Npool_litter_green_bg     =nbalance%Npool_litter_green_bg(kidx0:kidx1,:),      & 
    Npool_litter_wood_ag      =nbalance%Npool_litter_wood_ag(kidx0:kidx1,:),       &    
    Npool_litter_wood_bg      =nbalance%Npool_litter_wood_bg(kidx0:kidx1,:),       &      
    Npool_slow                =nbalance%Npool_slow(kidx0:kidx1,:),                 &     
    SMINN_pool                =nbalance%SMINN_pool(kidx0:kidx1,:),                 &    
    minNflux_litter_green_ag  =nbalance%minNflux_litter_green_ag(kidx0:kidx1,:),   &  
    minNflux_litter_green_bg  =nbalance%minNflux_litter_green_bg(kidx0:kidx1,:),   & 
    minNflux_litter_wood_ag   =nbalance%minNflux_litter_wood_ag(kidx0:kidx1,:),    &     
    minNflux_litter_wood_bg   =nbalance%minNflux_litter_wood_bg(kidx0:kidx1,:),    &   
    minNflux_slow             =nbalance%minNflux_slow(kidx0:kidx1,:),              &  
    Nplant_demand             =nbalance%Nplant_demand(kidx0:kidx1,:),              & 
    Nsoil_demand              =nbalance%Nsoil_demand(kidx0:kidx1,:),               &     
    Ntotal_demand             =nbalance%Ntotal_demand(kidx0:kidx1,:),              &    
    SMINN_herbivory           =nbalance%SMINN_herbivory(kidx0:kidx1,:),            &  
    N2O_emissions_mineraliz   =nbalance%N2O_emissions_mineraliz(kidx0:kidx1,:) ,   &
    N2O_emissions_slow        =nbalance%N2O_emissions_slow(kidx0:kidx1,:),         &
    N2O_emissions_grazing     =nbalance%N2O_emissions_grazing(kidx0:kidx1,:),      &
    Nflx_2_crop_harvest       =nbalance%Nflx_2_crop_harvest(kidx0:kidx1,:)         &
                               )
         ENDIF !(with_yasso)

      ELSE !! carbon-only mode

         IF (with_yasso) THEN
            CALL update_Cpools(cbalance_diag%LAI_yDayMean(kidx0:kidx1,:),          &
                               cbalance_diag%LAI_previousDayMean(kidx0:kidx1,:),   &
                               cbalance_diag%NPP_yDayMean(kidx0:kidx1,:),          &
                               SPREAD(cbalance_diag%topSoilTemp_yDayMean(kidx0:kidx1),DIM=2,NCOPIES=ntiles),  &
                               cbalance_diag%alpha_yDayMean(kidx0:kidx1,:),        &
                               frac_npp_2_woodPool(1:nidx,:),                      &
                               frac_npp_2_reservePool(1:nidx,:),                   &
                               frac_npp_2_exudates(1:nidx,:),                      &
                               frac_green_2_herbivory(1:nidx,:),                   &
                               tau_Cpool_litter_green(1:nidx,:),                   &
                               tau_Cpool_litter_wood(1:nidx,:),                    &
                               tau_Cpool_woods(1:nidx,:),                          &
                               LAI_shed_constant(1:nidx,:),                        &
                               frac_C_litter_green2atmos(1:nidx,:),                &
                               Max_C_content_woods(1:nidx,:),                      &
                               specificLeafArea_C(1:nidx,:),                       &
                               reserveC2leafC(1:nidx,:),                           &
                               LAI_max(1:nidx,:),                                  &
                               surface%is_vegetation(kidx0:kidx1,:),               &
                               surface%is_crop(kidx0:kidx1,:),                     &
                               cbalance%Cpool_green(kidx0:kidx1,:),                &
                               cbalance%Cpool_woods(kidx0:kidx1,:),                &
                               cbalance%Cpool_reserve(kidx0:kidx1,:),              &
                               cbalance%Cpool_crop_harvest(kidx0:kidx1,:),         &
                               cbalance%Cpool_litter_green_ag(kidx0:kidx1,:),      &
                               cbalance%Cpool_litter_green_bg(kidx0:kidx1,:),      &
                               cbalance%Cpool_litter_wood_ag(kidx0:kidx1,:),       &
                               cbalance%Cpool_litter_wood_bg(kidx0:kidx1,:),       &
                               cbalance%Cpool_slow(kidx0:kidx1,:),                 &
                               cbalance%soil_respiration(kidx0:kidx1,:),           &
                               cbalance%soil_respiration_pot(kidx0:kidx1,:),       &
                               cbalance%NPP_flux_correction(kidx0:kidx1,:),        &
                               cbalance%excess_NPP(kidx0:kidx1,:),                 &
                               cbalance%root_exudates(kidx0:kidx1,:),              &
                               cbalance_diag%litter_flux(kidx0:kidx1,:),           &
                               cbalance%Cflux_herbivory(kidx0:kidx1,:),            &
                               cbalance%Cflux_herbivory_LG(kidx0:kidx1,:),         &
                               cbalance%Cflux_herbivory_2_atm(kidx0:kidx1,:),      &
                               cbalance%Cflx_2_crop_harvest(kidx0:kidx1,:),        &
                               cbalance%Cflx_crop_harvest_2_atm(kidx0:kidx1,:),    &
                               cbalance%NPP_act_yDayMean(kidx0:kidx1,:),           &
                              cbalance_diag%leaf_shedding_debug(kidx0:kidx1,:),           & !HW
                              cbalance_diag%Cpool_green_minus_shed_leaves(kidx0:kidx1,:),           & !HW
                              cbalance_diag%excess_carbon_debug(kidx0:kidx1,:),           & !HW
                               cbalance_diag%litter_leaf(kidx0:kidx1,:),           & !HW
                              cbalance_diag%leaf_shedding_rate(kidx0:kidx1,:),           & !HW
    frac_litter_wood_new    =  cbalance%frac_litter_wood_new(kidx0:kidx1,:),       &
                     ! variables only needed with yasso   
    temp2_30d               =  SPREAD(cbalance_diag%pseudo_temp_yDay(kidx0:kidx1),DIM=2,NCOPIES=ntiles),    &
    precip_30d              =  SPREAD(cbalance_diag%pseudo_precip_yDay(kidx0:kidx1),DIM=2,NCOPIES=ntiles),  &
    YCpool_acid_ag1         =  cbalance%YCpool_acid_ag1(kidx0:kidx1,:),             &
    YCpool_water_ag1        =  cbalance%YCpool_water_ag1(kidx0:kidx1,:),            &
    YCpool_ethanol_ag1      =  cbalance%YCpool_ethanol_ag1(kidx0:kidx1,:),          &
    YCpool_nonsoluble_ag1   =  cbalance%YCpool_nonsoluble_ag1(kidx0:kidx1,:),       &
    YCpool_acid_bg1         =  cbalance%YCpool_acid_bg1(kidx0:kidx1,:),             &
    YCpool_water_bg1        =  cbalance%YCpool_water_bg1(kidx0:kidx1,:),            &
    YCpool_ethanol_bg1      =  cbalance%YCpool_ethanol_bg1(kidx0:kidx1,:),          &
    YCpool_nonsoluble_bg1   =  cbalance%YCpool_nonsoluble_bg1(kidx0:kidx1,:),       &
    YCpool_humus_1          =  cbalance%YCpool_humus_1(kidx0:kidx1,:),              &
    YCpool_acid_ag2         =  cbalance%YCpool_acid_ag2(kidx0:kidx1,:),             &
    YCpool_water_ag2        =  cbalance%YCpool_water_ag2(kidx0:kidx1,:),            &
    YCpool_ethanol_ag2      =  cbalance%YCpool_ethanol_ag2(kidx0:kidx1,:),          &
    YCpool_nonsoluble_ag2   =  cbalance%YCpool_nonsoluble_ag2(kidx0:kidx1,:),       &
    YCpool_acid_bg2         =  cbalance%YCpool_acid_bg2(kidx0:kidx1,:),             &
    YCpool_water_bg2        =  cbalance%YCpool_water_bg2(kidx0:kidx1,:),            &
    YCpool_ethanol_bg2      =  cbalance%YCpool_ethanol_bg2(kidx0:kidx1,:),          &
    YCpool_nonsoluble_bg2   =  cbalance%YCpool_nonsoluble_bg2(kidx0:kidx1,:),       &
    YCpool_humus_2          =  cbalance%YCpool_humus_2(kidx0:kidx1,:),              &
    LeafLit_coef1           =  LeafLit_coef(1:nidx,:,1),                           &         ! Yasso vectors
    LeafLit_coef2           =  LeafLit_coef(1:nidx,:,2),                           &
    LeafLit_coef3           =  LeafLit_coef(1:nidx,:,3),                           &
    LeafLit_coef4           =  LeafLit_coef(1:nidx,:,4),                           &
    LeafLit_coef5           =  LeafLit_coef(1:nidx,:,5),                           &
    WoodLit_coef1           =  WoodLit_coef(1:nidx,:,1),                           &
    WoodLit_coef2           =  WoodLit_coef(1:nidx,:,2),                           &
    WoodLit_coef3           =  WoodLit_coef(1:nidx,:,3),                           &
    WoodLit_coef4           =  WoodLit_coef(1:nidx,:,4),                           &
    WoodLit_coef5           =  WoodLit_coef(1:nidx,:,5),                           &     
    WoodLitterSize          =  WoodLitterSize(1:nidx,:)                            &
                               )
         ELSE !! without yasso
            CALL update_Cpools(cbalance_diag%LAI_yDayMean(kidx0:kidx1,:),          &
                               cbalance_diag%LAI_previousDayMean(kidx0:kidx1,:),   &
                               cbalance_diag%NPP_yDayMean(kidx0:kidx1,:),          &
                               SPREAD(cbalance_diag%topSoilTemp_yDayMean(kidx0:kidx1),DIM=2,NCOPIES=ntiles),  &
                               cbalance_diag%alpha_yDayMean(kidx0:kidx1,:),        &
                               frac_npp_2_woodPool(1:nidx,:),                      &
                               frac_npp_2_reservePool(1:nidx,:),                   &
                               frac_npp_2_exudates(1:nidx,:),                      &
                               frac_green_2_herbivory(1:nidx,:),                   &
                               tau_Cpool_litter_green(1:nidx,:),                   &
                               tau_Cpool_litter_wood(1:nidx,:),                    &
                               tau_Cpool_woods(1:nidx,:),                          &
                               LAI_shed_constant(1:nidx,:),                        &
                               frac_C_litter_green2atmos(1:nidx,:),                &
                               Max_C_content_woods(1:nidx,:),                      &
                               specificLeafArea_C(1:nidx,:),                       &
                               reserveC2leafC(1:nidx,:),                           &
                               LAI_max(1:nidx,:),                                  &
                               surface%is_vegetation(kidx0:kidx1,:),               &
                               surface%is_crop(kidx0:kidx1,:),                     &
                               cbalance%Cpool_green(kidx0:kidx1,:),                &
                               cbalance%Cpool_woods(kidx0:kidx1,:),                &
                               cbalance%Cpool_reserve(kidx0:kidx1,:),              &
                               cbalance%Cpool_crop_harvest(kidx0:kidx1,:),         &
                               cbalance%Cpool_litter_green_ag(kidx0:kidx1,:),      &
                               cbalance%Cpool_litter_green_bg(kidx0:kidx1,:),      &
                               cbalance%Cpool_litter_wood_ag(kidx0:kidx1,:),       &
                               cbalance%Cpool_litter_wood_bg(kidx0:kidx1,:),       &
                               cbalance%Cpool_slow(kidx0:kidx1,:),                 &
                               cbalance%soil_respiration(kidx0:kidx1,:),           &
                               cbalance%soil_respiration_pot(kidx0:kidx1,:),       &
                               cbalance%NPP_flux_correction(kidx0:kidx1,:),        &
                               cbalance%excess_NPP(kidx0:kidx1,:),                 &
                               cbalance%root_exudates(kidx0:kidx1,:),              &
                               cbalance_diag%litter_flux(kidx0:kidx1,:),           &
                               cbalance%Cflux_herbivory(kidx0:kidx1,:),            &
                               cbalance%Cflux_herbivory_LG(kidx0:kidx1,:),         &
                               cbalance%Cflux_herbivory_2_atm(kidx0:kidx1,:),      &
                               cbalance%Cflx_2_crop_harvest(kidx0:kidx1,:),        &
                               cbalance%Cflx_crop_harvest_2_atm(kidx0:kidx1,:),    &
                               cbalance%NPP_act_yDayMean(kidx0:kidx1,:),           &
                              cbalance_diag%leaf_shedding_debug(kidx0:kidx1,:),           & !HW
                              cbalance_diag%Cpool_green_minus_shed_leaves(kidx0:kidx1,:),           & !HW
                              cbalance_diag%excess_carbon_debug(kidx0:kidx1,:),           & !HW
                               cbalance_diag%litter_leaf(kidx0:kidx1,:),           & !HW
                              cbalance_diag%leaf_shedding_rate(kidx0:kidx1,:),           & !HW
       frac_litter_wood_new  = cbalance%frac_litter_wood_new(kidx0:kidx1,:)        &
                               )

         END IF
      ENDIF

      !! Update Carbon and Nitrogen conservation test (should be zero if everything is correct, except for green litter pool)

      IF (debug_Cconservation) THEN
         cbalance_diag%testCconserv(kidx0:kidx1,:) =                     &
            cbalance_diag%testCconserv(kidx0:kidx1,:)                    & !! Sum of pools at beginning of day 
            + cbalance%NPP_act_yDayMean(kidx0:kidx1,:)      * 86400._dp  & !!    plus NPP added 
            + cbalance%soil_respiration(kidx0:kidx1,:)      * 86400._dp  & !!    minus carbon lost to atmosph.
            + cbalance%Cflux_herbivory_2_atm(kidx0:kidx1,:) * 86400._dp  & !!
            - sumCnatural(cbalance, with_yasso)                            !!    minus the pools at end of day
      END IF                                                               !! should always give zero!

      IF (with_nitrogen .AND. debug_Nconservation) THEN
         nbalance_diag%testNconserv(kidx0:kidx1,:) =                         &
              nbalance_diag%testNconserv(kidx0:kidx1,:)                      &
              + nbalance%NetEcosyst_N_flux(kidx0:kidx1,:)       * 86400._dp  &
              - nbalance%N2O_emissions_mineraliz(kidx0:kidx1,:) * 86400._dp  &
              - nbalance%N2O_emissions_grazing(kidx0:kidx1,:)   * 86400._dp  &
              - sumNpools(nbalance)
      END IF

      ! Compute carbon contents per grid box area by weighting pools with fractions of grid box covered by vegetation 
      ! -------------------------------------------------------------------------------------------------------------
      
      cbalance_diag%boxC_green(kidx0:kidx1,:)           = cbalance%Cpool_green(kidx0:kidx1,:)        * areaWeightingFactor(:,:)
      cbalance_diag%boxC_reserve(kidx0:kidx1,:)         = cbalance%Cpool_reserve(kidx0:kidx1,:)      * areaWeightingFactor(:,:)
      cbalance_diag%boxC_woods(kidx0:kidx1,:)           = cbalance%Cpool_woods(kidx0:kidx1,:)        * areaWeightingFactor(:,:)
      cbalance_diag%boxC_crop_harvest(kidx0:kidx1,:)    = cbalance%Cpool_crop_harvest(kidx0:kidx1,:) * areaWeightingFactor(:,:)

      IF (.NOT. with_yasso) THEN
         cbalance_diag%boxC_litter_green_ag(kidx0:kidx1,:)= cbalance%Cpool_litter_green_ag(kidx0:kidx1,:) * areaWeightingFactor(:,:)
         cbalance_diag%boxC_litter_green_bg(kidx0:kidx1,:)= cbalance%Cpool_litter_green_bg(kidx0:kidx1,:) * areaWeightingFactor(:,:)
         cbalance_diag%boxC_litter_wood(kidx0:kidx1,:)    = &
            (cbalance%Cpool_litter_wood_ag(kidx0:kidx1,:) + cbalance%Cpool_litter_wood_bg(kidx0:kidx1,:)) * areaWeightingFactor(:,:)
         cbalance_diag%boxC_slow(kidx0:kidx1,:)           = cbalance%Cpool_slow(kidx0:kidx1,:)            * areaWeightingFactor(:,:)
      ELSE
         cbalance_diag%boxYC_acid_ag1(kidx0:kidx1,:)      = cbalance%YCpool_acid_ag1(kidx0:kidx1,:)       * areaWeightingFactor(:,:)
         cbalance_diag%boxYC_water_ag1(kidx0:kidx1,:)     = cbalance%YCpool_water_ag1(kidx0:kidx1,:)      * areaWeightingFactor(:,:)
         cbalance_diag%boxYC_ethanol_ag1(kidx0:kidx1,:)   = cbalance%YCpool_ethanol_ag1(kidx0:kidx1,:)    * areaWeightingFactor(:,:)
         cbalance_diag%boxYC_nonsoluble_ag1(kidx0:kidx1,:)= cbalance%YCpool_nonsoluble_ag1(kidx0:kidx1,:) * areaWeightingFactor(:,:)
         cbalance_diag%boxYC_acid_bg1(kidx0:kidx1,:)      = cbalance%YCpool_acid_bg1(kidx0:kidx1,:)       * areaWeightingFactor(:,:)
         cbalance_diag%boxYC_water_bg1(kidx0:kidx1,:)     = cbalance%YCpool_water_bg1(kidx0:kidx1,:)      * areaWeightingFactor(:,:)
         cbalance_diag%boxYC_ethanol_bg1(kidx0:kidx1,:)   = cbalance%YCpool_ethanol_bg1(kidx0:kidx1,:)    * areaWeightingFactor(:,:)
         cbalance_diag%boxYC_nonsoluble_bg1(kidx0:kidx1,:)= cbalance%YCpool_nonsoluble_bg1(kidx0:kidx1,:) * areaWeightingFactor(:,:)
         cbalance_diag%boxYC_humus_1(kidx0:kidx1,:)       = cbalance%YCpool_humus_1(kidx0:kidx1,:)        * areaWeightingFactor(:,:)
         cbalance_diag%boxYC_acid_ag2(kidx0:kidx1,:)      = cbalance%YCpool_acid_ag2(kidx0:kidx1,:)       * areaWeightingFactor(:,:)
         cbalance_diag%boxYC_water_ag2(kidx0:kidx1,:)     = cbalance%YCpool_water_ag2(kidx0:kidx1,:)      * areaWeightingFactor(:,:)
         cbalance_diag%boxYC_ethanol_ag2(kidx0:kidx1,:)   = cbalance%YCpool_ethanol_ag2(kidx0:kidx1,:)    * areaWeightingFactor(:,:)
         cbalance_diag%boxYC_nonsoluble_ag2(kidx0:kidx1,:)= cbalance%YCpool_nonsoluble_ag2(kidx0:kidx1,:) * areaWeightingFactor(:,:)
         cbalance_diag%boxYC_acid_bg2(kidx0:kidx1,:)      = cbalance%YCpool_acid_bg2(kidx0:kidx1,:)       * areaWeightingFactor(:,:)
         cbalance_diag%boxYC_water_bg2(kidx0:kidx1,:)     = cbalance%YCpool_water_bg2(kidx0:kidx1,:)      * areaWeightingFactor(:,:)
         cbalance_diag%boxYC_ethanol_bg2(kidx0:kidx1,:)   = cbalance%YCpool_ethanol_bg2(kidx0:kidx1,:)    * areaWeightingFactor(:,:)
         cbalance_diag%boxYC_nonsoluble_bg2(kidx0:kidx1,:)= cbalance%YCpool_nonsoluble_bg2(kidx0:kidx1,:) * areaWeightingFactor(:,:)
         cbalance_diag%boxYC_humus_2(kidx0:kidx1,:)       = cbalance%YCpool_humus_2(kidx0:kidx1,:)        * areaWeightingFactor(:,:)
      END IF

      cbalance_diag%box_Cpools_total(kidx0:kidx1) = SUM(sumCnatural(cbalance, with_yasso) * areaWeightingFactor(:,:), DIM=2)
      IF (lcc_scheme==2) THEN
         cbalance_diag%box_Cpools_total(kidx0:kidx1) = cbalance_diag%box_Cpools_total(kidx0:kidx1)     &
              + surface%veg_ratio_max(kidx0:kidx1) * sumCanthro(cbalance, with_landcover_transitions)
      END IF
      
      IF (lcc_scheme==2) THEN

         ! Compute averages of anthropogenic pools if necessary
         ! ------------------------------------------------------
         cbalance_diag%boxC_onSite_avg      (kidx0:kidx1) = cbalance_diag%boxC_onSite_avg      (kidx0:kidx1) &
            + cbalance%Cpool_onSite         (kidx0:kidx1) * 86400._dp  *  surface%veg_ratio_max(kidx0:kidx1)
         cbalance_diag%boxC_paper_avg       (kidx0:kidx1) = cbalance_diag%boxC_paper_avg       (kidx0:kidx1) &
            + cbalance%Cpool_paper          (kidx0:kidx1) * 86400._dp  *  surface%veg_ratio_max(kidx0:kidx1)
         cbalance_diag%boxC_construction_avg(kidx0:kidx1) = cbalance_diag%boxC_construction_avg(kidx0:kidx1) &
            + cbalance%Cpool_construction   (kidx0:kidx1) * 86400._dp  *  surface%veg_ratio_max(kidx0:kidx1)
         IF (with_landcover_transitions) THEN
            cbalance_diag%boxC_onSite_harvest_avg      (kidx0:kidx1) = cbalance_diag%boxC_onSite_harvest_avg      (kidx0:kidx1) &
               + cbalance%Cpool_onSite_harvest         (kidx0:kidx1) * 86400._dp  *  surface%veg_ratio_max        (kidx0:kidx1)
            cbalance_diag%boxC_paper_harvest_avg       (kidx0:kidx1) = cbalance_diag%boxC_paper_harvest_avg       (kidx0:kidx1) &
               + cbalance%Cpool_paper_harvest          (kidx0:kidx1) * 86400._dp  *  surface%veg_ratio_max        (kidx0:kidx1)
            cbalance_diag%boxC_construction_harvest_avg(kidx0:kidx1) = cbalance_diag%boxC_construction_harvest_avg(kidx0:kidx1) &
               + cbalance%Cpool_construction_harvest   (kidx0:kidx1) * 86400._dp  *  surface%veg_ratio_max        (kidx0:kidx1)
         ENDIF
      ENDIF

      ! Add FOM product pools to box_Cpools_total if required
      IF (with_forest_management) THEN
        cbalance_diag%box_Cpools_total(kidx0:kidx1) = cbalance_diag%box_Cpools_total(kidx0:kidx1) &
          & + cbalance%FOM_fm_Box_burning_pool(kidx0:kidx1) + cbalance%FOM_fm_Box_paper_pool(kidx0:kidx1) &
          & + cbalance%FOM_fm_Box_construction_pool(kidx0:kidx1) + cbalance%FOM_fm_Box_eternity_pool(kidx0:kidx1)
      ENDIF

      ! Other diagnostics
      ! -----------------
      cbalance_diag%box_soil_respiration(kidx0:kidx1,:)   = cbalance%soil_respiration(kidx0:kidx1,:)    * areaWeightingFactor(:,:)
      cbalance_diag%box_soil_respiration_pot(kidx0:kidx1,:) = cbalance%soil_respiration_pot(kidx0:kidx1,:) &
                                                              * areaWeightingFactor(:,:)
      cbalance_diag%box_NPP_yDayMean(kidx0:kidx1,:)       = cbalance_diag%NPP_yDayMean(kidx0:kidx1,:)   * areaWeightingFactor(:,:)
      cbalance_diag%box_NPP_act_yDayMean(kidx0:kidx1,:)   = cbalance%NPP_act_yDayMean(kidx0:kidx1,:)    * areaWeightingFactor(:,:)
      cbalance_diag%box_NPP_flux_correction(kidx0:kidx1,:)= cbalance%NPP_flux_correction(kidx0:kidx1,:) * areaWeightingFactor(:,:)
      cbalance_diag%box_GPP_yDayMean(kidx0:kidx1,:)       = cbalance_diag%GPP_yDayMean(kidx0:kidx1,:)   * areaWeightingFactor(:,:)
      cbalance_diag%box_litter_flux(kidx0:kidx1,:)        = cbalance_diag%litter_flux(kidx0:kidx1,:)    * areaWeightingFactor(:,:)
      cbalance_diag%box_litter_leaf(kidx0:kidx1,:)        = cbalance_diag%litter_leaf(kidx0:kidx1,:)    * areaWeightingFactor(:,:) ! HW
      cbalance_diag%box_root_exudates(kidx0:kidx1,:)      = cbalance%root_exudates(kidx0:kidx1,:)       * areaWeightingFactor(:,:)
      cbalance_diag%box_Cflux_herbivory(kidx0:kidx1,:)    = cbalance%Cflux_herbivory(kidx0:kidx1,:)     * areaWeightingFactor(:,:)
      cbalance_diag%box_Cflux_herbivory_2_atm(kidx0:kidx1,:) = cbalance%Cflux_herbivory_2_atm(kidx0:kidx1,:) &
                                                              * areaWeightingFactor(:,:)
      cbalance_diag%box_Cflx_2_crop_harvest(kidx0:kidx1)  = SUM(cbalance%Cflx_2_crop_harvest(kidx0:kidx1,:) &
                                                           * areaWeightingFactor(:,:), DIM=2)
      cbalance_diag%box_Cflx_crop_harvest_2_atm(kidx0:kidx1) = SUM(cbalance%Cflx_crop_harvest_2_atm(kidx0:kidx1,:) &
                                                              * areaWeightingFactor(:,:), DIM=2)

      IF (with_nitrogen) THEN

         ! Compute nitrogen contents per grid box area by weighting pools with fractions of grid box covered by vegetation
         ! -------------------------------------------------------------------------------------------------------------
         nbalance_diag%boxN_green(kidx0:kidx1,:)          = nbalance%Npool_green(kidx0:kidx1,:)          * areaWeightingFactor(:,:)
         nbalance_diag%boxN_mobile(kidx0:kidx1,:)         = nbalance%Npool_mobile(kidx0:kidx1,:)         * areaWeightingFactor(:,:)
         nbalance_diag%boxN_woods(kidx0:kidx1,:)          = nbalance%Npool_woods(kidx0:kidx1,:)          * areaWeightingFactor(:,:)
         nbalance_diag%boxN_crop_harvest(kidx0:kidx1,:)   = nbalance%Npool_crop_harvest(kidx0:kidx1,:)   * areaWeightingFactor(:,:)
         nbalance_diag%boxN_litter_green_ag(kidx0:kidx1,:)= nbalance%Npool_litter_green_ag(kidx0:kidx1,:)* areaWeightingFactor(:,:)
         nbalance_diag%boxN_litter_green_bg(kidx0:kidx1,:)= nbalance%Npool_litter_green_bg(kidx0:kidx1,:)* areaWeightingFactor(:,:)
         nbalance_diag%boxN_litter_wood_ag(kidx0:kidx1,:) = nbalance%Npool_litter_wood_ag(kidx0:kidx1,:) * areaWeightingFactor(:,:)
         nbalance_diag%boxN_litter_wood_bg(kidx0:kidx1,:) = nbalance%Npool_litter_wood_bg(kidx0:kidx1,:) * areaWeightingFactor(:,:)
         nbalance_diag%boxN_slow(kidx0:kidx1,:)           = nbalance%Npool_slow(kidx0:kidx1,:)           * areaWeightingFactor(:,:)
         nbalance_diag%boxN_SMINN(kidx0:kidx1,:)          = nbalance%SMINN_pool(kidx0:kidx1,:)           * areaWeightingFactor(:,:)

         nbalance_diag%box_N2O_total(kidx0:kidx1,:) = (nbalance%N2O_emissions_mineraliz(kidx0:kidx1,:) +  &
                                                       nbalance%N2O_emissions_depfix(kidx0:kidx1,:) +     &
                                                       nbalance%N2O_emissions_grazing(kidx0:kidx1,:) +    &
                                                       nbalance%N2O_emissions_nfert(kidx0:kidx1,:) )      &
                                                       * areaWeightingFactor(:,:) * 1.e06_dp
         ! factor 1.e06_dp: conversion from Mol to mikro-Mol to overcome Grib problem with numbers smaller 1.e-13

         nbalance_diag%box_N2O_emissions_depfix(kidx0:kidx1,:) = nbalance%N2O_emissions_depfix(kidx0:kidx1,:) &
                                                                 * areaWeightingFactor(:,:) * 1.e06_dp

         nbalance_diag%box_N2O_emissions_mineraliz(kidx0:kidx1,:) = nbalance%N2O_emissions_mineraliz(kidx0:kidx1,:) &
                                                                    * areaWeightingFactor(:,:) * 1.e06_dp

         nbalance_diag%box_minNflux_litter_green(kidx0:kidx1,:) = &
              (nbalance%minNflux_litter_green_ag(kidx0:kidx1,:) + nbalance%minNflux_litter_green_bg(kidx0:kidx1,:)) &
               * areaWeightingFactor(:,:)
         nbalance_diag%box_minNflux_litter_wood(kidx0:kidx1,:) = &
              (nbalance%minNflux_litter_wood_ag(kidx0:kidx1,:) + nbalance%minNflux_litter_wood_bg(kidx0:kidx1,:)) &
               * areaWeightingFactor(:,:)
         nbalance_diag%box_minNflux_slow(kidx0:kidx1,:)   = nbalance%minNflux_slow(kidx0:kidx1,:)   * areaWeightingFactor(:,:)
         nbalance_diag%box_Nplant_demand(kidx0:kidx1,:)   = nbalance%Nplant_demand(kidx0:kidx1,:)   * areaWeightingFactor(:,:)
         nbalance_diag%box_Nsoil_demand(kidx0:kidx1,:)    = nbalance%Nsoil_demand(kidx0:kidx1,:)    * areaWeightingFactor(:,:)
         nbalance_diag%box_Ntotal_demand(kidx0:kidx1,:)   = nbalance%Ntotal_demand(kidx0:kidx1,:)   * areaWeightingFactor(:,:)
         nbalance_diag%box_SMINN_herbivory(kidx0:kidx1,:) = nbalance%SMINN_herbivory(kidx0:kidx1,:) * areaWeightingFactor(:,:)
         nbalance_diag%box_nfix_to_sminn(kidx0:kidx1,:)   = nbalance%nfix_to_sminn(kidx0:kidx1,:)   * areaWeightingFactor(:,:)
         nbalance_diag%box_ndep_to_sminn(kidx0:kidx1,:)   = nbalance%ndep_to_sminn(kidx0:kidx1,:)   * areaWeightingFactor(:,:)
         nbalance_diag%box_nfert_to_sminn(kidx0:kidx1,:)  = nbalance%nfert_to_sminn(kidx0:kidx1,:)  * areaWeightingFactor(:,:)
         nbalance_diag%box_sminn_to_denit(kidx0:kidx1,:)  = nbalance%sminn_to_denit(kidx0:kidx1,:)  * areaWeightingFactor(:,:)
         nbalance_diag%box_sminn_leach(kidx0:kidx1,:)     = nbalance%sminn_leach(kidx0:kidx1,:)     * areaWeightingFactor(:,:)
         nbalance_diag%box_tot_sminn_loss(kidx0:kidx1,:)  = nbalance_diag%box_sminn_leach(kidx0:kidx1,:) &
                                                            + nbalance_diag%box_sminn_to_denit(kidx0:kidx1,:) &
                                                            + nbalance_diag%box_N2O_total(kidx0:kidx1,:) * 1e-06_dp &
                                                            + nbalance_diag%box_Nplant_demand(kidx0:kidx1,:) &
                                                            + nbalance_diag%box_Nsoil_demand(kidx0:kidx1,:)
         nbalance_diag%box_Nflx_2_crop_harvest(kidx0:kidx1) = SUM(nbalance%Nflx_2_crop_harvest(kidx0:kidx1,:) &
                                                             * areaWeightingFactor(:,:),  DIM=2)
         !! box_tot_sminn_gain does not yet include all gain fluxes,
         !! box_tot_sminn_gain may be better determined by summing the differences in
         !! nbalance_diag%boxN_SMINN between 2 time steps and nbalance_diag%box_tot_sminn_loss
         !! TR nbalance_diag%box_tot_sminn_gain(kidx0:kidx1,:)  = nbalance_diag%box_nfix_to_sminn(kidx0:kidx1,:)  & 
         !! TR                                                   + nbalance_diag%box_nfert_to_sminn(kidx0:kidx1,:) &
         !! TR                                                   + nbalance_diag%box_ndep_to_sminn(kidx0:kidx1,:)

         nbalance_diag%box_Npools_total(kidx0:kidx1) = SUM(sumNpools(nbalance) * areaWeightingFactor(:,:), DIM=2)

      END IF !! with nitrogen
   END IF

   ! Compute net CO2 fluxes exchanged with atmosphere at each time step
   !-----------------------------------------------------------------
   ! Note: carbon loss of biosphere means a positive CO2 flux to atmosphere (i.e. NEP and net CO2-flux have opposite signs)

   CO2_flx2atm_npp(1:nidx) = molarMassCO2_kg *                  & !! Conversion factor from mol to kg CO2
        SUM(- areaWeightingFactor(1:nidx,:)                     & !! Minus: atmosphere gain is positive
                * (cbalance%NPP_rate(kidx0:kidx1,:)             & !! current (not actual) NPP rate
                    - (cbalance_diag%NPP_yDayMean(kidx0:kidx1,:)-cbalance%NPP_act_yDayMean(kidx0:kidx1,:))) &
              , DIM=2)                                            !! corrected with yesterdays actual NPP defizit
   CO2_flx2atm_soilresp(1:nidx) = molarMassCO2_kg *             & !! Conversion factor from mol to kg CO2
        SUM(- areaWeightingFactor(1:nidx,:)                     & !! Minus: atmosphere gain is positive
                * cbalance%soil_respiration(kidx0:kidx1,:)      & !! .. soil respiration
              , DIM=2)
   CO2_flx2atm_herbivory(1:nidx) = molarMassCO2_kg *            & !! Conversion factor from mol to kg CO2
        SUM(- areaWeightingFactor(1:nidx,:)                     & !! Minus: atmosphere gain is positive
                * cbalance%Cflux_herbivory_2_atm(kidx0:kidx1,:) & !! .. herbivory
              , DIM=2)

   ! The carbon conservation test in check_Cconservation is based on yDayMean Fluxes. It cannot cope with the current
   ! NPP rate. Thus the diurnal cycle is ignored for the test. 
   CO2_flx2atm_npp_test(1:nidx) = molarMassCO2_kg *             & !! Conversion factor from mol to kg CO2
        SUM(- areaWeightingFactor(1:nidx,:)                     & !! Minus: atmosphere gain is positive
                * cbalance%NPP_act_yDayMean(kidx0:kidx1,:)      & !! yesterdays actual NPP rate
              , DIM=2)

   IF (with_nitrogen) THEN

      ! Compute net N2O release each time step
      ! positive flux to atmosphere: emissions of N2O

      N2O_flx2atm_mineraliz(1:nidx) = molarMassN2O_kg * & !! Conversion factor from mol to kg N2O
         SUM( areaWeightingFactor(1:nidx,:)             & !! N2O emissions are already positive
         *  nbalance%N2O_emissions_mineraliz(kidx0:kidx1,:), DIM=2)
      N2O_flx2atm_depfix(1:nidx) = molarMassN2O_kg *    & !! Conversion factor from mol to kg N2O
         SUM( areaWeightingFactor(1:nidx,:)             &
         *  nbalance%N2O_emissions_depfix(kidx0:kidx1,:), DIM=2)
      N2O_flx2atm_nfert(1:nidx) = molarMassN2O_kg *     & !! Conversion factor from mol to kg N2O
         SUM( areaWeightingFactor(1:nidx,:)             &
         *  nbalance%N2O_emissions_nfert(kidx0:kidx1,:), DIM=2)
      N2O_flx2atm_grazing(1:nidx) = molarMassN2O_kg *   & !! Conversion factor from mol to kg N2O
         SUM( areaWeightingFactor(1:nidx,:)             & !! N2O emissions are already positive
         *  nbalance%N2O_emissions_grazing(kidx0:kidx1,:), DIM=2)

      N2_flx2atm_ecosystem(1:nidx) = molarMassN2_kg *   & !! Conversion factor from mol to kg N2
         SUM(- areaWeightingFactor(1:nidx,:)            &
         *  nbalance%NetEcosyst_N_flux(kidx0:kidx1,:), DIM=2)
   END IF

 END SUBROUTINE update_cbalance_bethy


  ! --- NPP_rate_bethy() ----------------------------------------------------------------------------------------------------------

  ! Following BETHY the net primary productivity is estimated from dark respiration (R_d) and  gross primary production (GPP). 
  ! More precisely:  
  !
  !    (1) NPP = GPP - R_m - R_g,
  !
  ! where "R_g" is the growth respiration and "R_m" the maintenance respiration. Maintenance respiration can be estimated from
  ! dark respiration by
  !
  !                 R_d
  !    (2) R_m = ----------             (Eq. (128) in Knorr (in molar units))
  !              f_aut_leaf
  !
  ! where "f_aut_leaf" is the leaf fraction of plant-total (autotrophic) respiration. Note that "GPP" is here the Farquar
  ! productivity, where it is already accounted for photorespiration.  To estimate NPP from (1) it therefore remains to 
  ! determine the growth respiration. There is no growth respiration for NPP<0 and is otherwise a certain fraction of NPP:
  !
  !               / (cCost -1) NPP  for NPP>0
  !    (3) R_g = <
  !               \ 0               otherwise,
  !
  ! where "cCost" are the relative costs (measured in carbon) to produce 1 unit of carbon. Entering this into (1) one thus finds
  !
  !
  !                      cCost -1
  !    (4) R_g = MAX(0,  --------- (GPP -R_m)) 
  !                        cCost
  !
  ! See: W. Knorr, "Satellite Remote Sensing and Modelling of the Global CO2 Exchange of Land Vegetation: A Synthesis Study",
  !      Examensarbeit Nr. 49, (MPI for Meterology, Hamburg, 1998).
  !

  elemental pure function NPP_rate_bethy(grossAssimilation,darkRespiration)
    real(dp),intent(in) :: grossAssimilation ! Gross primary product (where photorespiration has already  ..
                                                                 ! .. been accounted for) [mol(CO2)/(m^2 s)]
    real(dp),intent(in) :: darkRespiration   ! dark respiration of leaves [mol(CO2)/(m^2 s)]
    real(dp)            :: NPP_rate_bethy    ! net primary production rate [mol(CO2)/(m^2 s)]

    ! locals

    REAL(dp), parameter :: f_aut_leaf             = 0.40_dp      ! leaf fraction of plant-total (autotrophic) respiration
    REAL(dp), parameter :: cCost                  = 1.25_dp      ! relative costs (measured in carbon) to produce 1 unit of carbon

    real(dp),parameter   :: hlp =  (cCost -1._dp)/cCost
    real(dp)             :: maintenanceRespiration ! R_m
    real(dp)             :: growthRespiration      ! R_g

    ! Go ...
    
    maintenanceRespiration = darkRespiration/f_aut_leaf
    growthRespiration = max(0._dp,hlp*(grossAssimilation - maintenanceRespiration))
    NPP_rate_bethy = grossAssimilation - maintenanceRespiration - growthRespiration

  end function NPP_rate_bethy


  !! --- read_nitrogen_deposition_fluxes() ----------------------------------------------------
  !!
  !! This routine should be called during initialization and at the beginning of the first time step of each month
  !! to read in the new nitrogen deposition fluxes [kg(N) m-2(grid-box) s-1].
  !!
  subroutine read_nitrogen_deposition_fluxes(current_year, current_month, grid, domain, nbalance)
    use mo_util_string,      only: int2string
    use mo_netcdf,           only: FILE_INFO, IO_inq_dimid, IO_inq_dimlen, IO_inq_varid, IO_get_var_double, IO_get_vara_double
    USE mo_tr_scatter,       ONLY: scatter_field
    USE mo_temp,             only: zreal2d,zzreal2d,zreal2d_ptr

    INTEGER, intent(in)                :: current_year
    INTEGER, intent(in)                :: current_month
    type(grid_type),   intent(in)      :: grid
    type(domain_type), intent(in)      :: domain
    TYPE(nbalance_type),intent(inout)  :: nbalance

    !! -- parameters

    character(len=*),parameter  :: varName="Ndepo"    !! Name of variable of the nitrogen deposition fluxes in netcdf input file

    !! -- locals

    character(len=1024)  :: filename   !! Name of netcdf input file from which nitrogen deposition shall be read
    integer              :: IO_file_id, IO_var_id, IO_dim_id
    type(FILE_INFO)      :: IO_file
    integer              :: znlon, znlat, zntime
    INTEGER              :: istart(3), icount(3)

    IF (current_year < 10000) THEN
       WRITE(filename,'(A,I4.4,A)') "Ndepo.", current_year, ".nc"
    ELSE
       CALL finish('read_nitrogen_deposition_fluxes','Only years between 0 and 9999 supported currently')
    END IF

    if (p_parallel_io) then
       ! Open ini file
       call message('read_nitrogen_deposition_fluxes()','Reading nitrogen deposition fluxes from '//trim(filename))
       IO_file%opened = .false.
       call IO_open(trim(filename), IO_file, IO_READ)
       IO_file_id = IO_file%file_id

       ! Check resolution
       call IO_inq_dimid  (IO_file_id, 'lat', IO_dim_id)
       call IO_inq_dimlen (IO_file_id, IO_dim_id, znlat)
       call IO_inq_dimid  (IO_file_id, 'lon', IO_dim_id)
       call IO_inq_dimlen (IO_file_id, IO_dim_id, znlon)
       call IO_inq_dimid  (IO_file_id, 'time', IO_dim_id)
       call IO_inq_dimlen (IO_file_id, IO_dim_id, zntime)

       if (znlon /= grid%nlon .or. znlat /= grid%nlat) then
          call finish('read_nitrogen_deposition_fluxes', 'Unexpected grid resolution:'//int2string(znlon)//'x'//int2string(znlat))
       endif


       !! read nitrogen deposition flux

       CALL IO_inq_varid(IO_file_id, TRIM(varName), IO_var_id)
       ALLOCATE(zreal2d(grid%nlon,grid%nlat))
       IF (zntime == 1) THEN
         ! annual data
         CALL IO_get_var_double(IO_file_id, IO_var_id, zreal2d)
       ELSE IF (zntime == 12) THEN
         ! monthly mean data
         istart = (/1, 1, current_month/)
         icount = (/grid%nlon, grid%nlat, 1/)
         CALL IO_get_vara_double(IO_file_id, IO_var_id, istart, icount, zreal2d)
       ELSE
         CALL finish('read_nitrogen_deposition_fluxes', &
              'Unexpected time dimension: one annual value or 12 monthly values expected')
       END IF
       CALL IO_close(IO_file)
    END IF !! end parallel_io

    !! Bring nitrogen fluxes to the other processors
    ALLOCATE(zzreal2d(domain%ndim,domain%nblocks))
    NULLIFY(zreal2d_ptr)
    IF (p_parallel_io) zreal2d_ptr => zreal2d(:,:)
    CALL scatter_field(zreal2d_ptr, zzreal2d)
    nbalance%Ndep_forc(:) = PACK(zzreal2d, MASK=domain%mask)
    IF (p_parallel_io) DEALLOCATE(zreal2d)
    DEALLOCATE(zzreal2d)

  END SUBROUTINE read_nitrogen_deposition_fluxes
  
  !! --- read_nitrogen_fertilizer_fluxes() ----------------------------------------------------
  !!
  !! This routine should be called during initialization and at the beginning of the first time step of each year
  !! to read in the new nitrogen fertilizer fluxes.
  !!
  subroutine read_nitrogen_fertilizer_fluxes(current_year, grid, domain, &
                                             nbalance, is_crop, ntiles)
    use mo_jsbach_grid,      only: grid_type,domain_type
    use mo_exception,        only: finish, message
    use mo_util_string,      only: int2string
    use mo_netcdf,           only: FILE_INFO,IO_inq_dimid,IO_inq_dimlen,IO_inq_varid,io_get_var_double
    use mo_tr_scatter,       only: scatter_field
    use mo_io,               only: IO_open, IO_READ, IO_close
    USE mo_temp,             only: zreal2d,zzreal2d,zreal2d_ptr

    INTEGER, intent(in)                :: current_year
    type(grid_type),   intent(in)      :: grid
    type(domain_type), intent(in)      :: domain
    TYPE(nbalance_type),intent(inout)  :: nbalance
    LOGICAL, intent(in)                :: is_crop(:,:)
    INTEGER, intent(in)                :: ntiles

    !! -- parameters
    character(len=*),parameter  :: varName="Nfertilizer"  !! Name of variable of the nitrogen fertilizer fluxes in netcdf input file

    !! -- locals

    character(len=1024)  :: filename   !! Name of netcdf input file from which the yearly nitrogen fertilizer fields shall be read
    integer              :: IO_file_id, IO_var_id, IO_dim_id
    type(FILE_INFO)      :: IO_file
    integer              :: znlon, znlat
    integer              :: status, itile
#if defined (__PGI)
    integer              :: i
#endif

    IF (current_year < 10000) THEN
       WRITE(filename,'(A,I4.4,A)') "Nfertilizer_", current_year, ".nc"
    ELSE
       CALL finish('read_nitrogen_fertilizer_fluxes','Only years between 0 and 9999 supported currently')
    END IF

    if (p_parallel_io) then
    
       ! Open ini file
       call message('read_nitrogen_fertilizer_fluxes()','Reading new nitrogen fertilizer fluxes from '//trim(filename))
       IO_file%opened = .false.
       call IO_open(trim(filename), IO_file, IO_READ)
       IO_file_id = IO_file%file_id

       ! Check resolution
       call IO_inq_dimid  (IO_file_id, 'lat', IO_dim_id)
       call IO_inq_dimlen (IO_file_id, IO_dim_id, znlat)
       call IO_inq_dimid  (IO_file_id, 'lon', IO_dim_id)
       call IO_inq_dimlen (IO_file_id, IO_dim_id, znlon)

       if (znlon /= grid%nlon .or. znlat /= grid%nlat) then
          call finish('read_nitrogen_fertilizer_fluxes()', 'Unexpected grid resolution:'//int2string(znlon)//'x'//int2string(znlat))
       endif

       !! allocate temporary memory

       allocate(zreal2d(grid%nlon,grid%nlat),STAT=status)
       if(status .ne. 0) call finish('read_nitrogen_fertilizer_fluxes()','Allocation failure (1)')

       !! read nitrogen fertilizer flux

       call IO_inq_varid(IO_file_id, trim(varName), IO_var_id)
       call IO_get_var_double(IO_file_id, IO_var_id, zreal2d)

       call IO_close(IO_file)
    endif !! end parallel_io

    !! Bring nitrogen fluxes to the other processors

    allocate(zzreal2d(domain%ndim,domain%nblocks),STAT=status)
    if(status .ne. 0) call finish('read_nitrogen_fertilizer_fluxes()','Allocation failure (3)')
    nullify(zreal2d_ptr)

    if (p_parallel_io) zreal2d_ptr => zreal2d(:,:)
    call scatter_field(zreal2d_ptr, zzreal2d)
    nbalance%nfert_forc(:) = pack(zzreal2d, MASK=domain%mask) 

#if defined (__PGI)

    do itile=1,ntiles
      do i = 1, size(is_crop,1)
        if (is_crop(i,itile)) then
          nbalance%nfert_forc_2d(i,itile) = nbalance%nfert_forc(i)
        else
          nbalance%nfert_forc_2d(i,itile) = 0.0_dp
        endif
      enddo
    enddo

#else

    do itile=1,ntiles
       where (is_crop(:,itile))
          nbalance%nfert_forc_2d(:,itile) = pack(zzreal2d, MASK=domain%mask)	  
       elsewhere 
          nbalance%nfert_forc_2d(:,itile) = 0.0_dp
       end where
    end do

#endif
    
    where ( nbalance%nfert_forc_2d(:,:) .lt. 0.0_dp )  
            nbalance%nfert_forc_2d(:,:) = 0.0_dp
    end where
    
!alternativ:
!    FORALL(n=1:nidx,itile=1:ntiles, is_Crop(n,itile))
!      nbalance%nfert_forc_2d(n,itile) = nbalance%nfert_forc(n)
!    end forall    

    !! Free temporary memory
    if (p_parallel_io) deallocate(zreal2d)
    deallocate(zzreal2d)

  end subroutine read_nitrogen_fertilizer_fluxes

  !----------------------------------------------------------------------------
  FUNCTION sumCnatural(cbalance, with_yasso)
  !----------------------------------------------------------------------------
    TYPE(cbalance_type), INTENT(in) :: cbalance
    LOGICAL,             INTENT(in) :: with_yasso

    REAL(dp) :: sumCnatural(kstart:kend,cbalance%ntiles)

    !
    ! sum up all natural carbon pools (i.e. C pools with tiles dimension)
    !
    
    IF (with_yasso) THEN
      sumCnatural(:,:) =   cbalance%Cpool_green(kstart:kend,:)                                                            &
                         + cbalance%Cpool_woods(kstart:kend,:)                                                            &
                         + cbalance%Cpool_reserve(kstart:kend,:)         + cbalance%Cpool_crop_harvest(kstart:kend,:)     &
                         + cbalance%YCpool_acid_ag1(kstart:kend,:)       + cbalance%YCpool_acid_ag2(kstart:kend,:)        &
                         + cbalance%YCpool_water_ag1(kstart:kend,:)      + cbalance%YCpool_water_ag2(kstart:kend,:)       &
                         + cbalance%YCpool_ethanol_ag1(kstart:kend,:)    + cbalance%YCpool_ethanol_ag2(kstart:kend,:)     &
                         + cbalance%YCpool_nonsoluble_ag1(kstart:kend,:) + cbalance%YCpool_nonsoluble_ag2(kstart:kend,:)  &
                         + cbalance%YCpool_acid_bg1(kstart:kend,:)       + cbalance%YCpool_acid_bg2(kstart:kend,:)        &
                         + cbalance%YCpool_water_bg1(kstart:kend,:)      + cbalance%YCpool_water_bg2(kstart:kend,:)       &
                         + cbalance%YCpool_ethanol_bg1(kstart:kend,:)    + cbalance%YCpool_ethanol_bg2(kstart:kend,:)     &
                         + cbalance%YCpool_nonsoluble_bg1(kstart:kend,:) + cbalance%YCpool_nonsoluble_bg2(kstart:kend,:)  &
                         + cbalance%YCpool_humus_1(kstart:kend,:)        + cbalance%YCpool_humus_2(kstart:kend,:)
    ELSE

      sumCnatural(:,:) =   cbalance%Cpool_green(kstart:kend,:)                                                            &
                         + cbalance%Cpool_woods(kstart:kend,:)                                                            &
                         + cbalance%Cpool_reserve(kstart:kend,:)         + cbalance%Cpool_crop_harvest(kstart:kend,:)     &
                         + cbalance%Cpool_litter_green_ag(kstart:kend,:) + cbalance%Cpool_litter_wood_ag(kstart:kend,:)   &
                         + cbalance%Cpool_litter_green_bg(kstart:kend,:) + cbalance%Cpool_litter_wood_bg(kstart:kend,:)   &
                         + cbalance%Cpool_slow(kstart:kend,:)
    END IF
  END FUNCTION sumCnatural

  !----------------------------------------------------------------------------
  FUNCTION sumCanthro(cbalance, with_harvest)
  !----------------------------------------------------------------------------
    TYPE(cbalance_type), INTENT(in) :: cbalance
    LOGICAL,             INTENT(in) :: with_harvest

    REAL(dp) :: sumCanthro(kstart:kend)

    sumCanthro(:) =   cbalance%Cpool_onSite(kstart:kend)         &
                    + cbalance%Cpool_paper(kstart:kend)          &
                    + cbalance%Cpool_construction(kstart:kend)

    IF (with_harvest) THEN  ! i.e. with_landuse_transitions
      sumCanthro(:) =   sumCanthro(kstart:kend)                           &
                      + cbalance%Cpool_onSite_harvest(kstart:kend)        &
                      + cbalance%Cpool_paper_harvest(kstart:kend)         &
                      + cbalance%Cpool_construction_harvest(kstart:kend)
    END IF

  END FUNCTION sumCanthro

  !----------------------------------------------------------------------------
  FUNCTION sumNpools(nbalance)
  !----------------------------------------------------------------------------
    TYPE(nbalance_type), INTENT(in) :: nbalance

    REAL(dp) :: sumNpools(kstart:kend,nbalance%ntiles)

    sumNpools(:,:) =   nbalance%Npool_green(kstart:kend,:)                                                           &
                     + nbalance%Npool_woods(kstart:kend,:)                                                           &
                     + nbalance%Npool_litter_green_ag(kstart:kend,:) + nbalance%Npool_litter_wood_ag(kstart:kend,:)  &
                     + nbalance%Npool_litter_green_bg(kstart:kend,:) + nbalance%Npool_litter_wood_bg(kstart:kend,:)  &
                     + nbalance%Npool_slow(kstart:kend,:)                                                            &
                     + nbalance%Npool_mobile(kstart:kend,:) + nbalance%Npool_crop_harvest(kstart:kend,:)             &
                     + nbalance%sminN_pool(kstart:kend,:)

  END FUNCTION sumNpools

  !----------------------------------------------------------------------------
  SUBROUTINE check_C_conservation (CO2_fluxes)
  !----------------------------------------------------------------------------

    ! Routine to check the carbon conservation within jsbach. There is a time 
    ! mismatch between the calculation of carbon fluxes and the update
    ! of the carbon pools.
    ! The current CO2-fluxes entering the routine are base on calculations at
    ! 00:00 (if new_day). They are accumulated for one day in this routine.
    ! They thus correspond to the change in the carbon pools from the day 
    ! before yesterday to yesterday. 

    USE mo_time_control,           ONLY: lstart
    USE mo_jsbach,                 ONLY: new_day

    REAL(dp), INTENT(in) :: CO2_fluxes(:) ! Sum of the carbon fluxes [kg(CO2) m-2 s-1]

    IF (.NOT. test_Cconservation) RETURN

    IF (lstart) THEN
       cbalance_diag%CO2_fluxes_accu(kstart:kend) = 0._dp
    END IF

    IF (new_day) THEN
       cbalance_diag%jsbachCconserv(kstart:kend) = cbalance_diag%box_cpools_total_old(kstart:kend) &
                                                 - cbalance_diag%box_cpools_total_old2(kstart:kend) &
                                                 + cbalance_diag%CO2_fluxes_accu(kstart:kend)
       cbalance_diag%CO2_fluxes_accu(kstart:kend) = 0._dp
       cbalance_diag%box_cpools_total_old2(kstart:kend) = cbalance_diag%box_cpools_total_old(kstart:kend)
       cbalance_diag%box_cpools_total_old(kstart:kend) = cbalance_diag%box_cpools_total(kstart:kend)
    END IF

    cbalance_diag%CO2_fluxes_accu(kstart:kend) = cbalance_diag%CO2_fluxes_accu(kstart:kend) &
            + CO2_fluxes(:) * delta_time / molarMassCO2_kg

    IF (lstart) THEN
       cbalance_diag%jsbachCconserv(kstart:kend) = 0._dp
    END IF

  END SUBROUTINE check_C_conservation

  !----------------------------------------------------------------------------
  SUBROUTINE check_N_conservation (N2_fluxes)
  !----------------------------------------------------------------------------

    ! Routine to check the nitrogen conservation within jsbach. There is a time 
    ! mismatch between the calculation of nitrogen fluxes and the update
    ! of the nitrogen pools.
    ! The current N2-fluxes entering the routine are base on calculations at
    ! 00:00 (if new_day). They are accumulated for one day in this routine.
    ! They thus correspond to the change in the carbon pools from the day 
    ! before yesterday to yesterday. 

    USE mo_time_control,           ONLY: lstart
    USE mo_jsbach,                 ONLY: new_day

    REAL(dp), INTENT(in) :: N2_fluxes(:) ! Sum of the nitrogen fluxes [kg(N2) m-2 s-1]

    IF (.NOT. test_Nconservation) RETURN

    IF (lstart) THEN
       nbalance_diag%N2_fluxes_accu(kstart:kend) = 0._dp
    END IF

    IF (new_day) THEN
       nbalance_diag%jsbachNconserv(kstart:kend) = nbalance_diag%box_npools_total_old(kstart:kend) &
                                                 - nbalance_diag%box_npools_total_old2(kstart:kend) &
                                                 + nbalance_diag%N2_fluxes_accu(kstart:kend)
       nbalance_diag%N2_fluxes_accu(kstart:kend) = 0._dp
       nbalance_diag%box_npools_total_old2(kstart:kend) = nbalance_diag%box_npools_total_old(kstart:kend)
       nbalance_diag%box_npools_total_old(kstart:kend) = nbalance_diag%box_npools_total(kstart:kend)
    END IF

    nbalance_diag%N2_fluxes_accu(kstart:kend) = nbalance_diag%N2_fluxes_accu(kstart:kend) &
            + N2_fluxes(:) * delta_time / molarMassN2_kg

    IF (lstart) THEN
       nbalance_diag%jsbachNconserv(kstart:kend) = 0._dp
    END IF

  END SUBROUTINE check_N_conservation

end module mo_cbal_bethy
