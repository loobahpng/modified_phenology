!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_jsbach_interface

  ! 
  !! Description: 
  !!   This module defines the interface between a main program (climate model
  !!   or driver program).
  !!   The interface is independent of the forcing data which are provided from
  !!   the main program. It is specific to any one land surface scheme (JSBACH
  !!   in this case) and does all the conversions from the general grid-based
  !!   forcing input to the variables (structure, units, etc.) as needed by
  !!   the land surface scheme. In particular:
  !!   - Input fields are gathered to keep just land points
  !!   - Call the routines that make up the JSBACH land surface scheme
  !!   - Scatter output fields to complete global grids
  ! 
  !! Current Code Owner: jsbach_admin
  ! 
  !! History: 
  !  
  !! Version   Date       Comment 
  !! -------   ----       ------- 
  !! 0.1       2001/06/28 Original code (based on J. Polcher's
  !!                      routine intersurf.f90). Reiner Schnur
  ! 
  Use mo_jsbach,        Only : options_type,                        &
                               debug, test_stream
  USE mo_soil,          ONLY : soil_type, soil_param_type, update_soil, soil_diagnostics, get_soil_diag
  USE mo_land_surface,  ONLY : land_surface_type, land_surface_diagnostics
  USE mo_jsbach_lctlib, ONLY : lctlib_type
  USE mo_jsbach_grid,   ONLY : grid_type, domain_type, kstart, kend, nidx, kindex
  USE mo_exception,     ONLY : finish, message
  USE mo_util_string,   ONLY : int2string
  USE mo_kind,          ONLY : dp
  USE mo_jsbach_veg,    ONLY : vegetation_type, veg_diagnostics
  USE mo_cbal_bethy,    ONLY : cbalance_type, nbalance_type, init_cbalance_bethy, update_cbalance_bethy, &
                               check_C_conservation, check_N_conservation
  USE mo_bethy,         ONLY : bethy_type, init_bethy, update_bethy, bethy_diagnostics
  USE mo_cbal_landcover_change, ONLY : landcover_change_type

  USE mo_mpi,           ONLY : p_io, p_parallel_io, p_parallel, p_bcast

  USE mo_linked_list, ONLY: t_stream
  
  USE mo_test, ONLY: init_test, init_test_2

  USE mo_time_control,    ONLY: lresume
  USE mo_climatology,     ONLY: get_from_climatology

  USE mo_dynveg,          ONLY: init_dynveg, update_dynveg
  USE mo_climbuf,         ONLY: climbuf_type, init_climbuf, update_climbuf
  USE mo_disturbance,     ONLY: init_disturbance, update_disturbance

#if defined (__SX__) && defined (_OPENMP)
  USE omp_lib,          ONLY: omp_get_thread_num, &
                              omp_get_num_threads
#endif

  IMPLICIT NONE

  PRIVATE

  ! Global (i.e. public) Declarations: 
  PUBLIC :: jsbach_interface, jsbach_inter_1d, update_current_call, stepon_jsbach, stepoff_jsbach, get_dates
  INTERFACE jsbach_interface
     MODULE PROCEDURE jsbach_inter_1d !, jsbach_inter_2d
  END INTERFACE


  TYPE land_type
     PRIVATE
     TYPE(grid_type)                :: Grid
     TYPE(domain_type)              :: Domain
     INTEGER                        :: nlct
     INTEGER                        :: npft
     INTEGER                        :: ntiles
     TYPE(land_surface_type)        :: Surface          ! State of land surface as seen by the atmosphere
     TYPE(lctlib_type)              :: LctLibrary       ! Lookup table of land cover characteristics (<nlct>)
     TYPE(soil_param_type)          :: SoilParam        ! Spatial distribution of soil characteristics
     TYPE(soil_type)                :: Soil             ! State of soil (<ntiles>)
     TYPE(vegetation_type)          :: Vegetation       ! State of vegetation (<ntiles>)
     TYPE(bethy_type)               :: Bethy            ! State of BETHY (including %Canopy for each canopy layer)
     TYPE(cbalance_type)            :: Cbalance         ! State of cbalance (NPP etc.)
     TYPE(nbalance_type)            :: Nbalance         ! State of cbalance (NPP etc.)
     TYPE(climbuf_type)             :: Climbuf          ! State of the long term climate buffer
     TYPE(landcover_change_type)    :: landcover_change ! State of landcover change and anthropogenic C-pools
     TYPE(t_stream), POINTER :: IO_stream
     TYPE(t_stream), POINTER :: IO_diag_stream
     TYPE(t_stream), POINTER :: IO_accw
     TYPE(t_stream), POINTER :: IO_veg
     TYPE(t_stream), POINTER :: IO_nitro
     TYPE(t_stream), POINTER :: IO_disturbance
     TYPE(t_stream), POINTER :: IO_yasso
  END TYPE land_type

  PUBLIC :: jsbach_init, jsbach_restart

  TYPE(options_type), SAVE       :: theOptions
  TYPE(land_type), SAVE          :: theLand
  PUBLIC :: theOptions, theLand

#ifdef STANDALONE
  EXTERNAL jsbalone_init_decomp
#endif

CONTAINS 
  
  !!+ JSBACH interface routine for 1-d call
  SUBROUTINE jsbach_inter_1d ( &
       ! Input arguments
       kdim, &                            !! Length of vectors
       kland, &                           !! Number of land points in vectors ( = SUM(lmask) )
       ! First level conditions
       wind, wind10, temp_air, qair,  &
       ! Rain, snow, radiation, surface pressure
       precip_rain, precip_snow, &
       lwdown, &
       sw_vis_net, &                      !! net solar radiation in the visible band
       sw_nir_net, &                      !! net solar radiation in the near infrared band
       sw_par_down, &                     !! downward PAR
       sw_par_frac_diffuse, &             !! fraction of diffuse solar radiation contained in sw_par_down
       pressure, &
       czenith, &
       ! CO2 concentration
       CO2_concentration, &
       ! Output: Fluxes
       evap_act, evap_pot, sensible, latent, &
       CO2_flux_npp, CO2_flux_soilresp, CO2_flux_herbivory, CO2_flux_dynveg, &
       CO2_emission_lcc, CO2_emission_harvest, &
       N2O_flux_mineraliz, N2O_flux_depfix, N2O_flux_nfert, N2O_flux_grazing, &
       ! Surface temperatures and atm. humidity
       tsoil_rad, temp_soil_new, qsurf, &
       ! Surface radiative properties
       albedo_vis, albedo_nir, emis, &
       ! roughness length and drag coef.
       z0h, z0m, cdrag, &
       ! Variables for implicit coupling
       etAcoef, etBcoef, eqAcoef, eqBcoef, cair, csat, &
       echam_zchl, &
       mask_land, land_fract, surf_dry_static_energy, &
       kblock, soil_wetness, snow_depth, &
       skin_res, veg_height, &
       tte_corr, &
       glacier_depth, &
       wsmax, &
       snow_melt &  ! for dry deposition modules
       )
    
    !! Description:

    USE mo_land_surface, ONLY: update_land_surface_fast, update_albedo, update_albedo_diag
    USE mo_land_boundary, ONLY: update_land_boundary_up
    USE mo_canopy, ONLY: unstressed_canopy_cond_par
    USE mo_jsbach_constants, ONLY: molarMassCO2_kg, molarMassN2O_kg, molarMassN2_kg, tmelt

    USE mo_time_control, ONLY : delta_time, lstart, lresume, l_trigrad, l_trigradm1
    USE mo_phenology,    ONLY : update_phenology
    USE mo_phenology_ccdas,    ONLY : update_phenology_ccdas
    USE mo_allocation,   ONLY : update_maxLai_according_to_funct_allometry
    USE mo_utils,        ONLY : average_tiles
    USE mo_cbal_landcover_change, ONLY: do_landcover_change, do_landuse_transitions
    USE mo_forest_management, ONLY: forest_management, update_forest_management
    USE mo_zenith, ONLY: cos_zenith
    USE mo_hydrology,    ONLY: hydrology_collect_land
    USE mo_test,         ONLY: write_interface_variables
    USE mo_input,        ONLY : InputUpdate

    ! Subroutine arguments
    INTEGER,            INTENT(in)    :: kdim                      !! Length of vectors (if call from ECHAM5 this
                                                                   !! will be nproma, or npromz for last block)
    INTEGER,  OPTIONAL, INTENT(in)    :: kblock                    !! Index of block  in domain that is to be
                                                                   !! processed. If missing: one call for whole domain
    INTEGER,            INTENT(in)    :: kland                     !! Number of land points in vectors
    REAL(dp),           INTENT(in)    :: wind(kdim)                !! Lowest level wind speed [m/s]
    REAL(dp),           INTENT(in)    :: wind10(kdim)              !! 10m wind speed [m/s] (for update_surface_down)
    REAL(dp),           INTENT(in)    :: temp_air(kdim)            !! Lowest level air temperature [Kelvin]
    REAL(dp),           INTENT(in)    :: qair(kdim)                !! Lowest level specific humidity
    REAL(dp),           INTENT(in)    :: precip_rain(kdim)         !! Precipitation as rain [kg/(m^2 s)]
    REAL(dp),           INTENT(in)    :: precip_snow(kdim)         !! Precipitation as snow [kg/(m^2 s)]
    REAL(dp),           INTENT(in)    :: lwdown(kdim)              !! Downward longwave flux
    REAL(dp),           INTENT(in)    :: sw_vis_net(kdim)          !! net surface visible radiation [W/m^2]
    REAL(dp),           INTENT(in)    :: sw_nir_net(kdim)          !! net surface NIR radiation [W/m^2]
    REAL(dp),           INTENT(in)    :: sw_par_down(kdim)         !! downward surface PAR [W/m^2]
    REAL(dp),           INTENT(in)    :: sw_par_frac_diffuse(kdim) !! fraction of diffuse radiation contained in sw_par_net
    REAL(dp),           INTENT(in)    :: pressure(kdim)            !! Surface pressure
    REAL(dp),           INTENT(in)    :: czenith(kdim)             !! Cosine of solar zenith angle
    REAL(dp),           INTENT(in)    :: CO2_concentration(kdim)   !! Atmospheric CO2 concentration [kg(CO2)/kg(air)]
    REAL(dp),           INTENT(in)    :: etAcoef(kdim)
    REAL(dp),           INTENT(in)    :: etBcoef(kdim)
    REAL(dp),           INTENT(in)    :: eqAcoef(kdim)
    REAL(dp),           INTENT(in)    :: eqBcoef(kdim)
    REAL(dp),           INTENT(in)    :: cdrag(kdim)               !! Surface drag
    REAL(dp),           INTENT(in)    :: echam_zchl(kdim)
    LOGICAL,  OPTIONAL, INTENT(in)    :: mask_land(kdim)           !! Land-sea mask (land includes glaciers)
    REAL(dp), OPTIONAL, INTENT(in)    :: land_fract(kdim)          !! Land fraction
    REAL(dp), OPTIONAL, INTENT(out)   :: evap_act(kdim)            !! Total of evaporation (actual) [kg/(m^2 s)]
    REAL(dp), OPTIONAL, INTENT(out)   :: evap_pot(kdim)            !! Potential evaporation [kg/(m^2 s)]
    REAL(dp), OPTIONAL, INTENT(out)   :: sensible(kdim)            !! Sensible heat flux [W/m^2]
    REAL(dp), OPTIONAL, INTENT(out)   :: latent(kdim)              !! Latent heat flux [W/m^2]
    REAL(dp), OPTIONAL, INTENT(out)   :: CO2_flux_npp(kdim)        !! CO2 flux due to NPP [kg(CO2)/m^2 s]
    REAL(dp), OPTIONAL, INTENT(out)   :: CO2_flux_soilresp(kdim)   !! CO2 flux due to soil respiration [kg(CO2)/m^2 s]
    REAL(dp), OPTIONAL, INTENT(out)   :: CO2_flux_herbivory(kdim)  !! CO2 flux due to herbivory [kg(CO2)/m^2 s]
    REAL(dp), OPTIONAL, INTENT(out)   :: CO2_flux_dynveg(kdim)     !! CO2 flux due to (natural) fires (dynveg) [kg(CO2)/m^2 s]
    REAL(dp), OPTIONAL, INTENT(out)   :: CO2_emission_lcc(kdim)    !! CO2 emission due to landcover change [kg(CO2)/m^2 s]
    REAL(dp), OPTIONAL, INTENT(out)   :: CO2_emission_harvest(kdim)!! CO2 emission due to harvest [kg(CO2)/m^2 s]
    REAL(dp), OPTIONAL, INTENT(out)   :: N2O_flux_mineraliz(kdim)  !! N2O emission from mineralisation [kg(N2O)/m^2 s]
    REAL(dp), OPTIONAL, INTENT(out)   :: N2O_flux_depfix(kdim)     !! N2O emission from deposition and fixation [kg(N2O)/m^2 s]
    REAL(dp), OPTIONAL, INTENT(out)   :: N2O_flux_nfert(kdim)      !! N2O emission from fertilization [kg(N2O)/m^2 s]
    REAL(dp), OPTIONAL, INTENT(out)   :: N2O_flux_grazing(kdim)    !! N2O emission from herbivores [kg(N2O)/m^2 s]
    REAL(dp), OPTIONAL, INTENT(out)   :: tsoil_rad(kdim)           !!
    REAL(dp), OPTIONAL, INTENT(out)   :: temp_soil_new(kdim)       !! New soil temperature
    REAL(dp), OPTIONAL, INTENT(out)   :: qsurf(kdim)               !! Surface specific humidity
    REAL(dp), OPTIONAL, INTENT(out)   :: albedo_vis(kdim)          !! Albedo of the visible range
    REAL(dp), OPTIONAL, INTENT(out)   :: albedo_nir(kdim)          !! Albedo of the NIR range
    REAL(dp), OPTIONAL, INTENT(out)   :: emis(kdim)                !! Emissivity
    REAL(dp), OPTIONAL, INTENT(out)   :: z0h(kdim)                 !! Surface roughness for heat
    REAL(dp), OPTIONAL, INTENT(out)   :: z0m(kdim)                 !! Surface roughness for momentum
    REAL(dp), OPTIONAL, INTENT(out)   :: skin_res(kdim)            !! Skin reservoir for g3b
    REAL(dp), OPTIONAL, INTENT(out)   :: veg_height(kdim)          !! average vegetation height above surface
    REAL(dp), OPTIONAL, INTENT(out)   :: tte_corr(kdim)            !! correction of tte without zdp (!)
    REAL(dp), OPTIONAL, INTENT(out)   :: snow_depth(kdim)          !! snow on soil acc. to sn in g3b
    REAL(dp), OPTIONAL, INTENT(out)   :: soil_wetness(kdim)        !! Soil moisture acc. to pws in echam for diagn - g3b
    REAL(dp), OPTIONAL, INTENT(out)   :: cair(kdim), csat(kdim)    !!
    REAL(dp), OPTIONAL, INTENT(out)   :: surf_dry_static_energy(kdim) !! qslnew in vdiff (zqsl*radiative_temp)
    ! The following variables are passed back to echam just to be put into the standard echam output streams.
    ! This should be revised once the echam and jsbach output is reconciled.
    REAL(dp), OPTIONAL, INTENT(out)   :: glacier_depth(kdim)
    REAL(dp), OPTIONAL, INTENT(out)   :: wsmax(kdim)
    REAL(dp), OPTIONAL, INTENT(out)   :: snow_melt(kdim)              ! for dry deposition modules
      
    !! Local declarations for packed (gathered) input fields
    REAL(dp), DIMENSION(kland) ::   zwind, zwind10, ztemp_air, zqair
    REAL(dp), DIMENSION(kland) ::   zetAcoef, zetBcoef, zeqAcoef, zeqBcoef
    REAL(dp), DIMENSION(kland) ::   zcdrag
    REAL(dp), DIMENSION(kland) ::   zprecip_rain, zprecip_snow         !! Packed precipitation
    REAL(dp), DIMENSION(kland) ::   zprecip_tot                        !! Packed total precipitation
    REAL(dp), DIMENSION(kland) ::   zCO2                               !! Packed CO2 concentration
    REAL(dp), DIMENSION(kland) ::   zlwdown                            !! Packed longwave radiation
    REAL(dp), DIMENSION(kland) ::   zsw_vis_net                        !! Packed radiation from visible band [W/m^2]
    REAL(dp), DIMENSION(kland) ::   zsw_nir_net                        !! Packed radiation from NIR band [W/m^2]
    REAL(dp), DIMENSION(kland) ::   zsw_par_down                       !! Packed surface downward PAR [W/m^2]
    REAL(dp), DIMENSION(kland) ::   zsw_par_frac_diffuse               !! Packed fraction of diffuse PAR
    REAL(dp), DIMENSION(kland) ::   zpressure                          !! Packed surface pressure
    REAL(dp), DIMENSION(kland) ::   zevap_act, zevap_pot               !! Packed actual evapotranspiration and potential evaporation
    REAL(dp), DIMENSION(kland) ::   zsensible, zlatent                 !! Packed surface fluxes
    REAL(dp), DIMENSION(kland) ::   zCO2_flux_npp                      !! CO2 flux from NPP [kg(CO2)/m^2 s]
    REAL(dp), DIMENSION(kland) ::   zCO2_flux_npp_test                 !! CO2 flux from NPP, only used for conservation test
    REAL(dp), DIMENSION(kland) ::   zCO2_flux_soilresp                 !! CO2 flux from soil respiration [kg(CO2)/m^2 s]
    REAL(dp), DIMENSION(kland) ::   zCO2_flux_herbivory                !! CO2 flux from herbivory [kg(CO2)/m^2 s]
    REAL(dp), DIMENSION(kland) ::   zCO2_flux_dynveg                   !! CO2 flux due to fires (dynveg)
    REAL(dp), DIMENSION(kland) ::   zCO2_emission_landcover_change     !! CO2 emission from landcover change
    REAL(dp), DIMENSION(kland) ::   zCO2_emission_harvest              !! CO2 emission from harvest
    REAL(dp), DIMENSION(kland) ::   zCO2_emission_forest_management    !! CO2 emission from forest management
    REAL(dp), DIMENSION(kland) ::   zN2O_flux_mineraliz                !! N2O emission from mineralisation [kg(N2O)/m^2 s]
    REAL(dp), DIMENSION(kland) ::   zN2O_flux_depfix                   !! N2O emission from deposition and fixation [kg(N2O)/m^2 s]
    REAL(dp), DIMENSION(kland) ::   zN2O_flux_nfert                    !! N2O emission from fertilization [kg(N2O)/m^2 s]
    REAL(dp), DIMENSION(kland) ::   zN2O_flux_grazing                  !! N2O emission from herbivores [kg(N2O)/m^2 s]
    REAL(dp), DIMENSION(kland) ::   zN2_flux_dynveg                    !! N2 emission from harvest
    REAL(dp), DIMENSION(kland) ::   zN2_emission_ecosystem             !! net N2 emission from ecosystem
    REAL(dp), DIMENSION(kland) ::   zN2_emission_landcover_change      !! N2 emission from landcover change
    REAL(dp), DIMENSION(kland) ::   zN2_emission_harvest               !! N2 emission from harvest
    REAL(dp), DIMENSION(kland) ::   ztsoil_rad, ztemp_soil_new, zqsurf
    REAL(dp), DIMENSION(kland) ::   zalbedo_vis                        !! Packed albedo of visible range
    REAL(dp), DIMENSION(kland) ::   zalbedo_nir                        !! Packed albedo of NIR range
    REAL(dp), DIMENSION(kland) ::   zalbedo_vis_nosno                  !! Packed albedo of visible range without snow
    REAL(dp), DIMENSION(kland) ::   zalbedo_nir_nosno                  !! Packed albedo of NIR range without snow
    REAL(dp), DIMENSION(kland) ::   zemis                              !! Packed emissivity
    REAL(dp), DIMENSION(kland) ::   zz0h                               !! Packed surface roughness for heat
    REAL(dp), DIMENSION(kland) ::   zz0m                               !! Packed surface roughness for momentum
    REAL(dp), DIMENSION(kland) ::   zsnow_depth                        !! Packed soil snow
    REAL(dp), DIMENSION(kland) ::   zsoil_wetness                      !! PAcked soil moisture (acc. to pws in old echam)
    REAL(dp), DIMENSION(kland) ::   zcair, zcsat
    REAL(dp), DIMENSION(kland) ::   p_echam_zchl
    REAL(dp), DIMENSION(kland) ::   zskin_res
    REAL(dp), DIMENSION(kland) ::   zveg_height

    !! Other local declarations
    LOGICAL,  DIMENSION(kdim)  ::   mask                               !! Land mask or all true's
    REAL(dp), DIMENSION(kland) ::   swnet                              !! Total net shortwave radiation [W/m^2]
    REAL(dp), DIMENSION(kland) ::   sw_par_fract_direct                !! fraction of direct radiation in PAR-band
    REAL(dp), DIMENSION(kland) ::   ztte_corr                          !! correction of tte multiplied by zqp (!)
    REAL(dp), DIMENSION(kland) ::   zdry_static_energy                 !! new surface dry static energy
    INTEGER,  DIMENSION(kland,theLand%ntiles)   :: pheno_type          !! Mapping of land cover types to phenology types
    INTEGER,  DIMENSION(kland,theLand%ntiles)   :: pft_type
    REAL(dp), DIMENSION(kland,theLand%ntiles)   :: canopy_conductance_unstressed
    REAL(dp), DIMENSION(kland,theLand%ntiles)   :: specificLeafArea_C 
    REAL(dp), DIMENSION(kland,theLand%ntiles)   :: cover_fract_old     !! cover fraction of previous timestep
    REAL(dp), DIMENSION(kland,theLand%ntiles,5) :: LeafLit_coef        !! factor to seperate non woody litterfall into yasso pools
    REAL(dp), DIMENSION(kland,theLand%ntiles,5) :: WoodLit_coef        !! factor to seperate non woody litterfall into yasso pools

    REAL(dp), DIMENSION(kland) :: zglac_runoff_evap
    REAL(dp), DIMENSION(kland) :: zsurf_runoff_hd  !! Surface runoff in [m] for whole time step
    REAL(dp), DIMENSION(kland) :: zdrainage_hd     !! Drainage in [m] for whole time step

    REAL(dp), DIMENSION(kland) :: relative_wilting_point   !! HW
    REAL(dp), DIMENSION(kland,theLand%ntiles) :: relative_extractable_water   !! HW REW
    REAL(dp), DIMENSION(kland,theLand%ntiles) :: leaf_shedding_rate   !! HW

    LOGICAL :: climbuf_init_running_means          !! climbuf option that is also needed in dynveg
    INTEGER :: ntiles, nsoil

    INTEGER                 ::   i, kidx0, kidx1, itile

#if defined (__SX__) && defined (_OPENMP)
    INTEGER :: tid          ! OpenMP thread number
#endif

#if defined (__SX__) && defined (_OPENMP)
    IF (debug) THEN
       tid = omp_get_thread_num()
       CALL message('jsbach_interface', 'OpenMP thread #'//int2string(tid))
       CALL message('                ', '   kdim    = '//int2string(kdim))
       CALL message('                ', '   kland   = '//int2string(kland))
       CALL message('                ', '   kblock  = '//int2string(kblock))
    END IF
#endif

    IF (PRESENT(mask_land)) THEN
       IF (COUNT(mask_land(1:kdim)) /= kland) THEN
          CALL message('jsbach_interface ','kland, COUNT(mask_land) = '// &
                            int2string(kland)//', '//int2string(COUNT(mask_land(1:kdim))))
          CALL finish('jsbach_interface ','wrong mask_land')
       END IF
    END IF

    IF (.NOT. PRESENT(kblock)) THEN
      CALL InputUpdate()
    ELSE
      IF (kblock == 1) CALL InputUpdate()
    ENDIF

    IF (kland == 0) THEN
       IF (PRESENT(evap_act)) evap_act                             = 0.0_dp
       IF (PRESENT(evap_pot)) evap_pot                             = 0.0_dp
       IF (PRESENT(sensible)) sensible                             = 0.0_dp
       IF (PRESENT(latent)) latent                                 = 0.0_dp
       IF (PRESENT(CO2_flux_npp)) CO2_flux_npp                     = 0.0_dp
       IF (PRESENT(CO2_flux_soilresp)) CO2_flux_soilresp           = 0.0_dp
       IF (PRESENT(CO2_flux_herbivory)) CO2_flux_herbivory         = 0.0_dp
       IF (PRESENT(CO2_flux_dynveg)) CO2_flux_dynveg               = 0.0_dp
       IF (PRESENT(CO2_emission_lcc)) CO2_emission_lcc             = 0.0_dp
       IF (PRESENT(CO2_emission_harvest)) CO2_emission_harvest     = 0.0_dp
       IF (PRESENT(N2O_flux_mineraliz)) N2O_flux_mineraliz         = 0.0_dp
       IF (PRESENT(N2O_flux_depfix)) N2O_flux_depfix               = 0.0_dp
       IF (PRESENT(N2O_flux_nfert)) N2O_flux_nfert                 = 0.0_dp
       IF (PRESENT(N2O_flux_grazing))  N2O_flux_grazing            = 0.0_dp
!!$       IF (PRESENT(tsoil_rad)) tsoil_rad                           = 280.0_dp
!!$       IF (PRESENT(temp_soil_new)) temp_soil_new                   = 280.0_dp
       IF (PRESENT(tsoil_rad)) tsoil_rad                           = 0.0_dp
       IF (PRESENT(temp_soil_new)) temp_soil_new                   = 0.0_dp
       IF (PRESENT(qsurf)) qsurf                                   = 0.0_dp
       IF (PRESENT(albedo_vis)) albedo_vis                         = 0.0_dp
       IF (PRESENT(albedo_nir)) albedo_nir                         = 0.0_dp
       IF (PRESENT(emis)) emis                                     = 0.0_dp
       IF (PRESENT(z0h)) z0h                                       = 0.0_dp
       IF (PRESENT(z0m)) z0m                                       = 0.0_dp
       IF (PRESENT(cair)) cair                                     = 0.0_dp
       IF (PRESENT(csat)) csat                                     = 0.0_dp
       IF (PRESENT(surf_dry_static_energy)) surf_dry_static_energy = 2.276e+05_dp
       IF (PRESENT(soil_wetness)) soil_wetness                     = 0.0_dp
       IF (PRESENT(snow_depth)) snow_depth                         = 0.0_dp
       IF (PRESENT(skin_res)) skin_res                             = 0.0_dp
       IF (PRESENT(veg_height)) veg_height                         = 0.0_dp
       IF (PRESENT(tte_corr)) tte_corr                             = 0.0_dp
       IF (PRESENT(glacier_depth)) glacier_depth                   = 0.0_dp
       IF (PRESENT(wsmax)) wsmax                                   = 0.0_dp
       IF (PRESENT(snow_melt)) snow_melt                           = 0.0_dp

       RETURN
    END IF

    CALL update_current_call(kdim, kland, kblock=kblock, mask=mask_land)

    IF (theOptions%WriteInterfaceVars) THEN
       CALL write_interface_variables(theLand%grid, theLand%domain, wind, wind10, temp_air, qair, &
            etAcoef, etBcoef, eqAcoef, eqBcoef, cdrag, precip_rain, precip_snow, CO2_concentration, lwdown, &
            sw_vis_net, sw_nir_net, sw_par_frac_diffuse, pressure, echam_zchl, &
            sw_par_down, czenith)
    END IF

    kidx0 = kstart
    kidx1 = kend

#if defined (__SX__) && defined (_OPENMP)
    IF (debug) THEN
       tid = omp_get_thread_num()
       CALL message('jsbach_interface', 'OpenMP thread #'//int2string(tid)//' '//int2string(kblock)//' '//int2string(kland)//&
            ' '//int2string(nidx)//' '//int2string(kidx0)//' '//int2string(kidx1))
    END IF
#endif

    IF (PRESENT(mask_land)) THEN
       mask = mask_land
    ELSE
       mask = .TRUE.
    ENDIF

    IF (PRESENT(land_fract)) THEN
      theLand%Surface%land_fract(kidx0:kidx1) = PACK(land_fract, MASK=mask)
    END IF

    ntiles = theLand%ntiles
    nsoil = theLand%Soil%nsoil

    !Generate local packed forcing arrays for each domain (processor)
    zwind        = PACK(wind, MASK=mask)
    zwind10      = PACK(wind10, MASK=mask)
    ztemp_air    = PACK(temp_air, MASK=mask)
    zqair        = PACK(qair, MASK=mask)
    zcdrag       = PACK(cdrag, MASK=mask)
    zetAcoef     = PACK(etAcoef, MASK=mask)
    zetBcoef     = PACK(etBcoef, MASK=mask)
    zeqAcoef     = PACK(eqAcoef, MASK=mask)
    zeqBcoef     = PACK(eqBcoef, MASK=mask)
    zprecip_rain = PACK(precip_rain, MASK=mask)
    zprecip_snow = PACK(precip_snow, MASK=mask)
    zlwdown      = PACK(lwdown, MASK=mask)
    zsw_vis_net  = PACK(sw_vis_net, MASK=mask)
    zsw_nir_net  = PACK(sw_nir_net, MASK=mask)
    cos_zenith(kidx0:kidx1) = PACK(czenith, MASK=mask)
    zsw_par_down = PACK(sw_par_down, MASK=mask)
    zsw_par_frac_diffuse = PACK(sw_par_frac_diffuse, MASK=mask)
    zpressure    = PACK(pressure, MASK=mask)
    zCO2         = PACK(CO2_concentration, MASK=mask)
    p_echam_zchl = PACK(echam_zchl, MASK=mask)

    zprecip_tot = zprecip_rain + zprecip_snow

    !! Fill global arrays from lctlib information
    !! TR: this should be moved to the initialization routine in mo_jsbach_veg.f90 and mo_cbal_bethy.f90
    DO itile=1,ntiles
       theLand%Soil%relative_moisture(kidx0:kidx1,itile) = &
            theLand%Soil%moisture(kidx0:kidx1,itile)/theLand%soilparam%MaxMoisture(kidx0:kidx1)
       specificLeafArea_C(1:nidx,itile) = &
            theLand%LctLibrary%specificLeafArea_C(theLand%Surface%cover_type(kidx0:kidx1,itile))
       LeafLit_coef(1:nidx,itile,:) = &
            theLand%LctLibrary%LeafLit_coef(theLand%Surface%cover_type(kidx0:kidx1,itile),:)
       WoodLit_coef(1:nidx,itile,:) = &
            theLand%LctLibrary%WoodLit_coef(theLand%Surface%cover_type(kidx0:kidx1,itile),:)
       theLand%Vegetation%veg_height(kidx0:kidx1,itile) = &
            theLand%LctLibrary%VegHeight(theLand%Surface%cover_type(kidx0:kidx1,itile))

    ! proof that the values in lctlib are consistent with YASSO
       DO i = 1,nidx 
          IF ( ABS(SUM(LeafLit_coef(i,itile,:)) - 1.0_dp) > EPSILON(1.0_dp) ) THEN
             WRITE(*,*) "Leaf litter coefficients have invalid values; sum must be 1, but is:" 
             WRITE(*,*) SUM(LeafLit_coef(i,itile,:)), "for tile=", itile
             CALL finish('jsbach_interfacei_1d','invalid sum of leaf litter coefficients')
          END IF        

          IF ( ABS(SUM(WoodLit_coef(i,itile,:)) - 1.0_dp) > EPSILON(1.0_dp) ) THEN
             WRITE(*,*) "Wood litter coefficients have invalid values; sum must be 1, but is:" 
             WRITE(*,*) SUM(WoodLit_coef(i,itile,:)), "for tile=", itile
             CALL finish('jsbach_interfacei_1d','invalid sum of wood litter coefficients')
          END IF        
          IF (theOptions%Withnitrogen) THEN
             IF ( ABS(LeafLit_coef(i,itile,5)) > EPSILON(1.0_dp) ) THEN
                WRITE(*,*) "Leaf litter coefficient has invalid value for use with nitrogen"
                WRITE(*,*) "Leaf litter coefficient must be 0 for humus, but is:" 
                WRITE(*,*) LeafLit_coef(i,itile,5), "for tile=", itile
                CALL finish('jsbach_interfacei_1d','non zero leaf litter coefficient for humus')
             ENDIF
             IF ( ABS(WoodLit_coef(i,itile,5)) > EPSILON(1.0_dp) ) THEN
                WRITE(*,*) "Wood litter coefficient has invalid value for use with nitrogen"
                WRITE(*,*) "Wood litter coefficient must be 0 for humus, but is:" 
                WRITE(*,*) WoodLit_coef(i,itile,5), "for tile=", itile
                CALL finish('jsbach_interfacei_1d','non zero wood litter coefficient for humus') 
             ENDIF
          ENDIF
       END DO
    END DO


    cover_fract_old(1:nidx,1:theLand%ntiles) = theLand%Surface%cover_fract(kidx0:kidx1,1:theLand%ntiles)

    !! Recompute natural vegetation cover

    IF (theOptions%UseDynveg .OR. theOptions%UseDisturbance) THEN
       CALL update_climbuf(nidx, kidx0, kidx1, ntiles, nsoil, &
            theLand%LctLibrary%gdd_base(:), theLand%LctLibrary%upper_tlim(1:ntiles), &
            ztemp_air - tmelt, &
            theLand%Soil%surface_temperature(kidx0:kidx1,:) - tmelt, &
            zwind10, &
            zprecip_rain(:) + zprecip_snow(:), &
            theLand%Soil%relative_moisture(kidx0:kidx1,1:ntiles), &
            theLand%cbalance%NPP_rate(kidx0:kidx1,1:ntiles), &
            theLand%Soil%relative_humidity_air(kidx0:kidx1), &
            Climbuf_Init_running_means,theLand%Climbuf)

       IF (theOptions%UseDynveg) THEN
          IF (debug) CALL message('jsbach_interface_1d','Calling update_dynveg')
          CALL update_dynveg (lstart, lresume, nidx, kidx0, kidx1, ntiles, &
               theOptions%UseDisturbance, theOptions%ReadCoverFract, theOptions%WithYasso, theOptions%WithNitrogen, &
               theOptions%LCC_scheme==2, theOptions%UseLanduseTransitions, &
               theLand%LctLibrary, theLand%Surface, theLand%Cbalance, theLand%Nbalance, theLand%Climbuf, &
               SpecificLeafArea_C(1:nidx,1:ntiles), &
               theLand%Vegetation%veg_fract_correction(kidx0:kidx1,1:ntiles),  &
               zCO2_flux_dynveg(1:nidx), zN2_flux_dynveg(1:nidx), &
               LeafLit_coef(1:nidx,1:ntiles,:), WoodLit_coef(1:nidx,1:ntiles,:))
          IF (debug) CALL message('jsbach_interface_1d','Returned from update_dynveg')
       ELSE
          CALL update_disturbance(nidx, kidx0, kidx1, theLand%climbuf, theLand%Surface, theLand%Cbalance, &
               theOptions%WithYasso, theOptions%WithNitrogen, theLand%Nbalance, &
               theLand%LctLibrary, theLand%Vegetation%veg_fract_correction(kidx0:kidx1,1:ntiles), &
               zCO2_flux_dynveg(1:nidx), zN2_flux_dynveg(1:nidx), &
               LeafLit_coef(1:nidx,1:ntiles,:), WoodLit_coef(1:nidx,1:ntiles,:))
       ENDIF

    ELSE
       zCO2_flux_dynveg = 0._dp
       zN2_flux_dynveg = 0._dp
    ENDIF

    !! --- Landcover change ---
    !! In case landcover change is driven by a sequence of external maps or landuse transitions, read in the appropriate 
    !! maps/transitions, update the landcover map (theLand%Surface%cover_fract), perform the necessary changes in the carbon 
    !! pools and compute the CO2 emissions from the change of landcover. 
    !! NOTE: landcover change has to be done before any other submodels of JSBACH are called to guarantee that during the
    !!       whole timestep the same landcover map is used! 
    !! NOTE: If the model is coupled to mpiom/hamocc, the prism setup starts the run one time step before the actual initial
    !!       date, usually the last time step in December of the previous year. Therefore, although do_landcover_change is
    !!       called, no land cover data are read (since current_date and previous_date belong to the same year). This is ok,
    !!       since for this first time step, the cover fractions from the ini file (in %cover_fract) can be used. However, it
    !!       seems to be important that the cover_fract in the ini file is the same as the one for the first year ... 
    !!       otherwise, the model crashes with negative CO2 concentrations, probably because the change in land cover is
    !!       too sudden, leading to strange effects in the CO2 transport.

    zCO2_emission_landcover_change(:) = 0._dp
    zCO2_emission_harvest(:) = 0._dp
    zN2_emission_landcover_change(:) = 0._dp
    zN2_emission_harvest(:) = 0._dp
    zN2_emission_ecosystem(:) = 0._dp

    IF (theOptions%UseExternalLandcoverMaps) THEN

       CALL do_landcover_change(theLand%ntiles,                                                        &
                                theLand%LctLibrary,                                                    &
                                theLand%Surface,                                                       &
                                theLand%Cbalance,                                                      & !! <-- inout
                                theLand%Nbalance,                                                      & !! <-- inout
                                theOptions%WithNitrogen,                                               &
                                theOptions%WithYasso,                                                  &
                                LeafLit_coef(1:nidx,1:theLand%ntiles,1:5),                             &
                                WoodLit_coef(1:nidx,1:theLand%ntiles,1:5),                             &
                                theLand%Surface%veg_ratio_max(kidx0:kidx1),                            &
                                theLand%Vegetation%veg_fract_correction(kidx0:kidx1,1:theLand%ntiles), &
                                theLand%Surface%cover_fract(kidx0:kidx1,1:theLand%ntiles),             & !! <-- inout
                                zCO2_emission_landcover_change(1:nidx),                                & !! <-- out
                                zN2_emission_landcover_change(1:nidx),                                 & !! <-- out
                                theLand%landcover_change,                                              & !! <-- inout
                                theOptions%lcc_scheme)
    ELSE IF (theOptions%UseLanduseTransitions) THEN

       CALL do_landuse_transitions(theLand%LctLibrary,theLand%Surface,theLand%ntiles,             &
                                   theOptions%UseNetLanduseTransitions,                           &
                                   theLand%Surface%is_naturalVeg(kidx0:kidx1,1:theLand%ntiles),   &
                                   theLand%Surface%is_C4vegetation(kidx0:kidx1,1:theLand%ntiles), &
                                   theLand%Surface%is_grass(kidx0:kidx1,1:theLand%ntiles),    &
                                   theLand%Surface%is_crop(kidx0:kidx1,1:theLand%ntiles),     &
                                   theLand%Surface%is_pasture(kidx0:kidx1,1:theLand%ntiles),  &
                                   theLand%Surface%is_vegetation(kidx0:kidx1,:),              &
                                   theLand%Surface%is_glacier(kidx0:kidx1,:),                 &
                                   theLand%Surface%veg_ratio_max(kidx0:kidx1),                            &
                                   theLand%Vegetation%veg_fract_correction(kidx0:kidx1,1:theLand%ntiles), &
                                   theLand%Surface%cover_fract_pot(kidx0:kidx1,1:theLand%ntiles), &
                                   theLand%Surface%cover_fract(kidx0:kidx1,1:theLand%ntiles), & !! <-- inout
                                   theLand%Cbalance,                                          & !! <-- inout
                                   zCO2_emission_landcover_change(1:nidx),                    & !! <-- out
                                   zCO2_emission_harvest(1:nidx),                             & !! <-- out
                                   zN2_emission_landcover_change(1:nidx),                     & !! <-- out
                                   zN2_emission_harvest(1:nidx),                              & !! <-- out
                                   theLand%landcover_change,                                  &
                                   LeafLit_coef(1:nidx,1:theLand%ntiles,1:5),                 &
                                   WoodLit_coef(1:nidx,1:theLand%ntiles,1:5),                 &
                                   theOptions%WithYasso,                                      &
                                   theLand%Nbalance,                                          & !! <-- inout
                                   theOptions%WithNitrogen,                                   &
                                   theOptions%lcc_scheme)

       IF (debug) CALL message('jsbach_interface_1d','Returned from do_landuse_transitions()')
    END IF

    IF (theOptions%UseForestManagement) THEN
       IF (debug) CALL message('jsbach_interface_1d','Calling update_forest_management')
       CALL update_forest_management (nidx, kidx0, kidx1, ntiles, theOptions%WithYasso, LeafLit_coef, WoodLit_coef, &
         & theOptions%UsePFTDiameterForFM, theOptions%UseDiameterCWoodDistr, theOptions%tau_fm_burning, &
         & theOptions%tau_fm_paper, theOptions%tau_fm_construction, theLand%Surface, theLand%LctLibrary, &
         & theLand%Vegetation%veg_fract_correction, theLand%Cbalance, zCO2_emission_forest_management(1:nidx),&
         & theOptions%UseRotationMap)

       IF (debug) CALL message('jsbach_interface_1d','Returned from update_forest_management')

       !Since we anyway currently only simulate forest management in cases without lcc, emissions from forest management are put
       !into zCO2_emission_harvest here, such that the carbon cycle can be closed but mo_surface etc do no need to be adapted
       !-> ASSERTION: this has the precondition, that lcc is really switched off, so assure this here again:
       IF(theOptions%UseLanduseTransitions .OR. theOptions%UseExternalLandcoverMaps) CALL finish('jsbach_inter_1d', &
         & 'FM not implemented to work with landuse change!')
       zCO2_emission_harvest = zCO2_emission_forest_management
    ENDIF

    IF (theOptions%UsePhenology) THEN
       DO itile=1,ntiles
          WHERE (theLand%Surface%is_vegetation(kidx0:kidx1,itile))
             pheno_type(1:nidx,itile) = theLand%LctLibrary%PhenologyType(theLand%Surface%cover_type(kidx0:kidx1,itile))
             pft_type(1:nidx,itile) = theLand%Surface%cover_type(kidx0:kidx1,itile)
          ELSEWHERE
             pheno_type(1:nidx,itile) = 0
             pft_type(1:nidx,itile) = 0
          END WHERE
       END DO
       canopy_conductance_unstressed=0._dp
    END IF


    ! If demanded, update max lai according to allometric relationships
    ! -> LAI_max, veg_fract_correction and the c densities might be adapted
    IF (theOptions%UseFunctAllometry) THEN

          IF (theOptions%UseDisturbance) THEN
             CALL finish('jsbach_config','using allometric relationships not tested (and without cohorts not logical) with disturbance')
          END IF
          IF (theOptions%UseDynveg) THEN
             CALL finish('jsbach_config','using allometric relationships not tested (and without cohorts not logical) with dynveg')
          END IF
          IF (theOptions%UseLanduseTransitions .OR. theOptions%UseExternalLandcoverMaps) THEN
             CALL finish('jsbach_config','using allometric relationships not tested (and without cohorts not logical) with landuse')
          END IF


        !ALLOCATE (stopCO2_flag(kland,theLand%ntiles))
        IF (debug) CALL message('jsbach_interface_1d','Entering allocation') !KN
        CALL update_maxLai_according_to_funct_allometry(nidx, kidx0, kidx1, ntiles, theLand%Surface, theLand%Cbalance,  &
                                                        theLand%LctLibrary,                                             &
                                                        theLand%Vegetation%LAI_max(kidx0:kidx1,1:ntiles),               &
                                                        theLand%Vegetation%veg_fract_correction(kidx0:kidx1,1:ntiles),  &
                                                        theOptions%UseLimitCO2Fert, theOptions%WithNitrogen, theOptions%WithYasso)
        IF (debug) CALL message('jsbach_interface_1d','Returning from allocation') !KN

        IF (.NOT. lstart) THEN
            ! Update of lai and box_veg_ratio (in case of lstart they are set below)
            DO i = 1,ntiles
                theLand%Vegetation%lai(kidx0:kidx1,i) = &
                     MIN(theLand%Vegetation%lai(kidx0:kidx1,i),theLand%Vegetation%lai_max(kidx0:kidx1,i))
                theLand%Surface%veg_ratio(kidx0:kidx1,i) = theLand%Surface%veg_ratio_max(kidx0:kidx1) &
                    * (1.0_dp - exp(-theLand%Vegetation%lai(kidx0:kidx1,i) &
                               / theLand%LctLibrary%clumpinessFactor(theLand%Surface%cover_type(kidx0:kidx1,i))))
            END DO
            theLand%Surface%box_veg_ratio(kidx0:kidx1,1:ntiles) = &
                theLand%Surface%veg_ratio(kidx0:kidx1,1:ntiles) * theLand%Surface%cover_fract(kidx0:kidx1,1:ntiles)
        ENDIF
    ENDIF

    ! If first time step initialize lai with climatological value

    IF (lstart) THEN
       IF (debug) CALL message('jsbach_interface','Get LAI from climatology at model start')

       CALL get_from_climatology(nidx, ntiles, theLand%Vegetation%lai_clim(kidx0:kidx1,1:ntiles,0:13), &
            theLand%Vegetation%lai(kidx0:kidx1,1:ntiles))
       IF (theOptions%UsePhenology) THEN
          DO i=1,ntiles
             theLand%Vegetation%lai(kidx0:kidx1,i) = &
                  MIN(theLand%Vegetation%lai(kidx0:kidx1,i),theLand%Vegetation%lai_max(kidx0:kidx1,i))
          END DO
          DO i = 1,ntiles
             theLand%Surface%veg_ratio(kidx0:kidx1,i) = theLand%Surface%veg_ratio_max(kidx0:kidx1) &
                  * (1.0_dp - exp(-theLand%Vegetation%lai(kidx0:kidx1,i) &
                               / theLand%LctLibrary%clumpinessFactor(theLand%Surface%cover_type(kidx0:kidx1,i))))
          END DO
       ELSE
          CALL get_from_climatology(nidx, ntiles, theLand%Surface%veg_fract(kidx0:kidx1,0:13), &
               theLand%Surface%veg_ratio(kidx0:kidx1,1:ntiles))
       END IF
       theLand%Surface%box_veg_ratio(kidx0:kidx1,1:ntiles) = &
            theLand%Surface%veg_ratio(kidx0:kidx1,1:ntiles) * theLand%Surface%cover_fract(kidx0:kidx1,1:ntiles)
    END IF

    ! Calculate albedo

    IF (.NOT. theOptions%UseAlbedo) THEN
      IF (lstart) THEN
        CALL update_land_surface_fast(nidx, theLand%LctLibrary, theLand%Surface, &
             theLand%Soil%surface_temperature(kidx0:kidx1,:),                    &
             theLand%Soil%snow_fract         (kidx0:kidx1,:),                    &
             theLand%Soil%albedo             (kidx0:kidx1,:),                    &
             theLand%Vegetation%snow_fract_canopy(kidx0:kidx1,:),                &
             theLand%Vegetation%lai          (kidx0:kidx1,:))
        theLand%Surface%albedo_vis(kidx0:kidx1,:) = theLand%Surface%albedo(kidx0:kidx1,:)
        theLand%Surface%albedo_nir(kidx0:kidx1,:) = theLand%Surface%albedo(kidx0:kidx1,:)
      END IF
    ELSE
       IF (l_trigrad .OR. theOptions%Standalone .OR. theOptions%InterfaceTest) THEN
          IF (debug) CALL message('jsbach_interface', &
                                  'update_albedo_diag called at radiation time step')
          CALL update_albedo_diag(              nidx,               &
               theLand%Surface%ntiles,                              &
               zsw_vis_net                         (1:nidx),        &
               zsw_nir_net                         (1:nidx),        &
               theLand%Surface%albedo_vis          (kidx0:kidx1,:), &
               theLand%Surface%albedo_nir          (kidx0:kidx1,:), &
               theLand%Surface%albedo              (kidx0:kidx1,:))
       END IF
    END IF

    ! Diagnose weighted albedo (of previous pre-radiation timestep)
    CALL average_tiles(theLand%Surface%albedo_vis(kidx0:kidx1,:), theLand%Surface%is_present(kidx0:kidx1,:), &
                       theLand%Surface%cover_fract(kidx0:kidx1,:), zalbedo_vis(:))
    CALL average_tiles(theLand%Surface%albedo_nir(kidx0:kidx1,:), theLand%Surface%is_present(kidx0:kidx1,:), &
                       theLand%Surface%cover_fract(kidx0:kidx1,:), zalbedo_nir(:))

    CALL average_tiles(theLand%Surface%albedo_vis_nosno_tiles(kidx0:kidx1,:), theLand%Surface%is_present(kidx0:kidx1,:), &
                       theLand%Surface%cover_fract(kidx0:kidx1,:), zalbedo_vis_nosno(:))
    CALL average_tiles(theLand%Surface%albedo_nir_nosno_tiles(kidx0:kidx1,:), theLand%Surface%is_present(kidx0:kidx1,:), &
                       theLand%Surface%cover_fract(kidx0:kidx1,:), zalbedo_nir_nosno(:))

    ! Determine total solar shortwave radiation

    swnet(:) = zsw_vis_net(:) + zsw_nir_net(:)

    ! Compute fraction of direct part in PAR 

    sw_par_fract_direct(:) = 1._dp - zsw_par_frac_diffuse(:)

    theLand%Surface%swdown_acc(kidx0:kidx1) = theLand%Surface%swdown_acc(kidx0:kidx1) + &
                                              (zsw_vis_net(:) / (1.0_dp - zalbedo_vis(:)) + &
                                              zsw_nir_net(:) / (1.0_dp - zalbedo_nir(:))) * delta_time
    theLand%Surface%swdown_reflect_acc(kidx0:kidx1) = theLand%Surface%swdown_reflect_acc(kidx0:kidx1) + &
                                              (zalbedo_vis(:) * zsw_vis_net(:) / (1.0_dp - zalbedo_vis(:)) + &
                                              zalbedo_nir(:) * zsw_nir_net(:) / (1.0_dp - zalbedo_nir(:))) * delta_time

    theLand%Surface%swdown_absorb_nosno_acc(kidx0:kidx1) = theLand%Surface%swdown_absorb_nosno_acc(kidx0:kidx1)           &
                                            +((1._dp - zalbedo_vis_nosno(:)) * zsw_vis_net(:) / (1.0_dp - zalbedo_vis(:)) &
                                            + (1._dp - zalbedo_nir_nosno(:)) * zsw_nir_net(:) / (1.0_dp - zalbedo_nir(:)))&
                                            * delta_time
    IF (theOptions%ls3m > 0) THEN
      theLand%Surface%swdown_vis_acc(kidx0:kidx1) = theLand%Surface%swdown_vis_acc(kidx0:kidx1) &
                                                  + zsw_vis_net(:) / (1.0_dp - zalbedo_vis(:)) * delta_time

      theLand%Surface%swdown_nir_acc(kidx0:kidx1) = theLand%Surface%swdown_nir_acc(kidx0:kidx1) &
                                                  + zsw_nir_net(:) / (1.0_dp - zalbedo_nir(:)) * delta_time

      theLand%Surface%sw_vis_absorb_nosno_acc(kidx0:kidx1) = theLand%Surface%sw_vis_absorb_nosno_acc(kidx0:kidx1) &
                              + ((1._dp - zalbedo_vis_nosno(:)) * zsw_vis_net(:) / (1.0_dp - zalbedo_vis(:))) * delta_time

      theLand%Surface%sw_nir_absorb_nosno_acc(kidx0:kidx1) = theLand%Surface%sw_nir_absorb_nosno_acc(kidx0:kidx1) &
                              + ((1._dp - zalbedo_nir_nosno(:)) * zsw_nir_net(:) / (1.0_dp - zalbedo_nir(:))) * delta_time
    ENDIF

    zCO2_flux_npp       = 0._dp
    zCO2_flux_npp_test  = 0._dp
    zCO2_flux_soilresp  = 0._dp
    zCO2_flux_herbivory = 0._dp

    zN2O_flux_mineraliz     = 0._dp
    zN2O_flux_depfix        = 0._dp
    zN2O_flux_nfert         = 0._dp
    zN2O_flux_grazing       = 0._dp
        
    ! Fist call to BETHY model (water unlimited case): 
    ! calculate unstressed stomatal conductance at prescribed leaf internal CO2 concentration
    IF (theOptions%UseBethy) THEN  

       IF (debug) CALL message('jsbach_interface_1d','Calling bethy (first time)')
       CALL update_bethy(nidx, theLand%Domain, theLand%Surface%is_vegetation(kidx0:kidx1,1:ntiles), & 
            theLand%Bethy, &
            theLand%Surface%cover_type(kidx0:kidx1,:), &
            theLand%LctLibrary, &
            theOptions%UseAlbedo, &                                    !! Flag (true if interactive albedo is used)
            .FALSE.,                                               &   !! Flag for water limitation in photosynthesis
            theLand%Vegetation%lai(kidx0:kidx1,1:ntiles),          &   !! Leaf area index
            zsw_par_down         (:),                              &   !! incoming shortwave radiation from PAR-band
            sw_par_fract_direct  (:),                              &   !! fraction of direct radiation in PAR-band
            zpressure            (:),                              &   !! Surface pressure
            ztemp_air            (:),                              &   !! canopy temperature (= air temperature ?)
            theLand%Soil%albedo(kidx0:kidx1,1:ntiles),             &   !! Soil albedo
            zCO2                 (:),                              &   !! CO2 concentration air
            theLand%Vegetation%canopy_conductance_bethy(kidx0:kidx1,1:ntiles), &   !! Unlimited canopy conductance (output)
            theLand%Vegetation%canopy_conductance_limited(kidx0:kidx1,1:ntiles), &   !! dummy
            !stopCO2_flag(kidx0:kidx1,1:ntiles), & !, stopCO2_flag_leaf &
            lstart, lresume, &
            theOptions%UseLimitCO2Fert)
       IF (debug) CALL message('jsbach_interface_1d','Leaving bethy (first time)')
    ELSE

       ! Compute canopy resistance without water stress

       IF (debug .AND. theOptions%UseEchamLand) CALL message('jsbach_interface_1d','Calling unstressed_canopy_cond_par')

       DO itile=1,ntiles
          IF (theOptions%UseEchamLand) THEN
             CALL unstressed_canopy_cond_par(theLand%Vegetation%lai(kidx0:kidx1,itile), zsw_vis_net(:), &
                  theLand%Vegetation%canopy_conductance_bethy(kidx0:kidx1,itile))
!!$          ELSE IF (theOptions%UseVic) THEN
!!$             CALL unstressed_canopy_resistance_rmin( &
!!$                  theLand%Vegetation% lai                (kidx0:kidx1,itile),                             &
!!$                  theLand%LctLibrary% CanopyResistanceMin(theLand%Surface%cover_type(kidx0:kidx1,itile)), &
!!$                  theLand%Vegetation% canopy_resistance  (kidx0:kidx1,itile))
          END IF
       END DO

       theLand%Vegetation%canopy_conductance_bethy(kidx0:kidx1,1:ntiles) = 1.e-5_dp   ! XXX

    END IF
    canopy_conductance_unstressed(1:nidx,1:ntiles)=theLand%Vegetation%canopy_conductance_bethy(kidx0:kidx1,1:ntiles) 
    
    ! Note: theLand%Vegetation%canopy_resistance now contains the canopy resistance for the case that there's no water stress.
    ! It will be used to compute potential evapotranspiration (but not in the ECHAM formulation). Afterwards, it will be 
    ! corrected for water stress and used for the computation of actual evapotranspiration from vegetation.

    ! --- Run the land surface scheme ---------------------------------------------------------------------------------------

    If (theOptions%Standalone) THEN
      CALL update_soil(nidx, theLand%Surface, theLand%Soil, theLand%SoilParam, &
                     theOptions%UseDynveg .or. theOptions%UseDisturbance, &
                     theOptions%Standalone, &
                     theLand%Vegetation%canopy_conductance_bethy(kidx0:kidx1,:), &
                     cover_fract_old(1:nidx,1:theLand%ntiles),                   &
                     theLand%Vegetation%lai(kidx0:kidx1,:), &
                     theLand%Vegetation%snow_depth_canopy(kidx0:kidx1,:), &
                     theLand%Vegetation%snow_fract_canopy(kidx0:kidx1,:), &
                     zcdrag, zetAcoef, zetBcoef, zeqAcoef, zeqBcoef, ztemp_air, zqair,  &
                     zpressure, zwind, zwind10, &
                     zlwdown, swnet, zprecip_rain, zprecip_snow, zcair, zcsat, &
                     p_echam_zchl, ztte_corr, &
                     theLand%Vegetation%canopy_conductance_limited(kidx0:kidx1,:), &
                     zglac_runoff_evap, zsurf_runoff_hd, zdrainage_hd, &
                     theOptions%HeightWind, theOptions%HeightHumidity, &
                     theOptions%CoefRilTM1, theOptions%CoefRilT, theOptions%CoefRilTP1,&
                     theLand%Soil%z0m(kidx0:kidx1), theLand%Soil%z0h(kidx0:kidx1))
    ELSE
      CALL update_soil(nidx, theLand%Surface, theLand%Soil, theLand%SoilParam, &
                     theOptions%UseDynveg .or. theOptions%UseDisturbance, &
                     theOptions%Standalone, &
                     theLand%Vegetation%canopy_conductance_bethy(kidx0:kidx1,:), &
                     cover_fract_old(1:nidx,1:theLand%ntiles),                   &
                     theLand%Vegetation%lai(kidx0:kidx1,:), &
                     theLand%Vegetation%snow_depth_canopy(kidx0:kidx1,:), &
                     theLand%Vegetation%snow_fract_canopy(kidx0:kidx1,:), &
                     zcdrag, zetAcoef, zetBcoef, zeqAcoef, zeqBcoef, ztemp_air, zqair,  &
                     zpressure, zwind, zwind10, &
                     zlwdown, swnet, zprecip_rain, zprecip_snow, zcair, zcsat, &
                     p_echam_zchl, ztte_corr, &
                     theLand%Vegetation%canopy_conductance_limited(kidx0:kidx1,:), &
                     zglac_runoff_evap, zsurf_runoff_hd, zdrainage_hd)
    END IF
!!$    ,  zlatent, zsensible, &
!!$                     zevap_pot, zevap_act)

!!$    WHERE (theLand%Vegetation%canopy_resistance(kidx0:kidx1,:) > 0.9_dp*HUGE(1._dp)) &
!!$         theLand%Vegetation%canopy_resistance(kidx0:kidx1,:) = 1.E20_dp

   IF (theOptions%WithHD) THEN
      CALL hydrology_collect_land(kidx0, kidx1, nidx, zsurf_runoff_hd, zdrainage_hd, zglac_runoff_evap)
   END IF

    ! Second call to BETHY model (water limited case): 
    ! calculate net assimilation at prescribed unstressed stomatal conductance

   IF (theOptions%UseBethy) THEN
      IF (debug) CALL message('jsbach_interface_1d','Calling bethy (second time)')
      CALL update_bethy(nidx, theLand%Domain, theLand%Surface%is_vegetation(kidx0:kidx1,1:ntiles), &
            theLand%Bethy, &
            theLand%Surface%cover_type(kidx0:kidx1,:), &
            theLand%LctLibrary, &
            theOptions%UseAlbedo, &                                    !! Flag (true if interactive albedo is used)
            .TRUE.,                                                &   !! Flag for water limited photosynthesis
            theLand%Vegetation%lai(kidx0:kidx1,1:ntiles),          &   !! Leaf area index
            zsw_par_down         (:),                              &   !! incoming shortwave radiation from PAR-band
            sw_par_fract_direct  (:),                              &   !! fraction of direct radiation in PAR-band
            zpressure            (:),                              &   !! Surface pressure
            ztemp_air            (:),                              &   !! canopy temperature (= lowest level air temperature)
            theLand%Soil%albedo(kidx0:kidx1,1:ntiles),             &   !! Soil albedo
            zCO2                 (:),                              &   !! CO2 concentration air
            theLand%Vegetation%canopy_conductance_bethy(kidx0:kidx1,1:ntiles), &   !! Canopy conductance (output, 
                                               !! should be the same as limited conductance input on following line)
            theLand%Vegetation%canopy_conductance_limited(kidx0:kidx1,1:ntiles), &   !! Limited canopy conductance (input)
            lstart, lresume, & 
            theOptions%UseLimitCO2Fert)
 
 DO itile=1,ntiles ! HW calculation of REW 
   relative_extractable_water(kidx0:kidx1,itile)=&
 (theLand%Soil%moisture(kidx0:kidx1,itile)/theLand%soilparam%RootDepth(kidx0:kidx1/itile)-theLand%soilparam%WiltingPoint(kidx0:kidx1))/&
 (theLand%soilparam%FieldCapacity(kidx0:kidx1)-theLand%soilparam%WiltingPoint(kidx0:kidx1)) 
 ! moisture in m; WiltingPoint and FieldCapacity in m/m
 END DO     
      CALL update_cbalance_bethy(   &
           theLand%Cbalance,        &
           theLand%Nbalance,        &
           theOptions%WithNitrogen, &
           theOptions%WithYasso, &
           theOptions%UseExternalLandcoverMaps, &
           theOptions%UseLandUseTransitions, &
           theOptions%UseForestManagement, &
           theLand%LctLibrary,      &
           theLand%Surface,         &
           theOptions%lcc_scheme,   &
           theLand%Bethy%gross_assimilation(kidx0:kidx1,1:ntiles),        & ! input: GPP rel. to ground area
           theLand%Bethy%dark_respiration(kidx0:kidx1,1:ntiles),          & ! input: dark respiration rel. to ground area
           theLand%Soil%soil_temperature(kidx0:kidx1,3,1),                & ! input (third layer, only one tile because temperature
                                                                            !        is identical on all tiles)
           theLand%Soil%relative_moisture(kidx0:kidx1,1:ntiles),          & ! input
           theLand%Vegetation%lai(kidx0:kidx1,1:ntiles),                  & ! Leaf area index (input)
           theLand%Vegetation%LAI_max(kidx0:kidx1,1:ntiles),              & ! Maximum value of leaf area index
           theLand%Vegetation%veg_fract_correction(kidx0:kidx1,1:ntiles), & ! correction factor 1-exp(-LAI_max/2)
           ztemp_air(:), zprecip_tot(:),                                  & ! variables needed with yasso
           zdrainage_hd(1:nidx)/delta_time,                               & ! drainage in [m/s] --> for N-leaching
           theLand%soilparam%MaxMoisture(kidx0:kidx1),                    & ! max water content in soil bucket [m]
           zCO2_flux_npp(1:nidx),                                         & ! output: [kg(CO2)/m^2(gridbox) s]
           zCO2_flux_soilresp(1:nidx),                                    & ! output: [kg(CO2)/m^2(gridbox) s]
           zCO2_flux_herbivory(1:nidx),                                   & ! output: [kg(CO2)/m^2(gridbox) s]
           zCO2_flux_npp_test(1:nidx),                                    & ! only used for conservation test
           zN2O_flux_mineraliz(1:nidx),                                   & ! N2O emission from mineralisation [kg(N2O)/m^2 s]
           zN2O_flux_depfix(1:nidx),                                      & ! N2O emission from deposition and fixation [kg(N2O)/m2 s]
           zN2O_flux_nfert(1:nidx),                                       & ! N2O emission from fertilization [kg(N2O)/m^2 s]
           zN2O_flux_grazing(1:nidx),                                     & ! N2O emission from herbivores [kg(N2O)/m^2 s]
           zN2_emission_ecosystem(1:nidx),                &                  ! net ecosystem N2 emissions [kg(N2)/m^2 s]
           leaf_shedding_rate(kidx0:kidx1,1:ntiles),& ! HW
           relative_extractable_water(kidx0:kidx1,1:ntiles))  ! HW update_cbalance_bethy passing REW to bethy

   ELSE
      zCO2_flux_npp(1:nidx) = 0._dp
      zCO2_flux_soilresp(1:nidx) = 0._dp
      zCO2_flux_herbivory(1:nidx) = 0._dp
      zCO2_flux_npp_test(1:nidx) = 0._dp
   END IF

   ! HW calculate relative wilting point
   relative_wilting_point(kidx0:kidx1)=theLand%soilparam%WiltingPoint(kidx0:kidx1)/theLand%soilparam%FieldCapacity(kidx0:kidx1)


   ! Update phenology, i.e. recompute LAI

   IF (theOptions%UsePhenology) THEN

      IF (debug) CALL message('jsbach_interface_1d','Calling update_phenology')
      SELECT CASE(theOptions%PhenoScheme)
      CASE('LOGROP')
         CALL update_phenology(nidx, &
              pheno_type(:, 1:ntiles), &
              theLand%Vegetation%lai(kidx0:kidx1, 1:ntiles), &
              theLand%Vegetation%lai_max(kidx0:kidx1, 1:ntiles), &
              (ztemp_air(:)-tmelt), &
              theLand%Soil%relative_moisture(kidx0:kidx1,1:ntiles), &
              theLand%cbalance%NPP_rate(kidx0:kidx1,1:ntiles), &
              specificLeafArea_C(1:nidx,1:ntiles), &
              theLand%Domain%lat(kidx0:kidx1), &
              relative_wilting_point(kidx0:kidx1),  &  ! HW
              leaf_shedding_rate(kidx0:kidx1,1:ntiles) & ! HW
              )
      CASE('KNORR')
         CALL update_phenology_ccdas(nidx, theLand%LctLibrary, &
              theLand%Vegetation%lai(kidx0:kidx1, 1:ntiles), &
              theLand%Soil%relative_moisture(kidx0:kidx1,1:ntiles), &
              pheno_type(:, 1:ntiles),&
              pft_type(:,1:ntiles),&
              ztemp_air(:),&                   
              theLand%Domain%sinlat(kidx0:kidx1), &
              theLand%Domain%coslat(kidx0:kidx1), &
              zlwdown(:),&
              swnet(:),&
         ! TR swnet*1.2 is a proxy for swdown and used to calculate potential transpiration in update_phenology_ccdas
         ! TR This should be replaced soon by potential transpiration calculated in update_soil
              (swnet(:) * 1.2_dp), &
              theLand%Soil%surface_temperature(kidx0:kidx1,:),&
              theLand%Bethy%apar_acc(kidx0:kidx1,:), &
              theLand%Bethy%par_acc(kidx0:kidx1,:),&
              zpressure(:),&
             ! zvapor_pressure(:),&
              zwind10(:),&
              theLand%Vegetation%veg_height(kidx0:kidx1, 1:ntiles), &
              canopy_conductance_unstressed(1:nidx,1:ntiles), &
              theLand%Vegetation%lai_max(kidx0:kidx1, 1:ntiles) ,&
              zqair(:) &
                ) 
      CASE default 
         CALL finish('jsbach_interfacei_1d','Unknown phenology scheme')
      END SELECT
      IF (debug) CALL message('jsbach_interface_1d','Returned from update_phenology')

      ! calculate the grid cell fraction currently covered with vegetation
      DO i = 1,ntiles
         theLand%Surface%veg_ratio(kidx0:kidx1,i) = theLand%Surface%veg_ratio_max(kidx0:kidx1) &
              * (1.0_dp - exp(-theLand%Vegetation%lai(kidx0:kidx1,i) &
                               / theLand%LctLibrary%clumpinessFactor(theLand%Surface%cover_type(kidx0:kidx1,i))))
         theLand%Surface%box_veg_ratio(kidx0:kidx1,i) = &
              theLand%Surface%veg_ratio(kidx0:kidx1,i) * theLand%Surface%cover_fract(kidx0:kidx1,i)
      END DO

   ELSE

      IF (debug) CALL message('jsbach_interface_1d','Get LAI from climatology')

      CALL get_from_climatology(nidx, ntiles, theLand%Vegetation%lai_clim(kidx0:kidx1,1:ntiles,0:13), &
           theLand%Vegetation%lai(kidx0:kidx1,1:ntiles))

      CALL get_from_climatology(nidx, ntiles, theLand%Surface%veg_fract(kidx0:kidx1,0:13), &
           theLand%Surface%veg_ratio(kidx0:kidx1,1:ntiles))
   END IF

   WHERE (theLand%Vegetation%lai(kidx0:kidx1,1:ntiles) < EPSILON(1._dp)) &
        theLand%Vegetation%lai(kidx0:kidx1,1:ntiles) = 0._dp

   ! Calculate albedo

   IF (.NOT. theOptions%UseAlbedo) THEN
      CALL update_land_surface_fast(nidx, theLand%LctLibrary, theLand%Surface, &
           theLand%Soil%surface_temperature(kidx0:kidx1,:),                    &
           theLand%Soil%snow_fract         (kidx0:kidx1,:),                    &
           theLand%Soil%albedo             (kidx0:kidx1,:),                    &
           theLand%Vegetation%snow_fract_canopy(kidx0:kidx1,:),                &
           theLand%Vegetation%lai          (kidx0:kidx1,:))
      theLand%Surface%albedo_vis(kidx0:kidx1,:) = theLand%Surface%albedo(kidx0:kidx1,:)
      theLand%Surface%albedo_nir(kidx0:kidx1,:) = theLand%Surface%albedo(kidx0:kidx1,:)
   ELSE
      IF (l_trigradm1 .OR. theOptions%Standalone .OR. theOptions%InterfaceTest) THEN
         IF (debug) CALL message('jsbach_interface', 'update_albedo called l_trigradm1')
         IF (theOptions%WithYasso) THEN
            CALL update_albedo(theLand%LctLibrary,                 &
              nidx,                                                &
              theLand%Surface%ntiles,                              &
              theLand%Surface%cover_type          (kidx0:kidx1,:), &
              theLand%Surface%cover_fract         (kidx0:kidx1,:), &
              theLand%Vegetation%veg_fract_correction(kidx0:kidx1,:), &
              theLand%Surface%is_glacier          (kidx0:kidx1,:), &
              theLand%Surface%is_forest           (kidx0:kidx1,:), &
              theLand%Surface%veg_ratio_max       (kidx0:kidx1),   &
              zsw_vis_net                         (1:nidx),        &
              zsw_nir_net                         (1:nidx),        &
              theLand%Soil%snow_age               (kidx0:kidx1),   &
              theLand%Soil%surface_temperature    (kidx0:kidx1,:), &
              theLand%Soil%snow_fract             (kidx0:kidx1,:), &
              theLand%Soil%albedo_soil_vis        (kidx0:kidx1,:), &
              theLand%Soil%albedo_soil_nir        (kidx0:kidx1,:), &
              theLand%Soil%albedo_vegetation_vis  (kidx0:kidx1,:), &
              theLand%Soil%albedo_vegetation_nir  (kidx0:kidx1,:), &
              theLand%Vegetation%lai              (kidx0:kidx1,:), &
              theLand%Vegetation%snow_fract_canopy(kidx0:kidx1,:), &
              theLand%Surface%albedo_vis          (kidx0:kidx1,:), &
              theLand%Surface%albedo_nir          (kidx0:kidx1,:), &
              theLand%Surface%albedo_vis_nosno_tiles(kidx0:kidx1,:), &
              theLand%Surface%albedo_nir_nosno_tiles(kidx0:kidx1,:), &
              theOptions%WithYasso,                                &
              YCpool_acid_ag2 = theLand%Cbalance%YCpool_acid_ag2(kidx0:kidx1,:), &
              YCpool_water_ag2 = theLand%Cbalance%YCpool_water_ag2(kidx0:kidx1,:), &
              YCpool_ethanol_ag2 = theLand%Cbalance%YCpool_ethanol_ag2(kidx0:kidx1,:), &
              YCpool_nonsoluble_ag2 = theLand%Cbalance%YCpool_nonsoluble_ag2(kidx0:kidx1,:), &
              YCpool_acid_bg1 = theLand%Cbalance%YCpool_acid_bg1(kidx0:kidx1,:), &
              YCpool_water_bg1 = theLand%Cbalance%YCpool_water_bg1(kidx0:kidx1,:), &
              YCpool_ethanol_bg1 = theLand%Cbalance%YCpool_ethanol_bg1(kidx0:kidx1,:), &
              YCpool_nonsoluble_bg1 = theLand%Cbalance%YCpool_nonsoluble_bg1(kidx0:kidx1,:), &
              YCpool_acid_bg2 = theLand%Cbalance%YCpool_acid_bg2(kidx0:kidx1,:), &
              YCpool_water_bg2 = theLand%Cbalance%YCpool_water_bg2(kidx0:kidx1,:), &
              YCpool_ethanol_bg2 = theLand%Cbalance%YCpool_ethanol_bg2(kidx0:kidx1,:), &
              YCpool_nonsoluble_bg2 = theLand%Cbalance%YCpool_nonsoluble_bg2(kidx0:kidx1,:), &
              YCpool_humus_1 = theLand%Cbalance%YCpool_humus_1(kidx0:kidx1,:), &
              YCpool_humus_2 = theLand%Cbalance%YCpool_humus_2(kidx0:kidx1,:))
         ELSE
            CALL update_albedo(theLand%LctLibrary,                 &
              nidx,                                                &
              theLand%Surface%ntiles,                              &
              theLand%Surface%cover_type          (kidx0:kidx1,:), &
              theLand%Surface%cover_fract         (kidx0:kidx1,:), &
              theLand%Vegetation%veg_fract_correction(kidx0:kidx1,:), &
              theLand%Surface%is_glacier          (kidx0:kidx1,:), &
              theLand%Surface%is_forest           (kidx0:kidx1,:), &
              theLand%Surface%veg_ratio_max       (kidx0:kidx1),   &
              zsw_vis_net                         (1:nidx),        &
              zsw_nir_net                         (1:nidx),        &
              theLand%Soil%snow_age               (kidx0:kidx1),   &
              theLand%Soil%surface_temperature    (kidx0:kidx1,:), &
              theLand%Soil%snow_fract             (kidx0:kidx1,:), &
              theLand%Soil%albedo_soil_vis        (kidx0:kidx1,:), &
              theLand%Soil%albedo_soil_nir        (kidx0:kidx1,:), &
              theLand%Soil%albedo_vegetation_vis  (kidx0:kidx1,:), &
              theLand%Soil%albedo_vegetation_nir  (kidx0:kidx1,:), &
              theLand%Vegetation%lai              (kidx0:kidx1,:), &
              theLand%Vegetation%snow_fract_canopy(kidx0:kidx1,:), &
              theLand%Surface%albedo_vis          (kidx0:kidx1,:), &
              theLand%Surface%albedo_nir          (kidx0:kidx1,:), &
              theLand%Surface%albedo_vis_nosno_tiles(kidx0:kidx1,:), &
              theLand%Surface%albedo_nir_nosno_tiles(kidx0:kidx1,:), &
              theOptions%WithYasso,                                &
              Cpool_litter_green_ag = theLand%Cbalance%Cpool_litter_green_ag(kidx0:kidx1,:), &
              Cpool_slow = theLand%Cbalance%Cpool_slow(kidx0:kidx1,:))
         END IF
      END IF
      CALL average_tiles(theLand%Surface%albedo_vis(kidx0:kidx1,:), theLand%Surface%is_present(kidx0:kidx1,:), &
           theLand%Surface%cover_fract(kidx0:kidx1,:), zalbedo_vis(:))
      CALL average_tiles(theLand%Surface%albedo_nir(kidx0:kidx1,:), theLand%Surface%is_present(kidx0:kidx1,:), &
           theLand%Surface%cover_fract(kidx0:kidx1,:), zalbedo_nir(:))
   END IF
   IF (debug) CALL message('jsbach_interface', 'return from update_albedo l_trigradm1') !KN 

   ! Compute diagnostic output

   CALL land_surface_diagnostics(theLand%Surface)
   CALL soil_diagnostics(theLand%Surface, theLand%Soil)
   CALL veg_diagnostics (theLand%Surface, theLand%Vegetation)
   IF (theOptions%UseBethy) THEN
      ! calculate net land co2 flux and convert to mol CO2 / m2 (gridbox) / s
      theLand%Bethy%CO2_flux_net_acc(kidx0:kidx1) = theLand%Bethy%CO2_flux_net_acc(kidx0:kidx1) &
           + (zCO2_flux_npp + zCO2_flux_soilresp + zCO2_flux_herbivory + &
              zCO2_emission_landcover_change + zCO2_emission_harvest + zCO2_flux_dynveg ) / molarMassCO2_kg * delta_time
      theLand%Bethy%CO2_flux_herbivory_acc(kidx0:kidx1) = theLand%Bethy%CO2_flux_herbivory_acc(kidx0:kidx1) &
           + zCO2_flux_herbivory / molarMassCO2_kg * delta_time
      theLand%Bethy%CO2_emission_landcover_change_acc(kidx0:kidx1) = theLand%Bethy%CO2_emission_landcover_change_acc(kidx0:kidx1) &
           + zCO2_emission_landcover_change / molarMassCO2_kg * delta_time
      theLand%Bethy%CO2_emission_harvest_acc(kidx0:kidx1) = theLand%Bethy%CO2_emission_harvest_acc(kidx0:kidx1) &
           + zCO2_emission_harvest / molarMassCO2_kg * delta_time
      theLand%Bethy%CO2_flux_dynveg_acc(kidx0:kidx1) = theLand%Bethy%CO2_flux_dynveg_acc(kidx0:kidx1) &
           + zCO2_flux_dynveg / molarMassCO2_kg * delta_time
      CALL bethy_diagnostics(theLand%Surface, theLand%Bethy, theLand%Vegetation%veg_fract_correction(kidx0:kidx1,1:ntiles))
   END IF
   IF (debug) CALL message('jsbach_interface', 'diagnostic output computed') !KN  

   ! Compute roughness length and grid box areal average values
   CALL update_land_boundary_up(theLand%LctLibrary, theLand%Surface, theLand%Soil, theLand%Vegetation, &
        theOptions%UseRoughnessLAI, theOptions%UseRoughnessOro, theOptions%Standalone,&
        theLand%Vegetation%lai(kidx0:kidx1,1:ntiles), theLand%SoilParam%Roughness(kidx0:kidx1), zz0h, zz0m, &
        ztemp_soil_new, zqsurf, zevap_act, zevap_pot,                            &
        zsensible, zlatent, ztsoil_rad, zdry_static_energy, zsoil_wetness, zsnow_depth,   &
        zskin_res, zveg_height)

   theLand%Soil%z0m(kidx0:kidx1) = zz0m
   theLand%Soil%z0h(kidx0:kidx1) = zz0h

   ! Global Carbon and Nitrogen conservation tests
   CALL check_C_conservation(zCO2_flux_npp_test(:) + zCO2_flux_soilresp(:) + zCO2_flux_herbivory(:) &
        + zCO2_flux_dynveg(:) + zCO2_emission_landcover_change(:) + zCO2_emission_harvest(:))
   CALL check_N_conservation( (zN2O_flux_mineraliz(:) + zN2O_flux_grazing(:)) * molarMassN2_kg/molarMassN2O_kg     &
        + zN2_emission_ecosystem + zN2_flux_dynveg(:) + zN2_emission_landcover_change(:) + zN2_emission_harvest(:))

   ! Unpack output fields
   IF (PRESENT(evap_act)) evap_act                             = UNPACK(zevap_act,mask,0.0_dp)
   IF (PRESENT(evap_pot)) evap_pot                             = UNPACK(zevap_pot,mask,0.0_dp)
   IF (PRESENT(sensible)) sensible                             = UNPACK(zsensible,mask,0.0_dp)
   IF (PRESENT(latent)) latent                                 = UNPACK(zlatent,mask,0.0_dp)
   IF (PRESENT(CO2_flux_npp)) CO2_flux_npp                     = UNPACK(zCO2_flux_npp, mask, 0.0_dp)
   IF (PRESENT(CO2_flux_soilresp)) CO2_flux_soilresp           = UNPACK(zCO2_flux_soilresp, mask, 0.0_dp)
   IF (PRESENT(CO2_flux_herbivory)) CO2_flux_herbivory         = UNPACK(zCO2_flux_herbivory, mask, 0.0_dp)
   IF (PRESENT(CO2_flux_dynveg)) CO2_flux_dynveg               = UNPACK(zCO2_flux_dynveg, mask, 0.0_dp)
   IF (PRESENT(CO2_emission_lcc)) CO2_emission_lcc             = UNPACK(zCO2_emission_landcover_change, mask, 0.0_dp)
   IF (PRESENT(CO2_emission_harvest)) CO2_emission_harvest     = UNPACK(zCO2_emission_harvest, mask, 0.0_dp)
   IF (PRESENT(N2O_flux_mineraliz)) N2O_flux_mineraliz         = UNPACK(zN2O_flux_mineraliz, mask, 0.0_dp)                !dk new
   IF (PRESENT(N2O_flux_depfix)) N2O_flux_depfix               = UNPACK(zN2O_flux_depfix, mask, 0.0_dp)
   IF (PRESENT(N2O_flux_nfert)) N2O_flux_nfert                 = UNPACK(zN2O_flux_nfert, mask, 0.0_dp)
   IF (PRESENT(N2O_flux_grazing)) N2O_flux_grazing             = UNPACK(zN2O_flux_grazing, mask, 0.0_dp)
   IF (PRESENT(tsoil_rad)) tsoil_rad                           = UNPACK(ztsoil_rad,mask,0.0_dp)
   IF (PRESENT(temp_soil_new)) temp_soil_new                   = UNPACK(ztemp_soil_new,mask,0.0_dp)
   IF (PRESENT(qsurf)) qsurf                                   = UNPACK(zqsurf,mask,0.0_dp)
   IF (PRESENT(albedo_vis)) albedo_vis                         = UNPACK(zalbedo_vis,mask,0.0_dp)
   IF (PRESENT(albedo_nir)) albedo_nir                         = UNPACK(zalbedo_nir,mask,0.0_dp)
   IF (PRESENT(emis)) emis                                     = UNPACK(zemis,mask,0.0_dp)
   IF (PRESENT(z0h)) z0h                                       = UNPACK(zz0h,mask,0.0_dp)
   IF (PRESENT(z0m)) z0m                                       = UNPACK(zz0m,mask,0.0_dp)
   IF (PRESENT(cair)) cair                                     = UNPACK(zcair,mask,0.0_dp)
   IF (PRESENT(csat)) csat                                     = UNPACK(zcsat,mask,0.0_dp)
   IF (PRESENT(surf_dry_static_energy)) surf_dry_static_energy = UNPACK(zdry_static_energy,mask,2.276e+05_dp)
   IF (PRESENT(soil_wetness)) soil_wetness                     = UNPACK(zsoil_wetness,mask,0.0_dp)
   IF (PRESENT(snow_depth)) snow_depth                         = UNPACK(zsnow_depth,mask,0.0_dp)
   IF (PRESENT(skin_res)) skin_res                             = UNPACK(zskin_res,mask,0.0_dp)
   IF (PRESENT(veg_height)) veg_height                         = UNPACK(zveg_height,mask,0.0_dp)
   IF (PRESENT(tte_corr)) tte_corr                             = UNPACK(ztte_corr,mask,0.0_dp)

   IF (PRESENT(glacier_depth)) glacier_depth = UNPACK(get_soil_diag(nidx, 'glacier_depth'),mask,0.0_dp)
   IF (PRESENT(snow_melt)) snow_melt = UNPACK(get_soil_diag(nidx, 'snow_melt_inst'),mask,0.0_dp)
   IF (PRESENT(wsmax)) wsmax =         UNPACK(theLand%SoilParam%MaxMoisture(kidx0:kidx1),mask,0.0_dp)
   
   IF (debug) CALL message('jsbach_interface', 'output fields unpacked') !KN  

   ! kindex is allocated at the beginning of the interface through call to update_current_call
   DEALLOCATE(kindex)
   IF (debug) CALL message('jsbach_interface', 'variables deallocated') !KN

   ! Do some things if current time step is complete
   !    IF (theLand%Domain%LastBlock) THEN
   !
   !    ENDIF

 END SUBROUTINE jsbach_inter_1d

  SUBROUTINE jsbach_init(grid, domain, options, ntile &
#ifdef STANDALONE
                        , nproca, nprocb, npedim, forcing, driving &                    
#endif
                            )

    ! Initialization of the global grid, local PE domains, and all submodels of JSBACH

    USE mo_jsbach,        ONLY: jsbach_config, nml_unit
    USE mo_jsbach_grid,   ONLY: init_grid, init_domain
    USE mo_land_surface,  ONLY: init_land_surface
    USE mo_jsbach_lctlib, ONLY: init_lctlib
    USE mo_jsbach_veg,    ONLY: init_vegetation
    USE mo_soil,          ONLY: init_soil
    USE mo_phenology,     ONLY: init_phenology
    USE mo_phenology_ccdas, ONLY: init_phenology_ccdas
    USE mo_allocation,    ONLY: init_funct_allometry
    USE mo_forest_management, ONLY: init_forest_management
    USE mo_time_control,  ONLY: lfirst_cycle, init_manager, init_times, print_events, &
                                ec_manager_init, lresume
    USE mo_linked_list,   ONLY: LAND, TILES, print_stream
    USE mo_output,        ONLY: init_output, open_output_streams, &
                                land_table, ldiag_table, ldiag2_table, veg_table, nitrogen_table, yasso_table
    USE mo_memory_base,   ONLY: new_stream, default_stream_setting, ostreams, nstreams
    USE mo_mvstream,         ONLY : init_mvstream, mvstream_update_cache
    USE mo_filename,      ONLY: standard_rerun_file
#ifdef STANDALONE
    USE mo_filename,      ONLY: GRIB
    USE mo_echam_convect_tables,ONLY: init_convect_tables
    USE mo_gaussgrid,     ONLY: inigau
    USE mo_time_control,  ONLY: init_events
    USE mo_jsbalone,      ONLY: init_driving, drive_type
    USE mo_jsbalone_forcing, ONLY: init_forcing, forcing_type
    USE mo_input_interface,  ONLY: InputInitialize, INPUT_MODEL_JSBACH
#else
    USE mo_input_interface,  ONLY: InputInitialize, INPUT_MODEL_ECHAM
    USE mo_echam_yaxt,       ONLY: add_yaxt_gp_nlevs
#endif
    USE mo_zenith,        ONLY: init_zenith
    USE mo_cbal_landcover_change, ONLY: init_landcover_change
    USE mo_sub_nml,     ONLY : set_stream_element_nml, set_stream_nml
    USE mo_jsbach,      ONLY: nml_unit, veg_putdata
    USE mo_hydrology,   ONLY: init_hydrology

    TYPE(grid_type),    INTENT(out), OPTIONAL :: grid
    TYPE(domain_type),  INTENT(out), OPTIONAL :: domain
    TYPE(options_type), INTENT(out), OPTIONAL :: options
    INTEGER ,           INTENT(out), OPTIONAL :: ntile
#ifdef STANDALONE
    INTEGER,            INTENT(in),    OPTIONAL :: nproca
    INTEGER,            INTENT(in),    OPTIONAL :: nprocb
    INTEGER,            INTENT(in),    OPTIONAL :: npedim
    TYPE(forcing_type), INTENT(inout), OPTIONAL :: forcing
    TYPE(drive_type),   INTENT(inout), OPTIONAL :: driving
#endif

    INTEGER :: IO_timestep, istep, i

    IF (debug) CALL message('jsbach_init','Begin')
    CALL jsbach_config(theOptions)
    standard_rerun_file='rerun_'//theOptions%Experiment

    IF (theOptions%Standalone) THEN
       CALL get_dates(IO_timestep, istep)           ! get timestep
       CALL ec_manager_init(IO_timestep, istep)     ! time manager initialization
       IF (lfirst_cycle) CALL init_manager
       CALL init_times                              ! initialize all dates
    END IF

    CALL init_grid(theLand%Grid, &                                      ! Set global dimensions
                   theOptions%Standalone, lresume, theOptions%GridFile)

#ifdef STANDALONE
    IF (theOptions%Standalone) THEN
       CALL jsbalone_init_decomp(theLand%Grid%nlon, theLand%Grid%nlat, & ! Initialize decomposition
            theLand%Grid%mask, nproca, nprocb, npedim)
       CALL inigau
       CALL init_convect_tables

       CALL init_events

       ! GRIB output does not work correctly for arrays with a tile dimension if
       ! npedim (=nproma) differs from -1
       IF (theOptions%FileType == GRIB .AND. npedim /= -1) THEN
          CALL finish('jsbach_init', 'nproma has to be -1 with GRIB output')
       END IF
       CALL InputInitialize(INPUT_MODEL_JSBACH,TRIM(theOptions%input_verbose),file_name=TRIM(theOptions%GridFile), &
                            nprocx=nproca,nprocy=nprocb)
    END IF
#else
    CALL InputInitialize(INPUT_MODEL_ECHAM,TRIM(theOptions%input_verbose),file_name=TRIM(theOptions%GridFile))
#endif

    ! Add JSBACH I/O stream
    CALL new_stream(theLand%IO_stream, 'jsbach', filetype=theOptions%FileType, ztype=theOptions%FileZtype, &
                    lpost=theOptions%OutputModelState, lrerun=.TRUE.)
    CALL default_stream_setting(theLand%IO_stream, lpost=theOptions%OutputModelState, lrerun=.TRUE., repr=LAND, table=land_table)

    CALL new_stream(theLand%IO_diag_stream, 'land', filetype=theOptions%FileType, ztype=theOptions%FileZtype, &
                    lpost=.TRUE., lrerun=.FALSE.)
    CALL default_stream_setting(theLand%IO_diag_stream, lpost=.TRUE., lrerun=.FALSE., repr=LAND, table=ldiag_table)

    ! Stream that contains accumulated values weighted by land fraction, necessary to compute water and energy 
    ! balance in ECHAM postprocessing
    CALL new_stream(theLand%IO_accw, 'accw', filetype=theOptions%FileType, ztype=theOptions%FileZtype, &
                    lpost=.NOT. theOptions%Standalone, lrerun=.NOT. theOptions%Standalone)
    CALL default_stream_setting(theLand%IO_accw, lpost=.NOT. theOptions%Standalone, lrerun=.NOT. theOptions%Standalone, &
                    repr=LAND, table=ldiag2_table)

    CALL new_stream(theLand%IO_veg, 'veg', filetype=theOptions%FileType, ztype=theOptions%FileZtype, &
                    lpost=.TRUE., lrerun=.TRUE., interval=veg_putdata)
    CALL default_stream_setting(theLand%IO_veg, lpost=.TRUE., lrerun=.TRUE., repr=LAND, table=veg_table)

    IF (theOptions%withNitrogen) THEN
       CALL new_stream(theLand%IO_nitro, 'nitro', filetype=theOptions%FileType, ztype=theOptions%FileZtype, &
                       lpost=.TRUE., lrerun=.TRUE., lcontnorest=.TRUE., interval=veg_putdata)
       CALL default_stream_setting(theLand%IO_nitro, lpost=.TRUE., lrerun=.TRUE., repr=LAND, table=nitrogen_table)
    ELSE
       theLand%IO_nitro => theLand%IO_veg ! To secure correct passing of parameters (no variables will be added)
    ENDIF
    
    IF (theOptions%withYasso) THEN  !DSG 2.6.2013
       CALL new_stream(theLand%IO_yasso, 'yasso', filetype=theOptions%FileType, ztype=theOptions%FileZtype, &
                       lpost=.TRUE., lrerun=.TRUE., interval=veg_putdata)
       CALL default_stream_setting(theLand%IO_yasso, lpost=.TRUE., lrerun=.TRUE., repr=LAND, table=yasso_table)
    ELSE
       theLand%IO_yasso => theLand%IO_veg ! To secure correct passing of parameters (no variables will be added)
    ENDIF

    CALL init_domain(theLand%Grid, theLand%Domain, theOptions%FileType, theOptions%FileZtype, theLand%IO_stream)

    CALL init_lctlib(theOptions%LctlibFile, theLand%LctLibrary)       ! Gets number of land cover types

    !------------------------------------------------------------------------------------------------------
    ! Determine total number of tiles
    theLand%ntiles = theOptions%ntiles

    if (present(ntile)) ntile =  theLand%ntiles

    ! Here, other subroutines could be called which introduce additional sub-tiling structures 
    ! (e.g. snowbands, distributed precipitation). In this case, theLand%ntiles, theLand%Surface%ntiles,
    ! theLand%Soil%ntiles and theLand%Vegetation%ntiles would be updated here.
    ! ...
    !
    ! The number of soil tiles is set to the total number of land cover types (times additional sub-tiling)
    theLand%Surface%ntiles = theLand%ntiles
    theLand%Soil%ntiles = theLand%ntiles
    theLand%Vegetation%ntiles = theLand%ntiles
    IF (theOptions%UseBethy) THEN
       theLand%Cbalance%ntiles = theLand%ntiles
       theLand%Nbalance%ntiles = theLand%ntiles
       theLand%Bethy%ntiles = theLand%ntiles
    ENDIF

    !------------------------------------------------------------------------------------------------------
    ! Now that the total number of tiles is determined, the land surface, vegetation and soil modules can be initialized.
    CALL init_land_surface(theLand%Grid, theLand%Domain, theLand%LctLibrary, theLand%Surface, theOptions%SurfFile,  &
                         theOptions%UseDynveg .or. theOptions%UseDisturbance, theOptions%ls3m, &
                         theOptions%UseLanduseTransitions, &
                         theOptions%Standalone, lresume, theOptions%ReadCoverFract, theOptions%WithHD, &
                         theOptions%FileType, theOptions%FileZtype, theLand%IO_diag_stream, theLand%IO_stream)
    theLand%Domain%elev => theLand%Surface%elevation   !! Note: when this is a restart run, theLand%Surface%elevation has been
                                                       !!       allocated but not initialized
    ! default_stream_settings needs to be called here again, as TILES (index into IO_dim_ids) has only been defined in
    ! the previous call to init_land_surface (through add_dim(..., levtyp=70, indx=TILES) in  land_surface_init_io).
    CALL default_stream_setting(theLand%IO_stream, leveltype=TILES)
    CALL default_stream_setting(theLand%IO_veg, leveltype=TILES)
    CALL default_stream_setting(theLand%IO_yasso, leveltype=TILES)

    IF (test_stream) CALL init_test(theLand%Grid, theLand%Domain, theLand%ntiles)

    CALL init_vegetation(theLand%Grid, theLand%Domain, theLand%Vegetation, &
         theOptions%VegFile, lresume, theOptions%FileType, theOptions%FileZtype, theLand%IO_diag_stream, theLand%IO_stream)

    CALL init_soil(theLand%Grid, theLand%Domain, theLand%SoilParam, theLand%Soil, theOptions%SoilFile, &
         theOptions%Standalone, lresume, theOptions%UseAlbedo, theOptions%UseDynveg .OR. theOptions%UseDisturbance, &
         theOptions%FileType, theOptions%FileZtype, theLand%IO_diag_stream, theLand%IO_accw, theLand%IO_stream)
#ifndef STANDALONE
    CALL add_yaxt_gp_nlevs((/theLand%ntiles, theLand%Soil%nsoil, theLand%Soil%ntsoil/))
#endif
    IF (theOptions%UsePhenology) THEN 
       SELECT CASE(theOptions%PhenoScheme)
       CASE ('LOGROP')
          IF (debug) CALL message('jsbach_init','LoGro-P phenology scheme')
          CALL init_phenology(theLand%Grid%nland, theLand%Domain%nland, theLand%ntiles,theLand%Soil%nsoil, &
               lresume, theOptions%FileType, theOptions%FileZtype, theLand%IO_veg)
       CASE('KNORR')
          IF (debug) CALL message('jsbach_init','KNORR phenology scheme')
          CALL init_phenology_ccdas(theLand%Grid%nland, theLand%Domain%nland, theLand%ntiles,lresume, &
               theOptions%filetype, theLand%IO_veg,theLand%LctLibrary%nlct)
       CASE default
          CALL finish('jsbach_init','Unknown phenology scheme')
       END SELECT
    ENDIF

    IF (theOptions%UseFunctAllometry) THEN
        CALL init_funct_allometry(theLand%Grid, theLand%Domain, theLand%ntiles, theOptions%FileType, &
            theOptions%FileZtype, theLand%IO_veg)
    ENDIF

    IF (theOptions%UseBethy) THEN
       CALL init_cbalance_bethy(theLand%Grid, theLand%Domain,theLand%Cbalance, theLand%Nbalance, &
                                theOptions%WithNitrogen, theOptions%lcc_scheme, theOptions%UseForestManagement, &
                                theOptions%UseLanduseTransitions, &
                                theOptions%FileType, theOptions%FileZtype, theLand%IO_veg, theLand%IO_nitro, &
                                theLand%IO_yasso, theOptions%WithYasso)
       CALL init_bethy(theLand%Grid, theLand%Domain, theLand%Bethy, theOptions%FileType, theOptions%FileZtype, &
                       theLand%IO_diag_stream, theLand%IO_stream)
    ENDIF

    IF (theOptions%UseExternalLandcoverMaps .or. theOptions%UseLanduseTransitions) &
       CALL init_landcover_change(theLand%Grid%nland, theLand%Domain%nland, theLand%ntiles, &
                                  lresume, theOptions%UseLanduseTransitions, theOptions%WithNitrogen, &
                                  theOptions%lcc_scheme,theLand%landcover_change, &
                                  theOptions%FileType, theOptions%FileZtype, theLand%IO_veg, theLand%IO_nitro)

    IF (theOptions%UseForestManagement) &
       CALL init_forest_management(theLand%Grid, theLand%Domain, theLand%ntiles, theOptions%FileType, theLand%IO_veg, &
                                  & theOptions%UseRotationMap)

    IF (theOptions%UseDynveg .OR. theOptions%UseDisturbance) &
       CALL init_climbuf(theLand%Grid, theLand%Domain, theLand%ntiles, lresume, theOptions%FileType, theOptions%FileZtype, &
            theLand%Climbuf, theLand%IO_veg)
    CALL init_disturbance(theLand%Grid, theLand%Domain, theLand%lctLibrary, theLand%surface, lresume, &
                          theOptions%UseDisturbance, theOptions%WithYasso, theOptions%WithNitrogen, &
                          nml_unit, theOptions%FileType, theOptions%FileZtype, theLand%IO_veg, theLand%IO_disturbance)
    IF (theOptions%UseDynveg) &
       CALL init_dynveg(theOptions%UseDisturbance, theLand%Grid, theLand%Domain, &
                        theLand%ntiles, lresume, theOptions%FileType, theOptions%FileZtype, theLand%IO_veg)

#ifdef STANDALONE
    IF (theOptions%Standalone) THEN
       CALL init_forcing(theLand%Grid, theLand%Domain, theOptions, forcing)
       CALL init_driving(theLand%Grid, theLand%Domain, theOptions, driving)
    END IF
#endif
    IF (.NOT. theOptions%Standalone) CALL init_zenith(theLand%Domain)

    ! initialize HD model
    IF (theOptions%WithHD) THEN
       CALL init_hydrology (theLand%Grid, theLand%Domain, theOptions%Standalone, &
                            theOptions%FileType, theOptions%FileZtype, theLand%IO_stream)
    ENDIF

    ! initialize mvstream routine
    IF (theOptions%Standalone) CALL init_mvstream('namelist.jsbach')

    !
    !------------------------------------------------------------------------------------------------------
    ! Read restart files if this is a restarted run
    IF (lresume) THEN
       CALL jsbach_restart
    ENDIF

    ! *stepon_jsbach* needs to be called the first time after the initial or restart files are
    ! read, and before the interface is called the first time.
    ! For coupled runs, stepon_jsbach is called in *control* after *ioinitial* or *iorestart*
    ! and before *init_surface*
    IF (theOptions%Standalone) CALL stepon_jsbach !! Initializes cover fractions and veg_ratio

    CALL set_stream_element_nml('namelist.jsbach')
    CALL set_stream_nml('namelist.jsbach')
    IF (theOptions%Standalone) CALL init_output
    IF (theOptions%Standalone) CALL open_output_streams
    IF (theOptions%Standalone) CALL mvstream_update_cache

    IF (theOptions%Standalone) THEN
       ! Print status of streams
       CALL message('','')
       DO i = 1, nstreams
          CALL print_stream (ostreams (i))
       ENDDO
       ! Print events for model run
       CALL print_events
    END IF

    IF (PRESENT(grid)) grid = theLand%Grid
    IF (PRESENT(domain)) domain = theLand%Domain
    IF (PRESENT(options)) options = theOptions

    IF (debug) CALL message('jsbach_init','End')

  END SUBROUTINE jsbach_init

  SUBROUTINE jsbach_restart

!LK - please, use your own timer:    USE mo_timer,            ONLY: timer_start, timer_stop, timer_restart
    USE mo_io,               ONLY: IO_read_streams

    ! Note: in coupled mode the restart streams are read from the GCM 
    !       (SUBROUTINE *iorestart* in ECHAM5, called from *control*)
    IF (.not. theOptions%Standalone) THEN !! case ECHAM + JSBACH
       ! Convert REAL netCDF representation of surface cover types to INTEGER
       theLand%Surface%cover_type = NINT(theLand%Surface%cover_type_real)

    else !! case JBACH ALONE

       !! -- Read restart files and set time
!LK - please, use your own timer:       IF (theOptions%Timer) CALL timer_start(timer_restart)
       CALL io_read_streams
!LK - please, use your own timer:       IF (theOptions%Timer) CALL timer_stop(timer_restart)

       ! Convert REAL netCDF representation of surface cover types to INTEGER
       theLand%Surface%cover_type = NINT(theLand%Surface%cover_type_real)

    END IF

  END SUBROUTINE jsbach_restart

  SUBROUTINE stepon_jsbach

    ! This subroutine is called from *stepon* (ECHAM) or *jsbalone_driver* (standalone JSBACH) at the 
    ! beginning of each time step. That means that in ECHAM it is called before first scans over
    ! gaussian latitudes.
    ! At the beginning of a run (lstart or lresume), *stepon_jsbach* is instead called from *control*
    ! (ECHAM) or *jsbach_init* (standalone JSBACH) before the JSBACH interface is called for the first
    ! time and after initial and restart files have been read.
    USE mo_land_surface,          ONLY: init_cover_fract, init_veg_ratio_max, init_cover_fract_pot, &
                                        define_masks, scale_cover_fract, &
                                        fract_small
    USE mo_cbal_landcover_change, ONLY: read_landcover_fractions, read_landuse_transitions, read_harvest
    USE mo_cbal_bethy,            ONLY: read_ndepo, read_nitrogen_deposition_fluxes, read_nfert, read_nitrogen_fertilizer_fluxes
    USE mo_jsbach,                ONLY: new_day, new_month, new_year
    USE mo_time_control,          ONLY: get_date_components, current_date, previous_date, lstart, lresume
    USE mo_jsbach_constants,      ONLY: zsigfac, zepsec, zqsncr
    USE mo_forest_management,     ONLY: read_rotation_length, read_age_map

    INTEGER :: current_year, previous_year, current_month, previous_month, current_day, previous_day
    INTEGER :: i

    ! Set logical flags new_day and new_year
    !     Comparing current_date and previous_date leads to new_day for the first time step after 00:00.
    !     It might be more correct to define new_day by comparing current_date and next_date. We stay
    !     however with the estabished version.
    CALL get_date_components (current_date, year=current_year, month=current_month, day=current_day)
    CALL get_date_components (previous_date,  year=previous_year, month=previous_month, day=previous_day)
    new_day   = current_day   /= previous_day
    new_month = current_month /= previous_month
    new_year  = current_year  /= previous_year

    ! Initialize cover fractions and glacier mask
    IF (lstart .OR. (lresume .AND. theOptions%ReadCoverFract)) THEN
       ! Initialize cover_fract
       theLand%Surface%cover_fract(:,:) = init_cover_fract(:,:)
       CALL define_masks(theLand%Surface, theLand%LctLibrary)  ! Generate masks from new cover fractions
       CALL scale_cover_fract(theLand%Domain%nland, theLand%ntiles, theLand%Surface%is_present(:,:), &
                              theLand%Surface%is_glacier(:,:), theLand%Surface%is_naturalVeg(:,:), &
                              theLand%Surface%cover_fract(:,:))
       CALL message('jsbach_interface','Using cover fractions from ini file')

       ! Initialize potential natural vegetation cover fractions
       IF (theOptions%UseDynveg .OR. theOptions%UseDisturbance .OR. theOptions%UseLanduseTransitions) THEN
          theLand%Surface%cover_fract_pot(:,:) = init_cover_fract_pot(:,:)
          CALL scale_cover_fract (theLand%Domain%nland, theLand%ntiles,  theLand%Surface%is_present(:,:), &
                                  theLand%Surface%is_glacier(:,:), theLand%Surface%is_naturalVeg(:,:), &
                                  theLand%Surface%cover_fract_pot(:,:))
          CALL message('jsbach_interface','Using potential natural vegetation from ini file')
       END IF

       ! Initialize veg_ratio_max
       theLand%Surface%veg_ratio_max(:) = init_veg_ratio_max(:)
       WHERE (ANY(theLand%Surface%is_glacier(:,:),DIM=2))
          theLand%Surface%veg_ratio_max(:) = 0._dp
       ELSEWHERE
          theLand%Surface%veg_ratio_max(:) = MAX(fract_small, theLand%Surface%veg_ratio_max(:))
       END WHERE
       CALL message('jsbach_interface','Using veg_ratio_max from ini file')

    ELSE
       CALL define_masks(theLand%Surface, theLand%LctLibrary)   ! Generate masks from restart cover fractions
    END IF

    IF (lstart .OR. lresume) THEN
      DO i=0,13
         theLand%Vegetation%lai_clim(:,:,i) = MERGE(theLand%Vegetation%lai_clim(:,:,i),0._dp,theLand%Surface%is_vegetation)
      END DO

      IF ( .NOT. lresume) THEN ! these variables are in the restart now -- because they are required when restarting with allometry
        ! Calculate lai_max and veg_fract_correction (representing vegetation gaps)
        DO i=1,theLand%ntiles
          WHERE (theLand%Surface%is_vegetation(:,i))
            theLand%Vegetation%lai_max(:,i) = theLand%LctLibrary%MaxLAI(theLand%Surface%cover_type(:,i))
            theLand%Vegetation%veg_fract_correction(:,i) = 1.0_dp - exp(-theLand%Vegetation%lai_max(:,i) &
                                                  / theLand%LctLibrary%clumpinessFactor(theLand%Surface%cover_type(:,i)))
          ELSEWHERE
            theLand%Vegetation%lai_max(:,i) = 0.0_dp
            theLand%Vegetation%veg_fract_correction(:,i) = 1.0_dp
          END WHERE
        END DO
      END IF
    END IF


    IF (lstart) THEN
       ! Roesch et al, 2002, Climate Dynamics (compare with snow_fract calculation in update_soil)
       DO i=1,theLand%ntiles
          WHERE (theLand%Surface%is_glacier(:,i))
             theLand%Soil%snow_fract(:,i) = 1.0_dp
          ELSEWHERE
             theLand%Soil%snow_fract(:,i) = zqsncr * TANH(theLand%Soil%snow(:,i) * 100._dp) *              &
                  SQRT(theLand%Soil%snow(:,i) * 1000._dp / (theLand%Soil%snow(:,i) * 1000._dp + zepsec +   &
                  zsigfac * theLand%Surface%oro_std_dev(:)))
          ENDWHERE
       END DO

    END IF

    !Read rotation length from maps 
    IF ((lstart .OR. lresume) .AND. theOptions%UseRotationMap) THEN
       CALL read_rotation_length(theLand%Grid,theLand%Domain,theLand%Surface)
    ENDIF
    
    !Read forest age from maps 
    IF (lstart .AND. theOptions%UseRotationMap) THEN
       CALL read_age_map(theLand%Grid,theLand%Domain,theLand%Surface)
    ENDIF

    !!
    !! === Do stuff at beginning of new year
    !!
    ! If at the first time step of a new year or first time step of new run ==> read in new landcover map
    !    (run might start in the middle of a year)
    IF (theOptions%UseExternalLandcoverMaps .AND. (new_year .OR. lstart)) THEN
       CALL read_landcover_fractions(current_year,theLand%Grid,theLand%Domain,theLand%Surface)
    END IF
    ! --- update Nitrogen deposition forcing
    IF (theOptions%WithNitrogen .AND. read_ndepo .AND. (lstart .OR. lresume .OR. new_month)) then
       call read_nitrogen_deposition_fluxes(current_year,current_month,theLand%Grid,theLand%Domain,theLand%Nbalance)
    end IF

    IF (theOptions%WithNitrogen .AND. read_nfert .AND. (lstart .OR. lresume .OR. new_year)) then
       call read_nitrogen_fertilizer_fluxes(current_year,theLand%Grid,theLand%Domain,theLand%Nbalance,          &
                                                     theLand%surface%is_crop(:,1:theLand%ntiles),theLand%ntiles)
    end IF
    ! === If first time step in new year read in 3x3 NHHLP landuse transitions and compute 4x4 landuse transitions
    IF (theOptions%UseLanduseTransitions .AND. (new_year .OR. lstart)) THEN
       call read_landuse_transitions(current_year,theLand%Grid,theLand%Domain)
       call read_harvest(current_year,theLand%Grid,theLand%Domain)
    END IF

  END SUBROUTINE stepon_jsbach

  SUBROUTINE stepoff_jsbach

    ! This subroutine is called from *stepon* (ECHAM) or *jsbalone_driver* (standalone JSBACH) at the
    ! end of each time step.

    USE mo_hydrology,             ONLY: hydrology_model, hydrology_restart
    USE mo_time_control,          ONLY: lstop, l_putrerun, l_puthd


    ! --- call hydrological discharge and glacier calving model

    IF (theOptions%WithHD .AND. l_puthd) THEN
       CALL hydrology_model(theLand%Domain%mask, theOptions%Standalone)
    END IF

    ! --- clean up at the end of a simulation
    IF (theOptions%WithHD .AND. (l_putrerun .OR. lstop)) THEN
       CALL hydrology_restart
    END IF

  END SUBROUTINE stepoff_jsbach

  SUBROUTINE get_dates(IO_timestep, istep)

    ! For standalone mode only.
    ! Corresponds to ECHAM5 subroutine IO_init in mo_io.f90
    ! This routine doesn't need to be called from ECHAM in coupled mode.

    USE mo_io, ONLY: IO_READ, IO_open, IO_close
    USE mo_netCDF, ONLY: FILE_INFO, NF_GLOBAL, io_get_att_int
    USE mo_time_control, ONLY: INIT_STEP, delta_time, dt_start, &
                               lresume, resume_date, start_date, inp_convert_date, write_date
    USE mo_time_conversion, ONLY : time_native, TC_set, TC_convert

    INTEGER, INTENT(out) :: IO_timestep, istep

    TYPE(FILE_INFO) :: IO_file
    INTEGER :: IO_file_id
    INTEGER :: forecast_date, verification_date   ! YYYYMMDD 
    INTEGER :: forecast_time, verification_time   ! HHMMSS
    TYPE(time_native) :: date_nat

    IF (debug) CALL message('jsbach_init_io','BEGIN')

    ! Initialize time
    IF (p_parallel_io) THEN
       IF (lresume) THEN
          ! Get time information from restart file
          IO_file%opened = .FALSE.
          CALL IO_open(trim(theOptions%RestartPrefix)//'_jsbach.nc', IO_file, IO_READ)
          IO_file_id = IO_file%file_id
          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'nstep', istep)
          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'timestep', IO_timestep)
          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'fdate', forecast_date)
          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'ftime', forecast_time)
          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'vdate', verification_date)
          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'vtime', verification_time)
          CALL IO_close(IO_file)
       ELSE
          istep = INIT_STEP
          IO_timestep = delta_time
       ENDIF
    ENDIF
       
    IF (p_parallel) THEN
       CALL p_bcast(istep, p_io)
       CALL p_bcast(IO_timestep, p_io)
       IF (lresume) THEN
          CALL p_bcast(forecast_date, p_io)
          CALL p_bcast(forecast_time, p_io)
          CALL p_bcast(verification_date, p_io)
          CALL p_bcast(verification_time, p_io)
       ENDIF
    ENDIF
       
    IF (lresume) THEN
       CALL inp_convert_date(forecast_date, forecast_time, start_date)
       CALL inp_convert_date(verification_date, verification_time, resume_date)
    ELSE
       IF (SUM(dt_start(:)) /= 0) THEN
          CALL TC_set(dt_start(1),dt_start(2),dt_start(3), &
               dt_start(4),dt_start(5),dt_start(6), date_nat)
          CALL TC_convert(date_nat, start_date)
          resume_date = start_date
       ELSE
          CALL finish('jsbach_init_io', 'Start date not set')
       ENDIF
    ENDIF

    IF (p_parallel_io) THEN
       CALL write_date(start_date, 'Start date (initial/restart): ')
       CALL write_date(resume_date, 'Resume date (initial/restart): ')
    ENDIF

    IF (debug) CALL message('get_dates','END')

  END SUBROUTINE get_dates

  SUBROUTINE update_current_call(kdim, kland, kblock, mask)

    USE mo_jsbach_grid,     ONLY: update_comm_to_echam5mods

    INTEGER,           INTENT(in) :: kdim
    INTEGER,           INTENT(in) :: kland
    INTEGER, OPTIONAL, INTENT(in) :: kblock
    LOGICAL, OPTIONAL, INTENT(in) :: mask(kdim)

!!$    INTEGER, ALLOCATABLE, TARGET, SAVE :: kindex_help(:)
    LOGICAL :: mask_help(kdim)    !! Land mask or all true's
    INTEGER,  ALLOCATABLE  :: zint1d(:)
    INTEGER,  ALLOCATABLE  :: zint2d(:,:)

    INTEGER :: i

#if defined (__SX__) && defined (_OPENMP)
    INTEGER :: tid
    tid = omp_get_thread_num()
#endif

    IF (PRESENT(kblock)) THEN
       theLand%Domain%kblock = kblock
    ELSE
       theLand%Domain%kblock = 1
       IF (kdim /= theLand%Domain%nlon * theLand%Domain%nlat .AND. kdim /= theLand%Domain%nland) THEN
          call message('jsbach_inter_1d: kdim, nlon, nlat, nland for domain: ',&
               int2string(kdim)//' '//int2string(theLand%Domain%nlon)//' '//int2string(theLand%Domain%nlat)//' '// &
                                                                    int2string(theLand%Domain%nland))
          CALL finish('jsbach_inter_1d','Dimension in call inconsistent with domain')
       ENDIF
    ENDIF

    IF (PRESENT(mask)) THEN
       mask_help = mask
    ELSE
       mask_help = .TRUE.
    ENDIF

    ! Index of current land boxes into packed vector of land boxes for each PE's domain
    nidx = kland
    ALLOCATE(kindex(kland))

    IF (ALLOCATED(zint1d)) CALL finish('jsbach_inter_1d', 'zint1d already allocated') 

    ALLOCATE(zint1d(theLand%Domain%nland)) ; 
    zint1d = (/ (i, i=1,theLand%Domain%nland) /)        ! Index of domain's land points
    IF (PRESENT(kblock)) THEN
       ALLOCATE(zint2d(theLand%Domain%ndim,theLand%Domain%nblocks))
       zint2d(:,:) = UNPACK(zint1d, MASK=theLand%Domain%mask, FIELD=0) 
                                              ! Location of all domain land points in (ndim,nblocks) grid
       kindex(1:kland) = PACK(zint2d(1:kdim,kblock), MASK=mask_help)
                                              ! Note: kdim can be smaller than grid%ndim for last block!
                                              ! Note: if no mask given, mask=.T. and kland=kdim
       DEALLOCATE(zint2d)
    ELSE   ! kblock not present, i.e. kdim = grid%ndim = grid%nlon * grid%nlat, i.e. whole domain in one call (block, nblocks=1)
       kindex(:) = zint1d(:)     ! Note: kland = grid%nland
    ENDIF
    kstart = kindex(1)
    kend = kindex(kland)
    IF (ANY(kindex == 0)) THEN
       CALL message('mo_jsbach_interface - update_current_call','nidx, kstart, kend = '// &
                                      int2string(nidx)//' '//int2string(kstart)//' '//int2string(kend))
       CALL finish('              -  ','Grid of calling programm inconsistent with JSBACH grid')
    ENDIF

    DEALLOCATE(zint1d)

    IF (test_stream) CALL init_test_2(kstart,kend, kblock , &
         theLand%Surface%cover_fract(kstart:kend,:), &
         theLand%Surface%is_present(kstart:kend,:) .AND. &
         .NOT. theLand%Surface%is_lake(kstart:kend,:), theLand%Domain)

#if defined (__SX__) && defined (_OPENMP)
    IF (tid == 0) &
#endif
         CALL update_comm_to_echam5mods(theLand%Grid, theLand%Domain)

  END SUBROUTINE update_current_call

END MODULE mo_jsbach_interface


!Local Variables:
!mode: f90
!fill-column: 100
!End:
