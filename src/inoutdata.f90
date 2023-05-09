module inoutdata
   USE nrtype
   character(len=13) :: model                                        ! Type of model
   character(len=13) :: dep_model                                    ! Depositional model
!     Constants
   real(dp), parameter ::  denatm = 1.225d0                            ! Air density
   real(dp), parameter ::  r_air = 287.d0							   ! Gas constant of air
   real(dp), parameter ::  g = 9.81d0                                  ! Gravity acceleration
   real(dp), parameter ::  df2 = 2.d0, df100 = 100.d0                     ! Reference minimum and maximum flow density
   real(dp), parameter ::  theta = 0.015d0!,pn=2.5d0                    ! Shield parameter, Rouse number at settling
   real(dp), parameter ::  dz = 0.01d0                                 ! Vertical space step for profiles
   real(dp), parameter ::  small = 1.d-15                              ! Small number
   real(dp), parameter ::  xacc = 1.d-15                               ! Accuracy level in rtflsp
   real(dp), parameter ::  kvk = 0.41d0                                ! Von Karman's constant
   real(dp), parameter ::  xmaxchi = 300.d0                          ! Maximum Chi value
   real(dp), parameter :: undefined = 9.87654321d31
   integer, parameter :: undefined_i = 987654321
   character, parameter :: undefined_c = ' '
   real(dp) :: mu, dengas                                            ! Gas viscosity and density
! input data common throughout the code
   real(dp) :: dens_ent, dm_ent                                       ! Entrained particle density and diameter
   real(dp), dimension(6) :: phi50, sorting, davgeqsph, d50mm, phi84, d84mm, phi16, d16mm           ! Median grainsize, sorting, average diameter of the equivalent sphere of particles in the median size class, median grainsize in mm
   real(dp) :: probt, zlam, zlams, c0, ks                                ! Significance level of the T test, layer thickness, sublayer thickness, reference particle concentration, substrate roughness
   logical :: distr1, distr2                                          ! Flags to decide if to read grainsize distributions of the two components from separate files
   real(dp) :: dummy1, dummy2, dummy3                                  ! Dummy variables read from the particle data files and ignored

   real(dp) :: pnsavgguess, pnsmaxguess, pnsminguess, zsfavgguess, z0minguess, z0maxguess, z0avgguess ! Starting values for profile calculations
   !FABIO: new guess variables for the new method
   real(dp) :: rhogavgguess, rhogmaxguess, rhogminguess, ztavgguess, ztmaxguess, ztminguess ! Starting values for profile calculations
   real(dp) :: rhogavg, rhogmax, rhogmin
   
   logical :: lEXISTS                                                ! Flag to check if input data files exist in the folder
! data transferred between routines
   real(dp) :: fx                                                    ! Cd/d calculated in different ways
   real(dp) :: alfa                                                  ! Two-tails t-Student significance level
   real(dp) :: usqnrm                                                ! Normalized squared shear velocity
   integer :: nfunc, nfunc1                                           ! Number of the function to be used in the routine "func" and "func1"
   real(dp) :: tcalc, ttab, nfreed, nclasst                              ! Calculated and tabulated t values, degrees of freedom and number of grainsize classes
   real(dp) :: dennrm, denmax, denmin, tauavg, taumax, taumin, ushavg, ushmax, ushmin !Medium,maximum,minimum flow density, shear stress and shear velocity
   real(dp) :: densp                                                 ! Particle density
   real(dp) :: ztot                                                  ! Total flow height
   real(dp) :: ztavg, ztmax, ztmin                                     ! Average, maximum, minimum total flow height
   real(dp) :: z0avg, z0max, z0min                                     ! Average, maximum, minimum z0 height in the Rouse equation
   real(dp) :: pnsavg, pnsmax, pnsmin                                  ! Average, maximum, minimum Rouse number
   real(dp) :: zsfavg, zsfmax, zsfmin                                  ! Average, maximum, minimum shear flow height
   real(dp) :: zttemp, z0temp, pnstemp, cgastemp, cairtemp			 ! Temporary ztot, z0, pns and gas and air concentration used to calculate the average temperature profile
   real(dp) :: p10av1, p10mx1, p10mn1, c2av1, c2max1, c2min1              ! Temporary average, maximum, minimum 10 m dynamic pressure and 2 m particle concentration
   real(dp) :: p10avg, p10max, p10min, c2avg, c2max, c2min                ! Average, maximum, minimum 10 m dynamic pressure and 2 m particle concentration
   real(dp), dimension(100) :: pzav1, pzmax1, pzmin1, czav1, czmax1, czmin1, tzav1, tzmax1, tzmin1 ! Temporary average, maximum, minimum dynamic pressure, particle concentration and temperature at user requested height
   real(dp), dimension(100) :: zdynpr, zc, zt                              ! User selected heights for dynamic pressure, particle concentration and temperature
   real(dp), dimension(100) :: pzavg, pzmax, pzmin, czavg, czmax, czmin, tzav, tzmax, tzmin    ! Average, maximum, minimum  dynamic pressure, particle concentration and temperature at user requested height
   integer :: ipr, ic, itemp                                           ! Number of user requested heights for dynamic pressure, particle concentration and flow temperature
   logical :: usr_z_dynpr, usr_z_c, usr_z_t                            ! User's choice on if to calculate dynamic pressure, particle concentration and temperature at different heights
   real(dp) :: d1m, d2m, cd1                                           ! Diameter of the first and second component (m) and average drag coefficient of the first component in the two components method
   logical :: checkzbr                                               ! Logical variable to check if the root finding routine zbrent converges
   real(dp) :: z0, zshr                                               ! Reference level thickness in the Rouse equation and shear flow thickness (Newton routine)
   real(dp) :: pns, pnstmp                                            ! Rouse number and temporary Rouse number (Newton routine)
   real(dp) :: den                                                   ! Flow density used (Newton routine)
   integer :: nnewt                                                 ! Control for the system of equations to be solved with the Newton's method
   logical :: usr_pcx_sol                                            ! User's choice on if to calculate the function values at any desired percentile
   real(dp), dimension(100) :: pcx                                   ! User requested percentiles
   real(dp) :: px                                                    ! User requested percentiles
   real(dp) :: sigsim, musim, mm                                       ! Standard deviation, median and simmetrization exponent of the symmetrized distribution
   real(dp) :: mudstr, mxdstr, mndstr                                  ! Median, maximum and minimum value of the original distribution results
   logical :: calc_t_mix
   character(len=10) :: input_weight                                  ! If "MASS" the weights in input are read as mass (gr), if WT as weight fraction (wt%); in the latter case, wt% is referred to the WHOLE samples (all the considered components)
! Variables for grainsize
   logical, dimension(6) :: dotestchi                                ! Flag to decide if to perform the Chi test for each component
   real(dp) :: fxdistr
   real(dp) :: nfreedchi, probchi, sensgrchi                           ! Degrees of freedom for the Chi test, Significance level of the Chi test
   real(dp) :: xnmax                                                 ! Real variable for calculating number of grainsize classes
   real(dp), dimension(6) :: siglevchi, senschi                       ! Significance level of the Chi test, sensitivity in the merging of grainsize classes weightfractions
   real(dp), dimension(6) :: dphi, phimin, phimax                      ! Phi stepsize, minimum and maximum dimension (phi) of the grainsize distribution
   integer, dimension(6) :: nclass                                   ! Number of classes in the grainsize distributions
   real(dp), dimension(6, 30) :: weight                               ! Weight of the grainsize classes(gr)
   real(dp), dimension(6, 0:30) :: rhos                               ! Density of grainsize classes (0 refer to the median grainsize)
! Variables for flow temperature
   real(dp) :: t_gas, t_particles, t_air							! Temperature of the magmatic gas, solid particles and atmosphere
   real(dp) :: rho_gas, rho_air, rho_particles						! Density of magmatic gas, air and solid particles
   real(dp) :: cp_air, cp_particles, cp_gas							! Specific heat at constant at constant pressure of atmosphere, particles and magmatic gas
   real(dp) :: r_gas												! Gas constant of magmatic gases
   real(dp) :: p_air												! Atmospheric pressure
   real(dp) :: t_mix_avg, t_mix_max, t_mix_min						! Mixture temperature (average, maximum and minimum)
   real(dp) :: c_gas_avg, c_gas_max, c_gas_min						! Gas concentration (average, maximum and minimum)
   real(dp) :: c_air_avg, c_air_max, c_air_min						! Air concentration (average, maximum and minimum)
   real(dp) :: r_mix                                          
! Deposition
   logical :: deprates, only_deprates                                 ! If true, the deposition rates and times are calculated for each solution; if the second is true, the code will only run the deposition rate model
   logical :: pn_cut                                                 ! If true, the code considers also classes for which Pn>5
   real(dp) :: wtot_sample, psreadtot                      ! Total weight and weight fraction read from the grainsize distributions in input.dat, minimum reference level in the Rouse equation
   real(dp) :: zlam_massive                                          ! Thickness of the fine massive layer
   integer :: ncomp                                                  ! Number of components in the deposit
   integer, dimension(6) :: nclass_dep                               ! To store temporary number of classes
   real(dp), dimension(6, 30) :: psread, dread, rhosread                ! Read weight fractions, read grainsize dimensions, read densities
   real(dp), dimension(6, 0:30) :: shpar1read, shpar2read, shpar3read     ! To store all shape parameters (worst case, three for some laws)
   character(len=10), dimension(6) :: rho_custom                     ! CONSTANT=constant size-independent value; VARIABLE=size dependent density (to be read from rhos(i,j)
   character(len=5) :: fluid                                          ! If water, Yalin law is used for the Shields parameter, otherwise the Miller law
   logical, dimension(6) :: merge_classes                            ! Flag to decide if to merge grainsize classes of each component
   real(dp), dimension(6) :: sensmerge                               ! Sensitivity of grainsize classes merging of each component
   integer :: n_solutions, kmax                                       ! Number of solutions for DEPRATES if only_deprates=.T., final number of considered solutions                     !
   real(dp), dimension(5) :: rho_flow, ztot_flow, ush_flow, pns_flow   !Flow density, thickness, shear velocity, thickness, suspension Rouse number,
   real(dp) :: rhoflow, ushearflow, ztotflow, ctotflow, pnsflow          ! flow density,shear velocity, shear stress, total flow thickness, flow concentration, population Rouse number
   real(dp) :: z0_x, pntemp                                           ! Temporary z0 and Rouse number
   real(dp) :: wt                                                    ! Terminal velocity calculated in the routine cdmodel
   real(dp), dimension(5) :: z0_final, zlam_final, zlam_susp, zlam_wash ! Final values of z0 and layer thickness after reworking
   real(dp), dimension(5) :: rtot_susp, rtot_wash, rtot_massive, ar, qtot, srw, sqratio ! Sedimentation rate from turbulent suspension, wash load and wash load forming the fine massive layer

   real(dp), dimension(5) :: ctot_flow, ctot_dep, ctot_susp, ctot_wash, ctot_massive        ! Total particle concentration in the turbulent suspension, at deposition and that remains in suspension, Total particle volumetric concentration from turbulent suspension, wash load
   real(dp), dimension(5) :: tdep_susp, tdep_wash, tdep_massive        ! Deposition time from turbulent suspension, wash load and wash load forming the fine massive layer
   real(dp), dimension(5) :: rtot, ctot, tdep                          ! Total deposition rate, total particle volumetric concentration, deposition time
   real(dp) :: rtot_max, rtot_min, rtot_avg, tdep_max, tdep_min, tdep_avg ! For probability function
   real(dp) :: slope_ground                                          ! Slope. If underfined, it is recalculated
   real(dp) :: dep_median, rhos_median                               ! Median of the total deposit (for bedload transportation calculation), density of the median
! Shape parameters for the different Cd laws
   integer :: ishape, jshape                                          ! Indexes for tansferring shape parameters data to cdmodel
   real(dp), dimension(6, 0:30) :: shpar1, shpar2, shpar3
!      logical, dimension(6,0:30) :: shpar4
   logical, dimension(6) :: shpar4
   real(dp) :: p11, p12, p13, p14                                       ! Different shape parameters
   logical :: p15                                                    ! isometric?
   real(dp), dimension(6, 0:30) :: sphericity, longspher, crossspher, circularity, voleqsphd,&
  &circeqard, corey, flatratio, shapefact, fractdim
!      logical, dimension(6,0:30) :: isometric
   logical, dimension(6) :: isometric
   logical, dimension(6) :: fractal
   character(len=10), dimension(6) :: cdlaw, rholaw

!******PARAMETERS OF THE DIOGUARDI & MELE DRAG LAW*****
   real(dp), parameter :: acd = 1.1883d0, bcd = 0.4826d0                  ! Fitting parameters
   real(dp), parameter :: exp1 = -0.23d0, exp2 = 0.05d0                   ! Exponent of the shape factor function
!******PARAMETERS OF THE DIOGUARDI et al. (2018) LAW*****
   real(dp), parameter :: exp_a = 0.25d0, exp_b = 0.08d0, exp_c = -5.05d0

   real(dp), parameter :: restart = 100.d0                             ! Starting Re number for the iterative Cd calculations
   real(dp), parameter :: tolcd = 1.d-10                                ! Tolerance in the iterative Cd calculation
!******PARAMETERS OF THE DIOGUARDI et al. (2017) LAW*****
   real(dp), parameter :: a_fd = 0.3492d0, b_fd = 1.3358d0                      ! Fitting parameters
   real(dp), parameter :: a_sph = 0.559d0, b_sph = 0.5134d0                      ! Fitting parameters
   real(dp), parameter :: exp1_fd = 1.62d0, exp2_fd = -0.13d0                   ! Exponents
   real(dp), parameter :: exp1_sph = 4.18d0, exp2_sph = -0.2d0                   ! Exponents

end module inoutdata
