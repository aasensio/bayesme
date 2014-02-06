!********************************************************
! Global variables
!********************************************************
module vars
implicit none

	real(kind=8), parameter :: PK = 1.3806503d-16, UMA = 1.66053873d-24, PC = 2.99792458d10
	real(kind=8), parameter :: PH = 6.62606876d-27, PHK = PH / PK, PHC = 2.d0 * PH / PC**2
	real(kind=8), parameter :: PME = 9.10938188d-28, PE = 4.8032d-10, PI = 3.14159265359d0
	real(kind=8), parameter :: PHC2 = 2.d0 * PH * PC**2, OPA = PI * PE**2 / (PME * PC)
	real(kind=8), parameter :: SQRTPI = 1.77245385091d0, EARTH_SUN_DISTANCE = 1.495979d13
	real(kind=8), parameter :: LARMOR = PE/(4.d0*PI*PME*PC), PMP = 1.67262158d-24
	real(kind=8), parameter :: NUCLEAR_MAGNETON = PE*PH/(2.d0*PI*PMP)
	real(kind=8), parameter :: BOHR_MAGNETON = PE*PH/(2.d0*PI*PME)
	real(kind=8), parameter :: PH2C3 = 2.d0 * PH**2 * PC**3, PERTURBATION = 0.005d0
		
	real(kind=8), allocatable :: zeeman_voigt(:,:), zeeman_faraday(:,:)
	
	integer :: n_lineas
	character(len=120) :: fich_malla, fich_lineas, fich_chain, file_stray_light
	
	type stokes_type
		integer :: nlambda
		real(kind=8), pointer :: stokes(:,:), lambda(:), sigma(:,:), variance_total(:,:)
		real(kind=8), pointer :: sigma_minimum(:)
		real(kind=8) :: step_lambda
	end type stokes_type
	
	type modelo_type
		logical :: stray_light_component
		integer :: stray_light_nlambda, npix_stray
		real(kind=8) :: Bfield, theta, chi, ff, vmac, damping, beta, macrot, mu
		real(kind=8) :: doppler, filling_factor
		real(kind=8), dimension(2) :: Bfield_range, theta_range, chi_range, ff_range
		real(kind=8), dimension(2) :: vmac_range, damping_range, beta_range, macrot_range
		real(kind=8), dimension(2) :: doppler_range, filling_factor_range
		real(kind=8), pointer :: kl_range(:,:) !,kl(:)
		real(kind=8), pointer :: straylight(:,:), straylight_cov(:,:)
		real(kind=8) :: kl(10)
	end type modelo_type
	
	type inversion_type
		logical :: free_filling_factors
		integer :: n_cycles, n_params, nparams_invert
		integer, pointer, dimension(:) :: which_to_invert, filling_factor_position,&
			dimension_space_B, which_parameter, which_component, which_line		
		real(kind=8), pointer, dimension(:,:) :: range, range_modified
		real(kind=8) :: log_prior_volume_inv, chi2_min
	end type inversion_type
	
	type statistics
		integer :: nparams, accepted_models
		integer, pointer :: is_accepted(:)
		real(kind=8), pointer :: covariance(:,:), parameters_proposed(:), &
			parameters_last(:), mean(:), parameters_MAP(:)
		real(kind=8), pointer :: eigenvectors(:,:), eigenvalues(:)
		real(kind=8) :: alpha
		logical :: calculated
	end type statistics
	
	type line_type
		character(len=2) :: theory
		real(kind=8) :: lambda_init, lambda_end, lambda_step, lambda
		real(kind=8) :: wave0, Jup, Jlow, gup, glow
		real(kind=8) :: Lu, Su, Au, Bu, Ll, Sl, Al, Bl, I
	end type line_type

	type posterior_function_type
		integer :: n_variables
		character(len=10), pointer :: variables(:)
		character(len=120) :: posterior_function
	end type posterior_function_type
	
	type chain_analysis_type
		integer :: nlength, nparams, n_posterior_combinations
		real(kind=8) :: evidence, avg_lnL, Dkl
		real(kind=8), pointer :: chain(:,:)
		type(posterior_function_type), pointer :: posterior_fun(:)
	end type chain_analysis_type

	type instrumental_profile_type
		integer :: nlambda_original, nlambda_forfilter, n2
		real(kind=8), pointer :: lambda_original(:), transmission_original(:)
		real(kind=8), pointer :: lambda_forfilter(:), transmission(:), data(:), respns(:)
		complex(kind=8), pointer :: transmission_fft(:)
		real(kind=8) :: step_lambda
	end type instrumental_profile_type
	
	type(stokes_type) :: Observation, Emergent
	type(stokes_type) :: Stokes_Syn
	
	type(modelo_type), pointer :: model(:)
	
	type(inversion_type) :: inversion
	
	type(line_type), pointer :: linea(:)
	
	type(chain_analysis_type) :: chain_analysis

	type(instrumental_profile_type) :: instrumental_profile
	
	integer :: what_to_do, number_of_components, idum, chain_maxlength, uniform_proposal_steps
	integer :: start_or_restart, uniform_B, uniform_thB, use_stray_light
	real(kind=8) :: desired_acceptance
	real(kind=8) :: stokes_weights(4)
	logical :: using_instrumental_profile

	character(len=10) :: parameters_name(8) = (/ 'B         ','thetaB    ','chiB      ',&
		'vdopp     ','v_mac     ','damp      ','beta      ','eta0      '/)
end module vars
