module markov_chain_montecarlo
use vars
use atomic_functions
use chain_analysis_module
! use l_bfgs_b
implicit none

	type(statistics) :: stat_nonlinear
	logical :: physical_nonlinear
	type(modelo_type), allocatable :: proposed_nonlinear(:)
	
contains

!------------------------------------------------------------
! Diagonalize the covariance matrix S
!  D=transpose(T)##S##T , where D is diagonal
!  T is saved in stat%eigenvectors
!------------------------------------------------------------
	subroutine diagonalize_covariance(stat)
	type(statistics) :: stat
	integer :: nrot, i
	real(kind=8) :: A(stat%nparams,stat%nparams)
		
		A = stat%covariance
!		print *, A
!		call jacobi(A, stat%nparams, stat%eigenvalues, stat%eigenvectors, nrot)
		
!		print *
!		print *, stat%eigenvalues
!		pause
		if (stat%calculated) then
			!print *, A
			!pause
			!call cholesky(A, stat%nparams, stat%eigenvalues, stat%eigenvectors)		
		
			!print *, stat%eigenvectors
			!print *, matmul(stat%eigenvectors,transpose(stat%eigenvectors))
			!pause
		endif
				
!		print *, stat%covariance
!		print *, 'eigen : ', stat%eigenvalues		
				
	end subroutine diagonalize_covariance
	
!------------------------------------------------------------
! Compute the statistics of the likelihood function
! It uses the new value of the parameters to update the mean and the covariance matrix
!------------------------------------------------------------
	subroutine compute_statistics(stat,inv,n,reset)	
	type(statistics) :: stat
	type(inversion_type) :: inv
	integer :: n, i, j
	real(kind=8) :: mold(stat%nparams)
	logical :: reset

! The first call, reset the mean and the covariance
		if (reset) then
			stat%mean = stat%parameters_last
			stat%covariance = 0.d0
			return
		endif
		
! Update the mean
		mold = stat%mean
		stat%mean = mold + (stat%parameters_last - mold) / (n+1.d0)
		
! Update the covariance matrix
		do i = 1, stat%nparams
			do j = 1, stat%nparams
				stat%covariance(i,j) = (n-1.d0)/n * stat%covariance(i,j) + &
					(stat%parameters_last(i)-mold(i))*(stat%parameters_last(j)-mold(j)) / &
					(n+1.d0)**2 + &
					(stat%parameters_last(i)-stat%mean(i))*&
					(stat%parameters_last(j)-stat%mean(j)) / (n*1.d0)
			enddo
		enddo
				
	end subroutine compute_statistics
	
!------------------------------------------------------------
! Verify if the model has physical meaning
!------------------------------------------------------------
	function test_physical(model)
	type(modelo_type) :: model(:)
	logical :: test_physical
	integer :: i, j
	real(kind=8) :: ff_total
	
		test_physical = .TRUE.
		ff_total = 0.d0
						
		do i = 1, number_of_components

			ff_total = ff_total + model(i)%filling_factor
									
			if (model(i)%filling_factor < model(i)%filling_factor_range(1) .or. &
					model(i)%filling_factor > model(i)%filling_factor_range(2)) then				
				test_physical = .FALSE.								
				return
			endif
						
			if (model(i)%bfield < model(i)%bfield_range(1) .or. &
					model(i)%bfield > model(i)%bfield_range(2)) then
				test_physical = .FALSE.					
				return
			endif
			if (model(i)%doppler < model(i)%doppler_range(1) .or. &
					model(i)%doppler > model(i)%doppler_range(2)) then
				test_physical = .FALSE.
				return
			endif
			if (model(i)%vmac < model(i)%vmac_range(1) .or. &
					model(i)%vmac > model(i)%vmac_range(2)) then
				test_physical = .FALSE.
				return
			endif
						
			if (model(i)%theta < model(i)%theta_range(1) .or. &
					model(i)%theta > model(i)%theta_range(2)) then
				test_physical = .FALSE.		
				return
			endif
			if (model(i)%chi < model(i)%chi_range(1) .or. &
					model(i)%chi > model(i)%chi_range(2)) then
				test_physical = .FALSE.
				return
			endif
			if (model(i)%damping < model(i)%damping_range(1) .or. &
					model(i)%damping > model(i)%damping_range(2)) then
				test_physical = .FALSE.
				return
			endif
			if (model(i)%beta < model(i)%beta_range(1) .or. &
					model(i)%beta > model(i)%beta_range(2)) then
				test_physical = .FALSE.
				return
			endif
			do j = 1, n_lineas
				if (model(i)%kl(j) < model(i)%kl_range(1,j) .or. &
						model(i)%kl(j) > model(i)%kl_range(2,j)) then
					test_physical = .FALSE.
					return
				endif
			enddo
		enddo
		
		if (abs(ff_total-1.d0) > 1.d-4) then
			test_physical = .FALSE.
		endif
		
	end function test_physical
		
!------------------------------------------------------------
! Locate the position in the parameter's list of the filling factors
! This is necessary for later fixing the sum of them to 1
!------------------------------------------------------------
	subroutine locate_filling_factors(inv)
	type(inversion_type) :: inv
	integer :: i, j, k
		
		allocate(inv%filling_factor_position(number_of_components))
		
		j = 1
		do i = 1, inv%nparams_invert
			k = mod((inv%which_to_invert(i)-1),(8+n_lineas)) + 1			
			if (k == 8+n_lineas) then
				inv%filling_factor_position(j) = i
				j = j + 1
			endif
		enddo		
												
	end subroutine locate_filling_factors
	
!------------------------------------------------------------
! Adjust the filling factors so that there sum is equal to 1 if the number of components
! is at least 2
!  It adjusts the filling factor of the last component so that the sum is
!  equal to 1. If the correction gives a negative filling factor, the model
!  will be rejected in the test_physical subroutine
!------------------------------------------------------------	
	subroutine adjust_filling_factors(stat,inv)
	type(inversion_type) :: inv
	type(statistics) :: stat	
	integer :: i, j
	real(kind=8) :: ff_total
		
		if (number_of_components > 1) then
			ff_total = 0.d0
			do i = 1, number_of_components
				j = inv%filling_factor_position(i)
				ff_total = ff_total + stat%parameters_proposed(j)
			enddo
			
			if (ff_total /= 1.d0) then
				j = inv%filling_factor_position(number_of_components)
				stat%parameters_proposed(j) = 1.d0 - ff_total + stat%parameters_proposed(j)
			endif
		endif			
												
	end subroutine adjust_filling_factors

!------------------------------------------------------------
! Fill a model with a new set of parameters
!------------------------------------------------------------
	subroutine fill_model(inv,stat,proposed)
	type(inversion_type) :: inv
	type(statistics) :: stat
	type(modelo_type) :: proposed(:)
	integer :: i, k, which_component, which_line
										
		do i = 1, inv%nparams_invert			
			k = inv%which_parameter(i)
			which_component = inv%which_component(i)
			if (k == 1) then				
				proposed(which_component)%bfield = stat%parameters_proposed(i)
    			cycle
			endif
			if (k == 2) then
				proposed(which_component)%theta = stat%parameters_proposed(i)
    			cycle
			endif
			if (k == 3) then
				proposed(which_component)%chi = stat%parameters_proposed(i)
    			cycle
			endif
			if (k == 4) then
				proposed(which_component)%doppler = stat%parameters_proposed(i)
    			cycle
			endif
			if (k == 5) then
				proposed(which_component)%vmac = stat%parameters_proposed(i)
    			cycle
			endif
			if (k == 6) then
				proposed(which_component)%damping = stat%parameters_proposed(i)
    			cycle
			endif
			if (k == 7) then
				proposed(which_component)%beta = stat%parameters_proposed(i)
    			cycle
			endif
			if (k == 8) then
				which_line = inv%which_line(i)
				proposed(which_component)%kl(which_line) = stat%parameters_proposed(i)
    			cycle
			endif
			if (k == 9) then
				proposed(which_component)%filling_factor = stat%parameters_proposed(i)
    			cycle
			endif	
		enddo

 	end subroutine fill_model
! 	
! !------------------------------------------------------------
! ! Fill a model with a new set of parameters
! ! B is proposed as B^3
! ! theta is proposed as cos(theta)
! !------------------------------------------------------------
! 	subroutine fill_model_uniformvector(inv,stat,proposed)
! 	type(inversion_type) :: inv
! 	type(statistics) :: stat
! 	type(modelo_type) :: proposed(:)
! 	integer :: i, k, which_component, which_line
! 	real(kind=8) :: expon
! 					
! ! Adjust the last filling factor so that they add up to 1		
! 		if (inv%free_filling_factors) call adjust_filling_factors(stat,inv)
! 						
! 		do i = 1, inv%nparams_invert
! 			k = mod((inv%which_to_invert(i)-1),(8+n_lineas)) + 1
! 			which_component = (inv%which_to_invert(i)-1) / (8+n_lineas) + 1
! 			which_line = k - 8			
! 			if (k == 1) then
! 				expon = inv%dimension_space_B(which_component)					
! 				proposed(which_component)%bfield = (stat%parameters_proposed(i))**(1.d0/expon)
!     			cycle
! 			endif
! 			if (k == 2) then
! ! Uniform sampling in cos(thetaB) or in thetaB
! 				if (uniform_thB == 1) then
! 					proposed(which_component)%theta = &
! 						acos(stat%parameters_proposed(i))*180.d0/PI
! 				else
! 					proposed(which_component)%theta = stat%parameters_proposed(i)
! 				endif
!     			cycle
! 			endif
! 			if (k == 3) then
! 				proposed(which_component)%chi = stat%parameters_proposed(i)
!     			cycle
! 			endif
! 			if (k == 4) then
! 				proposed(which_component)%doppler = stat%parameters_proposed(i)
!     			cycle
! 			endif
! 			if (k == 5) then
! 				proposed(which_component)%vmac = stat%parameters_proposed(i)
!     			cycle
! 			endif
! 			if (k == 6) then
! 				proposed(which_component)%damping = stat%parameters_proposed(i)
!     			cycle
! 			endif
! 			if (k == 7) then
! 				proposed(which_component)%beta = stat%parameters_proposed(i)
!     			cycle
! 			endif
! 			if (k >= 8 .and. k < 8+n_lineas) then
! 				proposed(which_component)%kl(which_line) = stat%parameters_proposed(i)
!     			cycle
! 			endif
! 			if (k == 8+n_lineas) then
! 				proposed(which_component)%filling_factor = stat%parameters_proposed(i)
!     			cycle
! 			endif	
! 		enddo		
! 		
! 	end subroutine fill_model_uniformvector

!------------------------------------------------------------
! Fill a model with a new set of parameters
! B is proposed as B^3
! theta is proposed as cos(theta)
!------------------------------------------------------------
	subroutine fill_model_uniformvector_new(inv,stat,proposed)
	type(inversion_type) :: inv
	type(statistics) :: stat
	type(modelo_type) :: proposed(:)
	integer :: i, k, which_component, which_line
	real(kind=8) :: expon
											
		do i = 1, inv%nparams_invert			
			k = inv%which_parameter(i)
			which_component = inv%which_component(i)
			if (k == 1) then
				expon = inv%dimension_space_B(which_component)					
				proposed(which_component)%bfield = (stat%parameters_proposed(i))**(1.d0/expon)
    			cycle
			endif
			if (k == 2) then
! Uniform sampling in cos(thetaB) or in thetaB
				if (uniform_thB == 1) then
					proposed(which_component)%theta = &
						acos(stat%parameters_proposed(i))*180.d0/PI
				else
					proposed(which_component)%theta = stat%parameters_proposed(i)
				endif
    			cycle
			endif
			if (k == 3) then
				proposed(which_component)%chi = stat%parameters_proposed(i)
    			cycle
			endif
			if (k == 4) then
				proposed(which_component)%doppler = stat%parameters_proposed(i)
    			cycle
			endif
			if (k == 5) then
				proposed(which_component)%vmac = stat%parameters_proposed(i)
    			cycle
			endif
			if (k == 6) then
				proposed(which_component)%damping = stat%parameters_proposed(i)
    			cycle
			endif
			if (k == 7) then
				proposed(which_component)%beta = stat%parameters_proposed(i)
    			cycle
			endif
			if (k == 8) then
				which_line = inv%which_line(i)
				proposed(which_component)%kl(which_line) = stat%parameters_proposed(i)
    			cycle
			endif
			if (k == 9) then
				proposed(which_component)%filling_factor = stat%parameters_proposed(i)
    			cycle
			endif	
		enddo
		
	end subroutine fill_model_uniformvector_new

			
!------------------------------------------------------------
! Calculate the likelihood, i.e., exp(-1/2*chi^2) function
!------------------------------------------------------------
	function compute_likelihood(model,linea,in_observation,in_inversion,in_synthetic,physical,chi2,reduc)
	type(stokes_type) :: in_observation, in_synthetic
	type(modelo_type) :: model(:)
	type(line_type) :: linea(:)
	type(inversion_type) :: in_inversion
	real(kind=8) :: compute_likelihood, weight, factor, sigma_min, chi2, reduc, stray_light_contamination
	real(kind=8) :: normalization_loglike
	integer :: nf, i, j
	logical :: physical
		
! Verify if the proposed model has physical sense		
 		physical = test_physical(model)

 		compute_likelihood = -1.d100
 		 		
		if (physical) then
			
			call synthesize(model,linea,in_observation,in_synthetic)
						
			compute_likelihood = 0.d0

! Compute the chisq using the appropriate weights for each cycle
			do i = 1, 4
				if (stokes_weights(i) /= 0.d0) then

! If no stray-light obtained from the observations is used, calculate
! standard chi^2 with noise variance
					if (use_stray_light == 0) then
					
						compute_likelihood = compute_likelihood + stokes_weights(i) * &
							sum(((in_observation%stokes(i,:)-in_synthetic%stokes(i,:))**2) / &
								(reduc*in_observation%sigma(i,:))**2)
					else
					
! If the stray-light contamination is obtained from the observations
! modify accordingly the chi^2 using the correct variance that includes
! the covariance between the profile and the stray light
						stray_light_contamination = model(use_stray_light)%filling_factor

						in_observation%variance_total(i,:) = (1.d0 + stray_light_contamination**2 / model(use_stray_light)%npix_stray) * in_observation%sigma(i,:)**2
						in_observation%variance_total(i,:) = in_observation%variance_total(i,:) - &
							2.d0 * stray_light_contamination * model(use_stray_light)%straylight_cov(i,:)

						compute_likelihood = compute_likelihood + stokes_weights(i) * &
							sum(((in_observation%stokes(i,:)-in_synthetic%stokes(i,:))**2) / in_observation%variance_total(i,:))
								
					endif

				endif
				
				sigma_min = sqrt(maxval((in_observation%stokes(i,:)-&
					in_synthetic%stokes(i,:))**2))
				if (in_observation%sigma_minimum(i) > sigma_min) then
					in_observation%sigma_minimum(i) = sigma_min
				endif
			enddo

! Normalization factors
			normalization_loglike = 0.5d0 * 4.d0 * in_observation%nlambda * log(2.d0*PI) + sum(log(in_observation%sigma))

! Include normalization of the likelihood and normalization of the priors
			chi2 = compute_likelihood
			compute_likelihood = -normalization_loglike - 0.5d0 * compute_likelihood + in_inversion%log_prior_volume_inv

! If the covariance is negative at some point, return zero likelihood
			if (minval(in_observation%variance_total) < 0) then
				chi2 = 1.d100
				physical = .FALSE.
				compute_likelihood = -1.d100
			endif

! Return likelihood=1. This is used for estimating the priors			
			if (start_or_restart == 2) compute_likelihood = 1.d0
		endif
				
	end function compute_likelihood
	
!------------------------------------------------------------
! Proposed a new set of parameters
!------------------------------------------------------------
	subroutine proposal(inv,stat,ind,uniform)
	type(inversion_type) :: inv
	type(statistics) :: stat
	real(kind=8) :: prop, prop_last
	real(kind=8), allocatable :: parameters_last(:)
	integer :: i, j, k, which_component, which_line, ind, uniform
			
	
		stat%parameters_proposed = mrandomn(idum,stat%parameters_last,stat%alpha*stat%covariance)
!  		print *, stat%parameters_proposed
!  		pause
		
		
! Propose changes to the model		
! 		do i = 1, inv%nparams_invert
! 		
! ! Propose with gaussians with a given covariance
! 			prop = stat%alpha * sqrt(stat%covariance(i,i)) * randomn(idum) + &
! 				stat%parameters_last(i)
! 										
! 			stat%parameters_proposed(i) = prop
! 		enddo
		
	end subroutine proposal	
	
!------------------------------------------------------------
! Print output in screen
!------------------------------------------------------------
	subroutine print_output(stat,inv,format_output)
	type(inversion_type) :: inv
	type(statistics) :: stat
	integer :: i, j, k, which_component, which_line, accepted
	real(kind=8), allocatable :: transformed_mean(:), transformed_sigma(:)
	real(kind=8) :: value, value_sigma, expon
	character(len=40) :: format_output
			
		allocate(transformed_mean(inv%nparams_invert))
		allocate(transformed_sigma(inv%nparams_invert))
	
		do i = 1, inv%nparams_invert
			k = mod((inv%which_to_invert(i)-1),(8+n_lineas)) + 1
			which_component = (inv%which_to_invert(i)-1) / (8+n_lineas) + 1
			which_line = k - 8
			value = stat%mean(i)
			value_sigma = sqrt(stat%covariance(i,i))
			select case(k)				
				case(2)					
					expon = inv%dimension_space_B(which_component)
					value = value**(1.d0/expon)
					if (expon > 1) then
						value_sigma = value_sigma / (expon*value**(expon-1.d0))
					endif						
				case(3)
! Uniform sampling in cos(thetaB) or in thetaB
					if (uniform_thB == 1) then
						value = acos(value)*180.d0/PI
						value_sigma = value_sigma / sin(value*PI/180.d0) * 180.d0 / PI
					else
						value = value
						value_sigma = value_sigma
					endif
			end select			
			transformed_mean(i) = value
			transformed_sigma(i) = value_sigma
		enddo		
				
		write(*,FMT=format_output) 'Avg:   ', (transformed_mean(j),j=1,stat%nparams)		
		write(*,FMT=format_output) 'Sigma: ', (transformed_sigma(j),j=1,stat%nparams)
						
		deallocate(transformed_mean)
		deallocate(transformed_sigma)
		
					
	end subroutine print_output
	
!------------------------------------------------------------
! Print output in screen
!------------------------------------------------------------
	subroutine write_step_chain(stat,inv,flast,chi2last,acceptance_rate)
	type(inversion_type) :: inv
	type(statistics) :: stat
	integer :: i, j, k, which_component, which_line, accepted
	real(kind=8), allocatable :: transformed_par(:)
	real(kind=8) :: value, value_sigma, flast, acceptance_rate, expon, chi2last
			
		allocate(transformed_par(inv%nparams_invert))
		
		do i = 1, inv%nparams_invert
			k = mod((inv%which_to_invert(i)-1),(8+n_lineas)) + 1
			which_component = (inv%which_to_invert(i)-1) / (8+n_lineas) + 1
			which_line = k - 8
			value = stat%parameters_last(i)		
			select case(k)				
				case(2)
					expon = inv%dimension_space_B(which_component)
					value = value**(1.d0/expon)
				case(3)
! Uniform sampling in cos(thetaB) or in thetaB
					if (uniform_thB == 1) then
						value = acos(value)*180.d0/PI
					else
						value = value
					endif						
			end select			
			transformed_par(i) = value
		enddo		
		
		write(12) transformed_par, flast, chi2last, acceptance_rate
		
		deallocate(transformed_par)		
					
	end subroutine write_step_chain

!------------------------------------------------------------
! Do initialization
!------------------------------------------------------------
! 	subroutine lbfgs_init(inv, Observation)
! 	type(inversion_type) :: inv
! 	type(stokes_type) :: observation
! 	integer :: nmax
! 	character (len=60) :: task
! 	integer            :: n, m, iprint, i
! 	integer, allocatable :: nbd(:)
! 	real(kind=8)       :: f, factr, pgtol, fnew, delta, chi2
! 	real(kind=8), allocatable :: x(:), l(:), u(:), g(:)
! 	type(statistics) :: stat
! 	logical :: physical
! 	type(modelo_type), allocatable :: proposed(:)
! 
! 		nmax = inv%nparams_invert
! 
! 		stat%nparams = inv%nparams_invert
! 		allocate(stat%parameters_proposed(stat%nparams))
! 		allocate(proposed(number_of_components))
! 		proposed = model
! 
! 		allocate(x(nmax),l(nmax),u(nmax),g(nmax),nbd(nmax))
! 
! 		iprint = 101
! 		factr = 1.d7
! 		pgtol = 1.d-5
! 
! 		n = nmax
! 		m = 5 + nmax / 5
! 
! 		do i = 1, n
! 			nbd(i) = 2			
! 			l(i) = inv%range(i,1)
! 			u(i) = inv%range(i,2)
! 			x(i) = 0.5d0*(l(i) + u(i))
! 		enddo
! 
! 		task = 'START'
! 
! 		delta = 1.d-5
! 
! ! Call L-BFGS-B
! 		do while (task /= 'STOP')
! 			CALL mainlb(n, m, x, l, u, nbd, f, g, factr, pgtol, task, iprint)
! 
! 			print *, 'TASK -> ', task
! 			pause
! 
! !        The minimization routine has returned to request the
! !        function f and gradient g values at the current x.
! 			if (task(1:2) == 'FG') then
! 			
! 				stat%parameters_proposed = x
! 				call fill_model_uniformvector(inversion,stat,proposed)
! 
! 				f = compute_likelihood(proposed,linea,Observation,Stokes_Syn(1),physical,chi2,1.d0)
! 				
! 				do i = 1, n
! 					stat%parameters_proposed = x
! 
! 					stat%parameters_proposed(i) = stat%parameters_proposed(i) + delta
! 					call fill_model_uniformvector(inversion,stat,proposed)
! 
! 					print *, '--------B--------', proposed(1)%bfield
! 					print *, '--------thB--------', proposed(1)%theta
! 					print *, '--------dop--------', proposed(1)%doppler
! 		
! 					fnew = compute_likelihood(proposed,linea,Observation,Stokes_Syn(1),physical,chi2,1.d0)
! 
! 					g(i) = (fnew-f) / delta
! 					print *, f, fnew
! 				enddo
! 
! 				f = -f
! 				g = -g
! ! 				print *, x
! ! 				print *, f
! ! 				print *, g
! ! 				pause
! 			endif
! 
! 		enddo
! 
! 		deallocate(stat%parameters_proposed)
! 		deallocate(proposed)
! 
! 		stop
! 
! 	end subroutine lbfgs_init
! 
! !------------------------------------------------------------
! ! Do initialization
! !------------------------------------------------------------
! 	subroutine funct(n, xc, fc, iuser, ruser)
! 	integer :: n, iuser(:)
! 	real(kind=8) :: fc, xc(n), ruser(:), chi2
! 	logical :: physical
! 
! 		stat_nonlinear%parameters_proposed(1:n) = xc
! ! 		xc(1) = 1000.d0
! ! 		xc(2) = 40.d0
! ! 		xc(3) = 0.2d-1
! 		call fill_model_uniformvector(inversion,stat_nonlinear,proposed_nonlinear)
! 
! 		fc = compute_likelihood(proposed_nonlinear,linea,Observation,Stokes_Syn(1),physical,chi2,1.d0)
! 
! 		fc = chi2 / (4.d0*Observation%nlambda)
! 
! 		print *, xc, fc
! 		
! 	end subroutine funct
! 
! 
! !------------------------------------------------------------
! ! Do initialization
! !------------------------------------------------------------
! 	subroutine nag_init
! 	real(kind=8), allocatable :: x(:), l(:), u(:), w(:), ruser(:)
! 	type(statistics) :: stat
! 	logical :: physical
! 	type(modelo_type), allocatable :: proposed(:)
! 	integer, allocatable :: iw(:), iuser(:)
! 	integer :: ifail, nmax, i, liw, lw
! 	real(kind=8) :: f
! 
! 		nmax = inversion%nparams_invert
! 
! 		stat_nonlinear%nparams = inversion%nparams_invert
! 		allocate(stat_nonlinear%parameters_proposed(stat_nonlinear%nparams))
! 		allocate(proposed_nonlinear(number_of_components))
! 		proposed_nonlinear = model
! 
! 		allocate(x(nmax),l(nmax),u(nmax))
! 
! 		do i = 1, nmax			
! 			l(i) = inversion%range(i,1)
! 			u(i) = inversion%range(i,2)
! 			x(i) = 0.5d0*(l(i) + u(i))
! 		enddo
! 
! 		liw = nmax + 3
! 		lw = max(nmax*(nmax-1)/2+12*nmax,13) + 2		
! 		allocate(iw(lw), w(lw), iuser(1), ruser(1))
! 		ifail = 0
! 
! 		call e04jyf(nmax, 0, funct, l, u, x, f, iw, liw, w, lw, iuser, ruser, ifail)
! 		print *, ifail
! 
! 		deallocate(stat_nonlinear%parameters_proposed)
! 		deallocate(proposed_nonlinear)
! 
! 		stop
! 
! 	end subroutine nag_init


!------------------------------------------------------------
! Carry out the Markov-chain Montecarlo
!------------------------------------------------------------
	subroutine do_mcmc(inv,Observation)
	real(kind=8) :: flast, fproposed, fmax, r, alpha, ran, acceptance_rate1, chi2, chi2_max
	real(kind=8) :: acceptance_rate2, sigma_factor, alpha_theory, chi2last
	real(kind=8) :: foriginal, chi2original
	real(kind=8), allocatable :: parameters_last(:), parameters_proposed(:)
	type(inversion_type) :: inv
	type(stokes_type) :: observation
	type(modelo_type), allocatable :: proposed(:), last(:)
	type(statistics) :: stat
	integer :: i, j, k, accepted, n
	logical :: physical, sigma_working
	character(len=40) :: format_output, nparams_str

		print *, 'Starting...'
		call locate_filling_factors(inv)

! 		call lbfgs_init(inv, Observation)
! 		call nag_init
		
		allocate(proposed(number_of_components))
		allocate(last(number_of_components))
		
		stat%nparams = inv%nparams_invert
		allocate(stat%parameters_proposed(stat%nparams))
		allocate(stat%parameters_last(stat%nparams))
		allocate(stat%parameters_MAP(stat%nparams))
		
		allocate(stat%is_accepted(uniform_proposal_steps))
		stat%is_accepted = 0
		
		do i = 1, stat%nparams
			ran = randomu()
			stat%parameters_last(i) = ran*inv%range_modified(i,1)+(1.d0-ran)*inv%range_modified(i,2)
			stat%parameters_proposed(i) = stat%parameters_last(i)
		enddo
		
		
		allocate(stat%covariance(stat%nparams,stat%nparams))
		allocate(stat%eigenvectors(stat%nparams,stat%nparams))
		allocate(stat%eigenvalues(stat%nparams))
		
		stat%covariance = 0.d0
		
! The first estimation of the covariance matrix is just the diagonal matrix
! with variances equal to a small portion of the interval where the prior is non-zero
		do i = 1, stat%nparams
			stat%covariance(i,i) = (inv%range_modified(i,2)-inv%range_modified(i,1))*0.1d0
		enddo
		
! Theoretical estimation of the alpha constant
		alpha_theory = 2.4d0**2 / stat%nparams
		
		allocate(stat%mean(stat%nparams))
		stat%mean = 0.d0
		stat%alpha = 1.d0
		stat%calculated = .FALSE.
		
		write(nparams_str,FMT='(I2)') stat%nparams
		format_output = '(4X,A7,'//trim(adjustl(nparams_str))//'(F9.4,1X))'
				
		accepted = 0
		n = 0
		proposed = model
		last = proposed
		stat%accepted_models = 0
		
		call locate_filling_factors(inv)
				
		flast = compute_likelihood(model,linea,Observation,inv,Stokes_Syn,physical,chi2,1.d0)

		fmax = flast
		chi2_max = chi2
		chi2last = chi2
		chi2original = chi2
		foriginal = flast

! If starting a new chain
		if (start_or_restart == 0 .or. start_or_restart == 2) then
			print *, 'Starting a new chain...'
			open(unit=12,file=trim(adjustl(fich_chain))//'.header',action='write',status='replace')
			write(12,*) stat%nparams, chain_maxlength-uniform_proposal_steps
			close(12)
			open(unit=12,file=trim(adjustl(fich_chain))//'.data',action='write',status='replace',form='unformatted')			
		else
			print *, 'Appending to an existing chain...'
! When continuing a chain, do not perform any step with a uniform proposal			
			uniform_proposal_steps = 0
			
! Read the previous header to get the number of steps of the chain
			open(unit=12,file=trim(adjustl(fich_chain))//'.header',action='read',status='old')
			read(12,*) i, j			
			close(12)
			
! Test if at least the number of parameters is the same			
			if (stat%nparams /= i) then
				print *, 'The chain has no the same number of parameters as the previous one'
				stop
			endif
			
! Write the new chain header
			open(unit=12,file=trim(adjustl(fich_chain))//'.header',action='write',status='replace')
			write(12,*) stat%nparams, j + chain_maxlength
			close(12)
			
			n = j

! Open the chain data file with the pointer ready to append data
			open(unit=12,file=trim(adjustl(fich_chain))//'.data',action='write',status='old',&
				form='unformatted',position='append')
				
! Read the covariance matrix and the averages from the previous chain
			open(unit=13,file=trim(adjustl(fich_chain))//'.statistics',action='read',status='old')
			read(13,*) i, stat%nparams, accepted
			call lb(13,1)
			read(13,*) stat%mean
			call lb(13,1)
			read(13,*) stat%covariance
			close(13)
		endif
				
		sigma_working = .TRUE.
! MAIN LOOP
		do i = 1, chain_maxlength
		
! Propose and fill the model
			call proposal(inv,stat,i,uniform_proposal_steps)			
			
! Propose with B vector uniform (uniform in B^3 and cos(theta))
			call fill_model_uniformvector_new(inv,stat,proposed)
												
! Calculate the likelihood
! This is a trick that increases the value of sigma in the first steps
! of the chain to allow a better mixing of the chain
			sigma_factor = 1.d0
			if (sigma_working) then
				sigma_factor = (10.d0-1.d0)*exp(-n/(0.03*chain_maxlength))+1.d0
			endif
			
			if (sigma_factor < 1.01d0 .and. sigma_working) then
				sigma_working = .FALSE.
				fmax = foriginal
				chi2_max = chi2original
			endif
				
			fproposed = compute_likelihood(proposed,linea,Observation,inv,&
				Stokes_Syn,physical,chi2,sigma_factor)
									
			stat%is_accepted(mod(n,uniform_proposal_steps)+1) = 0
! If the model has a physical meaning, carry out the MCMC acceptance			
			if (physical) then				
! 				r = fproposed / flast
				r = exp(fproposed - flast)
							
				alpha = min(1.d0,r)
			
				ran = randomu()
				if (ran < alpha) then
					last = proposed
					stat%parameters_last = stat%parameters_proposed
					flast = fproposed
					chi2last = chi2
					if (flast > fmax) then
						fmax = flast
						chi2_max = chi2last
						stat%parameters_MAP = stat%parameters_last
					endif
					accepted = accepted + 1
					stat%accepted_models = stat%accepted_models + 1
					stat%is_accepted(mod(n,uniform_proposal_steps)+1) = 1
				endif
			endif
																		
! Write the value of the parameters and the likelihood
			if (i > uniform_proposal_steps) then
				call write_step_chain(stat,inv,flast,chi2last,acceptance_rate1)
			endif
			

			if (n >= uniform_proposal_steps) then
				call compute_statistics(stat,inv,n,.FALSE.)
			endif
						
! Write some information if proposing with the gaussian proposal
			if (n > uniform_proposal_steps) then				
				
				if (n / 100 == n / 100.d0) then
																		
					acceptance_rate1 = 1.d0*accepted / n * 100.d0
					acceptance_rate2 = &
						sum(stat%is_accepted) / float(uniform_proposal_steps) * 100.d0		
			
					write(*,FMT='(I1,1X,I6,2X,A,F6.2,A,F6.2,A,F9.4,1X,F9.4,A,E11.4,A,E10.4)') &
						sigma_working, n, &
						'Acceptance rate: ', &
						acceptance_rate1, '%/',acceptance_rate2, '% -- (alpha/alpha_th): ', &
						stat%alpha, alpha_theory, ' -- Fmax: ', fmax,&
						' -- chi^2: ', chi2_max / observation%nlambda
					call print_output(stat,inv,format_output)
					
					open(unit=13,file=trim(adjustl(fich_chain))//'.statistics',action='write',&
						status='replace')
					write(13,*) n, stat%nparams, accepted
					write(13,*) 'Averages'
					write(13,*) stat%mean
					write(13,*) 'Covariance matrix'
					write(13,*) stat%covariance
					close(13)
						
! Correct the scaling of the covariance matrix to maintain the acceptance rate
! in the range between 20% and 30%
					if (desired_acceptance > 0.d0) then
						if (acceptance_rate2 > desired_acceptance + 5.d0) then
							stat%alpha = stat%alpha * 1.04d0
						else if (acceptance_rate2 < desired_acceptance - 5.d0) then
							stat%alpha = stat%alpha * 0.96d0
						endif
					endif
					if (stat%alpha < 0.001) then
						stat%alpha = alpha_theory
					endif
					
				endif
			endif
			
			n = n + 1
			
		enddo
		
		print *, 'End.'
		close(12)
						
		do i = 1, 4
			write(*,FMT='(A,I1,A,E12.4)') 'Maxim. diff Stokes ', &
				i,' = ',Observation%sigma_minimum(i)
		enddo
		
! If the global acceptance rate is very small, something has gone bad
! and the analysis of the chains will fail
		open(unit=12,file='finished.state',action='write',status='replace')
		if (acceptance_rate1 < 5.d0) then
			write(12,*) 0
			stop
		else
			write(12,*) 1
		endif
		close(12)

		stat%parameters_proposed = stat%parameters_MAP
! Analyze the resulting chains (1D histograms and confidence intervals)
		call analyze_chains_mcmc(chain_analysis,fich_chain,stat%parameters_proposed)
		
! Carry out a synthesis with the most probable values of the inverted parameters
		stat%parameters_proposed = stat%parameters_MAP
		call fill_model_uniformvector_new(inv,stat,proposed)

		call synthesize(proposed,linea,Observation,Stokes_Syn)
		
		open(unit=12,file=trim(adjustl(fich_chain))//'.stokes_map',&
			action='write',status='replace')
		do i = 1, Stokes_Syn%nlambda
			write(12,FMT='(F10.4,2X,4(E12.4))') Stokes_Syn%lambda(i),&
				Stokes_Syn%stokes(1,i), Stokes_Syn%stokes(2,i), &
				Stokes_Syn%stokes(3,i), Stokes_Syn%stokes(4,i)
		enddo
		close(12)
				
	end subroutine do_mcmc

end module markov_chain_montecarlo

