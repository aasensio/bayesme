module read_config_mod
use vars
use math_functions
implicit none
contains

!--------------------------------------------------
! Read the configuration file
!--------------------------------------------------
	subroutine read_config(modelo,linea,inv_parameter,observed)
	type(inversion_type) :: inv_parameter
	type(stokes_type) :: observed
	type(modelo_type), pointer :: modelo(:)
	type(line_type), pointer :: linea(:)
	character(len=80) :: temp, file, text_proposal
	integer :: i, j, index, k
	integer, allocatable :: indice_linea(:)
	
	real(kind=8) :: wl, elow, mult_low, Jlow, Alow, Blow, eup, mult_up, Jup, Aup, Bup, glow, gup, factor
	integer :: t1, t2, npossible_params
	real(kind=8) :: ff_total, mu, t3, t4, sigma_estimation, maxim, minim
	character(len=1) :: lab_low, Llow, lab_up, Lup
	integer, allocatable :: where_magnetic_field(:)
	integer :: ios, nargs
	integer :: ind_par
	character(len=132) :: config_file, file_stray_light
	logical :: exists

      use_stray_light = 0
	
! Verify if a configuration file is passed as an argument
! If not, use the default 'config' file
		nargs = iargc()
		config_file = 'config'
		if (nargs == 1) then
			call getarg(1,config_file)
		endif		
		
		open(unit=12,file=config_file,status='old',action='read',IOSTAT=ios)
		
! Is the configuration file present? If not, we cannot continue because we don't
! know what to do!!
		if (ios /= 0) then
			write(*,*) "Configuration file not found. Verify that a 'config' file is present."
			stop
		endif
		
		write(*,*) 'Using configuration file : ', trim(adjustl(config_file))

! Read whether to start a new chain or restart a previous chain
! 		call lb(12,2)
! 		read(12,*) start_or_restart

! Read the maximum length of the chain
! 		call lb(12,1)
! 		read(12,*) chain_maxlength
		
! Read the number of steps with uniform proposal
! 		call lb(12,1)		
! 		read(12,*) uniform_proposal_steps
		
! Read the desired acceptance
! 		call lb(12,1)		
! 		read(12,*) desired_acceptance
		
! Read if sampling is uniform in B and theta
		call lb(12,2)		
		read(12,*) uniform_B, uniform_thB

! Read if posteriors for combinations of parameters are desired
		call lb(12,1)
		read(12,*) chain_analysis%n_posterior_combinations

! Read functions to be calculated
		if (chain_analysis%n_posterior_combinations > 0) then
			allocate(chain_analysis%posterior_fun(chain_analysis%n_posterior_combinations))			
			call lb(12,1)
			write(*,*) 'Posteriors for the following functions are also calculated:'
			do i = 1, chain_analysis%n_posterior_combinations				
				read(12,*) chain_analysis%posterior_fun(i)%n_variables

				allocate(chain_analysis%posterior_fun(i)%variables(chain_analysis%posterior_fun(i)%n_variables))
				read(12,*) (chain_analysis%posterior_fun(i)%variables(j),&
					j=1,chain_analysis%posterior_fun(i)%n_variables)	
				read(12,*) chain_analysis%posterior_fun(i)%posterior_function
				write(*,*) i, ' -> ', trim(adjustl(chain_analysis%posterior_fun(i)%posterior_function))
				write(*,*) '    Variables : ', &
					(chain_analysis%posterior_fun(i)%variables(j),j=1,chain_analysis%posterior_fun(i)%n_variables)
			enddo
		endif
		
! Read the file with the observations
		call lb(12,3)
		read(12,*) file
		print *, 'Reading observed profiles...'
		open(unit=13,file=file,action='read',status='old')
		read(13,*)
		read(13,*) observed%nlambda
		print *, 'Number of wavelength points of observed profile : ', observed%nlambda
		
		allocate(observed%lambda(observed%nlambda))
		allocate(observed%stokes(4,observed%nlambda))
		allocate(observed%sigma(4,observed%nlambda))
		allocate(observed%variance_total(4,observed%nlambda))
		allocate(observed%sigma_minimum(4))
		
		observed%sigma_minimum = 1.d10
		
		do i = 1, observed%nlambda
			read(13,*) observed%lambda(i), (observed%stokes(j,i),j=1,4), (observed%sigma(j,i),j=1,4)
		enddo		
		close(13)
		observed%step_lambda = observed%lambda(2)-observed%lambda(1)

! Read the file with the spectral PSF
		using_instrumental_profile = .FALSE.
		inquire (file='instrumental_profile.dat', exist=exists)

! If including PSF
		if (exists) then
			using_instrumental_profile = .TRUE.
			print *, 'Using instrumental profile'

			open(unit=13,file='instrumental_profile.dat',action='read',status='old')
			read(13,*) instrumental_profile%nlambda_original
			allocate(instrumental_profile%lambda_original(instrumental_profile%nlambda_original))
			allocate(instrumental_profile%transmission_original(instrumental_profile%nlambda_original))			
			
			do i = 1, instrumental_profile%nlambda_original
				read(13,*) instrumental_profile%lambda_original(i), instrumental_profile%transmission_original(i)
				instrumental_profile%lambda_original(i) = instrumental_profile%lambda_original(i) * 1.d-3
			enddo

! Calculate lambda step
			instrumental_profile%step_lambda = instrumental_profile%lambda_original(2)-&
				instrumental_profile%lambda_original(1)

! Instrumental profile is sampled finer than the observation
! In this case, 
			if (instrumental_profile%step_lambda < observed%step_lambda) then

				print *, 'Lambda step in observation < lambda step in instrum. profile'

! Calculate new axis for synthesis
				maxim = maxval(observed%lambda)
				minim = minval(observed%lambda)

				instrumental_profile%nlambda_forfilter = (maxim-minim) / instrumental_profile%step_lambda + 1

				allocate(instrumental_profile%lambda_forfilter(instrumental_profile%nlambda_forfilter))

				print *, 'Changing wavelength axis to include instrumental profile : '
				print *, 'Going from ', observed%nlambda, ' to ', instrumental_profile%nlambda_forfilter

				do i = 1, instrumental_profile%nlambda_forfilter
					instrumental_profile%lambda_forfilter(i) = minim + (i-1.d0) * instrumental_profile%step_lambda
				enddo

! Set instrumental profile so that the FFT can be applied later on				
				allocate(instrumental_profile%transmission(instrumental_profile%nlambda_original))


! Find nearest power of 2
				i = 1
      		do while (2.d0**i < instrumental_profile%nlambda_forfilter)
 					i = i + 1
	 			enddo
      		instrumental_profile%n2 = 2.d0**i

      		print *, 'Sizes : ', instrumental_profile%nlambda_original, instrumental_profile%n2,&
      			instrumental_profile%nlambda_forfilter

      		allocate(instrumental_profile%data(instrumental_profile%n2))
      		allocate(instrumental_profile%respns(instrumental_profile%n2))

      		instrumental_profile%respns(1:instrumental_profile%nlambda_original) = &
      			instrumental_profile%transmission_original(1:instrumental_profile%nlambda_original)

      		print *, instrumental_profile%respns(1:instrumental_profile%nlambda_original)

      		instrumental_profile%respns = cshift(instrumental_profile%respns,&
					instrumental_profile%nlambda_original/2)

! Normalize profile
				instrumental_profile%respns = instrumental_profile%respns /&
					sum(instrumental_profile%respns)


			else
				print *, 'Lambda step in observation > lambda step in instrum. profile'
				print *, 'Resampling filter...'
				print *, 'Not implemented yet...'
				stop
			endif
			close(13)
		endif
		
		call lb(12,1)
		read(12,*) file_stray_light
		
		call lb(12,1)
		read(12,*) fich_chain
				
! 		call lb(12,1)
! 		read(12,*) sigma_estimation
! 
! 		if (sigma_estimation /= 0.d0) then
! 		
! 			do j = 1, 4
! ! 			observed%sigma(j,:) = maxval(abs(observed%stokes(j,:)))
! 				observed%sigma(j,:) = sigma_estimation
! 			enddo
! 		endif
		
		call lb(12,1)
		read(12,*) (stokes_weights(j),j=1,4)
! Now read the rest of the config file that defines the number of lines,
! components, which parameters to invert and their range of variation
		call lb(12,1)
		read(12,*) n_lineas
		allocate(indice_linea(n_lineas))
		allocate(linea(n_lineas))
		
		print *, 'Number of lines : ', n_lineas
		
		call lb(12,1)		
		read(12,*) fich_lineas
		
		call lb(12,1)
		do i = 1, n_lineas			
			read(12,*) linea(i)%theory
			read(12,*) indice_linea(i)
		enddo
		
		call lb(12,1)
		read(12,*) mu
		
		call lb(12,1)
		read(12,*) number_of_components
		
		print *, 'Number of components : ', number_of_components
		allocate(modelo(number_of_components))
		
		ff_total = 0.d0
		
		j = 0
		
		npossible_params = number_of_components*(8+n_lineas)
		allocate(inv_parameter%which_to_invert(npossible_params))
		allocate(inv_parameter%which_parameter(npossible_params))
		allocate(inv_parameter%which_component(npossible_params))
		allocate(inv_parameter%which_line(npossible_params))
		
		inv_parameter%which_to_invert = 0
		inv_parameter%which_parameter = 0
		inv_parameter%which_component = 0
		inv_parameter%which_line = 0
		inv_parameter%free_filling_factors = .FALSE.
		allocate(inv_parameter%range(npossible_params,2))
		allocate(inv_parameter%range_modified(npossible_params,2))
		allocate(inv_parameter%filling_factor_position(number_of_components))
		
		open(unit=14,file=trim(adjustl(fich_chain))//'.params_inverted',action='write',status='replace')
				
		allocate(inv_parameter%dimension_space_B(number_of_components))
		allocate(where_magnetic_field(number_of_components))
		
! Dimension of the space where the B field is moving		
		inv_parameter%dimension_space_B = 1
		where_magnetic_field = -99

		ind_par = 0

		inv_parameter%log_prior_volume_inv = 0.d0
		
		do i = 1, number_of_components
		
			!allocate(modelo(i)%kl(n_lineas))
			allocate(modelo(i)%kl_range(2,n_lineas))
			
			print *, '  ---------------------------------'
			print *, '   COMPONENT  ', i
			print *, '  ---------------------------------'
			call lb(12,5)
						
! MAGNETIC FIELD
! 			call lb(12,1)
			read(12,*) t2, modelo(i)%Bfield, t3, t4
			ind_par = ind_par + 1
			write(*,FMT='(A23,F10.5)') 'Magnetic field : ', modelo(i)%Bfield
			modelo(i)%Bfield_range(1) = t3
			modelo(i)%Bfield_range(2) = t4
			if (t2 == 1) then
				j = j + 1
				inv_parameter%which_to_invert(j) = ind_par
				inv_parameter%which_parameter(j) = 1
				inv_parameter%which_component(j) = i
				inv_parameter%range(j,1) = t3
				inv_parameter%range(j,2) = t4
				inv_parameter%range_modified(j,1) = t3
				inv_parameter%range_modified(j,2) = t4				
				where_magnetic_field(i) = j
				inv_parameter%log_prior_volume_inv = inv_parameter%log_prior_volume_inv - log((t4-t3))
				write(*,FMT='(20X,A20)') '--- will be inverted'
				write(14,FMT='(A1,I1)') 'B', i
			endif

			
! INCLINATION
			call lb(12,1)		
			read(12,*) t2, modelo(i)%theta, t3, t4
			ind_par = ind_par + 1
			write(*,FMT='(A23,F10.5)') 'Inclination : ', modelo(i)%theta
			modelo(i)%theta_range(1) = t3
			modelo(i)%theta_range(2) = t4
			if (t2 == 1) then
				j = j + 1
				inv_parameter%which_to_invert(j) = ind_par
				inv_parameter%which_parameter(j) = 2
				inv_parameter%which_component(j) = i
				inv_parameter%range(j,1) = t3
				inv_parameter%range(j,2) = t4
				
! Uniform sampling in cos(thetaB) or in thetaB
				if (uniform_thB == 1) then
					inv_parameter%range_modified(j,1) = cos(t3*PI/180.d0)
					inv_parameter%range_modified(j,2) = cos(t4*PI/180.d0)
				else
					inv_parameter%range_modified(j,1) = t3
					inv_parameter%range_modified(j,2) = t4
				endif
				
! Uniform in the modulues or in the square/cube of the modulus
				if (uniform_B == 1) then
					inv_parameter%dimension_space_B(i) = inv_parameter%dimension_space_B(i) + 1
				endif
				inv_parameter%log_prior_volume_inv = inv_parameter%log_prior_volume_inv - log((t4-t3))
				write(*,FMT='(20X,A20)') '--- will be inverted'
				write(14,FMT='(A5,I1)') 'theta', i
			endif
			
! AZIMUTH
			call lb(12,1)
			read(12,*) t2, modelo(i)%chi, t3, t4
			ind_par = ind_par + 1
			write(*,FMT='(A23,F10.5)') 'Azimuth : ', modelo(i)%chi
			modelo(i)%chi_range(1) = t3
			modelo(i)%chi_range(2) = t4
			if (t2 == 1) then
				j = j + 1
				inv_parameter%which_to_invert(j) = ind_par
				inv_parameter%which_parameter(j) = 3
				inv_parameter%which_component(j) = i
				inv_parameter%range(j,1) = t3
				inv_parameter%range(j,2) = t4
				inv_parameter%range_modified(j,1) = t3
				inv_parameter%range_modified(j,2) = t4
				if (uniform_B == 1) then
					inv_parameter%dimension_space_B(i) = inv_parameter%dimension_space_B(i) + 1
				endif
				inv_parameter%log_prior_volume_inv = inv_parameter%log_prior_volume_inv - log((t4-t3))
				write(*,FMT='(20X,A20)') '--- will be inverted'
				write(14,FMT='(A7,I1)') 'azimuth', i
			endif						

! DOPPLER WIDTH
			call lb(12,1)		
			read(12,*) t2, modelo(i)%doppler, t3, t4
			ind_par = ind_par + 1
			write(*,FMT='(A23,F10.5)') 'Doppler width : ', modelo(i)%doppler
			modelo(i)%doppler_range(1) = t3
			modelo(i)%doppler_range(2) = t4
			if (t2 == 1) then
				j = j + 1
				inv_parameter%which_to_invert(j) = ind_par
				inv_parameter%which_parameter(j) = 4
				inv_parameter%which_component(j) = i
				inv_parameter%range(j,1) = t3
				inv_parameter%range(j,2) = t4
				inv_parameter%range_modified(j,1) = t3
				inv_parameter%range_modified(j,2) = t4
				inv_parameter%log_prior_volume_inv = inv_parameter%log_prior_volume_inv - log((t4-t3))
				write(*,FMT='(20X,A20)') '--- will be inverted'
				write(14,FMT='(A3,I1)') 'vth', i
			endif
			

! MACROSCOPIC VELOCITY
			call lb(12,1)		
			read(12,*) t2, modelo(i)%vmac, t3, t4
			ind_par = ind_par + 1
			write(*,FMT='(A23,F10.5)') 'Macroscopic velocity : ', modelo(i)%vmac
			modelo(i)%vmac_range(1) = t3 !* 1.d5
			modelo(i)%vmac_range(2) = t4 !* 1.d5
			if (t2 == 1) then
				j = j + 1
				inv_parameter%which_to_invert(j) = ind_par
				inv_parameter%which_parameter(j) = 5
				inv_parameter%which_component(j) = i
				inv_parameter%range(j,1) = t3 !* 1.d5
				inv_parameter%range(j,2) = t4 !* 1.d5
				inv_parameter%range_modified(j,1) = t3
				inv_parameter%range_modified(j,2) = t4
				inv_parameter%log_prior_volume_inv = inv_parameter%log_prior_volume_inv - log((t4-t3))
				write(*,FMT='(20X,A20)') '--- will be inverted'
				write(14,FMT='(A4,I1)') 'vmac', i
			endif
			
			modelo(i)%vmac = modelo(i)%vmac ! * 1.d5

! DAMPING
			call lb(12,1)		
			read(12,*) t2, modelo(i)%damping, t3, t4
			ind_par = ind_par + 1
			write(*,FMT='(A23,F10.5)') 'Damping : ', modelo(i)%damping
			modelo(i)%damping_range(1) = t3
			modelo(i)%damping_range(2) = t4
			if (t2 == 1) then
				j = j + 1
				inv_parameter%which_to_invert(j) = ind_par
				inv_parameter%which_parameter(j) = 6
				inv_parameter%which_component(j) = i
				inv_parameter%range(j,1) = t3
				inv_parameter%range(j,2) = t4
				inv_parameter%range_modified(j,1) = t3
				inv_parameter%range_modified(j,2) = t4
				inv_parameter%log_prior_volume_inv = inv_parameter%log_prior_volume_inv - log((t4-t3))
				write(*,FMT='(20X,A20)') '--- will be inverted'
				write(14,FMT='(A4,I1)') 'damp', i
			endif
			

! BETA
			call lb(12,1)		
			read(12,*) t2, modelo(i)%beta, t3, t4
			ind_par = ind_par + 1
			write(*,FMT='(A23,F10.5)') 'B1/B0 : ', modelo(i)%beta
			modelo(i)%beta_range(1) = t3
			modelo(i)%beta_range(2) = t4
			if (t2 == 1) then
				j = j + 1
				inv_parameter%which_to_invert(j) = ind_par
				inv_parameter%which_parameter(j) = 7
				inv_parameter%which_component(j) = i
				inv_parameter%range(j,1) = t3
				inv_parameter%range(j,2) = t4
				inv_parameter%range_modified(j,1) = t3
				inv_parameter%range_modified(j,2) = t4
				inv_parameter%log_prior_volume_inv = inv_parameter%log_prior_volume_inv - log((t4-t3))
				write(*,FMT='(20X,A20)') '--- will be inverted'
				write(14,FMT='(A4,I1)') 'beta', i
			endif			
			
! LINE STRENGTH			
			call lb(12,1)		
			do k = 1, n_lineas
				read(12,*) t2, modelo(i)%kl(k), t3, t4
				ind_par = ind_par + 1
				write(*,FMT='(A18,I2,A3,F10.5)') 'kl(',k,'): ', modelo(i)%kl(k)
				modelo(i)%kl_range(1,k) = t3
				modelo(i)%kl_range(2,k) = t4
				if (t2 == 1) then
					j = j + 1
					inv_parameter%which_to_invert(j) = ind_par
					inv_parameter%which_parameter(j) = 8
					inv_parameter%which_component(j) = i
					inv_parameter%which_line(j) = k
					inv_parameter%range(j,1) = t3
					inv_parameter%range(j,2) = t4
					inv_parameter%range_modified(j,1) = t3
					inv_parameter%range_modified(j,2) = t4
					inv_parameter%log_prior_volume_inv = inv_parameter%log_prior_volume_inv - log((t4-t3))
					write(*,FMT='(20X,A20)') '--- will be inverted'
					write(14,FMT='(A3,I1,A1,I1)') 'eta', i, '_', k
				endif
			enddo

! FILLING FACTOR
			call lb(12,1)
			read(12,*) t2, modelo(i)%filling_factor, t3, t4
			ind_par = ind_par + 1
			
! If negative, then use a stray light component
			modelo%stray_light_component = .FALSE.
			if (modelo(i)%filling_factor < 0) then
				write(*,*) 'Using stray light component'
				modelo(i)%stray_light_component = .TRUE.
				modelo(i)%filling_factor = abs(modelo(i)%filling_factor)
				use_stray_light = i
				if (number_of_components < 2) then
					print *, 'Use of stray light needs at least two components'
					stop
				endif
			endif
			
			write(*,FMT='(A23,F10.5)') 'Filling factor : ', modelo(i)%filling_factor
			modelo(i)%filling_factor_range(1) = t3
			modelo(i)%filling_factor_range(2) = t4			
			if (t2 == 1) then
				inv_parameter%free_filling_factors = .TRUE. 
				j = j + 1
				inv_parameter%which_to_invert(j) = ind_par
				inv_parameter%which_parameter(j) = 9
				inv_parameter%which_component(j) = i
				inv_parameter%range(j,1) = t3
				inv_parameter%range(j,2) = t4
				inv_parameter%range_modified(j,1) = t3
				inv_parameter%range_modified(j,2) = t4
				inv_parameter%log_prior_volume_inv = inv_parameter%log_prior_volume_inv - log((t4-t3))
				write(*,FMT='(20X,A20)') '--- will be inverted'
				write(14,FMT='(A5,I1)') 'alpha', i
			endif
			inv_parameter%filling_factor_position(i) = j
										
			ff_total = ff_total + modelo(i)%filling_factor
			
		enddo
				
		close(14)

! Number of parameters to invert
		inv_parameter%nparams_invert = j
		write(*,FMT='(A,I2,A)') 'Inverting ', inv_parameter%nparams_invert, ' parameters'
		if (uniform_thB == 1) then
			write(*,*) 'Proposing uniformly in cos(thetaB)'
		else
			write(*,*) 'Proposing uniformly in thetaB'
		endif
		write(*,*) 'Dimensions of B : '
		
! Write how we propose values for the field
		do i = 1, number_of_components
			
! Modify the ranges to take into account the correct proposal
			text_proposal = ' - Not inverted'
			if (where_magnetic_field(i) /= -99) then
				inv_parameter%range_modified(where_magnetic_field(i),1) = &
					inv_parameter%range_modified(where_magnetic_field(i),1)**inv_parameter%dimension_space_B(i)
				inv_parameter%range_modified(where_magnetic_field(i),2) = &
					inv_parameter%range_modified(where_magnetic_field(i),2)**inv_parameter%dimension_space_B(i)			
				
				select case(inv_parameter%dimension_space_B(i))
					case(1) 
						text_proposal = ' - Proposing as d(r)'
					case(2) 
						text_proposal = ' - Proposing as d(r^2)'
					case(3) 
						text_proposal = ' - Proposing as d(r^3)'
				end select
			endif
			write(*,FMT='(A,I1,A,I2,A)') '   - Component ', i, ' : ',&
				inv_parameter%dimension_space_B(i),	text_proposal			
		enddo
		
		if (ff_total /= 1.d0) then
			print *, 'Wrong filling factors. They do not add up to 1'
			stop
		endif
								
		
		do i = 1, number_of_components
			modelo(i)%mu = mu
		enddo
		
		print *, 'mu : ', modelo(1)%mu
				
		close(12)
		
		deallocate(where_magnetic_field)
		
! Read stray light profile
		if (use_stray_light > 0) then
			open(unit=13,file=file_stray_light,action='read',status='old')
			read(13,*) modelo(use_stray_light)%stray_light_nlambda, modelo(use_stray_light)%npix_stray
			write(*,*) 'Reading stray-light profile'
			
			allocate(modelo(use_stray_light)%straylight(4,modelo(use_stray_light)%stray_light_nlambda))
			allocate(modelo(use_stray_light)%straylight_cov(4,modelo(use_stray_light)%stray_light_nlambda))
			
			do i = 1, modelo(use_stray_light)%stray_light_nlambda
				read(13,*) t1, (modelo(use_stray_light)%straylight(j,i),j=1,4), (modelo(use_stray_light)%straylight_cov(j,i),j=1,4)
			enddo		
			close(13)
		endif
						
! Read the file with the line data						
		do j = 1, n_lineas
			open(unit=12,file=fich_lineas,status='old',action='read')
			if (linea(j)%theory == 'ZE') then
				read(12,*)
				do i = 1, indice_linea(j)
					call lb(12,1)
				enddo
				read(12,*) index, linea(j)%wave0, linea(j)%glow, linea(j)%Jlow, linea(j)%gup, linea(j)%Jup,&
					linea(j)%lambda_init, linea(j)%lambda_end, linea(j)%lambda_step
			endif
		
		
			if (linea(j)%theory == 'HF') then
				factor = 0.001d0  ! Although the data of Pickering (1996, ApJS 107, 811) is given in mK, she
									! says that 1mK = 0.001 cm^-1, and not PK / (PH*PC)*1.d-3
				read(12,*)
				read(12,*)
				do i = 0, indice_linea(j)-1
					read(12,*)
				enddo
				read(12,fmt='(i2,1x,f9.3,3x,f9.3,3x,a1,1x,f3.1,1x,a1,1x,f3.1,2x,f6.2,2x,f4.1,3x,f9.3,3x,&
					&a1,1x,f3.1,1x,a1,1x,f3.1,2x,f6.2,2x,f4.1,3x,f5.3,3x,f5.3,3x,f3.1,2x,f9.3,2x,f9.3,2x,f5.3)') &
					index, wl, elow, lab_low, mult_low, llow, jlow, alow, blow, eup, lab_up, mult_up, lup, jup, &
					aup, bup, glow, gup, linea(j)%i, linea(j)%lambda_init, linea(j)%lambda_end, linea(j)%lambda_step
					
							
				select case(Llow)
					case('S') 
						linea(j)%Ll = 0
					case('P') 
						linea(j)%Ll = 1
					case('D') 
						linea(j)%Ll = 2
					case('F') 
						linea(j)%Ll = 3
					case('G') 
						linea(j)%Ll = 4
					case('H') 
						linea(j)%Ll = 5
					case('I') 
						linea(j)%Ll = 6
					case('J') 
						linea(j)%Ll = 7
				end select
				linea(j)%Su = (mult_up-1.d0) / 2.d0
				linea(j)%Jup = Jup
				linea(j)%Au = Aup * factor
				linea(j)%Bu = Bup * factor

				select case(Lup)
					case('S') 
						linea(j)%Lu = 0
					case('P') 
						linea(j)%Lu = 1
					case('D') 
						linea(j)%Lu = 2
					case('F') 
						linea(j)%Lu = 3
					case('G') 
						linea(j)%Lu = 4
					case('H') 
						linea(j)%Lu = 5
					case('I') 
						linea(j)%Lu = 6
					case('J') 
						linea(j)%Lu = 7
				end select
				linea(j)%Sl = (mult_low-1.d0) / 2.d0
				linea(j)%Jlow = Jlow
				linea(j)%Al = Alow * factor
				linea(j)%Bl = Blow * factor
				
				linea(j)%wave0 = wl
			
				print *, 'Line with hyperfine structure'
				write(*,FMT='(A,F3.1)') 'I : ', linea(j)%I
				write(*,FMT='(6(A,F3.1,2X))') 'Su : ', linea(j)%Su, 'Sl : ',linea(j)%Sl, 'Lu : ',linea(j)%Lu, &
					'Ll : ',linea(j)%Ll, 'Ju : ',linea(j)%Jup, 'Jl : ',linea(j)%Jlow
				write(*,FMT='(4(A,F6.2,2X))') 'Au : ',linea(j)%Au / factor, 'Bu : ',linea(j)%Bu / factor, &
					'Al : ',linea(j)%Al / factor, 'Bl : ', linea(j)%Bl / factor
		
			endif
			close(12)
		enddo
		
		deallocate(indice_linea)
		
	end subroutine read_config
		
end module read_config_mod
