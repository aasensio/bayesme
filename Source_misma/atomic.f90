module atomic_functions
use math_functions
use vars
use hyperfine
use convol
implicit none
contains

	
!-----------------------------------------------------------------
! Return the profiles weighted by the strength of the components for a given frequency
! It returns zeeman_profile(q,n_depths) with
!  q=1  Mlow=Mup-1  (sigma blue)
!  q=2  Mlow=Mup    (sigma pi)
!  q=3  Mlow=Mup+1  (sigma red)
!-----------------------------------------------------------------	
	subroutine zeeman_profile(Stokes_Syn,model,linea,zeeman_voigt,zeeman_faraday)
	type(stokes_type) :: Stokes_Syn
	type(modelo_type) :: model
	type(line_type) :: linea
	real(kind=8) :: zeeman_voigt(:,:), zeeman_faraday(:,:)
	integer :: n, nlow, nup, iup, ilow, i_pi, i_blue, i_red, cual
	real(kind=8) :: Mup, Mlow, strength, va, vb, splitting
	real(kind=8), allocatable :: profile(:,:), v(:)
	
		n = size(zeeman_voigt(1,:))		
		
		allocate(profile(2,n))
		allocate(v(n))
		
		nup = 2*linea%Jup+1
		nlow = 2*linea%Jlow+1
		
		zeeman_voigt = 0.d0
		zeeman_faraday = 0.d0
		
		i_red = 0
		i_pi = 0
		i_blue = 0
		
! Macroscopic velocity is in km/s		
		v = (Stokes_Syn%lambda-linea%wave0) / (model%doppler)
		va = linea%wave0 * 1.d5 * model%vmac / (PC*model%doppler)
		vb = linea%wave0**2 * model%Bfield * 4.6686d-13 / model%doppler
				
		do iup = 1, nup
			Mup = linea%Jup + 1 - iup  ! Mup=J...-J
			do ilow = 1, 3
				Mlow = Mup-2+ilow			! Mlow=Mup-1,Mup,Mup+1 (the allowed transitions)
				if (abs(Mlow) <= linea%Jlow) then
					
					if (ilow == 1) then
						i_blue = i_blue + 1
						cual = i_blue
					endif
					if (ilow == 2) then
						i_pi = i_pi + 1
						cual = i_pi
					endif
					if (ilow == 3) then
						i_red = i_red + 1
						cual = i_red
					endif
					
					strength = strength_zeeman(linea%Jup,linea%Jlow,Mup,Mlow)
					splitting = linea%gup*Mup - linea%glow*Mlow

					profile = fvoigt_zeeman(model%damping,v-va+vb*splitting)
					zeeman_voigt(ilow,:) = zeeman_voigt(ilow,:) + strength * profile(1,:) / sqrt(PI)
					zeeman_faraday(ilow,:) = zeeman_faraday(ilow,:) + strength * profile(2,:) / sqrt(PI)		
				endif
			enddo
		enddo	
		
		deallocate(profile)
		deallocate(v)
		
	end subroutine zeeman_profile
	
!-----------------------------------------------------------------
! Return the profiles weighted by the strength of the components for a given frequency for the hyperfine case
! It returns zeeman_profile(q,n_depths) with
!  q=1  Mlow=Mup-1  (sigma blue)
!  q=2  Mlow=Mup    (sigma pi)
!  q=3  Mlow=Mup+1  (sigma red)
!-----------------------------------------------------------------	
	subroutine zeeman_profile_hyperfine(Stokes_Syn,model,linea,zeeman_voigt,zeeman_faraday)
	type(stokes_type) :: Stokes_Syn
	type(modelo_type) :: model
	type(line_type) :: linea
	real(kind=8) :: zeeman_voigt(:,:), zeeman_faraday(:,:)
	integer :: n, nlow, nup, iup, ilow, i_pi, i_blue, i_red, cual
	real(kind=8) :: Mup, Mlow, strength, va, vb, splitting
	real(kind=8), allocatable :: profile(:,:), v(:)
	
	integer :: nfl_up, f2m_up, nfl_low, f2m_low, nPi, nRight, nLeft, loop
	integer, allocatable :: nflev_up(:), nflev_low(:)
	real(kind=8), allocatable :: energy_up(:,:), c_up(:,:,:), energy_low(:,:), c_low(:,:,:)
	real(kind=8), allocatable :: SplitPi(:), SplitRight(:), SplitLeft(:), StrPi(:), StrRight(:), StrLeft(:)
	
		n = size(zeeman_voigt(1,:))
		
		allocate(profile(2,n))
		allocate(v(n))
				
! Calculate the size of the energy and eigenvectors arrays
		call hyper_size(linea%Jup, linea%I, nfl_up, f2m_up)
		call hyper_size(linea%Jlow, linea%I, nfl_low, f2m_low)
						
		allocate(nflev_up(0:2*f2m_up))	
		allocate(energy_up(0:2*f2m_up,0:nfl_up-1))
		allocate(c_up(0:2*f2m_up,0:nfl_up-1,0:f2m_up))
		nflev_up = 0
		energy_up = 0.d0
		c_up = 0.d0
		
		allocate(nflev_low(0:2*f2m_low))
		allocate(energy_low(0:2*f2m_low,0:nfl_low-1))
		allocate(c_low(0:2*f2m_low,0:nfl_low-1,0:f2m_low))
		nflev_low = 0
		energy_low = 0.d0
		c_low = 0.d0
		
! Diagonalize the hyperfine+magnetic Hamiltonian of the upper and lower levels
		call hyper_diagonalize(linea%Lu, linea%Su, linea%Jup, linea%I, linea%Au, linea%Bu, model%Bfield, &
			nflev_up, energy_up, c_up, f2m_up, nfl_up)
		call hyper_diagonalize(linea%Ll, linea%Sl, linea%Jlow, linea%I, linea%Al, linea%Bl, model%Bfield, &
			nflev_low, energy_low, c_low, f2m_low, nfl_low)		
		
! Calculate the number of Zeeman components for the transition
		call hyper_components(linea%Jup, linea%Jlow, linea%I, f2m_up, f2m_low, nflev_up, nflev_low, nPi, nLeft, nRight)
				
		allocate(SplitPi(nPi))
		allocate(SplitRight(nRight))
		allocate(SplitLeft(nLeft))
		
		allocate(StrPi(nPi))
		allocate(StrRight(nRight))
		allocate(StrLeft(nLeft))

! Calculate the Zeeman pattern, returning the splitting and strength of each pi and sigma component		
		call hyper_pattern(linea%I, linea%Lu, linea%Su, linea%Jup, linea%Au, linea%Bu, linea%Ll, linea%Sl, linea%Jlow, &
			linea%Al, linea%Bl, model%Bfield, &
			nflev_up, energy_up, c_up, f2m_up, nfl_up, nflev_low, energy_low, c_low, f2m_low, nfl_low, &
			SplitPi, SplitRight, SplitLeft, StrPi, StrRight, StrLeft)
		
		
		
		zeeman_voigt = 0.d0
		zeeman_faraday = 0.d0
		
		
		v = (Stokes_Syn%lambda-linea%wave0) / model%doppler
		va = linea%wave0 * model%vmac * 1.d5 / (PC*model%doppler)
		vb = 1.d-8 * linea%wave0**2 / (model%doppler)
		
		do loop = 1, nLeft
			profile = fvoigt_zeeman(model%damping,v-va+vb*SplitLeft(loop))			
			strength = StrLeft(loop)
			zeeman_voigt(1,:) = zeeman_voigt(1,:) + strength * profile(1,:) / sqrt(PI)
			zeeman_faraday(1,:) = zeeman_faraday(1,:) + strength * profile(2,:) / sqrt(PI)
		enddo
		
		do loop = 1, nPi
			profile = fvoigt_zeeman(model%damping,v-va+vb*SplitPi(loop))
			strength = StrPi(loop)
			zeeman_voigt(2,:) = zeeman_voigt(2,:) + strength * profile(1,:) / sqrt(PI)
			zeeman_faraday(2,:) = zeeman_faraday(2,:) + strength * profile(2,:) / sqrt(PI)
		enddo
		
		do loop = 1, nRight
			profile = fvoigt_zeeman(model%damping,v-va+vb*SplitRight(loop))
			strength = StrRight(loop)
			zeeman_voigt(3,:) = zeeman_voigt(3,:) + strength * profile(1,:) / sqrt(PI)
			zeeman_faraday(3,:) = zeeman_faraday(3,:) + strength * profile(2,:) / sqrt(PI)
		enddo
		
		
						
		deallocate(profile)
		deallocate(v)
		
	end subroutine zeeman_profile_hyperfine	
	
!-----------------------------------------------------------------
! Return the seven independent elements of the absorption matrix
! Remember that, zeeman_voigt(q,:) and zeeman_faraday(q,:) have
!  q=1  Mlow=Mup-1  (sigma blue)
!  q=2  Mlow=Mup    (sigma pi)
!  q=3  Mlow=Mup+1  (sigma red)
!-----------------------------------------------------------------	
	subroutine zeeman_opacity(model,line,zeeman_voigt,zeeman_faraday,ki,kq,ku,kv,fq,fu,fv)
	type(modelo_type) :: model
	integer :: line
	real(kind=8) :: zeeman_voigt(:,:), zeeman_faraday(:,:)
	real(kind=8) :: ki(:), kq(:), ku(:), kv(:), fq(:), fu(:), fv(:)
	real(kind=8) :: sin_theta, cos_theta, sin_2chi, cos_2chi

		sin_theta = sin(model%theta * PI / 180.d0)
		cos_theta = cos(model%theta * PI / 180.d0)
		
		sin_2chi = sin(2.d0 * model%chi * PI / 180.d0)
		cos_2chi = cos(2.d0 * model%chi * PI / 180.d0)

! Classical absorption coefficients
		ki = 0.5d0 * (zeeman_voigt(2,:)*sin_theta**2 + 0.5d0*(zeeman_voigt(1,:)+zeeman_voigt(3,:))*(1.d0+cos_theta**2))  ! eta_I
		kq = 0.5d0 * (zeeman_voigt(2,:) - 0.5d0*(zeeman_voigt(1,:)+zeeman_voigt(3,:))) * sin_theta**2*cos_2chi  ! eta_Q
		ku = 0.5d0 * (zeeman_voigt(2,:) - 0.5d0*(zeeman_voigt(1,:)+zeeman_voigt(3,:))) * sin_theta**2*sin_2chi  ! eta_U
		kv = 0.5d0 * (zeeman_voigt(3,:)-zeeman_voigt(1,:)) * cos_theta  ! eta_V

! Magneto-optical coefficients		
		fq = 0.5d0 * (zeeman_faraday(2,:) - 0.5d0*(zeeman_faraday(1,:)+zeeman_faraday(3,:))) * sin_theta**2*cos_2chi  ! rho_Q
		fu = 0.5d0 * (zeeman_faraday(2,:) - 0.5d0*(zeeman_faraday(1,:)+zeeman_faraday(3,:))) * sin_theta**2*sin_2chi  ! rho_U
		fv = 0.5d0 * (zeeman_faraday(3,:)-zeeman_faraday(1,:)) * cos_theta  ! rho_V
		
! 		ki = model%kl(line) * ki
! 		kq = model%kl(line) * kq
! 		ku = model%kl(line) * ku
! 		kv = model%kl(line) * kv
! 		fq = model%kl(line) * fq
! 		fu = model%kl(line) * fu
! 		fv = model%kl(line) * fv
				
	end subroutine zeeman_opacity	
	
!-----------------------------------------------------------------
! Add the atomic opacity to the opacity including the effect of a magnetic field
!-----------------------------------------------------------------	
	subroutine ApplyFilter(Stokes_Syn, Filter)	
	type(stokes_type) :: Stokes_Syn
	type(instrumental_profile_type) :: Filter
	integer :: i, j

		do i = 1, 4

			call lin_interpol(Stokes_Syn%lambda(1:Stokes_Syn%nlambda), &
				Stokes_Syn%stokes(i,1:Stokes_Syn%nlambda),&
				Filter%lambda_forfilter(1:Filter%nlambda_forfilter),&
				Filter%data(1:Filter%nlambda_forfilter))

			Filter%data(Filter%nlambda_forfilter+1:Filter%n2) = Filter%data(Filter%nlambda_forfilter)

! 			do j = 1, Filter%n2
! 				write(25,*) Filter%lambda_forfilter(j), Filter%data(j), Filter%respns(j)
! 			enddo

			call convolve(Filter%data, Filter%respns)

! 			do j = 1, Stokes_Syn%nlambda
! 				write(26,*) Stokes_Syn%lambda(j), Stokes_Syn%stokes(i,j)
! 			enddo
			
! 			do j = 1, Filter%nlambda_forfilter
! 				write(27,*) Filter%lambda_forfilter(j), Filter%data(j), Filter%respns(j)
! 			enddo


			call lin_interpol(Filter%lambda_forfilter(1:Filter%nlambda_forfilter),&
				Filter%data(1:Filter%nlambda_forfilter),&
				Stokes_Syn%lambda(1:Stokes_Syn%nlambda), &
				Stokes_Syn%stokes(i,1:Stokes_Syn%nlambda))

! 			do j = 1, Stokes_Syn%nlambda
! 				write(27,*) Stokes_Syn%lambda(j), Stokes_Syn%stokes(i,j)
! 			enddo

! 			stop
		enddo
			
	end subroutine ApplyFilter

!-----------------------------------------------------------------
! Add the atomic opacity to the opacity including the effect of a magnetic field
!-----------------------------------------------------------------	
	subroutine synthesize(model,linea,Observation,Stokes_Syn)
	type(modelo_type) :: model(:)
	type(stokes_type) :: Stokes_Syn, Observation
	type(line_type) :: linea(:)
	integer :: i, j, k, n, n_lineas, nHR
	real(kind=8), allocatable :: ki(:), kq(:), ku(:), kv(:), fq(:), fu(:), fv(:), stokes(:,:), delta(:)
	real(kind=8), allocatable :: ki_partial(:), kq_partial(:), ku_partial(:), kv_partial(:)
	real(kind=8), allocatable :: fq_partial(:), fu_partial(:), fv_partial(:)
	real(kind=8) :: factor1 ,factor2, lmax, lmin, lstep, lshift
	character(len=1) :: str
		
		
		n_lineas = size(linea)
					
		n = Observation%nlambda
		nHR = filter_profiles%nlambda
		Stokes_Syn%nlambda = n
		Stokes_Syn%lambda = Observation%lambda
			
		if (.not.associated(Stokes_Syn%lambda)) allocate(Stokes_Syn%lambda(n))
		if (.not.associated(Stokes_Syn%stokes)) allocate(Stokes_Syn%stokes(4,n))
		Stokes_Syn%stokes = 0.d0

! If we want to do the synthesis with filters
		if (spectro_or_filter == 1) then
			if (.not.associated(Stokes_Syn_HR(1)%lambda)) allocate(Stokes_Syn_HR(1)%lambda(nHR))
			if (.not.associated(Stokes_Syn_HR(1)%stokes)) allocate(Stokes_Syn_HR(1)%stokes(4,nHR))
			Stokes_Syn_HR(1)%stokes = 0.d0
			Stokes_Syn_HR(1)%lambda = filter_profiles%lambda
			Stokes_Syn_HR%nlambda = nHR

			n = nHR
		endif

		allocate(zeeman_voigt(3,n))
		allocate(zeeman_faraday(3,n))
	
		allocate(ki_partial(n))
		allocate(kq_partial(n))
		allocate(ku_partial(n))
		allocate(kv_partial(n))
		allocate(fq_partial(n))
		allocate(fu_partial(n))
		allocate(fv_partial(n))
			
		allocate(ki(n))
		allocate(kq(n))
		allocate(ku(n))
		allocate(kv(n))
		allocate(fq(n))
		allocate(fu(n))
		allocate(fv(n))
						
		allocate(delta(n))
		allocate(stokes(4,n))
					
		do j = 1, number_of_components
		
! If it is a stray-light contamination
			if (model(j)%stray_light_component) then				
				if (model(j)%stray_light_nlambda /= Stokes_Syn%nlambda) then
					print *, 'Not consistent number of wavelengths in stray light profile...'
					stop
				endif

! Include possibility of a velocity shift
				stokes = model(j)%filling_factor * model(j)%straylight

! Wavelength shift
				lshift = linea(1)%wave0 * model(j)%vmac * 1.d5 / PC

! Transform the shift to pixels (taking into account the step in wavelength)
				lshift = lshift / Observation%step_lambda
				
				do i = 1, 4
					stokes(i,:) = fft_shift(stokes(i,:), lshift)
				enddo
				Stokes_Syn%stokes = Stokes_Syn%stokes + stokes
			else
! Or a standard model

				factor1 = 1.d0 / (1.d0 + model(j)%beta*model(j)%mu)
				factor2 = -model(j)%beta*model(j)%mu * factor1
				
				ki = 0.d0
				kq = 0.d0
				ku = 0.d0
				kv = 0.d0
				fq = 0.d0
				fu = 0.d0
				fv = 0.d0
			
				do i = 1, n_lineas
					if (linea(i)%theory == 'ZE') then
						if (spectro_or_filter == 1) then
							call zeeman_profile(Stokes_Syn_HR(1),model(j),linea(i),zeeman_voigt,zeeman_faraday)
						else
							call zeeman_profile(Stokes_Syn,model(j),linea(i),zeeman_voigt,zeeman_faraday)
						endif
					endif
					if (linea(i)%theory == 'HF') then
						if (spectro_or_filter == 1) then
							call zeeman_profile_hyperfine(Stokes_Syn_HR(1),model(j),linea(i),zeeman_voigt,zeeman_faraday)
						else
							call zeeman_profile_hyperfine(Stokes_Syn,model(j),linea(i),zeeman_voigt,zeeman_faraday)
						endif
					endif
				
					call zeeman_opacity(model(j),i,zeeman_voigt,zeeman_faraday,ki_partial,kq_partial,&
						ku_partial,kv_partial,fq_partial,fu_partial,fv_partial)
						
					ki = ki + ki_partial * model(j)%kl(i)
					kq = kq + kq_partial * model(j)%kl(i)
					ku = ku + ku_partial * model(j)%kl(i)
					kv = kv + kv_partial * model(j)%kl(i)
					fq = fq + fq_partial * model(j)%kl(i)
					fu = fu + fu_partial * model(j)%kl(i)
					fv = fv + fv_partial * model(j)%kl(i)
				enddo

! Do the optically thin addition (only valid for two components)
				do k = 1, number_of_components

					if (k /= j .and. model(k)%optically_thin == j) then

! Introduce weight in the component for optically thin adding a second component
! Add it to the emission coefficients
						ki = ki * (1.d0 - model(k)%filling_factor)
						kq = kq * (1.d0 - model(k)%filling_factor)
						ku = ku * (1.d0 - model(k)%filling_factor)
						kv = kv * (1.d0 - model(k)%filling_factor)
						fq = fq * (1.d0 - model(k)%filling_factor)
						fu = fu * (1.d0 - model(k)%filling_factor)
						fv = fv * (1.d0 - model(k)%filling_factor)

						do i = 1, n_lineas
							if (linea(i)%theory == 'ZE') then
								if (spectro_or_filter == 1) then
									call zeeman_profile(Stokes_Syn_HR(1),model(k),linea(i),zeeman_voigt,zeeman_faraday)
								else
									call zeeman_profile(Stokes_Syn,model(k),linea(i),zeeman_voigt,zeeman_faraday)
								endif
							endif
							if (linea(i)%theory == 'HF') then
								if (spectro_or_filter == 1) then
									call zeeman_profile_hyperfine(Stokes_Syn_HR(1),model(k),linea(i),zeeman_voigt,zeeman_faraday)
								else
									call zeeman_profile_hyperfine(Stokes_Syn,model(k),linea(i),zeeman_voigt,zeeman_faraday)
								endif
							endif

							call zeeman_opacity(model(k),i,zeeman_voigt,zeeman_faraday,ki_partial,kq_partial,&
								ku_partial,kv_partial,fq_partial,fu_partial,fv_partial)

							ki = ki + ki_partial * model(k)%filling_factor * model(j)%kl(i)
							kq = kq + kq_partial * model(k)%filling_factor * model(j)%kl(i)
							ku = ku + ku_partial * model(k)%filling_factor * model(j)%kl(i)
							kv = kv + kv_partial * model(k)%filling_factor * model(j)%kl(i)
							fq = fq + fq_partial * model(k)%filling_factor * model(j)%kl(i)
							fu = fu + fu_partial * model(k)%filling_factor * model(j)%kl(i)
							fv = fv + fv_partial * model(k)%filling_factor * model(j)%kl(i)
						enddo
					endif
				enddo

				delta = (1.d0+ki)**4 + (1.d0+ki)**2 * (fq**2+fu**2+fv**2-kq**2-ku**2-kv**2) - (kq*fq+ku*fu+kv*fv)**2
		
				stokes(1,:) = factor1 * (1.d0+model(j)%beta*model(j)%mu*(1.d0+ki) / delta * &
					((1.d0+ki)**2 + fq**2 + fu**2 + fv**2))
				stokes(2,:) = factor2 / delta * ((1.d0+ki)**2*kq - (1.d0+ki)*(ku*fv-kv*fu) + fq*(kq*fq+ku*fu+kv*fv))
				stokes(3,:) = factor2 / delta * ((1.d0+ki)**2*ku - (1.d0+ki)*(kv*fq-kq*fv) + fu*(kq*fq+ku*fu+kv*fv))
				stokes(4,:) = factor2 / delta * ((1.d0+ki)**2*kv - (1.d0+ki)*(kq*fu-ku*fq) + fv*(kq*fq+ku*fu+kv*fv))
					
				if (spectro_or_filter == 1) then
					Stokes_Syn_HR(1)%stokes = Stokes_Syn_HR(1)%stokes + model(j)%filling_factor * stokes
				else
					Stokes_Syn%stokes = Stokes_Syn%stokes + model(j)%filling_factor * stokes
				endif
				
			endif
			
		enddo

! Apply filters
		if (spectro_or_filter == 1) then
			do i = 1, 4
				do j = 1, filter_profiles%nfilters
					Stokes_Syn%stokes(i,j) = &
						int_tabulated(Stokes_Syn_HR(1)%lambda, filter_profiles%transmission(j,:)*Stokes_Syn_HR(1)%stokes(i,:))
				enddo
			enddo
		endif


		deallocate(zeeman_voigt)
		deallocate(zeeman_faraday)
		deallocate(ki)
		deallocate(kq)
		deallocate(ku)
		deallocate(kv)
		deallocate(fq)
		deallocate(fu)
		deallocate(fv)
		
		deallocate(ki_partial)
		deallocate(kq_partial)
		deallocate(ku_partial)
		deallocate(kv_partial)
		deallocate(fq_partial)
		deallocate(fu_partial)
		deallocate(fv_partial)
		
		deallocate(stokes)
		deallocate(delta)

		if (using_instrumental_profile) then
			call ApplyFilter(Stokes_Syn, instrumental_profile)			
		endif

	end subroutine synthesize
		
end module atomic_functions
