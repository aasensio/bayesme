module hyperfine
use math_functions
implicit none
contains

! ---------------------------------------------------------
! Calculate the sizes of the eigenvalues and eigenvectors
! ---------------------------------------------------------
	subroutine hyper_size(J, I, nflargest, f2max)
	real(kind=8) :: J, I
	integer :: nflargest, F2min, F2max, m2, f2a, f2b, nf
		
		F2min = abs(2*J-2*I)
	 	F2max = 2*J+2*I
		
		nflargest = 0
	 	do m2 = -f2max, f2max, 2
			f2a = abs(m2)
	 	  	if (f2a <= f2min) f2a = f2min			
	 	  	f2b = f2max
	 	  	nf = 1 + (f2b-f2a) / 2.d0
		  	if (nf > nflargest) then 
				nflargest = nf
			endif
		enddo
				
	end subroutine hyper_size
	
! ---------------------------------------------------------
! Calculate the sizes of the eigenvalues and eigenvectors
! ---------------------------------------------------------
	subroutine hyper_components(Ju, Jl, I, f2m_up, f2m_low, nflev_up, nflev_low, nPi, nLeft, nRight)
	real(kind=8) :: Ju, Jl, I
	integer :: nPi, nLeft, nRight, m2u, m2l, f2a_up, f2a_low, iu, il, f2m_up, f2m_low
	integer :: f2min_up, f2max_up, f2min_low, f2max_low
	integer :: nflev_up(0:2*f2m_up), nflev_low(0:2*f2m_low)
	real(kind=8) :: deltaM, q
		
		f2min_up = abs(2*Ju-2*I)
	 	f2max_up = 2*Ju+2*I
	 
	 	f2min_low = abs(2*Jl-2*I)
	 	f2max_low = 2*Jl+2*I
	 	 
	 
! First count the number of hyperfine components
		nPi = 0
		nRight = 0
		nLeft = 0
	 	do m2u = -f2max_up, f2max_up, 2
	 	  
		  f2a_up = abs(m2u)
	 	  if (f2a_up <= f2min_up) f2a_up = f2min_up
		  
	 	  do m2l = -f2max_low, f2max_low, 2
		   					
		   	f2a_low = abs(m2l)
	 	   	if (f2a_low <= f2min_low) f2a_low = f2min_low
		   	
		   	deltaM = (m2u - m2l) / 2.d0
				q = -deltaM
				if (abs(deltaM) <= 1) then					
					do iu = 0, nflev_up(f2max_up+m2u)-1
						do il = 0, nflev_low(f2max_low+m2l)-1
							if (q == 0) nPi = nPi + 1
							if (q == -1) nLeft = nLeft + 1
							if (q == 1) nRight = nRight + 1
					 	enddo
					enddo
				endif
			enddo
		enddo	
				
	end subroutine hyper_components
	
! ---------------------------------------------------------
! Calculate the hyperfine eigenvalues and eigenvectors
! ---------------------------------------------------------
	subroutine hyper_diagonalize(L, S, J, I, A, B, Bfield, nflev, energy, c, f2m, nfl)
	integer :: f2m, nfl
	real(kind=8) :: L, S, J, I, A, B, Bfield, energy(0:2*f2m,0:nfl-1), c(0:2*f2m,0:nfl-1,0:f2m)
	
	integer :: nflev(0:2*f2m), F2min, F2max, m2, f2a, f2b, nu, nl, loop, nf, l1, l2
	real(kind=8) :: gLS, mu0, Fu, Fl, t1, t1_diag, t2_diag, K, M
	real(kind=8), allocatable :: H(:,:)
	real(kind=8), allocatable :: d(:), ud(:), z(:,:)
		
		gLS = 1.d0 + 0.5d0*(J*(J+1)+S*(S+1)-L*(L+1)) / (J*(J+1))
	 
	 	F2min = abs(2*J-2*I)
	 	F2max = 2*J+2*I
		
	 	mu0 = 9.27d-21 / (6.626d-27*3.d10)
	 	 
		do m2 = -f2max, f2max, 2
			M = m2 / 2.d0
		  	f2a = abs(m2)
		  	if (f2a <= f2min) f2a = f2min
		  	f2b = f2max
		  	nf = 1 + (f2b-f2a) / 2.d0
		  	nflev(f2max+m2) = nf
			
		  	allocate(H(nf,nf))

	 	  	do nu = 1, nf
				Fu = (f2a + 2*(nu-1)) / 2.d0
				do nl = 1, nf
					Fl = (f2a + 2*(nl-1)) / 2.d0
					
! Non-diagonal part
					t1 = sqrt(J*(J+1)*(2*J+1)*(2*Fu+1)*(2*Fl+1)) * &
						w6js(int(2*Fl),int(2*Fu),2,int(2*J),int(2*J),int(2*I)) * &
						w3js(int(2*Fu),int(2*Fl),2,-int(2*M),int(2*M),0)

					H(nu,nl) = mu0*Bfield*gLS*(-1.d0)**(J+I-M) * t1
					 
! Diagonal part
	 	   		if (fu == fl) then
						K = Fu*(Fu+1) - I*(I+1) - J*(J+1)
						t1_diag = 0.5d0 * A * K
						t2_diag = 0.d0
						if (I /= 0 .and. I /= 0.5d0 .and. J /= 0 .and. J /= 0.5d0) then
							t2_diag = 0.5d0 * B * (0.75d0*K*(K+1)-I*(I+1)*J*(J+1)) / (I*(2*I-1)*J*(2*J-1))
						endif
						H(nu,nl) = H(nu,nl) + t1_diag + t2_diag
					endif

				enddo
			enddo
			
			

! The matrix is tridiagonal. We use an appropriate routine for the diagonalization
			allocate(d(nf))
			d = 0.d0
			allocate(ud(nf))
			ud = 0.d0
			allocate(z(nf,nf))
			
			z = 0.d0
			do loop = 1, nf
				z(loop,loop) = 1.d0
			enddo
			
			do loop = 1, nf
				d(loop) = H(loop,loop)
			enddo
			
			do loop = 1, nf-1
				ud(loop+1) = H(loop+1,loop)
			enddo
						
			call tqli(d,ud,nf,nf,z)
			
! Store the eigenvalues and eigenvectors		  
			do l1 = 0, nf-1
				energy(f2max+m2,l1) = d(l1+1)
				do l2 = 0, nf-1
					Fu = f2a + 2*l2
					c(f2max+m2,l1,Fu) = z(l2+1,l1+1)
				enddo
			enddo
			
			deallocate(d)
			deallocate(ud)
			deallocate(z)
			deallocate(H)
		enddo

	end subroutine hyper_diagonalize
	
! ---------------------------------------------------------
! Calculate the hyperfine pattern
! ---------------------------------------------------------
	subroutine hyper_pattern(I, Lu, Su, Ju, Au, Bu, Ll, Sl, Jl, Al, Bl, Bfield, &
		nflev_up, energy_up, c_up, f2m_up, nfl_up, nflev_low, energy_low, c_low, f2m_low, nfl_low, &
		SplitPi, SplitRight, SplitLeft, StrPi, StrRight, StrLeft)
		
	integer :: f2m_up, nfl_up, f2m_low, nfl_low
	real(kind=8) :: I, Lu, Su, Ju, Au, Bu, Ll, Sl, Jl, Al, Bl
	real(kind=8) :: Bfield
	real(kind=8) :: energy_up(0:2*f2m_up,0:nfl_up-1), energy_low(0:2*f2m_low,0:nfl_low-1)
	real(kind=8) :: c_up(0:2*f2m_up,0:nfl_up-1,0:f2m_up), c_low(0:2*f2m_low,0:nfl_low-1,0:f2m_low)
	integer :: nflev_up(0:2*f2m_up), nflev_low(0:2*f2m_low)
	
	integer :: f2min_up, f2max_up, f2min_low, f2max_low, ncomponents, m2u, m2l, f2a_up, f2a_low, iu, il
	integer :: nPi, nRight, nLeft, F2, Fp2, Fs2, Ft2
	real(kind=8) :: deltaM, q, splitting, sum, t1, t2, t3, t4
	
	real(kind=8) :: SplitPi(:), SplitRight(:), SplitLeft(:), StrPi(:), StrRight(:), StrLeft(:)
			
		
		f2min_up = abs(2*Ju-2*I)
	 	f2max_up = 2*Ju+2*I
	 
	 	f2min_low = abs(2*Jl-2*I)
	 	f2max_low = 2*Jl+2*I	 
		
		nPi = 0
		nRight = 0
		nLeft = 0
	 	do m2u = -f2max_up, f2max_up, 2
	 	  
		  f2a_up = abs(m2u)
	 	  if (f2a_up <= f2min_up) f2a_up = f2min_up
		  
	 	  do m2l = -f2max_low, f2max_low, 2
		   					
		   	f2a_low = abs(m2l)
	 	   	if (f2a_low <= f2min_low) f2a_low = f2min_low
		   	
		   	deltaM = (m2u - m2l) / 2.d0
				q = -deltaM
				if (abs(deltaM) <= 1) then					
					do iu = 0, nflev_up(f2max_up+m2u)-1
						do il = 0, nflev_low(f2max_low+m2l)-1
							
							splitting = energy_up(f2max_up+m2u,iu) - energy_low(f2max_low+m2l,il)
							
							sum = 0.d0
							
	 	   		 	   do F2 = f2a_low, f2max_low, 2
					 	   		 
								do Fs2 = f2a_low, f2max_low, 2
										  
									do Fp2 = f2a_up, f2max_up, 2
					 	   		 	   	
										do Ft2 = f2a_up, f2max_up, 2
													 
											t1 = c_low(f2max_low+m2l,il,F2) * c_low(f2max_low+m2l,il,Fs2) * &
												c_up(f2max_up+m2u,iu,Fp2) * c_up(f2max_up+m2u,iu,Ft2)
											t2 = 3.d0 / (2.d0*I+1.d0) * sqrt((F2+1.d0)*(Fp2+1.d0)*(Fs2+1.d0)*(Ft2+1.d0))
											t3 = w6js(int(2*Jl),int(2*Ju),2,int(Fp2),int(F2),int(2*I)) * &
												w6js(int(2*Jl),int(2*Ju),2,int(Ft2),int(Fs2),int(2*I))
											t4 = w3js(int(Fp2),int(F2),2,-int(m2u),int(m2l),-int(2*q)) * &
												w3js(int(Ft2),int(Fs2),2,-int(m2u),int(m2l),-int(2*q))
													 
											sum = sum + t1 * t2 * t3 * t4											
										enddo
									enddo
								enddo
							enddo							
							
							if (q == 0) then
								nPi = nPi + 1
								SplitPi(nPi) = splitting
								StrPi(nPi) = sum
							endif
							
							if (q == -1) then
								nLeft = nLeft + 1
								SplitLeft(nLeft) = splitting
								StrLeft(nLeft) = sum								
							endif
							
							if (q == 1) then
								nRight = nRight + 1
								SplitRight(nRight) = splitting
								StrRight(nRight) = sum								
							endif
							
							
					 	enddo
					enddo
				endif
			enddo
		enddo		
		
	end subroutine hyper_pattern

end module hyperfine
