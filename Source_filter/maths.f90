module math_vars
implicit none
	integer, parameter :: nfac = 301
	real(kind=8) :: fact(0:nfac)
	real(kind=8) :: identity(4,4)
	
	real(kind=8), parameter :: L_weideman = 2.2248028d0
	real(kind=8), parameter :: a_weideman(7) = (/1.0547115,0.60504614,0.20181603,0.0080914418,-0.021688188,-0.0038423489,0.0030903311/) 
end module math_vars


module math_functions
use math_vars
use vars
use Fast_Fourier, only : fft_forward, fft_backward
implicit none
contains

!-------------------------------------------------------------
! Carry out the Cholesky decomposition of a symmetric matrix
!-------------------------------------------------------------
	subroutine cholesky(a,n,p)
	integer :: n
	real(kind=8) :: a(n,n), p(n)
	integer :: i, j, k
	real(kind=8) :: sum
	
		do i = 1, n
			do j = i, n
				sum = a(i,j)
				do k = i-1, 1, -1
					sum = sum - a(i,k)*a(j,k)
				enddo
				if (i == j) then
					if (sum == 0.d0) then
						print *, 'Cholesky decomposition failed...'												
					endif
					p(i) = dsqrt(sum)
				else
					a(j,i) = sum / p(i)
				endif
			enddo
		enddo
						
	end subroutine cholesky

!-------------------------------------------------------------
! Diagonalizes a symmetric matrix
!-------------------------------------------------------------
	subroutine jacobi(a,n,d,v,nrot)
	integer :: n, nrot
	real(kind=8) :: a(n,n), d(n), v(n,n)      
	integer :: i,ip,iq,j
	real(kind=8) :: c,g,h,s,sm,t,tau,theta,tresh,b(n),z(n)
	
		do ip=1,n
			do iq=1,n
				v(ip,iq)=0.d0
			enddo			
        	v(ip,ip)=1.d0
		enddo
		
      do ip=1,n
			b(ip)=a(ip,ip)
        	d(ip)=b(ip)
        	z(ip)=0.d0
		enddo
		
      nrot = 0
      do i=1,50
			sm=0.d0
        	do ip=1,n-1
          	do iq=ip+1,n
            	sm=sm+abs(a(ip,iq))
				enddo
			enddo
        	if(sm == 0.d0) return
        	if(i < 4)then
          	tresh=0.2*sm/n**2
        	else
          	tresh=0.
        	endif
        	
        	do ip=1,n-1
          	do iq=ip+1,n
            	g=100.*abs(a(ip,iq))
            	if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
              		a(ip,iq)=0.
            	else if(abs(a(ip,iq)).gt.tresh) then
              		h=d(iq)-d(ip)
              	if(abs(h)+g.eq.abs(h))then
                	t=a(ip,iq)/h
              	else
                	theta=0.5*h/a(ip,iq)
                	t=1./(abs(theta)+sqrt(1.+theta**2))
                	if(theta.lt.0.)t=-t
              	endif
              	c=1./sqrt(1+t**2)
              	s=t*c
              	tau=s/(1.+c)
              	h=t*a(ip,iq)
              	z(ip)=z(ip)-h
              	z(iq)=z(iq)+h
              	d(ip)=d(ip)-h
              	d(iq)=d(iq)+h
              	a(ip,iq)=0.
              	do j=1,ip-1
                	g=a(j,ip)
                	h=a(j,iq)
                	a(j,ip)=g-s*(h+g*tau)
                	a(j,iq)=h+s*(g-h*tau)
					enddo
              	do j=ip+1,iq-1
                	g=a(ip,j)
                	h=a(j,iq)
                	a(ip,j)=g-s*(h+g*tau)
                	a(j,iq)=h+s*(g-h*tau)
					enddo
              	do j=iq+1,n
                	g=a(ip,j)
                	h=a(iq,j)
                	a(ip,j)=g-s*(h+g*tau)
                	a(iq,j)=h+s*(g-h*tau)
					enddo
              	do j=1,n
                	g=v(j,ip)
                	h=v(j,iq)
                	v(j,ip)=g-s*(h+g*tau)
                	v(j,iq)=h+s*(g-h*tau)
					enddo
              	nrot=nrot+1
            endif
			enddo
			enddo
        	do ip=1,n
          	b(ip)=b(ip)+z(ip)
          	d(ip)=b(ip)
          	z(ip)=0.
			enddo
		enddo
      print *, 'too many iterations in jacobi'
      stop 
      return
      end subroutine jacobi
      
      
!-------------------------------------------------------------
! Initialize the random number generator
!-------------------------------------------------------------
	subroutine init_random_seed()
	integer :: i, n, clock
   integer, dimension(:), allocatable :: seed
          
   	call random_seed(size = n)
      allocate(seed(n))
          
      call system_clock(count=clock)
          
      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      call random_seed(put = seed)
          
      deallocate(seed)
   end subroutine init_random_seed
          
!-------------------------------------------------------------
! Generates a random number following an uniform distribution in the interval [0,1]
! Call it with idum<0 for initialization
!-------------------------------------------------------------
	function randomu()
	real(kind=8) :: randomu
	
		call random_number(randomu)
							
	end function randomu
	
!-------------------------------------------------------------
! Generates a random number following an normal distribution with zero mean
! and unit variance
!-------------------------------------------------------------
	function randomn()

	real(kind=8) :: randomn

	real(kind=8) :: u, sum
	real(kind=8), save :: v, sln
	logical, save :: second = .false.
	real(kind=8), parameter :: one = 1.0, vsmall = tiny( one )

! if second, use the second random number generated on last call
	if (second) then

		second = .false.
  		randomn = v*sln
	else
! first call; generate a pair of random normals

  		second = .true.
  		do
    		call random_number( u )
    		call random_number( v )
    		u = scale( u, 1 ) - one
    		v = scale( v, 1 ) - one
    		sum = u*u + v*v + vsmall         ! vsmall added to prevent log(zero) / zero
    		if(sum < one) exit
  		end do
  		sln = sqrt(- scale( log(sum), 1 ) / sum)
  		randomn = u*sln
	end if

	return
	end function randomn
	
!-------------------------------------------------------------
! Generates a random number following an normal distribution with zero mean
! and unit variance
! Call it with idum<0 for initialization
!-------------------------------------------------------------
!	function randomn()
!	real(kind=8) :: randomn, gasdev, fac, gset, rsq, v1, v2, ran1
!	integer :: iset
!	save :: iset, gset		
!		if (iset == 0) then
!			rsq = 0.0
!			do while(rsq > 1 .or. rsq == 0.0)
!				v1 = 2.0*randomu()-1.0
!				v2 = 2.0*randomu()-1.0
!				rsq = v1**2+v2**2
!			enddo
!			fac = sqrt(-2.0*log(rsq)/rsq)
!			gset = v1*fac
!			randomn = v2*fac
!			iset = 1
!		else
!			randomn = gset
!			iset = 0
!		endif
!	end function randomn
	
!-------------------------------------------------------------
! Generates a multivariate normal random number with a given
! mean and covariance matrix
!-------------------------------------------------------------
	function mrandomn(idum,rmean,covar)
	integer :: idum, n, i, j
	real(kind=8) :: rmean(:), covar(:,:), mrandomn(size(rmean))
	real(kind=8) :: chol(size(rmean),size(rmean)), p(size(rmean)), eps(size(rmean))
	
		n = size(rmean)
				
		chol = covar
		
		do i = 1, n
			eps(i) = randomn()
			chol(i,i) = chol(i,i) + 1.d-7    ! Some regularization
		enddo
								
		call cholesky(chol,n,p)
										
		do j = 1, n			
			do i = j, n
				chol(j,i) = 0.d0
			enddo
			chol(j,j) = p(j)
		enddo
								
		mrandomn = matmul(chol,eps) + rmean			
				
		
	end function mrandomn

!-------------------------------------------------------------------
! Calculate the natural logarithm of the Gamma function
!-------------------------------------------------------------------
	function gammln(xx)
	real(kind=8) :: gammln, xx
	integer :: j 
	real(kind=8) :: ser,tmp,x,y
	real(kind=8) :: cof(6) = (/76.18009172947146d0,-86.50532032941677d0, 24.01409824083091d0,-1.231739572450155d0,&
		.1208650973866179d-2,-.5395239384953d-5/)
	real(kind=8) :: stp = 2.5066282746310005d0
		
		x = xx 
		y = x 
		tmp = x+5.5d0 
		tmp = (x+0.5d0)*log(tmp)-tmp 
		ser = 1.000000000190015d0 
		do j=1,6 
			y = y+1.d0 
			ser = ser + cof(j)/y 
		enddo 
		gammln = tmp + log(stp*ser/x)
	end function gammln

!-------------------------------------------------------------------
! Read n lines at the selected unit 
!-------------------------------------------------------------------
	subroutine lb(u, n)
	integer :: u, n
	integer :: i
		do i = 1, n
	 		read(u,*) 
		enddo
	end subroutine lb

!-----------------------------------------------------------------
! Returns the weights (w) and the abscissas (x) for a Gaussian integration using the 
! Gauss-Legendre formula, using n points
!-----------------------------------------------------------------
   subroutine gauleg(x1,x2,x,w,n)
   integer, INTENT(IN) :: n
   real(kind=8), INTENT(IN) :: x1,x2
   real(kind=8), INTENT(INOUT) :: x(n),w(n)
   real(kind=8), parameter :: eps = 3.d-14
   integer :: i,j,m
   real(kind=8) :: p1,p2,p3,pp,xl,xm,z,z1

   	m = (n+1)/2
		xm = 0.5d0*(x2+x1)
		xl = 0.5d0*(x2-x1)
		do i=1,m
   		z = dcos(PI*(i-.25d0)/(n+.5d0))
			z1 = z + 2.d0*EPS
			do while(abs(z-z1) > EPS)

   			p1 = 1.d0
   			p2 = 0.d0
   			do j = 1, n
					p3 = p2
   				p2 = p1
   				p1 = ((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
				enddo
   			pp = n*(z*p1-p2)/(z*z-1.d0)
   			z1 = z
   			z = z1-p1/pp
			enddo
   		x(i) = xm-xl*z
   		x(n+1-i) = xm+xl*z
   		w(i) = 2.d0*xl/((1.d0-z*z)*pp*pp)
   		w(n+1-i) = w(i)
   	enddo

   end subroutine gauleg
!-------------------------------------------------------------------
! Convert vacuum wavelengths to air wavelengths
!       Corrects for the index of refraction of air under standard conditions.  
!       Wavelength values below 2000 A will not be altered.  Accurate to 
!       about 0.005 A 
!       It uses an approximation to the 4th power of inverse wavenumber is used
!       See IUE Image Processing Manual   Page 6-15.
!	wave : wavelength in A
!-------------------------------------------------------------------
	function vactoair(wave)
  	real(kind=8) :: vactoair, wave, wave2, factor
		
		wave2 = wave*wave
  		factor = 1.d0 + 2.735182d-4 + 131.4182d0/wave2 + 2.76249d8/(wave2*wave2)
		if (wave <= 2000.d0) factor = 1.d0

  		vactoair = wave / factor

	end function vactoair
	
!-------------------------------------------------------------------
! Convert air wavelengths to vacuum wavelengths
!       Corrects for the index of refraction of air under standard conditions.  
!       Wavelength values below 2000 A will not be altered.  Accurate to 
!       about 0.005 A 
!       It uses an approximation to the 4th power of inverse wavenumber is used
!       See IUE Image Processing Manual   Page 6-15.
!		  Morton (ApJS 77, 119)
!	wave : wavelength in A
!-------------------------------------------------------------------
	function airtovac(wave)
  	real(kind=8) :: airtovac, wave, sigma2, factor
	
		sigma2 = (1.d4/wave)**2
		
  		factor = 1.d0 + 6.4328d-5 + 2.94981d-2/(146.d0-sigma2) + 2.5540d-4/(41.d0-sigma2)
		if (wave <= 2000.d0) factor = 1.d0

  		airtovac = wave * factor

	end function airtovac

!-------------------------------------------------------------------
! Returns which is the index of the array wave_total closer to wave
!-------------------------------------------------------------------
	function close_to(wave_total, wave)
	real(kind=8) :: wave_total(:), wave
	real(kind=8), allocatable :: diff(:)
	integer :: n, i, location(1)
	integer :: close_to
		
		n = size(wave_total)
		allocate(diff(n))
		diff = dabs(wave_total-wave)
		location = minloc(diff)
		deallocate(diff)
		close_to = location(1)
		
	end function close_to
	
!-------------------------------------------------------------------
! Returns which is the index of the array wave_total closer to wave
!-------------------------------------------------------------------
	function extremes(freq, line_frequency, delta)
	real(kind=8) :: freq(:), line_frequency, delta
	integer :: left, right, extremes(2)
		

		left = close_to(freq, line_frequency-delta)
		right = close_to(freq, line_frequency+delta)
		
		extremes(1) = min(left,right)
		extremes(2) = max(left,right)
				
	end function extremes	
	
!-------------------------------------------------------------------
! Returns true if the wavelength we give is inside the limits of the spectral regions used
!-------------------------------------------------------------------
	function in_range(initial_wavelength, final_wavelength, wave)
	real(kind=8) :: initial_wavelength(:), final_wavelength(:), wave
	integer :: n, i
	logical :: in_range
		in_range = .FALSE.
		n = size(initial_wavelength)
		i = 1
		do while(.not.in_range .and. i <= n)
			if (wave >= initial_wavelength(i) .and. wave <= final_wavelength(i)) then
				in_range = .TRUE.
			endif
			i = i + 1
		enddo
	end function in_range

! ---------------------------------------------------------
! Returns Planck's function for a frequency and a set of temperatures in cgs
! INPUT:
!     - nu : frequency
!     - T : temperature
! ---------------------------------------------------------               
   function planck_vect(nu,T)
	real*8, INTENT(IN) :: nu, T(:)
   real*8 :: planck_vect(size(T))
   real*8 :: temporal(size(T))

		where (T /= 0.d0)
			planck_vect = PHC * nu**3.d0 / (dexp(PHK * nu / T) - 1.d0)
		elsewhere
			planck_vect = 0.d0
		endwhere

   end function planck_vect

! ---------------------------------------------------------
! Returns Planck's function for a frequency and temperature in cgs
! INPUT:
!     - nu : frequency (Hz)
!     - T : temperature (K)
! ---------------------------------------------------------               
   function planck(nu,T)
   real*8 :: planck
   real*8, INTENT(IN) :: nu, T
   real*8 :: temporal
      if (T /= 0.d0) then
         planck = PHC * nu**3.d0 / (dexp(PHK * nu / T) - 1.d0)
      else
			planck = 0.d0
      endif   

   end function planck
	
! ---------------------------------------------------------
! Returns Planck's function for a wavelength and a set of temperatures in cgs
! INPUT:
!     - lambda : wavelength
!     - T : temperature
! ---------------------------------------------------------               
   function planck_vect_lambda(lambda,T)
	real*8, INTENT(IN) :: lambda, T(:)
   real*8 :: planck_vect_lambda(size(T))
   real*8 :: temporal(size(T))

		where (T /= 0.d0)
			planck_vect_lambda = PHC2 / lambda**5.d0 / (dexp(PHK * PC / (lambda * T)) - 1.d0)
		elsewhere
			planck_vect_lambda = 0.d0
		endwhere

   end function planck_vect_lambda

! ---------------------------------------------------------
! Returns Planck's function for a wavelength and temperature in cgs
! INPUT:
!     - lambda : wavelength (cm)
!     - T : temperature (K)
! ---------------------------------------------------------               
   function planck_lambda(lambda,T)
   real*8 :: planck_lambda
   real*8, INTENT(IN) :: lambda, T
   real*8 :: temporal

      if (T /= 0.d0) then
         temporal = dexp(PHK * PC / (lambda * T) - 1.d0)
         planck_lambda = PHC2 / lambda**5.d0 / temporal
      else
			planck_lambda = 0.d0
      endif   

   end function planck_lambda	
	
! ---------------------------------------------------------
! Returns the derivative of the Planck's function for a wavelength and a set of temperatures in cgs
! INPUT:
!     - lambda : wavelength
!     - T : temperature
! ---------------------------------------------------------               
   function planck_vect_lambda_deriv(lambda,T)
	real*8, INTENT(IN) :: lambda, T(:)
   real*8 :: planck_vect_lambda_deriv(size(T))
   real*8 :: temporal(size(T))

		where (T /= 0.d0)
			temporal = dexp(PHK * PC / (lambda * T))
			planck_vect_lambda_deriv = PH2C3 / (PK*T**2*lambda**6.d0) / (temporal-1.d0) * temporal
		elsewhere
			planck_vect_lambda_deriv = 0.d0
		endwhere

   end function planck_vect_lambda_deriv

! ---------------------------------------------------------
! Returns the derivative of the Planck's function for a wavelength and temperature in cgs
! INPUT:
!     - lambda : wavelength (cm)
!     - T : temperature (K)
! ---------------------------------------------------------               
   function planck_lambda_deriv(lambda,T)
   real*8 :: planck_lambda_deriv
   real*8, INTENT(IN) :: lambda, T
   real*8 :: temporal

      if (T /= 0.d0) then
         temporal = dexp(PHK * PC / (lambda * T))
         planck_lambda_deriv = PH2C3 / (PK*T**2*lambda**6.d0) / (temporal-1.d0) * temporal
      else
			planck_lambda_deriv = 0.d0
      endif   

   end function planck_lambda_deriv

! ---------------------------------------------------------
! Diagonalizes a tridiagonal matrix
! ---------------------------------------------------------
	subroutine tqli(d,e,n,np,z)
	integer :: n, np, i, l, m, iter, k
	real(kind=8) :: d(np),e(np),z(np,np), g, r, s, c, p, f, b, dd
      
		if (n.gt.1) then
      do 11 i=2,n
      e(i-1)=e(i)
 11   continue
      e(n)=0.d0
      do 15 l=1,n
      iter=0
 1    do 12 m=l,n-1
      dd=dabs(d(m))+dabs(d(m+1))
      if(dabs(e(m))+dd.eq.dd) go to 2
 12   continue
      m=n
 2    if(m.ne.l) then
      if(iter.eq.30) then
	   print*, 'too many iterations'
	   stop
      endif
      iter=iter+1
      g=(d(l+1)-d(l))/(2.d0*e(l))
      r=dsqrt(g**2+1.d0)
      g=d(m)-d(l)+e(l)/(g+dsign(r,g))
      s=1.d0
      c=1.d0
      p=0.d0
      do 14 i=m-1,l,-1
      f=s*e(i)
      b=c*e(i)
      if(dabs(f).ge.dabs(g)) then
      c=g/f
      r=dsqrt(c**2+1.d0)
      e(i+1)=f*r
      s=1.d0/r
      c=c*s
      else
      s=f/g
      r=dsqrt(s**2+1.d0)
      e(i+1)=g*r
      c=1.d0/r
      s=s*c
      endif
      g=d(i+1)-p
      r=(d(i)-g)*s+2.d0*c*b
      p=s*r
      d(i+1)=g+p
      g=c*r-b
      do 13 k=1,n
      f=z(k,i+1)
      z(k,i+1)=s*z(k,i)+c*f
      z(k,i)=c*z(k,i)-s*f
 13   continue
 14   continue
      d(l)=d(l)-p
      e(l)=g
      e(m)=0.d0
      go to 1
      endif
 15   continue
      endif
      return
      end subroutine tqli		
		
!----------------------------------------------------------------
! Calculates the factorial of the integers up to 301 and save it into the fact array
!----------------------------------------------------------------
   subroutine factrl
	integer :: i
      fact(0) = 1.d0
      do i=1,301
			fact(I) = fact(I-1) * dble(I)
		enddo
   END subroutine factrl

!----------------------------------------------------------------
! This function calculates the 3-j symbol
! J_i and M_i have to be twice the actual value of J and M
!----------------------------------------------------------------
   function w3js(j1,j2,j3,m1,m2,m3)
   integer :: m1, m2, m3, j1, j2, j3
	integer :: ia, ib, ic, id, ie, im, ig, ih, z, zmin, zmax, jsum
	real(kind=8) :: w3js, cc, denom, cc1, cc2

!		w3js = w3js_regge(j1/2,j2/2,j3/2,m1/2,m2/2,m3/2)
!		return
      w3js = 0.d0
      if (m1+m2+m3 /= 0) return
      ia = j1 + j2
      if (j3 > ia) return
      ib = j1 - j2
      if (j3 < abs(ib)) return
		if (abs(m1) > j1) return
		if (abs(m2) > j2) return
		if (abs(m3) > j3) return
		
      jsum = j3 + ia
      ic = j1 - m1
      id = j2 - m2
      
		ie = j3 - j2 + m1
		im = j3 - j1 - m2
		zmin = max0(0,-ie,-im)
		ig = ia - j3
		ih = j2 + m2
		zmax = min0(ig,ih,ic)
		cc = 0.d0

		do z = zmin, zmax, 2
			denom = fact(z/2)*fact((ig-z)/2)*fact((ic-z)/2)*fact((ih-z)/2)*&
				fact((ie+z)/2)*fact((im+z)/2)
      	if (mod(z,4) /= 0) denom = -denom
			cc = cc + 1.d0/denom
		enddo

		cc1 = fact(ig/2)*fact((j3+ib)/2)*fact((j3-ib)/2)/fact((jsum+2)/2)
      cc2 = fact((j1+m1)/2)*fact(ic/2)*fact(ih/2)*fact(id/2)*fact((j3-m3)/2)*fact((j3+m3)/2)
      cc = cc * sqrt(1.d0*cc1*cc2)
		if (mod(ib-m3,4) /= 0) cc = -cc
		w3js = cc
		if (abs(w3js) < 1.d-8) w3js = 0.d0		
1000 		return
   end function w3js
	
	
!----------------------------------------------------------------
! This function calculates the 6-j symbol
! J_i and M_i have to be twice the actual value of J and M
!----------------------------------------------------------------
	function w6js(j1,j2,j3,l1,l2,l3)
	integer :: j1,j2,j3,l1,l2,l3
	integer :: ia, ib, ic, id, ie, iif, ig, ih, sum1, sum2, sum3, sum4
	integer :: w, wmin, wmax, ii, ij, ik
   real(kind=8) :: w6js, omega, theta1, theta2, theta3, theta4, theta, denom

      w6js = 0.d0
		ia = j1 + j2
		if (ia < j3) return
		ib = j1 - j2
		if (abs(ib) > j3) return
		ic = j1 + l2
		if (ic < l3) return
		id = j1 - l2
		if (abs(id) > l3) return
		ie = l1 + j2
		if (ie < l3) return
		iif = l1 - j2
		if (abs(iif) > l3) return
		ig = l1 + l2
		if (ig < j3) return
		ih = l1 - l2
		if (abs(ih) > j3) return
      sum1=ia + j3
      sum2=ic + l3
      sum3=ie + l3
      sum4=ig + j3
		wmin = max0(sum1, sum2, sum3, sum4)
		ii = ia + ig
		ij = j2 + j3 + l2 + l3
		ik = j3 + j1 + l3 + l1
		wmax = min0(ii,ij,ik)
		omega = 0.d0
		do w = wmin, wmax, 2
	      denom = fact((w-sum1)/2)*fact((w-sum2)/2)*fact((w-sum3)/2)&
				*fact((w-sum4)/2)*fact((ii-w)/2)*fact((ij-w)/2)&
				*fact((ik-w)/2)
			if (mod(w,4) /= 0) denom = -denom
			omega = omega + fact(w/2+1) / denom
		enddo	  
      theta1 = fact((ia-j3)/2)*fact((j3+ib)/2)*fact((j3-ib)/2)/fact(sum1/2+1)
      theta2 = fact((ic-l3)/2)*fact((l3+id)/2)*fact((l3-id)/2)/fact(sum2/2+1)
      theta3 = fact((ie-l3)/2)*fact((l3+iif)/2)*fact((l3-iiF)/2)/fact(sum3/2+1)
      theta4 = fact((ig-j3)/2)*fact((j3+ih)/2)*fact((j3-ih)/2)/fact(sum4/2+1)
      theta = theta1 * theta2 * theta3 * theta4
		w6js = omega * sqrt(theta)
      if (abs(w6js) < 1.d-8) w6js = 0.d0
1000 		return
   end function w6js

!--------------------------------------------------------------
! Generates Voigt and anomalous dispersion profiles using a different damping for each frequency
! point
! See J.A.C. Weideman, SIAM J. Numer. Anal. Vol 31, pp 1497
! http://dip.sun.ac.za/~weideman/research/cef.html
!  - Values for N=7
!		L=2.2248028
!		a=[1.0547115,0.60504614,0.20181603,0.0080914418,-0.021688188,-0.0038423489,0.0030903311]
!--------------------------------------------------------------
	function voigt_weideman(da, dv)
	real(kind=8) :: da(:)
	real(kind=8) :: dv(:), voigt_weideman(size(dv))
	complex :: z(size(dv)), p(size(dv)), ii, z2(size(dv))
	integer :: i, n

		n = size(dv)
		ii = cmplx(0.d0,1.d0)
		z = cmplx(dv,da)
		
		z2 = (L_weideman+ii*z)/(L_weideman-ii*z)
		
		p = a_weideman(1)+a_weideman(2)*z2+a_weideman(3)*z2**2+a_weideman(4)*z2**3+a_weideman(5)*z2**4+a_weideman(6)*z2**5+a_weideman(7)*z2**6
		
		p = 2.d0*p/(L_weideman-ii*z)**2 + (1.d0/sqrt(PI))/(L_weideman-ii*z)

   	voigt_weideman = dble(p)
		
	end function voigt_weideman
			
!--------------------------------------------------------------
! Generates Voigt and anomalous dispersion profiles using a different damping for each frequency
! point
! See J.A.C. Weideman, SIAM J. Numer. Anal. Vol 31, pp 1497
! http://dip.sun.ac.za/~weideman/research/cef.html
!  - Values for N=7
!		L=2.2248028
!		a=[1.0547115,0.60504614,0.20181603,0.0080914418,-0.021688188,-0.0038423489,0.0030903311]
!--------------------------------------------------------------
	function voigt_weideman_zeeman(da, dv)
	real(kind=8) :: da(:)
	real(kind=8) :: dv(:), voigt_weideman_zeeman(2,size(dv))
	complex :: z(size(dv)), p(size(dv)), ii, z2(size(dv))
	integer :: i, n

		n = size(dv)
		ii = cmplx(0.d0,1.d0)
		z = cmplx(dv,da)
		
		z2 = (L_weideman+ii*z)/(L_weideman-ii*z)
		
		p = a_weideman(1)+a_weideman(2)*z2+a_weideman(3)*z2**2+a_weideman(4)*z2**3+a_weideman(5)*z2**4+a_weideman(6)*z2**5+a_weideman(7)*z2**6
		
		p = 2.d0*p/(L_weideman-ii*z)**2 + (1.d0/sqrt(PI))/(L_weideman-ii*z)

   	voigt_weideman_zeeman(1,:) = dble(p)
		voigt_weideman_zeeman(2,:) = aimag(p)

	end function voigt_weideman_zeeman
		
!--------------------------------------------------------------
! Generates Voigt and anomalous dispersion profiles
! See Humlicek (1982) JQSRT 27, 437
!--------------------------------------------------------------
	function fvoigt(da, dv)
	real*8 :: da
	real*8 :: dv(:), fvoigt(size(dv))
	complex :: w4, z, t, u, v4
	real*8 :: s
	integer :: i, n
		
		n = size(dv)
		do i = 1, n
			z = cmplx(dv(i), da)
			t = cmplx(da, -dv(i))
			s = dabs(dv(i)) + da
			u = t*t


			if (s >= 15.d0) then
				w4 = t * 0.5641896d0 / (0.5d0+u)
			elseif (s >= 5.5) then
				w4 = t*(1.410474d0+u*0.5641896d0)/(0.75d0+u*(3.d0+u))
			elseif (da >= 0.195d0*dabs(dv(i))-0.176d0) then
				w4 = (16.4955d0+t*(20.20933d0+t*(11.96482d0+t*(3.778987d0+t*0.5642236d0)))) / &
					(16.4955d0+t*(38.82363d0+t*(39.27121d0+t*(21.69274d0+t*(6.699398d0+t)))))
			else 
				w4 = t*(36183.31d0-u*(3321.9905d0-u*(1540.787d0-u*(219.0313d0-u*(35.76683d0-&
					u*(1.320522d0-u*0.56419d0))))))
				v4 = (32066.6d0-u*(24322.84d0-u*(9022.228d0-u*(2186.181d0-u*(364.2191d0-&
					u*(61.57037d0-u*(1.841439d0-u)))))))
				w4 = cexp(u) - w4/v4
			endif
			fvoigt(i) = dble(w4)
		enddo                        

end function fvoigt

!--------------------------------------------------------------
! Generates Voigt and anomalous dispersion profiles
! See Humlicek (1982) JQSRT 27, 437
!--------------------------------------------------------------
	function fvoigt_zeeman(da, dv)
	real*8 :: da
	real*8 :: dv(:), fvoigt_zeeman(2,size(dv))
	complex :: w4, z, t, u, v4
	real*8 :: s
	integer :: i, n
		
		n = size(dv)
		do i = 1, n
			z = cmplx(dv(i), da)
			t = cmplx(da, -dv(i))
			s = dabs(dv(i)) + da
			u = t*t


			if (s >= 15.d0) then
				w4 = t * 0.5641896d0 / (0.5d0+u)
			elseif (s >= 5.5) then
				w4 = t*(1.410474d0+u*0.5641896d0)/(0.75d0+u*(3.d0+u))
			elseif (da >= 0.195d0*dabs(dv(i))-0.176d0) then
				w4 = (16.4955d0+t*(20.20933d0+t*(11.96482d0+t*(3.778987d0+t*0.5642236d0)))) / &
					(16.4955d0+t*(38.82363d0+t*(39.27121d0+t*(21.69274d0+t*(6.699398d0+t)))))
			else 
				w4 = t*(36183.31d0-u*(3321.9905d0-u*(1540.787d0-u*(219.0313d0-u*(35.76683d0-&
					u*(1.320522d0-u*0.56419d0))))))
				v4 = (32066.6d0-u*(24322.84d0-u*(9022.228d0-u*(2186.181d0-u*(364.2191d0-&
					u*(61.57037d0-u*(1.841439d0-u)))))))
				w4 = cexp(u) - w4/v4
			endif
			fvoigt_zeeman(1,i) = dble(w4)
			fvoigt_zeeman(2,i) = aimag(w4)
		enddo                        

	end function fvoigt_zeeman		
	
!--------------------------------------------------------------
! Generates Voigt and anomalous dispersion profiles using a different damping for each frequency
! point
! See Humlicek (1982) JQSRT 27, 437
!--------------------------------------------------------------
	function vecfvoigt(da, dv)
	real*8 :: da(:)
	real*8 :: dv(:), vecfvoigt(size(dv))
	complex :: w4, z, t, u, v4
	real*8 :: s
	integer :: i, n
                
		n = size(dv)
		do i = 1, n
! There is a problem for very small dampings. We force it to be 0 using a threshold                     
			if (da(i) < 1.d-3) da(i) = 0.d0
			z = cmplx(dv(i), da(i))
			t = cmplx(da(i), -dv(i))
			s = dabs(dv(i)) + da(i)
			u = t*t


			if (s >= 15.d0) then
				w4 = t * 0.5641896d0 / (0.5d0+u)
			elseif (s >= 5.5) then
				w4 = t*(1.410474d0+u*0.5641896d0)/(0.75d0+u*(3.d0+u))
			elseif (da(i) >= 0.195d0*dabs(dv(i))-0.176d0) then
				w4 = (16.4955d0+t*(20.20933d0+t*(11.96482d0+t*(3.778987d0+t*0.5642236d0)))) / &
					(16.4955d0+t*(38.82363d0+t*(39.27121d0+t*(21.69274d0+t*(6.699398d0+t)))))
			else 
				w4 = t*(36183.31d0-u*(3321.9905d0-u*(1540.787d0-u*(219.0313d0-u*(35.76683d0-&
					u*(1.320522d0-u*0.56419d0))))))
				v4 = (32066.6d0-u*(24322.84d0-u*(9022.228d0-u*(2186.181d0-u*(364.2191d0-&
					u*(61.57037d0-u*(1.841439d0-u)))))))
				w4 = cexp(u) - w4/v4
			endif
			vecfvoigt(i) = dble(w4)
		enddo                        

	end function vecfvoigt

!--------------------------------------------------------------
! Generates Voigt and anomalous dispersion profiles using a different damping for each frequency
! point
! See Humlicek (1982) JQSRT 27, 437
!--------------------------------------------------------------
	function vecfvoigt_zeeman(da, dv)
	real*8 :: da(:)
	real*8 :: dv(:), vecfvoigt_zeeman(2,size(dv))
	complex :: w4, z, t, u, v4
	real*8 :: s
	integer :: i, n
                
		n = size(dv)
		do i = 1, n
! There is a problem for very small dampings. We force it to be 0 using a threshold                     
			if (da(i) < 1.d-3) da(i) = 0.d0
			z = cmplx(dv(i), da(i))
			t = cmplx(da(i), -dv(i))
			s = dabs(dv(i)) + da(i)
			u = t*t


			if (s >= 15.d0) then
				w4 = t * 0.5641896d0 / (0.5d0+u)
			elseif (s >= 5.5) then
				w4 = t*(1.410474d0+u*0.5641896d0)/(0.75d0+u*(3.d0+u))
			elseif (da(i) >= 0.195d0*dabs(dv(i))-0.176d0) then
				w4 = (16.4955d0+t*(20.20933d0+t*(11.96482d0+t*(3.778987d0+t*0.5642236d0)))) / &
					(16.4955d0+t*(38.82363d0+t*(39.27121d0+t*(21.69274d0+t*(6.699398d0+t)))))
			else 
				w4 = t*(36183.31d0-u*(3321.9905d0-u*(1540.787d0-u*(219.0313d0-u*(35.76683d0-&
					u*(1.320522d0-u*0.56419d0))))))
				v4 = (32066.6d0-u*(24322.84d0-u*(9022.228d0-u*(2186.181d0-u*(364.2191d0-&
					u*(61.57037d0-u*(1.841439d0-u)))))))
				w4 = cexp(u) - w4/v4
			endif
			vecfvoigt_zeeman(1,i) = dble(w4)
			vecfvoigt_zeeman(2,i) = aimag(w4)
		enddo                        							

	end function vecfvoigt_zeeman

!--------------------------------------------------------------
! Generates Voigt and anomalous dispersion profiles using a different damping for each frequency
! point
! See Humlicek (1982) JQSRT 27, 437
!--------------------------------------------------------------
	function vecfvoigt_zeeman2(da, dv)
	real*8 :: da(:)
	real*8 :: dv(:), vecfvoigt_zeeman2(2,size(dv))
	complex :: w4(size(dv)), z(size(dv)), t(size(dv)), u(size(dv)), v4(size(dv))
	real*8 :: s(size(dv))
	integer :: i, n
                
		n = size(dv)
		
!		allocate(w4(n))
!		allocate(z(n))
!		allocate(t(n))
!		allocate(u(n))
!		allocate(v4(n))
!		allocate(s(n))
		
		where(da < 1.d-3) 
			da = 0.d0
		endwhere
		
		z = cmplx(dv,da)
		t = cmplx(da, -dv)
		s = dabs(dv) + da
		u = t*t
		
		where(s >= 15.d0)
			w4 = t * 0.5641896d0 / (0.5d0+u)
		endwhere
		where(s < 15.d0 .and. s >= 5.5d0)
			w4 = t*(1.410474d0+u*0.5641896d0)/(0.75d0+u*(3.d0+u))
		endwhere
		where(s < 15.d0 .and. s < 5.5d0 .and. da >= 0.195d0*dabs(dv)-0.176d0)
			w4 = (16.4955d0+t*(20.20933d0+t*(11.96482d0+t*(3.778987d0+t*0.5642236d0)))) / &
					(16.4955d0+t*(38.82363d0+t*(39.27121d0+t*(21.69274d0+t*(6.699398d0+t)))))
		endwhere
		where(s < 15.d0 .and. s < 5.5d0 .and. da < 0.195d0*dabs(dv)-0.176d0)
			w4 = t*(36183.31d0-u*(3321.9905d0-u*(1540.787d0-u*(219.0313d0-u*(35.76683d0-&
					u*(1.320522d0-u*0.56419d0))))))
			v4 = (32066.6d0-u*(24322.84d0-u*(9022.228d0-u*(2186.181d0-u*(364.2191d0-&
					u*(61.57037d0-u*(1.841439d0-u)))))))
			w4 = cexp(u) - w4/v4
		endwhere
		
		vecfvoigt_zeeman2(1,:) = dble(w4)
		vecfvoigt_zeeman2(2,:) = aimag(w4)
		
!		deallocate(w4)
!		deallocate(z)
!		deallocate(t)
!		deallocate(u)
!		deallocate(v4)
!		deallocate(s)

	end function vecfvoigt_zeeman2

	
! ---------------------------------------------------------
! Returns the mean of a vector
! ---------------------------------------------------------
      function mean(v)
		real*8 :: mean, v(:)
		integer :: n
			n = size(v)
			mean = sum(v) / n
		end function mean

! ---------------------------------------------------------
! Returns a Voigt profile for the line
! ---------------------------------------------------------
      function voigt(a,vv,j)
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      complex :: z
      real*8 :: xdws(28),ydws(28)

      data a0,a1,a2,a3,a4,a5,a6,b0,b1,b2,b3,b4,b5,b6/&
      122.607931777104326d0,214.382388694706425d0,181.928533092181549d0,&
      93.155580458138441d0,30.180142196210589d0,5.912626209773153d0,&
      .564189583562615d0,122.60793177387535d0,352.730625110963558d0,&
      457.334478783897737d0,348.703917719495792d0,170.354001821091472d0,&
      53.992906912940207d0,10.479857114260399d0/

      data xdws/.1d0,.2d0,.3d0,.4d0,.5d0,.6d0,.7d0,.8d0,.9d0,1.d0,1.2d0,1.4d0,1.6d0,1.8d0,2.d0,&
      3.d0,4.d0,5.d0,6.d0,7.d0,8.d0,9.d0,10.d0,12.d0,14.d0,16.d0,18.d0,20.d0/,ydws/&
      9.9335991d-02,1.9475104d-01,2.8263167d-01,3.5994348d-01,&
      4.2443639d-01,4.7476321d-01,5.1050407d-01,5.3210169d-01,&
      5.4072434d-01,5.3807950d-01,5.0727350d-01,4.5650724d-01,&
      3.9993989d-01,3.4677279d-01,3.0134040d-01,1.7827103d-01,&
      1.2934799d-01,1.0213407d-01,8.4542692d-02,7.2180972d-02,&
      6.3000202d-02,5.5905048d-02,5.0253846d-02,4.1812878d-02,&
      3.5806101d-02,3.1311397d-02,2.7820844d-02,2.5031367d-02/

      v=dabs(vv)
      IF(A.NE.0) GOTO 1 
      IF(J.NE.0) GOTO 3 
      VOIGT=DEXP(-V*V)   
      RETURN
   3  IF(V.GT.XDWS(1)) GOTO 4   
      D=V*(1.-.66666667d0*V*V)
      GOTO 8
   4  IF(V.GT.XDWS(28)) GOTO 5  
      K=27  
      DO 7 I=2,27   
      IF(XDWS(I).LT.V) GOTO 7   
      K=I   
      GOTO 6
   7  CONTINUE  
   6  KK=K-1
      KKK=K+1   
      D1=V-XDWS(KK) 
      D2=V-XDWS(K)  
      D3=V-XDWS(KKK)
      D12=XDWS(KK)-XDWS(K)  
      D13=XDWS(KK)-XDWS(KKK)
      D23=XDWS(K)-XDWS(KKK) 
      D=YDWS(KK)*D2*D3/(D12*D13)-YDWS(K)*D1*D3/(D12*D23)+YDWS(KKK)*&
        D1*D2/(D13*D23)
      GOTO 8
   5  Y=.5/V
      D=Y*(1.+Y/V)  
   8  VOIGT=5.641895836d-1*D
   9  IF(VV.LT.0.) VOIGT=-VOIGT 
      RETURN
   1  Z=CMPLX(A,-V) 
      Z=((((((A6*Z+A5)*Z+A4)*Z+A3)*Z+A2)*Z+A1)*Z+A0)/&
      (((((((Z+B6)*Z+B5)*Z+B4)*Z+B3)*Z+B2)*Z+B1)*Z+B0)
      IF(J.NE.0) GOTO 2 
      VOIGT=REAL(Z) 
      RETURN
   2  VOIGT=.5d0*AIMAG(Z) 
      GOTO 9
      END function voigt				
				
! ---------------------------------------------------------
! Given x(:) and y(:) which tabulate a function and the derivative at the boundary points
! this function returns the second derivative of the spline at each point
! ---------------------------------------------------------
		subroutine splin1(x,y,yp1,ypn,y2)
		real*8, INTENT(IN) :: x(:), y(:), yp1, ypn
		real*8, INTENT(INOUT) :: y2(size(x))
		integer :: n, i, k
		real*8 :: p, qn, sig, un, u(size(x))

			n = size(x)
			
			if (yp1 > .99d30) then
				y2(1) = 0.d0
				u(1) = 0.d0
			else
				y2(1) = -0.5d0
				u(1) = (3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
			endif

			do i = 2, n-1
				sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))				
				p = sig * y2(i-1)+2.d0
				y2(i) = (sig-1.d0)/p
				u(i) = (6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/&
					(x(i+1)-x(i-1))-sig*u(i-1))/p
			enddo
			if (ypn > .99d30) then
				qn = 0.d0
				un = 0.d0
			else
				qn = 0.5d0
				un = (3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
			endif
			
			y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.d0)

			do k = n-1, 1, -1
				y2(k) = y2(k)*y2(k+1)+u(k)
			enddo

		end subroutine splin1

! ---------------------------------------------------------
! Given xa(:) and ya(:) which tabulate a function, returns the interpolation using
! splines of vector x(:) in y(:)
! ---------------------------------------------------------
		subroutine spline(xa,ya,x,y)
		real*8, INTENT(INOUT) :: y(:)
		real*8, INTENT(IN) :: xa(:), ya(:), x(:)
		real*8 :: y2a(size(xa))
		integer :: n_x, n, i, k, khi, klo
		real*8 :: a, b, h, extrap
			
			n = size(xa)
			n_x = size(x)
			call splin1(xa,ya,1.d30,1.d30,y2a)

			do i = 1, n_x					

! Downward extrapolation 
				if (x(i) < xa(1)) then
!					y(i) = ya(1)
					y(i) = ya(1) + (ya(1)-ya(2))/(xa(1)-xa(2)) * (xa(1) - x(i))
				else 

! Upward extrapolation
				if (x(i) > xa(n)) then
!					y(i) = ya(n)
					y(i) = ya(n) + (ya(n)-ya(n-1)) / (xa(n)-xa(n-1)) * (x(i) - xa(n))
				else
! In range
						klo = 1
						khi = n
1						if(khi-klo > 1) then
							k = (khi+klo)/2
							if (xa(k) > x(i)) then
								khi = k
							else
								klo = k
							endif					
							go to 1
						endif

						h = xa(khi)-xa(klo)

						if (h == 0.d0) then
							print *, 'bad xa input in spline'
							stop	
						endif
						a = (xa(khi)-x(i))/h
						b = (x(i)-xa(klo))/h

						y(i) = a*ya(klo)+b*ya(khi)+((a**3.d0-a)*y2a(klo)+(b**3.d0-b)*y2a(khi))*(h**2.d0)/6.d0		
					endif
				endif
			enddo

		end subroutine spline

! ---------------------------------------------------------
!	Given an array xx(1:n), and given a value x, returns a value jlo such that x is between
!	xx(jlo) and xx(jlo+1). xx(1:n) must be monotonic, either increasing or decreasing.
!	jlo=0 or jlo=n is returned to indicate that x is out of range. jlo on input is taken as
!	the initial guess for jlo on output.
! ---------------------------------------------------------
	subroutine hunt(xx,n,x,jlo)
	integer :: jlo,n
	real(kind=8) :: x,xx(n)
	integer :: inc,jhi,jm
	logical :: ascnd

		ascnd=xx(n).ge.xx(1)
		if (jlo.le.0.or.jlo.gt.n) then
			jlo=0
			jhi=n+1
			goto 3
		endif
		inc=1
		if (x.ge.xx(jlo).eqv.ascnd) then
1     	jhi=jlo+inc
			if (jhi.gt.n) then
				jhi=n+1
			else if (x.ge.xx(jhi).eqv.ascnd) then
				jlo=jhi
				inc=inc+inc
				goto 1
			endif
		else
			jhi=jlo
2     	jlo=jhi-inc
			if (jlo.lt.1) then
				jlo=0
			else if (x.lt.xx(jlo).eqv.ascnd) then
				jhi=jlo
				inc=inc+inc
				goto 2
			endif
		endif
3 		if (jhi-jlo.eq.1) then
			if(x.eq.xx(n)) jlo=n-1
			if(x.eq.xx(1)) jlo=1
			return
		endif
		jm = (jhi+jlo)/2
		if (x.ge.xx(jm).eqv.ascnd) then
			jlo=jm
		else
			jhi=jm
		endif
		goto 3
	end subroutine hunt

! ---------------------------------------------------------
! Given xa(:) and ya(:) which tabulate a function, returns the interpolation using
! linear interpolation of vector x(:) in y(:)
! ---------------------------------------------------------
	subroutine lin_interpol(xa,ya,x,y)
   real(kind=8), INTENT(IN) :: xa(:), ya(:), x(:)
   real(kind=8), INTENT(INOUT) :: y(:)
   integer :: i, n, na
   integer :: loc, jlo

   	n = size(x)
   	na = size(xa)

      do i = 1, n
      	call hunt(xa,na,x(i),jlo)
         if (jlo == 0) then
         	y(i) = ya(1)
         else if (jlo == na) then
         	y(i) = ya(na)
         else
            y(i) = (ya(jlo+1)-ya(jlo))/(xa(jlo+1)-xa(jlo)) * (x(i)-xa(jlo)) + ya(jlo)
         endif
      enddo

   end subroutine lin_interpol
		
! ---------------------------------------------------------
! Given xa(:) and ya(:) which tabulate a function z(nxmny), returns the interpolation using
! spline interpolation in 2D
! ---------------------------------------------------------		
		function interpol_2d(xa,ya,z,x0,y0)
		real(kind=8) :: interpol_2d
		real(kind=8) :: xa(:), ya(:), z(:,:), x0, y0, res(1), xin(1)
		real(kind=8), allocatable :: temp(:)
		integer :: nx, ny, j
			nx = size(xa)
			ny = size(ya)
			allocate(temp(nx))
			do j = 1, nx
				xin(1) = y0
				call spline(ya,z(j,:),xin,res)
				temp(j) = res(1)
			enddo
			
			xin(1) = x0
			call spline(xa,temp,xin,res)
			
			interpol_2d = res(1)
			
		end function interpol_2d
		
!--------------------------------------------------------------
! Inversion of a 4x4 matrix
!--------------------------------------------------------------
	subroutine invert(a)
	real*8 :: a(4,4)
	real*8 :: b(4,4), det, maxim, fabsmax
	! First some tests of singularity
		b = dabs(a)
		maxim = maxval(b)
		fabsmax = 1.d0 / maxim
		if (maxim == 0.d0) then
			print *, 'Singularity in the inversion'
!			stop
		endif

		a = a * fabsmax

   	b(1,1) = a(2,2) * a(3,3) * a(4,4) + a(2,3) * a(3,4) * a(4,2)&
      	+ a(2,4) * a(3,2) * a(4,3) - a(2,2) * a(3,4) * a(4,3)&
	   	- a(2,3) * a(3,2) * a(4,4) - a(2,4) * a(3,3) * a(4,2)
   	b(2,1) = a(2,3) * a(3,1) * a(4,4) + a(2,4) * a(3,3) * a(4,1)&
	   	+ a(2,1) * a(3,4) * a(4,3) - a(2,3) * a(3,4) * a(4,1)&
   		- a(2,4) * a(3,1) * a(4,3) - a(2,1) * a(3,3) * a(4,4)
   	b(3,1) = a(2,4) * a(3,1) * a(4,2) + a(2,1) * a(3,2) * a(4,4)&
	   	+ a(2,2) * a(3,4) * a(4,1) - a(2,4) * a(3,2) * a(4,1)&
   		- a(2,1) * a(3,4) * a(4,2) - a(2,2) * a(3,1) * a(4,4)
   	b(4,1) = a(2,1) * a(3,3) * a(4,2) + a(2,2) * a(3,1) * a(4,3)&
	   	+ a(2,3) * a(3,2) * a(4,1) - a(2,1) * a(3,2) * a(4,3)&
   		- a(2,2) * a(3,3) * a(4,1) - a(2,3) * a(3,1) * a(4,2)
   	b(1,2) = a(3,2) * a(4,4) * a(1,3) + a(3,3) * a(4,2) * a(1,4)&
	   	+ a(3,4) * a(4,3) * a(1,2) - a(3,2) * a(4,3) * a(1,4)&
   		- a(3,3) * a(4,4) * a(1,2) - a(3,4) * a(4,2) * a(1,3)
   	b(2,2) = a(3,3) * a(4,4) * a(1,1) + a(3,4) * a(4,1) * a(1,3)&
	   	+ a(3,1) * a(4,3) * a(1,4) - a(3,3) * a(4,1) * a(1,4)&
   		- a(3,4) * a(4,3) * a(1,1) - a(3,1) * a(4,4) * a(1,3)
   	b(3,2) = a(3,4) * a(4,2) * a(1,1) + a(3,1) * a(4,4) * a(1,2)&
	   	+ a(3,2) * a(4,1) * a(1,4) - a(3,4) * a(4,1) * a(1,2)&
   		- a(3,1) * a(4,2) * a(1,4) - a(3,2) * a(4,4) * a(1,1)
   	b(4,2) = a(3,1) * a(4,2) * a(1,3) + a(3,2) * a(4,3) * a(1,1)&
	   	+ a(3,3) * a(4,1) * a(1,2) - a(3,1) * a(4,3) * a(1,2)&
   		- a(3,2) * a(4,1) * a(1,3) - a(3,3) * a(4,2) * a(1,1)
   	b(1,3) = a(4,2) * a(1,3) * a(2,4) + a(4,3) * a(1,4) * a(2,2)&
   		+ a(4,4) * a(1,2) * a(2,3) - a(4,2) * a(1,4) * a(2,3)&
	   	- a(4,3) * a(1,2) * a(2,4) - a(4,4) * a(1,3) * a(2,2)
   	b(2,3) = a(4,3) * a(1,1) * a(2,4) + a(4,4) * a(1,3) * a(2,1)&
	   	+ a(4,1) * a(1,4) * a(2,3) - a(4,3) * a(1,4) * a(2,1)&
   		- a(4,4) * a(1,1) * a(2,3) - a(4,1) * a(1,3) * a(2,4)
   	b(3,3) = a(4,4) * a(1,1) * a(2,2) + a(4,1) * a(1,2) * a(2,4)&
	   	+ a(4,2) * a(1,4) * a(2,1) - a(4,4) * a(1,2) * a(2,1)&
   		- a(4,1) * a(1,4) * a(2,2) - a(4,2) * a(1,1) * a(2,4)
   	b(4,3) = a(4,1) * a(1,3) * a(2,2) + a(4,2) * a(1,1) * a(2,3)&
	   	+ a(4,3) * a(1,2) * a(2,1) - a(4,1) * a(1,2) * a(2,3)&
   		- a(4,2) * a(1,3) * a(2,1) - a(4,3) * a(1,1) * a(2,2)
   	b(1,4) = a(1,2) * a(2,4) * a(3,3) + a(1,3) * a(2,2) * a(3,4)&
	   	+ a(1,4) * a(2,3) * a(3,2) - a(1,2) * a(2,3) * a(3,4)&
   		- a(1,3) * a(2,4) * a(3,2) - a(1,4) * a(2,2) * a(3,3)
   	b(2,4) = a(1,3) * a(2,4) * a(3,1) + a(1,4) * a(2,1) * a(3,3)&
	   	+ a(1,1) * a(2,3) * a(3,4) - a(1,3) * a(2,1) * a(3,4)&
   		- a(1,4) * a(2,3) * a(3,1) - a(1,1) * a(2,4) * a(3,3)
   	b(3,4) = a(1,4) * a(2,2) * a(3,1) + a(1,1) * a(2,4) * a(3,2)&
	   	+ a(1,2) * a(2,1) * a(3,4) - a(1,4) * a(2,1) * a(3,2)&
   		- a(1,1) * a(2,2) * a(3,4) - a(1,2) * a(2,4) * a(3,1)
   	b(4,4) = a(1,1) * a(2,2) * a(3,3) + a(1,2) * a(2,3) * a(3,1)&
	   	+ a(1,3) * a(2,1) * a(3,2) - a(1,1) * a(2,3) * a(3,2)&
   		- a(1,2) * a(2,1) * a(3,3) - a(1,3) * a(2,2) * a(3,1)

		det = a(1,1) * b(1,1) + a(1,2) * b(2,1) + a(1,3) * b(3,1) + a(1,4) * b(4,1)			

		a = b * (fabsmax / det)

	end subroutine invert		
		
!----------------------------------------------------------------
! Given the values of the Zeeman opacities, fill the absorption matrix
! for each depth point
!----------------------------------------------------------------		
		subroutine fill_matrix(opacity,ab_matrix)
		real(kind=8) :: opacity(:,:), ab_matrix(:,:,:)
		
			ab_matrix(1,1,:) = opacity(1,:)   !eta_I
			ab_matrix(2,2,:) = opacity(1,:)   !eta_I
			ab_matrix(3,3,:) = opacity(1,:)   !eta_I
			ab_matrix(4,4,:) = opacity(1,:)   !eta_I
			ab_matrix(1,2,:) = opacity(2,:)   !eta_Q
			ab_matrix(2,1,:) = opacity(2,:)   !eta_Q			
			ab_matrix(1,3,:) = opacity(3,:)   !eta_U
			ab_matrix(3,1,:) = opacity(3,:)   !eta_U			
			ab_matrix(1,4,:) = opacity(4,:)   !eta_V
			ab_matrix(4,1,:) = opacity(4,:)   !eta_V			
			ab_matrix(2,3,:) = opacity(7,:)   !rho_V
			ab_matrix(3,2,:) = -opacity(7,:)  !-rho_V
			ab_matrix(2,4,:) = -opacity(6,:)  !-rho_U
			ab_matrix(4,2,:) = opacity(6,:)   !rho_U						
			ab_matrix(3,4,:) = opacity(5,:)   !rho_Q			
			ab_matrix(4,3,:) = -opacity(5,:)  !-rho_Q
			
!			ab_matrix(1,1,:) = opacity(1,:)   !eta_I
!			ab_matrix(2,2,:) = opacity(1,:)   !eta_I
!			ab_matrix(3,3,:) = opacity(1,:)   !eta_I
!			ab_matrix(4,4,:) = opacity(1,:)   !eta_I
!			ab_matrix(2,1,:) = opacity(2,:)   !eta_Q
!			ab_matrix(1,2,:) = opacity(2,:)   !eta_Q			
!			ab_matrix(3,1,:) = opacity(3,:)   !eta_U
!			ab_matrix(1,3,:) = opacity(3,:)   !eta_U			
!			ab_matrix(4,1,:) = opacity(4,:)   !eta_V
!			ab_matrix(1,4,:) = opacity(4,:)   !eta_V			
!			ab_matrix(3,2,:) = opacity(7,:)   !rho_V
!			ab_matrix(2,3,:) = -opacity(7,:)  !-rho_V
!			ab_matrix(4,2,:) = -opacity(6,:)  !-rho_U
!			ab_matrix(2,4,:) = opacity(6,:)   !rho_U						
!			ab_matrix(4,3,:) = opacity(5,:)   !rho_Q			
!			ab_matrix(3,4,:) = -opacity(5,:)  !-rho_Q
			
		end subroutine fill_matrix
		
!----------------------------------------------------------------
! Given the values of the Zeeman opacities, fill the absorption matrix
! for each depth point
!----------------------------------------------------------------		
		subroutine fill_matrix_one(opacity,ab_matrix)
		real(kind=8) :: opacity(:), ab_matrix(:,:)
		
			ab_matrix(1,1) = opacity(1)	!eta_I
			ab_matrix(2,2) = opacity(1)	!eta_I
			ab_matrix(3,3) = opacity(1)	!eta_I
			ab_matrix(4,4) = opacity(1)	!eta_I
			ab_matrix(1,2) = opacity(2)	!eta_Q
			ab_matrix(2,1) = opacity(2)	!eta_Q		  
			ab_matrix(1,3) = opacity(3)	!eta_U
			ab_matrix(3,1) = opacity(3)	!eta_U		  
			ab_matrix(1,4) = opacity(4)	!eta_V
			ab_matrix(4,1) = opacity(4)	!eta_V		  
			ab_matrix(2,3) = opacity(7)	!rho_V
			ab_matrix(3,2) = -opacity(7)  !-rho_V
			ab_matrix(2,4) = -opacity(6)  !-rho_U
			ab_matrix(4,2) = opacity(6)	!rho_U					  
			ab_matrix(3,4) = opacity(5)	!rho_Q		  
			ab_matrix(4,3) = -opacity(5)  !-rho_Q
			
!			ab_matrix(1,1,:) = opacity(1,:)   !eta_I
!			ab_matrix(2,2,:) = opacity(1,:)   !eta_I
!			ab_matrix(3,3,:) = opacity(1,:)   !eta_I
!			ab_matrix(4,4,:) = opacity(1,:)   !eta_I
!			ab_matrix(2,1,:) = opacity(2,:)   !eta_Q
!			ab_matrix(1,2,:) = opacity(2,:)   !eta_Q			
!			ab_matrix(3,1,:) = opacity(3,:)   !eta_U
!			ab_matrix(1,3,:) = opacity(3,:)   !eta_U			
!			ab_matrix(4,1,:) = opacity(4,:)   !eta_V
!			ab_matrix(1,4,:) = opacity(4,:)   !eta_V			
!			ab_matrix(3,2,:) = opacity(7,:)   !rho_V
!			ab_matrix(2,3,:) = -opacity(7,:)  !-rho_V
!			ab_matrix(4,2,:) = -opacity(6,:)  !-rho_U
!			ab_matrix(2,4,:) = opacity(6,:)   !rho_U						
!			ab_matrix(4,3,:) = opacity(5,:)   !rho_Q			
!			ab_matrix(3,4,:) = -opacity(5,:)  !-rho_Q
			
		end subroutine fill_matrix_one		

!----------------------------------------------------------------
! Given the values of the Zeeman opacities, fill the absorption matrix
! for each depth point and frequency
!----------------------------------------------------------------		
		subroutine fill_matrix_vec(opacity,ab_matrix)
		real(kind=8) :: opacity(:,:,:), ab_matrix(:,:,:,:)
		
			ab_matrix(1,1,:,:) = opacity(1,:,:)   !eta_I
			ab_matrix(2,2,:,:) = opacity(1,:,:)   !eta_I
			ab_matrix(3,3,:,:) = opacity(1,:,:)   !eta_I
			ab_matrix(4,4,:,:) = opacity(1,:,:)   !eta_I
			ab_matrix(1,2,:,:) = opacity(2,:,:)   !eta_Q
			ab_matrix(2,1,:,:) = opacity(2,:,:)   !eta_Q
			ab_matrix(1,3,:,:) = opacity(3,:,:)   !eta_U
			ab_matrix(3,1,:,:) = opacity(3,:,:)   !eta_U
			ab_matrix(1,4,:,:) = opacity(4,:,:)   !eta_V
			ab_matrix(4,1,:,:) = opacity(4,:,:)   !eta_V
			ab_matrix(2,3,:,:) = opacity(7,:,:)   !rho_V
			ab_matrix(3,2,:,:) = -opacity(7,:,:)  !-rho_V
			ab_matrix(2,4,:,:) = -opacity(6,:,:)  !-rho_U
			ab_matrix(4,2,:,:) = opacity(6,:,:)   !rho_U						
			ab_matrix(3,4,:,:) = opacity(5,:,:)   !rho_Q			
			ab_matrix(4,3,:,:) = -opacity(5,:,:)  !-rho_Q
			
		end subroutine fill_matrix_vec

!----------------------------------------------------------------
! This function returns the strength of the Zeeman components
!----------------------------------------------------------------				
	function strength_zeeman(J_up,J_low,M_up,M_low)
	real(kind=8) :: J_up, J_low, M_up, M_low, strength_zeeman

		strength_zeeman = 3.d0 * w3js(int(2.d0*J_up),int(2.d0*J_low),2,-int(2.d0*M_up),&
					int(2.d0*M_low),-int(2.d0*(M_low-M_up)))**2
	end function strength_zeeman

!----------------------------------------------------------------
! This function initializes some things
!----------------------------------------------------------------		
		subroutine init_maths
		integer :: i
			call factrl
			identity = 0.d0
			do i = 1, 4
				identity(i,i) = 1.d0
			enddo
      end subroutine init_maths
      
	subroutine qsortd(x,ind,n)
	
	implicit none
	integer, parameter  :: dp = 8
	
	real (kind=8), intent(in)  :: x(:)
	integer, intent(out)   :: ind(:)
	integer, intent(in)    :: n
	
	!***************************************************************************
	
	!                                                         robert renka
	!                                                 oak ridge natl. lab.
	
	!   this subroutine uses an order n*log(n) quick sort to sort a real (dp)
	! array x into increasing order.  the algorithm is as follows.  ind is
	! initialized to the ordered sequence of indices 1,...,n, and all interchanges
	! are applied to ind.  x is divided into two portions by picking a central
	! element t.  the first and last elements are compared with t, and
	! interchanges are applied as necessary so that the three values are in
	! ascending order.  interchanges are then applied so that all elements
	! greater than t are in the upper portion of the array and all elements
	! less than t are in the lower portion.  the upper and lower indices of one
	! of the portions are saved in local arrays, and the process is repeated
	! iteratively on the other portion.  when a portion is completely sorted,
	! the process begins again by retrieving the indices bounding another
	! unsorted portion.
	
	! input parameters -   n - length of the array x.
	
	!                      x - vector of length n to be sorted.
	
	!                    ind - vector of length >= n.
	
	! n and x are not altered by this routine.
	
	! output parameter - ind - sequence of indices 1,...,n permuted in the same
	!                          fashion as x would be.  thus, the ordering on
	!                          x is defined by y(i) = x(ind(i)).
	
	!*********************************************************************
	
	! note -- iu and il must be dimensioned >= log(n) where log has base 2.
	
	!*********************************************************************
	
	integer   :: iu(21), il(21)
	integer   :: m, i, j, k, l, ij, it, itt, indx
	real(kind=8)      :: r
	real(kind=8) :: t
	
	! local parameters -
	
	! iu,il =  temporary storage for the upper and lower
	!            indices of portions of the array x
	! m =      index for iu and il
	! i,j =    lower and upper indices of a portion of x
	! k,l =    indices in the range i,...,j
	! ij =     randomly chosen index between i and j
	! it,itt = temporary storage for interchanges in ind
	! indx =   temporary index for x
	! r =      pseudo random number for generating ij
	! t =      central element of x
	
	if (n <= 0) return
	
	! initialize ind, m, i, j, and r
	
	do  i = 1, n
	ind(i) = i
	end do
	
	m = 1
	i = 1
	j = n
	r = .375
	
	! top of loop
	
	20 if (i >= j) go to 70
	if (r <= .5898437) then
	r = r + .0390625
	else
	r = r - .21875
	end if
	
	! initialize k
	
	30 k = i
	
	! select a central element of x and save it in t
	
	ij = i + r*(j-i)
	it = ind(ij)
	t = x(it)
	
	! if the first element of the array is greater than t,
	!   interchange it with t
	
	indx = ind(i)
	if (x(indx) > t) then
	ind(ij) = indx
	ind(i) = it
	it = indx
	t = x(it)
	end if
	
	! initialize l
	
	l = j
	
	! if the last element of the array is less than t,
	!   interchange it with t
	
	indx = ind(j)
	if (x(indx) >= t) go to 50
	ind(ij) = indx
	ind(j) = it
	it = indx
	t = x(it)
	
	! if the first element of the array is greater than t,
	!   interchange it with t
	
	indx = ind(i)
	if (x(indx) <= t) go to 50
	ind(ij) = indx
	ind(i) = it
	it = indx
	t = x(it)
	go to 50
	
	! interchange elements k and l
	
	40 itt = ind(l)
	ind(l) = ind(k)
	ind(k) = itt
	
	! find an element in the upper part of the array which is
	!   not larger than t
	
	50 l = l - 1
	indx = ind(l)
	if (x(indx) > t) go to 50
	
	! find an element in the lower part of the array whcih is not smaller than t
	
	60 k = k + 1
	indx = ind(k)
	if (x(indx) < t) go to 60
	
	! if k <= l, interchange elements k and l
	
	if (k <= l) go to 40
	
	! save the upper and lower subscripts of the portion of the
	!   array yet to be sorted
	
	if (l-i > j-k) then
	il(m) = i
	iu(m) = l
	i = k
	m = m + 1
	go to 80
	end if
	
	il(m) = k
	iu(m) = j
	j = l
	m = m + 1
	go to 80
	
	! begin again on another unsorted portion of the array
	
	70 m = m - 1
	if (m == 0) return
	i = il(m)
	j = iu(m)
	
	80 if (j-i >= 11) go to 30
	if (i == 1) go to 20
	i = i - 1
	
	! sort elements i+1,...,j.  note that 1 <= i < j and j-i < 11.
	
	90 i = i + 1
	if (i == j) go to 70
	indx = ind(i+1)
	t = x(indx)
	it = indx
	indx = ind(i)
	if (x(indx) <= t) go to 90
	k = i
	
	100 ind(k+1) = ind(k)
	k = k - 1
	indx = ind(k)
	if (t < x(indx)) go to 100
	
	ind(k+1) = it
	go to 90
	end subroutine qsortd
	
	
	
	subroutine cubint ( ftab, xtab, ntab, ia, ib, result, error )
!
!***********************************************************************
!
!! CUBINT approximates an integral using cubic interpolation of data.
!
!
!  Discussion:
!
!    The integral to be approximated is
! 
!      INTEGRAL (XTAB(IB) to XTAB(IA)) F(X) DX
!
!    The routine estimates the error in integration.
!
!  Reference:
!
!    Philip Davis and Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Blaisdell Publishing, 1967.
!
!    P E Gill and G F Miller
!    An algorithm for the integration of unequally spaced data,
!    Comput J, Number 15, 1972, pages 80-83.
!
!  Modified:
!
!    30 October 2000
!
!  Parameters:
!
!    Input, real FTAB(NTAB), contains the tabulated function
!    values, FTAB(I) = F(XTAB(I)).
!
!    Input, real XTAB(NTAB), contains the points at which the
!    function was tabulated.  XTAB should contain distinct
!    values, given in ascending order.
!
!    Input, integer NTAB, the number of tabulated points.
!    NTAB must be at least 4.
!
!    Input, integer IA, the entry of XTAB at which integration
!    is to begin.  IA must be no less than 1 and no greater
!    than NTAB.
!
!    Input, integer IB, the entry of XTAB at which integration
!    is to end.  IB must be no less than 1 and no greater than
!    NTAB.
!
!    Output, real RESULT, the approximate value of the
!    integral from XTAB(IA) to XTAB(IB) of the function.
!
!    Output, real ERROR, an estimate of the error in
!    integration.
!
  implicit none
!
  integer ntab
!
  real(kind=8) c
  real(kind=8) d1
  real(kind=8) d2
  real(kind=8) d3
  real(kind=8) error
  real(kind=8) ftab(ntab)
  real(kind=8) h1
  real(kind=8) h2
  real(kind=8) h3
  real(kind=8) h4
  integer i
  integer ia
  integer ib
  integer ind
  integer it
  integer j
  integer k
  real(kind=8) r1
  real(kind=8) r2
  real(kind=8) r3
  real(kind=8) r4
  real(kind=8) result
  real(kind=8) s
  real(kind=8) term
  real(kind=8) xtab(ntab)
!
  result = 0.0E+00
  error = 0.0E+00
 
  if ( ia == ib ) then
    return
  end if
 
  if ( ntab < 4 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBINT - Fatal error!'
    write ( *, '(a,i6)' ) '  NTAB must be at least 4, but input NTAB = ',ntab
    stop
  end if
 
  if ( ia < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBINT - Fatal error!'
    write ( *, '(a,i6)' ) '  IA must be at least 1, but input IA = ',ia
    stop
  end if
 
  if ( ia > ntab ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBINT - Fatal error!'
    write ( *, '(a,i6)' ) '  IA must be <= NTAB, but input IA=',ia
    stop
  end if
 
  if ( ib < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBINT - Fatal error!'
    write ( *, '(a,i6)' ) '  IB must be at least 1, but input IB = ',ib
    stop
  end if
 
  if ( ib > ntab ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBINT - Fatal error!'
    write ( *, '(a,i6)' ) '  IB must be <= NTAB, but input IB=',ib
    stop
  end if
!
!  Temporarily switch IA and IB, and store minus sign in IND
!  so that, while integration is carried out from low X's
!  to high ones, the sense of the integral is preserved.
!
  if ( ia > ib ) then
    ind = -1
    it = ib
    ib = ia
    ia = it
  else
    ind = 1
  end if
 
  s = 0.0E+00
  c = 0.0E+00
  r4 = 0.0E+00
  j = ntab-2
  if ( ia < ntab-1 .or. ntab == 4 ) then
    j=max(3,ia)
  end if

  k = 4
  if ( ib > 2 .or. ntab == 4 ) then
    k=min(ntab,ib+2)-1
  end if
 
  do i = j, k
 
    if ( i <= j ) then
 
      h2 = xtab(j-1)-xtab(j-2)
      d3 = (ftab(j-1)-ftab(j-2)) / h2
      h3 = xtab(j)-xtab(j-1)
      d1 = (ftab(j)-ftab(j-1)) / h3
      h1 = h2+h3
      d2 = (d1-d3)/h1
      h4 = xtab(j+1)-xtab(j)
      r1 = (ftab(j+1)-ftab(j)) / h4
      r2 = (r1-d1) / (h4+h3)
      h1 = h1+h4
      r3 = (r2-d2) / h1
 
      if ( ia <= 1 ) then
        result = h2 * (ftab(1)+h2*(0.5*d3-h2*(d2/6.0-(h2+h3+h3)*r3/12.)))
        s = -h2**3 * (h2*(3.0*h2+5.0*h4)+10.0*h3*h1)/60.0
      end if
 
    else
 
      h4 = xtab(i+1)-xtab(i)
      r1 = (ftab(i+1)-ftab(i))/h4
      r4 = h4+h3
      r2 = (r1-d1)/r4
      r4 = r4+h2
      r3 = (r2-d2)/r4
      r4 = (r3-d3)/(r4+h1)
 
    end if
 
    if ( i > ia .and. i <= ib ) then
 
      term = h3*((ftab(i)+ftab(i-1))*0.5-h3*h3*(d2+r2+(h2-h4)*r3) / 12.0 )
      result = result+term
      c = h3**3*(2.0E+00 *h3*h3+5.*(h3*(h4+h2) + 2.0 * h2 * h4 ) ) / 120.0E+00
      error = error+(c+s)*r4
 
      if ( i /= j ) then
        s = c
      else
        s = s+c+c
      end if
 
    else
 
      error = error+r4*s
 
    end if
 
    if ( i >= k ) then
 
      if ( ib >= ntab ) then
        term = h4*(ftab(ntab) - h4*(0.5*r1+h4*(r2/6.0 +(h3+h3+h4)*r3/12.)))
        result = result + term
        error = error - h4**3 * r4 * &
          ( h4 * ( 3.0 * h4 + 5.0 * h2 ) &
          + 10.0 * h3 * ( h2 + h3 + h4 ) ) / 60.0E+00
      end if
 
      if ( ib >= ntab-1 ) error=error+s*r4
    else
      h1 = h2
      h2 = h3
      h3 = h4
      d1 = r1
      d2 = r2
      d3 = r3
    end if
 
  end do
!
!  Restore original values of IA and IB, reverse signs
!  of RESULT and ERROR, to account for integration
!  that proceeded from high X to low X.
!
  if ( ind /= 1 ) then
    it = ib
    ib = ia
    ia = it
    result = -result
    error = -error
  end if
 
  return
end subroutine 

!----------------------------------------------------------------
! This function integrates a tabulated function
!----------------------------------------------------------------		
		function int_tabulated(x, f)
		real(kind=8) :: x(:), f(:), int_tabulated, res, error_res
		integer :: n
			n = size(x)
			call cubint (f, x, n, 1, n, res, error_res)
			int_tabulated = res
		end function int_tabulated

! ---------------------------------------------------------
! Return the input array but shifted by a sub-pixel amount
! ---------------------------------------------------------
		function fft_shift(x, sh)
		real(kind=8) :: x(:), fft_shift(size(x)), sh
		complex(kind=8) :: k(size(x)), ff(size(x))
		integer :: nx, i, n21
		
			nx = size(x)

			do i = 1, nx
         	k(i) = i-1
         enddo
         n21 = nx/2+1
         do i = n21+1, nx
         	k(i) = n21-nx+(i-1-n21)
         enddo

         k = 2.d0*PI*k/nx
         k = k * cmplx(0.d0,1.d0)

         ff = x

         ff = fft_forward(ff, nx)

         fft_shift = fft_backward(ff * exp(k*sh), nx)
		
		end function fft_shift


end module math_functions
