module chain_analysis_module
use vars
use math_functions
use fparser, only : initf, parsef, evalf
implicit none

contains

!------------------------------------------------------------
! Read the chains
!------------------------------------------------------------
	subroutine read_chains(fich_chain, chain_analysis)
	type(chain_analysis_type) :: chain_analysis
	character(len=120) :: fich_chain
	integer :: nparams, chain_length, i
	
		open(unit=12,file=trim(adjustl(fich_chain))//'.header',action='read',status='old')
		read(12,*) chain_analysis%nparams, chain_analysis%nlength
		close(12)
		
		allocate(chain_analysis%chain(chain_analysis%nparams+3,chain_analysis%nlength))
		open(unit=12,file=trim(adjustl(fich_chain))//'.data',action='read',&
			status='old',form='unformatted')
		do i = 1, chain_analysis%nlength
			read(12) chain_analysis%chain(:,i)
		enddo
		close(12)
				
	end subroutine read_chains

!------------------------------------------------------------
! Read the chains
!------------------------------------------------------------
	subroutine read_chains_multinest(fich_chain, chain_analysis)
	type(chain_analysis_type) :: chain_analysis
	character(len=120) :: fich_chain
	integer :: nparams, chain_length, i, j

		i = 0
		open(unit=12,file=trim(adjustl(fich_chain))//'post_equal_weights.dat',action='read',status='old')
		do while (.true.)
			read(12,*,end=999)
			i = i + 1
		enddo
999   continue
		close(12)
		chain_analysis%nparams = inversion%nparams_invert
		
		chain_analysis%nlength = i
		print *, 'Length of posterior samples : ', chain_analysis%nlength
		print *, 'Number of parameters : ', chain_analysis%nparams
		
		allocate(chain_analysis%chain(chain_analysis%nparams+1,chain_analysis%nlength))
		open(unit=12,file=trim(adjustl(fich_chain))//'post_equal_weights.dat',action='read',&
			status='old')
		do i = 1, chain_analysis%nlength
			read(12,*) (chain_analysis%chain(j,i),j=1,chain_analysis%nparams+1)
		enddo
		close(12)

! Read evidence
		open(unit=12,file=trim(adjustl(fich_chain))//'stats.dat',action='read',status='old')
		read(12,FMT='(16X,E20.12)') chain_analysis%evidence
		close(12)
		
	end subroutine read_chains_multinest


!------------------------------------------------------------
! Remove burn-in
!------------------------------------------------------------
	subroutine remove_burnin(chain_analysis, percent)
	type(chain_analysis_type) :: chain_analysis
	real(kind=8), pointer :: chain_new(:,:)
	real(kind=8) :: percent
	integer :: burn_in
		
		burn_in = chain_analysis%nlength * percent / 100.d0		
		if (burn_in > chain_analysis%nlength) burn_in = chain_analysis%nlength
				
		write(*,*) 'Burn-in of ', percent, '%'
		write(*,*) 'Leaving a chain of length : ', chain_analysis%nlength-burn_in
		
		chain_analysis%nlength = chain_analysis%nlength-burn_in
		
		allocate(chain_new(chain_analysis%nparams+2,chain_analysis%nlength))
		chain_new = chain_analysis%chain(:,burn_in:)
		
		deallocate(chain_analysis%chain)
		allocate(chain_analysis%chain(chain_analysis%nparams+2,chain_analysis%nlength))
		chain_analysis%chain = chain_new
				
		deallocate(chain_new)
			
	end subroutine remove_burnin

!------------------------------------------------------------
! Get the optimal binning for building a histogram of a variable
! Ref: Freedman, D, & Diaconis, P. (1981). "On the histogram as a density
!      estimator: L2 theory". Zeitschrift für Wahrscheinlichkeitstheorie und
!      verwandte Gebiete 57 (4): 453–476
!  BIN=2*IQR(var)/n^(1/3)   with n the number of elements of the array 'var'
!------------------------------------------------------------
	function optbin(var)
	real(kind=8) :: optbin, var(:)	
	integer :: n, quart1, quart3
	integer, allocatable :: ind(:)
	
		n = size(var)
		allocate(ind(n))
		
		call qsortd(var,ind,n)
		
		quart1 = n*0.25
		quart3 = n*0.75
		
		optbin = 2.d0*(var(ind(quart3)) - var(ind(quart1))) / (n*1.d0)**(1.d0/3.d0)
		
	end function optbin

!------------------------------------------------------------
! Estimate 0.2%, 2.5%, 16%, 50%, 84%, 97.5%, 99.8% percentiles
! of a given chain
!------------------------------------------------------------
	function percentile_summary(var)
	real(kind=8) :: percentile_summary(7), var(:)
	integer :: n, loc
	integer, allocatable :: ind(:)
	
		n = size(var)
		allocate(ind(n))
		
		call qsortd(var,ind,n)

		loc = n * 0.2 / 100.0
		percentile_summary(1) = var(ind(loc))
		loc = n * 2.5 / 100.0
		percentile_summary(2) = var(ind(loc))
		loc = n * 16.0 / 100.0
		percentile_summary(3) = var(ind(loc))
		loc = n * 50.0 / 100.0
		percentile_summary(4) = var(ind(loc))
		loc = n * 84.0 / 100.0
		percentile_summary(5) = var(ind(loc))
		loc = n * 97.5 / 100.0
		percentile_summary(6) = var(ind(loc))
		loc = n * 99.8 / 100.0
		percentile_summary(7) = var(ind(loc))

		deallocate(ind)
		
	end function percentile_summary

!------------------------------------------------------------
! 1D histograms
!------------------------------------------------------------
	function oned_histogram_size(chain)
	real(kind=8) :: chain(:)
	real(kind=8) :: step, xmin, xmax
	integer :: oned_histogram_size

! Bins
		step = optbin(chain)

		xmin = minval(chain)
		xmax = maxval(chain)
		oned_histogram_size = (xmax-xmin) / step

		return

	end function oned_histogram_size

!------------------------------------------------------------
! 1D histograms
!------------------------------------------------------------
	subroutine oned_histogram(chain, x, yGauss, yStep)
	real(kind=8) :: chain(:)	
	real(kind=8) :: step, xmin, xmax
	real(kind=8) :: x(:), yGauss(:), yStep(:)
	real(kind=8), allocatable :: xx(:)
	integer :: i, n, nn, ncut
	real(kind=8) :: wei, norm, error_norm, sig
	
! Bins
		step = optbin(chain)
! 		write(*,*) '  Histogram step = ', step
		
		xmin = minval(chain)
		xmax = maxval(chain)
		n = (xmax-xmin) / step
		step = (xmax-xmin) / float(n)
		
		do i = 1, n
			x(i) = xmin + step * (i-1.d0)
		enddo
		
! Variables for the normalization
		if (n <= 50) then
			nn = 50
			allocate(xx(nn))
			do i = 1, nn
				xx(i) = xmin + (xmax-xmin)/nn * i
			enddo
		else
			allocate(xx(n))
			xx = x
		endif
				
! Doing the histogram
		sig = 1.2d0
		do i = 1, n
			ncut = count(chain >= x(i)-step/2.d0 .and. chain < x(i)+step/2.0)
			
! Smoothing kernel
			wei = sum(exp(-0.5d0 * (chain-x(i))**2 / (sig*step)**2.d0))
			norm = int_tabulated(xx, exp(-0.5d0 * (xx-x(i))**2 / (sig*step)**2.d0))			
			yGauss(i) = wei / norm
			yStep(i) = ncut
		enddo
		
		deallocate(xx)
		
		yGauss = yGauss / maxval(yGauss)
		yStep = yStep / maxval(yStep)
		
	end subroutine 

!------------------------------------------------------------
! Confidence intervals
! INPUT
!   xdata : x position of the histogram
!   ydata : y position of the histogram
! OUTPUT
!   est : estimated value of the parameter
! errup : value of x at the upper confidence interval
! errdown : value of x at the lower confidence interval
! conf : desired confidence level
!------------------------------------------------------------
	subroutine conf_limits_ccf(xdata, ydata, conf, est, errup, errdown)
	real(kind=8) :: xdata(:), ydata(:)
	real(kind=8) :: conf, est, errup, errdown
	real(kind=8) :: xpdf_low, xpdf_up, xpdf_mid, xmin, xmax, lower, upper, norm
	real(kind=8), allocatable :: xx(:), yy(:), x(:), y(:), xF(:), F(:)
	integer, allocatable :: ind(:)
	integer :: n, npto, loc(1), i
	
		n = size(xdata)
		
		xpdf_low = 0.5d0 - conf/200.
		xpdf_mid = 0.5
		xpdf_up  = 0.5 + conf/200.
		
		if (xpdf_low < 0.0 .or. xpdf_up > 1) then 
			print *,'wrong value for CONF'
			stop
		endif
				
		allocate(xx(n))
		allocate(yy(n))
		allocate(ind(n))
		
		xx = xdata
		yy = ydata / maxval(ydata)
						
! Sorting
		call qsortd(xx,ind,n)
		xx = xx(ind)
		yy = yy(ind)
		deallocate(ind)
						
! Interpolation by splines
		npto = 300
		xmin = xx(1)
		xmax = xx(n)
		
		if (npto > n) then
			allocate(x(npto))
			allocate(y(npto))
			
			do i = 1, npto
				x(i) = xmin + (xmax-xmin)/npto * i
			enddo
			call spline(xx,yy,x,y)
			n = npto
		else
			allocate(x(n))
			allocate(y(n))
			x = xx
			y = yy
		endif
		
! Peak normalization
		y = y / maxval(y)
		
! Computing the cumulative distribution function
		norm = int_tabulated(x,y)
		
		allocate(xF(n-1))
		allocate(F(n-1))
		xF = x(1:n-1)
				
		do i = 5, n
			F(i-1) = int_tabulated(x(1:i-1),y(1:i-1)) / norm
		enddo
				
		loc = minloc(abs(F - xpdf_low))		
		lower = xF(loc(1))
		
		loc = minloc(abs(F - xpdf_mid))		
		est = xF(loc(1))
		
		loc = minloc(abs(F - xpdf_up))
		upper = xF(loc(1))
		
		errup = upper - est
		errdown = est - lower

		deallocate(x,y)
		deallocate(xF,F)
		deallocate(xx,yy)
	
	end subroutine conf_limits_ccf

!------------------------------------------------------------
! Locate index for a give parameter name
!------------------------------------------------------------
	function locate_parameter_byname(chain_analysis, fich_chain, name)
	integer :: locate_parameter_byname, i
	type(chain_analysis_type) :: chain_analysis
	character(len=120) :: fich_chain, str_parameter
	character(len=10) :: name
	
		open(unit=17,file=trim(adjustl(fich_chain))//'.params_inverted',&
			action='read',status='old')

		do i = 1, chain_analysis%nparams
			read(17,*) str_parameter			
			if (trim(adjustl(str_parameter)) == trim(adjustl(name))) then
				locate_parameter_byname = i
				close(17)
				return
			endif
		enddo

		close(17)
		write(*,*) 'Parameter ', name, ' not present'
		stop		

	end function locate_parameter_byname

!------------------------------------------------------------
! Analyze the chains
!------------------------------------------------------------
	subroutine analyze_chains_mcmc(chain_analysis, fich_chain, most_probable)
	type(chain_analysis_type) :: chain_analysis
	real(kind=8) :: most_probable(:)
	character(len=120) :: fich_chain, str_parameter
	real(kind=8), pointer :: x(:), yGauss(:), yStep(:)
	integer :: i, j, n
	real(kind=8) :: est, errup, errdown, chi2_mean, chi2_max
	integer, pointer :: indx(:)
	real(kind=8), allocatable :: chain_function(:)
	
		call read_chains(fich_chain, chain_analysis)
		call remove_burnin(chain_analysis, 50.d0)
		
		open(unit=12,file=trim(adjustl(fich_chain))//'.hist1D',&
			action='write',status='replace',&
			form='unformatted')
		write(12) chain_analysis%nparams, chain_analysis%n_posterior_combinations
		
		open(unit=13,file=trim(adjustl(fich_chain))//'.confidence',&
			action='write',status='replace')
		
		open(unit=14,file=trim(adjustl(fich_chain))//'.params_inverted',&
			action='read',status='old')

! Do histograms for all inverted variables
		do i = 1, chain_analysis%nparams
			call oned_histogram(chain_analysis%chain(i,:), x, yGauss, yStep)
			call conf_limits_ccf(x, yStep, 95.449972d0, est, errup, errdown)			
			write(13,*) est, errdown, errup
			write(13,*) most_probable(i)
			read(14,*) str_parameter
			write(*,FMT='(A,A,F9.4,A,F9.4,A,F9.4)') trim(adjustl(str_parameter)), &
				' : E(x)=', est, '  +', errdown, '  -', errup
			n = size(x)
			write(12) i, n
			write(12) x, yGauss, yStep
			deallocate(x)
			deallocate(yGauss)
			deallocate(yStep)
		enddo

! And then do histograms for combinations of parameters as indicated in the input file

		if (chain_analysis%n_posterior_combinations > 0) then
! First call the initialization routine for the parser
			call initf (chain_analysis%n_posterior_combinations)

			write(*,*)
			write(*,*) 'Posterior for combinations of parameters'
			write(*,*)
			allocate(chain_function(chain_analysis%nlength))
				
			do i = 1, chain_analysis%n_posterior_combinations
			
				allocate(indx(chain_analysis%posterior_fun(i)%n_variables))

				write(*,FMT='(A,A)') ' Function => ', chain_analysis%posterior_fun(i)%posterior_function
			
! Locate the index of the variables that enter into the function
				do j = 1, chain_analysis%posterior_fun(i)%n_variables
					indx(j) = locate_parameter_byname(chain_analysis, fich_chain, &
						chain_analysis%posterior_fun(i)%variables(j))
					write(*,FMT='(A,I2,2X,A,A,I3)') '     Variable ', j,&
					 	trim(adjustl(chain_analysis%posterior_fun(i)%variables(j))), &
						' - index= ', indx(j)
				enddo

! Parse the function generting bytecode
				call parsef(i, chain_analysis%posterior_fun(i)%posterior_function, &
					chain_analysis%posterior_fun(i)%variables)
	
! Then evaluate the chain for the new parameter
				do j = 1, chain_analysis%nlength				
					chain_function(j) = evalf(i,chain_analysis%chain(indx,j))				
				enddo
	

				deallocate(indx)

! Write the posterior marginal distributions
				call oned_histogram(chain_function, x, yGauss, yStep)
				call conf_limits_ccf(x, yStep, 95.449972d0, est, errup, errdown)
				most_probable(i) = est
				write(13,*) est, errdown, errup			
				write(*,FMT='(A,F9.4,A,F9.4,A,F9.4)') ' E(x)=', est, '  +', errdown, '  -', errup
				n = size(x)
				write(12) i, n
				write(12) x, yGauss, yStep
				deallocate(x)
				deallocate(yGauss)
				deallocate(yStep)
			
			enddo

			deallocate(chain_function)
		endif
		
		chi2_mean = sum(chain_analysis%chain(chain_analysis%nparams+2,:)) / &
			max(1,size(chain_analysis%chain(chain_analysis%nparams+2,:)))
		chi2_max = minval(chain_analysis%chain(chain_analysis%nparams+2,:))
		
		write(*,'(A,F9.3)') 'Bayesian complexity = ', chi2_mean - chi2_max
		write(13,'(A,F9.3)') 'Bayesian complexity = ', chi2_mean - chi2_max
		
		close(12)
		close(13)
		
	end subroutine analyze_chains_mcmc


!------------------------------------------------------------
! Analyze the chains
!------------------------------------------------------------
	subroutine analyze_chains_multinest(chain_analysis, fich_chain, most_probable, chi2_min)
	type(chain_analysis_type) :: chain_analysis
	real(kind=8) :: most_probable(:), chi2_min
	character(len=120) :: fich_chain, str_parameter
	real(kind=8), pointer :: x(:), yGauss(:), yStep(:)
	integer :: i, j, n
	real(kind=8) :: est, errup, errdown, chi2_mean, chi2_max, percentiles(7)
	integer, pointer :: indx(:)
	real(kind=8), allocatable :: chain_function(:)
	
		call read_chains_multinest(fich_chain, chain_analysis)
		
		open(unit=12,file=trim(adjustl(fich_chain))//'.hist1D',&
			action='write',status='replace')
		write(12,*) chain_analysis%nparams, chain_analysis%n_posterior_combinations
		
		open(unit=13,file=trim(adjustl(fich_chain))//'.confidence',&
			action='write',status='replace')
		
		open(unit=14,file=trim(adjustl(fich_chain))//'.params_inverted',&
			action='read',status='old')

		open(unit=15,file=trim(adjustl(fich_chain))//'.percentiles',&
			action='write',status='replace')
		write(15,*) 'Percentiles'
		write(15,*) '     0.2%       2.5%        16%         50%         84%        97.5%      99.8%'

! Do histograms for all inverted variables
		do i = 1, chain_analysis%nparams

			read(14,*) str_parameter
			
! Print 1sigma values
			n = oned_histogram_size(chain_analysis%chain(i,:))
			allocate(x(n))
			allocate(yGauss(n))
			allocate(yStep(n))
			
			call oned_histogram(chain_analysis%chain(i,:), x, yGauss, yStep)
			call conf_limits_ccf(x, yStep, 68.268949d0, est, errup, errdown)
			write(13,*) trim(adjustl(str_parameter))
			write(13,*) est, errdown, errup
						
			write(*,FMT='(A,A,F9.4,A,F9.4,A,F9.4)') trim(adjustl(str_parameter)), &
				' (1s) : E(x)=', est, '  +', errdown, '  -', errup

! Print 2sigma values
			call conf_limits_ccf(x, yStep, 95.449972d0, est, errup, errdown)			
			write(13,*) est, errdown, errup
			write(*,FMT='(A,A,F9.4,A,F9.4,A,F9.4)') trim(adjustl(str_parameter)), &
				' (2s) : E(x)=', est, '  +', errdown, '  -', errup
			
! Print MAP values
			write(13,*) most_probable(i)			
			
			n = size(x)
			write(12,*) i, n
			do j = 1, n
				write(12,*) x(j), yGauss(j), yStep(j)
			enddo
			deallocate(x)
			deallocate(yGauss)
			deallocate(yStep)

! Calculate the percentiles
			percentiles = percentile_summary(chain_analysis%chain(i,:))
			write(15,FMT='(7(F12.4))') percentiles

		enddo

		close(14)

! And then do histograms for combinations of parameters as indicated in the input file

		if (chain_analysis%n_posterior_combinations > 0) then
! First call the initialization routine for the parser
			call initf (chain_analysis%n_posterior_combinations)

			write(*,*)
			write(*,*) 'Posterior for combinations of parameters'
			write(*,*)
			allocate(chain_function(chain_analysis%nlength))
				
			do i = 1, chain_analysis%n_posterior_combinations
			
				allocate(indx(chain_analysis%posterior_fun(i)%n_variables))

				write(*,FMT='(A,A)') ' Function => ', chain_analysis%posterior_fun(i)%posterior_function
			
! Locate the index of the variables that enter into the function
				do j = 1, chain_analysis%posterior_fun(i)%n_variables
					indx(j) = locate_parameter_byname(chain_analysis, fich_chain, &
						chain_analysis%posterior_fun(i)%variables(j))
					write(*,FMT='(A,I2,2X,A,A,I3)') '     Variable ', j,&
					 	trim(adjustl(chain_analysis%posterior_fun(i)%variables(j))), &
						' - index= ', indx(j)
				enddo

! Parse the function generting bytecode
				call parsef(i, chain_analysis%posterior_fun(i)%posterior_function, &
					chain_analysis%posterior_fun(i)%variables)
	
! Then evaluate the chain for the new parameter
				do j = 1, chain_analysis%nlength				
					chain_function(j) = evalf(i,chain_analysis%chain(indx,j))				
				enddo
	

				deallocate(indx)

				n = oned_histogram_size(chain_function)

! Variables
				allocate(x(n))
				allocate(yGauss(n))
				allocate(yStep(n))

! Write the posterior marginal distributions
				call oned_histogram(chain_function, x, yGauss, yStep)
				call conf_limits_ccf(x, yStep, 68.268949d0, est, errup, errdown)
! 				most_probable(i) = est
				write(13,*) trim(adjustl(chain_analysis%posterior_fun(i)%posterior_function))
				write(13,*) est, errdown, errup
				write(*,FMT='(A,I1,A,F9.4,A,F9.4,A,F9.4,A)') 'Combination ', i, ' (1s) : E(x)=', &
					est, '  +', errdown, '  -', errup
				call conf_limits_ccf(x, yStep, 95.449972d0, est, errup, errdown)
				write(13,*) est, errdown, errup
				write(*,FMT='(A,I1,A,F9.4,A,F9.4,A,F9.4,A)') 'Combination ', i, ' (2s) : E(x)=', &
					est, '  +', errdown, '  -', errup
				n = size(x)
				write(12,*) i, n
				do j = 1, n
					write(12,*) x(j), yGauss(j), yStep(j)
				enddo
				deallocate(x)
				deallocate(yGauss)
				deallocate(yStep)

! Calculate the percentiles
				percentiles = percentile_summary(chain_function)
				write(15,FMT='(7(F12.4))') percentiles
			
			enddo

			deallocate(chain_function)
		endif

! Calculate average ln L for calculating Kullback-Leibler distance
		chain_analysis%avg_lnL = sum(chain_analysis%chain(chain_analysis%nparams+1,:)) / chain_analysis%nlength
		
		write(13,*) 'Evidence, <ln L>, Kullback-Leibler divergence, chi^2(min)'
		write(13,*) chain_analysis%evidence, chain_analysis%avg_lnL, &
			-chain_analysis%evidence + chain_analysis%avg_lnL, chi2_min

		write(*,*) 'Evidence, <ln L>, Kullback-Leibler divergence, chi^2(min)'
		write(*,*) chain_analysis%evidence, chain_analysis%avg_lnL, &
			-chain_analysis%evidence + chain_analysis%avg_lnL, chi2_min
			
			
		close(12)
		close(13)
		
	end subroutine analyze_chains_multinest

end module chain_analysis_module
