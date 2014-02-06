! Parameters for the multinest algorithm
module multinest_mod
use vars
use Nested
use markov_chain_montecarlo, only : fill_model_uniformvector_new, compute_likelihood, locate_filling_factors,&
	fill_model
use atomic_functions, only : synthesize
use chain_analysis_module
implicit none

!dimensionality
	integer :: sdim
      
      
!priors on the parameters
!uniform priors (-6,6) are used for all dimensions & are set in main.f90
   real(kind=8), allocatable :: prior_range(:,:) !spriorran(sdim,2)

   integer, allocatable :: nest_pWrap(:) ! Wraparound parameters
      
! Parameters for Nested Sampler
	
!whether to do multimodal sampling
	logical :: nest_mmodal
 	parameter(nest_mmodal=.true.)
	
!max no. of live points
   integer :: nest_nlive
	parameter(nest_nlive=300)
      
!tot no. of parameters, should be sdim in most cases but if you need to
!store some additional parameters with the actual parameters then
!you need to pass them through the likelihood routine
	integer :: nest_nPar

	logical :: nest_ceff
	parameter(nest_ceff = .TRUE.)
      
!seed for nested sampler, -ve means take it from sys clock
	integer :: nest_rseed
	parameter(nest_rseed=-1)
      
!evidence tolerance factor
   real(kind=8) :: nest_tol
   parameter(nest_tol=0.1)
      
!enlargement factor reduction parameter
	real(kind=8) :: nest_efr
   parameter(nest_efr=0.3d0)
      
!root for saving posterior files
	character*100 nest_root	
	
!no. of iterations after which the ouput files should be updated
	integer nest_updInt
	parameter(nest_updInt=100)

!null evidence (set it to very high negative no. if null evidence is unknown)
	real*8 nest_Ztol
	parameter(nest_Ztol=-1.d90)
      
!max modes expected, for memory allocation
  	integer nest_maxModes
  	parameter(nest_maxModes=20)
      
!no. of parameters to cluster (for mode detection)
  	integer nest_nClsPar
!   	parameter(nest_nClsPar=10)
      
!whether to resume from a previous run
  	logical nest_resume
  	parameter(nest_resume=.false.)
      
!feedback on the sampling progress?
  	logical nest_fb
  	parameter(nest_fb=.true.)

	type(statistics) :: stat
	type(modelo_type), pointer :: proposed(:)

	real(kind=8) :: maxlhood, chi2_min

contains

!----------------------------------------------------------------------
	subroutine nest_Sample
   integer :: nclusters,context !total number of clusters found
   integer :: maxNode !variables used by the posterior routine
   
   	call nestRun(nest_mmodal,nest_ceff,nest_nlive,nest_tol,nest_efr,sdim,nest_nPar, &
   		nest_nClsPar,nest_maxModes,nest_updInt,nest_Ztol,nest_root,nest_rseed, nest_pWrap, &
   		nest_fb,nest_resume,getLogLike,dumper,context)

	end subroutine nest_Sample
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! Dumper is called after updInt*10 iterations
!----------------------------------------------------------------------
	subroutine dumper(nSamples, nlive, nPar, physLive, posterior, paramConstr, maxLogLike, logZ)

        implicit none

        integer :: nSamples                                ! number of samples in posterior array
        integer :: nlive                                   ! number of live points
        integer :: nPar                                    ! number of parameters saved (physical plus derived)
        real(kind=8), pointer :: physLive(:,:)      ! array containing the last set of live points
        real(kind=8), pointer :: posterior(:,:)     ! array with the posterior distribution
        real(kind=8), pointer :: paramConstr(:)     ! array with mean, sigmas, maxlike & MAP parameters
        real(kind=8) :: maxLogLike                     ! max loglikelihood value
        real(kind=8) :: logZ                           ! log evidence

	end subroutine dumper

!----------------------------------------------------------------------
! Wrapper around Likelihood Function
! Cube(1:n_dim) has nonphysical parameters
! scale Cube(1:n_dim) & return the scaled parameters in Cube(1:n_dim) &
! additional parameters in Cube(n_dim+1:nPar)
! return the log-likelihood in lnew
!----------------------------------------------------------------------
	subroutine getLogLike(Cube,n_dim,nPar,lnew,context)

   integer :: n_dim,nPar,context
   real(kind=8) :: lnew,Cube(nPar)
   
   !call your loglike function here   
   !lnew=loglike(likeindx,Cube)
   	call slikelihood(Cube,lnew)

	end subroutine getLogLike

!----------------------------------------------------------------------
! Evaluate likelihood
!----------------------------------------------------------------------
	subroutine slikelihood(Cube,slhood)              
	real(kind=8) :: Cube(nest_nPar),slhood, chi2
	real(kind=8) :: temp(sdim),dist,loclik, ff_total
	integer :: i,j
	logical :: physical
	

		Cube(1:sdim) = (prior_range(1:sdim,2)-prior_range(1:sdim,1))*Cube(1:sdim) + prior_range(1:sdim,1)

! If more than one component, the last filling factor is not a free
! parameter but fulfills ff(last) = 1-sum(ff(1:last-1))

		if (number_of_components > 1) then

! Calculate the value of the last filling factor
! and put it into Cube(nest_npar)
			ff_total = 0.d0
			do i = 1, number_of_components
				j = inversion%filling_factor_position(i)
				ff_total = ff_total + Cube(j)
			enddo
			
			if (ff_total /= 1.d0) then
				j = inversion%filling_factor_position(number_of_components)
				Cube(j) = 1.d0 - ff_total + Cube(j)
			endif
		endif

		stat%parameters_proposed = Cube
		
		call fill_model_uniformvector_new(inversion,stat,proposed)
		
		slhood = compute_likelihood(proposed,linea,Observation,inversion,Stokes_Syn,physical,chi2,1.d0)

		if (slhood > maxlhood) then
			maxlhood = slhood
			chi2_min = chi2
			stat%parameters_MAP = stat%parameters_proposed
		endif

		return
			
	end subroutine slikelihood


!----------------------------------------------------------------------
! Do MULTINEST MCMC
!----------------------------------------------------------------------
	subroutine do_multinest(inv,Observation)
	type(inversion_type) :: inv
	type(stokes_type) :: observation
	real(kind=8) :: chi2
	logical :: physical
	integer :: i

		nest_nPar = inv%nparams_invert
		sdim = nest_npar		

! If there are several components, one of the filling factors is not a
! free parameter but is related to the sum of the others
		if (number_of_components > 1) then
			sdim = nest_npar - 1
		endif

! Number of parameters to cluster
		nest_nClsPar = sdim

		nest_root = fich_chain

		maxlhood = -1.d100

! Allocate memory for the prior ranges
		allocate(prior_range(nest_nPar,2))

! For periodic boundary conditions
		allocate(nest_pWrap(sdim))
		nest_pWrap = 0

! Set them equal to the prior ranges read before
		prior_range(:,1) = inv%range(1:inv%nparams_invert,1)
		prior_range(:,2) = inv%range(1:inv%nparams_invert,2)

! Set number of parameters to invert
		stat%nparams = inv%nparams_invert

! Get some memory
		allocate(stat%parameters_proposed(stat%nparams))
		allocate(stat%parameters_MAP(stat%nparams))
		allocate(proposed(number_of_components))

! Identify filling factors
		proposed = model
! 		call locate_filling_factors(inv)

! Sample
 		call nest_Sample

! Save MAP profile
		stat%parameters_proposed = stat%parameters_MAP
 		call fill_model_uniformvector_new(inv,stat,proposed)
 
 		call synthesize(proposed,linea,Observation,Stokes_Syn)
 		
 		open(unit=12,file=trim(adjustl(fich_chain))//'.stokes_map',&
 			action='write',status='replace')
 		write(12,*) maxlhood, chi2_min
 		write(12,*) Stokes_Syn%nlambda
 		do i = 1, Stokes_Syn%nlambda
 			write(12,FMT='(F10.4,2X,4(E12.4))') Stokes_Syn%lambda(i),&
 				Stokes_Syn%stokes(1,i), Stokes_Syn%stokes(2,i), &
 				Stokes_Syn%stokes(3,i), Stokes_Syn%stokes(4,i)
 		enddo
 		close(12)

! Analyze chains
		call analyze_chains_multinest(chain_analysis,fich_chain,stat%parameters_MAP,chi2_min)

		write(*,*) 'MAP values'
		write(*,*) stat%parameters_MAP
		write(*,*)

 		print *, 'Maximum posterior value = ', maxlhood
 		print *, 'Minimum chi^2 = ', chi2_min
 		print *, 'Reduced minimum chi^2 = ', chi2_min / (4.d0*Stokes_Syn%nlambda)

 		inv%chi2_min = chi2_min

		deallocate(prior_range)		
		
	end subroutine do_multinest

end module multinest_mod
