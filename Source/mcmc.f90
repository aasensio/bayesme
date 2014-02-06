program mcmc
use vars
use read_config_mod
use markov_chain_montecarlo
use multinest_mod

integer :: k
real(kind=4) :: dum
	
! 	idum = -240
! 	dum = randomu(idum)
	
	call init_random_seed()
	
	call factrl
	
	call read_config(model,linea,inversion,Observation)
	
	!call read_config_inversion(inversion,Observation)
		
! 	if (start_or_restart <= 2) then
! 		call do_mcmc(inversion,Observation)
! 	else
	call do_multinest(inversion,Observation)
! 	endif
					
end program mcmc
