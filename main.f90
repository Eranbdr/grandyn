

PROGRAM GranLayer
	use xml_data_config_template
	use mycommons
	IMPLICIT NONE

	!include "mpif.h"
	

!    *******************************************************************
!    ** FORTRAN GRANULAR DYANMICS PROGRAM*         
!    **  
!    **   BUILDS AND/OR SHEARS A LAYER OF CIRCULAR GRAINS
!    **     (NOTE NOT SET UP FOR MULTI-PARTICLE GRAINS)
!    **     PERIODIC IN THE HORIZONTAL DIRECTION
!    **     VARIOUS CONDITIONS ON THE TOP AND BOTTOM WALLS
!    **     C    *******************************************************************
!    ** PRINCIPAL VARIABLES:                   
!    **                                             
!    **  N		NO. OF GRAINS
!    **  RX,RY	GRAIN POSITIONS
!    **  VX,VY	GRAIN VELOCITIES
!    **  FX,FY	GRAIN FORCES
!    **  RADIUS GRAIN RADII


!    ** BOUNDARIES ARE LISTED IN THE ORDER: BOTTOM,TOP,LEFT,RIGHT
!    **  NBOUND(4)	NO. OF GRAINS IN EACH BOUNDARY
!    **  IB(4)	TYPE OF EACH BOUNDARY

!    **  HEAD	HEAD OF CHAIN FOR EACH CELL
!    **  LIST	LINKED LIST OF GRAINS
!    **  MAP		LIST OF NEIGHBORING CELLS (ONLY UPPER-RIGHT HALF)
!    **  MX 		TOTAL NUMBER OF CELLS IN X-DIRECTION
!    **  MY 		TOTAL NUMBER OF CELLS IN Y-DIRECTION

!    **  SIG		AVERAGE GRAIN DIAMETER (SET TO 1.)

!    ** **INPUT PARAMETERS**
!    **  BOXX, BOXY  SIZE OF BOX IN X AND Y (RELATIVE TO AVERAGE
!    **	    GRAIN DIAMETER)
!    **  SIGB	AVERAGE BOUNDARY GRAIN DIAMETER 
!    **  shear_spring_coef     SPRING CONSTANT FOR LEAF SPRING
!    **  GAMMA 	DAMPING PARAMETER FOR NORMAL VELOCITY
!    **  GAMSHEAR    DAMPING PARAMETER FOR TANGENTIAL VELOCITY
!    **  FRICTION   COEFFICIENT OF FRICTION


!  DWS 11/9 
!  ^^^^non-dimensionalization^^^^
! choosing values for :
!  k = stiffness between particles 
!  m = mass of standard particle
!  length scale , x0, is one standard grain diameter

! the following scales apply:
!  velocity scale,v0 = x0 * sqrt(k/m)
!  time scale, t0 = x0/v0
!  resulting force scale, f0 = k*x0
!  resulting energy scale, e0 = k*x0*x0

!  non-dim equations become
!    f = dabs(overlap) - gamma*(velocity difference) 
!                     - rho*(unit vertical vector)
!
!    x(t+dt) = x(t) + v(t)*dt  + .5*dt*dt*f(t)
!
!    v(t+dt) = v(t) + .5*dt*(f(t) + f(t+dt))
!
!   where
!   damping parameter gamma = gamma0/sqrt(k*m)
!			    = gamma0/grainradius/sqrt(density*pi*k)
!   buoyancy parameter rho  = g*m/(k*x0)
!                           = g *(grainradius)**2/(acoustic velocity)**2/x0
!
!   g=acceleration of gravity
!   gamma0 = some inherent damping(kg-m/s)
!
!      note: rho (input as g(1) and g(2), will be small (order 10-3 or less),
!            gamma will be of order 1 for well-damped systems
!       

!  critical time step = grain radius/acoustic velocity (from Mora and Place)
!                     ~ r/( sqrt(9/8 * k/m) *r)
!
!  critical time step scaled by t0 = sqrt(8/9)
!  Mora and Place recommend a time step of 0.2*critical, or about .15
!

!  DWS 11/9 

	
        
!    *******************************************************************
	

	integer i, j, k
	real*8 rsmall
	real rand1
	integer :: values(1:8)
!    ** READ IN INITIAL PARAMETERS **
	
		  
	print*, ' **  PROGRAM START **'		
	degree_to_radian	= pi/180d0
	step = 1
	output_burst_step = 0
	if (iargc() .eq. 0) then
		config_file = default_config_file
	else
		call getarg(1, config_file)
	endif
	
	call read_xml_file_config_template(config_file , 0, config_status )
	if ( config_status) then
		WRITE(*,*) 'failed to read ', config_file
		stop
	endif
	WRITE(*,*) 'loaded ', config_file
	! call write_xml_file_config_template( 'config/out_config.xml', 2000 )
	
	! we only deal with a single population
	frac1 = 1
	write(fnn,fmt='(i12.12)') 0
	
	call read_init_phase_params()                
	write(*,*) 'got init params, set phase ', currphase
	write(*,*) 'writing output to ' , TRIM(output_directory)
	call write_xml_file_config_template( TRIM(output_directory)//'/config.xml', 20000 )
	!open_file( len(trim(get_output_directory), randno_file, .true.)
	!open_file( len(trim(get_output_directory), randno_file, .false.)
!   generate random number seed, or use one previous stored in file 'randno'
	if (irandfile) then
		open(unit=randno_file,file=TRIM(output_directory)//'/randno')
		read(randno_file,*) randum
		close(randno_file)
	else 
		call date_and_time(values=values)
		randum = -values(8)
		open(unit=randno_file,file=TRIM(output_directory)//'/randno', status='replace')
		write(randno_file,*) randum
		close(randno_file)
	endif
      
    call srand(randum) ! seed for the rand  
	amp_rand = (rand(0) - 0.5) / 4	
	random_flag = .true.
! presy0 is the config input
! presy is used for: a. crushing, b. periodic change of external stress
 
	presy0=presy
	last_contact_step = 0
!   start at confining pressure of 1.0e-5, increase pressure gradually
 
	if (phigoal .gt. 0.) then 
		presyfinal=presy ! take phigoal?
		presy=1.0d-5
	endif
	
	
! &&&&&INIT&&&&&INIT&&&&&INIT&&&&&INIT&&&&&INIT&&&&&INIT&&&&&INIT
!  INITIALIZATION SECTION: BUILD OR READ IN GRAINS
           
! read a restart file
	if (irstart) then
		write(*,*) 'restart flag on'
		call read_restart()
		call dumpgrains(afrac1)

		!  adjust number of cells to account for final box size and grain sizes
		bg=MAXVAL(radius)*2.
		!     ******************************************************************
		boxx=boxx*(1-presy)

		! maximum number of cells allowed for neighbour calculation
		mxmax=int(xright-xleft)/bg-1
		mymax=max(1,int((ytop-ybot)/bg))
		print*, '***  read restart file   ****'
		print*, 'biggest grain =',bg
		write(luout, *) 'read restart file, biggest grain =',bg
		my=min(my,mymax)
		mx=min(mx,mxmax)

		print*, 'mx and my =', mx,my

		write(*,*) 'n,nbound(4)',n,nbound(4)
		write(*,*) 'position of boundaries:',ybot,ytop,xleft,xright
		write(*,*) 'boxx,boxy,mx,my',boxx,boxy,mx,my
		write(*,*) 'tau',tau
		
		! zero out tau at beginning of run
		if ( .not. iztau) tau=0.
		
		if (.not. izmom) then
			!  zero out the initial momentum, for a dead restart
			v(1,1:n)=0.
			v(2,1:n)=0.
			w(1:n)=0.
		endif

!  zero out the slip info, to avoid crazy springs during
!   no friction relaxation
		if ( .not. izslip) contft(1:contactknt)=0.

	else ! no restart file, build all the grains
	 
		tau=0.
		
		v(1,nbound(1):n)=fb(2)
		v(2,nbound(1):n)=0.
		w(nbound(1):n)=0.
		
		contft(1:contactknt)=0.
		
		grainarea=0d0
		call initbound()
		boxvol=boxx*boxfy
		write (*,*) ' partgen params = ', phigoal,boxvol,grainarea,alap
		call partgen()

		print*, 'finished generating ',n,' grains'
		call dumpgrains(afrac1)
		! no fluid in init step
	endif

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! AFTER GENERATING SYSTEM OR READING RSTART FILE:

!  set up factor for getting 1/mass of a particle (particle of
!    radius=0.5 has mass=1
       
	radinv2(1:n)=(.5/radius(1:n))**3
	mass(1:n)=4.*(radius(1:n)**2)

!  calculate the area taken up by particles,
	grainarea=0.
	do i=1,n
		if (bdgrain(i) .eq. 0) then
			grainarea = grainarea+radius(i)**2
		else if ((ib(2) .ge. 0) .or. (bdgrain(i) .eq. 1)) then
			grainarea = grainarea + .5*radius(i)**2
		endif
	enddo
	grainarea = grainarea*pi
	write(luout,*) 'initial grain area in main=',grainarea
	!print*, 'initial porosity = ', 1. - grainarea/((xright-xleft)*(ytop-ybot))	  
  
! find smallest grain for calculating time step
	rsmall=MINVAL(radius(1:n))

! Calculate dt
	print*, 'smallest radius before dt calculation: ',rsmall
	dt=dt*((2*rsmall)**1.5)
	dt2 = 0.5*dt  
	dtsq2 = dt*dt2
	WRITE(luout,*) 'rsmall ', rsmall
	WRITE(luout,*) '*********************************'
	WRITE(LUOUT, *) 'NUMBER OF GRAINS = ',N
	WRITE(luout, *)' NUMBER OF BOUNDARY GRAINS:' 
	WRITE(luout, *)'  BOTTOM: ',NBOUND(1)
	WRITE(luout, *)'     TOP: ',NBOUND(2)-NBOUND(1)
	WRITE(luout, *)' NO. OF CELLS =', MX,' x ', MY
	WRITE(luout, *)' '
	WRITE(luout, *)'OUTPUT FREQUENCY = ',  iprint
	WRITE(luout, *)' '
	WRITE(luout, *)'NORMAL STIFFNESS = 1.'
	WRITE(luout, *)'TANGENTIAL STIFFNESS = ', shear_spring_coef
	WRITE(luout, *)'NORMAL DAMPING COEF= ', gamma
	WRITE(luout, *)'FRICTION        = ', friction
	WRITE(luout, *)' '
	WRITE(luout, *)' TIME step        = ', dt
	WRITE(luout, *)' EXTERNAL LOADING       = ', presy
	WRITE(luout, *)' BODY FORCES       = ', g(1),g(2)
	WRITE(luout, *)' TOP BOTTOM LEFT RIGHT       = ', ytop, ybot,xleft,xright


!   index of first interior grain
	indint=nbound(4)+1 
	
! END OF BUILDING AND READING
! &&&&&INIT&&&&&INIT&&&&&INIT&&&&&INIT&&&&&INIT&&&&&INIT&&&&&INIT


!    ** CALCULATE INITIAL VALUES, NEED TO FORM INITIAL  LIST **
!    ** AND GET FORCES AT T=0 **

	if ( iprint .le. 0 ) iprint = nstep
	
!      ** INITIALIZE THE MAP OF NEIGHBOR CELLS
	call maps
	call links
	call dumplinklist


	
! if freezing(or resetting) boundaries, also need to freeze velocities,
!  otherwise, there will be slight movement during first 
!  time step at the previous velocity
	if (ib(1).eq.3) v(2,1:nbound(1))=fb(1)
	if (ib(2).eq.3) v(2,nbound(1)+1:nbound(2))=fb(2)
	if (ib(1) .eq. 4) then
		v(1,1:nbound(1))=fb(1)
		v(1,nbound(1)+1:nbound(2))=fb(2)
	endif
       
!  top wall pulled by shear spring
!   initialize the distance moved by the wall, dwall
!   and the pulling spring, dspring 
	if (ib(2) .eq. 5) then
		if ( .not. iztau) then
			dwall=0.
			dspring=0.
		endif
		uspring=fb(2)
	endif

!  treat the top row as a rigid block of some weight:
!   applied boundary force (presy)/gravity (g(2)) gives the
!   height of the overlying mass in particles. Give the
!   boundary particles effectively the mass of the overlying
!  column
	if (itmass .eq. 1) then
		topmass=0.00463/kspring
	else
		topmass=1.
	endif

! &&&&&BOUNDARY&&&&BOUNDARY&&&&&BOUNDARY&&&&&BOUNDARY&&&&&BOUNDARY&&&&&


! 11/30/01 - forces are now read in from restart file
! GET INITIAL FORCES AND PLOT
	print*, 'calling first force'
	!call links
	print*, 'number of contacts =',contactknt
	call getneighbors()
	if (ifluid) then
	   toprint = 0
	   if (mod(step,iprint).eq.0) toprint = 1
	  
	   call smooth(toprint)
	   do i = 1,TRATIO
		  call updatepressure(toprint,i)
	   enddo
	endif
    

!  open monitor files 
	open(unit=coordination_file,file=TRIM(output_directory)//'/coordination')
	open(unit=porosity_file,file=TRIM(output_directory)//'/porosity')
	open(unit=stresses_file,file=TRIM(output_directory)//'/stresses')
	open(unit=momentum_file,file=TRIM(output_directory)//'/momentum')
	open(unit=energy_file,file=TRIM(output_directory)//'/energy')
	if (ib(2) .eq. 5) open(unit=spring_file,file=TRIM(output_directory)//'/spring')
	open(unit=strain_file,file=TRIM(output_directory)//'/strain')
	open(unit=topfriction_file,file=TRIM(output_directory)//'/topfriction')
	open(unit = dilatation_file, file=TRIM(output_directory)//'/dilatationa')
	open(unit=flux_file,file=TRIM(output_directory)//'/flux_'//TRIM(output_directory)//'.csv')
	open(unit=vol_flux_file,file=TRIM(output_directory)//'/vol_flux_'//TRIM(output_directory)//'.csv')
	open(unit=grains_per_layer_file,file=TRIM(output_directory)//'/grains_per_layer_'//TRIM(output_directory)//'.csv')
	open(unit=velocity_pro_file,file=TRIM(output_directory)//'/velocity_profile_'//TRIM(output_directory)//'.csv')
	write(grains_per_layer_file, '(a)', advance='no') ' , internal_top-ybot'
	write(flux_file, '(a)', advance='no') ' , average ek'
	write(vol_flux_file, '(a)', advance='no') ' , average ek'
	write(velocity_pro_file, '(a)', advance='no') ' , max comp'
	do i=1,number_of_layers 
		write(grains_per_layer_file, '(", grains_flux_calc_per_layer (over flush time) layer ",i2)', advance='no') i
		write(grains_per_layer_file, '(", grains_per_layer (now) layer ",i2)', advance='no') i
		write(flux_file, '(", flux ",i2)', advance='no') i
		write(vol_flux_file, '(", flux ",i2)', advance='no') i
		write(velocity_pro_file, '(", vel ",i2)', advance='no') i
	enddo
	write(grains_per_layer_file, '(a)') ', sum'
	write(flux_file,'(a)') ', total flux'
	write(vol_flux_file,'(a)') ', total flux'
	write(velocity_pro_file,'(a)') ', ave vel'
	nbcheck = 20

!    *******MAIN*******MAIN*******MAIN*******MAIN*******MAIN*******MAIN
!    ** MAIN LOOP BEGINS                                            **
!    *******************************************************************
       
	WRITE(*,'(//1X,''**** START OF DYNAMICS ****'')')
	
	! main loop for init 
	stop_condition_reached = 0
	phase_step = 1
	! don't run init phase if from restart file
	if ( .not. irstart) then 
		write(*,*) 'setting init forces'
		
		curpresy = presy
		call force()
		call dump_forces()
		DO  step = 1, nstep	
			write(fnn,fmt='(i12.12)') step
			call iteration()
			if (stop_condition_reached .eq. stop_condition_limit) then
				WRITE(*,*) '******   stop condition reached in init*******'
				exit
			endif
			phase_step = phase_step + 1 
		enddo
		if (stop_condition_reached .eq. stop_condition_limit) step = step + 1 ! exit skips the advance
		if (stop_condition_reached .lt. stop_condition_limit &
				.and. terminate_on_condition_not_met) then
			WRITE(*,*) '**** END OF PHASE, stop condition not met in init, terminating   **** '
			stop
		endif
		
		 ! find the top of the grains
		if (ib(2) .ge. 0 .and. phigoal .gt. 0.) then
			topg=r(2,1)+radius(1)
			do i=nbound(2)+1,n
				topg=max(topg,r(2,i)+radius(i))
			enddo
			area = (xright-xleft)*(topg-ybot)
			phi = 1.-grainarea/area
			!          print*, 'porosity=',phi
			wmax=0.
			do i=nbound(1)+1,nbound(2)
				wmax=max(wmax,radius(i))
			enddo
			do i=nbound(1)+1,nbound(2)
				r(2,i)=topg+wmax
			enddo
			presy=presyfinal
			phigoal=0.
			ytop = r(2,nbound(2))
		endif
	endif

    
		
	WRITE(*,'(/1X,''**** END OF INITIALIZATION **** ''//)')

	call pscriptplot(0)
	call write_restart()
    	internal_top = maxval(r(2,(nbound(2)+1):n)) + maxval(radius((nbound(2)+1):n))
	write(*,*) 'start with phase of phases ' ,currphase,  size(phase)
    
! loop through phases	 
	do currphase=currphase, size(phase)  
		stop_condition_reached = 0
		phase_start_height = internal_top - ybot
		write (*,*) 'reading phase ', currphase, ' step ' , step
		write(*,*) 'phase height ', phase_start_height
		write(*,*) 'solid fraction = ', grainarea/((xright-xleft)*phase_start_height)
		write(*,*) 'l = ', xleft
		write(*,*) 'r = ', xright
		phase_step = 1
		call read_phase_parameters(currphase)
		
		!open(unit=layer_change_file,file=TRIM(output_directory)//'/layer_change'//TRIM(file_postfix)//'.csv')
		!write(layer_change_file,*) ', grain, start, end'
		write (*,*) 'g(2) = ', g(2)
		write (*,*) 'g(1) = ', g(1)
		write (*,*) 'ib 1 ', ib(1)
		write (*,*) 'ib 2 ', ib(2)
		write (*,*) 'ib 3 ', ib(3)
		write (*,*) 'ib 4 ', ib(4)
		write (*,*) 'mx ', mx
		write (*,*) 'my ', my
		write (*,*) 'dt ', dt
		

		! main loop for phase, step is total steps from init
		!write(flux_file,*) file_postfix
		DO  step=step, nstep
			
			write(fnn,fmt='(i12.12)') step
			call iteration()
			if (stop_condition_reached .eq. stop_condition_limit) then
				WRITE(*,*) '******   stop condition reached phase ', currphase, phase_step,' phase_step   **** '
				

				exit
			endif
			if (skip_phase) then
				WRITE(*,*) '******   skip phase ', currphase, phase_step,' phase_step   **** '
				OPEN (UNIT=555, FILE=TRIM(output_directory)//'/skip', STATUS="OLD")  
				CLOSE (UNIT=555, STATUS="DELETE")
				skip_phase = .false.
				exit
			endif
			if (output_burst) then
				OPEN (UNIT=111, FILE=TRIM(output_directory)//'/output', STATUS="OLD")  
				read(111, *) output_burst_step, output_burst_interval, output_burst_id, output_burst_id_step
				CLOSE (UNIT=111, STATUS="DELETE")
				WRITE(*,*) '******   output burst ', output_burst_step, output_burst_id, output_burst_id_step
				output_burst = .false.
			endif
			phase_step = phase_step + 1
		enddo
		if (stop_condition_reached .eq. stop_condition_limit) step = step + 1 ! exit skips the advance
		if (stop_condition_reached .lt. stop_condition_limit &
				.and. terminate_on_condition_not_met) then
			WRITE(*,*) '**** END OF PHASE, terminating, stop condition not met: phase ', currphase, ', &
			 stop_condition_reached ',stop_condition_reached, ' phase_step ',phase_step
			stop
		endif
		WRITE(*,*) '**** END OF PHASE ', currphase, '   **** '
		call pscriptplot(0)
		call write_restart()
		write (*,*) 'back from ps and restart'
		!close(layer_change_file)
	enddo
	WRITE(*,'(/1X,''**** END OF DYNAMICS **** ''//)')
	close(coordination_file)
	close(porosity_file)
	close(stresses_file)
	close(momentum_file)
	close(energy_file)
	if (ib(2) .eq. 5) close(spring_file)
	close(strain_file)
	close(topfriction_file)
	close(dilatation_file)
	close(grains_per_layer_file)
	close(flux_file)
	close(vol_flux_file)
	close(velocity_pro_file)
	
	STOP

END

	real function rand1(idum)
		integer idum,ia,im,iq,ir,ntab,ndiv
		real am,eps,rnmx
		parameter(ia=16807,im=2147483647,am=1./im,iq=127773,ir=2836,ntab=32)
		parameter (ndiv=1+(im-1)/ntab,eps=1.2e-7,rnmx=1.-eps)
		
		integer j,k,iv(ntab),iy
		save iv,iy
		data iv /ntab*0/, iy /0/

		if (idum .le. 0 .or. iy .eq. 0) then
			idum=max(-idum,1)
			do j=ntab+8,1,-1
				k=idum/iq
				idum=ia*(idum-k*iq)-ir*k
				if (idum .lt. 0) idum=idum+im
				if (j .le. ntab) iv(j)=idum
			enddo
			iy=iv(1)
		endif

		k=idum/iq
		idum=ia*(idum-k*iq)-ir*k
		if (idum .lt. 0) idum=idum+im
		j=1+iy/ndiv
		iy=iv(j)
		iv(j)=idum
		rand1=min(am*iy,rnmx)
		return
	end

	subroutine read_init_phase_params()
		use xml_data_config_template
		use mycommons
		!	init only params
		! restart
		irstart = init_phase%restart_file
		
		if (irstart) then
			restart_file_name = init_phase%restart_file_name
			iztau  = init_phase%use_restart_tau
			izmom  = init_phase%use_restart_velocity
			izslip  = init_phase%use_restart_slip			
		else 
			iztau  = .false.
			izmom  = .false.
			izslip  = .false.	
		endif
		
		ipulse = 0
		fstress = 0.
		fbc = 0
		ifluid = .false.
		  
		dt = init_phase%dt
		irandfile = init_phase%use_rand_file
		currphase = init_phase%start_from_index
		write(*,*) '**************    ', currphase
		! dimensions
		boxx = init_phase%dimensions%horizontal
		boxy = init_phase%dimensions%vertical
		mx = init_phase%dimensions%horizontal_cells
		my = init_phase%dimensions%vertical_cells
		boxfy = init_phase%dimensions%goal_vertical
		phigoal = init_phase%dimensions%goal_porosity
		
		! grain population size
		mean1 = init_phase%grain_population_size%mean_diameter
		std1 = init_phase%grain_population_size%standard_deviation
		logr = init_phase%grain_population_size%smallest
		higr = init_phase%grain_population_size%largest
		sigb = (logr + higr) / 2 
		
		! grain properties
		normal_spring_coef = init_phase%grain_properties%normal_spring
		shear_spring_coef = init_phase%grain_properties%shear_spring
		friction = init_phase%grain_properties%friction
		gamma = init_phase%grain_properties%normal_damping_coefficient
		gamshear = init_phase%grain_properties%shear_damping_coefficient
		iroll = init_phase%grain_properties%roll
		hertz_mindlin = init_phase%grain_properties%hertz_mindlin
		
		! output
		iprint = init_phase%output%plot_every
		gen_restart	= init_phase%output%restart_file_every
		iprintf = init_phase%output%dump_system_data_every
		igif = init_phase%output%generate_jpgs
		file_postfix = init_phase%output%postfix
		number_of_layers = init_phase%output%number_of_layers
		output_directory = init_phase%output_dir
		CALL system('mkdir -p '//TRIM(output_directory))
		
		! forces
		presx = init_phase%forces%horizontal_stress
		presy = init_phase%forces%vertical_stress
		g(1) = init_phase%forces%horizontal_gravity
		g(2) = init_phase%forces%vertical_gravity
		
		
		
		! progress
		nstep = init_phase%progress%steps
		target_ek = init_phase%progress%stop_when_kinetic_energy
		target_flux_change = init_phase%progress%stop_when_flux_change
		stop_condition_limit = init_phase%progress%number_of_consecutive_checks
		terminate_on_condition_not_met = init_phase%progress%kill_if_condtion_not_met
		check_condition_every = init_phase%progress%check_condition_every
		
		! ib = phase(iphase_index)%walls
		ib(1) = init_phase%walls%bottom
		ib(2)  = init_phase%walls%top
		ib(3) = init_phase%walls%left
		ib(4) = init_phase%walls%right
		
		write (*,*) 'g(1) = ', g(1)
		write (*,*) 'g(2) = ', g(2)
		write (*,*) 'ib 1 ', ib(1)
		write (*,*) 'ib 2 ', ib(2)
		write (*,*) 'ib 3 ', ib(3)
		write (*,*) 'ib 4 ', ib(4)
		write (*,*) 'mean1 ', mean1	
		write (*,*) 'std1 ', std1
		write (*,*) 'logr ', logr
		write (*,*) 'higr ', higr
		write (*,*) 'boxx ', boxx
		write (*,*) 'boxy ', boxy
		write (*,*) 'mx ', mx
		write (*,*) 'my ', my
		write (*,*) 'dt ', dt
	
	end subroutine
	
	subroutine read_phase_parameters(iphase_order_index)
	! the input xml may not be sorted by order_index
	! when we execute the phase, it must be in the specified order
	! init is order_index ==
		use xml_data_config_template
		use mycommons
		integer i
		real*8 tilt
		
		integer iphase_order_index, iphase_index
		iphase_index = -1
		! start with init : i =0
		write(*,*) "Reading configuration for phase order index ", iphase_order_index
		

		do i = 1,size(phase)
			if ( phase(i)%order_index == iphase_order_index ) then
				iphase_index = i
				exit
			endif
		enddo
		if (iphase_index .lt. 0) then
			write(*,*) "Failed to load configuration for actual phase order index ", iphase_order_index
			stop
		endif
		
		! output
		iprint = phase(iphase_index)%output%plot_every
		gen_restart	= phase(iphase_index)%output%restart_file_every 
		iprintf = phase(iphase_index)%output%dump_system_data_every
		igif = phase(iphase_index)%output%generate_jpgs
		number_of_layers = phase(iphase_index)%output%number_of_layers
		
		
		! forces
		presx = phase(iphase_index)%forces%horizontal_stress
		presy = phase(iphase_index)%forces%vertical_stress
		tilt = phase(iphase_index)%forces%tilt
		if (tilt .ne. 0) then 
			write(*,*) "Reading Gs - vertical gravity base " , phase(iphase_index)%forces%vertical_gravity
			

			g(1) = -phase(iphase_index)%forces%vertical_gravity*dsin(tilt*degree_to_radian)
			write(*,*) "g1 " , g(1)
			

			g(2) = phase(iphase_index)%forces%vertical_gravity*dcos(tilt*degree_to_radian)
			write(*,*) "g2 " , g(2)
			

			!write( file_postfix, '(f4.1)' ) tilt
			!file_postfix = '_tilt_'//trim(adjustl(file_postfix))//phase(iphase_index)%output%postfix
			file_postfix = phase(iphase_index)%output%postfix
			write(*,*) 'postfix' , file_postfix
			

		else
			write(*,*) 'else'
			

			g(1) = phase(iphase_index)%forces%horizontal_gravity
			write(*,*) "g1 " , g(1)
			

			g(2) = phase(iphase_index)%forces%vertical_gravity
			write(*,*) "g2 " , g(2)
			

			file_postfix = phase(iphase_index)%output%postfix
		endif
		acceleration_amplitude(1) = phase(iphase_index)%forces_modification%horizontal_amplitude
		acceleration_amplitude(2) = phase(iphase_index)%forces_modification%vertical_amplitude
		acceleration_function(1) = phase(iphase_index)%forces_modification%horizontal_function
		acceleration_function(2) = phase(iphase_index)%forces_modification%vertical_function
		acceleration_decay_type(1) = phase(iphase_index)%forces_modification%horizontal_decay_type
		acceleration_decay_type(2) = phase(iphase_index)%forces_modification%vertical_decay_type
		acceleration_decay_limit(1) = phase(iphase_index)%forces_modification%horizontal_decay_limit
		acceleration_decay_limit(2) = phase(iphase_index)%forces_modification%vertical_decay_limit
		
		wave_length(1) = phase(iphase_index)%forces_modification%horizontal_wave_length
		wave_length(2) = phase(iphase_index)%forces_modification%vertical_wave_length
		no_modification_from_last_contact = phase(iphase_index)%forces_modification%no_mod_from_last_contact
		do i = 1,2
			if (acceleration_amplitude(i) .ge. 0.0001 .and. tilt .ne. 0) then
				! amplitude is given in multiplications of g
				acceleration_amplitude(i) = -acceleration_amplitude(i) * phase(iphase_index)%forces%vertical_gravity
			endif
			
			if (wave_length(i) .eq. 0) then 
				angle_change_rate(i) = 0d0
			else
				angle_change_rate(i) = 360d0/wave_length(i)
			endif
			
			write(*,*) i, ' acceleration_amplitude ', acceleration_amplitude(i), ' wave_length ',wave_length(i), ' angle_change_rate ',angle_change_rate(i)
			

			curr_angle(i) = 0d0
		enddo
		
		
		if (acceleration_function(2) .eq. 8) then
			! vertical sine
			acceleration_amplitude(1) = -acceleration_amplitude(2)*dsin(tilt*degree_to_radian)
			acceleration_amplitude(2) = acceleration_amplitude(2)*dcos(tilt*degree_to_radian)
			acceleration_function(1) = 8
			write(*,*) 'sound acceleration_amplitude ', acceleration_amplitude(1), acceleration_amplitude(2)
		endif
		
		
		
		! progress
		nstep = step + phase(iphase_index)%progress%steps - 1
		target_ek = phase(iphase_index)%progress%stop_when_kinetic_energy
		target_flux_change = phase(iphase_index)%progress%stop_when_flux_change
		stop_condition_limit = phase(iphase_index)%progress%number_of_consecutive_checks
		terminate_on_condition_not_met = phase(iphase_index)%progress%kill_if_condtion_not_met
		check_condition_every = phase(iphase_index)%progress%check_condition_every
		! is assigning the xml structure directly wrong?
		! ib = phase(iphase_index)%walls
        ib(1) = phase(iphase_index)%walls%bottom
		ib(2)  = phase(iphase_index)%walls%top
        ib(3) = phase(iphase_index)%walls%left
        ib(4) = phase(iphase_index)%walls%right
	
		! fb = phase(iphase_index)%walls_velocity
		fb(1) =  phase(iphase_index)%walls_velocity%bottom
		fb(2) =  phase(iphase_index)%walls_velocity%top
		fb(3) =  phase(iphase_index)%walls_velocity%left
		fb(4) =  phase(iphase_index)%walls_velocity%right
		
		!<walls_properties> - not supported yet
			 !<spring_coefficient>1.000</spring_coefficient> kspring = phase(iphase_index)%walls_properties%spring_coefficient
	!         <amplitude>1.000</amplitude> amp from wave input = phase(iphase_index)%walls_properties%
			 !<period>1.000</period> period from wave input = phase(iphase_index)%walls_properties%
			 
			 
! fbc are boundary conditions for the fluid pressure at the top and bottom of the box. possibilities:
! 1 = Dirichlet, 2 = Neumann with derivative = 0, 
! OBSELETE 3 = Dirichlet with permeability test i.e. not moving the grains in response to fluid pressure.
! 4 = Dirichlet with moving grains
! 5 = Mixed boundary conditions. Constant pressure gradient at the bottom and constant pressure at top
! 6 = Sin/Cos wave - time dependent pressure on same b.c.
! 7 = Undrained BC that are implemented 

		ifluid = phase(iphase_index)%fluid%with_fluid
		if (ifluid) then
			fstress = phase(iphase_index)%fluid%pressure_ratio
			fbc = phase(iphase_index)%fluid%fluid_boundary_conditions
		else 
			fstress = 0.
			fbc = 0
		endif
		
		WRITE(luout,*) '********************************* phase ', iphase_index
	WRITE(LUOUT, *) 'NUMBER OF GRAINS = ',N
	WRITE(luout, *)' NUMBER OF BOUNDARY GRAINS:' 
	WRITE(luout, *)'  BOTTOM: ',NBOUND(1)
	WRITE(luout, *)'     TOP: ',NBOUND(2)-NBOUND(1)
	WRITE(luout, *)' NO. OF CELLS =', MX,' x ', MY
	WRITE(luout, *)' '
	WRITE(luout, *)'OUTPUT FREQUENCY = ',  iprint
	WRITE(luout, *)' Compression model - Hertz_Mindlin ' , hertz_mindlin
	WRITE(luout, *)'NORMAL STIFFNESS = ', normal_spring_coef
	WRITE(luout, *)'TANGENTIAL STIFFNESS = ', shear_spring_coef
	WRITE(luout, *)'NORMAL DAMPING COEF= ', gamma
	WRITE(luout, *)'FRICTION        = ', friction
	write (*,*) 'ib 1 ', ib(1)
	write (*,*) 'ib 2 ', ib(2)
	write (*,*) 'ib 3 ', ib(3)
	write (*,*) 'ib 4 ', ib(4)
	WRITE(luout, *)' '
	WRITE(luout, *)' TIME step        = ', dt
	WRITE(luout, *)' EXTERNAL LOADING       = ', presy
	WRITE(luout, *)' BODY FORCES       = ', g(1),g(2)
	WRITE(luout, *)' TOP BOTTOM LEFT RIGHT       = ', ytop, ybot,xleft,xright
	!
	end
	
	
	
	subroutine dumpgrains(afrac)
		use mycommons
		
		real*8 a1sum,a2sum,r1sum,r2sum,s1sum,s2sum,afrac
		real*8 amean1,amean2,astd1,astd2,alogr2,ahigr2,alogr,ahigr
		integer n1sum,n2sum, i
		real*8 area1,area2,tmean,tstd

		open(unit=grainsizes_file, file=TRIM(output_directory)//'/grainsizes')
	! set Young modulus for the grains
		do i=1,n
			emod(i)=1.
			write(grainsizes_file,*) gtype(i),2*radius(i),bdgrain(i)
		enddo
		close(grainsizes_file)

	!  get stats about grains
		a1sum=0.
		a2sum=0.
		n1sum=0.
		n2sum=0.
		r1sum=0.
		r2sum=0.
		s1sum=0.
		s2sum=0.
		alogr=higr
		ahigr=logr
		alogr2=higr2
		ahigr2=logr2
		do i=1,n
			if (gtype(i) .eq. 1) then
				r1sum=r1sum+radius(i)
				area1=pi*radius(i)*radius(i)
				if (bdgrain(i) .eq. 0) then
					a1sum=a1sum+area1
				else
					a1sum=a1sum+0.5*area1
				endif
				n1sum=n1sum+1
				alogr=min(alogr,2*radius(i))
				ahigr=max(ahigr,2*radius(i))
			else
				r2sum=r2sum+radius(i)
				area2=pi*radius(i)*radius(i)
				if (bdgrain(i) .eq. 0) then
					a2sum=a2sum+area2
				else
					a2sum=a2sum+0.5*area2
				endif
				n2sum=n2sum+1
				alogr2=min(alogr2,2*radius(i))
				ahigr2=max(ahigr2,2*radius(i))
			endif
		enddo
		amean1=2.*r1sum/n1sum
		if (n2sum .gt. 0.) amean2=2.*r2sum/n2sum

		do i=1,n
			if (gtype(i) .eq. 1) then
				s1sum=s1sum+(2.*radius(i)-amean1)**2
			else
				s2sum=s2sum+(2.*radius(i)-amean2)**2
			endif
		enddo
		astd1=s1sum/n1sum
		astd2=s2sum/max(1,n2sum)
	 
	405 format(a,2(F10.3,2x))
		open(unit=graininfo_file,file=TRIM(output_directory)//'/graininfo')
		write(graininfo_file,*) 'DISTRIBUTION 1: '
		write(graininfo_file,*) '   total grains:', n1sum
		write(graininfo_file,*) '   total area:', a1sum
		write(graininfo_file,*)   '            target values     actual values'
		write(graininfo_file,405) '   mean     ', mean1,amean1
		write(graininfo_file,405) '   std. dev.', std1,astd1
		write(graininfo_file,405) '   smallest ', logr,alogr
		write(graininfo_file,405) '   largest  ', higr,ahigr
		write(graininfo_file,*) ''
		write(graininfo_file,*) 'DISTRIBUTION 2: '
		write(graininfo_file,*) '   total grains:', n2sum
		write(graininfo_file,*) '   total area:', a2sum
		write(graininfo_file,*) '              target values     actual values'
		write(graininfo_file,405) '   mean     ', mean2,amean2
		write(graininfo_file,405) '   std. dev.', std2,astd2
		write(graininfo_file,405) '   smallest ', logr2,alogr2
		write(graininfo_file,405) '   largest  ', higr2,ahigr2
		write(graininfo_file,*) ''

		afrac=a1sum/(a1sum+a2sum)
		write(graininfo_file,*) 'fractional area of Dist. 1:', afrac

		write(graininfo_file,*) ''
		tmean = (amean1*n1sum+amean2*n2sum)/dble(n)
		write(graininfo_file,*) 'mean diameter of entire assemblage:', tmean
		
		do i=1,n
			tstd=tstd+(2.*radius(i)-tmean)**2
		enddo
		tstd=tstd/dble(n)
		write(graininfo_file,*) 'std. dev. of entire assemblage:', tstd
		close(graininfo_file)

	end
	
	subroutine iteration()
		use mycommons
		integer i
		real*8 vm,vmax
		internal_top = maxval(r(2,(nbound(2)+1):n)) + maxval(radius((nbound(2)+1):n))
		!call dump_forces
		if (mod(step, iprintf) .eq. 0) then
			write (*,*) 'dump ', step, iprintf
			

			dump_flux_file = .true.
		else
			dump_flux_file = .false.
		endif
		!  GET NEW POSITION AND VELOCITY AND ROTATION
		CALL MOVEA ()
		!     COMMENT NEXT LINE TO STOP ROTATION
		if (iroll) CALL ROTATEA
		if (ib(2) .eq. 5) then
			dwall = dwall+dt*v(1,nbound(1)+1)+dtsq2*f(1,nbound(1)+1)
			! spring-driven run

			dspring=dspring + dt*uspring
			fspring=kspring*(dspring-dwall)
		else
			fspring=0.
		endif

			


	! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
			  
			  
	!  IF NECESSARY BUILD NEIGHBOR LIST
	!    INCLUDES ALL PARTICLES THAT ARE distint FROM OVERLAP OR LESS.
	!   the maximum number of time steps between rebuilding the list (nbcheck) is
	!   calculated based on the time step and particle velocities. 
	!   however (based on tests with similar sized grains),
	!   once nbcheck gets over 100 it has minimal effect on 
	!   running time, and at 500 the solution begins to change.
			 
		if ( mod( step, nbcheck) .eq. 0 ) then
			
			call links
			!	   call dumplinklist
			call getneighbors()
			
			vmax=0d0
			do i=1,n
				vm=v(1,i)**2+v(2,i)**2
				vmax=max(vmax,vm)
			enddo
			vmax=dsqrt(vmax)
			! at very low speeds nbcheck overflows, 
			! but at very high speeds we want to check more frequently
			nbcheck=ceiling(distint/vmax/dt/5)
			!           print*, "check neighbor list at step", step
			if (nbcheck .le. 0) then
			!	overflow
				nbcheck = 100
			else 			
				! run nbcheck at least every 100 steps
				nbcheck=min(100,nbcheck)
			endif
			
			!           nbcheck=min(1,nbcheck)
		endif

			   
				
		if (ifluid) then
			toprint = 0
			if (mod(step,iprint).eq.0) toprint = 1
			call smooth(toprint)
			do i = 1,TRATIO
				call updatepressure(toprint,i)
			enddo
		endif
		curpresy = presy
		if (IPULSE .eq. 1) then  
			curpresy = presy + amp*dcos(omega*tau)
			if (mod(step,iprint).eq.0) then 
				write(61,*) tau, curpresy ! ERAN - file?
			endif
		endif
		CALL FORCE ( )
		!write (*,*) 'after force'
			   
		   
	!  UPDATE VELOCITIES
	if (fbc .ne.3) then
		
		CALL MOVEB ()
		!write (*,*) 'after moveb'
		if (iroll) CALL ROTATEB
		!write (*,*) 'after rotate b'
		
	endif



	
			   
					 

	! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
	!    IF SPRING-DRIVEN, UPDATE WALL POSITION
	!          if (ib(2) .eq. 5) then
	!	   dwall = dwall+dt*v(1,nbound(1)+1)+dtsq2*f(1,nbound(1)+1)
	!	  endif
	! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

		tau = tau + dt

	  

	!  ##### OUTPUT ##### OUTPUT ##### OUTPUT ##### OUTPUT ##### OUTPUT 
	! ****************************************************************
	! &&&  FREQUENT OUTPUT &&&&&&&  FREQUENT OUTPUT &&&&
		IF ( MOD( step, iprintf ) .EQ. 0 ) THEN
			write(*,*) 'step, iprintf ' , step, iprintf
			

	!  * coordination number
	!		coordave=0.
	!		coordave=coordave+sum(coordn(indint:n))
	!		write(coordination_file,fmt='(e14.7,1x,e14.7)') tau, coordave

			!call dumplinklist
			!call dump_forces
	! height between the two walls
			write(dilatation_file,*) tau,(ytop-ybot)


	!  * porosity
			if (phigoal .gt. 0.) then
				topg=r(2,1)+radius(1)
				do i=nbound(2)+1,n
					topg=max(topg,r(2,i)+radius(i))
				enddo
				area = (xright-xleft)*(topg-ybot)
			else
				area = (xright-xleft)*(ytop-ybot)
			endif

			phi = 1.-grainarea/area
			if (ib(2) .lt. 0) then
				phi = 1.-grainarea/((xright-xleft)*(internal_top-ybot))
			endif
			write(porosity_file,fmt='(e14.7,1x,e14.7,1x,e14.7)') tau,phi

			 
	!  * stresses
		!	fshear = fxbot - fxtop
		!	fnormal = fytop - fybot
			
	! these are the average forces on each grain. 
	!  convert to stresses
			sxt=fxtop*(nbound(2)-nbound(1))/(xright-xleft)              
			syt=fytop*(nbound(2)-nbound(1))/(xright-xleft)              
			sxb=fxbot*nbound(1)/(xright-xleft)              
			syb=fybot*nbound(1)/(xright-xleft)  
			if (ifluid) then
				write(stresses_file,716) tau,sxt,syt,sxb,syb,fgrain
			else
				write(stresses_file,716) tau,sxt,syt,sxb,syb
			endif
		716 format(6(e14.7,1x))              

					
		!  waves at walls
			if (ipulse .eq. 1 .and. ib(2) .eq. 5) then
				write(39,714) tau, fspring,fxbot,fybot ! ERAN file?
			endif
			if (ib(1) .eq. 6 .or. ib(1) .eq. 3) then
				write(49,714) tau, r(1,1) ! ERAN file?
			endif


	!  spring-driven wall
	!  if with fluid top friction also contains the ratio shear stress/effective normal
	!  when fgrain contains the effective normal stress
			if (ib(2) .eq. 5) then
				if (ifluid) then
					write(topfriction_file,714) tau, fspring/presy, fspring/fgrain, v(1,nbound(1)+1)/fb(2)
				else
					write(topfriction_file,714) tau, fspring/presy, v(1,nbound(1)+1)/fb(2)
				endif
				write(spring_file,714) tau, f(1,nbound(1)+1),dwall,dspring
			else
				if (ifluid) then
					write(topfriction_file,714) tau, -sxt/presy, -sxt/fgrain
				else
					write(topfriction_file,714) tau, -sxt/presy
				endif
			endif


	!  * kinetic energy and momentum
			ek = 0.0
			ekr = 0.0
			mvx=0.
			mvy=0.
			angmom=0.
			do i=indint,n
				
				ek = ek + mass(i)*(v(1,i)*v(1,i) + v(2,i)*v(2,i)) 
				ekr = ekr + 0.5*mass(i)*radius(i)*radius(i)*w(i)*w(i)
				mvx = mvx + mass(i)*v(1,i) 
				mvy = mvy + mass(i)*v(2,i) 
				angmom = angmom + 0.5*mass(i)*radius(i)*radius(i)*w(i)
			enddo
			ek = 0.5  * ek
			ekr = 0.5  * ekr 
			
	!  boundary dissipation
			write(momentum_file,714) tau, mvx,mvy,angmom
			write(energy_file,716) tau, ek, ekr, ev ! Eran ev?
	714 	format(5(e14.7,1x))
		ENDIF
	! &&&  FREQUENT OUTPUT &&&&&&&  FREQUENT OUTPUT &&&&

			   
				   
	! &&&  INFREQUENT OUTPUT &&&&&&&  INFREQUENT OUTPUT &&&&
	!  &&&&&&& SNAPSHOTS &&&&&&& SNAPSHOTS &&&&&&& SNAPSHOTS &&&&&&& SNAPSHOTS
		if ( mod( step, iprint ) .eq. 0 ) then
			print*, 'writing postscript step = ', step
			call pscriptplot(0)
			call pscriptplot(99999)
		endif
		if (output_burst_step * output_burst_id_step .gt. output_burst_id .and. mod(step, output_burst_interval) .eq. 0) then
			call pscriptplot(output_burst_id)
			output_burst_id = output_burst_id + output_burst_id_step
			if (output_burst_step * output_burst_id_step .le. output_burst_id) then
				output_burst_step = 0;
				WRITE(*,*) '******   output burst finished ', output_burst_id, phase_step,' phase_step   **** '
			endif
		endif
		if ( mod( step, gen_restart ) .eq. 0 ) then
			write(*,*) 'gen restart ', step, gen_restart
			

			call write_restart()           
		endif
	! &&&  INFREQUENT OUTPUT &&&&&&&  INFREQUENT OUTPUT &&&&
	!  ##### OUTPUT ##### OUTPUT ##### OUTPUT ##### OUTPUT ##### OUTPUT 

		if (mod(step, check_condition_every) .eq. 0) then
			write(*,*) 'check stop ', step, check_condition_every
			

			call check_stop_condition()
			INQUIRE(FILE=TRIM(output_directory)//'/skip', EXIST=skip_phase)
			if (output_burst_step .eq. 0) then
				INQUIRE(FILE=TRIM(output_directory)//'/output', EXIST=output_burst)
			endif
		endif
		
		
	end

	
subroutine check_stop_condition()
	use mycommons
	integer i
	real*8 curr_ave_flux
	! average Ek stop condition
	ek = 0.0
	grains_per_layer = 0
	velocities = 0d0
	do i=indint,n	
		ek = ek + mass(i)*(v(1,i)*v(1,i) + v(2,i)*v(2,i)) 
		velocities(grain_output_layer(i)) = velocities(grain_output_layer(i)) + v(1,i)
		grains_per_layer(grain_output_layer(i)) = grains_per_layer(grain_output_layer(i)) + 1
	enddo
	!write(*,*) 'v ', velocities(grain_output_layer(1000)), ' gpl ', grains_per_layer(grain_output_layer(1000)), 
	ek = 0.5  * ek
	average_ek = ek / dble(n)
	
	write (*,*) ' phase_step ', phase_step, ', average Ek: ',average_ek,' target ', target_ek
	if (average_ek < target_ek) then
		stop_condition_reached = stop_condition_reached + 1
	else
		stop_condition_reached = 0
	endif
	
	curr_ave_flux = sum(flux) / dble(check_condition_every)
	! flux change stop condition
	if (dabs(curr_ave_flux - average_flux_for_condition) .lt. target_flux_change) then
		stop_condition_reached = stop_condition_reached + 1
		average_flux_for_condition = (average_flux_for_condition + curr_ave_flux)/stop_condition_reached
		write (*,*) 'met step ', step, ', average flux: ',average_flux_for_condition,' curr_ave_flux ', curr_ave_flux
	else
		!stop_condition_reached = 0
		write (*,*) 'not met step ', step, ', average flux: ',average_flux_for_condition,' curr_ave_flux ', curr_ave_flux
		average_flux_for_condition = curr_ave_flux
	
	endif
	write(grains_per_layer_file,'(I12)', advance='no') step
	write(flux_file,'(I12)', advance='no') step
	write(vol_flux_file,'(I12)', advance='no') step
	write(velocity_pro_file,'(I12)', advance='no') step
	!write(flux_file,'(",",E21.15)', advance='no') curr_angle(1)
	!write(flux_file,'(",",E21.15)', advance='no') curr_angle(2)
	!write(flux_file,'(",",E21.15)', advance='no') tau
	!write(flux_file,'(",",E21.15)', advance='no') phi
	write(grains_per_layer_file,'(",",E21.15)', advance='no') internal_top - ybot
	write(flux_file,'(",",E21.15)', advance='no') average_ek
	write(vol_flux_file,'(",",E21.15)', advance='no') average_ek
	write(velocity_pro_file,'(",",E21.15)', advance='no') maxval(compression)
	
	do i=1,number_of_layers !int(internal_top - ybot) / grains_per_layer
		write(grains_per_layer_file,'(",",I10)', advance='no') grains_per_column_per_layer(i)
		write(grains_per_layer_file,'(",",I5)', advance='no') grains_per_layer(i)
		write(flux_file,'(",",E21.15)', advance='no') flux(i)
		write(vol_flux_file,'(",",E21.15)', advance='no') vol_flux(i)
		write(velocity_pro_file,'(",",E21.15)', advance='no') velocities(i)/grains_per_layer(i)
	enddo
	write(grains_per_layer_file, '(",",I10)') sum(grains_per_layer)
	write(flux_file,'(",",E21.15)') sum(flux)
	write(vol_flux_file,'(",",E21.15)') sum(vol_flux)
	write(velocity_pro_file, '(",",E21.15)') sum(velocities)/sum(grains_per_layer)
	
	grains_per_column_per_layer = 0
	flux = 0
	vol_flux = 0
	
	write (*,*) 'step ', step, 'stop_condition_reached ', stop_condition_reached
end

