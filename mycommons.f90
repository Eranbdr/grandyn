MODULE mycommons

	integer maxn,maxcell,maxk,maxd, tratio, max_flux_layers
	parameter (MAXN=40000)
	parameter (MAXK=30*MAXN)
	parameter (MAXCELL = 9025)
	parameter (MAXD = 55)
	parameter (TRATIO = 1)
	parameter (max_flux_layers = 100)
	
! CURRENT TALLY OF ARRAYS AND SIZES
!  SIZE MAXN: 29 real*8 + 5 integer
!  SIZE MAXK=MAXN*25: 4 real*8 + 3 integer
!  SIZE MAXCELL: 2
!  TOTAL: 1352*MAXN bytes

!  ATOM ARRAYS
	REAL*8      r(2,MAXN)
	REAL*8      v(2,MAXN)
	REAL*8      f(2,MAXN)
	REAL*8		accel(2,MAXN)
	REAL*8 	    orientation(MAXN),w(MAXN),tq(MAXN)
	REAL*8	    radius(MAXN),radinv2(MAXN), mass(MAXN), velocities(max_flux_layers)
	real*8		flux(max_flux_layers), ave_flux(max_flux_layers), vol_flux(max_flux_layers), ave_vol_flux(max_flux_layers)
	!REAL*8      coordn(MAXN)
	real*8 	    emod(maxn), amp_rand, grain_rand(2,MAXN), phase_start_height
	INTEGER     n,indint, no_modification_from_last_contact
	INTEGER     color(MAXN), last_contact_step(MAXN), grains_per_layer(max_flux_layers), grains_per_column_per_layer(max_flux_layers)
	integer	    gtype(maxn),bdgrain(maxn), grain_output_layer(maxn)
	logical 	first_active_layer(maxn), random_flag
!	**  N		NO. OF ATOMS
!	**  INDINT	INDEX OF FIRST NON-BOUNDARY ATOM
!	**  NINTN	NO. OF ATOMS (NOT IN BOUNDARY)
!	**  RX,RY	ATOMIC POSITIONS
!	**  VX,VY	ATOMIC VELOCITIES
!	**  FX,FY	ATOMIC FORCES
!	**  orientation		ORIENTATION OF THE ATOM (DEGREES CCW FROM EAST)
!	**  W		ANGULAR VELOCITY OF ATOM
!	**  TQ		TORQUE ON ATOM
!	**  RADIUS 	ATOMIC RADIUS
!	**  RADINV2	INVERSE SQUARE OF RADIUS
!       **  ONBD	TRUE IF ATOM IS TOUCHING A BOUNDARY ATOM
!       **  SX,SY,SXY   NORMAL AND SHEAR STRESSES IN X-Y SYSTEM
!       **  S1,S2       PRINCIPLE STRESSES
!       **  THETA       ORIENTATION OF PRINCIPLE COMPRESSIVE STRESS
!       **  COORDN      COORDINATION NUMBER OF GRAIN
!       **  DEDGE	DISTANCE TO PERIODIC BOUNDARY
!	**  COLOR	COLOR OF ATOM IN LEFT-HAND PICTURE
!       **  GTYPE	INDEX OF THE DISTRIBUTION (GRAIN TYPE)
!       **  EMOD        YOUNG'S MODULUS OF THE GRAIN

!  BOUNDARY CONDITION ARRAYS
	REAL*8	    fb(4)
	INTEGER	    nbound(4), ib(4)
! 	** BOUNDARIES ARE LISTED IN THE ORDER: BOTTOM,TOP,LEFT,RIGHT
!	**  NBOUND(4)	NO. OF ATOMS IN EACH BOUNDARY
!	**  IB(4)	TYPE OF EACH BOUNDARY

!   LINK-LIST ARRAYS
	INTEGER	    list(2*MAXN),head(MAXCELL),map(4*MAXCELL)
	INTEGER	    mx,my
!	**  HEAD	HEAD OF CHAIN FOR EACH CELL
!	**  LIST	LINKED LIST OF GRAINS
! 	**  MAP		LIST OF NEIGHBORING CELLS (ONLY UPPER-RIGHT HALF)
!	**  MX 		TOTAL NUMBER OF CELLS IN X-DIRECTION
!	**  MY 		TOTAL NUMBER OF CELLS IN Y-DIRECTION

! CONTACT ARRAYS
	INTEGER contacti(maxk),contactj(maxk),contactknt
	REAL*8  contfn(maxk),contft(maxk), compression(maxk)

!  SCALAR PARAMETERS
	REAL*8 boxx,boxy,sigb,xleft,xright,ytop,ybot, normal_spring_coef,gamma,tau, &
		 shear_spring_coef,gamshear,friction,dt
	INTEGER natoms
!	**  SIG		AVERAGE ATOM DIAMETER (DEFINED AS 1.)
!       ** **INPUT PARAMETERS**
!	**  BOXX, BOXY  SIZE OF BOX IN X AND Y (RELATIVE TO AVERAGE
!	**	ATOM DIAMETER
!	**  SIGB	AVERAGE BOUNDARY ATOM DIAMETER 
!	**  DMOL	SETS DISTANCE FROM MOLECULE CENTER TO CENTER OF ATOMS
!	**		1.=ATOMS ARE TANGENT, < 1. ATOMS OVERLAP
!	**  NATOMS	NUMBER OF ATOMS PER MOLECULE
!	**  normal_spring_coef	SPRING CONSTANT FOR INTRAMOLECULAR FORCES
!	**  shear_spring_coef     SPRING CONSTANT FOR LEAF SPRING
!       **  GAMMA 	DAMPING PARAMETER FOR NORMAL VELOCITY
!       **  GAMMAMOL    DAMPING PARAMETER FOR INTERACTIONS BETWEEN ATOMS
!       **              IN A MOLECULE 
!       **  GAMSHEAR    DAMPING PARAMETER FOR TANGENTIAL VELOCITY
!	**  TAU		CURRENT TIME
!	**  TAUADD	TIME AT WHICH NEXT PARTICLE ADDED (FOR GROWING
!			SANDPILE PROBLEM)
!
!	** **CALCULATED PARAMETERS**
!	** XLEFT, XRIGHT, YTOP,YBOT
! 	** 		CURRENT COORDINATES FOR BOX EDGES


! WITH THE ADDITION OF FLUID
	real*8 		dxc,dyc
	real*8 		phimat(MAXD+1,MAXD+1)
 	real*8 		perm(MAXD+1, MAXD+1)
	real*8 		sitevx(MAXD+1, MAXD+1)
	real*8		sitevy(MAXD+1,MAXD+1)
	real*8 		divvel(MAXD+1,MAXD+1)
	real*8		press(MAXD+1,MAXD+1)
	real*8		gradpx(MAXD+1,MAXD+1)
	real*8		gradpy(MAXD+1,MAXD+1)
	real*8		solfra(MAXD+1,MAXD+1)
	real*8		numpar(MAXD+1,MAXD+1)
	real*8 		vc(MAXD+1,MAXD+1)
	real*8 		siterd2(MAXD+1, MAXD+1)
	real*8		pbot,ptop,fluxbot
 	real*8 		famp, fomega,pout
	real*8		distint
	
	data	distint /0.1d0/
	!data	nbcheck /100/
	
!   originally from main

	integer*8     nstep, step, iprint,toprint, gen_restart
	real*8        G(2)
	real*8        FYTOP, FYBOT
	real*8        mvx, mvy
	
 
	real*8 ev,ek,ekr, average_ek, target_ek
	real*8 area,phi,fxbot,fxtop
	real*8 boxfy 
	
	character*12 fnn
	integer stop_condition_reached
	
	real*8 acceleration_amplitude(2),acceleration_decay_limit(2),angle_change_rate(2), curr_angle(2)
	integer randomize_acceleration, acceleration_function(2), acceleration_decay_type(2), check_condition_every
	logical  wait_for_condition
	
	real*8 intensity(2)
	integer wave_length(2), currphase, phase_step, number_of_layers, output_burst_step, output_burst_interval, output_burst_id, output_burst_id_step
	
	real*8 average_flux_for_condition, curr_flux, target_flux_change, internal_top
	
	
	
	logical dump_flux_file, skip_phase, output_burst
	
	
	integer LUOUT, coordination_file, porosity_file, stresses_file, momentum_file, energy_file, spring_file, waveinput_file, &
			strain_file, topfriction_file, dilatation_file, periodicf_file, randno_file, grainsizes_file, graininfo_file, &
			link_list_dump_file, contact_info_file, phi_perm_file, areas_file, velocity_file, num_velocity_file, velocity_pro_file, &
			pressure_file, pressure_grad_file, psunit, restart_file, debug_file, forces_file, flux_file, vol_flux_file, layer_change_file, &
			contact_info_zero_file, contact_info_mid_file, contact_info_high_file, grains_per_layer_file
	
	character*17 :: default_config_file = 'config/config.xml'
	character*26 :: config_template_file = 'config/config_template.xml'
	character*64 :: config_file 
	data randno_file /2/
	data restart_file /3/
	data luout  /6/
	data waveinput_file /7/
	data phi_perm_file /9/
	data areas_file /10/
	data velocity_file /11/
	data velocity_pro_file /12/
	data link_list_dump_file /13/
	data coordination_file /14/
	data porosity_file /15/
	data stresses_file /16/
	data momentum_file /17/
	data energy_file /18/
	data grainsizes_file /19/
	data graininfo_file /20/
	data flux_file /21/
	data vol_flux_file /22/
	data grains_per_layer_file /23/
	data strain_file /29/
	data contact_info_file /31/
	data spring_file /33/
	data dilatation_file /34/
	data contact_info_zero_file /35/
	data contact_info_mid_file /36/
	data contact_info_high_file /37/
	data psunit /45/
	data debug_file /46/
	data forces_file /47/
	data topfriction_file /57/   
	data pressure_file /58/
	data pressure_grad_file /59/
	data layer_change_file /66/
	data num_velocity_file /78/
	data periodicf_file /81/


	real*8 boxvol,alap,angmom,presy,presx,phigoal !, coordave,fnormal,fshear
	logical irstart, iroll,iztau,igif,izmom,izslip, irandfile, hertz_mindlin
	character*256 restart_file_name
	character*256 output_directory
	character*128 file_postfix
	real*8 topmass
	

	integer randum
	
	real*8 pi
	data pi /3.1415926535897932d0/
	real*8 degree_to_radian
	
    
	real*8 dspring,dwall,kspring,fspring,uspring
	integer itmass

	data itmass /0/
	
	real*8 frac1,mean1,std1,mean2,std2,logr,higr,logr2,higr2,e1,e2
	real*8 presyfinal
	integer nbcheck
	real*8 bg,dtsave
	real*8 tau00,wavelength,period,omega,twave,amp,presy0,cs
	real*8 es,ctime,rybot,rytop,gwide
	integer ipulse, stop_condition_limit
	real*8 wmax,topg
	real*8 afrac1,ftsum,prop1,prop2,diffsum1,diffsum2,diffrac1,diffrac2
	real*8 sxt,syt,sxb,syb
    logical terminate_on_condition_not_met
	integer iprintf
	real*8 fstress, dtsq2, dt2
	integer fbc
	integer mxmax,mymax
	real*8 curpresy, grainarea
	real*8 fgrain, fgrain0
	
 
	logical config_status, ifluid
        
	integer STATUS
!   ***********  util functions    *********************
!	contains
	!subroutine use_a_flexible_string( length )
		!
		!character(len=length) :: out_str
		!print*, 'in_str ', output_directory, '111'
		!out_str = trim(output_directory)
		!print*, 'out_str ', out_str, '111'
!
!
	!end subroutine
	
	
	!function get_output_directory( length , is_init )
	!	
	!	character(len=length) :: out_str
	!	logical is_init
	!	save out_str
	!	if (is_init) out_str = trim(output_directory)
	!	get_output_directory = out_str
	!	return

	!end function
	
!	subroutine open_file( output_length , file_handle, is_init )
		!
		!integer output_length
		!character(len=output_length) :: trimmed_output
		!logical is_init
		!integer file_handle
		!!save trimmed_output
		!print*, 'trimmed_output : ', trimmed_output, 'more'
		!if (is_init) trimmed_output = trim(output_directory)
		!!print*, 'trimmed_output//files(file_handle) : ', trimmed_output//files(file_handle), 'more'
		!!open(unit=file_handle,file=trimmed_output//files(file_handle))
		!return
!
	!end subroutine
	contains
	real*8 function acceleration_modification(direction, i)
	!use mycommons
	integer i, direction, grain_layer, contact
	real*8 relative_location, currand
	real*8 actual_modification_fraction
	
	acceleration_modification = 0.
	
	if ((no_modification_from_last_contact .ne. 0) .and. (step - last_contact_step(i) .gt. no_modification_from_last_contact)) then
		!write(*,*) step ,' grain ', i, ' last contact ', last_contact_step(i)
		return
	endif
	
!	acceleration decay type -
!	0 - no decay  
!	1 - depending on distance from top or left wall
!		decaying linearly so that force modification is 0 when 
!		r(direction,i)/rmax(direction)) <= acceleration_decay(direction)
!	2 - as above, from bottom of right wall
!	3 - number of layers from top or left that feel the modification, linear decay
!	4 - as above, from bottom of right wall
!	5 - only grains in contact with bottom wall
!	6 - only grains in specific layer feel the modification (each layer has the width of the largest grain radius)
!	7 - linear decay from top + random choice 
!	8 - random, no decay
!	9 - quadratic bottom - not relative!
	
	SELECT CASE (acceleration_decay_type(direction))  
		CASE (0)
			actual_modification_fraction = 1.
		CASE (1)
			if (direction .eq. 1) then ! horizontal
				relative_location = (r(direction,i)- xleft)/boxx
			else ! vertical
				relative_location = (r(direction,i) - ybot)/(internal_top - ybot)
			endif
			
			if (relative_location .gt. acceleration_decay_limit(direction)) then ! grain feels modification
				actual_modification_fraction = 1-(1-relative_location)/(1-acceleration_decay_limit(direction))
			else ! grain does not feel modification
				actual_modification_fraction = 0.
			endif
		CASE (2) 
			if (direction .eq. 1) then ! horizontal
				relative_location = (r(direction,i)- xleft)/boxx
			else ! vertical
				relative_location = (r(direction,i) - ybot)/(internal_top - ybot)
			endif
			
			if (relative_location .lt. acceleration_decay_limit(direction)) then ! grain feels modification
				actual_modification_fraction = 1-relative_location/acceleration_decay_limit(direction)
			else ! grain does not feel modification
				actual_modification_fraction = 0.
			endif
		CASE(3)
			
			!grains bucketed to bg/2 depth (width) rows (columns)
			if (direction .eq. 1) then ! horizontal
				grain_layer = (r(direction,i)- xleft - radius(i))/( bg / 2)
			else ! vertical
				grain_layer = (r(direction,i)- ybot - radius(i))/( bg / 2)
			endif
			
			if (grain_layer .gt. acceleration_decay_limit(direction)) then ! grain feels modification
				actual_modification_fraction = 1-(1-relative_location)/(1-acceleration_decay_limit(direction))
			else ! grain does not feel modification
				actual_modification_fraction = 0.
			endif
		CASE(4)
			
			!grains bucketed to bg/2 depth (width) rows (columns)
			if (direction .eq. 1) then ! horizontal
				grain_layer = (r(direction,i)- xleft - radius(i))/( bg / 2)
			else ! vertical
				grain_layer = (r(direction,i)- ybot - radius(i))/( bg / 2)
			endif
			
			if (grain_layer .lt. acceleration_decay_limit(direction)) then ! grain feels modification
				actual_modification_fraction = 1-relative_location/acceleration_decay_limit(direction)
			else ! grain does not feel modification
				actual_modification_fraction = 0.
			endif
		CASE(5)
			if (.not. first_active_layer(i)) then
				actual_modification_fraction = 0.
			else
				!write(*,*) 'acceleration ', acceleration_modification , ' on grain ', i
				!if (mod(step, 100000) .eq. 0) write(*,*) 'decay type 5, grain ', i
				actual_modification_fraction = 1.
			endif
		CASE(6)
			
			!grains bucketed to bg/2 depth (width) 
			grain_layer = (r(2,i)- ybot - radius(i))/( bg / 2)
			
			
			if (grain_layer .eq. acceleration_decay_limit(2)) then ! grain feels modification
				actual_modification_fraction = 1.
				if (mod(step, 100000) .eq. 0) write(*,*) 'decay type 6, grain ', i
			else ! grain does not feel modification
				actual_modification_fraction = 0.
			endif		
		CASE (7)
			
			relative_location = (r(2,i) - ybot)/(internal_top - ybot)
			
			!if (rand(0) .le. (1d0 / wave_length(direction))) then			
				if (relative_location .gt. acceleration_decay_limit(2)) then ! grain feels modification
					actual_modification_fraction = 1-(1-relative_location)/(1-acceleration_decay_limit(direction))
					if (rand(0) .gt. actual_modification_fraction / wave_length(direction)) then
						actual_modification_fraction = 0.
						!write(*,*) 'NO ******** acceleration fraction ', actual_modification_fraction , ' on grain ', i, 'direction ', direction
					else
						!write(*,*) 'acceleration fraction ', actual_modification_fraction , ' on grain ', i, 'direction ', direction
					endif
				else ! grain does not feel modification
					actual_modification_fraction = 0.
				endif
			!endif	
		
		CASE (8)
			if (rand(0) .le. (1d0 / wave_length(direction))) then	
				actual_modification_fraction = 1.
			endif
		CASE(9)
			actual_modification_fraction = 1d0/(r(2,i) - ybot)**2
		CASE DEFAULT
			write(*,*) 'unknown decay type ', acceleration_decay_type(direction)
			stop
	END SELECT
	
	
	!	acceleration_function types
!	0 - none
!	1 - perturbation sin (use given amplitude)
!	2 - perturbation cos (use given amplitude)
!	3 - sin, only first half of the wave (0-180 degrees)
!	4 - cos, only first half of the wave (0-180 degrees)
!	5 - random -amp to amplitude every wave_length steps
!	6 - random negative amplitude to amplitude every wave_length steps, 
!		negative values correlate to 0 acceleration
!	7 - modification randomly set horizontally -amp to amp, vertically -0.25 amp to amp
!	8 - vertical white noise (tilt adjusted in main)
!	9 - perturbation sin negative values correlate to 0 acceleration
!	10 - perturbation cos negative values correlate to 0 acceleration
!	11 - modification randomly set each step -amp to amp
!	12 - 9 + 5/2
!	13 - 12 full sin + 45 deg phase
!	14 - 12 full sin
	if (actual_modification_fraction .gt. 0) then

		SELECT CASE (acceleration_function(direction))  
			CASE (0)
				acceleration_modification = 0.
				return
			CASE (1)
				acceleration_modification = acceleration_amplitude(direction)*dsin(curr_angle(direction)*degree_to_radian)
			CASE (2)
				acceleration_modification = acceleration_amplitude(direction)*dcos(curr_angle(direction)*degree_to_radian)
			CASE (3)
				acceleration_modification = acceleration_amplitude(direction)*abs(dsin(curr_angle(direction)*degree_to_radian))
			CASE (4)
				acceleration_modification = acceleration_amplitude(direction)*abs(dcos(curr_angle(direction)*degree_to_radian))
			CASE(5)
				if (mod(int(step * sqrt(abs(grain_rand(direction,i)))), wave_length(direction)) .eq. 0) then
					
					grain_rand(direction,i) = 1d0 - rand(0) * 2d0
					!if (mod(i, 500) .eq. 0) then
					!	write(*,*) step, i, direction, grain_rand(direction,i)
					!endif
				endif
				acceleration_modification = acceleration_amplitude(direction) * dble(grain_rand(direction,i))
				!write(*,*) 'acceleration_amplitude ', acceleration_amplitude(direction), dble(amp_rand), acceleration_modification
			CASE(6)
				if (mod(step, wave_length(direction)) .eq. 0) then
					if (random_flag) then 
						amp_rand = 1d0 - 2d0 * rand(0) 
						!write(*,*) 'rand ', amp_rand
						random_flag = .false.
						
					endif
				else
					random_flag = .true.
				endif
				acceleration_modification = max(0d0, acceleration_amplitude(direction) * dble(amp_rand))
				!write(*,*) 'acceleration_amplitude ', acceleration_amplitude(direction), dble(amp_rand), acceleration_modification
			
			CASE(7)
				if (direction .eq. 1) then 	
					acceleration_modification = acceleration_amplitude(direction) * (1d0 - 2d0 * dble(rand(0)))
				else	
					acceleration_modification = acceleration_amplitude(direction) * (-0.25 + 1.25*dble(rand(0)))
				endif
				!write(*,*) 'actual_modification_fraction, modification ', actual_modification_fraction, acceleration_modification
			CASE(8)
				if (mod(step, wave_length(direction))/4 .eq. 0) then
					if (random_flag) then 
						amp_rand = 0.75 + rand(0) * 0.5
						!write(*,*) 'rand ', amp_rand
						random_flag = .false.
						
					endif
				else
					random_flag = .true.
				endif
				!acceleration_modification = acceleration_amplitude(direction)*abs(dsin(curr_angle(direction)*degree_to_radian))* dble(amp_rand)
				acceleration_modification = acceleration_amplitude(direction)*max(0d0, dsin(curr_angle(direction)*degree_to_radian))* dble(amp_rand)
			CASE (9)
				acceleration_modification = acceleration_amplitude(direction)*max(0d0, dsin(curr_angle(direction)*degree_to_radian))
			CASE (10)
				acceleration_modification = acceleration_amplitude(direction)*max(0d0, dcos(curr_angle(direction)*degree_to_radian))
			CASE(11)	
				acceleration_modification = acceleration_amplitude(direction) * (1d0 - 2d0 * dble(rand(0)))
			CASE (12)
				 
				if (mod(int(step * sqrt(abs(grain_rand(direction,i)))),wave_length(direction)) .eq. 0) then
					
					grain_rand(direction,i) = 1d0 - rand(0) * 2d0
					!if (mod(i, 500) .eq. 0) then
					!	write(*,*) step, i, direction, grain_rand(direction,i)
					!endif
				endif
				acceleration_modification = acceleration_amplitude(direction)*max(0d0, dsin(curr_angle(direction)*degree_to_radian)) + acceleration_amplitude(direction) * dble(grain_rand(direction,i))/2
				
			CASE (13)
			 
				if (mod(int(step * sqrt(abs(grain_rand(direction,i)))),wave_length(direction)) .eq. 0) then
					
					grain_rand(direction,i) = 1d0 - rand(0) * 2d0
					!if (mod(i, 500) .eq. 0) then
					!	write(*,*) step, i, direction, grain_rand(direction,i)
					!endif
				endif
				acceleration_modification = acceleration_amplitude(direction)*dsin((45 + curr_angle(direction))*degree_to_radian) + acceleration_amplitude(direction) * dble(grain_rand(direction,i))/2
			 
			CASE (14)
				if (mod(int(step * sqrt(abs(grain_rand(direction,i)))),wave_length(direction)) .eq. 0) then
					
					grain_rand(direction,i) = 1d0 - rand(0) * 2d0
					!if (mod(i, 500) .eq. 0) then
					!	write(*,*) step, i, direction, grain_rand(direction,i)
					!endif
				endif
				acceleration_modification = acceleration_amplitude(direction)*dsin((curr_angle(direction))*degree_to_radian) + acceleration_amplitude(direction) * dble(grain_rand(direction,i))/2
			
			CASE DEFAULT
				acceleration_modification = 0.
				write(*,*) 'unknown acceleration function ', acceleration_function(direction)
				stop
		END SELECT
		acceleration_modification = acceleration_modification * actual_modification_fraction
	else
		acceleration_modification = 0.
	endif
		
end function
END MODULE
