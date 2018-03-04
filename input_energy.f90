PROGRAM find_accel
	IMPLICIT NONE
	integer grains, steps, numps
	parameter (steps = 1000000)
	parameter (grains = 100)
	parameter (numps = 4)
	
	real, dimension(7,numps,grains) :: grain_rand, act_accel, abs_accel, act_total_accel, &
										abs_total_accel, grain_vel, grain_vel_sq
	integer :: wls(7), amps(numps), sizes(4), wl, amp, grain, changes_file, totals, step, averages, avechange
	integer :: changes(7,numps,grains), step_changes(7,numps,grains), acc_type, run
	real :: absvrage, average, avel, avel_sq, step_factor!,  cont_changes(7,numps,grains, 600,2)
	character (len=31) :: file_name
	wls= (/10000, 12000, 15000 ,5000, 8000 ,18000, 20000 /)
	!wls= (/5000, 8000 ,10000, 12000, 15000 ,18000, 20000 /)
	amps = (/11, 19, 25, 31 /)
	acc_type = 12
	step_factor = 0.05
	run = 1
	!amps = (/11 /)
	sizes = (/ 7, 3 , 3, 2 /)
	grain_rand = 0
	act_accel = 0
	abs_accel = 0
	act_total_accel = 0
	abs_total_accel = 0
	changes = 0
	grain_vel_sq = 0
	grain_vel = 0
	step_changes = 0
	write(file_name,fmt='("_type_",I2,"_factor_",e4.2,"_run_",i1,".csv")') acc_type, step_factor, run
	open(unit=changes_file,file='changes'//file_name)
	write(changes_file, '(a)') 'grain wl amp, steps, wl, amp, act_accel '
	
	do amp=1, numps
		do wl=1, sizes(amp)
			do grain=1, grains
				write(changes_file, '(",",I4, "_", I6,"_",I3)' , advance='no') grain, wls(wl), amps(amp)
			enddo
		enddo
	enddo
	write(changes_file, '(a)') ''
	do step=1, steps
	!	if (mod(step, 1000) .eq. 0) then
	!		write(changes_file, '(I10)' , advance='no') step
	!	endif
		do amp=1, numps
			do wl=1, sizes(amp)
				do grain=1, grains
				!decay	!if (mod(int(step * sqrt(abs(grain_rand(wl,amp,grain)))), wls(wl)) .eq. 0) then
				!func	!grain_rand(wl,amp,grain) = 1d0 - rand(0) * 2d0
				!endif
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
	
		SELECT CASE (acc_type)
			CASE (0)
				act_accel(wl,amp,grain) = 0.
				return
			CASE (1)
				act_accel(wl,amp,grain) = amps(amp)*dsin(dble(step) / wls(wl))
			CASE (2)
				act_accel(wl,amp,grain) = amps(amp)*dcos(dble(step) / wls(wl))
			CASE (3)
				act_accel(wl,amp,grain) = amps(amp)*abs(dsin(dble(step) / wls(wl)))
			CASE (4)
				act_accel(wl,amp,grain) = amps(amp)*abs(dcos(dble(step) / wls(wl)))
			CASE(5)
				if (mod(int(step * sqrt(abs(grain_rand(wl,amp,grain)))), wls(wl)) .eq. 0) then
					grain_rand(wl,amp,grain) = 1d0 - rand(0) * 2d0
				endif
				act_accel(wl,amp,grain) = amps(amp) * grain_rand(wl,amp,grain)
				!write(*,*) 'acceleration_amplitude ', amps(amp), dble(amp_rand), act_accel(wl,amp,grain)
			CASE(6)
				if (mod(step, wls(wl)) .eq. 0) then
				!	if (random_flag) then 
				!		amp_rand = 1d0 - 2d0 * rand(0) 
						!write(*,*) 'rand ', amp_rand
				!		random_flag = .false.
						
				!	endif
				!else
				!	random_flag = .true.
				!endif
					act_accel(wl,amp,grain) = max(0d0, amps(amp) * (1d0 - rand(0) * 2d0))
				!write(*,*) 'acceleration_amplitude ', amps(amp), dble(amp_rand), act_accel(wl,amp,grain)
				endif
			CASE(7)
				!if (direction .eq. 1) then 	
					act_accel(wl,amp,grain) = amps(amp) * (1d0 - 2d0 * dble(rand(0)))
				!else	
					act_accel(wl,amp,grain) = amps(amp) * (-0.25 + 1.25*dble(rand(0)))
				!endif
				!write(*,*) 'actual_modification_fraction, modification ', actual_modification_fraction, act_accel(wl,amp,grain)
			CASE(8)
				if (mod(step, wls(wl)/4) .eq. 0) then
					!if (random_flag) then 
					grain_rand(wl,amp,grain) = 0.75 + rand(0) * 0.5
						!write(*,*) 'rand ', amp_rand
					!	random_flag = .false.
						
					endif
				!else
					!random_flag = .true.
					
				!endif
				!act_accel(wl,amp,grain) = amps(amp)*abs(dsin(step / wls(wl)))* dble(amp_rand)
				act_accel(wl,amp,grain) = amps(amp)*max(0d0, dsin(dble(step) / wls(wl)))* grain_rand(wl,amp,grain)
			CASE (9)
				act_accel(wl,amp,grain) = amps(amp)*max(0d0, dsin(dble(step) / wls(wl)))
			CASE (10)
				act_accel(wl,amp,grain) = amps(amp)*max(0d0, dcos(dble(step) / wls(wl)))
			CASE(11)	
				act_accel(wl,amp,grain) = amps(amp) * (1d0 - 2d0 * dble(rand(0)))
			CASE (12)
				if (mod(int(step * sqrt(abs(grain_rand(wl,amp,grain)))),wls(wl)) .eq. 0) then
					
					grain_rand(wl,amp,grain) = 1d0 - rand(0) * 2d0
					!if (mod(i, 500) .eq. 0) then
					!	write(*,*) step, i, direction, grain_rand(direction,i)
					!endif
				endif
				act_accel(wl,amp,grain) = amps(amp) * max(0d0, dsin(dble(step) / wls(wl))) + & 
										amps(amp) * grain_rand(wl,amp,grain)/2	
			CASE DEFAULT
				act_accel(wl,amp,grain) = 0.
				stop
		END SELECT
						! draw a new factor
						
						act_accel(wl,amp,grain) =  dble(grain_rand(wl,amp,grain)) * amps(amp)
						abs_accel(wl,amp,grain) = abs(act_accel(wl,amp,grain))
						!write(changes_file, '(I4, "_", I6,"_",I3, ",",I10,",",I5,",", I3,",",E21.15 )') grain, &
						!		 wls(wl), amps(amp),  step - step_changes(wl, amp, grain), wls(wl), &
						!		amps(amp), act_accel(wl,amp,grain)
						changes(wl,amp,grain) = changes(wl,amp,grain) + 1
						step_changes(wl, amp, grain)= step
						!cont_changes(wl, amp, grain, changes(wl,amp,grain), 1) = step
						!cont_changes(wl, amp, grain, changes(wl,amp,grain), 2) = act_accel(wl,amp,grain)
					!endif
					act_total_accel(wl,amp,grain) = act_total_accel(wl,amp,grain) + act_accel(wl,amp,grain)
					abs_total_accel(wl,amp,grain) = abs_total_accel(wl,amp,grain) + abs_accel(wl,amp,grain)
					grain_vel(wl,amp,grain) = grain_vel(wl,amp,grain) + act_accel(wl,amp,grain)
					grain_vel_sq(wl,amp,grain) = grain_vel(wl,amp,grain)**2
					!if (mod(step, 1000) .eq. 0) then
					!	write(changes_file, '(",", E21.15)' , advance='no') grain_rand(wl,amp,grain)
					!endif
				enddo
			enddo	
		enddo
	!	if (mod(step, 1000) .eq. 0) then		
	!		write(changes_file, '(a)') ''		
	!	endif
	enddo
	!close(changes_file)
	!call sleep(1)
	open(averages, file='averages_averages.csv')
	!call sleep(1)	
	!open(totals, file='totals'//file_name)
	!call sleep(1)
	write(averages, '(a)') 'wl, amp, act, abs, changes,V, V_sq'
	!write(totals, '(a)') 'grain wl amp, wl, amp, act_accel, abs_acce, changes, V, V_sq'
	do amp=1, size(amps)
		do wl=1, sizes(amp)
			average = 0
			absvrage = 0
			avechange = 0
			avel_sq = 0
			avel = 0
			
			do grain=1, grains		
	!			write(totals,'(I4,"_",I6,"_", I3,",",I6, ",", I3, ",", E21.15, ",", E21.15, "," , I10,",",E21.15,",",E21.15)') & 
	!				grain, wls(wl), amps(amp), wls(wl), amps(amp), act_total_accel(wl,amp,grain), &
	!				abs_total_accel(wl,amp,grain), changes(wl,amp,grain), grain_vel(wl,amp,grain), grain_vel_sq(wl,amp,grain)
					absvrage = absvrage + abs_total_accel(wl,amp,grain)
					average = average + act_total_accel(wl,amp,grain)
					avechange = avechange + changes(wl,amp,grain)
					avel_sq = avel_sq + grain_vel_sq(wl,amp,grain)
					avel = avel + grain_vel(wl,amp,grain)
			enddo
			write(averages, '(I6, ",", I3, "," , E21.15, ",", E21.15, ",", I10,",",E21.15,",",E21.15)') wls(wl), amps(amp), average/grains, &
					absvrage/grains, int(avechange/grains), avel/grains, avel_sq/grains
			
		enddo	
	enddo	
	!call sleep(1)
	!close(totals)
	!call sleep(1)
	close(averages)
end	
	
