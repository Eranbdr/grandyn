PROGRAM find_accel
	IMPLICIT NONE
	integer grains, steps, numps
	parameter (steps = 2000000)
	parameter (grains = 500)
	parameter (numps = 4)
	real, dimension(7,numps,grains) :: grain_rand, act_accel, abs_accel, act_total_accel, &
										abs_total_accel, grain_vel, grain_vel_sq
	integer :: wls(7), amps(numps), sizes(4), wl, amp, grain, changes_file, totals, step, averages, avechange
	integer :: changes(7,numps,grains), step_changes(7,numps,grains)
	real :: absvrage, average, avel, avel_sq!,  cont_changes(7,numps,grains, 600,2)
	wls= (/10000, 12000, 15000 ,5000, 8000 ,18000, 20000 /)
	!wls= (/5000, 8000 ,10000, 12000, 15000 ,18000, 20000 /)
	amps = (/11, 19, 25, 31 /)
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
	!open(unit=changes_file,file='changes__NN_19.csv')
	!write(changes_file, '(a)') 'grain wl amp, steps, wl, amp, act_accel '
	
	do amp=1, numps
		do wl=1, sizes(amp)
			do grain=1, grains
				!write(changes_file, '(",",I4, "_", I6,"_",I3)' , advance='no') grain, wls(wl), amps(amp)
			enddo
		enddo
	enddo
	!write(changes_file, '(a)') ''
	do step=1, steps
	!	if (mod(step, 1000) .eq. 0) then
	!		write(changes_file, '(I10)' , advance='no') step
	!	endif
		do amp=1, numps
			do wl=1, sizes(amp)
				do grain=1, grains
					if (mod(int(step * sqrt(abs(grain_rand(wl,amp,grain)))), wls(wl)) .eq. 0) then
						! draw a new factor
						grain_rand(wl,amp,grain) = 1d0 - rand(0) * 2d0
						act_accel(wl,amp,grain) =  dble(grain_rand(wl,amp,grain)) * amps(amp)
						abs_accel(wl,amp,grain) = abs(act_accel(wl,amp,grain))
						!write(changes_file, '(I4, "_", I6,"_",I3, ",",I10,",",I5,",", I3,",",E21.15 )') grain, &
						!		 wls(wl), amps(amp),  step - step_changes(wl, amp, grain), wls(wl), &
						!		amps(amp), act_accel(wl,amp,grain)
						changes(wl,amp,grain) = changes(wl,amp,grain) + 1
						step_changes(wl, amp, grain)= step
						!cont_changes(wl, amp, grain, changes(wl,amp,grain), 1) = step
						!cont_changes(wl, amp, grain, changes(wl,amp,grain), 2) = act_accel(wl,amp,grain)
					endif
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
	open(averages, file='averages__NN_23.csv')
	!call sleep(1)
	!open(totals, file='totals__NN_19.csv')
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
	
