!***************************************************************
!**  VELOCITY VERSION OF VERLET ALGORITHM                           **
!***************************************************************

!    *******************************************************************
!    ** TWO ROUTINES THAT TOGETHER IMPLEMENT VELOCITY VERLET METHOD.  **
!    **                                                               **
!    ** REFERENCE:                                                    **
!    **                                                               **
!    ** SWOPE ET AL., J. CHEM. PHYS. 76, 637, 1982.                   **
!    **                                                               **
!    ** ROUTINES SUPPLIED:                                            **
!    **                                                               **
!    ** SUBROUTINE MOVEA ( DT, N, NTOT )                                    **
!    **    MOVES POSITIONS AND PARTIALLY UPDATES VELOCITIES.          **
!    ** SUBROUTINE MOVEB ( DT, N, NTOT)                                **
!    **    COMPLETES VELOCITY MOVE AND CALCULATES KINETIC ENERGY.     **
!    **                                                               **
!    ** PRINCIPAL VARIABLES:                                          **
!    **                                                               **
!    ** INTEGER N                   NUMBER OF MOLECULES               **
!    ** REAL*8    DT                  TIMESTEP                          **
!    ** REAL*8    M                   ATOMIC MASS                       **
!    ** REAL*8    RX(N),RY(N)   POSITIONS                         **
!    ** REAL*8    VX(N),VY(N)   VELOCITIES                        **
!    ** REAL*8    FX(N),FY(N)   FORCES                            **
!    **                                                               **
!    ** USAGE:                                                        **
!    **                                                               **
!    ** AT THE START OF A TIMESTEP, MOVEA IS CALLED TO ADVANCE THE    **
!    ** POSITIONS AND 'HALF-ADVANCE' THE VELOCITIES.  THEN THE FORCE  **
!    ** ROUTINE IS CALLED, AND THIS IS FOLLOWED BY MOVEB WHICH        **
!    ** COMPLETES THE ADVANCEMENT OF VELOCITIES.                      **
!    *******************************************************************



        SUBROUTINE MOVEA ()


!    *******************************************************************
!    ** FIRST PART OF VELOCITY VERLET ALGORITHM                       **
!    **                                                               **
!    ** USAGE:                                                        **
!    **                                                               **
!    ** THE FIRST PART OF THE ALGORITHM IS A TAYLOR SERIES WHICH      **
!    ** ADVANCES POSITIONS FROM T TO T + DT AND VELOCITIES FROM       **
!    ** T TO T + DT/2.  AFTER THIS, THE FORCE ROUTINE IS CALLED.      **
!    *******************************************************************

        use mycommons

        INTEGER     i, layer
	real*8 	 tmass1, rel_ry

	real*8 xl, prev_rx, drx, curr_height, prev_distance_from_center, curr_distance_from_center
	real*8 prev_central_angle, curr_central_angle, grain_flux, grain_vol_flux

!    *******************************************************************
   


	tmass1=1./topmass
        
		        
!  bottom boundary
!    ib=0, fixed wall
	if (ib(1) .eq. 0) then
           do i=1,nbound(1)
              v(1,i)=0.
              v(2,i)=0.
           enddo
        
	else if (ib(1) .eq. 1) then
!          applied vertical force, free horizontal
	do i=1,nbound(1)
		      r(1,i) = r(1,i) + dt * v(1,i) + dtsq2*f(1,i) 
		      r(2,i) = r(2,i) + dt * v(2,i) + dtsq2*fybot
		      v(1,i) = v(1,i) + dt2* f(1,i) 
		      v(2,i) = v(2,i) + dt2* fybot
 	enddo
	else if (ib(1) .eq. 2) then
!          applied vertical force, applied horizontal velocity
	do i=1,nbound(1)
		      r(1,i) = r(1,i) + dt * v(1,i)
		      r(2,i) = r(2,i) + dt * v(2,i) + dtsq2*fybot
		      v(1,i) = fb(1) 
		      v(2,i) = v(2,i) + dt2* fybot
 	enddo
	else if (ib(1) .eq. 3) then
!          applied vertical velocity, free horizontal 
	do i=1,nbound(1)
		      r(1,i) = r(1,i) + dt * v(1,i) + dtsq2*f(1,i)
		      r(2,i) = r(2,i) + dt * v(2,i) 
		      v(1,i) = v(1,i) + dt2* f(1,i)
		      v(2,i) = fb(1)
 	enddo
	else if (ib(1) .eq. 4) then
!          zero vertical velocity, applied horizontal velocity
	do i=1,nbound(1)
		      r(1,i) = r(1,i) + dt * v(1,i)
		      v(1,i) = fb(1) 
		      v(2,i) = 0. 
 	enddo
	else if (ib(1) .eq. 5) then
!          applied vertical velocity, zero horizontal velocity
	do i=1,nbound(1)
		      r(2,i) = r(2,i) + dt * v(2,i)
		      v(2,i) = fb(1) 
		      v(1,i) = 0. 
 	enddo


	else if (ib(1) .eq. 6) then
!          zero vertical velocity, applied horizontal force
	do i=1,nbound(1)
		      r(1,i) = r(1,i) + dt * v(1,i) + dtsq2*fxbot 
		      v(1,i) = v(1,i) + dt2* fxbot 
		      v(2,i) = 0.
 	enddo
	endif

	

!  top boundary
	if (ib(2) .eq. 0.) then
	 do i=nbound(1)+1,nbound(2)
	  v(1,i)=0.
	  v(2,i)=0.
	 enddo
	else if (ib(2) .eq. 1 .or. ib(2) .eq. 5) then
	  do i=nbound(1)+1,nbound(2)
	      r(1,i) = r(1,i) + dt * v(1,i) + dtsq2*f(1,i)*tmass1
	      r(2,i) = r(2,i) + dt * v(2,i) + dtsq2*fytop*tmass1
	      v(1,i) = v(1,i) + dt2* f(1,i)*tmass1 
	      v(2,i) = v(2,i) + dt2* fytop*tmass1
 	  enddo
	else if (ib(2) .eq. 2) then
	  do i=nbound(1)+1,nbound(2)
	      r(1,i) = r(1,i) + dt * v(1,i)
	      r(2,i) = r(2,i) + dt * v(2,i) + dtsq2*fytop*tmass1
	      v(1,i) = fb(2) 
	      v(2,i) = v(2,i) + dt2* fytop*tmass1
 	  enddo
	else if (ib(2) .eq. 3) then
	  do i=nbound(1)+1,nbound(2)
	      r(1,i) = r(1,i) + dt * v(1,i) + dtsq2*f(1,i)*tmass1
	      r(2,i) = r(2,i) + dt * v(2,i) 
	      v(1,i) = v(1,i) + dt2* f(1,i)*tmass1
	      v(2,i) = fb(2)
 	  enddo
	else if (ib(2) .eq. 4) then
	  do i=nbound(1)+1,nbound(2)
	      r(1,i) = r(1,i) + dt * v(1,i)
	      v(1,i) = fb(2) 
 	  enddo
	endif
!	print*, 'vel on top particle=',v(2,nbound(1)+1)
!	print*, 'vel on bottom particle=',v(2,nbound(1))
	 
! interior grains
	
	curr_height = internal_top - ybot	
    do  i=nbound(2)+1,n
		 
		 !layer = int((r(2,i) - ybot) / number_of_layers) + 1
		 layer = floor(number_of_layers*((r(2,i) - ybot)/phase_start_height)) + 1 
		 if (layer .gt. number_of_layers) then
			! higher than original number of layers
			!write(*, *) 'grain ', i, ' wanted to be in layer ', layer, ' set to ', number_of_layers
			layer = number_of_layers
		endif
		if (layer .ne. grain_output_layer(i)) then
			!write(layer_change_file,*) step, ',', i, ',' , grain_output_layer(i), ',', layer
			grain_output_layer(i) = layer
			
		endif
		 ! calculate horizontal flux as area across the rx = 0 line
		 prev_rx = r(1,i)
		 drx = dt * v(1,i) + dtsq2 * f(1,i) *radinv2(i)
		 r(1,i) = r(1,i) + drx
		 if (dabs(prev_rx) .lt. radius(i) .or. dabs(r(1,i)) .lt. radius(i)) then
			
			! grain moved across the x == 0 line
			! we need central angles, that's why we multiply by 2
			prev_central_angle = 2.*dacos(min(1d0,dabs(prev_rx)/radius(i)))
			curr_central_angle = 2.*dacos(min(1d0,dabs(r(1,i))/radius(i)))
			if (prev_rx * r(1,i) .ge. 0) then ! grain center didn't cross x == 0
				! calculate area of delta circle segment
				! dA = |A2 - A1| = 0.5R^2|(alfa2 - sin(alfa2)) - (alfa1 - sin(alfa1))|
				grain_flux = dabs((curr_central_angle - dsin(curr_central_angle) - &
								  (prev_central_angle - dsin(prev_central_angle)))*0.5*radius(i)**2)
								  
				! calculate volume of delta sphere segment
				! http://www.ambrsoft.com/Equations/Circle/EQ-Circle.htm
				! dV = pi(R^2|a-b|-|a^3-b^3|/3)
				grain_vol_flux = pi*((radius(i)**2)*abs(drx)) - abs((prev_rx**3) - r(1,i)**3)/3d0 
				grains_per_column_per_layer(layer) = grains_per_column_per_layer(layer) + 1
			else 
				! grain center moved across x == 0
				! dA = A - (A2 + A1)
				grain_flux = (pi - 0.5*(curr_central_angle - dsin(curr_central_angle) + &
							  prev_central_angle - dsin(prev_central_angle)))*radius(i)**2	
				
				! dV = Vsphere - (Vcap1 + Vcap2) = 4pi/3*R^3 - pi/3*(R-a)^2(3R-(R-a)) - pi/3*(R-b)^2(2R+b)
				grain_vol_flux = (4*(radius(i)**3) - ((radius(i) - prev_rx)**2)*(2*radius(i) + prev_rx) - &
								((radius(i) - r(1,i))**2)*(2*radius(i) + r(1,i)))*pi/3
				grains_per_column_per_layer(layer) = grains_per_column_per_layer(layer) + 1
			endif
			! set flux direction
			grain_flux = dsign(grain_flux, drx)
			flux(layer) = flux(layer) + grain_flux
			grain_vol_flux = dsign(grain_vol_flux, drx)
			vol_flux(layer) = vol_flux(layer) + grain_vol_flux
			
		 endif
		 r(2,i) = r(2,i) + dt * v(2,i) + dtsq2 * f(2,i) *radinv2(i)
		 v(1,i) = v(1,i) + dt2 * f(1,i) *radinv2(i)
		 v(2,i) = v(2,i) + dt2 * f(2,i) *radinv2(i)
    enddo
	 
	
           
        

!  ugly bit to take care of grain wraparound
!  in case of periodic boundaries
        if (ib(3) .eq. -1) then
          xl=xright-xleft
     
		do  i=1,n
		  if (r(1,i) .ge. xright) then 
		     r(1,i)=r(1,i)-xl
		  else if (r(1,i) .lt. xleft) then
		     r(1,i)=r(1,i)+xl
		  endif
		enddo	           

        endif
        
        ybot=r(2,2)
        ytop=r(2,nbound(1)+2)
        if (ib(3) .ne. -1) then
         xleft=r(1,1)
         xright=r(1,nbound(1))
        endif
         

        RETURN
        END



        SUBROUTINE MOVEB ()

        use mycommons

!    *******************************************************************
!    ** SECOND PART OF VELOCITY VERLET ALGORITHM                      **
!    **                                                               **
!    ** USAGE:                                                        **
!    **                                                               **
!    ** THE SECOND PART OF THE ALGORITHM ADVANCES VELOCITIES FROM     **
!    ** T + DT/2 TO T + DT. THIS ASSUMES THAT FORCES HAVE BEEN        **
!    ** COMPUTED IN THE FORCE ROUTINE AND STORED IN FX, FY, FZ.       **
!    *******************************************************************
 

        INTEGER     i
	real*8 tmass1
!    *******************************************************************
      


	tmass1=1./topmass



!  bottom boundary
	if (ib(1) .eq. 1) then
	do i=1,nbound(1)
		      v(1,i) = v(1,i) + dt2* f(1,i) 
		      v(2,i) = v(2,i) + dt2* fybot
 	enddo
	else if (ib(1) .eq. 2) then
	do i=1,nbound(1)
		      v(1,i) = fb(1) 
		      v(2,i) = v(2,i) + dt2* fybot
 	enddo
	else if (ib(1) .eq. 3) then
	do i=1,nbound(1)
		      v(1,i) = v(1,i) + dt2* f(1,i)
		      v(2,i) = fb(1)
 	enddo
	else if (ib(1) .eq. 4) then
	do i=1,nbound(1)
		      v(1,I) = fb(1) 
 	enddo
	else if (ib(1) .eq. 5) then
!          applied vertical velocity, zero horizontal velocity
	do i=1,nbound(1)
		      v(2,I) = fb(1) 
 	enddo
	else if (ib(1) .eq. 6) then
!          zero vertical velocity, applied horizontal force
	do i=1,nbound(1)
		      v(1,i) = v(1,i) + dt2* fxbot 
		      v(2,i) = 0.
 	enddo
	endif

	

!  top boundary
	if (ib(2) .eq. 1 .or. ib(2) .eq. 5) then
	  do i=nbound(1)+1,nbound(2)
	      v(1,i) = v(1,i) + dt2* f(1,i)*tmass1 
	      v(2,i) = v(2,i) + dt2* fytop*tmass1
 	  enddo
	else if (ib(2) .eq. 2) then
	  do i=nbound(1)+1,nbound(2)
	      v(1,i) = fb(2) 
	      v(2,i) = v(2,i) + dt2* fytop*tmass1
 	  enddo
	else if (ib(2) .eq. 3) then
	  do i=nbound(1)+1,nbound(2)
	      v(1,i) = v(1,i) + dt2* f(1,i)*tmass1
	      v(2,i) = fb(2)
 	  enddo
	else if (ib(2) .eq. 4) then
	  do i=nbound(1)+1,nbound(2)
	      v(1,i) = fb(2) 
 	  enddo
	endif
	
	 
! interior grains
	do i=nbound(2)+1,n
		v(1,i) = v(1,i) + dt2 * f(1,i) * radinv2(i)
		v(2,i) = v(2,i) + dt2 * f(2,i) * radinv2(i)
	enddo 
           

        RETURN
        END


