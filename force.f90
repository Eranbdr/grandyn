   

        subroutine force ()
     
        use mycommons

!    *******************************************************************
!    ** ROUTINE TO COMPUTE FORCES AND POTENTIAL USING A LINK LIST     **
!    **                                                               **
!    ** PRINCIPAL VARIABLES:                                          **
!    **                                                               **
!    ** INTEGER NGRAINS            NUMBER OF GRAINS                   **
!    ** INTEGER NGLIST(NGRAINS)    POINTER TO 1ST ATOM IN EACH GRAIN  **    
!    ** INTEGER MX,MY              NUMBER OF CELLS IN EACH DIRECTION  **
!    ** INTEGER NCELL              NUMBER OF SMALL CELLS (M**3)       **
!    ** INTEGER LIST(N)            THE LINKED LIST                    **
!    ** INTEGER HEAD(NCELL)        THE HEAD OF CHAIN ARRAY            **
!    ** INTEGER MAP(MAPSIZ)        LIST OF NEIGHBOURING CELLS         **
!    ** REAL*8    RX(N),RY(N)        POSITIONS                          **
!    ** REAL*8    FX(N),FY(N)        FORCES                             **
!    ** REAL*8    SIGMA              THE LJ LENGTH PARAMETER            **
!    ** REAL*8    RCUT               THE CUT-OFF DISTANCE               **
!    **                                                               **
!    **                                                               **
!    ** PROCEDURE:						      **
!    **   LOOP OVER THE CELLS  					      **
!    **     WITHIN THIS CELL				              **
!    **         DO INTRAMOLECULAR FORCES IN GRAINS		      **
!    **         LOOP THROUGH THE LINK LIST FOR GRAINS		      **
!    **	    IN NEIGHBORING CELLS				      **
!    **         LOOP THROUGH THE LINK LIST FOR GRAINS		      **
!    **                                                               **
!    **                                                               **
!    **  16/09/09 - Liran Goren.                                      **
!    **  Modified from force15(h).f contains changes with the way the   **
!    **  the fluid force and the confining force are calculated on    **
!    **  the walls.       
!    **  23/09/2014 - Eran BenDror
!    **  merged with force17.f                                            **   
!    *******************************************************************


      
      integer i,j,k, interacting, slipping
      integer ipf
      real*8 xl
      integer ixcell,jxcell

!    *******************************************************************

!    ** ZERO FORCES AND POTENTIAL **
!       print*, 'IN FORCE'

! gravity forces proportional to mass of particle, 
!  assuming equal density, proportional to radius^2
	
      call initf()

      ev = 0.0


  !    do i=1,n
  !        coordn(i) =0
   !   enddo
   


      xl=xright-xleft
      dxc=xl/dble(mx)
	interacting = 0
	slipping = 0
!  loop over each particle
!   loop over neighbor list for that particle
    do k=1,contactknt
		i=contacti(k)
		j=contactj(k)
		ipf=0
		if (mx .eq. 2) then 
		   ipf=2
		else
		   ixcell=1+int((r(1,i)-xleft)/dxc)
		   jxcell=1+int((r(1,j)-xleft)/dxc)
		   if ((ixcell .eq. 1 .and. jxcell .eq. mx) .or. &
				(ixcell .eq. mx .and. jxcell .eq. 1)) ipf = 1 
		endif
	 
     	!print*, i,j,ipf
		call interactions(i,j,k,ipf, slipping)
		if (contfn(k) .ne. 0 .or. contft(k) .ne. 0.) interacting = interacting + 1
    enddo
!     add the force due to the fluid gradient (liran)
      if (mod(step, check_condition_every) .eq. 0) then 
		print*,' number of interactions / contacts / slipping   = ',interacting, contactknt, slipping
	endif
      if (ifluid) call fluidf()
      


!        call zerocont2
        call bsum()
        
        RETURN
        END subroutine



       
                       
        
      subroutine initf()
!*****************************************************
! set up initial forces on particles
!  * body forces on internal particles
!  * applied forces on boundary particles  
      
      
      use mycommons
      integer i
      
      real*8 wallv
      integer ix,iy 
      real*8 yb1,xb1,dxr,dyr
      
      fgrain = 0.
      fgrain0 = 0.
      wallv = 0.
! curpresy, presy0 are applied stresses
! convert to force on each grain
      do i = 1,nbound(1)
         fgrain0 = fgrain0 + presy0*pi*radius(i)**2
      enddo
      fgrain0 = fgrain0/dble(nbound(1))
      do i = nbound(1)+1,nbound(2)
         fgrain = fgrain + curpresy*pi*radius(i)**2
      enddo
      fgrain = fgrain/dble(nbound(2) - nbound(1))
	  

 
! bottom boundary
      do i = 1, nbound(1)
           f(1,i) = 0.0
	   if (ib(1) .eq. 6 .or. ib(1) .eq. 3) f(1,i)=fb(1)
	   if (ib(1).eq.1 .or. ib(1).eq.2) then
		  f(2,i)=fgrain0 
	   else
		  f(2,i)=0.
	   endif
      enddo
      
! top boundary
      do i =nbound(1)+1,nbound(2)
           if (ib(2) .eq. 5) then
             f(1,i) = fspring  ! + 0.001*(fb(2) - v(1,i))
           else
             f(1,i)=0.
           endif
           if (ib(2).eq.1 .or. ib(2).eq.2 .or. ib(2) .eq. 5) then
             f(2,i) = -fgrain-wallv*v(2,nbound(2))
           else
              f(2,i)=0.
           endif
      enddo
!	calculate current acceleration modification angle
	curr_angle(1) = mod(angle_change_rate(1) * phase_step, 360d0)
	curr_angle(2) = mod(angle_change_rate(2) * phase_step, 360d0)
	
	!write(*,*) 'step ', step,' curr_angle(1) ',curr_angle(1),' curr_angle(2) ',curr_angle(2)
! internal particles
	do i = nbound(2)+1,n
		accel(1,i) = g(1) + acceleration_modification(1,i)
		accel(2,i) = g(2) + acceleration_modification(2,i)
		
		f(1,i) = accel(1,i)/radinv2(i)
		f(2,i) = accel(2,i)/radinv2(i)
		tq(i)=0.
	enddo
	
!	  real*8 acceleration_amplitude(2), acceleration_decay(2),angle_change_rate,final_angle,curr_angle
!	integer acceleration_function(2), check_condition_every
!	logical random_changes_per_grain,replace_original_forces,wait_for_condition_between_changes
	
!	print*, nbound(2),g(2),f(2,55)
      
!	print*, 'old force on top particle=',f(2,nbound(1)+1)
      return
      end subroutine

subroutine dump_forces()
	use mycommons
	integer i
	open(unit=forces_file,file=TRIM(output_directory)//'/forces'//TRIM(file_postfix)//fnn//'.csv')
	write(forces_file, *) 'step ', step, ' dt ',dt,' g(1) ',g(1),' g(2) ',g(2),' fgrain ',fgrain
	write(forces_file, *) 'angle ',curr_angle(2),' contacts ', contactknt
	write(forces_file, *) 'grain	', 'f(1)	','f(2)	'
	do i=1,n
		write(forces_file, *) i,',',gtype(i),',',f(1,i),',',f(2,i),',',v(1,i),',',v(2,i)
	enddo
	close(forces_file)

end subroutine

        subroutine bsum()
			use mycommons
!*************************************************************        
! get the sums of stresses on the rigid boundaries
	
	integer i
	real*8 xw
	
    
        
        xw=xright-xleft
! find average normal and shear forces on boundaries
        fytop=0.
        fybot=0.
        fxtop=0.
        fxbot=0.
        do i=1,nbound(1)
	     fybot=fybot + f(2,i)
	     fxbot=fxbot + f(1,i)
	enddo
        do i=nbound(1)+1,nbound(2)
	     fytop=fytop + f(2,i)
	     fxtop=fxtop + f(1,i)
	enddo
	fybot=fybot/dble(nbound(1))
	fxbot=fxbot/dble(nbound(1))
	fytop=fytop/dble(nbound(2)-nbound(1))
	fxtop=fxtop/dble(nbound(2)-nbound(1))

! to keep rigid, apply the average force to all boundary particles
        do i=1,nbound(1)
	     f(2,i)=fybot
	     f(1,i)=fxbot
	enddo
        do i=nbound(1)+1,nbound(2)
	     f(2,i)=fytop
	     f(1,i)=fxtop
	enddo
!	print*, 'new force on bottom particle=',f(2,nbound(1))
!	print*, 'new force on top particle=',f(2,nbound(1)+1)
	

 100    format(2(e14.5))
	

        RETURN
        END subroutine



        subroutine interactions(i,j,knt,ipf, slipping)
!**************************************************************
!   calculate if two particles are interacting, and the
!   interaction forces
!c
!   fnij - normal force
!   ftij - tangential force
!
!   forces on paritcle imn are postive if outward or counterclockwise
!   relative velocities of paritcle imx wrt. imn are 
!     postive if outward or counterclockwise

        use mycommons
	real*8 vxij,vyij,vnij,vthij,sk2,rij,rijsq,sk1,comp,vtij0
	real*8 fyij,fxij,fnij,vtij,ftmax,ftij,sk
	real*8 xl,cohesion,ryij,rxij,rmaxsq,eqlength, rxij2
	integer i,j,imx,imn, slipping
    real*8 slipij
	integer ipf
	integer knt, interacting
	real*8 gammax,g1,g2
	real*8 rbar, ar
    real*8 nu
!        real*8 gammat

 
         cohesion = 0.0
		 interacting = 0
         sk=1.  ! normal spring???
         xl=xright-xleft

! denote the two interacting grains as imn, imx
	      imn=min(i,j)
		  imx=max(i,j)

		  eqlength = radius(imn)+radius(imx)
		  rmaxsq=eqlength*eqlength

	      ryij  = r(2,imx) - r(2,imn)
	      rxij  = r(1,imx) - r(1,imn)
   	      if (ipf .eq. 1) then
               rxij  =rxij-dsign(xl,rxij)
   	      else if (ipf .eq. 2) then
               rxij2 = rxij-dsign(xl,rxij)
               if (dabs(rxij2) .lt. dabs(rxij)) rxij=rxij2
		  endif
		  rijsq = (rxij*rxij + ryij*ryij)

! particles are interacting
		if (rijsq .lt. rmaxsq) then
		
			last_contact_step(i) = step
			last_contact_step(j) = step
         !       coordn(i)=coordn(i)+1
          !      coordn(j)=coordn(j)+1
!	         print*, imn,imx
                rij=dsqrt(rijsq)
                comp=eqlength-rij

! stiffness for two different grains, length dependent
!	         sk1=emod(i)*emod(j)/
!     1             (radius(j)*emod(i)+radius(i)*emod(j))

!  LINEAR SPRING CONSTANTS
!  stiffness for E=1 for both grains, but length dependent
			sk1 = normal_spring_coef/eqlength !sk1=sk/eqlength
			sk2=shear_spring_coef/eqlength
			
!     potential energy from overlap  -  ERAN no use 
			ev = ev + 0.5 * sk1*comp*comp
!      get relative translational velocities
			 vxij=v(1,imx)-v(1,imn)
			 vyij=v(2,imx)-v(2,imn)

!      get relative rotations
			 vthij=w(imx)*radius(imx)+w(imn)*radius(imn)

!      resolve relative velocties into normal and tangential parts
			 vnij = (vxij*rxij + vyij*ryij)/rij
			 vtij0 = (-vxij*ryij + vyij*rxij)/rij
			 vtij = vtij0 - vthij

! following Potopov, keep track of tangential force on contact

!      update tangential force on contact
			 slipij=contft(knt) 
!      calculate the maximum allowable damping factor, gamma (size dependent)
			 g1=radius(imn)+radius(imx)
			 g2=radius(imn)*radius(imn)+radius(imx)*radius(imx)
			 !gammax=4.*radius(imn)*radius(imx)/dsqrt(g1+g2)

			if (.not. hertz_mindlin) then 
!      calculate normal and tangential forces
				fnij = -sk1*comp + gamma*vnij
				slipij = slipij + sk2*dt*vtij
		 
			else
!  HERTZ-MINDLIN CONTACT MODEL
!     harmonic mean of grain radii, rbar
!     contact radius,ar
				rbar=2.*radius(imn)*radius(imx)/(radius(imn)+radius(imx))
				ar=dsqrt(rbar*comp)

! calculate Poisson's ratio from shear_spring_coef
				nu = (6.-2.*shear_spring_coef)/(6.-shear_spring_coef)

!      calculate normal and tangential forces
	!			 fnij = -ar*comp + gamma*gammax*vnij
	!			 fnij = -ar*comp + gamma*ar*vnij
				fnij = -(dsqrt(2d0)/3./(1.-nu*nu))*ar*comp + gamma*ar*vnij
							  ! FOR MINDLIN TANGENTIAL CONTACT
! sk2=ratio of tangential stiffness to normal stiffness
!         =6*(1-nu)/(2-nu)  , where nu=Poisson's ratio
!                  slipij = slipij + shear_spring_coef*ar*dt*vtij
                  slipij = slipij + (2.*dsqrt(2d0)/(2.-nu)/(1.+nu))*ar*dt*vtij
		 
			endif
		
			 
			 ftmax = friction*fnij

			 if (dabs(slipij) .gt. dabs(ftmax)) then
!                     print*, 'hey Im slipping!'
				 slipping = slipping + 1
				 slipij=dsign(ftmax,slipij)
			 endif
			ftij=slipij
!	        if (vtij*slipij .gt. 0.) slipij=0.
                       

                contfn(knt)=fnij
                contft(knt)=slipij
				compression(knt) = comp
				
 		 
!      resolve forces back into x and y and add to sums for atom imn
                 fxij= (rxij*fnij-ryij*ftij)/rij
                 fyij= (ryij*fnij+rxij*ftij)/rij
                 f(1,imn)   = f(1,imn) + fxij
                 f(2,imn)   = f(2,imn) + fyij
                 tq(imn)   = tq(imn) + ftij


  205       format(7(1x,e12.5))
!   resolve opposite forces back into x and y and add to sums for imx
                 f(1,imx)   = f(1,imx) - fxij
                 f(2,imx)   = f(2,imx) - fyij
                 tq(imx)   = tq(imx) + ftij


	       else
		
! particles are not interacting
                contfn(knt)=0.
                contft(knt)=0.
            endif
        return
        end subroutine




! the fluid pressure force is gradP*m/(rho * (1-phi)) or in a different way: gradP*V/N where V is a cell volume and N in 
! the number of particle in the cell given in the array numpar
! after using non-dimensioalization we get that the fluid force is grad'P'*m'/(1-phi) * (3*pi/32)

      subroutine fluidf()
      use mycommons
      integer ix,iy,k 
      real*8 yb1,xb1,dxr,dyr, addx, addy,vol
      real*8 prefac
      

      prefac = pi/6.
      dxc=(xright-xleft)/mx
      dyc=(ytop-ybot)/my
      vol =dxc*dyc
! inner grains
      do k = nbound(2) + 1, n
         if (gtype(k) .ne. 0) then
            if (r(1,k) .lt. xleft .or. r(1,k) .gt. xright .or. r(2,k) .lt. ybot .or. r(2,k) .gt. ytop) then
               print*, k, r(1,k),r(2,k),xleft,xright,ybot,ytop
               stop
            endif         
            ix = int((r(1,k) - xleft)/dxc) + 1
            iy = int((r(2,k) - ybot)/dyc) + 1
            
            xb1=xleft+(ix-1)*dxc
            yb1=ybot+(iy-1)*dyc
            dxr=(r(1,k)-xb1)/dxc
            dyr=(r(2,k)-yb1)/dyc

            addx = (1.-dyr)*(1.-dxr)*gradpx(ix,iy)/solfra(ix,iy) + dyr*(1.-dxr)*gradpx(ix,iy+1)/solfra(ix,iy+1) + &
                   dyr*dxr*gradpx(ix+1,iy+1)/solfra(ix+1,iy+1) + (1.-dyr)*dxr*gradpx(ix+1,iy)/solfra(ix+1,iy)
            
            
            addx = addx*prefac*8*radius(k)**3
            f(1,k) = f(1,k) - addx 
            
            
            addy = (1.-dyr)*(1-dxr)*gradpy(ix,iy)/solfra(ix,iy) + dyr*(1.-dxr)*gradpy(ix,iy+1)/solfra(ix,iy+1) + &
                   dyr*dxr*gradpy(ix+1,iy+1)/solfra(ix+1,iy+1) + (1.-dyr)*dxr*gradpy(ix+1,iy)/solfra(ix+1,iy)
 

            addy = addy*prefac*8*radius(k)**3 
            f(2,k) = f(2,k) - addy
         endif
      enddo
      if (fbc .eq. 1 .or. fbc .eq. 5) then         
!     bottom grains
         if (fbc .eq. 1) then
            do k = 1, nbound(1)         
               ix = int((r(1,k) - xleft)/dxc) + 1
               iy = 1            
               xb1=xleft+(ix-1)*dxc
                yb1=ybot
                dxr=(r(1,k)-xb1)/dxc
                dyr= 0.
                addx = 1.*(1.-dxr)*gradpx(ix,iy)/solfra(ix,iy) + 1.*dxr*gradpx(ix+1,iy)/solfra(ix+1,iy)
!     mass of boundary grains is half
                addx = addx*prefac*8*radius(k)**3/2 
                f(1,k) = f(1,k) - addx 
                
                
                addy = 1.*(1-dxr)*gradpy(ix,iy)/solfra(ix,iy) + 1.*dxr*gradpy(ix+1,iy)/solfra(ix+1,iy)
                
                addy = addy*prefac*8*radius(k)**3/2 
                f(2,k) = f(2,k) - addy 
             enddo
          endif
!     top grains
          do k = nbound(1) + 1,nbound(2)
             ix = int((r(1,k) - xleft)/dxc) + 1
             iy = my+1            
             xb1=xleft+(ix-1)*dxc
             yb1=ytop
             dxr=(r(1,k)-xb1)/dxc
             dyr= 0.
             addx = 1.*(1.-dxr)*gradpx(ix,iy)/solfra(ix,iy) + 1.*dxr*gradpx(ix+1,iy)/solfra(ix+1,iy)
!     mass of boundary grains is half
             addx = addx*prefac*8*radius(k)**3/2 
             f(1,k) = f(1,k) - addx             
             addy = 1.*(1-dxr)*gradpy(ix,iy)/solfra(ix,iy) + 1.*dxr*gradpy(ix+1,iy)/solfra(ix+1,iy)
             
             addy = addy*prefac*8*radius(k)**3/2 
             f(2,k) = f(2,k) - addy 
          enddo
       endif
       if (fbc .eq. 2) then
!  bottom grains
          do k = 1,nbound(1)
             ix = int((r(1,k) - xleft)/dxc) + 1          
             xb1=xleft+(ix-1)*dxc            
             dxr=(r(1,k)-xb1)/dxc
!     when calculating the force as pressure* surface area
             addy = ((1-dxr)*press(ix,1) + dxr*press(ix+1,1))*pi*radius(k)**2
             f(2,k) = f(2,k) - addy
          enddo
! top grains
          do k = nbound(1)+1,nbound(2)
            ix = int((r(1,k) - xleft)/dxc) + 1          
            iy = int((ytop - ybot)/dyc) + 1
            xb1=xleft+(ix-1)*dxc            
            dxr=(r(1,k)-xb1)/dxc
!     when calculating the force as pressure* surface area
            addy = ((1-dxr)*press(ix,iy) + dxr*press(ix+1,iy))*pi*radius(k)**2
!     here we use '+' sign because positive PP sould act to push the top wall up  
            f(2,k) = f(2,k) + addy
         enddo
      endif
      return 
      end subroutine


