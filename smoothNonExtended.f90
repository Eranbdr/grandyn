!cccccccccccccccc 15/2/09 cccccccccccccccccccccccccccccccc


!	Here we assign values to grid points representing some continous properties of the fluid+grain system.
!       the values are assign using the halo function which is a 2D linear interpolation fuction which allows
!       calculating the relative weight of each coordinate in the system to the grid point arounf it.
!
!       We use staggared grid  so that 
!       P, k and phi are given at grid point (i,j)
!       ux is between x grid points (i+0.5,j)
!       uz is between y grid points (i,j+0.5)
!       
! 
!       vc -  is a 2D array of the area that each grain contributes to the grid point. On grid
!       siterd2 -  is a 2D array of the average squre of radiuses. On grid. 
!                  This is not a conservative quantity because of division by numpar 
!       sitevx and sitevy - are 2D arrays of grains velocity in each grid point. Staggared
!       divvel - is a 2D array of the divergence of the velocity. On grid
!	solfra - is a 2D array of the bulk density of system grains area/total area. On grid
!                Not a conserved quantity because of boundary effect.
!	perm - is a 2D array of the permeability. On grid  
!	phimat - is a 2D array of the porosity. On grid
!       numpar - is a 2D array of the number of particles in each grid point, used for averaging radius^2. On grid 
!       numvx - is a 2D array of the number of particles in each v(1) grid point, used for averaging
!       numvy - is a 2D array of the number of particles in each v(2) grid point, used for averaging

!     Note that there are two more lines of grid point for sitevy and numvy at the top and at the bottom
!     and one more lines of grid point for phimat,perm 
!
!     This file should go together with updatePressure4
!     It was contracted from smoothExtended by removing the code of 
!     the top and bottom extra grains layers.
	
	subroutine smooth()
        use mycommons
        real*8 gs
        real*8 dxr,dyr
        real*8 xb1,yb1
        integer ix,iy,k
!        real*8 vc(MAXD+1,MAXD+1)
!	real*8 siterd2(MAXD+1, MAXD+1)
	real*8 numvx(MAXD+1, MAXD+1)
	real*8 numvy(MAXD+1, MAXD+1)
	real*8 rsmin
!       minimum solfra is rsmin because below that the Carman-Kozeny permeability law fails.
	parameter (rsmin=0.25)
	real*8 maxsiterd2
	integer ixvxh,iyvyh
	real*8 xbvx1,ybvy1,dxrvx,dyrvy
	
     	
        dxc=(xright-xleft)/mx
        dyc=(ytop-ybot)/my
!	print*,'ytop = ',ytop
!	print*,'yboy = ',ybot
!	print*,'dyc = ',dyc

! zeroing out the arrays
! my - # of cells in the y direction
! my+1 # of grid points lines inside the boc
! my+2 # number of grid point lines for y velocity
			

	do ix=1,mx+1
	   do iy=1,my+2
	      vc(ix,iy)=0.
	      siterd2(ix,iy)=0.
	      sitevx(ix,iy)=0.
	      sitevy(ix,iy)=0.
	      divvel(ix,iy)=0.
	      solfra(ix,iy)=0.
	      perm(ix,iy) = 1.			                
	      phimat(ix,iy) = 0.
	      numpar(ix,iy)=0.
	      numvx(ix,iy)=0.
	      numvy(ix,iy)=0.
              sitevy(ix,iy)=0.
              numvy(ix,iy) = 0.
           enddo
	enddo
		   
! ix,iy for regular grid
! ixvxh for staggared v(1)
! iyvyh for staggared v(2)
!	print*, 'n = ',n
	do k = 1,n
!	   print*, 'k = ', k
	   gs=pi*radius(k)**2
	   if (bdgrain(k) .gt. 0) gs=0.5*gs
	   ix = int((r(1,k) - xleft)/dxc) + 1
	   if (r(1,k) .lt. xleft + 0.5*dxc) then
	      ixvxh = mx
	   else
	      ixvxh = int((r(1,k)-(xleft + 0.5*dxc))/dxc)+1
	   endif
	   if (k .lt. nbound(1)+1 .or. k .gt. nbound(2))then       
	      iy = int((r(2,k) - ybot)/dyc) + 1
!              print*, '**********'
!              print*, k,  nbound(1),  nbound(2), r(2,k),  ybot, iy
	   else
	      iy = int((r(2,k) - ybot)/dyc) + 1
!              print*,  '**********'
!              print*, k, r(2,k),  my, iy, ybot, ytop, dyc
	   endif
	   if (r(2,k) .lt. ybot + 0.5*dyc) then
	      iyvyh = 1
	   else
	      if (r(2,k) .eq. ytop) then
		 iyvyh = my+1
	      else
		 iyvyh=int((r(2,k)-(ybot + 0.5*dyc))/dyc)+2
	      endif
	   endif
	   if (ixvxh .eq. mx +1) then
	      print*, '^^^',k
	      print*,'***',r(1,k),ix,dxc 
	      print*,'***',r(2,k)
	      print*,'^^^'
	   endif
	   if (ix .lt.1 .or. ix .gt. mx+1 .or. iy .lt. 1 .or. iy .gt. my+1) then
	      print*,'***',k
	      print*,'***',r(1,k),ix,dxc
	      print*,'***',r(2,k),iy,dyc
	      print*,'---',ybot,ytop
	   endif
! xb1, yb1 - coordinates of grid
! xbvx1,ybvy1 - coordinate of the relevant staggared 

	   xb1=xleft+(ix-1)*dxc
	   yb1=ybot+(iy-1)*dyc
           if (iy .eq. 1) yb1 = ybot
           if (iy .eq. my+1) yb1 = ytop
	   if (ixvxh .eq. mx) then
	      xbvx1 = xright - 0.5*dxc
	   else
	      xbvx1 = xleft+(ixvxh-0.5)*dxc 
	   endif
	   ybvy1 = ybot + (iyvyh - 1.5)*dyc 
	   if(iyvyh .eq. 1) ybvy1 = ybot - 0.5*dyc          
	   if(iyvyh .eq. my+2) ybvy1 = ytop + 0.5*dyc


! scaled distance of grain from grid point
		
!	   dxr=(r(1,k)-xb1)/dxc		    
!	   if (dxr .lt. 0.) dxr = 0. 
!	   dyr=(r(2,k)-yb1)/dyc
!	   if (dyr .gt. 1.) dyr = 1.		   
!	   if (dyr .lt. 0.) dyr = 0.
	   dxr=(r(1,k)-xb1)/dxc		    
	   if (dxr .lt. 0.) print*, 'dxr: ',dxr
	   if (dxr .gt. 1.) print*, 'dxr: ',dxr
	   dyr=(r(2,k)-yb1)/dyc
	   if (dyr .gt. 1.) print*, 'dyr: ',dyr, k		   
	   if (dyr .lt. 0.) print*, 'dyr: ',dyr, k

	   dxrvx = (r(1,k) - xbvx1)/dxc
	   if (r(1,k) .lt. xleft + 0.5*dxc) then 
	      dxrvx = (r(1,k) - xleft + 0.5*dxc)/dxc
	   endif
	   dyrvy = (r(2,k) - ybvy1)/dyc
           
!     debugging	   
	   if (dxrvx .gt. 1. .or. dxrvx .lt. 0.) then
!	      print*, 'dxrvx is : ',dxrvx
	   endif
	   if (dyrvy .gt. 1. .or. dyrvy .lt. 0. ) then
	      print*, 'dyrvy is : ',dyrvy
	      print*, r(2,k),ybvy1,dyc,(r(2,k) - ybvy1)
	   endif
!          starting to calculate relevant quantities
!          bootm left
	   vc(ix,iy)=vc(ix,iy)+((1.-dyr)*(1.-dxr)*gs)	
	   siterd2(ix,iy)=siterd2(ix,iy)+(1.-dyr)*(1.-dxr)*(radius(k)**2)
	   if (gtype(k) .ne. 0) then
	      numpar(ix,iy)=numpar(ix,iy)+(1.-dyr)*(1.-dxr) 
	   endif
	   sitevx(ixvxh,iy)=sitevx(ixvxh,iy)+v(1,k)*(1.-dyr)*(1.-dxrvx)
	   numvx(ixvxh,iy)=numvx(ixvxh,iy)+(1.-dyr)*(1.-dxrvx)
	   sitevy(ix,iyvyh)=sitevy(ix,iyvyh)+v(2,k)*(1.-dyrvy)*(1.-dxr)
	   numvy(ix,iyvyh)=numvy(ix,iyvyh)+(1.-dyrvy)*(1.-dxr)
	   		
!         top left	  
           vc(ix,iy+1)=vc(ix,iy+1)+(dyr*(1.-dxr)*gs)
           siterd2(ix,iy+1) = siterd2(ix,iy+1) + dyr*(1.-dxr)*(radius(k)**2)
           if (gtype(k) .ne. 0) then
              numpar(ix,iy+1) =numpar(ix,iy+1)+dyr*(1.-dxr) 
           endif
!           if (gtype(k) .eq. 0) print*, 'gtype(',k,')=0'
           sitevx(ixvxh,iy+1) = sitevx(ixvxh,iy+1) + v(1,k)*dyr*(1.-dxrvx)
           numvx(ixvxh,iy+1) = numvx(ixvxh,iy+1) + dyr*(1.-dxrvx)
       
           sitevy(ix,iyvyh+1) = sitevy(ix,iyvyh+1) + v(2,k)*dyrvy*(1.-dxr)
           numvy(ix,iyvyh+1) = numvy(ix,iyvyh+1) + dyrvy*(1.-dxr)
	  
!          bottom right           
	   if (ix .lt. mx) then
	      vc(ix+1,iy)=vc(ix+1,iy)+((1.-dyr)*dxr*gs)
	      siterd2(ix+1,iy) = siterd2(ix+1,iy) + (1.-dyr)*dxr*(radius(k)**2)
	      if (gtype(k) .ne. 0) then
		 numpar(ix+1,iy) =numpar(ix+1,iy) + (1.-dyr)*dxr
	      endif
	      sitevy(ix+1,iyvyh) = sitevy(ix+1,iyvyh) + v(2,k)*(1.-dyrvy)*dxr
	      numvy(ix+1,iyvyh) = numvy(ix+1,iyvyh) + (1.-dyrvy)*dxr
!          top right
              vc(ix+1,iy+1)=vc(ix+1,iy+1)+(dyr*dxr*gs) 
              siterd2(ix+1,iy+1) = siterd2(ix+1,iy+1) + dyr*dxr*(radius(k)**2)
              if (gtype(k) .ne. 0) then
                 numpar(ix+1,iy+1)=numpar(ix+1,iy+1)+dyr*dxr
              endif
              
              sitevy(ix+1,iyvyh+1) = sitevy(ix+1,iyvyh+1)+v(2,k)*dyrvy*dxr
              numvy(ix+1,iyvyh+1) = numvy(ix+1,iyvyh+1)+dyrvy*dxr              
!C if ix = mx
	   else
!          bottom right
	      vc(1,iy)=vc(1,iy)+((1.-dyr)*dxr*gs)
	      siterd2(1,iy) = siterd2(1,iy) + (1.-dyr)*dxr*(radius(k)**2)
	      if (gtype(k) .ne. 0) then
		 numpar(1,iy) =numpar(1,iy) + (1.-dyr)*dxr
	      endif
	      sitevy(1,iyvyh) = sitevy(1,iyvyh) + v(2,k)*(1.-dyrvy)*dxr
	      numvy(1,iyvyh)=numvy(1,iyvyh)+(1.-dyrvy)*dxr
!           top right	      
              vc(1,iy+1)=vc(1,iy+1)+(dyr*dxr*gs) 
              siterd2(1,iy+1) = siterd2(1,iy+1) + dyr*dxr*(radius(k)**2)
              if (gtype(k) .ne. 0) then
                 numpar(1,iy+1)=numpar(1,iy+1)+dyr*dxr 
              endif
              sitevy(1,iyvyh+1) = sitevy(1,iyvyh+1) + v(2,k)*dyrvy*dxr
              numvy(1,iyvyh+1)=numvy(1,iyvyh+1)+dyrvy*dxr
	   endif 

	   if (ixvxh .lt. mx) then
	      sitevx(ixvxh+1,iy) = sitevx(ixvxh+1,iy) + v(1,k)*(1.-dyr)*dxrvx
	      numvx(ixvxh+1,iy) = numvx(ixvxh+1,iy) + (1.-dyr)*dxrvx
              sitevx(ixvxh+1,iy+1) = sitevx(ixvxh+1,iy+1)+v(1,k)*dyr*dxrvx
              numvx(ixvxh+1,iy+1) = numvx(ixvxh+1,iy+1) + dyr*dxrvx
	   else
	      sitevx(1,iy)=sitevx(1,iy)+v(1,k)*(1.-dyr)*dxrvx
	      numvx(1,iy) = numvx(1,iy) + (1.-dyr)*dxrvx
              sitevx(1,iy+1)=sitevx(1,iy+1)+v(1,k)*dyr*dxrvx
              numvx(1,iy+1) = numvx(1,iy+1) + dyr*dxrvx
	   endif

!	   if (k .eq. 231) then
!	      print*, 'k is ',k
!	      print*, ybot,ytop,r(2,k)
!	      print*, iyvyh,ybvy1,dyrvy 
!	      print*, v(2,k)
!	      print*, '***************'
!	   endif
	enddo 
        





        do iy=1,my+1
	   vc(mx+1,iy)=vc(1,iy) 
	enddo


	do iy=1,my+1
	   do ix =1,mx+1
	      solfra(ix,iy) = vc(ix,iy)/(dxc*dyc)
	      if (solfra(ix,iy) .gt. 1.) print*, 'solfra(',ix,iy,')=',solfra(ix,iy),vc(ix,iy),dxc,dyc,my
	      if (solfra(ix,iy) .lt. rsmin) then
!		 if (toprint .gt. 0 )  print*,'solfra(',ix,iy,')=',solfra(ix,iy) 
		 solfra(ix,iy) = rsmin
	      endif
	   enddo
	enddo	
   

	do iy=1,my+1
	   do ix=1,mx
	      phimat(ix,iy)=1.- solfra(ix,iy)
	      if (phimat(ix,iy) .lt. 0.)then
!		 print*,'phimat(',ix,iy,')=',phimat(ix,iy)  
		 phimat(ix,iy) = 0. 
		 
	      endif
	   enddo
	enddo
	
	do iy=1,my+2
	   phimat(mx+1,iy)=phimat(1,iy)
	   siterd2(mx+1,iy)=siterd2(1,iy)
	   numpar(mx+1,iy)=numpar(1,iy)
	   sitevy(mx+1,iy) = sitevy(1,iy) 
	   numvy(mx+1,iy) = numvy(1,iy) 
      	enddo

	


! CALCULATE THE PERMEABILITY ACCORDING TO CARMAN-KOZENEY RELATION:
! PERM(PHI) = d^2*PHI^3/(180*(1-PHI)^2)


! the permeability is calculated in this way because rho_s{3D} = (2/3)rho_s{2D}
! the coefficient which enters the scaling of fluid equation is alpha*d^2
! where alpha = 1/45/12
	

       

	do ix=1,mx+1	   
	   do iy=1,my+1
	      if (numpar(ix,iy) .ne. 0.)then
		 siterd2(ix,iy)=siterd2(ix,iy)/numpar(ix,iy)		
	      else
		 siterd2(ix,iy) = 0.
	      endif
	   enddo
!           siterd2(ix,2) = 0.5*(siterd2(ix,2)+0.25)
!           siterd2(ix,my+2) = 0.5*(siterd2(ix,my+2)+0.25)
	enddo



	maxsiterd2 = 0.
	do ix=1,mx+1	   
	   do iy=1,my+1
	      if (siterd2(ix,iy) .gt. maxsiterd2) maxsiterd2 = siterd2(ix,iy)
	   enddo
	enddo

	do ix=1,mx+1	   
	   do iy=1,my+1
!              perm(ix,iy) = 0.25*
!     +                ((1+2*phimat(ix,iy))**3)/
!     +		      (solfra(ix,iy)**2)

	      if (siterd2(ix,iy) .gt. 0.) then
		 perm(ix,iy) = 0.001*siterd2(ix,iy)*((1+2*phimat(ix,iy))**3)/(solfra(ix,iy)**2)
	      else
		 perm(ix,iy) = 0.001*maxsiterd2*((1+2*phimat(ix,iy))**3)/(solfra(ix,iy)**2)
!		 perm(ix,iy) = 0.5*((1+2*phimat(ix,iy))**3)/
!     +		      (solfra(ix,iy)**2)
	      endif	      
                 
!	      perm(ix,iy) = (phimat(ix,iy)**3)/(180*(solfra(ix,iy))**2)
	      
	   enddo
	enddo


	do ix=1,mx	   
	   do iy=1,my+1
	      if (numvx(ix,iy) .lt. 0.) then
		 print*,'numvx is: ',numvx(ix,iy), ix, iy
	      endif
	      if (numvx(ix,iy) .gt. 0.)then		 
		 sitevx(ix,iy) = sitevx(ix,iy)/numvx(ix,iy)		 
	      else	
		 sitevx(ix,iy) = 0.	
	      endif
	   enddo
	enddo
	do ix=1,mx+1	   
	   do iy=1,my+2
	      if (numvy(ix,iy) .lt. 0.) then 
		 print*,'numvy is: ',numvy(ix,iy), ix, iy
	      endif
	      if (numvy(ix,iy) .gt. 0.)then		 
		 sitevy(ix,iy) = sitevy(ix,iy)/numvy(ix,iy)		 
	      else	
		 sitevy(ix,iy) = 0.	
	      endif
	   enddo
	enddo

!       * CALCULATING THE DIVERGENCE OF THE VELOCITY *
!       INNER GRID

	do ix=2,mx
	   do iy=1,my+1
	      divvel(ix,iy) = (sitevx(ix,iy) - sitevx(ix-1,iy))/(dxc) + (sitevy(ix,iy+1) - sitevy(ix,iy))/(dyc)
	   enddo
	enddo
!       RIGHT AND LEFT
	
	do iy=1,my+1
	   divvel(1,iy) = (sitevx(1,iy) - sitevx(mx,iy))/(dxc) + (sitevy(1,iy+1) - sitevy(1,iy))/(dyc)
	   divvel(mx+1,iy) = divvel(1,iy)
	enddo



	do ix = 2,mx
	   divvel(ix,1) = (sitevx(ix,1) - sitevx(ix-1,1))/(dxc) + (sitevy(ix,2) - sitevy(ix,1)) /(dyc)
	   divvel(ix,my+1)=(sitevx(ix,my+1) - sitevx(ix-1,my+1))/(dxc) + (sitevy(ix,my+2) - sitevy(ix,my+1))/(dyc)
	enddo


!       CORNERS
	divvel(1,1) = (sitevx(1,1) - sitevx(mx,1))/(dxc) + (sitevy(1,2) - sitevy(1,1))/(dyc)
	divvel(mx+1,1)=divvel(1,1)
	divvel(1,my+1)=(sitevx(1,my+1) - sitevx(mx,my+1))/(dxc) + (sitevy(1,my+2) - sitevy(1,my+1))/(dyc)
	divvel(mx+1,my+1) = divvel(1,my+1)


	
  
! files:
! phipermfile - porosity | permeability
! areasfile - total grain area | weighted radiuses^2 | volume fraction (area/cellsize) | number grains in a cell
! velocityfile - v(1) | v(2) | div(v)
! numvelofile - numvx | numvy

	if ( toprint .gt. 0 ) then
	   open(unit=phi_perm_file,file=TRIM(output_directory)//'/phipermfile'//TRIM(file_postfix)//fnn)
	   open(unit=areas_file,file=TRIM(output_directory)//'/areasfile'//TRIM(file_postfix)//fnn) 
	   open(unit=velocity_file,file=TRIM(output_directory)//'/velocityfile'//TRIM(file_postfix)//fnn)
	   open(unit=num_velocity_file,file=TRIM(output_directory)//'/numvelofile'//TRIM(file_postfix)//fnn)
	   write(phi_perm_file,*) mx+1,mx+1 
	   write(areas_file,*) mx+1,mx+1,mx+1,mx+1
	   write(velocity_file,*) mx+1,mx+1,mx+1
	   write(num_velocity_file,*) mx+1,mx+1
	   write(phi_perm_file,*) my + 1,my + 1
	   write(areas_file,*) my + 1, my + 1,my + 1,my + 1
	   write(velocity_file,*) my + 2, my + 2,my + 2
	   write(num_velocity_file,*) my + 2, my + 2
	   do iy=1,my+1
	      do ix = 1,mx+1
		 write(phi_perm_file,*)  phimat(ix,iy),perm(ix,iy)
		 write(areas_file,*) vc(ix,iy),siterd2(ix,iy),solfra(ix,iy),numpar(ix,iy)
		 write(velocity_file,*) sitevx(ix,iy),sitevy(ix,iy),divvel(ix,iy)
		 write(num_velocity_file,*) numvx(ix,iy),numvy(ix,iy)
	      enddo
	   enddo
	   do ix = 1,mx+1
	      write(velocity_file,*) sitevx(ix,my+2),sitevy(ix,my+2),divvel(ix,my+2)
	      write(num_velocity_file,*) numvx(ix,my+2),numvy(ix,my+2)
	   enddo
	   close(phi_perm_file)
	   close(areas_file) 
	   close(velocity_file)
	endif
	
 104	format(2(E20.15))
 105	format(2(I4))
	return
	end






