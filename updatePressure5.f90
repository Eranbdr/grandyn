!ccccccccccccccc 15/2/09
!   File should go with SmoothNonExtended
!   Modified from updatePressure3, but with the removal of the extra layers 
!   at the top and bottom of the visible domain.  

!cccccccccccccccc 16/9/09
!   Modified from updatePressure4.f but instead of regular mean for the permeability 
!   between grid points it uses harmonic mean.
! ***************************************************************************************
!       FROM NUMERICAL RECIPES P. 43

	subroutine tridag(a,b,c,rr,u,n)
	integer n,NMAX
	real*8 a(n), b(n), c(n), rr(n), u(n)    
	PARAMETER (NMAX=500)
	real*8 bet,GAM(NMAX)
	integer j

	
	if(b(1).eq.0.) stop 'pb pivot in tridag 1'
	bet=b(1)
	u(1)=rr(1)/bet
	do j=2,n
	   gam(j)=c(j-1)/bet
	   bet=b(j)-a(j)*gam(j)
	   if(bet.eq.0.)then 
	      print*, 'stopped at ',j
	      stop 'pb pivot in tridag'
	   endif
	   u(j)=(rr(j)-a(j)*u(j-1))/bet
	enddo
	do j=n-1,1,-1
	   u(j)=u(j)-gam(j+1)*u(j+1)
	enddo
	return
	end 



! **************************************************************************************
!       FROM NUMERICAL RECIPES P. 68
	subroutine cyclic(a,b,c,alpha,beta,rr,x,n)
	integer n,NMAX
	real*8 alpha, beta, a(n), b(n), c(n), rr(n), x(n)    
	PARAMETER (NMAX=500)
	integer i
	real*8 fact,gamma,bb(NMAX),u(NMAX),z(NMAX)
	if(n .le. 2)stop 'n too small in cyclic'
	if(n .gt. NMAX) stop 'NMAX too samll in cyclic'
	gamma = -b(1)
	bb(1)=b(1)-gamma
	bb(n)=b(n)-alpha*beta/gamma
	do i=2,n-1
	   bb(i) = b(i)
	enddo
	call tridag(a,bb,c,rr,x,n)
	u(1)=gamma
	u(n)=alpha
	do i=2,n-1
	   u(i)=0
	enddo
	call tridag(a,bb,c,u,z,n)
	fact=(x(1)+beta*x(n)/gamma)/(1.+z(1)+beta*z(n)/gamma)
	do i=1,n
	   x(i)=x(i)-fact*z(i)
	enddo
	return
	end


! ************************************************************************************
! dont forget to initialize the pressure matrix *************************************
	subroutine updatepressure(numit)
! When we introduce the fluid and calculate the fluid pressure, we must introduce scaling factors.    
! Note that in the dry granular everything is really non dimensional and only when presenting the result for 
! a physcial system we need to scale everything up.
!
! The equation to solve is: beta*Phi*(dp/dt) + nabla u_s - nabla[kappa/mu nabla P] = 0
! we use: P0 = k/x0
!         V0 = x0*dsqrt(k/m)
!         t0 = dsqrt(m/k)
!         kappa0 = alpha*x0^2, where alpha = 1/(45*12) = 1/540, 
!                  and the permeability is calculated with a radius and not a diamter
!         x0 = d, grain diameter 
! The non dim fluid equation is then:
!  (dp/dt) + (x0/k*beta)*(1/phi)*nabla u_s - dsqrt(m/k)*(alpha/beta*mu)*(1/phi) nabla[kappa nabla P] = 0
!  
!  assign values: beta = 4.5e-10 1/Pa, E (bulk modulus of grains) = 8e10 Pa, x0 = d = 1 mm = 1e-3 m 
!  therefore k = 8e10*1e-3 = 8e7 Pa m
!  m = pi*(4/3)*(d/2)^3*rho_g = 1.3823e-6 kg
!  (rho_g = 2640 kg/m^3)
! Rge resulting coefficients are: (x0/k*beta) = 2.78e-2 and dsqrt(m/k)*(alpha/beta*mu) = 540.94



	use mycommons

                                     
	real*8 acoef
! acoef = 1/E*fluid compessibility
!	parameter (acoef=2.78e-2)
	parameter (acoef=2.777777778e-2)


	real*8 bcoef
! bcoef is Pe^-1 = k/(mu * fluid compessibility *d *velocity)
!	parameter (bcoef= 540.94)
	parameter (bcoef= 540.9411345)


	real*8 a(MAXD+1),b(MAXD+1),c(MAXD+1)
	real*8 rr(MAXD+1),tempP(MAXD+1)
	real*8 halfP(MAXD+1,MAXD+1)
	
      
	integer i,j,debug
	real*8 csqx,csqy,cdiv
	real*8 plhli,plhlj,mnhli,mnhlj,coef2
	real*8 alpha,beta	

	integer bcl
	integer numit
	real*8  locdt
	integer stepprint
	integer locprint
	real*8 wavepress
	integer locdiv
	real*8 maxperm
	real*8 maxdim
	real*8 current
	integer n1,n2,n3,n4
	character a1,a2,a3,a4
	character*4 lfilnum
	
	


	debug = 0
	maxperm = 0
	do i = 1,mx+1
	   do j = 1,my+1
	      if (perm(i,j) .gt. maxperm)maxperm=perm(i,j)
	   enddo
	enddo
	if (dxc .gt. dyc) then
	   maxdim = dxc
	else
	   maxdim = dyc
	endif
	
	      
	
	
	bcl = fbc

	if(bcl .eq. 3 .or. bcl .eq. 4 .or. bcl .eq. 6)bcl=1
	locdt = dt/real(TRATIO)

	current = maxperm*locdt/(maxdim*maxdim)

	
	stepprint = TRATIO

! OUTPUT ONLY IF EXTERNAL TOPRINT IS 1 AND IF STEPPRINT == CURRENTINTERATION.
! NORMALLY STEPPRINT = TRATIO SO ONLY ONE FILES WILL BE PRINTED FOR TRATIO ITERATIONS.
! FOR STEPPRINT = TRATIO/2 WE WILL HAVE TWO PRINTTED FILES FOR TRATIO ITERATION.

	
	locprint = 0

	if (toprint .gt. 0) then
	   if (mod(numit,stepprint) .eq. 0) then
	      locprint = 1
	      write(lfilnum, fmt='(i4.4)') numit/stepprint
	      print*,'current # = ', current
	   endif
	endif
	
	


	

!       first part of ADI update for time n+1/2, for each of the rows.


	csqx = (locdt*bcoef)/(2.*dxc*dxc)
	csqy = (locdt*bcoef)/(2.*dyc*dyc)
	cdiv = (locdt*acoef)/2.
!	print*, dxc,dyc,locdt,cdivmore velo


! %%%%%%%%%%%%%%%%%%%% NEUMANN BOUNDARY CONDITION WITH DERIVATIVE = 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


       	if (bcl .eq. 2) then
	   do i=1,mx
	      do j=1,my+1
		 halfP(i,j)=0.
	      enddo
	      tempP(i) = 0.
	      a(i) = 0.
	      b(i) = 0.
	      c(i) = 0.
	      rr(i) = 0.
	   enddo
	   do j=1,my+1
	      do i = 1,mx
		 coef2 = phimat(i,j)
		 if (j .ne. 1 .and. j .ne. my+1) then
!		    plhlj = (perm(i,j+1)+perm(i,j))/2.
!		    mnhlj = (perm(i,j-1)+perm(i,j))/2.
		    plhlj =2.*perm(i,j+1)*perm(i,j)/(perm(i,j+1)+perm(i,j))
		    mnhlj = 2.*perm(i,j-1)*perm(i,j)/(perm(i,j-1)+perm(i,j))
		    rr(i)= press(i,j-1)*(csqy*mnhlj)+press(i,j)*(coef2-csqy*(plhlj+mnhlj))+press(i,j+1)*(csqy*plhlj) - cdiv*divvel(i,j)		    
		 else 
		    
!       BOTTOM BOUNDARY 
		    if (j .eq. 1) then
!		       plhlj = (perm(i,j+1)+perm(i,j))/2.
		       plhlj = 2.*perm(i,j+1)*perm(i,j)/(perm(i,j+1)+perm(i,j))
		       mnhlj = 0
		       rr(i) = 0. + press(i,j)*(coef2-csqy*plhlj)+press(i,j+1)*(csqy*plhlj) - cdiv*divvel(i,j)
!       TOP BOUNDARY
		    else
		       plhlj = 0.
!		       mnhlj = (perm(i,j-1)+perm(i,j))/2.
		       mnhlj = 2.*perm(i,j-1)*perm(i,j)/(perm(i,j-1)+perm(i,j))
		       rr(i)= press(i,j-1)*(csqy*mnhlj)+press(i,j)*(coef2-csqy*mnhlj)+ 0. - cdiv*divvel(i,j)

		    endif
		 endif
		 
		 if (i .ne. 1 .and. i .ne. mx) then
!		    plhli = (perm(i+1,j)+perm(i,j))/2.
!		    mnhli = (perm(i-1,j)+perm(i,j))/2.	
		    plhli = 2.*perm(i+1,j)*perm(i,j)/(perm(i+1,j)+perm(i,j))
		    mnhli = 2.*perm(i-1,j)*perm(i,j)/(perm(i-1,j)+perm(i,j))
		    a(i) =  - (csqx*mnhli)
		    b(i) = coef2 + csqx*(plhli + mnhli)
		    c(i) =  - (csqx*plhli)
		 else
!       LEFT BOUNDARY	       
		    if (i .eq. 1)then 
!		       plhli = (perm(i+1,j)+perm(i,j))/2.
!		       mnhli = (perm(mx,j)+perm(i,j))/2.
		       plhli = 2.*perm(i+1,j)*perm(i,j)/(perm(i+1,j)+perm(i,j))
		       mnhli = 2.*perm(mx,j)*perm(i,j)/(perm(mx,j)+perm(i,j))
		       a(i) = 0.
		       beta =  - (csqx*mnhli)
		       b(i) = coef2 + csqx*(plhli + mnhli)
		       c(1) =  - (csqx*plhli)
!       RIGHT BOUNDARY
		    else
!		       plhli = (perm(1,j)+perm(i,j))/2.
!		       mnhli = (perm(i-1,j)+perm(i,j))/2.
		       plhli = 2.*perm(1,j)*perm(i,j)/(perm(1,j)+perm(i,j))
		       mnhli = 2.*perm(i-1,j)*perm(i,j)/(perm(i-1,j)+perm(i,j))
		       a(i) = - (csqx*mnhli)
		       b(i) = coef2 + csqx*(plhli + mnhli)
		       c(i) = 0.
		       alpha = - (csqx*plhli)
		    endif
		 endif
	      enddo      
	      if (debug .eq. 1 .and. j .eq. 4) then
		 do i = 1,mx
		    print*, 'i,j = ',i,j
		    print*,'a(i,j) =',a(i)
		    print*,'b(i,j) =',b(i)
		    print*,'c(i,j) =',c(i)
		    print*,'rr(i,j) =',rr(i)
		 enddo
		 print*, 'alpha',alpha
		 print*, 'beta',beta
		 print*,'calling cyclic' 
	      endif
	      
	      call cyclic(a,b,c,alpha,beta,rr,tempP,mx)
	      if (debug .eq. 1 .and. j .eq. 4) print*, 'tempP = ('
	      do i=1,mx
		 if (debug .eq. 1  .and. j .eq. 4) print*,i,',',j,') = ', tempP(i)
		 halfP(i,j)=tempP(i)
	      enddo   
	      halfP(mx+1,j) = halfP(1,j)
	   enddo

!       second part of ADI update for the next time step n+1
	   do j=1,my+1
	      tempP(j) = 0.
	      a(j) = 0.
	      b(j) = 0.
	      c(j) = 0.
	      rr(j) = 0.
	   enddo
	    
	   
	   do i = 1,mx
	      do j=1,my+1
		 coef2 = phimat(i,j)
		 if (i .ne. 1 .and. i .ne.mx) then
!		    plhli = (perm(i+1,j)+perm(i,j))/2.
!		    mnhli = (perm(i-1,j)+perm(i,j))/2.
		    plhli = 2.*perm(i+1,j)*perm(i,j)/(perm(i+1,j)+perm(i,j))
		    mnhli = 2.*perm(i-1,j)*perm(i,j)/(perm(i-1,j)+perm(i,j))
		    rr(j) = halfP(i-1,j)*(csqx* mnhli) + halfP(i,j)*(coef2-csqx*(plhli+mnhli))+halfP(i+1,j)*(csqx* plhli) - cdiv*divvel(i,j)
		 else
		    if (i .eq. 1)then
!		       plhli = (perm(i+1,j)+perm(i,j))/2.
!		       mnhli = (perm(mx,j)+perm(i,j))/2.
		       plhli = 2.*perm(i+1,j)*perm(i,j)/(perm(i+1,j)+perm(i,j))
		       mnhli = 2.*perm(mx,j)*perm(i,j)/(perm(mx,j)+perm(i,j))
		       rr(j) = halfP(mx,j)*( csqx* mnhli) + halfP(i,j)*(coef2-csqx*(plhli+mnhli))+halfP(i+1,j)*(csqx* plhli)-cdiv*divvel(i,j)
		    else
!		       plhli = (perm(1,j)+perm(i,j))/2.
!		       mnhli = (perm(i-1,j)+perm(i,j))/2.
		       plhli = 2.*perm(1,j)*perm(i,j)/(perm(1,j)+perm(i,j))
		       mnhli = 2.*perm(i-1,j)*perm(i,j)/(perm(i-1,j)+perm(i,j))
		       rr(j) = halfP(i-1,j)*( csqx* mnhli) + halfP(i,j)*(coef2-csqx*(plhli+mnhli))+halfP(1,j)*( csqx* plhli)-cdiv*divvel(i,j) 
		    endif
		 endif
		 if(j .ne. 1 .and. j .ne. my+1) then
!		    plhlj = (perm(i,j+1)+perm(i,j))/2.
!		    mnhlj = (perm(i,j-1)+perm(i,j))/2.
		    plhlj = 2.*perm(i,j+1)*perm(i,j)/(perm(i,j+1)+perm(i,j))
     		    mnhlj = 2.*perm(i,j-1)*perm(i,j)/(perm(i,j-1)+perm(i,j))
		    a(j) =  - (csqy*mnhlj)
		    b(j) = coef2 + csqy*(plhlj + mnhlj)
		    c(j) =  - (csqy*plhlj)
		 else
		    if (j .eq. 1)then
!		       plhlj = (perm(i,j+1)+perm(i,j))/2.
                       plhlj = 2.*perm(i,j+1)*perm(i,j)/(perm(i,j+1)+perm(i,j))
		       a(j)= 0.
		       b(j)= coef2 + csqy*plhlj
		       c(j) = - csqy*plhlj
		    else
!		       mnhlj = (perm(i,j-1)+perm(i,j))/2.
                       mnhlj = 2.*perm(i,j-1)*perm(i,j)/(perm(i,j-1)+perm(i,j))
		       a(j) = - csqy*mnhlj
		       b(j) = coef2 + csqy*mnhlj
		       c(j) = 0.
		    endif

		 endif
	      enddo
	    	      
	      call tridag(a,b,c,rr,tempP,my+1)
	      do j=1,my+1
		 press(i,j)=tempP(j)
	      enddo
	   enddo
	   
	   do j = 1,my+1
	      press(mx+1,j)=press(1,j)
	   enddo
       	endif

! %%%%%%%%%%%%%%%% DIRICHLET BOUNDARY CONDITION %%%%%%%%%%%%%%%%%%%%%%%%%%%

	if (fbc .eq. 1 .or. fbc .eq. 3 .or. fbc .eq. 4) then
	   do i=1,mx
	      do j=2,my
		 halfP(i,j)=0.
	      enddo
	      tempP(i) = 0.
	      a(i) = 0.
	      b(i) = 0.
	      c(i) = 0.
	      rr(i) = 0.
	      press(i,1) = pbot
	      press(i,my+1) = ptop
	      halfP(i,1)=press(i,1)
	      halfP(i,my+1)=press(i,my+1)
	   enddo
	   
	   do j=2,my
	      do i = 1,mx
		 coef2 = phimat(i,j)
!		 plhlj = (perm(i,j+1)+perm(i,j))/2.
!		 mnhlj = (perm(i,j-1)+perm(i,j))/2.
                 plhlj = 2.*perm(i,j+1)*perm(i,j)/(perm(i,j+1)+perm(i,j))
                 mnhlj = 2.*perm(i,j-1)*perm(i,j)/(perm(i,j-1)+perm(i,j))
		 if (debug .eq. 1) then
!                   do l = 1,mx
!		      print*, 'i,j = ',l,j
!		      print*,'press(i,j-1) =',press(l,j-1)
!		      print*,'press(i,j) =',press(l,j)
!		      print*,'press(i,j+1) =',press(l,j+1)
!		      enddo
		 endif
		 rr(i) =  press(i,j-1)*(csqy*mnhlj) + press(i,j)*(coef2 - csqy*(plhlj+mnhlj)) + press(i,j+1)*(csqy*plhlj) - cdiv*divvel(i,j)
		 if (i .ne. 1 .and. i .ne. mx) then
!		    plhli = (perm(i+1,j)+perm(i,j))/2.
!		    mnhli = (perm(i-1,j)+perm(i,j))/2.
                    plhli = 2.*perm(i+1,j)*perm(i,j)/(perm(i+1,j)+perm(i,j))
                    mnhli = 2.*perm(i-1,j)*perm(i,j)/(perm(i-1,j)+perm(i,j))
		    a(i) =  - (csqx*mnhli)
		    b(i) = coef2 + csqx*(plhli + mnhli)
		    c(i) =  - (csqx*plhli)
		 else
!       LEFT BOUNDARY	       
		    if (i .eq. 1)then 
!		       plhli = (perm(i+1,j)+perm(i,j))/2.
!		       mnhli = (perm(mx,j)+perm(i,j))/2.
                       plhli = 2.*perm(i+1,j)*perm(i,j)/(perm(i+1,j)+perm(i,j))
		       mnhli = 2.*perm(mx,j)*perm(i,j)/(perm(mx,j)+perm(i,j))
		       a(i) = 0.
		       beta =  - (csqx*mnhli)
		       b(i) = coef2 + csqx*(plhli + mnhli)
		       c(1) =  - (csqx*plhli)
!       RIGHT BOUNDARY
		    else
!		       plhli = (perm(1,j)+perm(i,j))/2.
!		       mnhli = (perm(i-1,j)+perm(i,j))/2.
                       plhli = 2.*perm(1,j)*perm(i,j)/(perm(1,j)+perm(i,j))
		       mnhli = 2.*perm(i-1,j)*perm(i,j)/(perm(i-1,j)+perm(i,j))
		       a(i) = - (csqx*mnhli)
		       b(i) = coef2 + csqx*(plhli + mnhli)
		       c(i) = 0.
		       alpha = - (csqx*plhli)
		    endif
		 endif
	      enddo      

	    
	      call cyclic(a,b,c,alpha,beta,rr,tempP,mx)
	      do i=1,mx		 
		 halfP(i,j)=tempP(i)
	      enddo
	      halfP(mx+1,j)=halfP(1,j)        
	   enddo

	   do j=1,my+1
	      tempP(j) = 0.
	      a(j) = 0.
	      b(j) = 0.
	      c(j) = 0.
	      rr(j) = 0.
	   enddo	
	   do i = 1,mx 
	      do j=2,my
		 coef2 = phimat(i,j)
		 if (i .ne. 1 .and. i .ne.mx) then
!		    plhli = (perm(i+1,j)+perm(i,j))/2.
!		    mnhli = (perm(i-1,j)+perm(i,j))/2.
                    plhli = 2.*perm(i+1,j)*perm(i,j)/(perm(i+1,j)+perm(i,j))
		    mnhli = 2.*perm(i-1,j)*perm(i,j)/(perm(i-1,j)+perm(i,j))
		    rr(j) = halfP(i-1,j)*(csqx* mnhli) + halfP(i,j)*(coef2-csqx*(plhli+mnhli))+halfP(i+1,j)*( csqx* plhli) - cdiv*divvel(i,j)
		 else
		    if (i .eq. 1)then
!		       plhli = (perm(i+1,j)+perm(i,j))/2.
!		       mnhli = (perm(mx,j)+perm(i,j))/2.
		       plhli = 2.*perm(i+1,j)*perm(i,j)/(perm(i+1,j)+perm(i,j))
		       mnhli = 2.*perm(mx,j)*perm(i,j)/(perm(mx,j)+perm(i,j))
		       rr(j) = halfP(mx,j)*( csqx* mnhli)+halfP(i,j)*(coef2 - csqx*(plhli+mnhli))+halfP(i+1,j)*( csqx* plhli)-cdiv*divvel(i,j)
		    else
!		       plhli = (perm(1,j)+perm(i,j))/2.
!		       mnhli = (perm(i-1,j)+perm(i,j))/2.
		       plhli = 2.*perm(1,j)*perm(i,j)/(perm(1,j)+perm(i,j))
		       mnhli = 2.*perm(i-1,j)*perm(i,j)/(perm(i-1,j)+perm(i,j))
		       rr(j) = halfP(i-1,j)*( csqx* mnhli)+halfP(i,j)*(coef2-csqx*(plhli+mnhli))+halfP(1,j)*( csqx* plhli)-cdiv*divvel(i,j) 
		    endif
		 endif
		 
!		 plhlj = (perm(i,j+1)+perm(i,j))/2.
!		 mnhlj = (perm(i,j-1)+perm(i,j))/2.
		 plhlj = 2.*perm(i,j+1)*perm(i,j)/(perm(i,j+1)+perm(i,j))
		 mnhlj = 2.*perm(i,j-1)*perm(i,j)/(perm(i,j-1)+perm(i,j))
		 a(j) =  - (csqy*mnhlj)
		 b(j) = coef2 + csqy*(plhlj + mnhlj)
		 c(j) =  - (csqy*plhlj)
	      enddo
	      rr(1) = halfP(i,1)
	      a(1) = 0
	      b(1) = 1
	      c(1) = 0
	      rr(my+1) = halfP(i,my+1)
	      a(my+1) = 0
	      b(my+1) = 1
	      c(my+1) = 0
	      call tridag(a,b,c,rr,tempP,my+1)
	      do j=1,my+1
		 press(i,j)=tempP(j)
	      enddo
	   enddo
	   
	   do j = 1,my+1
	      press(mx+1,j)=press(1,j)
	   enddo
	endif





! %%%%%%%%%%%%%%%% MIXED B.C. DIRICHLET TOP AND NEUMANN BOTTOM %%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%  at the bottom 0 pressure gradient and at the top 0 pressure %%%%%%%%%%%%%%%%%%%%

	if (fbc .eq. 5) then
	   do i=1,mx
	      do j=1,my
		 halfP(i,j)=0.
	      enddo
	      tempP(i) = 0.
	      a(i) = 0.
	      b(i) = 0.
	      c(i) = 0.
	      rr(i) = 0.
	      press(i,my+1) = 0.
	      halfP(i,my+1)=press(i,my+1)
	   enddo
	   do j=1,my
	      do i = 1,mx
		 coef2 = phimat(i,j)
		 if (j .ne. 1) then
!		    plhlj = (perm(i,j+1)+perm(i,j))/2.
!		    mnhlj = (perm(i,j-1)+perm(i,j))/2.
		    plhlj = 2.*perm(i,j+1)*perm(i,j)/(perm(i,j+1)+perm(i,j))
		    mnhlj = 2.*perm(i,j-1)*perm(i,j)/(perm(i,j-1)+perm(i,j))

		    rr(i)= press(i,j-1)*(csqy*mnhlj)+press(i,j)*(coef2-csqy*(plhlj+mnhlj))+press(i,j+1)*(csqy*plhlj) - cdiv*divvel(i,j)		    
		 else 		    
!       BOTTOM BOUNDARY 		    
!		    plhlj = (perm(i,j+1)+perm(i,j))/2.
                    plhlj = 2.*perm(i,j+1)*perm(i,j)/(perm(i,j+1)+perm(i,j))
		    mnhlj = 0
		    rr(i) = 0. + press(i,j)*(coef2-csqy*plhlj)+press(i,j+1)*(csqy*plhlj) - cdiv*divvel(i,j)
		 endif		 		 
		 if (i .ne. 1 .and. i .ne. mx) then
!		    plhli = (perm(i+1,j)+perm(i,j))/2.
!		    mnhli = (perm(i-1,j)+perm(i,j))/2.	    
		    plhli = 2.*perm(i+1,j)*perm(i,j)/(perm(i+1,j)+perm(i,j))
		    mnhli = 2.*perm(i-1,j)*perm(i,j)/(perm(i-1,j)+perm(i,j))	    

		    a(i) =  - (csqx*mnhli)
		    b(i) = coef2 + csqx*(plhli + mnhli)
		    c(i) =  - (csqx*plhli)
		 else
!       LEFT BOUNDARY	       
		    if (i .eq. 1)then 
!		       plhli = (perm(i+1,j)+perm(i,j))/2.
!		       mnhli = (perm(mx,j)+perm(i,j))/2.
		       plhli = 2.*perm(i+1,j)*perm(i,j)/(perm(i+1,j)+perm(i,j))
		       mnhli = 2.*perm(mx,j)*perm(i,j)/(perm(mx,j)+perm(i,j))
		       a(i) = 0.
		       beta =  - (csqx*mnhli)
		       b(i) = coef2 + csqx*(plhli + mnhli)
		       c(1) =  - (csqx*plhli)
!       RIGHT BOUNDARY
		    else
!		       plhli = (perm(1,j)+perm(i,j))/2.
!		       mnhli = (perm(i-1,j)+perm(i,j))/2.
		       plhli = 2.*perm(1,j)*perm(i,j)/(perm(1,j)+perm(i,j))
		       mnhli = 2.*perm(i-1,j)*perm(i,j)/(perm(i-1,j)+perm(i,j))
		       a(i) = - (csqx*mnhli)
		       b(i) = coef2 + csqx*(plhli + mnhli)
		       c(i) = 0.
		       alpha = - (csqx*plhli)
		    endif
		 endif
	      enddo      
	      
	      
	      call cyclic(a,b,c,alpha,beta,rr,tempP,mx)
	      if (debug .eq. 1 .and. j .eq. 4) print*, 'tempP = ('
	      do i=1,mx
		 if (debug .eq. 1  .and. j .eq. 4) print*,i,',',j,') = ', tempP(i)
		 halfP(i,j)=tempP(i)
	      enddo   
	      halfP(mx+1,j) = halfP(1,j)
	   enddo

!       second part of AFI update for the next time step n+1
	   do j=1,my+1
	      tempP(j) = 0.
	      a(j) = 0.
	      b(j) = 0.
	      c(j) = 0.
	      rr(j) = 0.
	   enddo
	    
	   
	   do i = 1,mx
	      do j=1,my
		 coef2 = phimat(i,j)
		 if (i .ne. 1 .and. i .ne.mx) then
!		    plhli = (perm(i+1,j)+perm(i,j))/2.
!		    mnhli = (perm(i-1,j)+perm(i,j))/2.
		    plhli = 2.*perm(i+1,j)*perm(i,j)/(perm(i+1,j)+perm(i,j))
		    mnhli = 2.*perm(i-1,j)*perm(i,j)/(perm(i-1,j)+perm(i,j))
		    rr(j) = halfP(i-1,j)*(csqx* mnhli) + halfP(i,j)*(coef2-csqx*(plhli+mnhli))+halfP(i+1,j)*(csqx* plhli) - cdiv*divvel(i,j)
		 else
		    if (i .eq. 1)then
!		       plhli = (perm(i+1,j)+perm(i,j))/2.
!		       mnhli = (perm(mx,j)+perm(i,j))/2.
		       plhli = 2.*perm(i+1,j)*perm(i,j)/(perm(i+1,j)+perm(i,j))
		       mnhli = 2.*perm(mx,j)*perm(i,j)/(perm(mx,j)+perm(i,j))
		       rr(j) = halfP(mx,j)*( csqx* mnhli) + halfP(i,j)*(coef2-csqx*(plhli+mnhli))+halfP(i+1,j)*(csqx* plhli)-cdiv*divvel(i,j)
		    else
!		       plhli = (perm(1,j)+perm(i,j))/2.
!		       mnhli = (perm(i-1,j)+perm(i,j))/2.
                       plhli = 2.*perm(1,j)*perm(i,j)/(perm(1,j)+perm(i,j))
		       mnhli = 2.*perm(i-1,j)*perm(i,j)/(perm(i-1,j)+perm(i,j))
		       rr(j) = halfP(i-1,j)*( csqx* mnhli) + halfP(i,j)*(coef2-csqx*(plhli+mnhli))+halfP(1,j)*( csqx* plhli)-cdiv*divvel(i,j) 
		    endif
		 endif
		 if(j .ne. 1) then
!		    plhlj = (perm(i,j+1)+perm(i,j))/2.
!		    mnhlj = (perm(i,j-1)+perm(i,j))/2.
                    plhlj = 2.*perm(i,j+1)*perm(i,j)/(perm(i,j+1)+perm(i,j))
		    mnhlj = 2.*perm(i,j-1)*perm(i,j)/(perm(i,j-1)+perm(i,j))
		    a(j) =  - (csqy*mnhlj)
		    b(j) = coef2 + csqy*(plhlj + mnhlj)
		    c(j) =  - (csqy*plhlj)
		 else
! bottom boundary
!		    plhlj = (perm(i,j+1)+perm(i,j))/2.
                    plhlj = 2.*perm(i,j+1)*perm(i,j)/(perm(i,j+1)+perm(i,j))
		    a(j)= 0.
		    b(j)= coef2 + csqy*plhlj
		    c(j) = - csqy*plhlj		   
		 endif
	      enddo
	      rr(my+1) = halfP(i,my+1)
	      a(my+1) = 0
	      b(my+1) = 1
	      c(my+1) = 0
	      call tridag(a,b,c,rr,tempP,my+1)
	      do j=1,my+1
		 press(i,j)=tempP(j)
	      enddo
	   enddo
	   
	   do j = 1,my+1
	      press(mx+1,j)=press(1,j)
	   enddo
       	endif

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	 
	do j = 1,my+1
	   do i = 2,mx
	      gradpx(i,j) = (press(i+1,j) - press(i-1,j))/(2.*dxc)
	   enddo
	   gradpx(1,j) = (press(2,j) - press(mx,j))/(2.*dxc)
	   gradpx(mx+1,j) = gradpx(1,j)
	enddo
	
	do i = 1,mx + 1
	   do j = 2,my
	      gradpy(i,j) = (press(i,j+1)-press(i,j-1))/(2.*dyc)
	   enddo
!       neumann = zero flux b.c.
	   if (bcl .eq. 2) then
              gradpy(i,1) = (press(i,2) - press(i,1))/(dyc)
              gradpy(i,my+1)=(press(i,my+1)-press(i,my))/(dyc)	     
	   else 
	      if (bcl .eq. 1) then
!       dirichlet
		 gradpy(i,1) =  (press(i,2) - press(i,1))/(dyc)
		 gradpy(i,my+1) = (press(i,my+1) - press(i,my))/(dyc)
	      else
!       mixed
		 gradpy(i,1) = (press(i,2) - press(i,1))/(dyc)
		 gradpy(i,my+1) = (press(i,my+1) - press(i,my))/(dyc)
	      endif
	   endif	   
	enddo


	   
	      


	      
	
	
! pressfile: halfP | pressure 
! pgradfile:  gradpx | gradpy	
	
	if ( locprint .gt. 0 ) then
	   if (TRATIO .ne. 1) then
	      open(unit=pressure_file,file=TRIM(output_directory)//'/pressfile'//TRIM(file_postfix)//fnn//'.'//lfilnum)
	      open(unit=pressure_grad_file,file=TRIM(output_directory)//'/pgradfile'//TRIM(file_postfix)//fnn//'.'//lfilnum)
	   else
	      open(unit=pressure_file,file=TRIM(output_directory)//'/pressfile'//TRIM(file_postfix)//fnn)
	      open(unit=pressure_grad_file,file=TRIM(output_directory)//'/pgradfile'//TRIM(file_postfix)//fnn)
	   endif
	   write(pressure_file,*) mx,mx 
	   write(pressure_file,*) my + 1, my + 1
	   write(pressure_grad_file,*) mx, mx
	   write(pressure_grad_file,*) my + 1, my + 1
	   do j=1,my+1
	      do i = 1,mx
		 write(pressure_file,*) halfP(i,j),press(i,j)
		 write(pressure_grad_file,*) gradpx(i,j),gradpy(i,j)
	      enddo
	   enddo
	   close(pressure_file)
	   close(pressure_grad_file)
	endif
	
 104	format(2(E20.15))
 105	format(2(I4))
	return
	end

