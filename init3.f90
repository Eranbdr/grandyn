! ******************************************************************
! ** ROUTINES TO INITIATE VELOCITY AND POSITION VECTORS, **
! ******************************************************************

	SUBROUTINE INITBOUND()

	use mycommons

	real*8 dist,deviate,value,logg


	INTEGER   I

	real*8 xtogo
	real gauss
	real rand1
	real*8 rbb
	integer nt
	real*8 rxg,ryg,rad0
	integer dflag
	
!     ******************************************************************

	xright  = boxx*.5
	xleft   = -xright
	ytop    = boxy*.5
	ybot    = -ytop
	rbb = sigb

! BOUNDARY GRAINS ARE RAND0M-SIZE SINGLE ATOMS            
            
! for simplicity, make them doubly-periodic
!  fill in from left to right and top to bottom
	print*, 'initialize boundaries'
	print*, 'init xright,xleft,ytop,ybot,rbb,boxx,boxy ', xright,xleft,ytop,ybot,rbb,boxx,boxy
	print*, 'bottom: '

	if (ib(1) .ge. 0) then
		!     ** FORM THE BOTTOM BOUNDARY ***
		rxg = xleft
		ryg = ybot
		! find grain size from set of gaussian distributions
		dist=rand1(randum)
		if (dist .le. frac1) then 
			1          deviate=gauss(randum)
			deviate=std1*deviate
			value=mean1 + deviate 
			dflag=1
			if (value .lt. logr .or. value .gt. higr) goto 1
		else 
		2          deviate=gauss(randum)
			deviate=std2*deviate
			value=mean2 + deviate 
			dflag=2
			if (value .lt. logr2 .or. value .gt. higr2) goto 2
		endif
		rad0=rbb*0.5*value
		call addgrain2(rxg,ryg,rad0,1)
		gtype(n)=dflag
		xtogo=xright-2*radius(1)

		DO 200 i = 2, maxn
			dist=rand1(randum)
			if (dist .le. frac1) then 
				dflag=1
				logg=logr
			3          deviate=gauss(randum)
				deviate=std1*deviate
				value=mean1 + deviate 
				if (value .lt. logr .or. value .gt. higr) goto 3
			else 
				dflag=2
				logg=logr2
				4          deviate=gauss(randum)
				deviate=std2*deviate
				value=mean2 + deviate 
				if (value .lt. logr2 .or. value .gt. higr2) goto 4
			endif
			rad0=rbb*0.5*value
			rxg = r(1,i-1)+radius(i-1)+rad0
			
			
			if (2*rad0 .gt. xtogo) then
				if (xtogo .lt. logg) then
					radius(i-1)=radius(i-1)+.5*xtogo
					r(1,i-1) = r(1,i-1)+.5*xtogo
					nbound(1)=i-1
				else
					rad0=xtogo/2.
					rxg = r(1,i-1)+radius(i-1)+rad0
					nbound(1)=i
					call addgrain2(rxg,ryg,rad0,1)
					gtype(n)=dflag
				endif
				goto 201
			else
				call addgrain2(rxg,ryg,rad0,1)
				gtype(n)=dflag
				xtogo = xright-radius(1)-r(1,i)-radius(i)
			endif
	200       CONTINUE
	201       continue


	endif
	   
	print*, 'top: '
	nt=nbound(1)+1
	rxg = xleft
	ryg = ytop
	! find grain size from set of gaussian distributions
	dist=rand1(randum)
	if (dist .le. frac1) then 
		dflag=1
		logg=logr
		5          deviate=gauss(randum)
		deviate=std1*deviate
		value=mean1 + deviate 
		if (value .lt. logr .or. value .gt. higr) goto 5
	else 
		dflag=2
		logg=logr2
		6          deviate=gauss(randum)
		deviate=std2*deviate
		value=mean2 + deviate 
		if (value .lt. logr2 .or. value .gt. higr2) goto 6
	endif
	rad0=rbb*0.5*value
	call addgrain2(rxg,ryg,rad0,2)
	gtype(n)=dflag
	xtogo=xright-2*radius(nt)
	DO 300 i = nt+1, maxn
	dist=rand1(randum)
	if (dist .le. frac1) then 
		7          deviate=gauss(randum)
		deviate=std1*deviate
		value=mean1 + deviate 
		dflag=1
		if (value .lt. logr .or. value .gt. higr) goto 7
	else 
		8          deviate=gauss(randum)
		deviate=std2*deviate
		value=mean2 + deviate 
		dflag=2
		if (value .lt. logr2 .or. value .gt. higr2) goto 8
	endif
	rad0=rbb*0.5*value
	rxg = r(1,i-1)+radius(i-1)+rad0
	if (2*rad0 .gt. xtogo) then
		if (xtogo .lt. logg) then
			radius(i-1)=radius(i-1)+.5*xtogo
			r(1,i-1) = r(1,i-1)+.5*xtogo
			nbound(2)=i-1
		else
			rad0=xtogo/2.
			rxg = r(1,i-1)+radius(i-1)+rad0
			nbound(2)=i
			call addgrain2(rxg,ryg,rad0,2)
			gtype(n)=dflag
		endif
	goto 301
	else
		call addgrain2(rxg,ryg,rad0,2)
		gtype(n)=dflag
		xtogo = xright-radius(nt)-r(1,i)-radius(i)
	endif
	300       CONTINUE
	301       continue

	   do i=1,nbound(2)
             grainarea=grainarea+.5*pi*radius(i)*radius(i)
	   enddo

	nbound(3)=nbound(2)
	nbound(4)=nbound(3)


        RETURN
        END


