
	subroutine partgen()
! generate a randomly-filled box of particles, starting in the lower
! left corner and adding particles along the bottom, then starting
! a new row. Put the particles as close together as possible without
! overlapping them. By putting particles closer (reduce the factor
! gfac), reduces initial pore space but there will be overlap.)
 
		use mycommons
		real gauss,rand1
		real*8 yedge(maxn),yedgepos(maxn),ynewpos(maxn),ynewedge(maxn)
		real*8 gfac,rad0,rxg,grad0,ryg,xedge
		real*8 yt,solidvol,xr
		real*8 newarea
		real*8 highest
		real*8 dist,deviate,value
		real*8 averad
		integer nranatom,nlastrow,l,krow,nrow,nrowmax,np
		integer i
		integer dflag
		integer kl0,kl1

		 
		print*, 'generating random grid of grains...'

!  add grains, until their volume reaches the desired porosity
!  at the eventual box size
        solidvol = (1.-phigoal)*boxvol
                
        nlastrow=nbound(1)

! create a surface that drapes over the last row of particles, defined
!  by the points(yedgepos,yedge) in (x,y)
! add new grains on top of this surface, starting at left wall and working
! toward right wall. Then drape a new surface over the top of this and
! start over. 
!  rad0=atom radius in this grain
!  grad0=effective grain radius(some factor*rad0)
!  nlastrow = number of points that make up the underlying surface
!    (given by the number of grains placed into the previous row)
!  krow = counter for number of grains in this row
!  nrow = counter for number of rows
        xr = xright-0.5
        yt = ytop-0.5

        do l=1,nlastrow
			yedge(l)=ybot+radius(l)
			!         yedgepos(l)=xleft+
			!     &      (xright-xleft)*dble(l-1)/(dble(nlastrow)-1.)
			yedgepos(l)=r(1,l)
        enddo
        yedge(nlastrow+1)=yedge(1)
        yedgepos(nlastrow+1)=xr
        xedge=xleft
                
        krow=0
        nrow=1
        nrowmax=150
        kl0=1
		!kl1=1


        do 100 np=1,maxn

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
				if (value .lt. logr2 .or. value .gt. higr2) goto 2
					dflag=2
				endif
				rad0=0.5*value

				nranatom=1
				gfac=1.
				newarea = grainarea+rad0*rad0*nranatom*(pi-alap)
				if (grainarea.gt.solidvol) goto 101

				grad0=gfac*rad0
				rxg=xedge+grad0
				xedge=rxg+grad0          
          
!  make a new row
		if (xedge .gt. xright) then
			209		format(2x,f7.3,2x,f7.3)	      
			nlastrow=krow
			krow=0
			nrow=nrow+1
			kl0=1
			rxg=xleft+grad0
			xedge = rxg+grad0
			do l=1,nlastrow
				yedgepos(l)=ynewpos(l)
				yedge(l)=ynewedge(l)
			enddo
			yedge(nlastrow+1)=yedge(nlastrow)
			yedgepos(nlastrow+1)=xr
			!              print*, nrow-1,nlastrow
			if (nrow .gt. nrowmax) goto 101

		endif

!           loop over lastrow to find height of pile
   99        krow=krow+1
			do l=2,nlastrow+1
				if ((rxg-grad0) .lt. yedgepos(l)) then
					kl0=l-1
					goto 31
				endif
			enddo
   31        continue
		do l=kl0,nlastrow+1
			if ((rxg+grad0) .lt. yedgepos(l)) then
				kl1=l
				goto 32
			endif
		enddo
   32        continue
              highest=yedge(kl0)
              do l=kl0+1,kl1
				highest=max(highest,yedge(l))
              enddo
	      ryg=highest + grad0 
	     
	     ynewedge(krow)=ryg+grad0
	     ynewpos(krow)=rxg
!	     ynewpos(krow)=rxg+grad0

!  if this particle runs into the top boundary, quit adding particles
			if (ynewedge(krow) .gt. (yt)) then
				rxg=rxg+grad0
				xedge=rxg+grad0
				if (xedge .gt. xr) then
					nrowmax = nrow
					print*, 'box too short, increase boxy'
					stop
				else
					goto 99
				endif
			endif
		rxg=rxg+(ryg-ybot)*.001
		call addgrain2(rxg,ryg,rad0,0)
		grainarea=newarea
		gtype(n)=dflag
		emod(n)=e1
	      
             
 100        continue
 101        continue
          print*, 'grainarea, solidvol',grainarea, solidvol
		
		averad=sum(radius)/dble(n)
		print*, "average particle diameter =", 2*averad
 
		return
 	 end 
              

	subroutine addgrain2(rxg,ryg,rad0,isbound)

        use mycommons
        real*8 rxg,ryg,rad0
        integer isbound
        
        n=n+1
		
		if (n .gt. maxn) then
			print*, 'TOO MANY GRAINS!'
			stop
		endif
        
        bdgrain(n) = isbound
	 
! check array size
		if (n+1 .gt. maxn) then
			write(*,*) 'n = ', n
			print*, 'TOO MANY PARTICLES. INCREASE N IN GRANDAM.PAR'
			!          call writcn(tau)
			stop
		endif 
          
! set the positions of the particles
		r(1,n)=rxg
		r(2,n)=ryg
		radius(n) = rad0
		radinv2(n) = .125/(radius(n)*radius(n)*radius(n))
		orientation(n)=0.
		v(1,n)=0.
		v(2,n)=0.
		w(n)=0.
		gtype(n)=1
!	write(1,101) ngrains,'(',numatoms,')'
  100   format(i4,$)
  101   format(i5,a,i1,a,$)

		return
	end

 
	function gauss(idum)
  
		integer idum
		real gauss
		integer iset
		real fac,gset, rsq,v1,v2,rand1
		save iset,gset
		data iset /0/

		if (iset .eq. 0) then
	  1      v1=2*rand1(idum)-1.
			 v2=2*rand1(idum)-1.
			 rsq=v1*v1+v2*v2
			 if (rsq .ge. 1. .or. rsq .eq. 0.) goto 1
			 fac=sqrt(-2*log(rsq)/rsq)
			 gset=v1*fac
			 gauss=v2*fac
			 iset=1
		else
			 gauss=gset
			 iset=0
		endif
        return
	end
              


	subroutine shearstep(gstep)        

		use mycommons
		real*8 dxs,gstep
		integer i

		do i=1,n
			dxs=gstep*(r(2,i)-ybot) 
			r(1,i)=r(1,i)+dxs
		enddo

		return
	end

