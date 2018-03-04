   
	subroutine getneighbors()
 
	use mycommons

!*********************************


        
		INTEGER     ICELL, JCELL0, JCELL, L, NABOR
		integer ncell
		integer ixcell,iycell,ixcell2,ipf,iycell2
		integer i,j,k
		real*8 temp(maxk)
		integer kindex(maxk),kk,kindex2,kntold,kmin,kmax,kspan
!    *******************************************************************
		
		ncell=mx*my

!  
! save the shear forces on the contacts
		do k=1,contactknt
			kindex(k)=(contacti(k)-1)*n+contactj(k)
			temp(k)=contft(k)
			!         print*, contacti(k),contactj(k)
		enddo
		kntold=contactknt


!  zero out the neighbor lists
		do k=1,contactknt
			contacti(k)=0
			contactj(k)=0
			contft(k)=0.
		enddo

	first_active_layer(1:maxn) = .false.
	
! **LOOP OVER ALL GRAINS IN THE CELL **
		icell=0
		contactknt=0
		DO 8000 iycell=1,my
			do 8000 ixcell=1,mx
				icell=icell+1 
				K=HEAD(ICELL)
				7800  IF (K .GT. 0) THEN

				!  ** LOOP OVER ALL GRAINS BELOW K IN THE CURRENT CELL **
					L=LIST(K)
				7700    IF (L .GT. 0) THEN
					if ( (bdgrain(k).ne.bdgrain(l)).or.(bdgrain(k)*bdgrain(l).eq.0) ) then
						if (gtype(k)*gtype(l) .ne. 0) then
							call setneighbor(k,l,0)
						endif
					endif
					L = LIST(L)
					GO TO 7700
				ENDIF			!end 7700 if 

				!  ** LOOP OVER NEIGHBOURING CELLS **
				JCELL0 = 4 * (ICELL - 1)
				DO 7600 NABOR = 1, 4
				JCELL = MAP ( JCELL0 + NABOR )
				ipf=0
				if (mx .eq. 2) then 
					ipf=2
				else
					if (ixcell .eq. 1) then 
						iycell2=(jcell-1)/mx
						ixcell2=jcell-(iycell2)*mx
						if (ixcell2 .eq. mx) ipf=1
					else if (ixcell .eq. mx) then
						iycell2=(jcell-1)/mx
						ixcell2=jcell-(iycell2)*mx
						if (ixcell2 .eq. 1) ipf=1
					endif           
				endif 
				!           print*, k,l,icell,ixcell,jcell,ixcell2,ipf          

				!     **  LOOP OVER ALL GRAINS IN NEIGHBOURING CELLS **
				if (JCELL .eq. 0) goto 7600
				L = HEAD(JCELL)
				7650         IF ( L .NE. 0 ) THEN

					if ( (bdgrain(k).ne.bdgrain(l)).or.(bdgrain(k)*bdgrain(l).eq.0) ) then
						if (gtype(k)*gtype(l) .ne. 0) then
							call setneighbor(k,l,ipf)
						endif
					endif
				L = LIST(L)
				GO TO 7650
				ENDIF   		!end 7650 if

			7600       CONTINUE  		!end NABOR loop

			K = LIST(K)
			GO TO 7800

			ENDIF     		!end 7800 if

		8000  CONTINUE   		!end ICELL loop

!  now reset the shear force array
!	do k=1,contactknt
!	 contft(k)=temp(contacti(k),contactj(k))
!	enddo

!        write(69,*) 'old # = ',kntold,' new # = ',contactknt
	kspan=3
	do 2 k=1,contactknt
! index hasn't changed
	 kindex2=(contacti(k)-1)*n+contactj(k)
           if (kindex2 .eq. kindex(k)) then
            contft(k)=temp(k)
!            write(69,*) contacti(k),contactj(k),'  old k =',k,
!     1  '  new k =',k
            goto 2
           endif
! index is nearby in the list
         kmin=max(1,k-kspan)
         kmax=min(kntold,k+kspan)
         do kk=kmin,kmax
           if (kindex2 .eq. kindex(kk)) then
            contft(k)=temp(kk)
!            write(69,*) contacti(k),contactj(k),'  old k =',kk,
!     1  '  new k =',k
            goto 2
           endif
         enddo
! contact has moved a lot
         do kk=kmin-1,1,-1
           if (kindex2 .eq. kindex(kk)) then
            contft(k)=temp(kk)
!            write(69,*) contacti(k),contactj(k),'  old k =',kk,
!     1  '  new k =',k
            goto 2
           endif
         enddo
         do kk=kmax+1,kntold
           if (kindex2 .eq. kindex(kk)) then
            contft(k)=temp(kk)
!            write(69,*) contacti(k),contactj(k),'  old k =',kk,
!     1  '  new k =',k
            goto 2
           endif
         enddo
   2    continue
        


!	nbknt=0
!	print*, 'number of distances checked=',nknt
!	do i=1,n
!         nbknt=nbknt+numnbr(i)
!	enddo
!	print*, 'number of neighbors=',nbknt
!	do i=1,n
!	 do k=1,numnbr(i)
!	  j=neighbor(i,k)
!	  write(num_velocity_file,*) i,j
!         enddo
!        enddo
        
        RETURN
        END


        subroutine setneighbor(i,j,ipf)
!**************************************************************
!   check the neighbor list, add contact if particles are
!   within distint of each other

        use mycommons
		real*8 rijsq,rxij2
		real*8 xl,ryij,rxij,rmaxsq,rmax
		integer i,imx,imn,j, ipf
		xl=xright-xleft
		imn=min(i,j)
		imx=max(i,j)

		rmax = radius(imn)+radius(imx)+distint
		rmaxsq=rmax*rmax
		
		ryij  = r(2,imx) - r(2,imn)
		rxij  = r(1,imx) - r(1,imn)
		if (ipf .eq. 1) then
			rxij=rxij-dsign(xl,rxij)
		else if (ipf .eq. 2) then
			rxij2=rxij-dsign(xl,rxij)
			if (dabs(rxij2) .lt. dabs(rxij)) rxij=rxij2
		endif
		rijsq = (rxij*rxij + ryij*ryij)

		!	      write(77,*) i,j,rmax,dsqrt(rijsq)
		! particles are interacting
		if (rijsq .lt. rmaxsq) then
		!                neighbor(imn,numnbr(imn)+1)=imx
		!                numnbr(imn)=numnbr(imn)+1
			contactknt=contactknt+1
			contacti(contactknt)=imn
			contactj(contactknt)=imx
			if ((bdgrain(imn) .eq. 1) .and. (bdgrain(imx) .eq. 0)) then
				first_active_layer(imx) = .true.
			endif
		endif              



		7   format(i4,i4,6(1x,e12.5))


		return
	end
