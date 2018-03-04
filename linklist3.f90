!    *******************************************************************
!    ** CONSTRUCTION OF CELL LINKED-LISTS AND USE IN FORCE ROUTINE.   **
!    **                                                               **
!    ** REFERENCES:                                                   **
!    **                                                               **
!    ** QUENTREC AND BROT, J. COMPUT. PHYS. 13, 430, 1975.            **
!    ** HOCKNEY AND EASTWOOD, COMPUTER SIMULATION USING PARTICLES,    **
!    **    MCGRAW HILL, 1981.                                         **
!    **                                                               **
!    ** ROUTINES SUPPLIED:                                            **
!    **                                                               **
!    ** SUBROUTINE MAPS                                               **
!    **    SETS UP MAP OF CELL STRUCTURE FOR USE IN FORCE             **
!    ** SUBROUTINE LINKS ( RCUT )                                     **
!    **    SETS UP HEAD OF CHAIN ARRAY AND LINKED LIST                **
!    ** SUBROUTINE FORCE ( SIGMA, RCUT, V, W )                        **
!    **    CALCULATES FORCES USING A LINKED LIST                      **
!    **                                                               **
!    ** USAGE:                                                        **
!    **                                                               **
!    ** SUBROUTINE MAPS IS CALLED ONCE AT THE START OF A SIMULATION   **
!    ** TO ESTABLISH CELL NEIGHBOUR IDENTITIES.  AT EACH TIMESTEP,    **
!    ** SUBROUTINE LINKS IS CALLED TO SET UP THE LINKED LIST AND THIS **
!    ** IS IMMEDIATELY USED BY SUBROUTINE FORCE.                      **
!    *******************************************************************



        SUBROUTINE MAPS

        use mycommons
	
!    *******************************************************************
!    ** ROUTINE TO SET UP A LIST OF NEIGHBOURING CELLS                **
!    **                                                               **
!    ** PRINCIPAL VARIABLES:                                          **
!    **                                                               **
!    ** INTEGER M                  NUMBER OF CELLS IN EACH DIRECTION  **
!    ** INTEGER MAPSIZ             SIZE OF CELL-CELL MAP              **
!    ** INTEGER MAP(MAPSIZ)        LIST OF NEIGHBOURING CELLS         **
!    **                                                               **
!    ** USAGE:                                                        **
!    **                                                               **
!    ** THIS SUBROUTINE SETS UP A LIST OF THE THIRTEEN NEIGHBOURING   **
!    ** CELLS OF EACH OF THE SMALL CELLS IN THE CENTRAL BOX. THE      **
!    ** EFFECTS OF THE PERIODIC BOUNDARY CONDITIONS ARE INCLUDED.     **
!    ** THE SUBROUTINE IS CALLED ONCE AT THE BEGINNING OF THE         **
!    ** SIMULATION AND THE MAP IS USED IN THE FORCE SUBROUTINE        **
!    *******************************************************************

        INTEGER     IX, IY, IMAP, ICELL

!    *******************************************************************

!    ** STATEMENT FUNCTION TO GIVE CELL INDEX **

        ICELL ( IX, IY) = 1 + MOD ( IX - 1 + mx, mx ) + MOD ( IY - 1 + my, my ) * mx

!    ** FIND HALF THE NEAREST NEIGHBOURS OF EACH CELL **


           DO 40 IY = 1, my-1
             if (mx .eq. 1) then
                 IMAP = ( ICELL ( 1, IY ) - 1 ) * 4
                 MAP( IMAP + 1 ) = 0
                 MAP( IMAP + 2 ) = 0
                 MAP( IMAP + 3 ) = ICELL( 1    , IY + 1)
                 MAP( IMAP + 4 ) = 0
             else if (mx .eq. 2) then
                 IMAP = ( ICELL ( 1, IY ) - 1 ) * 4
                 MAP( IMAP + 1 ) = ICELL( 2, IY    )
                 MAP( IMAP + 2 ) = ICELL( 2, IY + 1)
                 MAP( IMAP + 3 ) = ICELL( 1    , IY + 1)
                 MAP( IMAP + 4 ) = 0

                 IMAP = ( ICELL ( 2, IY ) - 1 ) * 4
                 MAP( IMAP + 1 ) = 0
                 MAP( IMAP + 2 ) = ICELL(1     , IY + 1)
                 MAP( IMAP + 3 ) = ICELL( 2    , IY + 1)
                 MAP( IMAP + 4 ) = 0
             else

              DO 30 IX = 1, mx-1

                 IMAP = ( ICELL ( IX, IY ) - 1 ) * 4
                 MAP( IMAP + 1 ) = ICELL( IX + 1, IY    )
                 MAP( IMAP + 2 ) = ICELL( IX + 1, IY + 1)
                 MAP( IMAP + 3 ) = ICELL( IX    , IY + 1)
                 MAP( IMAP + 4 ) = ICELL( IX - 1, IY + 1)

30            CONTINUE
                 IMAP = ( ICELL ( mx, IY ) - 1 ) * 4
                 MAP( IMAP + 1 ) = ICELL( 1, IY    )
                 MAP( IMAP + 2 ) = ICELL( 1, IY + 1)
                 MAP( IMAP + 3 ) = ICELL( mx    , IY + 1)
                 MAP( IMAP + 4 ) = ICELL( mx - 1, IY + 1)

             endif

40         CONTINUE


!  top boundary
	iy=my
             if (mx .eq. 1) then
                 IMAP = ( ICELL ( 1, IY ) - 1 ) * 4
                 MAP( IMAP + 1 ) = 0
                 MAP( IMAP + 2 ) = 0
                 MAP( IMAP + 3 ) = 0
                 MAP( IMAP + 4 ) = 0
             else if (mx .eq. 2) then
                 IMAP = ( ICELL ( 1, IY ) - 1 ) * 4
                 MAP( IMAP + 1 ) = ICELL( 2, IY    )
                 MAP( IMAP + 2 ) = 0
                 MAP( IMAP + 3 ) = 0
                 MAP( IMAP + 4 ) = 0

                 IMAP = ( ICELL ( 2, IY ) - 1 ) * 4
                 MAP( IMAP + 1 ) = 0
                 MAP( IMAP + 2 ) = 0
                 MAP( IMAP + 3 ) = 0
                 MAP( IMAP + 4 ) = 0
             else

              DO 50 IX = 1, mx-1

                 IMAP = ( ICELL ( IX, IY ) - 1 ) * 4
                 MAP( IMAP + 1 ) = ICELL( IX + 1, IY    )
                 MAP( IMAP + 2 ) = 0
                 MAP( IMAP + 3 ) = 0
                 MAP( IMAP + 4 ) = 0

50            CONTINUE
                 IMAP = ( ICELL ( mx, IY ) - 1 ) * 4
                 MAP( IMAP + 1 ) = ICELL( 1, IY    )
                 MAP( IMAP + 2 ) = 0
                 MAP( IMAP + 3 ) = 0
                 MAP( IMAP + 4 ) = 0
             endif

        RETURN
        END



        SUBROUTINE LINKS 

        use mycommons

!    *******************************************************************
!    ** ROUTINE TO SET UP LINKED LIST AND THE HEAD OF CHAIN ARRAYS    **
!    **                                                               **
!    ** PRINCIPAL VARIABLES:                                          **
!    **                                                               **
!    ** INTEGER N                  NUMBER OF ATOMS                    **
!    ** INTEGER M                  NUMBER OF CELLS IN EACH DIRECTION  **
!    ** INTEGER NCELL              TOTAL NUMBER OF CELLS (M**3)       **
!    ** INTEGER LIST(NGRAINS+NBOUND)LINKED LIST OF ATOMS              **
!    ** INTEGER HEAD(NCELL)        HEAD OF CHAIN FOR EACH CELL        **
!    ** REAL*8  RX(N),RY(N),RZ(N)  POSITIONS                          **
!    ** REAL*8  RCUT               THE CUTOFF DISTANCE FOR THE FORCE  **
!    **                                                               **
!    ** USAGE:                                                        **
!    **                                                               **
!    ** EACH GRAIN IS SORTED INTO ONE OF THE MX*MY SMALL CELLS.       **
!    ** THE FIRST GRAIN IN EACH CELL IS PLACED IN THE HEAD ARRAY.     **
!    ** SUBSEQUENT GRAINS ARE PLACED IN THE LINKED LIST ARRAY.        **
!    ** THEN EACH OF THE BOUNDARY ATOMS IS PLACED IN THE LISTS        **
!    ** THE ROUTINE IS CALLED EVERY TIMESTEP BEFORE THE FORCE ROUTINE.**
!    *******************************************************************


	integer k,i,icell
	real*8 dxi,dyi,bx,by
	integer getcell

!    *******************************************************************
	bx=(xright-xleft)
	by=(ytop-ybot)
    dxi = dble(mx)/bx
    dyi = dble(my)/by
        

!    ** ZERO HEAD OF CHAIN ARRAY **
        DO 10 ICELL = 1, mx*my

           HEAD(ICELL) = 0

10      CONTINUE

!    ** SORT ALL GRAINS **  
           ICELL = getcell(r(1,1),r(2,1),bdgrain(1),dxi,dyi)
           LIST(1)     = HEAD(ICELL)
           HEAD(ICELL) = 1
        
        DO 20 i = 2, n
        
          if (bdgrain(i) .lt. 10) then
 
           ICELL = getcell(r(1,i),r(2,i),bdgrain(i),dxi,dyi)
		   !if (icell .lt. 1 .or. icell .gt. maxcell) then
			!write(*,*) 'i, icell, bx, by, dxi, dyi, mx, my v(2,i) ', i, icell, bx, by, dxi, dyi, mx, my , v(2,i)
           !write(*,*) 'xleft,xright,ybot,ytop,r(1,i),r(2,i),bdgrain(i) ',  xleft,xright,ybot,ytop,r(1,i),r(2,i),bdgrain(i)
		   !endif
           LIST(i)     = HEAD(ICELL)
           HEAD(ICELL) = i
!           print*, i,bdgrain(i),icell
           
          endif

20      CONTINUE
        


        RETURN
        END

           function getcell(x,y,ibd,dxi,dyi)
! **************************************************************
        use mycommons
	integer ibd
	integer getcell	
	real*8 dxi,dyi,x,y

           if (ibd .eq. 2) then
             getcell = 1 + INT ( ( x - xleft) * dxi ) + (my-1) * mx
             if (x .ge. xright) getcell = mx*my
           
           else if (ibd .eq. 4) then
             getcell = mx *(1+ INT ( ( y - ybot) * dyi ) )

           else
             getcell = 1 + INT ( ( x - xleft) * dxi ) + INT ( ( y - ybot ) * dyi ) * mx
     
           endif

        RETURN 
        END



                 subroutine dumplinklist
!  for debugging, dump a readable copy of the linklist
!    ***dont' want to do this every time step!***

!  output looks like 
! 	icell =   20
!  	grains:   10_b  11_b  14  16
! 	 neighbor(  21):  13_b  12_b  17
! 	 neighbor(  61):  22  24
! 	 neighbor(  60):
! 	 neighbor(  59):
!  so cell no. 20 contains the interior grains 14 and 16, and the
!   boundary grains 10 and 11
!  the 4 neighbor cells (on the row above, or to the immediate right
!    of cell 20) are also listed.
!  neighboring cells 21 and 61 also contain grains 
                 
      use mycommons
      integer k,l,ncell,icell,ivb,jcell,jcell0,nabor

      
      ivb=1
      if (ivb .eq. 1) then
      open(unit=link_list_dump_file,file=TRIM(output_directory)//'/celllist'//TRIM(file_postfix)//fnn)
      NCELL = mx*my
      do 4999 icell =1,ncell
        k=head(icell)
101     if (k .gt. 0) then
          write(link_list_dump_file,47) 'icell = ',icell
          if (bdgrain(k) .eq. 0) then
            write(link_list_dump_file,48) '  grains: ',k
           else
            write(link_list_dump_file,44) '  grains: ',k,'_b'
          endif
          l=list(k)
103       if (l .gt. 0) then
            if (bdgrain(l) .eq. 0) then
              write(link_list_dump_file,49) l 
            else
              write(link_list_dump_file,46) l,'_b'
            endif  
            l=list(l)
            goto 103
          endif
          
          jcell0 = 4*(icell-1)
          do 105 nabor =1,4
            jcell = map(jcell0+nabor)
            if (jcell .eq. 0) goto 105
            l=head(jcell)
            write(link_list_dump_file,43) '  neighbor(',jcell,'):'
 108        if (l .ne. 0) then
               if (bdgrain(l) .eq. 0) then
                 write(link_list_dump_file,49) l 
               else
                 write(link_list_dump_file,46) l,'_b'
               endif  
               l=list(l)
               goto 108
            endif
 105      continue
        endif
        
 4999 continue
  43  format(/a,i6,a,X,$)
  44  format(a,i6,a,X,$)
  46  format(i6,a,X,$)
  47  format(//a,X,i6)
  48  format(a,i6,X,$)
  49  format(i6,X,$)
  
      close(link_list_dump_file)
      endif
      
!      stop

      return
      end
      
! ******************************************************************

                 subroutine dumplinklist2()
!  for debugging, dump a readable copy of the linklist
!    ***dont' want to do this every time step!***

!  output looks like 
! 	icell =   20
!  	grains:   10_b  11_b  14  16
! 	 neighbor(  21):  13_b  12_b  17
! 	 neighbor(  61):  22  24
! 	 neighbor(  60):
! 	 neighbor(  59):
!  so cell no. 20 contains the interior grains 14 and 16, and the
!   boundary grains 10 and 11
!  the 4 neighbor cells (on the row above, or to the immediate right
!    of cell 20) are also listed.
!  neighboring cells 21 and 61 also contain grains 
                 
      use mycommons
      integer k,l,ncell,icell,ivb,jcell,jcell0,nabor
      
      ivb=1
      if (ivb .eq. 1) then
      open(unit=link_list_dump_file,file=TRIM(output_directory)//'/celllist'//TRIM(file_postfix)//fnn)
      NCELL = mx*my
      do 4999 icell =1,ncell
        k=head(icell)
101     if (k .gt. 0) then
          write(link_list_dump_file,47) 'icell = ',icell
          if (bdgrain(k) .eq. 0) then
            write(link_list_dump_file,48) '  grains: ',k
           else
            write(link_list_dump_file,44) '  grains: ',k,'_b'
          endif
          l=list(k)
103       if (l .gt. 0) then
            if (bdgrain(l) .eq. 0) then
              write(link_list_dump_file,49) l 
            else
              write(link_list_dump_file,46) l,'_b'
            endif  
            l=list(l)
            goto 103
          endif
          
          jcell0 = 4*(icell-1)
          do 105 nabor =1,4
            jcell = map(jcell0+nabor)
            if (jcell .eq. 0) goto 105
            l=head(jcell)
            write(link_list_dump_file,43) '  neighbor(',jcell,'):'
 108        if (l .ne. 0) then
               if (bdgrain(l) .eq. 0) then
                 write(link_list_dump_file,49) l 
               else
                 write(link_list_dump_file,46) l,'_b'
               endif  
               l=list(l)
               goto 108
            endif
 105      continue
        endif
        
 4999 continue
  43  format(/a,i6,a,X,$)
  44  format(a,i6,a,X,$)
  46  format(i6,a,X,$)
  47  format(//a,i6,X)
  48  format(a,i6,X,$)
  49  format(i6,X,$)
  
      close(link_list_dump_file)
      endif
      
!      stop

      return
      end
