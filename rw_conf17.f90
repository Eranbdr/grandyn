

	SUBROUTINE read_restart() 

		use mycommons

		integer i,k
		real*8 tauadd
		character*7 ident
		write(*,*) 'opening restart file ', TRIM(restart_file_name)
		OPEN ( UNIT = restart_file, FILE = TRIM(restart_file_name),FORM = 'UNFORMATTED' )
		read(restart_file) ident
		print*, ident
		read(restart_file) n,(nbound(i),i=1,4)
		read(restart_file) ybot,ytop,xleft,xright
		read(restart_file) tau,tauadd
		write(*,*) 'n,(nbound(i),i=1,4):',n,(nbound(i),i=1,4)
		write(*,*) 'boundaries:',ybot,ytop,xleft,xright
		write(*,*) 'tau',tau,tauadd
		if (ident .eq. 'flforce') then
			read(restart_file) (r(1,i),r(2,i),radius(i),orientation(i),v(1,i),v(2,i),w(i),f(1,i),f(2,i),tq(i),bdgrain(i),i=1,n)
		else
			read(restart_file) (r(1,i),r(2,i),radius(i),orientation(i),v(1,i),v(2,i),w(i),f(1,i),f(2,i),tq(i),bdgrain(i),gtype(i),i=1,n)
		endif
		read(restart_file) (color(i),i=1,n)

		read(restart_file) contactknt
		print*, contactknt
		do k=1,contactknt
			read(restart_file) contacti(k),contactj(k),contfn(k),contft(k)
		enddo

		do k = 1,MAXD+1 
			read(restart_file) (press(i,k),i=1,MAXD+1)
			read(restart_file) (gradpx(i,k),i=1,MAXD+1)
			read(restart_file) (gradpy(i,k),i=1,MAXD+1)
		enddo
		if (ib(2) .eq. 5) read(restart_file,end=10) dspring,dwall,fspring
		10	CLOSE ( restart_file )

		RETURN
	END


        
		
		
		SUBROUTINE write_restart()

        use mycommons
        real*8 tauadd
        character*7 ident
		
        
        integer i,k

!     ******************************************************************
		write(*,*) 'Writing restart file to ' , TRIM(output_directory)//'/restart'//TRIM(file_postfix)//fnn
		

        OPEN ( UNIT = restart_file, FILE = TRIM(output_directory)//'/restart'//TRIM(file_postfix)//fnn,FORM = 'UNFORMATTED' )
        OPEN ( UNIT = contact_info_file, FILE = TRIM(output_directory)//'/contactinfo'//TRIM(file_postfix)//fnn//'.csv')
		
     
          ident="gtype"
          write(restart_file) ident
	  write(restart_file) n,(nbound(i),i=1,4)
	  write(restart_file) ybot,ytop,xleft,xright
	  write(restart_file) tau,tauadd
      write(restart_file) (r(1,i),r(2,i),radius(i),orientation(i),v(1,i),v(2,i),w(i),f(1,i),f(2,i),tq(i),bdgrain(i),gtype(i),i=1,n)
	  
      write(restart_file) (color(i),i=1,n)
	  write(restart_file) contactknt
	  
	  !write(contact_info_file,'(a)') 'i , j, normal, tangential'
	  do k=1,contactknt

	    write(restart_file) contacti(k),contactj(k),contfn(k),contft(k)
	    write(contact_info_file,'(i6,",",i6,",",E21.15,",",E21.15)') contacti(k),contactj(k),contfn(k),contft(k)

	  enddo
        
	  do k = 1,MAXD+1 
		 write(restart_file) (press(i,k),i=1,MAXD+1)
		 write(restart_file) (gradpx(i,k),i=1,MAXD+1)
		 write(restart_file) (gradpy(i,k),i=1,MAXD+1)
	  enddo
	  if (ib(2) .eq. 5) write(restart_file) dspring,dwall,fspring
		
        CLOSE ( UNIT = restart_file )
        CLOSE (contact_info_file)

 951    format(i6,i6,5(e12.4,1x),1x,i1)


		write(*,*) 'done with restart'
		

        RETURN
        END


	subroutine pscriptplot(id_prefix)
!*******************************************************************
!*******************************************************************
!   DRAW A POSTSCRIPT PAGE !!!!
!*******************************************************************
!*******************************************************************

	use mycommons

	real*8  rd(256,2),bl(256,2),gr(256,2)
	save rd,bl,gr
 	integer logging, id_prefix
	character*6 id_prefix_str
    
	real*8 rshade(30),bshade(30),gshade(30)

	integer i,j,k,index
	real*8 gs,red,blue,green
	character*256 com1
	character*128 pspage
	character*128 gfile 
	character*5 xsize
	character*5 ysize
	
	real*8 p1
	integer icnt
	real*8 olap(MAXN*10),maxlap,minlap,wscale,width
	real*8 maxpr,pr(maxn)
	real*8 xl,tforce, height_, pr_i
        
	red = 0.
	blue = 0.
	green = 0.
	write(id_prefix_str,fmt='(i6.6)') id_prefix
!    ****************************************************************
! set up 2 color maps for different populations
!	if (bl(1,1) .eq. 0.) then
!	  do i=1,128
!	    rd(i,2) = 1.+.0*dble(i-1)/127.
!	    bl(i,2) = .7-.7*dble(i-1)/127.
!            gr(i,2)=1.-1.*dble(i-1)/127.
!	    bl(i,1) = 0.9+.0*dble(i-1)/127.
!	    rd(i,1) = .7-.6*dble(i-1)/127.
!            gr(i,1)=1.-0.9*dble(i-1)/127.
!	  enddo
!	endif

! set up red-blue color map 
	if (bl(1,1) .eq. 0.) then
	  do i=1,128
	    rd(i,1) = .3+.7*dble(i-1)/127.
	    bl(i,1) = .3+.7*(127.-dble(i-1))/127.
	    rd(i,2) = .3+.7*dble(i-1)/127.
	    bl(i,2) = .3+.7*(127.-dble(i-1))/127.
	  enddo
	endif
	

    pspage=TRIM(output_directory)//'/'//id_prefix_str//'_bpage'//TRIM(file_postfix)//fnn
		
! start new postscript page and plot box  
	open(unit=psunit,file=pspage)

	if (ib(2) .eq. -2 ) then
		call newpsdoc(internal_top,0,tau)
		if (id_prefix .eq. 0) print*, 'height = ',internal_top-ybot
		write(ysize, '(i5)') int(internal_top-ybot) * 10
	 else
		call newpsdoc(ytop,0,tau)
		print*, 'height = ',ytop-ybot
		if (id_prefix .eq. 0) write(ysize, '(i5)') int(ytop-ybot) * 20
     endif
    
	    
	xl = xright-xleft
	write(xsize, '(i5)') int(xl) * 20
 

!   plot pressures and stress directions
             
!  loop through all the contacts, find the overlaps

	  maxlap=0.
	  minlap=0.
	  maxpr=0.
	  pr=0.
	  !OPEN ( UNIT = contact_info_zero_file, FILE = TRIM(output_directory)//'/contactinfo_zero'//TRIM(file_postfix)//fnn//'.csv')
		!OPEN ( UNIT = contact_info_mid_file, FILE = TRIM(output_directory)//'/contactinfo_mid'//TRIM(file_postfix)//fnn//'.csv')
		!OPEN ( UNIT = contact_info_high_file, FILE = TRIM(output_directory)//'/contactinfo_high'//TRIM(file_postfix)//fnn//'.csv')
	  do k=1,contactknt
	     i=contacti(k)
	     j=contactj(k)
		 pr(i)=pr(i)-contfn(k)
		 pr(j)=pr(j)-contfn(k)
	     tforce=contfn(k)*contfn(k)+contft(k)*contft(k)
		 olap(k)=dsqrt(tforce)
	
		!if (contfn(k) .eq. 0) then
			!write(contact_info_zero_file,'(i6,",",i6,",",E21.15,",",E21.15)') contacti(k),contactj(k),contfn(k),contft(k)
		!else if (abs(contfn(k)) .lt. (internal_top - r(2,i)) * abs(g(2))) then
			!write(contact_info_mid_file,'(i6,",",i6,",",E21.15,",",E21.15)') contacti(k),contactj(k),contfn(k),contft(k)
		!else
			!write(contact_info_high_file,'(i6,",",i6,",",E21.15,",",E21.15)') contacti(k),contactj(k),contfn(k),contft(k)
		!endif
!  FOR LINES REPRESENTING FRICTION
!             if (contfn(k) .eq. 0.) then
!               olap(k)=0.
!             else
!               olap(k)=dabs(contft(k)/contfn(k))
!             endif

		 maxlap=max(maxlap,olap(k))
	  enddo
	  
	 ! close(contact_info_zero_file)
!		close(contact_info_mid_file)
		!close(contact_info_high_file)
		
		do i=1,n
		 
	!                pr(i)=max(0.,f(1,i))
			maxpr=max(maxpr,pr(i))
	!                maxpr=max(maxpr,pr(i))
		enddo
          
!            print*, maxlap , maxpr
!            do k=1,contactknt
!               print*, olap(k)/maxlap
!               enddo
!            do i=1,n
!               print*, pr(i)/maxpr
!               enddo
		
		do 605 i=1,n
		   if (gtype(i) .ne. 0) then
		   	!!	"normal" coloring
			  !if (maxpr .gt. 0.) then
				! p1=pr(i)/maxpr
			  !else
				! p1=0.
			  !endif

			  !if (p1 .gt. 1.) then
				!red=1.
				!blue=.6
				!green=.6
!	 !           print*, 'grain ',i, ' is pink'
!			  else
!				index =1+ max(0,nint(127*p1))
!				red = rd(index,gtype(i))
!				blue = bl(index,gtype(i))
!	! 		     green = gr(index,gtype(i))
!				gs=blue-1.8*(1.-blue)
!				green = max(0.0d0,gs)
!!	            print*, 'grain ',i, ' is ', index,gtype(i),red,green,blue
	!		  endif

	!! "chain" coloring
			height_ = internal_top - r(2,i)
			pr_i = abs(pr(i))
			if (id_prefix .eq. 0) then 
				if (pr_i .gt. height_ * 3. *abs(g(2))) then
				! less normal force than over burden
					red = 1.
					blue = 0.3
					green = 0.3
					
				else
					blue = 1- pr_i / (height_ * 4. *abs(g(2)))
					red = 1 - pr_i / (height_ * 4. *abs(g(2)))
					green = 1 - pr_i / (height_ * 4. *abs(g(2)))
					
				endif
			else
				blue = 1
				red = 1
				green = 1
				if (v(1,i) > 3e-5) then
					blue = 0
					red = 1
					green = 0.2
				else if (v(1,i) > 2e-5) then
					blue = 3 - v(1,i) * 1e5
					red = 1
					green = 0.3
				else if (v(1,i) > 1e-5) then
				! blue
					blue = 1
					red = 0.3
					green = 2 - v(1,i) * 1e5
				else if (v(1,i) > 0) then
				! yel(low) to green
					blue = 0
					red = 1 - v(1,i) * 1e5
					green = 1
				else 
					blue = 1
					red = 1
					green = 1
				endif
			endif
				
			
			!print*, 'grain ',i, ' is ', pr(i),g(2),red,green,blue
			write(psunit,28)  radius(i),r(1,i),r(2,i),red,green,blue,' tbgrain'

!            write(psunit,*) 'gsave'
!            dr=0.5*radius(i)
!            write(psunit,22) r(1,i)-dr,r(2,i)-dr,' moveto'
!            write(psunit,*) '0.05 0.05 scale 0 0 0 setrgbcolor'
!            write(psunit,*) '/Helvetica findfont 8 scalefont setfont'     
!	    write(psunit,contact_info_file) '(',i,') show'
!            write(psunit,*) 'grestore'

			if (ib(3) .eq. -1) then
!      duplicate atoms wrapping around x direction
				if (r(1,i) .gt. xright-radius(i)) then
					write(psunit,28)  radius(i),r(1,i)-xl,r(2,i),red,green,blue,' tbgrain'
				else if (r(1,i) .lt. xleft+radius(i)) then
					write(psunit,28)  radius(i),r(1,i)+xl,r(2,i),red,green,blue,' tbgrain'
				endif
			endif

		   endif
 605     continue
          icnt=0
! draw the network of contacts
          write(psunit,*) '0 0 0 setrgbcolor' 

           wscale=0.25

! set logging=1 for run with gravity, takes the log of
! overlaps when drawing force lines
           logging=0
           if (logging .eq. 1) then
             minlap=dlog10(minlap)
             maxlap=dlog10(maxlap)
           endif


 	    width=1.
          do k=1,contactknt
           if (olap(k) .gt. 0.) then
           
             icnt=icnt+1
	     i=contacti(k)
	     j=contactj(k)

              
	   if (maxlap .ne. minlap) then
            if (logging .eq. 1) then
              width=dlog10(olap(k))/maxlap
            else
               width=olap(k)/maxlap
            endif
           endif
           
           if (width .lt. 0.) width=0.
              write(psunit,29) wscale*width,' setlinewidth'

	      if (dabs(r(1,j)-r(1,i)) .gt. .5*xl) then

	       if (r(1,i) .lt. 0.) then
   	         write(psunit,21) r(1,i),r(2,i),' moveto'
	         write(psunit,22) r(1,j)-xl,r(2,j),' lineto stroke'
  	         write(psunit,21) r(1,i)+xl,r(2,i),' moveto'
	         write(psunit,22) r(1,j),r(2,j),' lineto stroke'
	        else if (r(1,i) .gt. 0.) then
  	         write(psunit,21) r(1,i)-xl,r(2,i),' moveto'
	         write(psunit,22) r(1,j),r(2,j),' lineto stroke'
  	         write(psunit,21) r(1,i),r(2,i),' moveto'
	         write(psunit,22) r(1,j)+xl,r(2,j),' lineto stroke'
	        endif

	      else

	       write(psunit,21) r(1,i),r(2,i),' moveto'
	       write(psunit,22) r(1,j),r(2,j),' lineto stroke'
	      endif
	      	      
	   endif
	  enddo

		


           if (id_prefix .eq. 0) print*,' number of interactions / contacts   = ',icnt, contactknt
   21	   format(e12.4,1x,e12.4,a,$)
   22	   format(e12.4,1x,e12.4,a)
   23	   format(a,e12.4,1x,e12.4,a)
   25      format(7e12.4/6e12.4/6e12.4,a)
   26      format(7e12.4,a)
   27      format(7e12.4/2e12.4,a)
   28      format(6e12.4,a)
   29      format(e12.4,a)
   31      format(a,i3,a)


 
 

        call endpspage()
        close(psunit)
		
		if (igif) then
          
!         gfile = ' trigrainxxx.gif'
!         gfile(10:12) = pspage(7:9)
			gfile = ' > '//TRIM(output_directory)//'/jpegjunk &'
			com1='./pstojpeg.dws.sh -xsize '//xsize//' -rmorig '//TRIM(pspage)//TRIM(gfile)
			if (id_prefix .eq. 0)  print*, com1
			call system(com1)
			

        endif
        RETURN
        END
        
