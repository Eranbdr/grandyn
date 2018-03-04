        subroutine rotatea
        use mycommons
	integer i
        

! do atom rotations
        do 100 i=1,n
          if (gtype(i) .ne. 0) then 
         if (bdgrain(i) .eq. 0) then 
            orientation(i)=orientation(i)+dt*w(i) + dt*dt*tq(i)*radinv2(i)/radius(i)
            orientation(i)=mod(orientation(i),2*pi)
            w(i) = w(i) + dt*(tq(i)*radinv2(i)/radius(i))
         else
          orientation(i)=0.
          w(i)=0.
         endif
         endif
  100   continue
   
		        
  	return
	end


        subroutine rotateb

        use mycommons
	integer i

     

        
       
        do 100 i=1,n
           if (gtype(i) .ne. 0) then 
         if (bdgrain(i) .eq. 0) then 
            w(i) = w(i) + dt*(tq(i)*radinv2(i)/radius(i))
         else
          w(i)=0.
         endif
      endif
  100   continue
!        write(*,*) 'point 300=',orientation(300),w(300),tq(300)

  	return
	end

