c	INTERPOLATODOS interpola todos
c_______________________________________________________________________
c
        subroutine interpolatodos(ngrado,
     &      ntauold,tauold,told,pold,micold,hold,vold,gamold,fiold,
     &      ntaunew,taunew,tnew,pnew,micnew,hnew,vnew,gamnew,finew)

	implicit real*4 (a-h,o-z)
	real*4 tauold(*),told(*),pold(*),micold(*)
        real*4 hold(*),vold(*),gamold(*),fiold(*)
        real*4 taunew(*),tnew(*),pnew(*),micnew(*)
        real*4 hnew(*),vnew(*),gamnew(*),finew(*)
	real*4 xa(11),ya(11)
      
        if(ngrado.gt.10)then
           print*,'el maximo grado del polinomio interpolador es 10'
           print*,'Stop en interpolatp'
           STOP
        end if

        n2=int(ngrado/2)
        do i=1,ntaunew
           taunewi=taunew(i)
           call locate(tauold,ntauold,taunewi,j)

	   n3=j-n2-1

           if(n3.lt.0.)n3=0

           if(n3+ngrado+1.gt.ntauold)n3=ntauold-ngrado-1

           if(n3.lt.0.)n3=0

	   do k=1,ngrado+1
	      xa(k)=tauold(n3+k)
	      ya(k)=told(n3+k)
	   end do
	   CALL POLINT(xa,ya,ngrado+1,taunewi,ynewi,ERROR)
           tnew(i)=ynewi

	   do k=1,ngrado+1
	      ya(k)=alog(pold(n3+k))
	   end do
	   CALL POLINT(xa,ya,ngrado+1,taunewi,ynewi,ERROR)
           pnew(i)=exp(ynewi)

	   do k=1,ngrado+1
	      ya(k)=micold(n3+k)
	   end do
	   CALL POLINT(xa,ya,ngrado+1,taunewi,ynewi,ERROR)
           micnew(i)=ynewi

	   do k=1,ngrado+1
	      ya(k)=hold(n3+k)
	   end do
	   CALL POLINT(xa,ya,ngrado+1,taunewi,ynewi,ERROR)
           hnew(i)=ynewi

	   do k=1,ngrado+1
	      ya(k)=vold(n3+k)
	   end do
	   CALL POLINT(xa,ya,ngrado+1,taunewi,ynewi,ERROR)
           vnew(i)=ynewi

	   do k=1,ngrado+1
	      ya(k)=gamold(n3+k)
	   end do
	   CALL POLINT(xa,ya,ngrado+1,taunewi,ynewi,ERROR)
           gamnew(i)=ynewi

	   do k=1,ngrado+1
	      ya(k)=fiold(n3+k)
	   end do
	   CALL POLINT(xa,ya,ngrado+1,taunewi,ynewi,ERROR)
           finew(i)=ynewi





       end do

       end

   
