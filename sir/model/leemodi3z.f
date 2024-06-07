

c leemodi3z
c iesc=0 :lee una atmosfera modelo contando puntos en tau
c iesc=1 :escribe
c ntau : numero de puntos en tau
c tau  : log10(tau)
	subroutine leemodi3z(iesc,model,vmac,fill,stray,atmos,ntau)

	include 'PARAMETER'
	real*4 atmos(*)
	real*4 fill,vmac,stray,a(8)
        character*20 model

	if(iesc.eq.0)then

	open(22,file=model)
	ntau=0
	read(22,*)a(1),a(2)
	do while(ntau.lt.kt)
	   ntau=ntau+1
	   read(22,*,end=221)a
	end do
221	ntau=ntau-1
        print*,'tengo',ntau,'puntos'
	close(22)
	open(22,file=model)
	read(22,*)vmac,fill,stray
	do i=1,ntau
	   read(22,*)(atmos(i+j*ntau),j=0,7)
	end do
	close(22)
	else
	open(22,file=model)
	write(22,*)vmac,fill,stray
	do i=1,ntau
	   write(22,100)(atmos(i+j*ntau),j=0,7)
	end do
	close(22)
	end if


100     format(1x,f11.4,1x,f11.3,1x,1pe12.5,1x,e12.5,1x,4(e12.5,1x))
	return
	end
