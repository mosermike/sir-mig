c LEEMODi1
c iesc=0 :Lee una atmosfera modelo contando puntos en tau
c iesc=1 :escribe
c ntau : numero de puntos en tau
c tau  : log10(tau)
	subroutine leemodi1(iesc,MODEL,VMAC,FILL,stray,atmos,ntau)

	include 'PARAMETER'
	real*4 atmos(*)
	real*4 FILL,vmac,a(8)
        character*(*) model

	IF(IESC.EQ.0)THEN

	open(22,file=model,status='old',err=991)
	ntau=0
	read(22,*)a(1),a(2)
	do while(ntau.lt.kt)
	   ntau=ntau+1
	   read(22,*,end=221)a
	end do
221	ntau=ntau-1
	close(22)
	open(22,file=model)
	read(22,*)VMAC,FILL,stray
	do i=1,ntau
	   read(22,*)(atmos(i+j*ntau),j=0,7)
	end do
	close(22)
	ELSE
	open(22,file=model)
	write(22,*)VMAC,FILL,stray
	do i=1,ntau
	   write(22,100)(SNGL(atmos(i+j*ntau)),j=0,7)
	end do
	close(22)
	END IF

	if(atmos(2).gt.atmos(1))then
	   print*,'da la vuelta al modelo (necesito tau decreciente)'
	   STOP
	end if

100     format(1x,f6.3,1x,f11.3,1x,1pe12.5,1x,5(e12.5,1x))
	return

991	print*,'el fichero del modelo ',model,' no existe'
	stop

	end
