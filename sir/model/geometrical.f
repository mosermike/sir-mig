c	zeta
c_______________________________________________________________________
c
	implicit real*4 (a-h,o-z)


	include 'PARAMETER'   !por kt
        real*4 tau(kt),t(kt),p(kt),atmos(8*kt),taue(kt)
        real*4 km(kt),ro(kt),pg(kt),y(kt),error(kt),z1(kt),stray
        character*20 modelin,modelout
	character*100 fichabun
	common/mu/cth
        common/constantes/g,avog	!gravedad,n. avogadro/pmu
        common/fichabun/fichabun


	print*,'Abundance file? : '
	read(*,'(a)')fichabun

	print*,'Input model atmosphere? : '
	read(*,'(a)')modelin
	print*,'Output model atmosphere?: '
	read(*,'(a)')modelout
	cth=1.

c leemos los modelos
	tau0=1.
	call leemodi3z(0,modelin,vmac,fill,stray,atmos,ntau)
	do i=1,ntau
	   tau(i)=atmos(i)
	   t(i)=atmos(i+ntau)
	   p(i)=atmos(i+2*ntau)
           if(abs(tau(i)).lt.tau0)then
              imin=i
              tau0=tau(i)
           end if
	end do

c calculamos la presion gaseosa y el kappa por gramo
	call subzeta(ntau,tau,t,p,km,pg,ro)
c        delta=(tau(1)-tau(2 ))*alog(10.)

c calculamos ro (cgs) 
        do i=1,ntau
c           y(i)=1.0/km(i)/ro(i) !para integrar con aformal2i
c           y(i)=(10.**(tau(i))) / (km(i)*ro(i))
            taue(i)=10.**(tau(i))
            y(i)=1. / (km(i)*ro(i))

        end do


c integramos z (la rutina calcula entre i y ntau)
c        call aformal2i(ntau,tau,y,z1,error)
        z1(1)=0.
        do i=2,ntau
c           z1(i)=z1(i-1)+alog(10.)*(tau(i-1)-tau(i))*(y(i)+y(i-1))/2.
            z1(i)=z1(i-1)+(taue(i-1)-taue(i))*(y(i-1)+y(i))/2.
c           print*,'z1(',i,')=',z1(i)
        end do
        
	do i=1,ntau
	   atmos(i)=(z1(i)-z1(imin))*1.e-5 
c           print*,'altura(',i,')=',atmos(i),pg(i),p(i),ro(i),km(i),y(i) 
	   atmos(i+3*ntau)=pg(i)
	   atmos(i+6*ntau)=ro(i)
	   atmos(i+7*ntau)=km(i)
	end do
	call leemodi3z(1,modelout,vmac,fill,stray,atmos,ntau)

	end
c..............................................................................

