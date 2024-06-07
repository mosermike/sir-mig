c	EQUILIBRIUM
c_______________________________________________________________________
c
	implicit real*4 (a-h,o-z)
	parameter (kt=10000)
        real*4 atmos(8*kt)
	real*4 tau(kt),t(kt),p(kt),mic(kt),h(kt),v(kt),gam(kt),fi(kt)   
        real*4 taunew(kt),tnew(kt),pnew(kt),micnew(kt),hnew(kt),vnew(kt)
        real*4 gamnew(kt),finew(kt)
        real*4 kmnew(kt),ronew(kt),pgnew(kt),z1new(kt)
        real*4 km(kt),ro(kt),pg(kt),y(kt),z1(kt),stray
        character*100 modelin,modelout,modelouttau
	character*100 fichabun
	common/mu/cth
        common/preciso/prec      
        common/fichabun/fichabun

	print*,'Abundance file? : '
	read(*,'(a)')fichabun
	print*,'Input model atmosphere:                        :'
	read(*,'(a)')modelin

        print*,'OPTIONS:'
c        print*,'(1) entra en tau y sale en tau en eq. hidrostatico'
c        print*,'(2) entra en tau y sale en z   en eq. hidrostatico'
c        print*,'(3) entra en tau y sale en tau y z en eq. hidrostatico'
c        print*,'(4) entra en tau y sale en z (NO equilibrio)'
c        print*,'(11) entra en tau y sale en tau en eq. dinamico'
c        print*,'(12) entra en tau y sale en z   en eq. dinamico'
c        print*,'(13) entra en tau y sale en tau y z en eq. dinamico'

        print*,'(1) Optical depth --> optical depth in hydrostatic equilibrium'
        print*,'(2) Optical depth --> geometrical height in hydrostatic equilibrium'
        print*,'(4) Optical depth --> geometrical height (no hydrostatic equilibrium)'

        read*,opc

	if(opc.eq.1.or.opc.eq.3.or.opc.eq.11.or.opc.eq.13)then
           print*,'Output model atmosphere: '
	   read(*,'(a)')modelouttau
        end if
 	if(opc.ne.1.and.opc.ne.11)then
	   print*,'Output file (geometrical heights): '
	   read(*,'(a)')modelout
        end if

        print*,'Cosine of theta (mu)='
        read*,cth
        print*,'Spacing in log tau='
        read*,taumax
	
        prec=1.e-8
        tau0=1.
        imin=1
c	taumax=.05
c leemos los modelos
	call leemodi1(0,modelin,vmac,fill,stray,atmos,ntau)
	do i=1,ntau
	   tau(i)=atmos(i)
	   t(i)  =atmos(i+ntau)
	   p(i)  =atmos(i+2*ntau)
	   mic(i)=atmos(i+3*ntau)
	   h(i)  =atmos(i+4*ntau)  !debe contener dv/dt en el caso dinamico
	   v(i)  =atmos(i+5*ntau)
	   gam(i)=atmos(i+6*ntau)
	   fi(i) =atmos(i+7*ntau)
	end do

c los interpolamos a una red de paso    (si es menor que el paso)
        ntaunew=1.+int((tau(1)-tau(ntau))/taumax+.6)
        if(ntaunew.gt.kt)then
           print*,' There are more than ',kt,' points in the new depth grid!! ' 
           stop
        end if
        if(ntaunew.gt.ntau)then
           ngrado=2
           pasomedio=2.30259*taumax/2.
           do i=1,ntaunew
              taunew(i)=tau(1)-(i-1)*taumax
              if(abs(taunew(i)).lt.tau0)then
                 imin=i
                 tau0=taunew(i)
              end if
           end do 

	   call interpolatodos(ngrado,
     &          ntau,tau,t,p,mic,h,v,gam,fi,
     &          ntaunew,taunew,tnew,pnew,micnew,hnew,vnew,gamnew,finew)


c los ponemos en equilibrio hidrostatico
	   if(opc.eq.1.or.opc.eq.2.or.opc.eq.3)call 
     &                equisubmu(ntaunew,taunew,tnew,pnew)	   		
c	   if(opc.eq.11.or.opc.eq.12.or.opc.eq.13)call 
c     &                equisubmu2(ntaunew,taunew,tnew,pnew,hnew)

	   if(opc.eq.1.or.opc.eq.3.or.opc.eq.11.or.opc.eq.13)then
	      do i=1,ntaunew
	         atmos(i)=taunew(i)
	         atmos(i+ntaunew)=tnew(i)
	         atmos(i+2*ntaunew)=pnew(i)
	         atmos(i+3*ntaunew)=micnew(i)
	         atmos(i+4*ntaunew)=hnew(i)
	         atmos(i+5*ntaunew)=vnew(i)
	         atmos(i+6*ntaunew)=gamnew(i)
	         atmos(i+7*ntaunew)=finew(i)
	      end do
              call leemodi1(1,modelouttau,vmac,fill,stray,atmos,ntaunew)
              if(opc.eq.1.or.opc.eq.11)stop
           end if

c calculamos la presion gaseosa y el kappa por gramo
	   call subzeta(ntaunew,taunew,tnew,pnew,kmnew,pgnew,ronew)

           do i=1,ntaunew
              y(i)=10.**(taunew(i)) / (kmnew(i)*ronew(i))
           end do

           z1new(1)=0.
           do i=2,ntaunew
              z1new(i)=z1new(i-1)+pasomedio*(y(i-1)+y(i))
           end do
        
	   do i=1,ntaunew
	      atmos(i)=(z1new(i)-z1new(imin))*1.e-5 
	      atmos(i+ntaunew)=tnew(i)
	      atmos(i+2*ntaunew)=pnew(i)
	      atmos(i+3*ntaunew)=pgnew(i)
	      atmos(i+4*ntaunew)=hnew(i)
	      atmos(i+5*ntaunew)=vnew(i)
	      atmos(i+6*ntaunew)=ronew(i)
	      atmos(i+7*ntaunew)=kmnew(i)
	   end do
	   call leemodi3z(1,modelout,vmac,fill,stray,atmos,ntaunew)

        else
c los ponemos en equilibrio hidrostatico

	   if(opc.eq.1.or.opc.eq.2.or.opc.eq.3)call 
     &                equisubmu(ntau,tau,t,p)	   		
c	   if(opc.eq.11.or.opc.eq.12.or.opc.eq.13)call 
c     &                equisubmu2(ntau,tau,t,p,h)

	   if(opc.eq.1.or.opc.eq.3.or.opc.eq.11.or.opc.eq.13)then

	      do i=1,ntau
	         atmos(i+2*ntau)=p(i)
	      end do
              call leemodi1(1,modelouttau,vmac,fill,stray,atmos,ntau)
              if(opc.eq.1.or.opc.eq.11)stop
           end if

	   call subzeta(ntau,tau,t,p,km,pg,ro)

           do i=1,ntau
c             taue(i)=10.**(tau(i))
c             y(i)=1. / (km(i)*ro(i))
              y(i)=10.**(tau(i)) / (km(i)*ro(i))

              if(abs(tau(i)).lt.tau0)then
                 imin=i
                 tau0=tau(i)
              end if
           end do

           z1(1)=0.
           do i=2,ntau
c              z1(i)=z1(i-1)+(taue(i-1)-taue(i))*(y(i-1)+y(i))/2.
              z1(i)=z1(i-1)+(tau(i-1)-tau(i))*2.30259*(y(i-1)+y(i))/2.
           end do
        
	   do i=1,ntau
	      atmos(i)=(z1(i)-z1(imin))*1.e-5 
	      atmos(i+ntau)=t(i)
	      atmos(i+2*ntau)=p(i)
	      atmos(i+3*ntau)=pg(i)
	      atmos(i+4*ntau)=h(i)
	      atmos(i+5*ntau)=v(i)
	      atmos(i+6*ntau)=ro(i)
	      atmos(i+7*ntau)=km(i)
	   end do
	   call leemodi3z(1,modelout,vmac,fill,stray,atmos,ntau)

        end if

	end

   
