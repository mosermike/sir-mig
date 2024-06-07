c	optical.f pasa un fichero en zetas a taus
c se sobreentiende que el fichero tiene una fila de cabecera (vmac,ff)
c y al menos 4 columnas:
c     la primera es z en km (positiva hacia arriba). 
c     la segunda es t en k 
c     la cuarta  es pg en dinas/cm^2
c_______________________________________________________________________
c
	implicit real*4 (a-h,o-z)

	include 'PARAMETER'
        real*4 tau(kt),t(kt),p(kt),f(kt),atmos(8*kt),h(kt),vz(kt)
        real*4 pg(kt),z(kt),d2,d3,kac,kap(kt),ro(kt)
        real*4 wgt,abu,ei1,ei2
        character*20 modelin,modelout
	real tsi,psi,psg,pp(10),d1(10)
	character*100 fichabun

        common/fichabun/fichabun

        common/preciso/prec

        prec=1.e-6
        bol=1.380662e-16      !cte Boltzmann cgs
	avog=6.023e23
	tausup=1.
	

	print*,'Abundance file? : '
	read(*,'(a)')fichabun

	print*,'Input model atmosphere? : '
	read(*,'(a)')modelin
	print*,'Output model atmosphere?: '
	read(*,'(a)')modelout

c Calculamos el peso molecular medio pmu
	pmu=0.0
        asum=0.0
	do i=1,92
  	   ii=i
	   call neldatb(ii,0.,wgt,abu,ei1,ei2)
	   pmu=pmu+wgt*abu
           asum=asum+abu
	end do

c leemos los modelos
        open(1,file=modelin,err=999)
	i=1
	read(1,*)vmac,fill,stray
        do while(i.lt.kt)
           read(1,*,end=911)z(i),t(i),p(i),pg(i),h(i),vz(i)
	   z(i)=z(i)*1.e5	!z se lee en km y se pasa a cm
           i=i+1
        end do
	
911	ntau=i-1 


        print*,'Give the logarithm of the optical depth at the surface:'
        read*,tausup
	tausup=10.**(tausup)	!prof. optica en la superficie
        i0=ntau

	if(p(1).lt.1.)then	!si no conocemos la presion electronica 
             print*,'Electron pressure is being computed from gas pressure' 
             pe0=1.e4	!inicializacion
             do i=1,ntau
           	if(i.gt.1)then
             	   p(i)=p(i-1)
           	else
              	   p(1)=pe0
                end if
                call pefrompg1(t(i),pg(i),p(i))
             end do
	     pnew=p(1)
	     itera=0
	     corr=1.
             do while(corr.gt.1.e-3.and.itera.lt.50) 
		itera=itera+1
		do i=1,ntau
                   call pefrompg1(t(i),pg(i),p(i))
		enddo
	   	corr=abs(p(1)-pnew)/pnew
		pnew=p(1)
             end do
	     if(itera.ge.50)print*,'Precision lower than ',corr*100,'% in the calculation of electron pressures'
	else
c             print*,'se supone conocida la presion electronica'
	end if

	do i=1,ntau
	   tsi=t(i) ! temperature
	   psi=p(i) ! electronic pressure
	   call gasc(tsi,psi,psg,pp)
	   pesomedio=pmu/(asum+pp(8))
           ro(i)=pesomedio*pg(i)/bol/t(i)/avog
	   call kappach(5.e-5,tsi,psi,pp,d1,d1,kac,d2,d3)
	   kap(i)=1.0081*kac*avog*ro(i)/pmu	!por unidad de longitud (cm2/cm3)
	   				        ! kac is abs. coef. per H particle
						! kap is per cm^-1
	end do

c integramos
	f(ntau)=0.
	do i=ntau-1,1,-1
          aa=(kap(i)+kap(i+1))*(z(i)-z(i+1))/2.
          f(i)=f(i+1)+aa
c	  print*,i,f(i)
        end do
           
c	call integral(ntau,z,kap,f) !f(i) es la integral entre 1 e i
c	print*,i0,f(i0),f(i0)+tausup

c	do i=1,ntau
c	   bb=f(i0)+tausup-f(i)
c           print*,'tau(sup)+f(sup)-f(',i,')=',bb
c        end do

	do i=1,ntau
	   bb=f(i0)+tausup-f(i)
	   if (bb.le.0.) then
	      print*,'tau(',i,')=',bb,' cannot be nor 0 neither negative!!!!!'
	      stop
	   end if
           tau(i)=alog10(bb)
           atmos(i)=tau(i)
           atmos(i+ntau)=t(i)
           atmos(i+2*ntau)=p(i)
           atmos(i+3*ntau)=0.
           atmos(i+4*ntau)=h(i)
           atmos(i+5*ntau)=vz(i)

           do j=6,7
              atmos(i+j*ntau)=0.
	   end do

	enddo

        call leemodi1(1,modelout,vmac,fill,stray,atmos,ntau)

	stop
999	print*,'The file ',modelin,' does not exist'	

	end

***********************************************************************

c integral: calcula la integral 
c usa la rutina "integran" entre x(1) y x(j),(d01gaf de la nag)
c dir(j) es la integral entre 1 y j

	subroutine integral(n,x,y,dir)

	implicit real*4 (a-h,o-z)
	real*4 x(*),y(*),dir(*)

	call integran(x,y,n,dir,error,ifail)

	if(ifail.ne.0)then
	   print*,'ifail en integral =',ifail
	   stop
	end if

	dir(1)=0.0
	do i=2,4
	     xpp=x(i)-x(i-1)
	     dir(i)=dir(i-1)+xpp*(y(i)+y(i-1))/2.0
	end do

	return
	end


