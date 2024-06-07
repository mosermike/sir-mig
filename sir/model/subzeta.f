c..............................................................................

c subzeta rutina que evalua la presion 

	subroutine subzeta(ntau,tau,t,pe,kap,pg,ro)

c ntau tau t pe are input 
c kap,pg,ro are output
c I think IT does not use equilibrio hidrostatico as it is said above (vmp)

	implicit real*4 (a-h,o-z)
	include 'PARAMETER'  !por kt

	real*4 tau(*),t(*),pe(*),kap(kt),pg(*),ro(*),pp(10),kac,d1(10),d2,d3
	real wgt,abu,ei1,ei2
        common/constantes/g,avog	!gravedad,n. avogadro/pmu
	common/mu/cth
      
	g=2.7414e+4*cth		!gravedad cm/s^2 en fotosfera solar
	avog=6.023e23
        bol=1.380662e-16      !cte boltzmann cgs

        do i=1,10
           d1(i)=0
        end do
c calculamos el peso molecular medio pmu
	pmusum=0.0
        asum=0.0
	do i=1,92
  	   ii=i
	   call neldatb(ii,0.,wgt,abu,ei1,ei2)
           asum=asum+abu		
	   pmusum=pmusum+wgt*abu
	end do

	do i=1,ntau
           tsi=t(i)
	   psi=pe(i)
           call gasc(tsi,psi,psg,pp)
           pg(i)=psg

	   call kappach(5.e-5,tsi,psi,pp,d1,d1,kac,d2,d3)

c  pp(8) es la abundancia electronica =pe/p(h') 
           pesomedio=pmusum/(asum+pp(8)) !peso molec. medio

           kap(i)=1.0081*kac*avog/pmusum ! abs. coef. per gram
           ro(i)=pg(i)*pesomedio/bol/t(i)/avog 
c  	   print*,i,pg(i),t(i),ro(i),pe(i),kap(i)
	end do
	return
	end


