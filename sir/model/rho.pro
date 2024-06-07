pro rho
;
close,/all
dat=fltarr(11,55)
openr,1,'hot11.mod'
readf,1,a,b,c
readf,1,dat
close,1
;
r=dat(10,*)
z=dat(8,*)
;
p0=where(abs(z) eq min(abs(z)))
p0=fix(p0(0))
r0=r(p0)
;
rf=fltarr(55,100)
hr=findgen(100)+50.
ch2=fltarr(100)

for j=0,99 do begin
   for k=0,54 do begin
      rf(k,j)=r0*exp(-z(k)/hr(j))
   endfor
   ch2(j)=total((r(p0:p0+10)-rf(p0:p0+10,j))^2.)
endfor
   
v=where(ch2 eq min(ch2))
v=fix(v(0))

set_plot,'ps'
device,filename='rho.ps'

plot,z,r*1e7,charsize=2,ytitle='!7q!5 [10!E7!N g cm!E-3!N]',xtitle='z [km]' $
     ,yrange=[0,10],/xsty,/ysty,xrange=[-150,400],thick=4,xthick=4,ythick=4 $
     ,xticklen=0.05,yticklen=0.05,charthick=3
oplot,z,rf(*,v)*1e7,line=2,thick=4
oplot,[100,160],[7,7],line=0,thick=4
oplot,[100,160],[8,8],line=2,thick=4
xyouts,180,7,'!5Cool umbral model !C (Collados et al. 1994)',CHARTHICK=2
xyouts,180,8,'!5Fit with Density scale !C height of 135 km',CHARTHICK=2

;
zz=z(p0:p0+8)
rr=r(p0:p0+8)
lrr=alog(rr)
cc=poly_fit(zz,lrr,1,yfit,/double)
plot,zz,rr*1e7,charsize=2,ytitle='!7q!5 [10!E7!N g cm!E-3!N]',xtitle='z [km]' $
     ,yrange=[3,6],/xsty,/ysty,xrange=[0,75],thick=4,xthick=4,ythick=4 $
     ,xticklen=0.05,yticklen=0.05,charthick=3
oplot,zz,(2.71829^yfit)*1e7,line=2,thick=4
oplot,[5,15],[3.5,3.5],line=0,thick=4
oplot,[5,15],[4,4],line=2,thick=4
xyouts,17,3.5,'!5Cool umbral model !C (Collados et al. 1994)',CHARTHICK=2
xyouts,17,4,'!5Fit with Density scale !C height of 135 km',CHARTHICK=2
device,/close
set_plot,'x'


stop
end
