pro nodosfit
;
model=dblarr(8,53)
;
openr,1,'penumjti.mod'
readf,1,a,b,c
readf,1,model
close,1
;
tau=reform(model(0,*))
tem=reform(model(1,*))
;
a=poly_fit(tau,tem,4,/double)





stop
end
