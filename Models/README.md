# Models for SIR

Here is a collection of models


Penumbral models:
	penumjti.mod .. Del Toro Iniesta J.C., Tarbell T.D., Ruiz Cobo B., 
        	        1994, ApJ 436,400 
 			(penumbral model)
	penumjti-inter.mod .. Interpolated penumbral model to 75 depth points
	penumjti-inter-short.mod .. Interpolated penumbral model to 61 depth points
	penumjti-inter-all.mod .. Interpolated penumbral model to 61 depth points with 11 columns
One-component quiet Sun models:
	atmcl1.mod
        atmcl1-B2000.mod .. Quiet sun model with 2 kG and constant in B,vlos and inc. This can be used as a base model.
        atmcl1-B2000-easy.mod .. Quiet sun model with 2 kG and constant in B,vlos and inc without z, rho and Pg. This can be used as a base model.
	atmcl1-5-2500.mod .. Quiet sun model from log tau 0.5 to -2.5
	atmcl1-500-1.mod .. Linear decreasing magnetic field and constant vlos
	model-B2000.mod .. Model with B = 2000, const. in B,vlos,phi,gamma and from log tau 1 to -2.9
	hsra.mod   ...  Harvard Smithsonian Reference Atmosphere (Gingerich O.,Noyes R.W., Kalkofen W., & Cuny Y., 1971. Sol. Phys. 18, 347)
	

Umbra models:
	atmcl10.mod .. Umbra models with constant B, inc, azi, vlos
	atmcl10-linear.mod .. Umbra models with constant vlos and linear decreasing B
	cool.mod .. Collados M., Martínez Pillet V., Ruiz Cobo B., 
                  Del Toro Iniesta J.C., & Vázquez M. 1994 A&A 291 622 
                  (Umbral model for a big spot)
	cool11.mod .. Collados M., Martínez Pillet V., Ruiz Cobo B., 
                  Del Toro Iniesta J.C., & Vázquez M. 1994 A&A 291 622 
                  (Umbral model for a big spot)
Other Models:	
