function [Vunbal,Cunbal] = unbal_cal(Va,Vb,Vc,Sa,Sb,Sc)

Ia = abs(Sa)./Va;
Ib = abs(Sb)./Vb;
Ic = abs(Sc)./Vc;

V=[Va Vb Vc];
I=[Ia Ib Ic];

Vavg = (Va+Vb+Vc)./3;
Iavg=(Ia+Ib+Ic)./3;

Vdiff = abs(V-Vavg);
Idiff = abs(I-Iavg);

Vunbal = (max(Vdiff)/Vavg)*100;
Cunbal = (max(Idiff)/Iavg)*100;
