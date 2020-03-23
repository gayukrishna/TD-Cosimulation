function [J1,J2,J3,J4]=Jac_calculation(v,Yth,VaT,VbT,VcT)

m = conj(0.0026 - 0.581*1i)/(2*abs(1));
%% thevenin voltages
V_ath=VaT;
V_bth=VbT;
V_cth=VcT;

a = ((2*abs(V_ath)*conj(Yth(1)))+(V_ath*conj(V_bth)*conj(Yth(2)))/abs(V_ath)+(V_ath*conj(V_cth)*conj(Yth(3)))/abs(V_ath));
b = (V_ath*conj(V_bth)*conj(Yth(2)))/abs(V_bth);
c = (V_ath*conj(V_cth)*conj(Yth(3)))/abs(V_cth);
d = (V_bth*conj(V_ath)*conj(Yth(4)))/abs(V_ath);
e = ((V_bth*conj(V_ath)*conj(Yth(4)))/abs(V_bth)+(2*abs(V_bth)*conj(Yth(5)))+(V_bth*conj(V_cth)*conj(Yth(6)))/abs(V_bth));
f = (V_bth*conj(V_cth)*conj(Yth(6)))/abs(V_cth);
g = (V_cth*conj(V_ath)*conj(Yth(7)))/abs(V_ath);
h = (V_cth*conj(V_bth)*conj(Yth(8)))/abs(V_bth);
j2 =((V_cth*conj(V_ath)*conj(Yth(7)))/abs(V_cth)+(V_cth*conj(V_bth)*conj(Yth(8)))/abs(V_cth)+(2*abs(V_cth)*conj(Yth(9))));

J1 = [1   0  
      1   0  
      1   0  
      0   1  
      0   1 
      0   1];  
      
 J2 = [-real(a) -real(b)  -real(c);
        -real(d) -real(e)  -real(f);
        -real(g) -real(h)  -real(j2);
        -imag(a) -imag(b)  -imag(c);
        -imag(d) -imag(e)  -imag(f);
        -imag(g) -imag(h)  -imag(j2);];
  
 J3 = [  -m 0 0;    
          0 -m 0;
          0 0 -m];
       
 J4 = [(1/3)    (1/3)      (1/3);
       (1/3)    ((1/3)*v)  ((1/3*v^2));
       (1/3)    ((1/3)*v^2) ((1/3)*v);];