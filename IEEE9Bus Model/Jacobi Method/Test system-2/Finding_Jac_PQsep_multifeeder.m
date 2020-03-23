 
function[Jacobian,J1,J2,J3,J4]=Finding_Jac_PQsep_multifeeder(v,VaT6,VbT6,VcT6)
% - Start OpenDSS
DSSObj=actxserver('OpenDSSEngine.DSS');
if ~DSSObj.Start(0)
    disp('Unable to start openDSS Engine');
    return
end

DSSText=DSSObj.Text;
DSSCircuit=DSSObj.ActiveCircuit;
DSSBus = DSSCircuit.ActiveBus;

DSSText.Command='Compile (C:\Users\user\Dropbox\personal\MATLABCODES\cosim_codes\Codes\Jacobi method\Test system-2\Ckt24\master_ckt24_cosim.dss)';
%DSSText.Command = 'Redirect PV_combined100per.txt';
%DSSText.command='BatchEdit Load..* yearly=default';
DSSText.Command='solve loadmult=2.5';
DSSText.Command='solve mode=Faultstudy';
DSSCircuit.SetActiveBus('SOURCEBUS');
DSSCircuit.ActiveBus.ZscRefresh
ZscMatrix= DSSCircuit.ActiveBus.ZscMatrix;
ThVoltage = DSSBus.puVmagAngle;
%ZscMatrix

%% findind the complex conj of the matrix
Zaa=(ZscMatrix(1,1)+(ZscMatrix(1,2)*sqrt(-1)));
Zab=(ZscMatrix(1,3)+(ZscMatrix(1,4)*sqrt(-1)));
Zac=(ZscMatrix(1,5)+(ZscMatrix(1,6)*sqrt(-1)));
Zba=(ZscMatrix(1,7)+(ZscMatrix(1,8)*sqrt(-1)));
Zbb=(ZscMatrix(1,9)+(ZscMatrix(1,10)*sqrt(-1)));
Zbc=(ZscMatrix(1,11)+(ZscMatrix(1,12)*sqrt(-1)));
Zca=(ZscMatrix(1,13)+(ZscMatrix(1,14)*sqrt(-1)));
Zcb=(ZscMatrix(1,15)+(ZscMatrix(1,16)*sqrt(-1)));
Zcc=(ZscMatrix(1,17)+(ZscMatrix(1,18)*sqrt(-1)));

Z = [Zaa,Zab,Zac;Zba,Zbb,Zbc;Zca,Zcb,Zcc];
Y = inv(Z);
Yaa= Y(1,1);
Yab= Y(1,2);
Yac=Y(1,3);
Yba=Y(2,1);
Ybb=Y(2,2);
Ybc=Y(2,3);
Yca=Y(3,1);
Ycb=Y(3,2);
Ycc=Y(3,3);


%% thevenin voltages
V_ath=VaT6;
V_bth=VbT6;
V_cth=VcT6;

%v=-0.5+0.8660*(sqrt(-1));

%% finding the jacobain matrix
%m0 = conj(0.0138 + 0.0947i)/(2*abs(V_1bus6));
m = conj(0.0026 - 0.581*1i)/(2*abs(1));

% a = ((2*V_ath*conj(Yaa))+((v)*V_bth*conj(Yba))+(v^2*V_cth*conj(Yca)));
% b = (V_ath*v*conj(Yba));
% c = (V_ath*(v^2)*conj(Yca));
% d = (V_bth*(v^2)*conj(Yab));
% e = ((v^2*V_ath*conj(Yab))+(2*V_bth*conj(Ybb))+((v)*V_cth*conj(Ycb)));
% f = (V_bth*v*conj(Ycb));
% g = (V_cth*v*conj(Yac));
% h = (V_cth*(v^2)*conj(Ybc)) ;
% j2 = ((v*V_ath*conj(Yac))+(v^2*V_bth*conj(Ybc))+(2*V_cth*conj(Ycc)));
% k = real(m);
% l = imag(m);

a = ((2*abs(V_ath)*conj(Yaa))+(V_ath*conj(V_bth)*conj(Yba))/abs(V_ath)+(V_ath*conj(V_cth)*conj(Yca))/abs(V_ath));
b = (V_ath*conj(V_bth)*conj(Yba))/abs(V_bth);
c = (V_ath*conj(V_cth)*conj(Yca))/abs(V_cth);
d = (V_bth*conj(V_ath)*conj(Yab))/abs(V_ath);
e = ((V_bth*conj(V_ath)*conj(Yab))/abs(V_bth)+(2*abs(V_bth)*conj(Ybb))+(V_bth*conj(V_cth)*conj(Ycb))/abs(V_bth));
f = (V_bth*conj(V_cth)*conj(Ycb))/abs(V_cth);
g = (V_cth*conj(V_ath)*conj(Yac))/abs(V_ath);
h = (V_cth*conj(V_bth)*conj(Ybc))/abs(V_bth);
j2 =((V_cth*conj(V_ath)*conj(Yac))/abs(V_cth)+(V_cth*conj(V_bth)*conj(Ybc))/abs(V_cth)+(2*abs(V_cth)*conj(Ycc)));
k = real(m);
l = imag(m);


Jacobian=[1   0  -real(a) -real(b)  -real(c);
          1   0  -real(d) -real(e)  -real(f);
          1   0  -real(g) -real(h)  -real(j2);
          0   1  -imag(a) -imag(b)  -imag(c);
          0   1  -imag(d) -imag(e)  -imag(f);
          0   1  -imag(g) -imag(h)  -imag(j2);
          0   0   (1/3)    (1/3)      (1/3);
          -k -l   (1/3)    ((1/3)*v)  ((1/3*v^2));
          0   0   (1/3)    ((1/3)*v^2) ((1/3)*v);];
%       

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


release(DSSObj)
%disp(m);
