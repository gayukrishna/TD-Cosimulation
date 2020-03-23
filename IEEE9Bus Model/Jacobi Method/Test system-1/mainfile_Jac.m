clc;
clear all;

%% defining transformation matrices and their inverse
v=-0.5+0.8660*(sqrt(-1));
T1= (1/3)*[1 1 1; 1 v (v^2); 1 (v^2) v;];
T2=[1 1 1; 1 (v^2) v; 1 v (v^2);];
T3= (1/(sqrt(3)))*[1 1 1; 1 v (v^2); 1 (v^2) v;]; %Power invariant sequence transformation matrix

%% start of the entire simulation
error_cosim=1; 
count_cosim=0;
V_asub=1.05;
V_bsub=1.05;
V_csub=1.05;
V_aangsub=0;
V_bangsub=240;
V_cangsub=120;
Pl_update=[0.178; 0.178; 0.175];
Ql_update=[0.051; 0.052; 0.05];
%%
[Y1,Y0,P_bus6,Q_bus6,S_bus6,P_seqbus6,Q_seqbus6,VaT6,VbT6,VcT6,V_0bus6,V_1bus6,V_2bus6,V_aangbus6,V_bangbus6,V_cangbus6]=Trans_loadflow(Pl_update,Ql_update);
V_abus6 = abs(VaT6);
V_bbus6 = abs(VbT6);
V_cbus6 = abs(VcT6);
S1 = P_seqbus6(2,1)+1i*Q_seqbus6(2,1);

%% solving distribution load flow
[S_phckt24sub,P_phckt24sub,Q_phckt24sub,V_0sub,V_1sub,V_2sub]=Dist_loadflow_Jac(T1,T3,V_asub,V_bsub,V_csub,V_aangsub,V_bangsub,V_cangsub);

    %% finding jacobian
[Jacobian,J1,J2,J3,J4]=Finding_Jac_PQsep(v,VaT6,VbT6,VcT6);

clc
while error_cosim>0.0001
%% solve transmission load flow
[Y1,Y0,P_bus6,Q_bus6,S_bus6,P_seqbus6,Q_seqbus6,VaT6,VbT6,VcT6,V_0bus6,V_1bus6,V_2bus6,V_aangbus6,V_bangbus6,V_cangbus6]=Trans_loadflow(Pl_update,Ql_update);
V_abus6 = abs(VaT6);
V_bbus6 = abs(VbT6);
V_cbus6 = abs(VcT6);
S1 = P_seqbus6(2,1)+1i*Q_seqbus6(2,1);

%% solving distribution load flow
[S_phckt24sub,P_phckt24sub,Q_phckt24sub,V_0sub,V_1sub,V_2sub]=Dist_loadflow_Jac(T1,T3,V_asub,V_bsub,V_csub,V_aangsub,V_bangsub,V_cangsub);

%% computing the difference and updates      
 D_matrix=-[real(S_bus6)-real(S_phckt24sub);
     imag(S_bus6)-imag(S_phckt24sub);
          (V_0sub-V_0bus6);
          (V_1sub-V_1bus6);
          (V_2sub-V_2bus6);];

%% Voltage and positive sequence power update to the transmissionsystem
diffV1 = inv(J4)*D_matrix(7:9);
diffV = diffV1;
Vabc_Dtemp = T2*([V_0sub;V_1sub;V_2sub]);
Vabc_update = Vabc_Dtemp+diffV;
Z = 0.0141 + 0.0958i;
m1 = conj(Z/Vabc_update(1,1));
m2 = conj(Z/Vabc_update(2,1));
m3 = conj(Z/Vabc_update(3,1));
J3 = [m1,m2,m3; 
      m1,m2./(v^2),m3./(v);
      m1,m2./(v),m3./(v^2);];
Sl_update1 =0;

for i = 1:10
Vabc_update1 = Vabc_update;
diffPQabc = eye(6,6)*(D_matrix(1:6)-J2*(abs(Vabc_update)-abs(Vabc_Dtemp)));
alpha = 1;
Sa_update = S_bus6(1)+alpha.*(diffPQabc(1,1)+1i*diffPQabc(4,1));
Sb_update = S_bus6(2)+alpha.*(diffPQabc(2,1)+1i*diffPQabc(5,1));
Sc_update = S_bus6(3)+alpha.*(diffPQabc(3,1)+1i*diffPQabc(6,1));
Sl_update = [Sa_update;Sb_update;Sc_update];
%%
diffV1 = inv(J4)*(D_matrix(7:9)-J3*(Sl_update-S_bus6));
Vabc_update = Vabc_Dtemp+diffV1;
if max(abs(Sl_update-Sl_update1))<0.00001
    break;
end
Sl_update1 = Sl_update;
end 
%%
Pl_update = real(Sl_update);
Ql_update = imag(Sl_update);

%% Voltage phase update to the distribution system
V_asub=abs(Vabc_update(1));
V_bsub=abs(Vabc_update(2));
V_csub=abs(Vabc_update(3));

%% voltage angle update to the distribution system
V_aangsub=angle(Vabc_update(1))*180/pi;
V_bangsub=angle(Vabc_update(2))*180/pi;
V_cangsub=angle(Vabc_update(3))*180/pi;

%% output the cosim_iteration values and error
error_cosim = max(abs(D_matrix)); 
count_cosim = count_cosim + 1;
disp(V_abus6);
disp(V_bbus6);
disp(V_cbus6);
disp('The co-simulation ran with ');
disp('Error_cosim=');
disp(error_cosim);
disp('iterations_cosim= ');
disp(count_cosim);
disp('..........................................................');
end

