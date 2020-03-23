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
% Pl_update=[0.1748; 0.1727; 0.1736];
% Ql_update=[0.0398; 0.0399; 0.0379];
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

%%
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
Vabc_Dtemp_old = T2*([V_0sub;V_1sub;V_2sub]);
%Vabc_seq_old = T1*Vabc_Dtemp_old;
Iasptemp_old=conj(S_phckt24sub(1,1)/Vabc_Dtemp_old(1,1));
Ibsptemp_old=conj(S_phckt24sub(2,1)/Vabc_Dtemp_old(2,1));
Icsptemp_old=conj(S_phckt24sub(3,1)/Vabc_Dtemp_old(3,1));
Itemp_old = T1*[Iasptemp_old;Ibsptemp_old;Icsptemp_old];

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
Vabc_update1 = Vabc_update;
diffPQ = pinv(J1)*(D_matrix(1:6)-J2*(abs(Vabc_update)-abs(Vabc_Dtemp)));
alpha = 2.5;
P1_update=P_seqbus6(2,1)+alpha*(diffPQ(1));
Q1_update=Q_seqbus6(2,1)+alpha*(diffPQ(2));
Vabc_seq = T1*Vabc_update;
 Z = 0.0141 - 0.5148i;
 m = conj(Z)/Vabc_seq(2);
%%
Vabc_seq = T1*Vabc_update;
Iasptemp=conj(S_phckt24sub(1,1)/Vabc_update(1,1));
Ibsptemp=conj(S_phckt24sub(2,1)/Vabc_update(2,1));
Icsptemp=conj(S_phckt24sub(3,1)/Vabc_update(3,1));
Itemp = T1*[Iasptemp;Ibsptemp;Icsptemp];
dI0 = Itemp(1) - Itemp_old(1);
dS1 = diffPQ(1)+1i*diffPQ(2);
dI2 = Itemp(3) - Itemp_old(3);
%Z = 0.0141 - 0.5148i;
Z = 0.0141 - 0.0958i;
m = conj(Z)/Vabc_seq(2); 
J3 = [-m,0,0;0,m,0;0,0,-m];
diffV1 = inv(J4)*(D_matrix(7:9)-(J3*[dI0;dS1;dI2]));
Vabc_update = Vabc_Dtemp + diffV1;
Vabc_seq = T1*Vabc_update;

% diffPQ = pinv(J1)*(D_matrix(1:6)-J2*(abs(Vabc_update)-abs(Vabc_Dtemp)));
% P1_update=P_seqbus6(2,1)+alpha*(diffPQ(1));
% Q1_update=Q_seqbus6(2,1)+alpha*(diffPQ(2));
%S0_update=(P_seqbus6(1,1)+1i*Q_seqbus6(1,1))+[1/3*(diffV(1,1)+v^2*diffV(2,1)+v*diffV(3,1))-D_matrix(7,1)]/m;
%S2_update=(P_seqbus6(3,1)+1i*Q_seqbus6(3,1))+[1/3*(diffV(1,1)+v*diffV(2,1)+v^2*diffV(3,1))-D_matrix(9,1)]/m;

S1_update = P1_update + 1i*Q1_update;
Iasptemp=conj(S_phckt24sub(1,1)/Vabc_update(1,1));
Ibsptemp=conj(S_phckt24sub(2,1)/Vabc_update(2,1));
Icsptemp=conj(S_phckt24sub(3,1)/Vabc_update(3,1));
Itemp = T1*[Iasptemp;Ibsptemp;Icsptemp];

%%
I0=Itemp(1,1);
I1=conj(S1_update/(3*Vabc_seq(2,1)));
I2=Itemp(3,1);
Itemp_new = T2*[I0;I1;I2];
Sl_update = Vabc_update.*conj(Itemp_new);
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
abs(D_matrix)
error_cosim = max(abs(D_matrix)); 
count_cosim = count_cosim + 1;
disp('The co-simulation ran with ');
disp('Error_cosim=');
disp(error_cosim);
disp('iterations_cosim= ');
disp(count_cosim);
disp('..........................................................');
end

