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

V_asub_bus6=1.05;
V_bsub_bus6=1.05;
V_csub_bus6=1.05;

V_asub_bus5=1.05;
V_bsub_bus5=1.05;
V_csub_bus5=1.05;

V_asub_bus8=1.05;
V_bsub_bus8=1.05;
V_csub_bus8=1.05;

V_aangsub_bus6=0;
V_bangsub_bus6=240;
V_cangsub_bus6=120;

V_aangsub_bus5=0;
V_bangsub_bus5=240;
V_cangsub_bus5=120;

V_aangsub_bus8=0;
V_bangsub_bus8=240;
V_cangsub_bus8=120;

Pl_update_bus5=[0.178; 0.178; 0.175];
Ql_update_bus5=[0.051; 0.052; 0.05];

Pl_update_bus6=[0.178; 0.178; 0.175];
Ql_update_bus6=[0.051; 0.052; 0.05];

Pl_update_bus8=[0.178; 0.178; 0.175];
Ql_update_bus8=[0.051; 0.052; 0.05];

Sl_bus6_update = Pl_update_bus6 + 1i*Ql_update_bus6;

Sl_bus5_update = Pl_update_bus5 + 1i*Ql_update_bus5;

Sl_bus8_update = Pl_update_bus8 + 1i*Ql_update_bus8;

count = 1;
%S_iter(count) = [Sl_bus6_update,Sl_bus6_update,Sl_bus6_update];


while error_cosim>0.0001
%% solve transmission load flow
[S_bus6,VaT6,VbT6,VcT6,S_bus5,VaT5,VbT5,VcT5,S_bus8,VaT8,VbT8,VcT8]=Trans_loadflow_multifeeder(Sl_bus6_update,Sl_bus5_update,Sl_bus8_update);
%[P_bus6,Q_bus6,V_abus6,V_bbus6,V_cbus6]=NR_tlf(Pl_update,Ql_update);

Strans_iter(:,:,count) = [S_bus6,S_bus5,S_bus8];
Vtrans_iter(:,:,count) = [[VaT6;VbT6;VcT6],[VaT5;VbT5;VcT5],[VaT8;VbT8;VcT8]];

V_abus6 = abs(VaT6);
V_bbus6 = abs(VbT6);
V_cbus6 = abs(VcT6);

V_aangbus6 = angle(VaT6)*180/pi;
V_bangbus6 = angle(VbT6)*180/pi;
V_cangbus6 = angle(VcT6)*180/pi;

V_abus5 = abs(VaT5);
V_bbus5 = abs(VbT5);
V_cbus5 = abs(VcT5);

V_aangbus5 = angle(VaT5)*180/pi;
V_bangbus5 = angle(VbT5)*180/pi;
V_cangbus5 = angle(VcT5)*180/pi;

V_abus8 = abs(VaT8);
V_bbus8 = abs(VbT8);
V_cbus8 = abs(VcT8);

V_aangbus8 = angle(VaT8)*180/pi;
V_bangbus8 = angle(VbT8)*180/pi;
V_cangbus8 = angle(VcT8)*180/pi;

%% solving distribution load flow 
Cktnum = 6;
[S_phckt24sub_bus6]=Dist_loadflow_multifeeder(Cktnum,T1,V_asub_bus6,V_bsub_bus6,V_csub_bus6,V_aangsub_bus6,V_bangsub_bus6,V_cangsub_bus6);

%% solving distribution load flow 
Cktnum = 5;
[S_phckt24sub_bus5]=Dist_loadflow_multifeeder(Cktnum,T1,V_asub_bus5,V_bsub_bus5,V_csub_bus5,V_aangsub_bus5,V_bangsub_bus5,V_cangsub_bus5);

%% solving distribution load flow 
Cktnum = 8;
[S_phckt24sub_bus8]=Dist_loadflow_multifeeder(Cktnum,T1,V_asub_bus8,V_bsub_bus8,V_csub_bus8,V_aangsub_bus8,V_bangsub_bus8,V_cangsub_bus8);

%%

Sdist_iter(:,:,count) = [S_phckt24sub_bus6,S_phckt24sub_bus5,S_phckt24sub_bus8];
Vdist_iter(:,:,count) = [[V_asub_bus6*(cos(V_aangsub_bus6)+1i*sin(V_aangsub_bus6));V_bsub_bus6*(cos(V_bangsub_bus6)+1i*sin(V_bangsub_bus6));V_csub_bus6*(cos(V_cangsub_bus6)+1i*sin(V_cangsub_bus6))],[V_asub_bus5*(cos(V_aangsub_bus5)+1i*sin(V_aangsub_bus5));V_bsub_bus5*(cos(V_bangsub_bus5)+1i*sin(V_bangsub_bus5));V_csub_bus5*(cos(V_cangsub_bus5)+1i*sin(V_cangsub_bus5))],[V_asub_bus8*(cos(V_aangsub_bus8)+1i*sin(V_aangsub_bus8));V_bsub_bus8*(cos(V_bangsub_bus8)+1i*sin(V_bangsub_bus8));V_csub_bus6*(cos(V_cangsub_bus6)+1i*sin(V_cangsub_bus8))]];

%% computing the difference and updates
% 
R_matrix=[real(S_phckt24sub_bus6)-real(S_bus6);
          imag(S_phckt24sub_bus6)-imag(S_bus6);
          real(S_phckt24sub_bus5)-real(S_bus5);
          imag(S_phckt24sub_bus5)-imag(S_bus5);
          real(S_phckt24sub_bus8)-real(S_bus8);
          imag(S_phckt24sub_bus8)-imag(S_bus8);
          V_abus6-V_asub_bus6;
          V_bbus6-V_bsub_bus6;
          V_cbus6-V_csub_bus6;
          V_abus5-V_asub_bus5;
          V_bbus5-V_bsub_bus5;
          V_cbus5-V_csub_bus5;
          V_abus8-V_asub_bus8;
          V_bbus8-V_bsub_bus8;
          V_cbus8-V_csub_bus8;];
    
%% give input updates for next iteration in cosim

P_diff_bus6=real(S_phckt24sub_bus6)-real(S_bus6);
Q_diff_bus6=imag(S_phckt24sub_bus6)-imag(S_bus6);



Pl_update_bus6=real(S_bus6)+P_diff_bus6;
Ql_update_bus6=imag(S_bus6)+Q_diff_bus6;
V_asub_bus6=V_asub_bus6+(V_abus6-V_asub_bus6);
V_bsub_bus6=V_bsub_bus6+(V_bbus6-V_bsub_bus6);
V_csub_bus6=V_csub_bus6+(V_cbus6-V_csub_bus6);
V_aangsub_bus6=V_aangbus6;
V_bangsub_bus6=V_bangbus6;
V_cangsub_bus6=V_cangbus6;
 
Sl_bus6_update = Pl_update_bus6 + 1i*Ql_update_bus6;

%%
P_diff_bus5=real(S_phckt24sub_bus5)-real(S_bus5);
Q_diff_bus5=imag(S_phckt24sub_bus5)-imag(S_bus5);
Pl_update_bus5=real(S_bus5)+P_diff_bus5;
Ql_update_bus5=imag(S_bus5)+Q_diff_bus5;
V_asub_bus5=V_asub_bus5+(V_abus5-V_asub_bus5);
V_bsub_bus5=V_bsub_bus5+(V_bbus5-V_bsub_bus5);
V_csub_bus5=V_csub_bus5+(V_cbus5-V_csub_bus5);
V_aangsub_bus5=V_aangbus5;
V_bangsub_bus5=V_bangbus5;
V_cangsub_bus5=V_cangbus5;
 
Sl_bus5_update = Pl_update_bus5 + 1i*Ql_update_bus5;
%%
P_diff_bus8=real(S_phckt24sub_bus8)-real(S_bus8);
Q_diff_bus8=imag(S_phckt24sub_bus8)-imag(S_bus8);
Pl_update_bus8=real(S_bus8)+P_diff_bus8;
Ql_update_bus8=imag(S_bus8)+Q_diff_bus8;
V_asub_bus8=V_asub_bus8+(V_abus8-V_asub_bus8);
V_bsub_bus8=V_bsub_bus8+(V_bbus8-V_bsub_bus8);
V_csub_bus8=V_csub_bus8+(V_cbus8-V_csub_bus8);
V_aangsub_bus8=V_aangbus8;
V_bangsub_bus8=V_bangbus8;
V_cangsub_bus8=V_cangbus8;
 
Sl_bus8_update = Pl_update_bus8 + 1i*Ql_update_bus8;
count=count+1;
%% output the cosim_iteration values and error
disp('The residual vector is=');
disp(R_matrix);   

error_cosim = max(abs(R_matrix)); 
count_cosim = count_cosim + 1;

disp('The co-simulation ran with ');
disp('Error_cosim=');
disp(error_cosim);
disp('iterations_cosim= ');
disp(count_cosim);
disp('..........................................................');
end
