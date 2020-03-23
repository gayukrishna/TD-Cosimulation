clc;
clear all;

%% defining transformation matrices and their inverse
v=-0.5+0.8660*(sqrt(-1));
T1= (1/3)*[1 1 1; 1 v (v^2); 1 (v^2) v;];
T2=[1 1 1; 1 (v^2) v; 1 v (v^2);];
T3= (1/(sqrt(3)))*[1 1 1; 1 v (v^2); 1 (v^2) v;]; %Power invariant sequence transformation matrix

%% start of the entire simulation
a_bus25=[];
a_bus26=[];
a_bus28=[];
r_bus25=[];
r_bus26=[];
r_bus28=[];
error_cosim=1; 
count_cosim=0;

V_asub_bus25=1.05;
V_bsub_bus25=1.05;
V_csub_bus25=1.05;

V_asub_bus26=1.05;
V_bsub_bus26=1.05;
V_csub_bus26=1.05;

V_asub_bus28=1.05;
V_bsub_bus28=1.05;
V_csub_bus28=1.05;

V_aangsub_bus25=0;
V_bangsub_bus25=240;
V_cangsub_bus25=120;

V_aangsub_bus26=0;
V_bangsub_bus26=240;
V_cangsub_bus26=120;

V_aangsub_bus28=0;
V_bangsub_bus28=240;
V_cangsub_bus28=120;

P_asub=0.1748;
P_bsub=0.1727;
P_csub=0.1736;

Q_asub=0.0398;
Q_bsub=0.0399;
Q_csub=0.0379;

Pl_update_bus25=[P_asub; P_bsub; P_csub];%.*10;
Ql_update_bus25=[Q_asub; Q_bsub; Q_csub];%.*25;

Pl_update_bus26=[P_asub; P_bsub; P_csub];%.*10;
Ql_update_bus26=[Q_asub; Q_bsub; Q_csub];%.*25;

Pl_update_bus28=[P_asub; P_bsub; P_csub];%.*5;
Ql_update_bus28=[Q_asub; Q_bsub; Q_csub];%.*25;

Sl_bus25_update = Pl_update_bus25 + 1i*Ql_update_bus25;

Sl_bus26_update = Pl_update_bus26 + 1i*Ql_update_bus26;

Sl_bus28_update = Pl_update_bus28 + 1i*Ql_update_bus28;

while error_cosim>0.0001
%% solve transmission load flow
[S_bus26,VaT26,VbT26,VcT26,V1T26,S_bus25,VaT25,VbT25,VcT25,V1T25,S_bus28,VaT28,VbT28,VcT28,V1T28]=Trans_loadflow_multifeeder(Sl_bus25_update,Sl_bus26_update,Sl_bus28_update);

V_abus25 = abs(VaT25);
V_bbus25 = abs(VbT25);
V_cbus25 = abs(VcT25);
V_aangbus25 = angle(VaT25)*180/pi;
V_bangbus25 = angle(VbT25)*180/pi;
V_cangbus25 = angle(VcT25)*180/pi;

V_abus26 = abs(VaT26);
V_bbus26 = abs(VbT26);
V_cbus26 = abs(VcT26);
V_aangbus26 = angle(VaT26)*180/pi;
V_bangbus26 = angle(VbT26)*180/pi;
V_cangbus26 = angle(VcT26)*180/pi;

V_abus28 = abs(VaT28);
V_bbus28 = abs(VbT28);
V_cbus28 = abs(VcT28);
V_aangbus28 = angle(VaT28)*180/pi;
V_bangbus28 = angle(VbT28)*180/pi;
V_cangbus28 = angle(VcT28)*180/pi;

%% solving distribution load flow 
Cktnum = 25;
[S_phckt24sub_bus25]=Dist_loadflow_multifeeder(Cktnum,T1,V_asub_bus25,V_bsub_bus25,V_csub_bus25,V_aangsub_bus25,V_bangsub_bus25,V_cangsub_bus25);

%% solving distribution load flow 
Cktnum = 26;
[S_phckt24sub_bus26]=Dist_loadflow_multifeeder(Cktnum,T1,V_asub_bus26,V_bsub_bus26,V_csub_bus26,V_aangsub_bus26,V_bangsub_bus26,V_cangsub_bus26);

%% solving distribution load flow 
Cktnum = 28;
[S_phckt24sub_bus28]=Dist_loadflow_multifeeder(Cktnum,T1,V_asub_bus28,V_bsub_bus28,V_csub_bus28,V_aangsub_bus28,V_bangsub_bus28,V_cangsub_bus28);

%% computing the difference and updates
R_matrix=[real(S_phckt24sub_bus26)-real(S_bus26);
          imag(S_phckt24sub_bus26)-imag(S_bus26);
          real(S_phckt24sub_bus25)-real(S_bus25);
          imag(S_phckt24sub_bus25)-imag(S_bus25);
          real(S_phckt24sub_bus28)-real(S_bus28);
          imag(S_phckt24sub_bus28)-imag(S_bus28);
          V_abus26-V_asub_bus26;
          V_bbus26-V_bsub_bus26;
          V_cbus26-V_csub_bus26;
          V_abus25-V_asub_bus25;
          V_bbus25-V_bsub_bus25;
          V_cbus25-V_csub_bus25;
          V_abus28-V_asub_bus28;
          V_bbus28-V_bsub_bus28;
          V_cbus28-V_csub_bus28;];
    
%% give input updates for next iteration in cosim
% values at bus 26
P_diff_bus26=real(S_phckt24sub_bus26)-real(S_bus26);
Q_diff_bus26=imag(S_phckt24sub_bus26)-imag(S_bus26);
Pl_update_bus26=real(S_bus26)+P_diff_bus26;
Ql_update_bus26=imag(S_bus26)+Q_diff_bus26;
V_asub_bus26=V_asub_bus26+(V_abus26-V_asub_bus26);
V_bsub_bus26=V_bsub_bus26+(V_bbus26-V_bsub_bus26);
V_csub_bus26=V_csub_bus26+(V_cbus26-V_csub_bus26);
V_aangsub_bus26=V_aangbus26;
V_bangsub_bus26=V_bangbus26;
V_cangsub_bus26=V_cangbus26;
V_bus26 =[(V_asub_bus26); (V_bsub_bus26); (V_csub_bus26)]; 
[vunb_bus26]= unbal_cal(V_bus26);
% Vseq6=T1*[(V_asub_bus6); (V_bsub_bus6); (V_csub_bus6)];
% V0T6=Vseq6(1,1);
% V1T6=Vseq6(2,1);
% V2T6=Vseq6(3,1);
Sl_bus26_update = Pl_update_bus26 + 1i*Ql_update_bus26;

% values at bus 25
P_diff_bus25=real(S_phckt24sub_bus25)-real(S_bus25);
Q_diff_bus25=imag(S_phckt24sub_bus25)-imag(S_bus25);
Pl_update_bus25=real(S_bus25)+P_diff_bus25;
Ql_update_bus25=imag(S_bus25)+Q_diff_bus25;
V_asub_bus25=V_asub_bus25+(V_abus25-V_asub_bus25);
V_bsub_bus25=V_bsub_bus25+(V_bbus25-V_bsub_bus25);
V_csub_bus25=V_csub_bus25+(V_cbus25-V_csub_bus25);
V_aangsub_bus25=V_aangbus25;
V_bangsub_bus25=V_bangbus25;
V_cangsub_bus25=V_cangbus25;
V_bus25 =[(V_asub_bus25); (V_bsub_bus25); (V_csub_bus25)]; 
[vunb_bus25]= unbal_cal(V_bus25);
% Vseq5=T1*[(V_asub_bus5); (V_bsub_bus5); (V_csub_bus5)];
% V0T5=Vseq5(1,1);
% V1T5=Vseq5(2,1);
% V2T5=Vseq5(3,1);
Sl_bus25_update = Pl_update_bus25 + 1i*Ql_update_bus25;

% values at bus 28
P_diff_bus28=real(S_phckt24sub_bus28)-real(S_bus28);
Q_diff_bus28=imag(S_phckt24sub_bus28)-imag(S_bus28);
Pl_update_bus28=real(S_bus28)+P_diff_bus28;
Ql_update_bus28=imag(S_bus28)+Q_diff_bus28;
V_asub_bus28=V_asub_bus28+(V_abus28-V_asub_bus28);
V_bsub_bus28=V_bsub_bus28+(V_bbus28-V_bsub_bus28);
V_csub_bus28=V_csub_bus28+(V_cbus28-V_csub_bus28);
V_aangsub_bus28=V_aangbus28;
V_bangsub_bus28=V_bangbus28;
V_cangsub_bus28=V_cangbus28;
V_bus28 =[(V_asub_bus28); (V_bsub_bus28); (V_csub_bus28)]; 
[vunb_bus28]= unbal_cal(V_bus28);
% Vseq8=T1*[(V_asub_bus8); (V_bsub_bus8); (V_csub_bus8)];
% V0T8=Vseq8(1,:);
% V1T8=Vseq8(2,:);
% V2T8=Vseq8(3,:);
Sl_bus28_update = Pl_update_bus28 + 1i*Ql_update_bus28;

%% output the cosim_iteration values and error
% disp('The residual vector is=');
% disp(R_matrix);   
error_cosim = max(abs(R_matrix)); 
count_cosim = count_cosim + 1;

disp('The co-simulation ran with ');
disp('Error_cosim=');
disp(error_cosim);
disp('iterations_cosim= ');
disp(count_cosim);

a_bus26=[a_bus26;V_asub_bus26 V_bsub_bus26 V_csub_bus26 V_aangsub_bus26 V_bangsub_bus26 V_cangsub_bus26 real(S_phckt24sub_bus26(1,1)) real(S_phckt24sub_bus26(2,1)) real(S_phckt24sub_bus26(3,1)) imag(S_phckt24sub_bus26(1,1)) imag(S_phckt24sub_bus26(2,1)) imag(S_phckt24sub_bus26(3,1))];
a_bus25=[a_bus25;V_asub_bus25 V_bsub_bus25 V_csub_bus25 V_aangsub_bus25 V_bangsub_bus25 V_cangsub_bus25 real(S_phckt24sub_bus25(1,1)) real(S_phckt24sub_bus25(2,1)) real(S_phckt24sub_bus25(3,1)) imag(S_phckt24sub_bus25(1,1)) imag(S_phckt24sub_bus25(2,1)) imag(S_phckt24sub_bus25(3,1))];
a_bus28=[a_bus28;V_asub_bus28 V_bsub_bus28 V_csub_bus28 V_aangsub_bus28 V_bangsub_bus28 V_cangsub_bus28 real(S_phckt24sub_bus28(1,1)) real(S_phckt24sub_bus28(2,1)) real(S_phckt24sub_bus28(3,1)) imag(S_phckt24sub_bus28(1,1)) imag(S_phckt24sub_bus28(2,1)) imag(S_phckt24sub_bus28(3,1))];

r_bus26=abs([r_bus26;R_matrix(1,1) R_matrix(2,1) R_matrix(3,1) R_matrix(4,1) R_matrix(5,1) R_matrix(6,1) R_matrix(19,1) R_matrix(20,1) R_matrix(21,1)]); 
r_bus25=abs([r_bus25;R_matrix(7,1) R_matrix(8,1) R_matrix(9,1) R_matrix(10,1) R_matrix(11,1) R_matrix(12,1) R_matrix(22,1) R_matrix(23,1) R_matrix(24,1)]); 
r_bus28=abs([r_bus28;R_matrix(13,1) R_matrix(14,1) R_matrix(15,1) R_matrix(16,1) R_matrix(17,1) R_matrix(18,1) R_matrix(25,1) R_matrix(26,1) R_matrix(27,1)]);
disp('..........................................................');
end

disp ('output at bus 26')
disp ('The converged power values at bus 26 are=')
disp (real(S_phckt24sub_bus26));
disp (imag(S_phckt24sub_bus26));
disp('Sequence Voltages at the interface at bus 26=');
%disp(V0T6);
disp(V1T26);
%disp(V2T6);
disp('Phase Voltages at the interface at bus 26=');
disp(V_asub_bus26);
disp(V_bsub_bus26);
disp(V_csub_bus26);
disp('Phase angles at the interface at bus 26=');
disp(V_aangsub_bus26);
disp(V_bangsub_bus26);
disp(V_cangsub_bus26);

disp('..........................................................');

disp ('output at bus 25')
disp ('The converged power values at bus 25 are=')
disp (real(S_phckt24sub_bus25));
disp (imag(S_phckt24sub_bus25));
disp('Sequence Voltages at the interface at bus 25=');
%disp(V0T5);
disp(V1T25);
%disp(V2T5);
disp('Phase Voltages at the interface at bus 25=');
disp(V_asub_bus25);
disp(V_bsub_bus25);
disp(V_csub_bus25);
disp('Phase angles at the interface at bus 25=');
disp(V_aangsub_bus25);
disp(V_bangsub_bus25);
disp(V_cangsub_bus25);

disp('..........................................................');

disp ('output at bus 28')
disp ('The converged power values at bus 28 are=')
disp (real(S_phckt24sub_bus28));
disp (imag(S_phckt24sub_bus28));
disp('Sequence Voltages at the interface at bus 28=');
%disp(V0T8);
disp(V1T28);
%disp(V2T8);
disp('Phase Voltages at the interface at bus 28=');
disp(V_asub_bus28);
disp(V_bsub_bus28);
disp(V_csub_bus28);
disp('Phase angles at the interface at bus 28=');
disp(V_aangsub_bus28);
disp(V_bangsub_bus28);
disp(V_cangsub_bus28);

disp('The unbalances are=')
disp(vunb_bus26);
disp(vunb_bus25);
disp(vunb_bus28);
