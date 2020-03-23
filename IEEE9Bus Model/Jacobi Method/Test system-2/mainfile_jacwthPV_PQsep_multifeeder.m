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

Pl_update_bus5=[0.1748; 0.1727; 0.1736];
Ql_update_bus5=[0.0398; 0.0399; 0.0379];

Pl_update_bus6=[0.1748; 0.1727; 0.1736];
Ql_update_bus6=[0.0398; 0.0399; 0.0379];

Pl_update_bus8=[0.1748; 0.1727; 0.1736];
Ql_update_bus8=[0.0398; 0.0399; 0.0379];

Sl_update_bus6 = Pl_update_bus6 + 1i*Ql_update_bus6;

Sl_update_bus5 = Pl_update_bus5 + 1i*Ql_update_bus5;

Sl_update_bus8 = Pl_update_bus8 + 1i*Ql_update_bus8;

%%
[S_bus6,P_seqbus6,Q_seqbus6,VaT6,VbT6,VcT6,S_bus5,P_seqbus5,Q_seqbus5,VaT5,VbT5,VcT5,S_bus8,P_seqbus8,Q_seqbus8,VaT8,VbT8,VcT8]=Trans_loadflow_multifeeder_Jacobian(Sl_update_bus6,Sl_update_bus5,Sl_update_bus8);
%[P_bus6,Q_bus6,V_abus6,V_bbus6,V_cbus6]=NR_tlf(Pl_update,Ql_update);
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

%% solving distribution load flow and %% finding jacobian
Cktnum = 6;
[S_phckt24sub_bus6,V_seqsub_bus6]=Dist_loadflow_multifeeder_Jac(Cktnum,T1,V_asub_bus6,V_bsub_bus6,V_csub_bus6,V_aangsub_bus6,V_bangsub_bus6,V_cangsub_bus6);
[Jacobian_bus6,J1_bus6,J2_bus6,J3_bus6,J4_bus6]=Finding_Jac_PQsep_multifeeder(v,VaT6,VbT6,VcT6);
%% solving distribution load flow and %% finding jacobian
Cktnum = 5;
[S_phckt24sub_bus5,V_seqsub_bus5]=Dist_loadflow_multifeeder_Jac(Cktnum,T1,V_asub_bus5,V_bsub_bus5,V_csub_bus5,V_aangsub_bus5,V_bangsub_bus5,V_cangsub_bus5);
[Jacobian_bus5,J1_bus5,J2_bus5,J3_bus5,J4_bus5]=Finding_Jac_PQsep_multifeeder(v,VaT5,VbT5,VcT5);
%% solving distribution load flow and %% finding jacobian
Cktnum = 8;
[S_phckt24sub_bus8,V_seqsub_bus8]=Dist_loadflow_multifeeder_Jac(Cktnum,T1,V_asub_bus8,V_bsub_bus8,V_csub_bus8,V_aangsub_bus8,V_bangsub_bus8,V_cangsub_bus8);
[Jacobian_bus8,J1_bus8,J2_bus8,J3_bus8,J4_bus8]=Finding_Jac_PQsep_multifeeder(v,VaT8,VbT8,VcT8);
    

%%
while error_cosim>0.0001
%% solve transmission load flow
[S_bus6,P_seqbus6,Q_seqbus6,VaT6,VbT6,VcT6,S_bus5,P_seqbus5,Q_seqbus5,VaT5,VbT5,VcT5,S_bus8,P_seqbus8,Q_seqbus8,VaT8,VbT8,VcT8]=Trans_loadflow_multifeeder_Jacobian(Sl_update_bus6,Sl_update_bus5,Sl_update_bus8);
%[P_bus6,Q_bus6,V_abus6,V_bbus6,V_cbus6]=NR_tlf(Pl_update,Ql_update);
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
%% solving distribution load flow and %% finding jacobian
Cktnum = 6;
[S_phckt24sub_bus6,V_seqsub_bus6]=Dist_loadflow_multifeeder_Jac(Cktnum,T1,V_asub_bus6,V_bsub_bus6,V_csub_bus6,V_aangsub_bus6,V_bangsub_bus6,V_cangsub_bus6);
%[Jacobian_bus6,J1_bus6,J2_bus6,J3_bus6,J4_bus6]=Finding_Jac_PQsep_multifeeder(v,VaT6,VbT6,VcT6);
%% solving distribution load flow and %% finding jacobian
Cktnum = 5;
[S_phckt24sub_bus5,V_seqsub_bus5]=Dist_loadflow_multifeeder_Jac(Cktnum,T1,V_asub_bus5,V_bsub_bus5,V_csub_bus5,V_aangsub_bus5,V_bangsub_bus5,V_cangsub_bus5);
%[Jacobian_bus5,J1_bus5,J2_bus5,J3_bus5,J4_bus5]=Finding_Jac_PQsep_multifeeder(v,VaT5,VbT5,VcT5);
%% solving distribution load flow and %% finding jacobian
Cktnum = 8;
[S_phckt24sub_bus8,V_seqsub_bus8]=Dist_loadflow_multifeeder_Jac(Cktnum,T1,V_asub_bus8,V_bsub_bus8,V_csub_bus8,V_aangsub_bus8,V_bangsub_bus8,V_cangsub_bus8);
%[Jacobian_bus8,J1_bus8,J2_bus8,J3_bus8,J4_bus8]=Finding_Jac_PQsep_multifeeder(v,VaT8,VbT8,VcT8);

%% computing the difference and updates      

V_T012_bus6 = T1*[VaT6;VbT6;VcT6;];
V_T012_bus5 = T1*[VaT5;VbT5;VcT5;];
V_T012_bus8 = T1*[VaT8;VbT8;VcT8;];

V_D012_bus6 = V_seqsub_bus6;
V_D012_bus5 = V_seqsub_bus5;
V_D012_bus8 = V_seqsub_bus8;

D_matrix=-[real(S_bus6)-real(S_phckt24sub_bus6);
            imag(S_bus6)-imag(S_phckt24sub_bus6);
            real(S_bus5)-real(S_phckt24sub_bus5);
            imag(S_bus5)-imag(S_phckt24sub_bus5);
            real(S_bus8)-real(S_phckt24sub_bus8);
            imag(S_bus8)-imag(S_phckt24sub_bus8);
          V_D012_bus6-V_T012_bus6;
          V_D012_bus5-V_T012_bus5;
          V_D012_bus8-V_T012_bus8;];

%% Voltage and positive sequence power update to the transmissionsystem

diffV_bus6 = inv(J4_bus6)*D_matrix(19:21);
diffV_bus5 = inv(J4_bus6)*D_matrix(22:24);
diffV_bus8 = inv(J4_bus6)*D_matrix(25:27);

Vabc_Dtemp_bus6 = T2*([ V_D012_bus6(1); V_D012_bus6(2); V_D012_bus6(3)]);
Vabc_update_bus6 = Vabc_Dtemp_bus6+diffV_bus6;

Vabc_Dtemp_bus5 = T2*([ V_D012_bus5(1); V_D012_bus5(2); V_D012_bus5(3)]);
Vabc_update_bus5 = Vabc_Dtemp_bus5+diffV_bus5;

Vabc_Dtemp_bus8 = T2*([ V_D012_bus8(1); V_D012_bus8(2); V_D012_bus8(3)]);
Vabc_update_bus8 = Vabc_Dtemp_bus8+diffV_bus8;

Vabc_update = [Vabc_update_bus6;Vabc_update_bus5;Vabc_update_bus8];
Vabc_Dtemp = [Vabc_Dtemp_bus6;Vabc_Dtemp_bus5;Vabc_Dtemp_bus8];

J1 = [J1_bus6,zeros(6,2),zeros(6,2);zeros(6,2),J1_bus5,zeros(6,2);zeros(6,2),zeros(6,2),J1_bus8];
J2 = [J2_bus6,zeros(6,3),zeros(6,3);zeros(6,3),J2_bus5,zeros(6,3);zeros(6,3),zeros(6,3),J2_bus8];

diffPQ = pinv(J1)*(D_matrix(1:18)-J2*(abs(Vabc_update)-abs(Vabc_Dtemp)));


alpha = 2.5;

Pseq1_update_bus6=P_seqbus6(2)+alpha*(diffPQ(1));
Qseq1_update_bus6=Q_seqbus6(2)+alpha*(diffPQ(2));

Pseq1_update_bus5=P_seqbus5(2)+alpha*(diffPQ(3));
Qseq1_update_bus5=Q_seqbus5(2)+alpha*(diffPQ(4));

Pseq1_update_bus8=P_seqbus8(2)+alpha*(diffPQ(5));
Qseq1_update_bus8=Q_seqbus8(2)+alpha*(diffPQ(6));

Vabc_seq_bus6 = T1*Vabc_update_bus6;

Vabc_seq_bus5 = T1*Vabc_update_bus5;

Vabc_seq_bus8 = T1*Vabc_update_bus8;


%%

Sseq1_update_bus6 = Pseq1_update_bus6 + 1i*Qseq1_update_bus6;
Sseq1_update_bus5 = Pseq1_update_bus5 + 1i*Qseq1_update_bus5;
Sseq1_update_bus8 = Pseq1_update_bus8 + 1i*Qseq1_update_bus8;

Iasptemp_bus6=conj(S_phckt24sub_bus6(1,1)/Vabc_update_bus6(1,1));
Ibsptemp_bus6=conj(S_phckt24sub_bus6(2,1)/Vabc_update_bus6(2,1));
Icsptemp_bus6=conj(S_phckt24sub_bus6(3,1)/Vabc_update_bus6(3,1));

Iasptemp_bus5=conj(S_phckt24sub_bus5(1,1)/Vabc_update_bus5(1,1));
Ibsptemp_bus5=conj(S_phckt24sub_bus5(2,1)/Vabc_update_bus5(2,1));
Icsptemp_bus5=conj(S_phckt24sub_bus5(3,1)/Vabc_update_bus5(3,1));

Iasptemp_bus8=conj(S_phckt24sub_bus8(1,1)/Vabc_update_bus8(1,1));
Ibsptemp_bus8=conj(S_phckt24sub_bus8(2,1)/Vabc_update_bus8(2,1));
Icsptemp_bus8=conj(S_phckt24sub_bus8(3,1)/Vabc_update_bus8(3,1));

Itemp_bus6 = T1*[Iasptemp_bus6;Ibsptemp_bus6;Icsptemp_bus6];
Itemp_bus5 = T1*[Iasptemp_bus5;Ibsptemp_bus5;Icsptemp_bus5];
Itemp_bus8 = T1*[Iasptemp_bus8;Ibsptemp_bus8;Icsptemp_bus8];

%S0_update = 3*V_0bus6*conj(Itemp(1,1));
%S2_update = 3*V_2bus6*conj(Itemp(3,1));

%%
I0_bus6=Itemp_bus6(1,1);
I1_bus6=conj(Sseq1_update_bus6/(3*Vabc_seq_bus6(2,1)));
I2_bus6=Itemp_bus6(3,1);

Itemp_new_bus6 = T2*[I0_bus6;I1_bus6;I2_bus6];

Sl_update_bus6 = Vabc_update_bus6.*conj(Itemp_new_bus6);
Pl_update_bus6 = real(Sl_update_bus6);
Ql_update_bus6 = imag(Sl_update_bus6);


I0_bus5=Itemp_bus5(1,1);
I1_bus5=conj(Sseq1_update_bus5/(3*Vabc_seq_bus5(2,1)));
I2_bus5=Itemp_bus5(3,1);

Itemp_new_bus5 = T2*[I0_bus5;I1_bus5;I2_bus5];

Sl_update_bus5 = Vabc_update_bus5.*conj(Itemp_new_bus5);
Pl_update_bus5 = real(Sl_update_bus5);
Ql_update_bus5 = imag(Sl_update_bus5);


I0_bus8=Itemp_bus8(1,1);
I1_bus8=conj(Sseq1_update_bus8/(3*Vabc_seq_bus8(2,1)));
I2_bus8=Itemp_bus8(3,1);

Itemp_new_bus8 = T2*[I0_bus8;I1_bus8;I2_bus8];

Sl_update_bus8 = Vabc_update_bus8.*conj(Itemp_new_bus8);
Pl_update_bus8 = real(Sl_update_bus8);
Ql_update_bus8 = imag(Sl_update_bus8);

%% Voltage phase update to the distribution system

V_asub_bus6=abs(Vabc_update_bus6(1));
V_bsub_bus6=abs(Vabc_update_bus6(2));
V_csub_bus6=abs(Vabc_update_bus6(3));

V_asub_bus5=abs(Vabc_update_bus5(1));
V_bsub_bus5=abs(Vabc_update_bus5(2));
V_csub_bus5=abs(Vabc_update_bus5(3));

V_asub_bus8=abs(Vabc_update_bus8(1));
V_bsub_bus8=abs(Vabc_update_bus8(2));
V_csub_bus8=abs(Vabc_update_bus8(3));

%% voltage angle update to the distribution system
V_aangsub_bus6=angle(Vabc_update_bus6(1))*180/pi;
V_bangsub_bus6=angle(Vabc_update_bus6(2))*180/pi;
V_cangsub_bus6=angle(Vabc_update_bus6(3))*180/pi;

V_aangsub_bus5=angle(Vabc_update_bus5(1))*180/pi;
V_bangsub_bus5=angle(Vabc_update_bus5(2))*180/pi;
V_cangsub_bus5=angle(Vabc_update_bus5(3))*180/pi;

V_aangsub_bus8=angle(Vabc_update_bus8(1))*180/pi;
V_bangsub_bus8=angle(Vabc_update_bus8(2))*180/pi;
V_cangsub_bus8=angle(Vabc_update_bus8(3))*180/pi;

%% output the cosim_iteration values and error
% disp('The residual vector is=');
% disp(R_matrix);   
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

