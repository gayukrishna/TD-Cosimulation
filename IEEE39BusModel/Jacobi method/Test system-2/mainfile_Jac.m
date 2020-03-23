clc;
clear all;

%% defining transformation matrices and their inverse
v=-0.5+0.8660*(sqrt(-1));
T1= (1/3)*[1 1 1; 1 v (v^2); 1 (v^2) v;];
T2=[1 1 1; 1 (v^2) v; 1 v (v^2);];
T3= (1/(sqrt(3)))*[1 1 1; 1 v (v^2); 1 (v^2) v;]; %Power invariant sequence transformation matrix

%% Initialization and variable definitions
Loads=[25; 26; 28; 32; 33; 34; 35; 36; 37; 38;];
Nload=length(Loads);
error_cosim=1; 
count_cosim=0;

Jacobian=cell(Nload,1);
J1=cell(Nload,1);
J2=cell(Nload,1);
J3=cell(Nload,1);
J4=cell(Nload,1);
Yth=cell(Nload,1);
% defining initial voltages, angles and power for distribution and
% transmission respectively
for k=1:Nload
    V_asub(k,1)=1.05;
    V_bsub(k,1)=1.05;
    V_csub(k,1)=1.05;
    V_aangsub(k,1)=0;
    V_bangsub(k,1)=240;
    V_cangsub(k,1)=120;
end

P_asub=0.1748;
P_bsub=0.1727;
P_csub=0.1736;
Q_asub=0.0398;
Q_bsub=0.0399;
Q_csub=0.0379;

Pl_update=[P_asub; P_bsub; P_csub];
Ql_update=[Q_asub; Q_bsub; Q_csub];
Sl_update=complex(Pl_update, Ql_update);

for k=1:Nload
Sl_aupdate(k)=Sl_update(1);
Sl_bupdate(k)=Sl_update(2);
Sl_cupdate(k)=Sl_update(3);
end

%% Initializing the T&D systems
%%% solving transmission load flow
[S_a,S_b,S_c,P_0seq,P_1seq,P_2seq,Q_0seq,Q_1seq,Q_2seq,VaT,VbT,VcT,V1T]=Trans_loadflow_multifeeder_Jac(Sl_aupdate,Sl_bupdate,Sl_cupdate);

V_a = abs(VaT)';
V_b = abs(VbT)';
V_c = abs(VcT)';
V_aang= angle(VaT)'.*180/pi;
V_bang= angle(VbT)'.*180/pi;
V_cang= angle(VcT)'.*180/pi;

%%% solving distribution load flow and %% finding thevenin equivalent
for x=1:Nload
    [S_phckt24sub,V_seqsub]=Dist_loadflow_multifeeder_Jac(x,T1,V_asub(x),V_bsub(x),V_csub(x),V_aangsub(x),V_bangsub(x),V_cangsub(x));
    S_ackt24sub(x,1)=S_phckt24sub(1);
    S_bckt24sub(x,1)=S_phckt24sub(2);
    S_cckt24sub(x,1)=S_phckt24sub(3);
    V_0sub(x,1)=V_seqsub(1);
    V_1sub(x,1)=V_seqsub(2);
    V_2sub(x,1)=V_seqsub(3);
    
    [Yth_f]=Finding_Jac();
    Yth{x}=Yth_f;
end

%% start of the simulation
time_start = tic;
while error_cosim>0.0001
%% solve transmission load flow
[S_a,S_b,S_c,P_0seq,P_1seq,P_2seq,Q_0seq,Q_1seq,Q_2seq,VaT,VbT,VcT,V0T,V1T,V2T]=Trans_loadflow_multifeeder_Jac(Sl_aupdate,Sl_bupdate,Sl_cupdate);

V_a = abs(VaT)';
V_b = abs(VbT)';
V_c = abs(VcT)';
V_aang= angle(VaT)'.*180/pi;
V_bang= angle(VbT)'.*180/pi;
V_cang= angle(VcT)'.*180/pi;
Sabc_T=cell(Nload,1);

for x=1:Nload
    Sabc_T{x}=[S_a(x); S_b(x); S_c(x)];
end

%% solving distribution load flow and %% finding jacobian
for x=1:Nload
    [S_phckt24sub,V_seqsub]=Dist_loadflow_multifeeder_Jac(x,T1,V_asub(x),V_bsub(x),V_csub(x),V_aangsub(x),V_bangsub(x),V_cangsub(x));
    S_ackt24sub(x,1)=S_phckt24sub(1);
    S_bckt24sub(x,1)=S_phckt24sub(2);
    S_cckt24sub(x,1)=S_phckt24sub(3);
    V_0sub(x,1)=V_seqsub(1);
    V_1sub(x,1)=V_seqsub(2);
    V_2sub(x,1)=V_seqsub(3);
end

for x=1:Nload
    [J1_mf,J2_mf,J3_mf,J4_mf]=Jac_calculation(v,Yth{x},V_asub(x),V_bsub(x),V_csub(x));
    J1{x}=J1_mf;
    J2{x}=J2_mf;
    J3{x}=J3_mf;
    J4{x}=J4_mf;
end

%% computing the difference and updates   

D_matrix=[real(S_a)'-real(S_ackt24sub);
           imag(S_a)'-imag(S_ackt24sub);
           real(S_b)'-real(S_bckt24sub);
           imag(S_b)'-imag(S_bckt24sub);
           real(S_c)'-real(S_cckt24sub);
           imag(S_c)'-imag(S_cckt24sub);
           V_0sub-V0T';
           V_1sub-V1T';
           V_2sub-V2T'];

%% Voltage and positive sequence power update to the transmission system

%%% voltage update
V_diff=cell(Nload,1);
Vabc_Dtemp=cell(Nload,1);
Vabc_update=cell(Nload,1);
for x=1:Nload
    V_diff{x} = inv(J4{x})*D_matrix([60+x 70+x 80+x]);  %% since V_0 of all distribtuion loads are in D_matrix 60's row, similarly for V_1 and V_2
end
for x=1:Nload
Vabc_Dtemp{x} = T2*([V_0sub(x,1); V_1sub(x,1); V_2sub(x,1)]);
Vabc_update{x} = Vabc_Dtemp{x}-V_diff{x};               %%% now updating the difference.
end

%%% sequence power update
Sl_update1=zeros(3*Nload,1);

for i = 1:20   %%% iterating the calculations of V_update and S_update with jacobians
Vabc_update1 = Vabc_update;

J1_new = blkdiag(J1{1},J1{2},J1{3},J1{4},J1{5},J1{6},J1{7},J1{8},J1{9},J1{10});
J2_new = blkdiag(J2{1},J2{2},J2{3},J2{4},J2{5},J2{6},J2{7},J2{8},J2{9},J2{10});

diffPQabc = eye(60,60)*(D_matrix(1:60)-J2_new*(abs(cell2mat(Vabc_update))-abs(cell2mat(Vabc_Dtemp))));

alpha = 0.8;
for x=1:Nload
Sa_update(x,1) = S_a(x)-alpha.*(diffPQabc(x,1)+1i*diffPQabc(10+x,1));
Sb_update(x,1) = S_b(x)-alpha.*(diffPQabc(20+x,1)+1i*diffPQabc(30+x,1));
Sc_update(x,1) = S_c(x)-alpha.*(diffPQabc(40+x,1)+1i*diffPQabc(50+x,1));
end

Sl_update= [Sa_update;Sb_update;Sc_update];
Sl_updC=cell(Nload,1);
for x=1:Nload
    Sl_updC{x}=[Sa_update(x,1);Sb_update(x,1);Sc_update(x,1)];
end
%% re updating the voltages using the full equation
% Z = 0.0141 + 0.0958i;
% V_upd=cell2mat(Vabc_update);
% y=1;
% for x=1:Nload
%     m1(x) = conj(Z/V_upd(y,1));
%     m2(x) = conj(Z/V_upd(y+1,1));
%     m3(x) = conj(Z/V_upd(y+2,1));
%     y=y+3;
% end

%Z = 0.0141 + 0.0958i;
V_upd=cell2mat(Vabc_update);
y=1;
for x=1:Nload
    m1(x) = conj(S_ackt24sub(x,1)/V_upd(y,1));
    m2(x) = conj(S_bckt24sub(x,1)/V_upd(y+1,1));
    m3(x) = conj(S_cckt24sub(x,1)/V_upd(y+2,1));
    y=y+3;
end

J3_new=cell(Nload,1);
for x=1:Nload
J3_new{x}= [m1(x),m2(x),m3(x); 
            m1(x),m2(x)./(v^2),m3(x)./(v);
            m1(x),m2(x)./(v),m3(x)./(v^2);];
end  

diff_Vabc=cell(Nload,1);
for x=1:Nload
    diff_Vabc{x}=inv(J4{x})*(D_matrix([60+x 70+x 80+x])+J3_new{x}*(Sl_updC{x}-Sabc_T{x}));
end

for x=1:Nload
Vabc_Dtemp{x} = T2*([V_0sub(x,1); V_1sub(x,1); V_2sub(x,1)]);
Vabc_update{x} = Vabc_Dtemp{x}-diff_Vabc{x};
end

if max(abs(Sl_update-Sl_update1))<0.00001
    break;
end

Sl_update1= Sl_update;
end

%% updates to the transmission system
%%% complex power updates to the transmission system
Sl_aupdate=Sa_update;
Sl_bupdate=Sb_update;
Sl_cupdate=Sc_update;

%% Updates to the distribution system
%%% Voltage phase update to the distribution system

V_update=cell2mat(Vabc_update);
y=1;
for k=1:Nload
    V_asub(k)=abs(V_update(y,1));
    V_bsub(k)=abs(V_update(y+1,1));
    V_csub(k)=abs(V_update(y+2,1));
    y=y+3;
end

%%% voltage angle update to the distribution system
y=1;
for k=1:Nload
    V_aangsub(k)=angle(V_update(y,1))*180/pi;
    V_bangsub(k)=angle(V_update(y+1,1))*180/pi;
    V_cangsub(k)=angle(V_update(y+2,1))*180/pi;
    y=y+3;
end

%% output the cosim_iteration values and error
% disp('The residual vector is=');
% disp(R_matrix);   
% abs(D_matrix)
error_cosim = max(abs(D_matrix)); 
count_cosim = count_cosim + 1;


disp('The co-simulation ran with ');
disp('Error_cosim=');
disp(error_cosim);
disp('iterations_cosim= ');
disp(count_cosim);
disp('..........................................................');

end
time_elapsed = toc(time_start);

V =[(V_asub) (V_bsub) (V_csub)];
S =[(Sl_aupdate) (Sl_bupdate) (Sl_cupdate)];

for k=1:Nload
[vunb,cunb]= unbal_cal(V(k,1), V(k,2), V(k,3), S(k,1), S(k,2), S(k,3));
Vunbal(k,1)=vunb;
Cunbal(k,1)=cunb;
end

disp ('output at buses')

disp ('The converged power values are [Pa Pb Pc],[Qa Qb, Qc]=')
disp ([real(S_ackt24sub) real(S_bckt24sub) real(S_bckt24sub)]);
disp ([imag(S_ackt24sub) imag(S_bckt24sub) imag(S_cckt24sub)]);

% disp('Positive sequence voltage at the interface [V1]=');
% disp(V1T);

disp('Phase Voltages at the interface [Va Vb Vc]=');
disp([V_asub V_bsub V_csub]);

disp('Phase angles at the interface =');
disp([V_aangsub V_bangsub V_cangsub]);

disp('The unbalances are [Vunb Cunb]=')
disp([Vunbal Cunbal]);

disp('..........................................................');

disp('Time taken for total simulation = ');
disp(time_elapsed)
disp('************************************************************');
