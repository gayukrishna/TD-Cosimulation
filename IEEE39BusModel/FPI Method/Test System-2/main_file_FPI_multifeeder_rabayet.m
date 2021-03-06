clc;
clear all;

%% defining transformation matrices and their inverse
v=-0.5+0.8660*(sqrt(-1));
T1= (1/3)*[1 1 1; 1 v (v^2); 1 (v^2) v;];
T2=[1 1 1; 1 (v^2) v; 1 v (v^2);];
T3= (1/(sqrt(3)))*[1 1 1; 1 v (v^2); 1 (v^2) v;]; %Power invariant sequence transformation matrix

%% start of the entire simulation
Loads=[25; 26; 28; 32; 33; 34; 35; 36; 37; 38;];
Nload=length(Loads);
error_cosim=1;
count_cosim=0;
% defining initial voltages, angles and power for distribution and
% transmission respectively

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

time_start = tic;

%% solve transmission load flow to initialize -
[S_a,S_b,S_c,VaT,VbT,VcT,V1T]=Trans_loadflow_multifeeder(Sl_aupdate,Sl_bupdate,Sl_cupdate);

V_a = abs(VaT)';
V_b = abs(VbT)';
V_c = abs(VcT)';
V_aang= angle(VaT)'.*180/pi;
V_bang= angle(VbT)'.*180/pi;
V_cang= angle(VcT)'.*180/pi;


while error_cosim > 0.001
    fclose('all');
    %% solve transmission load flow
    
    %giving input to distribution system-
    
    V_asub = V_a;
    V_bsub = V_b;
    V_csub = V_c;
    V_aangsub = V_aang;
    V_bangsub = V_bang;
    V_cangsub = V_cang;
    
    %% solving distribution load flow
    for x=1:Nload
        [S_phckt24sub]=Dist_loadflow_multifeeder(x,T1,V_asub(x),V_bsub(x),V_csub(x),V_aangsub(x),V_bangsub(x),V_cangsub(x));
        S_ackt24sub(x,1)=S_phckt24sub(1);
        S_bckt24sub(x,1)=S_phckt24sub(2);
        S_cckt24sub(x,1)=S_phckt24sub(3);
    end
    
    %% give input updates for next iteration in cosim
    Sl_aupdate = S_ackt24sub;
    Sl_bupdate = S_bckt24sub;
    Sl_cupdate = S_cckt24sub;
    
    
    %% solve transmission load flow
    [S_a,S_b,S_c,VaT,VbT,VcT,V1T]=Trans_loadflow_multifeeder(Sl_aupdate,Sl_bupdate,Sl_cupdate);
    
    V_a = abs(VaT)';
    V_b = abs(VbT)';
    V_c = abs(VcT)';
    V_aang= angle(VaT)'.*180/pi;
    V_bang= angle(VbT)'.*180/pi;
    V_cang= angle(VcT)'.*180/pi;
    
    %% computing the difference and updates
    R_matrix=[real(S_ackt24sub)-real(S_a)';
        imag(S_ackt24sub)-imag(S_a)';
        real(S_bckt24sub)-real(S_b)';
        imag(S_bckt24sub)-imag(S_b)';
        real(S_cckt24sub)-real(S_c)';
        imag(S_cckt24sub)-imag(S_c)';
        V_a-V_asub;
        V_b-V_bsub;
        V_c-V_csub;];
    
    V =[(V_asub) (V_bsub) (V_csub)];
    S =[(Sl_aupdate) (Sl_bupdate) (Sl_cupdate)];
    
    for k=1:Nload
        [vunb,cunb]= unbal_cal(V(k,1), V(k,2), V(k,3), S(k,1), S(k,2), S(k,3));
        Vunbal(k,1)=vunb;
        Cunbal(k,1)=cunb;
    end
    
    %% output the cosim_iteration values and error
    
    error_cosim = max(abs(R_matrix));
    count_cosim = count_cosim + 1;
    
    disp('The co-simulation ran with ');
    disp('Error_cosim=');
    disp(error_cosim);
    
    disp('..........................................................');
end

time_elapsed = toc(time_start);
% disp ('output at buses')
%
% disp ('The converged power values are [Pa Pb Pc]=')
% disp ([real(S_ackt24sub) real(S_bckt24sub) real(S_bckt24sub)]);
% disp ([imag(S_ackt24sub) imag(S_bckt24sub) imag(S_cckt24sub)]);
%
% disp('Positive sequence voltage at the interface [V1]=');
% disp(V1T);
%
% disp('Phase Voltages at the interface [Va Vb Vc]=');
% disp([V_asub V_bsub V_csub]);
%
% disp('Phase angles at the interface =');
% disp([V_aangsub V_bangsub V_cangsub]);
%
% disp('The unbalances are=')
% disp([Vunbal Cunbal]);

disp('************************************************************');

disp('iterations_cosim= ');
disp(count_cosim);
disp('Time taken for total simulation = ');
disp(time_elapsed)
disp('************************************************************');
