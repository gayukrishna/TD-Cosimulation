clc;
clear all;

%% defining transformation matrices and their inverse
v=-0.5+0.8660*(sqrt(-1));
T1= (1/3)*[1 1 1; 1 v (v^2); 1 (v^2) v;];
T2=[1 1 1; 1 (v^2) v; 1 v (v^2);];
T3= (1/(sqrt(3)))*[1 1 1; 1 v (v^2); 1 (v^2) v;]; %Power invariant sequence transformation matrix

%% start of the entire simulation
a=[];
r=[];
error_cosim=1; 
count_cosim=0;
V_asub=1.05;
V_bsub=1.05;
V_csub=1.05;
V_aangsub=0;
V_bangsub=240;
V_cangsub=120;
Pl_update=[0.1748; 0.1727; 0.1736];
Ql_update=[0.0398; 0.0399; 0.0379];
% Pl_update=[0.25; 0.55; 0.34];
% Ql_update=[0.27; 0.11; 0.15];

while error_cosim>0.0001
%% solve transmission load flow
[P_bus6,Q_bus6,S_bus6,P_seqbus6,Q_seqbus6,V_abus6,V_bbus6,V_cbus6,V_0bus6,V_1bus6,V_2bus6,V_aangbus6,V_bangbus6,V_cangbus6]=Trans_loadflow(Pl_update,Ql_update);

%% solving distribution load flow
[P_seqckt24sub,Q_seqckt24sub,S_seqckt24sub,S_phckt24sub,P_phckt24sub,Q_phckt24sub,V_0sub,V_1sub,V_2sub,V_ackt24sub,V_bckt24sub,V_cckt24sub,Vseq_sub]=Dist_loadflow(T1,T3,V_asub,V_bsub,V_csub,V_aangsub,V_bangsub,V_cangsub);

%% computing the difference and updates

R_matrix=[P_phckt24sub-P_bus6
          Q_phckt24sub-Q_bus6;
          V_abus6-V_asub;
          V_bbus6-V_bsub;
          V_cbus6-V_csub;];
      
%% give input updates for next iteration in cosim

P_diff=P_phckt24sub-P_bus6;
Q_diff=Q_phckt24sub-Q_bus6;
Pl_update=P_bus6+P_diff;
Ql_update=Q_bus6+Q_diff;
V_asub=V_asub+(V_abus6-V_asub);
V_bsub=V_bsub+(V_bbus6-V_bsub);
V_csub=V_csub+(V_cbus6-V_csub);
V_aangsub=V_aangbus6;
V_bangsub=V_bangbus6;
V_cangsub=V_cangbus6;

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

%a=[a;V_asub V_bsub V_csub V_aangsub V_bangsub V_cangsub P_phckt24sub(1,1) P_phckt24sub(2,1) P_phckt24sub(3,1) Q_phckt24sub(1,1) Q_phckt24sub(2,1) Q_phckt24sub(3,1)]
r=abs([r;R_matrix(1,1) R_matrix(2,1) R_matrix(3,1) R_matrix(4,1) R_matrix(5,1) R_matrix(6,1) R_matrix(7,1) R_matrix(8,1) R_matrix(9,1)]) 
disp('..........................................................');
end

% Vmag.a=V_abus6;
% Vmag.b=V_bbus6;
% Vmag.c=V_cbus6;
% Vang.a=V_aangbus6;
% Vang.b=V_bangbus6;
% Vang.c=V_cangbus6;
% Power.a=P_bus6(1,1);
% Power.b=P_bus6(2,1);
% Power.c=P_bus6(3,1);
% QPower.a=Q_bus6(1,1);
% QPower.b=Q_bus6(2,1);
% QPower.c=Q_bus6(3,1);
% save('Vgemagnitude','-struct','Vmag');
% save('Vgeangle','-struct','Vang');
% save('activepower','-struct','Power');
% save('reactivepower','-struct','QPower');

disp ('The converged power values are is=')
disp (P_phckt24sub);
disp (Q_phckt24sub);
  
disp('Phase Voltages at the interface=');
disp(V_asub);
disp(V_bsub);
disp(V_csub);

disp('Phase angles at the interface=');
disp(V_aangsub);
disp(V_bangsub);
disp(V_cangsub);

