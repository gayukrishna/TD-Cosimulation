function [S_phckt24sub,V_seqsub]=Dist_loadflow_multifeeder_Jac(Cktnum,T1,V_asub,V_bsub,V_csub,V_aangsub,V_bangsub,V_cangsub)
%% send bus 6 voltages to distribution substation 
string1 = 'New Circuit.ckt24 bus1=SourceBus.1 pu=';
string2='angle=';
string3='New Vsource.Source_2 bus1=SourceBus.2 pu=';
string4='angle=';
string5='New Vsource.Source_3 bus1=SourceBus.3  pu=';
string6='angle=';
string7 = 'basekV=132.79 phase=1 R1=0.63 X1=6.72 R0=4.07 X0=15.55';
%string4 ='basekV=66.395 phase=1 r1=0  x1=0.001  r0=0  x0=0.001';

if Cktnum == 6
fid=fopen('C:\Users\user\Dropbox\personal\MATLABCODES\cosim_codes\Codes\Jacobi method\Test system-2\Ckt24\Interfacevoltage_ckt6.txt','w');
fprintf(fid, '%s %f %s %f %s \r\n', string1, V_asub, string2, V_aangsub, string7, string3, V_bsub, string4, V_bangsub, string7, string5, V_csub, string6, V_cangsub, string7);
%fprintf(fid, '%s %f %s \r\n', string1, V_asub, string4, string2, V_bsub, string4, string3, V_csub, string4);
fclose(fid);
elseif Cktnum == 5
    fid=fopen('C:\Users\user\Dropbox\personal\MATLABCODES\cosim_codes\Codes\Jacobi method\Test system-2\Ckt24\Interfacevoltage_ckt5.txt','w');
fprintf(fid, '%s %f %s %f %s \r\n', string1, V_asub, string2, V_aangsub, string7, string3, V_bsub, string4, V_bangsub, string7, string5, V_csub, string6, V_cangsub, string7);
%fprintf(fid, '%s %f %s \r\n', string1, V_asub, string4, string2, V_bsub, string4, string3, V_csub, string4);
fclose(fid);
else
    fid=fopen('C:\Users\user\Dropbox\personal\MATLABCODES\cosim_codes\Codes\Jacobi method\Test system-2\Ckt24\Interfacevoltage_ckt8.txt','w');
fprintf(fid, '%s %f %s %f %s \r\n', string1, V_asub, string2, V_aangsub, string7, string3, V_bsub, string4, V_bangsub, string7, string5, V_csub, string6, V_cangsub, string7);
%fprintf(fid, '%s %f %s \r\n', string1, V_asub, string4, string2, V_bsub, string4, string3, V_csub, string4);
fclose(fid);
end
% - Start OpenDSS
DSSObj=actxserver('OpenDSSEngine.DSS');
if ~DSSObj.Start(0)
    disp('Unable to start openDSS Engine');
    return
end

DSSText=DSSObj.Text;
DSSCircuit=DSSObj.ActiveCircuit;
DSSBus = DSSCircuit.ActiveBus;

%% run load flow in snap shot mode

if Cktnum == 6
DSSText.Command='Compile (C:\Users\user\Dropbox\personal\MATLABCODES\cosim_codes\Codes\Jacobi method\Test system-2\Ckt24\master_ckt24_cosim_ckt6.dss)';
DSSText.Command='solve mode=snap loadmult = 2.5';
elseif Cktnum == 5
DSSText.Command='Compile (C:\Users\user\Dropbox\personal\MATLABCODES\cosim_codes\Codes\Jacobi method\Test system-2\Ckt24\master_ckt24_cosim_ckt5.dss)';
DSSText.Command='solve mode=snap loadmult = 2.5';
else
DSSText.Command='Compile (C:\Users\user\Dropbox\personal\MATLABCODES\cosim_codes\Codes\Jacobi method\Test system-2\Ckt24\master_ckt24_cosim_ckt8.dss)';
DSSText.Command='solve mode=snap loadmult = 2.5';
end
    
DSSCircuit.SetActiveBus('SOURCEBUS');
DSSCircuit.SetActiveElement('Transformer.SUBXFMR');

%% calculating input sequence voltages
V_asub = complex(V_asub*cos(V_aangsub*pi/180), V_asub*sin(V_aangsub*pi/180));
V_bsub = complex(V_bsub*cos(V_bangsub*pi/180), V_bsub*sin(V_bangsub*pi/180));
V_csub = complex(V_csub*cos(V_cangsub*pi/180), V_csub*sin(V_cangsub*pi/180));

V_seqsub=T1*[V_asub; V_bsub; V_csub;];

V_0sub=V_seqsub(1,1);
V_1sub=V_seqsub(2,1);
V_2sub=V_seqsub(3,1);
%% output 
%Power = DSSCircuit.TotalPower;
Voltage = DSSBus.puVmagAngle;
Current=DSSCircuit.ActiveCktElement.Currents;
Power=DSSCircuit.ActiveCktElement.Powers;

%% phase powers from substation
P_ackt24sub = Power(1,1)/(1000*100);
Q_ackt24sub = Power(1,2)/(1000*100);
P_bckt24sub = Power(1,3)/(1000*100);
Q_bckt24sub = Power(1,4)/(1000*100);
P_cckt24sub = Power(1,5)/(1000*100);
Q_cckt24sub = Power(1,6)/(1000*100);

S_ackt24sub=complex(Power(1,1),Power(1,2))/(1000*100);
S_bckt24sub=complex(Power(1,3),Power(1,4))/(1000*100);
S_cckt24sub=complex(Power(1,5),Power(1,6))/(1000*100);

S_phckt24sub=[S_ackt24sub; S_bckt24sub; S_cckt24sub];
P_phckt24sub=[P_ackt24sub; P_bckt24sub; P_cckt24sub];
Q_phckt24sub=[Q_ackt24sub; Q_bckt24sub; Q_cckt24sub];

% %% calculating sequence voltages from phase
% 
% [a1, b1]=pol2cart((pi*Voltage(1,2))/180,Voltage(1,1));
%    V_ackt24sub = complex(a1,b1);
% [a2, b2]=pol2cart((pi*Voltage(1,4))/180,Voltage(1,3));
%    V_bckt24sub = complex(a2,b2);    
% [a3, b3]=pol2cart((pi*Voltage(1,6))/180,Voltage(1,5));
%    V_cckt24sub = complex(a3,b3);     
% Vseq_sub=T1*[V_ackt24sub; V_bckt24sub; V_cckt24sub];
%  
% %% calculating seq currents from phase
% 
% I_a=complex(Current(1,1),Current(1,2));
% I_b=complex(Current(1,3),Current(1,4));
% I_c=complex(Current(1,5),Current(1,6));
% [c1, d1]=cart2pol(Current(1,1),Current(1,2)); 
% [c2, d2]=cart2pol(Current(1,3),Current(1,4));
% [c3, d3]=cart2pol(Current(1,5),Current(1,6));
% I_base=135.5524;
% MVA_base=18.000;
% % I_base=205.78;
% % MVA_base=2.566;
% d1=d1/I_base;
% d2=d2/I_base;
% d3=d3/I_base;
% [h1, g1]=pol2cart(c1,d1); 
% I_ackt24sub = complex(h1,g1);
% [h2, g2]=pol2cart(c2,d2); 
% I_bckt24sub = complex(h2,g2);
% [h3, g3]=pol2cart(c3,d3); 
% I_cckt24sub = complex(h3,g3);
% Iseq_sub=T1*[I_ackt24sub; I_bckt24sub; I_cckt24sub];
% 
% %% calculating seq powers using calculated Vseq and Iseq
% 
% S_0ckt24sub=3*Vseq_sub(1,1)*conj(Iseq_sub(1,1));
% S_1ckt24sub=3*Vseq_sub(2,1)*conj(Iseq_sub(2,1));
% S_2ckt24sub=3*Vseq_sub(3,1)*conj(Iseq_sub(3,1));
% 
% S_seqckt24sub=[S_0ckt24sub; S_1ckt24sub; S_2ckt24sub;];
% 
% P_1ckt24sub=(real(S_1ckt24sub)*MVA_base)/100;
% Q_1ckt24sub=(imag(S_1ckt24sub)*MVA_base)/100;
% P_0ckt24sub=(real(S_0ckt24sub)*MVA_base)/100;
% Q_0ckt24sub=(imag(S_0ckt24sub)*MVA_base)/100;
% P_2ckt24sub=(real(S_2ckt24sub)*MVA_base)/100;
% Q_2ckt24sub=(imag(S_2ckt24sub)*MVA_base)/100;
% 
% P_seqckt24sub=[P_0ckt24sub; P_1ckt24sub; P_2ckt24sub;];
% Q_seqckt24sub=[Q_0ckt24sub; Q_1ckt24sub; Q_2ckt24sub;];

% disp('Phase Voltages given to the distribution substation=');
% disp(abs(V_asub));
% disp(abs(V_bsub));
% disp(abs(V_csub));

% disp('Sequence Voltages given to the distribution substation=');
% disp(V_0sub);
% disp(V_1sub);
% disp(V_2sub);

% disp('Phase Voltages from the distribution substation=');
% disp(abs(V_ackt24sub));
% disp(abs(V_bckt24sub));
% disp(abs(V_cckt24sub));

% disp('Sequence Voltages from the distribution substation=');
% disp(abs(Vseq_sub));
% 
% disp('phase powers from the distribution substation=');
%  disp(P_phckt24sub);
%  disp(Q_phckt24sub);
% % 
%  disp('sequence powers from the distribution substation=');
%  disp(P_seqckt24sub);
%  disp(Q_seqckt24sub);
 
release(DSSObj)

