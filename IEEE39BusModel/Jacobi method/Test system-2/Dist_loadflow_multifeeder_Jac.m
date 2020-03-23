function [S_phckt24sub,V_seqsub]=Dist_loadflow_multifeeder_Jac(x,T1,V_asub,V_bsub,V_csub,V_aangsub,V_bangsub,V_cangsub)

Loads=[25; 26; 28; 32; 33; 34; 35; 36; 37; 38;];
Nload=length(Loads);

%% send transmission bus voltages to distribution substation 
string1 = 'New Circuit.ckt24 bus1=SourceBus.1 pu=';
string2='angle=';
string3='New Vsource.Source_2 bus1=SourceBus.2 pu=';
string4='angle=';
string5='New Vsource.Source_3 bus1=SourceBus.3  pu=';
string6='angle=';
string7 = 'basekV=132.79 phase=1 R1=0.63 X1=6.72 R0=4.07 X0=15.55';

%fid=fopen('C:\Users\user\Dropbox\personal\MATLABCODES\cosim_codes\Codes\bigger_model\Jacobi method\Test system-2\ckt24_multifeeder\Interfacevoltage.txt','w');
fid=fopen('C:\Users\user\Dropbox\personal\MATLABCODES\cosim_codes\Codes\bigger_model\Jacobi method\Test system-2\ckt24_unbalance\Interfacevoltage.txt','w');
fprintf(fid, '%s %f %s %f %s \r\n', string1, V_asub, string2, V_aangsub, string7, string3, V_bsub, string4, V_bangsub, string7, string5, V_csub, string6, V_cangsub, string7);
fclose(fid);

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

%DSSText.Command='Compile (C:\Users\user\Dropbox\personal\MATLABCODES\cosim_codes\Codes\bigger_model\Jacobi method\Test system-2\ckt24_multifeeder\master_ckt24_cosim.dss)';
DSSText.Command='Compile (C:\Users\user\Dropbox\personal\MATLABCODES\cosim_codes\Codes\bigger_model\Jacobi method\Test system-2\ckt24_unbalance\master_ckt24_cosim.dss)';
DSSText.Command='solve mode=snap loadmult = 1';
    
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

S_ackt24sub=complex(P_ackt24sub,Q_ackt24sub);
S_bckt24sub=complex(P_bckt24sub,Q_bckt24sub);
S_cckt24sub=complex(P_cckt24sub,Q_cckt24sub);

S_phckt24sub=[S_ackt24sub; S_bckt24sub; S_cckt24sub];
P_phckt24sub=[P_ackt24sub; P_bckt24sub; P_cckt24sub];
Q_phckt24sub=[Q_ackt24sub; Q_bckt24sub; Q_cckt24sub];
 
release(DSSObj)

