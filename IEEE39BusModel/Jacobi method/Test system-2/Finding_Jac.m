 
function[Yth]=Finding_Jac()
% - Start OpenDSS
DSSObj=actxserver('OpenDSSEngine.DSS');
if ~DSSObj.Start(0)
    disp('Unable to start openDSS Engine');
    return
end

DSSText=DSSObj.Text;
DSSCircuit=DSSObj.ActiveCircuit;
DSSBus = DSSCircuit.ActiveBus;

%DSSText.Command='Compile (C:\Users\user\Dropbox\personal\MATLABCODES\cosim_codes\Codes\bigger_model\Jacobi method\Test system-2\ckt24_multifeeder\master_ckt24_cosim.dss)';
DSSText.Command='Compile (C:\Users\user\Dropbox\personal\MATLABCODES\cosim_codes\Codes\bigger_model\Jacobi method\Test system-2\ckt24_unbalance\master_ckt24_cosim.dss)';
%DSSText.Command = 'Redirect PV_combined100per.txt';
%DSSText.command='BatchEdit Load..* yearly=default';
%DSSText.Command='solve loadmult=1.0';
DSSText.Command='solve mode=Faultstudy';
DSSCircuit.SetActiveBus('SOURCEBUS');
YscMatrix= DSSCircuit.ActiveBus.YscMatrix;

Yaa=(YscMatrix(1,1)+(YscMatrix(1,2)*sqrt(-1)));
Yab=(YscMatrix(1,3)+(YscMatrix(1,4)*sqrt(-1)));
Yac=(YscMatrix(1,5)+(YscMatrix(1,6)*sqrt(-1)));
Yba=(YscMatrix(1,7)+(YscMatrix(1,8)*sqrt(-1)));
Ybb=(YscMatrix(1,9)+(YscMatrix(1,10)*sqrt(-1)));
Ybc=(YscMatrix(1,11)+(YscMatrix(1,12)*sqrt(-1)));
Yca=(YscMatrix(1,13)+(YscMatrix(1,14)*sqrt(-1)));
Ycb=(YscMatrix(1,15)+(YscMatrix(1,16)*sqrt(-1)));
Ycc=(YscMatrix(1,17)+(YscMatrix(1,18)*sqrt(-1)));
Yth=[Yaa;Yab;Yac;Yba;Ybb;Ybc;Yca;Ycb;Ycc];


release(DSSObj)
