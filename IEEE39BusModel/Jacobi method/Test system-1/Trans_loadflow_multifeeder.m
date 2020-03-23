function [S_bus26,VaT26,VbT26,VcT26,V1T26,S_bus25,VaT25,VbT25,VcT25,V1T25,S_bus28,VaT28,VbT28,VcT28,V1T28]=Trans_loadflow_multifeeder(Sl_bus25_update,Sl_bus26_update,Sl_bus28_update)
%% data input

Busdata=[
% Busno BaseKV Type  Vsp    theta       Pg     Qg       Pl       Ql     Qmax   Qmin    
1	3	0.982	0	521.207	198.319     9.2     4.6     9900	-9999
2	2	0.9831	0	650     205.154     0       0       9900	-9999
3	2	0.9972	0	632     109.911     0       0       9900	-9999
4	2	1.0123	0	508     165.765     0       0       9900	-9999
5	2	1.0493	0	650     212.417     0       0       9900	-9999
6	2	1.0635	0	560     101.178     0       0       9900	-9999
7	2	1.0278	0	540     0.445       0       0       9900	-9999
8	2	1.0265	0	830     22.844      0       0       9900	-9999
9	2	1.03	0	1000	88.281      1104	250     9900	-9999
10	2	1.0475	0	250     146.154     0       0       9900	-9999
11	1	1       0	0       0           0       0       0           0
12	1	1       0	0       0           0       0       0           0
13	1	1       0	0       0           322     2.4     0           0
14	1	1       0	0       0           500     184     0           0
15	1	1       0	0       0           0       0       0           0
16	1	1       0	0       0           0       0       0           0
17	1	1       0	0       0           233.8	84      0           0
18	1	1       0	0       0           522     176     0           0
19	1	1       0	0       0           0       0       0           0
20	1	1       0	0       0           0       0       0           0
21	1	1       0	0       0           0       0       0           0
22	1	1       0	0       0           7.5     88      0           0
23	1	1       0	0       0           0       0       0           0
24	1	1       0	0       0           0       0       0           0
25	1	1       0	0       0           320     153     0           0
26	1	1       0	0       0           329.4	32.3	0           0
27	1	1       0	0       0           0       0       0           0
28	1	1       0	0       0           158     30      0           0
29	1	1       0	0       0           0       0       0           0
30	1	1       0	0       0           628     103     0           0
31	1	1       0	0       0           274     115     0           0
32	1	1       0	0       0           0       0       0           0
33	1	1       0	0       0           247.5	84.6	0           0
34	1	1       0	0       0           308.6	-92.2	0           0
35	1	1       0	0       0           224     47.2	0           0
36	1	1       0	0       0           139     17      0           0
37	1	1       0	0       0           281     75.5	0           0
38	1	1       0	0       0           206     27.6	0           0
39	1	1       0	0       0           283.5	26.9	0           0];


Branchdata=[
%From  To    R0      X0          B0         R1          X1          B1          R2          X2          B2      Tap  baseMVA    NominalVoltage  branch 
1	16	0	0.025	0	0	0.025	0	0	0.025	0	1	100	22	0
12	10	0	0.0181	0	0	0.0181	0	0	0.0181	0	1.025	100	22	0
20	2	0	0.02	0	0	0.02	0	0	0.02	0	1.07	100	22	0
22	21	0.0016	0.0435	0	0.0016	0.0435	0	0.0016	0.0435	0	1.006	100	345	0
22	23	0.0016	0.0435	0	0.0016	0.0435	0	0.0016	0.0435	0	1.006	100	345	0
29	3	0.0007	0.0142	0	0.0007	0.0142	0	0.0007	0.0142	0	1.07	100	22	0
29	30	0.0007	0.0138	0	0.0007	0.0138	0	0.0007	0.0138	0	1.06	100	345	0
30	4	0.0009	0.018	0	0.0009	0.018	0	0.0009	0.018	0	1.009	100	22	0
32	5	0	0.0143	0	0	0.0143	0	0	0.0143	0	1.025	100	22	0
33	6	0.0005	0.0272	0	0.0005	0.0272	0	0.0005	0.0272	0	1	100	22	0
35	7	0.0006	0.0232	0	0.0006	0.0232	0	0.0006	0.0232	0	1.025	100	22	0
39	8	0.0008	0.0156	0	0.0008	0.0156	0	0.0008	0.0156	0	1.025	100	22	0
11	12	0.0035	0.0411	0.6987	0.0035	0.0411	0.6987	0.0035	0.0411	0.6987	0	100	345	1
11	9	0.001	0.025	0.75	0.001	0.025	0.75	0.001	0.025	0.75	0	100	345	1
12	13	0.0013	0.0151	0.2572	0.0013	0.0151	0.2572	0.0013	0.0151	0.2572	0	100	345	1
12	35	0.007	0.0086	0.146	0.007	0.0086	0.146	0.007	0.0086	0.146	0	100	345	1
13	14	0.0013	0.0213	0.2214	0.0013	0.0213	0.2214	0.0013	0.0213	0.2214	0	100	345	1
13	28	0.0011	0.0133	0.2138	0.0011	0.0133	0.2138	0.0011	0.0133	0.2138	0	100	345	1
14	15	0.0008	0.0128	0.1342	0.0008	0.0128	0.1342	0.0008	0.0128	0.1342	0	100	345	1
14	24	0.0008	0.0129	0.1382	0.0008	0.0129	0.1382	0.0008	0.0129	0.1382	0	100	345	1
15	18	0.0008	0.0112	0.1476	0.0008	0.0112	0.1476	0.0008	0.0112	0.1476	0	100	345	1
16	15	0.0002	0.0026	0.0434	0.0002	0.0026	0.0434	0.0002	0.0026	0.0434	0	100	345	1
16	17	0.0006	0.0092	0.113	0.0006	0.0092	0.113	0.0006	0.0092	0.113	0	100	345	1
16	21	0.0007	0.0082	0.1389	0.0007	0.0082	0.1389	0.0007	0.0082	0.1389	0	100	345	1
17	18	0.0004	0.0046	0.078	0.0004	0.0046	0.078	0.0004	0.0046	0.078	0	100	345	1
18	19	0.0023	0.0363	0.3804	0.0023	0.0363	0.3804	0.0023	0.0363	0.3804	0	100	345	1
19	9	0.001	0.025	1.2	0.001	0.025	1.2	0.001	0.025	1.2	0	100	345	1
20	21	0.0004	0.0043	0.0729	0.0004	0.0043	0.0729	0.0004	0.0043	0.0729	0	100	345	1
20	23	0.0004	0.0043	0.0729	0.0004	0.0043	0.0729	0.0004	0.0043	0.0729	0	100	345	1
23	24	0.0009	0.0101	0.1723	0.0009	0.0101	0.1723	0.0009	0.0101	0.1723	0	100	345	1
24	25	0.0018	0.0217	0.366	0.0018	0.0217	0.366	0.0018	0.0217	0.366	0	100	345	1
25	26	0.0009	0.0094	0.171	0.0009	0.0094	0.171	0.0009	0.0094	0.171	0	100	345	1
26	27	0.0007	0.0089	0.1342	0.0007	0.0089	0.1342	0.0007	0.0089	0.1342	0	100	345	1
26	29	0.0016	0.0195	0.304	0.0016	0.0195	0.304	0.0016	0.0195	0.304	0	100	345	1
26	31	0.0008	0.0135	0.2548	0.0008	0.0135	0.2548	0.0008	0.0135	0.2548	0	100	345	1
26	34	0.0003	0.0059	0.068	0.0003	0.0059	0.068	0.0003	0.0059	0.068	0	100	345	1
27	28	0.0007	0.0082	0.1319	0.0007	0.0082	0.1319	0.0007	0.0082	0.1319	0	100	345	1
27	37	0.0013	0.0173	0.3216	0.0013	0.0173	0.3216	0.0013	0.0173	0.3216	0	100	345	1
31	32	0.0008	0.014	0.2565	0.0008	0.014	0.2565	0.0008	0.014	0.2565	0	100	345	1
32	33	0.0006	0.0096	0.1846	0.0006	0.0096	0.1846	0.0006	0.0096	0.1846	0	100	345	1
33	34	0.0022	0.035	0.361	0.0022	0.035	0.361	0.0022	0.035	0.361	0	100	345	1
35	36	0.0032	0.0323	0.513	0.0032	0.0323	0.513	0.0032	0.0323	0.513	0	100	345	1
36	37	0.0014	0.0147	0.2396	0.0014	0.0147	0.2396	0.0014	0.0147	0.2396	0	100	345	1
36	38	0.0043	0.0474	0.7802	0.0043	0.0474	0.7802	0.0043	0.0474	0.7802	0	100	345	1
36	39	0.0057	0.0625	1.029	0.0057	0.0625	1.029	0.0057	0.0625	1.029	0	100	345	1
38	39	0.0014	0.0151	0.249	0.0014	0.0151	0.249	0.0014	0.0151	0.249	0	100	345	1];
    
baseMVA = 100;
%% extracting data from the input file

Nbus = max(Busdata(:,1));                   %Number of buses
NGenbus = (length(find(Busdata(:,2)==2))) + (length(find(Busdata(:,2)==3)));
NPQbus = (length(find(Busdata(:,2)==1)));
busPQ = (find(Busdata(:,2)==1));            % PQ bus array
busPV= find(Busdata(:,2)==2);
Slackbus = find(Busdata(:,2)==3);
V = Busdata(:,3);                           % voltage magnitude
theta_1 = Busdata(:,4)*pi/180;                % voltage angle in radians
Pg=Busdata(:,5)./baseMVA;                   % P power generated in p.u
Qg=Busdata(:,6)./baseMVA;                   % Q power generated in p.u
Pd=Busdata(:,7)./(baseMVA*3);                   % P load in p.u
Qd=Busdata(:,8)./(baseMVA*3);                   % Q load in p.u

frombus = Branchdata(:,1);                  % branch from
tobus = Branchdata(:,2);                    % branch to
Nlines=length(find(Branchdata(:,15)==1));   % number of transmission lines (% 1-line, 0-transformer)
Nbranch =length(frombus);                   % Number of branches
R0=Branchdata(:,3);
X0=Branchdata(:,4);
B0=sqrt(-1)*Branchdata(:,5);
R1=Branchdata(:,6);
X1=Branchdata(:,7);
B1=sqrt(-1)*Branchdata(:,8);
R2=Branchdata(:,9);
X2=Branchdata(:,10);
B2=sqrt(-1)*Branchdata(:,11);

z0=R0+X0*sqrt(-1);
z1=R1+X1*sqrt(-1);
z2=R2+X2*sqrt(-1);
y0= 1./z0; 
y1= 1./z1;  
y2= 1./z2;  
Y0= zeros(Nbus,Nbus); 
Y1= zeros(Nbus,Nbus); 
Y2= zeros(Nbus,Nbus); 

%making tap at each branch. 1 means no effect of tap
tap=Branchdata(:,12);
for i=1:Nbranch
    if tap(i)==0
        tap(i)=1;
    end
end

%% formation of Ybus
% Y0 -  zero sequence nodal admittance matrix                 

% %Formation of the off diagonal elements of the bus admittance matrix
% for i= 1:Nbranch
%     Y0(frombus(i),tobus(i)) = Y0(frombus(i),tobus(i))-y0(i)/tap(i);
%     Y0(tobus(i),frombus(i)) = Y0(frombus(i),tobus(i));
% end

%Formation of the diagonal elements of the matrix
for i= 1:Nbranch
    for j = 1:Nbus
        if frombus(i) == j
            Y0(j,j) = Y0(j,j) + y0(i)/((tap(i))^2) + (B0(i)/2);
        end
            
    end
end

for i= 1:Nbranch
    for j = 1:Nbus
        if(tobus(i) == j)
            Y0(j,j) = Y0(j,j) + y0(i) + (B0(i)/2);
        end
            
    end
end

% Y1 -  positive sequence nodal admittance matrix 

%Formation of the off diagonal elements of the bus admittance matrix
for i= 1:Nbranch
    Y1(frombus(i),tobus(i)) = Y1(frombus(i),tobus(i))-y1(i)/tap(i);
    Y1(tobus(i),frombus(i)) = Y1(frombus(i),tobus(i));
end

%Formation of the diagonal elements of the matrix
for i= 1:Nbranch
    for j = 1:Nbus
        if frombus(i) == j
            Y1(j,j) = Y1(j,j) + y1(i)/((tap(i))^2) + (B1(i)/2);
        end
            
    end
end

for i= 1:Nbranch
    for j = 1:Nbus
        if(tobus(i) == j)
            Y1(j,j) = Y1(j,j) + y1(i) + (B1(i)/2);
        end
            
    end
end

% Y2 -  negative sequence nodal admittance matrix  

% Formation of the off diagonal elements of the bus admittance matrix
% for i= 1:Nbranch
%     Y2(frombus(i),tobus(i)) = Y2(frombus(i),tobus(i))-y2(i)/tap(i);
%     Y2(tobus(i),frombus(i)) = Y2(frombus(i),tobus(i));
% end

%Formation of the diagonal elements of the matrix
for i= 1:Nbranch
    for j = 1:Nbus
        if frombus(i) == j
            Y2(j,j) = Y2(j,j) + y2(i)/((tap(i))^2) + (B2(i)/2);
        end
            
    end
end

for i= 1:Nbranch
    for j = 1:Nbus
        if(tobus(i) == j)
            Y2(j,j) = Y2(j,j) + y2(i) + (B2(i)/2);
        end
            
    end
end
%% defining transformation matrices and their inverse

v=-0.5+0.8660*(sqrt(-1));
T1= (1/3)*[1 1 1; 1 v (v^2); 1 (v^2) v;];
T2=[1 1 1; 1 (v^2) v; 1 v (v^2);];
T3= (1/(sqrt(3)))*[1 1 1; 1 v (v^2); 1 (v^2) v;]; %Power invariant sequence transformation matrix
T4=(1/(sqrt(3)))*[1 1 1; 1 (v^2) v; 1 v (v^2);];  %Power invariant sequence transformation matrix

%% initial voltages in three phases
Va=V.*(v^3);
Vb=V.*(v^2);
Vc=V.*(v);

%% sequence voltages on all buses
for k=1:Nbus
Vseq(:,k)=T1*[Va(k,1); Vb(k,1);Vc(k,1);];
end
V_0=(Vseq(1,:)).';
V_1=(Vseq(2,:)).';
V_2=(Vseq(3,:)).';

for k=1:Nbus
Vseq_PV(:,k)=T3*[Va(k,1); Vb(k,1);Vc(k,1);];
end

V_0_PV=(Vseq_PV(1,:)).';
V_1_PV=(Vseq_PV(2,:)).';
V_2_PV=(Vseq_PV(3,:)).';

%% line series admittance data
Z_Seq=cell(Nbranch,1);  
Y_Seq=cell(Nbranch,1);

for ln=1:Nbranch
    if (Branchdata(ln,15)==1) %% transmission lines
        Z_Seq{ln}=  [z0(ln)  0  0; 0  z1(ln) 0; 0   0   z2(ln)];
        Y_Seq{ln} = inv(Z_Seq{ln});
    end 
end

%% line shunt admittance data
Y_Shunt=cell(Nbranch,1);

for ln=1:Nbranch
    if (Branchdata(ln,15)==1) %% transmission lines
        Y_Shunt{ln}=  [B0(ln)/2  0  0; 0  B1(ln)/2 0; 0   0   B2(ln)/2];
    end 
end

%% specified powers at load busbars
Pl_a=Pd; Ql_a=Qd;
Pl_b=Pd; Ql_b=Qd;
Pl_c=Pd; Ql_c=Qd;

%% update powers for next cosim iteration in all phases

Pl_a(25,1)=real(Sl_bus25_update(1,1));
Ql_a(25,1)=imag(Sl_bus25_update(1,1));
Pl_b(25,1)=real(Sl_bus25_update(2,1));
Ql_b(25,1)=imag(Sl_bus25_update(2,1));
Pl_c(25,1)=real(Sl_bus25_update(3,1));
Ql_c(25,1)=imag(Sl_bus25_update(3,1));
%%
Pl_a(26,1)=real(Sl_bus26_update(1,1));
Ql_a(26,1)=imag(Sl_bus26_update(1,1));
Pl_b(26,1)=real(Sl_bus26_update(2,1));
Ql_b(26,1)=imag(Sl_bus26_update(2,1));
Pl_c(26,1)=real(Sl_bus26_update(3,1));
Ql_c(26,1)=imag(Sl_bus26_update(3,1));
%%
Pl_a(28,1)=real(Sl_bus28_update(1,1));
Ql_a(28,1)=imag(Sl_bus28_update(1,1));
Pl_b(28,1)=real(Sl_bus28_update(2,1));
Ql_b(28,1)=imag(Sl_bus28_update(2,1));
Pl_c(28,1)=real(Sl_bus28_update(3,1));
Ql_c(28,1)=imag(Sl_bus28_update(3,1));

Ssp_a=Pl_a+Ql_a*(sqrt(-1));  
Ssp_b=Pl_b+Ql_b*(sqrt(-1));  
Ssp_c=Pl_c+Ql_c*(sqrt(-1));  

%% Convergence loop
error = 1;
count = 0;
theta_1=zeros(Nbus,1);

while error>0.00001
      
for t=1:Nbus
Iasp(t,1)=conj(Ssp_a(t,1)/Va(t,1));
Ibsp(t,1)=conj(Ssp_b(t,1)/Vb(t,1));
Icsp(t,1)=conj(Ssp_c(t,1)/Vc(t,1));
end
 Vabs_old = [abs(Va), abs(Vb),  abs(Vc)];
 
% sequence current injections
for h=1:Nbus
    Iseq(:,h)=T1*[Iasp(h,1); Ibsp(h,1); Icsp(h,1)];
end
for h=1:Nbus
    Iseq_PV(:,h)=T3*[Iasp(h,1); Ibsp(h,1); Icsp(h,1)];
end

Iload_0sp=(Iseq(1,:)).';
Iload_1=(Iseq(2,:)).';
Iload_2sp=(Iseq(3,:)).';

Iload_0sp_PV=(Iseq_PV(1,:)).';
Iload_1_PV=(Iseq_PV(2,:)).';
Iload_2sp_PV=(Iseq_PV(3,:)).';

if (count == 0)
    V_1_PV = (V_1_PV)*1;
    V_2_PV = (V_2_PV)*1;
    V_0_PV = (V_0_PV)*1;
else
    V_1_PV = (V_1)*sqrt(3);
    V_2_PV = (V_2)*sqrt(3);
    V_0_PV = (V_0)*sqrt(3);
end


for f=1:Nbus
Sload_1sp(f,1)=V_1_PV(f,1)*conj(Iload_1_PV(f,1));
Sload_0sp(f,1)=V_0_PV(f,1)*conj(Iload_0sp_PV(f,1));
Sload_2sp(f,1)=V_2_PV(f,1)*conj(Iload_2sp_PV(f,1));
end


Pload_1sp=real(Sload_1sp);
Qload_1sp=imag(Sload_1sp);
Pload_0sp=real(Sload_0sp);
Qload_0sp=imag(Sload_0sp);
Pload_2sp=real(Sload_2sp);
Qload_2sp=imag(Sload_2sp);

 %% off-diagonal current injection in series & shunt admittance matrix in untransposed lines
 
 dIz0=zeros(Nbranch,1);
 dIz1=zeros(Nbranch,1);
 dIz2=zeros(Nbranch,1);
 dIs0=zeros(Nbus,1);
 dIs1=zeros(Nbus,1);
 dIs2=zeros(Nbus,1);
 dI_0 = zeros(Nbus,1);
 dI_1 = zeros(Nbus,1);
 dI_2 = zeros(Nbus,1);
 
 % off-diagonal current injection in the series admittance matrix
for k=1:Nbranch
    if (Branchdata(k,15)==1)
         dIz0(k,1)=(Y_Seq{k,1}(1,2)*(V_1(frombus(k),1)-V_1(tobus(k),1)))+(Y_Seq{k,1}(1,3)*(V_2(frombus(k),1)-V_2(tobus(k),1)));
         dIz1(k,1)=(Y_Seq{k,1}(2,1)*(V_0(frombus(k),1)-V_0(tobus(k),1)))+(Y_Seq{k,1}(2,3)*(V_2(frombus(k),1)-V_2(tobus(k),1)));
         dIz2(k,1)=(Y_Seq{k,1}(3,1)*(V_0(frombus(k),1)-V_0(tobus(k),1)))+(Y_Seq{k,1}(3,2)*(V_1(frombus(k),1)-V_1(tobus(k),1)));
    end
end 

 % off-diagonal current injection in the shunt admittance matrix
 for m=1:Nbranch
     for n=1:Nbus
       if ((Branchdata(m,15)==1) && (frombus(m) == n))  
         dIs0(n,1)=dIs0(n,1)+(Y_Shunt{m,1}(1,2)*V_1(n,1))+(Y_Shunt{m,1}(1,3)*V_2(n,1));
         dIs1(n,1)=dIs0(n,1)+(Y_Shunt{m,1}(2,1)*V_0(n,1))+(Y_Shunt{m,1}(2,3)*V_2(n,1));
         dIs2(n,1)=dIs0(n,1)+(Y_Shunt{m,1}(3,1)*V_0(n,1))+(Y_Shunt{m,1}(3,2)*V_1(n,1)); 
       end
    end
 end
 
  for m=1:Nbranch
     for n=1:Nbus
       if ((Branchdata(m,15)==1) && (tobus(m) == n))  
         dIs0(n,1)=dIs0(n,1)+(Y_Shunt{m,1}(1,2)*V_1(n,1))+(Y_Shunt{m,1}(1,3)*V_2(n,1));
         dIs1(n,1)=dIs0(n,1)+(Y_Shunt{m,1}(2,1)*V_0(n,1))+(Y_Shunt{m,1}(2,3)*V_2(n,1));
         dIs2(n,1)=dIs0(n,1)+(Y_Shunt{m,1}(3,1)*V_0(n,1))+(Y_Shunt{m,1}(3,2)*V_1(n,1)); 
       end
    end
 end
 
%% sequence current  injections at each bus
% these represent coupling in the untransposed transmission lines

for i=1:Nbranch
    for j=1:Nbus
        if ((Branchdata(i,15)==1) && (frombus(i) == j))
            dI_0(frombus(i),1)=dI_0(frombus(i),1)+dIs0(frombus(i),1)+dIz0(i,1);
            dI_1(frombus(i),1)=dI_1(frombus(i),1)+dIs1(frombus(i),1)+dIz1(i,1);
            dI_2(frombus(i),1)=dI_2(frombus(i),1)+dIs2(frombus(i),1)+dIz2(i,1);
        end
    end
end

for i=1:Nbranch
    for j=1:Nbus
        if ((Branchdata(i,15)==1) && (tobus(i) == j))
            dI_0(tobus(i),1)=dI_0(tobus(i),1)+dIs0(tobus(i),1)-dIz0(i,1);
            dI_1(tobus(i),1)=dI_1(tobus(i),1)+dIs1(tobus(i),1)-dIz1(i,1);
            dI_2(tobus(i),1)=dI_2(tobus(i),1)+dIs2(tobus(i),1)-dIz2(i,1);
        end
    end
end

% positive sequence power injection
dI_1conj=conj(dI_1);

for q=1:Nbus
dS_1(q,1)=V_1_PV(q,1)*dI_1conj(q,1);
end

dP_1=real(dS_1);
dQ_1=imag(dS_1); 

%% total sequence specifications

P_1sp = zeros(Nbus,1);
Q_1sp = zeros(Nbus,1);
I_2sp = zeros(Nbus,1);
I_0sp = zeros(Nbus,1);

for i=1:Nbus 
     %if (Busdata(i,2)==1)
     P_1sp(i,1)=-Pload_1sp(i)+(dP_1(i)); 
     Q_1sp(i,1)=-Qload_1sp(i)+(dQ_1(i));
     I_2sp(i,1)=-Iload_2sp(i)-dI_2(i);
     I_0sp(i,1)=-Iload_0sp(i)-dI_0(i);
     %end
end

for i=1:Nbus 
     if (Busdata(i,2)==2)
     P_1sp(i,1)=P_1sp(i,1)+Pg(i,1); 
     Q_1sp(i,1)=Q_1sp(i,1)+Qg(i,1);
     end
end  

%% linear formulation for negative and zero sequence voltages
V_0=inv(Y0)*I_0sp;
V_2=inv(Y2)*I_2sp;

%% Positive sequence Newton raphson load flow
V_1=abs(V_1);

%Calculate the active and reactive powers according to initial values
        P_1=zeros(Nbus,1);            % Pcalculate
        Q_1=zeros(Nbus,1);            % Qcalculate
        
        for i=1:Nbus
            for j=1:Nbus
                if j~=i
                P_1(i)=P_1(i)+V_1(i)*V_1(j)*abs(Y1(i,j))*cos(angle(Y1(i,j))+theta_1(j)-theta_1(i));
                Q_1(i)=Q_1(i)+V_1(i)*V_1(j)*abs(Y1(i,j))*sin(angle(Y1(i,j))+theta_1(j)-theta_1(i));
                end
            end
             P_1(i)=P_1(i)+V_1(i)*V_1(i)*real(Y1(i,i));
             Q_1(i)=-Q_1(i)-V_1(i)*V_1(i)*imag(Y1(i,i));
        end
          
%Calculation of the mismatch matrix
deltaP = zeros(Nbus,1);
deltaQ = zeros(Nbus-NGenbus,1);

M = zeros(2*Nbus-NGenbus-1,1);

ctr=1;
for i = 2:Nbus
    deltaP(i) =P_1sp(i)- P_1(i) ;
    M(ctr,1) = deltaP(i);
    ctr = ctr +1;
end

for i= 2:Nbus
    %only done for PQ buses 
    if(Busdata(i,2)==1)
        deltaQ(i) = Q_1sp(i)-Q_1(i);
        M(ctr,1) = deltaQ(i);
        ctr = ctr +1 ;
    end
end


%Formation of The Jacobian Matrix
%Initialize The Jacobian Matrix to all Zeroes
J = zeros(((2*Nbus)-1-NGenbus),((2*Nbus)-1-NGenbus));

% Jacobian
% Jptheta - Derivative of Real Power Injection with Angles

    Jptheta = zeros(Nbus-1,Nbus-1);
    for i = 1:(Nbus-1)
        n=i+1;
        for k = 1:(Nbus-1)
            m=k+1;
            if m ==n
                %diagonal elements of Jptheta
                for n = 1:Nbus
                    if m~=n                        
                        Jptheta(i,k) = Jptheta(i,k) + V_1(m)*V_1(n)*abs(Y1(m,n))*sin(theta_1(n)-theta_1(m)+angle(Y1(m,n)));                 
                    end
                end
            else
                  %off diagonal elements of Jptheta                
                    Jptheta(k,i)=-V_1(m)*V_1(n)*abs(Y1(m,n))*sin(theta_1(n)-theta_1(m)+angle(Y1(m,n)));
                    Jptheta(i,k)=-V_1(m)*V_1(n)*abs(Y1(m,n))*sin(theta_1(m)-theta_1(n)+angle(Y1(m,n)));
            end
        end
    end
    
    
    % Jpv - Derivative of Real Power Injection with Voltage 
    
     Jpv = zeros(Nbus-1,NPQbus);
    for i = 1:(Nbus-1)
        m = i+1;
        for k = 1:NPQbus
            n = busPQ(k);
            if n == m
                for n = 1:Nbus
                    if m~=n
                    Jpv(i,k) = Jpv(i,k) + V_1(n)*abs(Y1(m,n))*cos(angle(Y1(m,n)+theta_1(n)-theta_1(m)));
                    end
                end
                Jpv(i,k) = Jpv(i,k) + 2*V_1(m)*(real(Y1(m,m)));
            else
                Jpv(i,k) = V_1(m)*abs(Y1(m,n))*cos(theta_1(n)-theta_1(m)+angle(Y1(m,n)));                
            end
        end
    end  
    
    
% Jqtheta - Derivative of Reactive Power Injection with Angles
    Jqtheta = zeros(NPQbus,Nbus-1);
    for i = 1:NPQbus
        m = busPQ(i);
        for k = 1:(Nbus-1)
            n = k+1;
            if n == m
                for n = 1:Nbus
                    if n~=m
                        Jqtheta(i,k) = Jqtheta(i,k) + V_1(m)* V_1(n)*abs(Y1(m,n))*cos(theta_1(n)-theta_1(m)+angle(Y1(m,n)));
                    end
                end                       
            else
                Jqtheta(i,k) = -V_1(m)* V_1(n)*abs(Y1(m,n))*cos(theta_1(n)-theta_1(m)+angle(Y1(m,n)));
            end
        end
    end
    
 % J4 - Derivative of Reactive Power Injections with Voltage
    Jqv = zeros(NPQbus,NPQbus);
    for i = 1:NPQbus
        m = busPQ(i);
        for k = 1:NPQbus
            n = busPQ(k);
            if n == m
                for n = 1:Nbus
                    if m~=n
                        Jqv(i,k) = Jqv(i,k) -V_1(n)*abs(Y1(m,n))*sin(theta_1(n)-theta_1(m)+angle(Y1(m,n)));
                    end
                end
                Jqv(i,k) = Jqv(i,k) - 2*V_1(m)*(imag(Y1(m,m)));
            else
                Jqv(i,k) = -V_1(m)*abs(Y1(m,n))*sin(theta_1(n)-theta_1(m)+angle(Y1(m,n)));
                Jqv(k,i) = -V_1(n)*abs(Y1(m,n))*sin(theta_1(m)-theta_1(n)+angle(Y1(m,n)));
            end
        end
    end
  
%The Complete Jacobian Matrix can now be formed
J = [Jptheta Jpv; Jqtheta Jqv];

D = inv(J)*M;

%Change V and theta matrix according to the new value of D
for i=2:Nbus
        theta_1(i) = theta_1(i) + D(i-1);
end

%first thirty eight values of D are taken by theta
deltaV = 39;
for i= 2:Nbus
    % updated for PQ buses
    if(Busdata(i,2)==1)
            V_1(i) = V_1(i) + D(deltaV);
            deltaV = deltaV + 1;           
    end
end

 for i = 1:Nbus
    V_1(i) = complex(V_1(i)*cos(theta_1(i)), V_1(i)*sin(theta_1(i)));
end

 for i=1:Nbus
      Vph(:,i)=T2*[(V_0(i,1)); (V_1(i,1)); (V_2(i,1))];
      
  end
  Va=(Vph(1,:)).';
  Vb=(Vph(2,:)).';
  Vc=(Vph(3,:)).';
 
Vabs_new = [abs(Va), abs(Vb),  abs(Vc)];
error = max(max(abs(Vabs_new -Vabs_old)));
count = count+1;

% disp('Error=');
% disp(error);
% disp('Number of Iterations To Converge= ');
% disp(count);

end

% disp('..........................................................')
% disp('Final Converged values, [Va Vb Vc]=');
% disp([abs(Va) abs(Vb) abs(Vc)]);
% disp('Final Converged values, [V0 V1 V2]=');
% disp([abs(V_0) abs(V_1) abs(V_2)]);
% disp('..........................................................')

%% calculating powers and voltages at distribution connected load buses
% phase powers at node 26
    S_abus26=Va(26,1)*conj(Iasp(26,1));
    P_abus26=real(S_abus26); Q_abus26=imag(S_abus26);
    S_bbus26=Vb(26,1)*conj(Ibsp(26,1));
    P_bbus26=real(S_bbus26); Q_bbus26=imag(S_bbus26);
    S_cbus26=Vc(26,1)*conj(Icsp(26,1));
    P_cbus26=real(S_cbus26); Q_cbus26=imag(S_cbus26);
    S_bus26=[S_abus26; S_bbus26; S_cbus26];
    %disp(S_bus26);
    P_bus26=[P_abus26; P_bbus26; P_cbus26];
    Q_bus26=[Q_abus26; Q_bbus26; Q_cbus26];
    P_seqbus26=[abs(Pload_0sp(26,1)); abs(Pload_1sp(26,1)); abs(Pload_2sp(26,1));];
    Q_seqbus26=[abs(Qload_0sp(26,1)); abs(Qload_1sp(26,1)); abs(Qload_2sp(26,1));];

% voltages and angles at bus 26
    VaT26=Va(26,1);
    VbT26=Vb(26,1);
    VcT26=Vc(26,1);
    Va_bus26=abs(VaT26);
    Vb_bus26=abs(VbT26);
    Vc_bus26=abs(VcT26);
    V_aangbus26=angle(Va(26,1))*180/pi;
    V_bangbus26=angle(Vb(26,1))*180/pi;
    V_cangbus26=angle(Vc(26,1))*180/pi; 

%% phase powers at node 25
    S_abus25=Va(25,1)*conj(Iasp(25,1));
    P_abus25=real(S_abus25); Q_abus25=imag(S_abus25);
    S_bbus25=Vb(25,1)*conj(Ibsp(25,1));
    P_bbus25=real(S_bbus25); Q_bbus25=imag(S_bbus25);
    S_cbus25=Vc(25,1)*conj(Icsp(25,1));
    P_cbus25=real(S_cbus25); Q_cbus25=imag(S_cbus25);
    S_bus25=[S_abus25; S_bbus25; S_cbus25];
    %disp(S_bus25);
    P_bus25=[P_abus25; P_bbus25; P_cbus25];
    Q_bus25=[Q_abus25; Q_bbus25; Q_cbus25];
    P_seqbus25=[abs(Pload_0sp(25,1)); abs(Pload_1sp(25,1)); abs(Pload_2sp(25,1));];
    Q_seqbus25=[abs(Qload_0sp(25,1)); abs(Qload_1sp(25,1)); abs(Qload_2sp(25,1));];

% voltages and angles at bus 26
    VaT25=Va(5,1);
    VbT25=Vb(5,1);
    VcT25=Vc(5,1);
    Va_bus25=abs(VaT25);
    Vb_bus25=abs(VbT25);
    Vc_bus25=abs(VcT25);
    V_aangbus25=angle(Va(25,1))*180/pi;
    V_bangbus25=angle(Vb(25,1))*180/pi;
    V_cangbus25=angle(Vc(25,1))*180/pi; 
  
%% phase powers at node 28
    S_abus28=Va(28,1)*conj(Iasp(28,1));
    P_abus28=real(S_abus28); Q_abus28=imag(S_abus28);
    S_bbus28=Vb(28,1)*conj(Ibsp(28,1));
    P_bbus28=real(S_bbus28); Q_bbus28=imag(S_bbus28);
    S_cbus28=Vc(28,1)*conj(Icsp(28,1));
    P_cbus28=real(S_cbus28); Q_cbus28=imag(S_cbus28);
    S_bus28=[S_abus28; S_bbus28; S_cbus28];
    %disp(S_bus28);
    P_bus28=[P_abus28; P_bbus28; P_cbus28];
    Q_bus28=[Q_abus28; Q_bbus28; Q_cbus28];
    P_seqbus28=[abs(Pload_0sp(28,1)); abs(Pload_1sp(28,1)); abs(Pload_2sp(28,1));];
    Q_seqbus28=[abs(Qload_0sp(28,1)); abs(Qload_1sp(28,1)); abs(Qload_2sp(28,1));];

% voltages and angles at bus 26
    VaT28=Va(28,1);
    VbT28=Vb(28,1);
    VcT28=Vc(28,1);
    Va_bus28=abs(VaT28);
    Vb_bus28=abs(VbT28);
    Vc_bus28=abs(VcT28);
    V_aangbus28=angle(Va(28,1))*180/pi;
    V_bangbus28=angle(Vb(28,1))*180/pi;
    V_cangbus28=angle(Vc(28,1))*180/pi; 
  
%% sequence voltages at bus 25
    V_0bus25=V_0(25,1);
    V_1bus25=V_1(25,1);
    V_2bus25=V_2(25,1);
    V0T25=abs(V_0bus25);
    V1T25=abs(V_1bus25);
    V2T25=abs(V_2bus25);
%% sequence voltages at bus 26
    V_0bus26=V_0(26,1);
    V_1bus26=V_1(26,1);
    V_2bus26=V_2(26,1);
    V0T26=abs(V_0bus26);
    V1T26=abs(V_1bus26);
    V2T26=abs(V_2bus26);
%% sequence voltages at bus 28
    V_0bus28=V_0(28,1);
    V_1bus28=V_1(28,1);
    V_2bus28=V_2(28,1);
    V0T28=abs(V_0bus28);
    V1T28=abs(V_1bus28);
    V2T28=abs(V_2bus28);
  
% disp ('The input power to the transmission load flow is=')
% disp (P_bus6);
% disp (Q_bus6);
% %   
%  disp('Phase Voltages from the load node at the interface=');
% disp(VaT6);
% disp(VbT6);
% disp(VcT6);
% 
% disp('Phase angles from the load node at the interface=');
% disp(V_aangbus6);
% disp(V_bangbus6);
% disp(V_cangbus6);
% 
% disp('Sequence Voltages from the load node at the interface=');
% disp(V_0bus6);
% disp(V_1bus6);
% disp(V_2bus6);

% 
% disp('The input positive Sequence powers in the load node=');
% disp(P_seqbus6);
% disp(Q_seqbus6);

% disp('angles');
% disp(del_1*180/pi);