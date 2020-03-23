function [P_bus6,Q_bus6,S_bus6,P_seqbus6,Q_seqbus6,V_abus6,V_bbus6,V_cbus6,V_0bus6,V_1bus6,V_2bus6,V_aangbus6,V_bangbus6,V_cangbus6]=Trans_loadflow(Pl_update,Ql_update)
 BUS_DATA=[
 % Busno  Type  Vsp   theta   Pg     Qg       Pl       Ql       Qmin     Qmax
   1      3     1.04   0.0     0.0    0.0       0.0      0.0       0.0      0.0;   
   2      2     1.025  0.0   163.0    5.0       0.0      0.0    -300.0    300.0;   
   3      2     1.025  0.0    85.0  -11.0       0.0      0.0    -300.0    300.0;   
   4      1     1      0.0     0.0    0.0       0.0      0.0       0.0      0.0;  
   5      1     1      0.0     0.0    0.0     125.0     50.0       0.0      0.0;  
   6      1     1      0.0     0.0    0.0     90.00     30.0       0.0      0.0;  
   7      1     1      0.0     0.0    0.0       0.0      0.0       0.0      0.0;  
   8      1     1      0.0     0.0    0.0     100.0     35.0       0.0      0.0;  
   9      1     1      0.0     0.0    0.0       0.0      0.0       0.0      0.0; ];

%% BRANCH DATA FOLLOWS   
BRANCH_DATA=[
% From  To     Tap   R0            X0         b0/2           R1           X1        b1/2       R2           X2        b2/2 
    1	4	    1    0            0.0576      0              0            0.0576      0      0            0.0576       0;
    2	7     	1    0            0.0625      0              0            0.0625      0      0            0.0625       0; 
    3	9	  	1    0            0.0586      0              0            0.0586      0      0            0.0586       0;  
	4	6		1    0.0268	      0.1869    0.0718           0.0268	      0.1869    0.13675  0.0268	      0.1869     0.13675;
	6	9		1    0.0268	      0.1869    0.0718           0.0268	      0.1869    0.13675  0.0268	      0.1869     0.13675;                                                          
	8	9		1    0.0268	      0.1869    0.0718           0.0268	      0.1869    0.13675  0.0268	      0.1869     0.13675;
	7	8		1    0.0268	      0.1869    0.0718           0.0268	      0.1869    0.13675  0.0268	      0.1869     0.13675;
	7	5	    1    0.0268	      0.1869    0.0718           0.0268	      0.1869    0.13675  0.0268	      0.1869     0.13675;
	5	4		1    0.0268	      0.1869    0.0718           0.0268	      0.1869    0.13675  0.0268	      0.1869     0.13675;];

BMva=100;   
BranchFrom=BRANCH_DATA(:,1);    % From bus number...
BranchTo=BRANCH_DATA(:,2);      % To bus number... 
nl = length(BranchFrom);        % No. of branches...
bus = BUS_DATA(:,1);            % Bus Number..
type = BUS_DATA(:,2);           % Type of Bus 3-Slack, 2-PV, 1-PQ.
pv = find(type == 2 );          % PV Buses..
pq = find(type == 1);           % PQ Buses..
Pg = BUS_DATA(:,5)/(BMva);    % PGi..
Qg = BUS_DATA(:,6)/(BMva);    % QGi..
Pl= BUS_DATA(:,7)/(BMva*3);         % PLi..
Ql = BUS_DATA(:,8)/(BMva*3);        % QLi..
nbus=length(bus);               % No. of buses
npv = 3;                        % No. of PV buses..
npq = length(pq);               % No. of PQ buses..
V=BUS_DATA(:,3);                % Voltage magnitude
del = BUS_DATA(:,4)*pi/180;     % Voltage Angle..
a = BRANCH_DATA(:,3);           % transformer taps
R0=BRANCH_DATA(:,4);
X0=BRANCH_DATA(:,5);
R1=BRANCH_DATA(:,7);
X1=BRANCH_DATA(:,8);
R2=BRANCH_DATA(:,10);
X2=BRANCH_DATA(:,11);
b0=sqrt(-1)*BRANCH_DATA(:,6);
b1=sqrt(-1)*BRANCH_DATA(:,9);
b2=sqrt(-1)*BRANCH_DATA(:,12);
z0=R0+X0*sqrt(-1);
z1=R1+X1*sqrt(-1);
z2=R2+X2*sqrt(-1);
y0= 1./z0; 
y1= 1./z1;  
y2= 1./z2;  
Y0= zeros(nbus,nbus); 
Y1= zeros(nbus,nbus); 
Y2= zeros(nbus,nbus); 

%% sequence admittance matrices (YBus)
% The y bus matrix is formed considering only the diagonal elements of the
% line series and shunt impedance data ( the mutal coupling is included as
% current injections)

% Y0
 % Formation of the Off Diagonal Elements...
 %for k = 1:nl
     %Y0(BranchFrom(k),BranchTo(k)) = Y0(BranchFrom(k),BranchTo(k)) - y0(k)/a(k);
     %Y0(BranchTo(k),BranchFrom(k)) = Y0(BranchFrom(k),BranchTo(k));
 %end
 
 % Formation of Diagonal Elements....
 for m = 1:nbus
     for n = 1:nl
         if BranchFrom(n) == m
             Y0(m,m) = Y0(m,m) + y0(n)/(a(n)^2) + b0(n);
         end
     end
 end
 for m = 1:nbus
     for n = 1:nl
         if BranchTo(n) == m
             Y0(m,m) = Y0(m,m) + y0(n)/(a(n)^2) + b0(n);
         end
     end
 end
% Y1 
 % Formation of the Off Diagonal Elements...
 for k = 1:nl
     Y1(BranchFrom(k),BranchTo(k)) = Y1(BranchFrom(k),BranchTo(k)) - y1(k)/a(k);
     Y1(BranchTo(k),BranchFrom(k)) = Y1(BranchFrom(k),BranchTo(k));
 end
 
 % Formation of Diagonal Elements....
 for m = 1:nbus
     for n = 1:nl
         if BranchFrom(n) == m
             Y1(m,m) = Y1(m,m) + y1(n)/(a(n)^2) + b1(n);
         end
     end
 end
 for m = 1:nbus
     for n = 1:nl
         if BranchTo(n) == m
             Y1(m,m) = Y1(m,m) + y1(n)/(a(n)^2) + b1(n);
         end
     end
 end
% Y2 
  % Formation of the Off Diagonal Elements...
 %for k = 1:nl
     %Y2(BranchFrom(k),BranchTo(k)) = Y2(BranchFrom(k),BranchTo(k)) - y2(k)/a(k);
     %Y2(BranchTo(k),BranchFrom(k)) = Y2(BranchFrom(k),BranchTo(k));
 %end
 
 % Formation of Diagonal Elements....
 for m = 1:nbus
     for n = 1:nl
         if BranchFrom(n) == m
             Y2(m,m) = Y2(m,m) + y2(n)/(a(n)^2) + b2(n);
         end
     end
 end
 for m = 1:nbus
     for n = 1:nl
         if BranchTo(n) == m
             Y2(m,m) = Y2(m,m) + y2(n)/(a(n)^2) + b2(n);
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
Va=V*(v^3);
Vb=V*(v^2);
Vc=V*(v);

 %% sequence voltages on all buses
for k=1:nbus
Vseq(:,k)=T1*[Va(k,1); Vb(k,1);Vc(k,1);];
end
V_0=(Vseq(1,:)).';
V_1=(Vseq(2,:)).';
V_2=(Vseq(3,:)).';
for k=1:nbus
Vseq_PV(:,k)=T3*[Va(k,1); Vb(k,1);Vc(k,1);];
end
V_0_PV=(Vseq_PV(1,:)).';
V_1_PV=(Vseq_PV(2,:)).';
V_2_PV=(Vseq_PV(3,:)).';

%% line series admittance data
% Z012=[0.0268+0.1869i   0.0018-0.0011i   -0.0018-0.0011i;
%       -0.0018-0.0011i   0.0052+0.0574i   -0.0037+0.0021i;
%       0.0018-0.0011i   0.0037+0.0021i   0.0052+0.0574i;];
%   Yz012=1./Z012;
Zabc = [0.0124 + 0.1020i   0.0072 + 0.0445i   0.0072 + 0.0445i 
   0.0072 + 0.0445i   0.0124 + 0.1020i    0.0072 + 0.0445i 
   0.0072 + 0.0445i     0.0072 + 0.0445i   0.0124 + 0.1020i];
Yabc = inv(Zabc); 
Yz012 = T1*Yabc*T2;

%%  line shunt admittance data
Ysabc=[ 0.2273i -0.0518i -0.0263i;
       -0.0518i  0.2361i -0.0518i;
       -0.0263i -0.0518i  0.2273i;];
SSA= T1*Ysabc*T2;
Ys012=(SSA/2);
error = 1;
count = 0;
del_1=zeros(nbus,1);

 %% off-diagonal current injection in series & shunt admittance matrix in untransposed lines
 dIz0=zeros(nl,nl);
 dIz1=zeros(nl,nl);
 dIz2=zeros(nl,nl);
 dIs0=zeros(1,nl);
 dIs1=zeros(1,nl);
 dIs2=zeros(1,nl);
for k=4:nl
         dIz0(BranchFrom(k),BranchTo(k))=(Yz012(1,2)*(V_1(BranchFrom(k),1)-V_1(BranchTo(k),1)))+(Yz012(1,3)*(V_2(BranchFrom(k),1)-V_2(BranchTo(k),1)));
         dIz1(BranchFrom(k),BranchTo(k))=(Yz012(2,1)*(V_0(BranchFrom(k),1)-V_0(BranchTo(k),1)))+(Yz012(2,3)*(V_2(BranchFrom(k),1)-V_2(BranchTo(k),1)));
         dIz2(BranchFrom(k),BranchTo(k))=(Yz012(3,1)*(V_0(BranchFrom(k),1)-V_0(BranchTo(k),1)))+(Yz012(3,2)*(V_1(BranchFrom(k),1)-V_1(BranchTo(k),1)));
end

 for m=4:nl
         dIs0(m)=(Ys012(1,2)*V_1(m,1))+(Ys012(1,3)*V_2(m,1));
         dIs1(m)=(Ys012(2,1)*V_0(m,1))+(Ys012(2,3)*V_2(m,1));
         dIs2(m)=(Ys012(3,1)*V_0(m,1))+(Ys012(3,2)*V_1(m,1));  
end
     dIs0=(dIs0).';
         dIs1=(dIs1).';
         dIs2=(dIs2).';
%% sequence current  injections at each bus
% these represent coupling in the untransposed transmission lines
dI_0 = zeros(nl,1);
dI_1 = zeros(nl,1);
dI_2 = zeros(nl,1);
for i=1:nbus
            dI_0(BranchFrom(i))=dI_0(BranchFrom(i))+dIs0(BranchFrom(i),1)+dIz0(BranchFrom(i),BranchTo(i));
            dI_0(BranchTo(i))=dI_0(BranchTo(i))+dIs0(BranchTo(i),1)-dIz0(BranchFrom(i),BranchTo(i));
            dI_1(BranchFrom(i))=dI_1(BranchFrom(i))+dIs1(BranchFrom(i),1)+dIz1(BranchFrom(i),BranchTo(i));
            dI_1(BranchTo(i))=dI_1(BranchTo(i))+dIs1(BranchTo(i),1)-dIz2(BranchFrom(i),BranchTo(i));
            dI_2(BranchFrom(i))=dI_2(BranchFrom(i))+dIs2(BranchFrom(i),1)+dIz2(BranchFrom(i),BranchTo(i));
            dI_2(BranchTo(i))=dI_2(BranchTo(i))+dIs2(BranchTo(i),1)-dIz2(BranchFrom(i),BranchTo(i));
end
dI_2=[0+0i; 0+0i; 0+0i; -0.0172-0.01i; -0.0172-0.01i; -0.0172-0.01i; -0.0172-0.01i; -0.0172-0.01i; -0.0172-0.01i]; 
dI_0=[0+0i; 0+0i; 0+0i; -0.0048+0.0028i; -0.0048+0.0028i; -0.0048+0.0028i; -0.0048+0.0028i; -0.0048+0.0028i; -0.0048+0.0028i];

% positive sequence power injection
dI_1conj=conj(dI_1);
for q=1:nbus
dS_1(q,1)=3*V_1_PV(q,1)*dI_1conj(q,1);
end
dP_1=real(dS_1);
dQ_1=imag(dS_1);

% specified powers at load busbars
Pl_a=Pl; Ql_a=Ql;
Pl_b=Pl; Ql_b=Ql;
Pl_c=Pl; Ql_c=Ql;

%% update powers for next cosim iteration in all phases
Pl_a(6,1)=Pl_update(1,1);
Ql_a(6,1)=Ql_update(1,1);
Pl_b(6,1)=Pl_update(2,1);
Ql_b(6,1)=Ql_update(2,1);
Pl_c(6,1)=Pl_update(3,1);
Ql_c(6,1)=Ql_update(3,1);

Ssp_a=Pl_a+Ql_a*(sqrt(-1));  
Ssp_b=Pl_b+Ql_b*(sqrt(-1));  
Ssp_c=Pl_c+Ql_c*(sqrt(-1));  

% phase components injected currents
while error>0.0001
for t=1:nbus
Iasp(t,1)=conj(Ssp_a(t,1)/Va(t,1));
Ibsp(t,1)=conj(Ssp_b(t,1)/Vb(t,1));
Icsp(t,1)=conj(Ssp_c(t,1)/Vc(t,1));
end

% sequence current injections
for h=1:nbus
    Iseq(:,h)=T1*[Iasp(h,1); Ibsp(h,1); Icsp(h,1);];
end
for h=1:nbus
    Iseq_PV(:,h)=T3*[Iasp(h,1); Ibsp(h,1); Icsp(h,1);];
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
for f=1:nbus
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
%% final sequence specifications at any busbar i,

P_1sp = zeros(nl,1);
Q_1sp = zeros(nl,1);
I_2sp = zeros(nl,1);
I_0sp = zeros(nl,1);
% total sequence specifications
for i=1:nbus 
     if(BUS_DATA(i,2)==1)
     P_1sp(i,1)=-(Pload_1sp(i)+(dP_1(i))); 
     Q_1sp(i,1)=-(Qload_1sp(i)+(dQ_1(i)));
     I_2sp(i,1)=(-Iload_2sp(i)-(dI_2(i)));
     I_0sp(i,1)=(-Iload_0sp(i)-(dI_0(i)));
     end
end
for i=1:nbus 
     if(BUS_DATA(i,2)==2)
     P_1sp(i,1)=Pg(i,1); 
     Q_1sp(i,1)=Qg(i,1);
     
     end
end  
% linear formulation for negative and zero sequence voltages
V_0=inv(Y0)*I_0sp;
V_2=inv(Y2)*I_2sp;
%% NR positive sequence

V_1=abs(V_1);
%while error >0.0001
%% calculated powers for NR
        P_1=zeros(nbus,1);            % Pcalculate
        Q_1=zeros(nbus,1);            % Qcalculate
        
        for i=1:nbus
            for j=1:nbus
                if j~=i
                P_1(i)=P_1(i)+V_1(i)*V_1(j)*abs(Y1(i,j))*cos(angle(Y1(i,j))+del_1(j)-del_1(i));
                Q_1(i)=Q_1(i)+V_1(i)*V_1(j)*abs(Y1(i,j))*sin(angle(Y1(i,j))+del_1(j)-del_1(i));
                end
            end
             P_1(i)=P_1(i)+V_1(i)*V_1(i)*real(Y1(i,i));
             Q_1(i)=-Q_1(i)-V_1(i)*V_1(i)*imag(Y1(i,i));
        end
    
%Calculation of the mismatch matrix
deltaP = zeros(nbus,1);
deltaQ = zeros(nbus-npv,1);

%Mismatch Matrix
M = zeros(2*nbus-npv-1,1);

ctr=1;
for i = 2:nbus
    deltaP(i,1) =P_1sp(i,1)- P_1(i,1);
    M(ctr,1) = deltaP(i,1);
    ctr = ctr +1;
end

for i= 2:nbus
    if(BUS_DATA(i,2)==1)
        deltaQ(i,1) = Q_1sp(i,1)-Q_1(i,1);
        M(ctr,1) = deltaQ(i,1);
        ctr = ctr +1 ;
    end
end

%Initialize The Jacobian Matrix to all Zeroes
J = zeros(((2*nbus)-1-npv),((2*nbus)-1-npv));

% Jptheta - Derivative of Real Power Injections with Angles..

    Jptheta = zeros(nbus-1,nbus-1);
    for i = 1:(nbus-1)
        n=i+1;
        for k = 1:(nbus-1)
            m=k+1;
            if m ==n
                %diagonal elements of Jptheta
                for n = 1:nbus
                    if m~=n                        
                        Jptheta(i,k) = Jptheta(i,k) + V_1(m)*V_1(n)*abs(Y1(m,n))*sin(del_1(n)-del_1(m)+angle(Y1(m,n)));                 
                    end
                end
            else
                  %off diagonal elements of Jptheta                
                    Jptheta(k,i)=-V_1(m)*V_1(n)*abs(Y1(m,n))*sin(del_1(n)-del_1(m)+angle(Y1(m,n)));
                    Jptheta(i,k)=-V_1(m)*V_1(n)*abs(Y1(m,n))*sin(del_1(m)-del_1(n)+angle(Y1(m,n)));
            end
        end
    end
    
    
    % Jpv - Derivative of Real Power Injections with V..
    
     Jpv = zeros(nbus-1,npq);
    for i = 1:(nbus-1)
        m = i+1;
        for k = 1:npq
            n = pq(k);
            if n == m
                for n = 1:nbus
                    if m~=n
                    Jpv(i,k) = Jpv(i,k) + V_1(n)*abs(Y1(m,n))*cos(angle(Y1(m,n))+del_1(n)-del_1(m));
                    end
                end
                Jpv(i,k) = Jpv(i,k) + 2*V_1(m)*(real(Y1(m,m)));
            else
                Jpv(i,k) = V_1(m)*abs(Y1(m,n))*cos(del_1(n)-del_1(m)+angle(Y1(m,n)));                
            end
        end
    end  
    
    
% Jqtheta - Derivative of Reactive Power Injections with Angles..
    Jqtheta = zeros(npq,nbus-1);
    for i = 1:npq
        m = pq(i);
        for k = 1:(nbus-1)
            n = k+1;
            if n == m
                for n = 1:nbus
                    if n~=m
                        Jqtheta(i,k) = Jqtheta(i,k) + V_1(m)* V_1(n)*abs(Y1(m,n))*cos(del_1(n)-del_1(m)+angle(Y1(m,n)));
                    end
                end                       
            else
                Jqtheta(i,k) = -V_1(m)* V_1(n)*abs(Y1(m,n))*cos(del_1(n)-del_1(m)+angle(Y1(m,n)));
            end
        end
    end
    
 % J4 - Derivative of Reactive Power Injections with V..
    Jqv = zeros(npq,npq);
    for i = 1:npq
        m = pq(i);
        for k = 1:npq
            n = pq(k);
            if n == m
                for n = 1:nbus
                    if m~=n
                        Jqv(i,k) = Jqv(i,k) -V_1(n)*abs(Y1(m,n))*sin(del_1(n)-del_1(m)+angle(Y1(m,n)));
                    end
                end
                Jqv(i,k) = Jqv(i,k) - 2*V_1(m)*(imag(Y1(m,m)));
            else
                Jqv(i,k) = -V_1(m)*abs(Y1(m,n))*sin(del_1(n)-del_1(m)+angle(Y1(m,n)));
                Jqv(k,i) = -V_1(n)*abs(Y1(m,n))*sin(del_1(m)-del_1(n)+angle(Y1(m,n)));
            end
        end
    end
  
%The Complete Jacobian Matrix can now be formed
J = [Jptheta Jpv; Jqtheta Jqv];
D= inv(J)*M;

for i=2:nbus
    del_1(i) = del_1(i) + D(i-1);
end

dV = 9;
for i= 2:nbus
    %Only should update the non generator buses
    if(BUS_DATA(i,2)==1)
            V_1(i) = V_1(i) + D(dV);
            dV = dV + 1;    
    end
end
 
 count = count + 1;
 error = max(abs(M));
 for i = 1:nbus
    V_1(i) = complex(V_1(i)*cos(del_1(i)), V_1(i)*sin(del_1(i)));
end

 for i=1:nbus
      Vph(:,i)=T2*[(V_0(i,1)); (V_1(i,1)); ((V_2(i,1)))];
      
  end
  Va=(Vph(1,:)).';
  Vb=(Vph(2,:)).';
  Vc=(Vph(3,:)).'; 
  
%   disp(abs(Va));
%   disp(abs(Vb));
%   disp(abs(Vc));
%     disp(abs(Va));
%   disp(abs(Vb));
%   disp(abs(Vc));
%  V_aangbus=angle(Va)*180/pi;
% V_bangbus=angle(Vb)*180/pi;
% V_cangbus=angle(Vc)*180/pi; 
%   disp(V_aangbus);
%   disp(V_bangbus);
%   disp(V_cangbus);
%   disp(count);
% disp(error);
end
%% phase powers at node 6
S_abus6=Va(6,1)*conj(Iasp(6,1));
P_abus6=real(S_abus6); Q_abus6=imag(S_abus6);
S_bbus6=Vb(6,1)*conj(Ibsp(6,1));
P_bbus6=real(S_bbus6); Q_bbus6=imag(S_bbus6);
S_cbus6=Vc(6,1)*conj(Icsp(6,1));
P_cbus6=real(S_cbus6); Q_cbus6=imag(S_cbus6);
S_bus6=[S_abus6; S_bbus6; S_cbus6];
%disp(S_bus6);
P_bus6=[P_abus6; P_bbus6; P_cbus6];
Q_bus6=[Q_abus6; Q_bbus6; Q_cbus6];
P_seqbus6=[abs(Pload_0sp(6,1)); abs(Pload_1sp(6,1)); abs(Pload_2sp(6,1));];
Q_seqbus6=[abs(Qload_0sp(6,1)); abs(Qload_1sp(6,1)); abs(Qload_2sp(6,1));];

% voltages and angles at bus 6
  V_abus6=abs(Va(6,1));
  V_bbus6=abs(Vb(6,1));
  V_cbus6=abs(Vc(6,1));

  V_aangbus6=angle(Va(6,1))*180/pi;
  V_bangbus6=angle(Vb(6,1))*180/pi;
  V_cangbus6=angle(Vc(6,1))*180/pi; 
  
  %% sequence voltages at bus 6
  V_0bus6=abs(V_0(6,1));
  V_1bus6=abs(V_1(6,1));
  V_2bus6=abs(V_2(6,1));
  
% disp ('The input power to the transmission load flow is=')
% disp (P_bus6);
% disp (Q_bus6);
%   
% disp('Phase Voltages from the load node at the interface=');
% disp(V_abus6);
% disp(V_bbus6);
% disp(V_cbus6);
% 
% disp('Phase angles from the load node at the interface=');
% disp(V_aangbus6);
% disp(V_bangbus6);
% disp(V_cangbus6);
% 
% % disp('Sequence Voltages from the load node at the interface=');
% % disp(V_0bus6);
% % disp(V_1bus6);
% % disp(V_2bus6);
% 
% % 
% % disp('The input positive Sequence powers in the load node=');
% % disp(P_seqbus6);
% % disp(Q_seqbus6);
% 
% % disp('angles');
% % disp(del_1*180/pi);