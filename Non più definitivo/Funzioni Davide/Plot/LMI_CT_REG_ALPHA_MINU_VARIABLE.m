function [K4,rho4,feas4,kappaL]=LMI_CT_REG_ALPHA_MINU_VARIABLE(A,B,C,N,ContStruc, alpha)
% Computes, using LMIs, the distributed "state feedback" control law for the continuous-time system, with reference to the control
% information structure specified by 'ContStruc'.
%
% Inputs:
% - A: system matrix.
% - B: input matrices (i.e., B{1},..., B{N} are the input matrices of the decomposed system, one for each channel).
% - C: output matrices  (i.e., C{1},..., C{N} are the output matrices of the decomposed system, one for each channel, where [Cdec{1}',...,
% Cdec{N}']=I).
% - N: number of subsystems.
% - ContStruc: NxN matrix that specifies the information structure
% constraints (ContStruc(i,j)=1 if communication is allowed between channel
% j to channel i, ContStruc(i,j)=0 otherwise).
%
% Output:
% - K: structured control gain
% - rho: spectral abscissa of matrix (A+B*K) - note that [C{1}',...,
% C{N}']=I
% - feas: feasibility of the LMI problem (=0 if yes)

Btot=[];
for i=1:N
    m(i)=size(B{i},2);
    n(i)=size(C{i},1);
    Btot=[Btot,B{i}];
end
ntot=size(A,1);
mtot=sum(m);

yalmip clear

if ContStruc==ones(N,N)
    % Centralized design
    Y=sdpvar(ntot);
    L=sdpvar(mtot,ntot);
else
    % Decentralized/distributed design
    Y=[];
    L=sdpvar(mtot,ntot);
    minc=0;
    for i=1:N
        Y=blkdiag(Y,sdpvar(n(i)));
        ninc=0;
        for j=1:N
            if ContStruc(i,j)==0
                L(minc+1:minc+m(i),ninc+1:ninc+n(j))=zeros(m(i),n(j));
            end
            ninc=ninc+n(j);
            
        end
        minc=minc+m(i);
    end  
end

theta=pi/6;
kL=sdpvar(1,1);
kY=sdpvar(1,1);

LMIconstr=[[sin(theta)*((A*Y+Y*A')+(Btot*L+L'*Btot'))    cos(theta)*(A*Y-Y*A'+Btot*L-L'*Btot');
            cos(theta)*(-A*Y+Y*A'-Btot*L+L'*Btot')      sin(theta)*(A*Y+Y*A'+Btot*L+L'*Btot')]<=(-1e-2*eye(2*ntot))];

LMIconstr=LMIconstr+[Y>=1e-2*eye(ntot)];


LMIconstr=LMIconstr+[Y*A'+A*Y+Btot*L+L'*Btot'+2*alpha*Y<=-1e-2*eye(ntot)];

LMIconstr=LMIconstr+[[kL*eye(ntot)   L'; 
                        L             eye(mtot)]>=1e-2*eye(ntot+mtot)];

LMIconstr=LMIconstr+[[kY*eye(ntot)  eye(ntot); 
                        eye(ntot)           Y]>=1e-2*eye(ntot*2)];


Cost_funct=0.01*kL+10*kY;  
options=sdpsettings('solver','sedumi');
J=optimize(LMIconstr,Cost_funct,options);   
feas4=J.problem;
L=double(L);
Y=double(Y);

K4=L/Y;
rho4=max(real(eig(A+Btot*K4)));
kappaL=double(kL);
