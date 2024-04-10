function [K4,rho4,feas4]=LMI_CT_REGION(A,B,C,N,ContStruc)
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
% N=[]
theta=pi/6;

% LMIconstr=[[sin(theta)*(A*Y+Y*A'+B*L+L'*B')  cos(theta)*(A*Y-Y*A'+B*L-L'*B');
%            cos(theta)*(-A*Y+Y*A'-B*L+L'*B')  sin(theta)*(A*Y+Y*A'+B*L+L'*B')]<=-1e-2*eye(ntot*2)]; %ntot?

LMIconstr=[[sin(theta)*(A*Y+Y*A'+Btot*L+L'*Btot')    cos(theta)*(A*Y-Y*A'+Btot*L-L'*Btot');
            cos(theta)*(-A*Y+Y*A'-Btot*L+L'*Btot')      sin(theta)*(A*Y+Y*A'+Btot*L+L'*Btot')]<=(-1e-2*eye(ntot*2))]+[Y>=1e-2*eye(ntot)];


options=sdpsettings('solver','sdpt3');
%Obj=norm(L,2);   %largest singular value of matrix L=K*Y
J=optimize(LMIconstr,[],options);   %don't really understand why it crashes if I specify that I want to minimize smth about L, as in his examples
feas4=J.problem;
L=double(L);
Y=double(Y);

K4=L/Y;
rho4=max(real(eig(A+Btot*K4)));
