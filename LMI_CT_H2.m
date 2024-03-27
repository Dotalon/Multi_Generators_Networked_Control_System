function [K2,rho2,feas2]=LMI_CT_H2(A,B,Ch,N,ContStruc)
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
    n(i)=size(Ch{i},1);
    Btot=[Btot,B{i}];
end
ntot=size(A,1);        %20 total number of states
mtot=sum(m);           %5 total number of inputs

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

Bw=rand([20,5]);  % 20x5 random matrix of noises, multiplies the 5x20 vector w (noise)

S=sdpvar(ntot+mtot);

Ch=[eye(ntot);zeros(mtot,ntot)]; % 25x20 matrix, first 20x20 is I (Q=I), the rest is zeros
Dh=[zeros(ntot,mtot);eye(mtot)]; % 25x5 matrix, first 20x5 is zeros, the rest is a 5x5 identity (R=I)

LMIconstr=[Y*A'+A*Y+Btot*L+L'*Btot'+Bw*Bw'<=-1e-2*eye(ntot)]+[Y>=1e-2*eye(ntot)]+[[S            Ch*Y+Dh*L;  L'*Dh'+Y*Ch'  Y]>=1e-2*eye(ntot+mtot*5)];
options=sdpsettings('solver','sedumi');
Obj=trace(S);
J=optimize(LMIconstr,Obj,options);   
feas2=J.problem;
L=double(L);
Y=double(Y);

K2=L/Y;
rho2=max(real(eig(A+Btot*K2)));