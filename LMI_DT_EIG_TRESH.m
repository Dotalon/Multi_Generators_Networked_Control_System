function [K2,rho2,feas2]=LMI_DT_EIG_CIRCLE(F,G,H,N,ContStruc,rho)
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

Gtot=[];
for i=1:N
    m(i)=size(G{i},2);
    n(i)=size(H{i},1);
    Gtot=[Gtot,G{i}];
end
ntot=size(F,1);
mtot=sum(m);

yalmip clear

if ContStruc==ones(N,N)
    % Centralized design
    P=sdpvar(ntot);
    L=sdpvar(mtot,ntot);
else
    % Decentralized/distributed design
    P=[];
    L=sdpvar(mtot,ntot);
    minc=0;
    for i=1:N
        P=blkdiag(P,sdpvar(n(i)));
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
LMIconstr=[[rho^2*P-F*P*F'-F*L'*Gtot'-Gtot*L*F'     Gtot*L;
            L'*Gtot'                               P]>=1e-2*eye(ntot*2)]+[P>=1e-2*eye(ntot)];
options=sdpsettings('solver','sedumi');
% Obj=norm(L,2)   %largest singular value of matrix L=K*Y
J=optimize(LMIconstr,[],options);   %don't really understand why it crashes if I specify that I want to minimize smth about L, as in his examples
feas2=J.problem;
L=double(L);
P=double(P);

K2=L/P;
rho2=max(real(eig(F+Gtot*K2)));
