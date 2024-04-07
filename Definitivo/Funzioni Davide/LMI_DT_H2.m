function [K2,rho2,feas2]=LMI_DT_H2(F,G,H,N,ContStruc)
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
ntot=size(F,1);        %20
mtot=sum(m);           %5

yalmip clear

    Gw=eye(20);  

    q = [3600 0 0 0; 0 3600 0 0; 0 0 0.01^2 0; 0 0 0 0.01^2];
    Q = blkdiag(q, q, q, q, q);
    Hq = sqrt(Q);
    % Hq = eye(20)
    r = 0.5;

if ContStruc==ones(N,N)
    % Centralized design
    P=sdpvar(ntot);
    L=sdpvar(mtot,ntot);
    
    S=sdpvar(ntot+mtot);

    Hh=[Hq;zeros(mtot,ntot)]; % 25x20 matrix, first 20x20 is I (Q=I), the rest is zeros
    Nh=r*[zeros(ntot,mtot);eye(mtot)]; % 25x5 matrix, first 20x5 is zeros, the rest is a 5x5 identity (R=I)

    %Ch=[sdpvar(ntot);zeros(mtot,ntot)];
    %Dh=[zeros(ntot,mtot);sdpvar(mtot)];

    %Q=sdpvar(ntot);
    %R=sdpvar(mtot);
    %     
else
    % Decentralized/distributed design
    P=[];
    L=sdpvar(mtot,ntot);

    S=sdpvar(ntot+mtot);
    
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
    Hh=[Hq;zeros(mtot,ntot)]; % 25x20 matrix, first 20x20 is I (Q=I), the rest is zeros
    Nh=r*[zeros(ntot,mtot);eye(mtot)]; % 25x5 matrix, first 20x5 is zeros, the rest is a 5x5 identity (R=I)
end
% N=[]


LMIconstr=[[P-F*P*F'-F*L'*Gtot'-Gtot*L*F'-Gw*Gw'   Gtot*L;
            L'*Gtot'                           P]>=1e-2*eye(ntot*2)]+[P>=1e-2*eye(ntot)]+[[S Hh*P+Nh*L;  L'*Nh'+P*Hh'  P]>=1e-2*eye(ntot+mtot*5)];

options=sdpsettings('solver','sedumi');
Obj=trace(S);
J=optimize(LMIconstr,Obj,options);   
feas2=J.problem;
L=double(L);
P=double(P);

K2=L/P;
rho2=max(real(eig(F+Gtot*K2)));