function [K2,rho2,feas2]=LMI_DT_H2_VARIABLE(F,G,H,N,ContStruc, v, r)
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

if ContStruc==ones(N,N)
    % Centralized design
    P=sdpvar(ntot);
    L=sdpvar(mtot,ntot);

    Gw=eye(20);  
    S=sdpvar(ntot+mtot); %25x25
    
    q = [3600 0 0 0; 0 3600 0 0; 0 0 0.01^2 0; 0 0 0 0.01^2];
    Q = blkdiag(q, q, q, q, q);
    R = (0.01^2)*eye(5);
    if v == 1
        Hq = sqrt(Q);
        Nq = sqrt(R);
    end
    if v == 2 
        Hq = eye(20)
        Nq = eye(5);
    end

    Hh=[Hq; zeros(mtot,ntot)]; % 25x20 matrix, first 20x20 is I (Q=I), the rest is zeros
    Dh=r*[zeros(ntot,mtot);Nq]; % 25x5 matrix, first 20x5 is zeros, the rest is a 5x5 identity (R=I)

else
    % Decentralized/distributed design
    P=[];
    L=sdpvar(mtot,ntot);

    Gw=eye(20); %hopefully same thing as before

    S=sdpvar(ntot+mtot); %25x25

    q = [3600 0 0 0; 0 3600 0 0; 0 0 0.01^2 0; 0 0 0 0.01^2];
    Q = blkdiag(q, q, q, q, q);
    R = (0.01^2)*eye(5);
    if v == 1
        Hq = sqrt(Q);
        Nq = sqrt(R);
    end
    if v == 2 
        Hq = eye(20)
        Nq = eye(5);
    end

    Hh=[Hq; zeros(mtot,ntot)]; % 25x20 matrix, first 20x20 is I (Q=I), the rest is zeros
    Dh=r*[zeros(ntot,mtot);Nq]; % 25x5 matrix, first 20x5 is zeros, the rest is a 5x5 identity (R=I)

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


LMIconstr=[[P-F*P*F'-F*L'*Gtot'-Gtot*L*F'-Gw*Gw'   Gtot*L;
            L'*Gtot'                           P]>=1e-2*eye(ntot*2)]+[P>=1e-2*eye(ntot)]+[[S Hh*P+Dh*L;  L'*Dh'+P*Hh'  P]>=1e-2*eye(ntot+mtot*5)];
options=sdpsettings('solver','sdpt3');
Obj=trace(S);
J=optimize(LMIconstr,Obj,options);   
feas2=J.problem;
L=double(L);
P=double(P);

K2=L/P;
rho2=max(real(eig(F+Gtot*K2)));