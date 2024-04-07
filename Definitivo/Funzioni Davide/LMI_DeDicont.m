function [K,rho,feas]=LMI_DeDicont(A,B,C,N,ContStruc, Mode, alpha, radius)
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
% - Mode: CT for Continuos Time, DT for Discrete Time
% - Alpha: Center of the disk region. If there is no input for alpha, it is set to 0. If there is no input for radius 
%   alpha is the limit value for the spectral abscissa instead of the center of a disk
% - Radius: Radius of of the disk region.
% Output:
% - K: structured control gain
% - rho: spectral abscissa of matrix (A+B*K) - note that [C{1}',...,
% C{N}']=I
% - feas: feasibility of the LMI problem (=0 if yes)

if Mode == 'CT'
    if nargin < 7
     alpha = 0;
    end

    if nargin < 8
     radius = -1;
    end
    
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
   if radius==-1
    LMIconstr=[P*A'+A*P+Btot*L+L'*Btot'+ 2*alpha*P<=-1e-2*eye(ntot)]+[P>=1e-2*eye(ntot)];
   else
    LMIconstr=[[(radius^2-alpha^2)*P-A*P*A'-A*L'*Btot'-Btot*L*A'-alpha*(P*A'+A*P+L'*Btot'+Btot*L), Btot*L;
                 L'*Btot', P]>=1e-2*eye(ntot*2)];
   end
    options=sdpsettings('solver','sedumi');
    J=optimize(LMIconstr,[],options);
    feas=J.problem;
    L=double(L);
    P=double(P);

    K=L/P;
    rho=max(real(eig(A+Btot*K)));
elseif Mode == 'DT'
    if nargin < 7
     alpha = 0;
     radius = 1;
    end

    if nargin < 8
     radius = 1;
    end
    
    F = A;
    G = B;
    H = C;

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
      % Dentralized/distributed design
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

if radius==-1
    LMIconstr=[[P-F*P*F'-F*L'*Gtot'-Gtot*L*F' Gtot*L;
               L'*Gtot' P]>=1e-2*eye(ntot*2)];
   else
    LMIconstr=[[(radius^2-alpha^2)*P-F*P*F'-F*L'*Gtot'-Gtot*L*F'-alpha*(P*F'+F*P+L'*Gtot'+Gtot*L), Gtot*L;
                 L'*Gtot', P]>=1e-2*eye(ntot*2)];
   end
    options=sdpsettings('solver','sedumi');
    J=optimize(LMIconstr,[],options);
    feas=J.problem;
    L=double(L);
    P=double(P);

    K=L/P;
    rho=max(abs(eig(F+Gtot*K)));
else print('Mode not vailid, please use CT for continuous time or DT for Discrete time')
end
end