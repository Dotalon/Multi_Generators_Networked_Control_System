function [K,rho,feas, gamma_x]=LMI_CT_DeDicont_Hinf_2(A,B,C,N,ContStruc)
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

Bw = rand([20,20]);

Ch = [eye(ntot); zeros(mtot, ntot)];
Du = [zeros(ntot,mtot);eye(mtot)];
Dw = rand([25,20]);

yalmip clear



if ContStruc==ones(N,N)
    % Centralized design
    gamma = sdpvar(1);
    Y=sdpvar(ntot);
    L=sdpvar(mtot,ntot);
else
    % Decentralized/distributed design
    Y=[];
    L=sdpvar(mtot,ntot);
    gamma = sdpvar(1);
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

LMIconstr=[[Y*A'+A*Y+Btot*L+L'*Btot',  Bw,  Y*Ch'+L'*Du';
            Bw',  -gamma*eye(ntot),  Dw';
            Ch*Y+Du*L,  Dw,  -gamma*eye(ntot+mtot)]<=1e-2*eye(3*ntot+mtot)]+[Y>=1e-2*eye(ntot)]+[gamma>=1e-3];
options=sdpsettings('solver','sdpt3');
sol = optimize(LMIconstr,gamma,options);
 
feas=sol.problem;
L=double(L);
Y=double(Y);

K=L/Y;
rho=max(real(eig(A+Btot*K)));
gamma_x = double(gamma);
