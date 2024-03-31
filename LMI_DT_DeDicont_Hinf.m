function [K,rho,feas]=LMI_DT_DeDicont(F,G,H,N,ContStruc)
% Computes, using LMIs, the distributed "state feedback" control law for the discrete-time system, with reference to the control
% information structure specified by 'ContStruc'.
%
% Inputs:
% - F: system matrix.
% - G: input matrices (i.e., G{1},..., G{N} are the input matrices of the decomposed system, one for each channel).
% - H: output matrices  (i.e., H{1},..., H{N} are the output matrices of the decomposed system, one for each channel, where [Hdec{1}',...,
% Hdec{N}']=I).
% - N: number of subsystems.
% - ContStruc: NxN matrix that specifies the information structure
% constraints (ContStruc(i,j)=1 if communication is allowed between channel
% j to channel i, ContStruc(i,j)=0 otherwise).
%
% Output:
% - K: structured control gain
% - rho: spectral radius of matrix (F+G*K) - note that [H{1}',...,
% H{N}']=I
% - feas: feasibility of the LMI problem (=0 if yes)

Gtot=[];
for i=1:N
    m(i)=size(G{i},2);
    n(i)=size(H{i},1);
    Gtot=[Gtot,G{i}];
end
ntot=size(F,1);
mtot=sum(m);

gamma_max = 100;
gamma_min = 0.1;

Gw = rand([20,20]);

Hh = [eye(ntot); zeros(mtot, ntot)];
Nu = [zeros(ntot,mtot);eye(mtot)];
Nw = rand([25,20]);

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

while(true)
    if [(gamma_max-gamma_min)/gamma_min]<1e-6
        gamma_test = (gamma_max+gamma_min/2);
        gamma_test
        break
    end
   gamma_test = gamma_max + gamma_min/2;
   LMIconstr = [[P  F*P+Gtot*L  Gw  zeros(20,25);
           P*F'+L'*Gtot'  P  zeros(20,20)  P*Hh'+L'*Nu';
           Gw'  zeros(20)  gamma_test*eye(ntot)  Nw';
           zeros(25, 20)  Hh*P+Nu*L  Nw   gamma_test*eye(ntot+mtot)]>=1e-6*eye(ntot*4+mtot)];
   ops = sdpsettings('solver', 'sedumi');
   sol = optimize(LMIconstr, [], ops);

   if (sol.problem==1)
       disp('Infeasible');
       gamma_min = gamma_test;
   else
       disp('Feasible')
       gamma_max = gamma_test;
   end
end

feas=sol.problem;
L=double(L);
P=double(P);

K=L/P;
rho=max(abs(eig(F+Gtot*K)));
