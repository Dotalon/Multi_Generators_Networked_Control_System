function [K,rho,feas]=LMI_CT_DeDicont_Hinf(A,B,C,N,ContStruc)
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

gamma_max = 100;
gamma_min = 0.1;

Bw = rand([20,20]);

Ch = [eye(ntot); zeros(mtot, ntot)];
Du = [zeros(ntot,mtot);eye(mtot)];
Dw = rand([25,20]);

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

while(true)                                          % while is used to find the minimum gamma (https://www.youtube.com/watch?v=ah4FIabnTzg&t=1192s&ab_channel=ArtCell)
    if [(gamma_max-gamma_min)/gamma_min]<1e-2
        gamma_test = (gamma_max+gamma_min/2);
        gamma_test
        break
    end
   gamma_test = gamma_max + gamma_min/2;
LMIconstr=[[Y*A'+A*Y+Btot*L+L'*Btot',  Bw,  Y*Ch'+L'*Du';
            Bw',  -gamma_test*eye(ntot),  Dw';
            Ch*Y+Du*L,  Dw,  -gamma_test*eye(25)]<=1e-6*eye(65)]+[Y>=1e-6*eye(ntot)];
options=sdpsettings('solver','sdpt3');
sol = optimize(LMIconstr,[],options);
 
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
Y=double(Y);

K=L/Y;
rho=max(real(eig(A+Btot*K)));
