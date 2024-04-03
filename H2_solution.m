function H2_solution(A,Bdec,Cdec,F,Gdec,Hdec,N,rounding_n,Bw)
%% Solution attempt with SIMPLY STABILIZING control law
% n= rounding decimal
[ContStruc,K,rho,feas]=Optimize_ContStruc(A,Bdec,Cdec,N,'LMI_CT_H2',Bw);
CFM_H2=di_fixed_modes(A,Bdec,Cdec,N,ContStruc, rounding_n)
DFM_H2=di_fixed_modes(F,Gdec,Hdec,N,ContStruc, rounding_n)

K
rho
feas
ContStruc
