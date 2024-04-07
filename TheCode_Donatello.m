% Power Network system
% Marcello Farina, 19/12/2018
% SISOgna_il_NW, --/04/2024

clc
clear all

%%
H1=12;      H2=10;      H3=8;       H4=8;       H5=10;
R1=0.05;    R2=0.0625;  R3=0.08;    R4=0.08;    R5=0.05;
D1=0.7;     D2=0.9;     D3=0.9;     D4=0.7;     D5=0.86;
Tt1=0.65;   Tt2=0.4;    Tt3=0.3;    Tt4=0.6;    Tt5=0.8;
Tg1=0.1;    Tg2=0.1;    Tg3=0.1;    Tg4=0.1;    Tg5=0.15;

P12=0.001*4;  P21=P12;
P23=0.001*2;  P32=P23;
P34=0.001*2;  P43=P34;
P45=0.001*3;  P54=P45;
P25=0.001*3;  P52=P25;


A11=[   0                               1               0               0;
        -(1/(2*H1))*P12                 -D1/(2*H1)      1/(2*H1)        0;
        0                               0               -1/Tt1          1/Tt1;
        0                               -1/(R1*Tg1)     0               -1/Tg1];
    
A22=[   0                               1               0               0;
        -(1/(2*H2))*(P21+P25+P23)       -D2/(2*H2)      1/(2*H2)        0;
        0                               0               -1/Tt2          1/Tt2;
        0                               -1/(R2*Tg2)     0               -1/Tg2];
    
A33=[   0                               1               0               0;
        -(1/(2*H3))*(P32+P34)           -D3/(2*H3)      1/(2*H3)        0;
        0                               0               -1/Tt3          1/Tt3;
        0                               -1/(R3*Tg3)     0               -1/Tg3];
    
A44=[   0                               1               0               0;
        -(1/(2*H4))*(P43+P45)           -D4/(2*H4)      1/(2*H4)        0;
        0                               0               -1/Tt4          1/Tt4;
        0                               -1/(R4*Tg4)     0               -1/Tg4];
    
A55=[   0                               1               0               0;
        -(1/(2*H5))*(P52+P54)           -D5/(2*H5)      1/(2*H5)        0;
        0                               0               -1/Tt5          1/Tt5;
        0                               -1/(R5*Tg5)     0               -1/Tg5];

B1=[0 0 0 1/Tg1]';
B2=[0 0 0 1/Tg2]';
B3=[0 0 0 1/Tg3]';
B4=[0 0 0 1/Tg4]';
B5=[0 0 0 1/Tg5]';

A12=[   0   0   0   0;  P12/(2*H1)  0   0   0;  0   0   0   0;  0   0   0   0];
A13=zeros(4,4);
A14=zeros(4,4);
A15=zeros(4,4);
A21=[   0   0   0   0;  P21/(2*H2)  0   0   0;  0   0   0   0;  0   0   0   0];
A23=[   0   0   0   0;  P23/(2*H2)  0   0   0;  0   0   0   0;  0   0   0   0];
A24=zeros(4,4);
A25=[   0   0   0   0;  P25/(2*H2)  0   0   0;  0   0   0   0;  0   0   0   0];
A31=zeros(4,4);
A32=[   0   0   0   0;  P32/(2*H3)  0   0   0;  0   0   0   0;  0   0   0   0];
A34=[   0   0   0   0;  P34/(2*H3)  0   0   0;  0   0   0   0;  0   0   0   0];
A35=zeros(4,4);
A41=zeros(4,4);
A42=zeros(4,4);
A43=[   0   0   0   0;  P43/(2*H4)  0   0   0;  0   0   0   0;  0   0   0   0];
A45=[   0   0   0   0;  P45/(2*H4)  0   0   0;  0   0   0   0;  0   0   0   0];
A51=zeros(4,4);
A52=[   0   0   0   0;  P52/(2*H5)  0   0   0;  0   0   0   0;  0   0   0   0];
A53=zeros(4,4);
A54=[   0   0   0   0;  P54/(2*H5)  0   0   0;  0   0   0   0;  0   0   0   0];

C1=eye(4); 
C2=eye(4);
C3=eye(4);
C4=eye(4);
C5=eye(4);
%Ctot is an identity! all states are measurable!

% centralized model matrices:
    
A=[ A11 A12 A13 A14 A15;
    A21 A22 A23 A24 A25;
    A31 A32 A33 A34 A35;
    A41 A42 A43 A44 A45;
    A51 A52 A53 A54 A55];

B=blkdiag(B1,B2,B3,B4,B5);
C=blkdiag(C1,C2,C3,C4,C5);

%% Decoupled and discretized model
Bdec=[];
Cdec=[];
Hdec=[];
Gdec=[];
%in teoria la decompongo guardando la matrice A: osservo le interazioni fra
%i vari sistemi in maniera tale da decidere come decomporre X vettore degli
%stati n=4. 
h=1;
N=5;
n_states=4;

Bw=rand([20,5]); %NOISEs TO BE USED IN H2 and Hinf
Gw=rand([20,5]);

[F,G,H,L,h]=ssdata(c2d(ss(A,B,C,[]),h));


for i=1:N
    Bdec{i}=B(:,i);
    Cdec{i}=C(n_states*(i-1)+1:n_states*i,:);
    Gdec{i}=G(:,i);
    Hdec{i}=H(n_states*(i-1)+1:n_states*i,:);
end

rounding_n=3;


%% spectral abscissa and radius
eigenvaluesCT = eig(A)
spectral_abscissa = round(max(real(eig(A))),10)

eigenvaluesDT = eig(F)
spectral_radius = round(max(abs(eig(F))),10)

%the open loop system is simply stable as we have only one eig on the border (zero for CT , 1 for DT),
%in both cases

%% state-feedback centralized control structures
rounding_n=3;
ContStruc=ones(N,N);
%[ 1 1 1 1 1
%  1 1 1 1 1
%  1 1 1 1 1
%  1 1 1 1 1];  %centralized structure

CFM=di_fixed_modes(A,Bdec,Cdec,N,ContStruc, rounding_n)
%no fixed modes in Continuous
DFM=di_fixed_modes(F,Gdec,Hdec,N,ContStruc, rounding_n)
%no fixed modes in Discrete

%Continuous Time stabilizing gains
[K_cent,rho_cent,feas_cent]=LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStruc);
rho_cent
feas_cent
%Discrete Time stabilizing gains
[K_dt_cent,rho_dt_cent,feas_dt_cent]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStruc);
rho_dt_cent
feas_dt_cent

eigenvaluesCT = eig(A+B*K_cent)
spectral_abscissa = round(max(real(eig(A))),10)

eigenvaluesDT = eig(F+G*K_dt_cent)
spectral_radius = round(max(abs(eig(F))),10)

%CT and DT H2 gains
[K_cent2,rho_cent2,feas_cent2]=LMI_CT_H2(A,Bdec,Cdec,N,ContStruc,Bw);
rho_cent2
feas_cent2
[K_dt_cent2,rho_dt_cent2,feas_dt_cent2]=LMI_DT_H2(F,Gdec,Hdec,N,ContStruc,Gw);
rho_dt_cent2
feas_dt_cent2

eigenvaluesCT = eig(A+B*K_cent2)
spectral_abscissa = round(max(real(eig(A))),10)

eigenvaluesDT = eig(F+G*K_dt_cent2)
spectral_radius = round(max(abs(eig(F))),10)

%% state-feedback Decentralized control structure
rounding_n=3;
ContStruc_dec=diag(ones(N,1));
%[ 1 0 0 0 0
%  0 1 0 0 0
%  0 0 1 0 0
%  0 0 0 1 0
%  0 0 0 0 1];  %decentralized structure

CFM_dec=di_fixed_modes(A,Bdec,Cdec,N,ContStruc_dec, rounding_n)
%no fixed modes in Continuous
DFM_dec=di_fixed_modes(F,Gdec,Hdec,N,ContStruc_dec, rounding_n)
%no fixed modes in Discrete

%Continuous Time gains
[K_dec,rho_dec,feas_dec]=LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStruc_dec);
K_dec
rho_dec
feas_dec

%Discrete Time gains
[K_dt_dec,rho_dt_dec,feas_dt_dec]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStruc_dec);
K_dt_dec
rho_dt_dec
feas_dt_dec

%% state-feedback distributed control structures
% In case one wants to create a sort of a hierarchy with a master subsystem
% e.g. the one that produces most, or whose production is most important.
% (in this case it's subsystem 1)

%Bidirectional star topology
ContStruc_BiStar= [ 1 1 1 1 1
                    1 1 0 0 0
                    1 0 1 0 0
                    1 0 0 1 0
                    1 0 0 0 1];

CFM_BiStar=di_fixed_modes(A,Bdec,Cdec,N,ContStruc_BiStar, rounding_n)
%no fixed modes in Continuous
DFM_BiStar=di_fixed_modes(F,Gdec,Hdec,N,ContStruc_BiStar, rounding_n)
%no fixed modes in Discrete

%Continuous Time gains
[K_BiStar,rho_BiStar,feas_BiStar]=LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStruc_BiStar);
K_BiStar
rho_BiStar
feas_BiStar

%Discrete Time gains
[K_dt_BiStar,rho_dt_BiStar,feas_dt_BiStar]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStruc_BiStar);
K_dt_BiStar
rho_dt_BiStar
feas_dt_BiStar

%% Reinforced Bi-Directional Star
% Like Bi Star but we add extra connections to try improve the overall
% stability.
%Reinforced Bidirectional Star topology
ContStruc_REBiStar= [ 1 1 1 1 1
                      1 1 0 0 1
                      1 1 1 0 0
                      1 0 1 1 0
                      1 0 0 1 1];

CFM_REBiStar=di_fixed_modes(A,Bdec,Cdec,N,ContStruc_REBiStar, rounding_n)
%no fixed modes in Continuous
DFM_REBiStar=di_fixed_modes(F,Gdec,Hdec,N,ContStruc_REBiStar, rounding_n)
%no fixed modes in Discrete

%Continuous Time stabilizing gains
[K_RE,rho_RE,feas_RE]=LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStruc_REBiStar);
rho_RE
feas_RE
%Discrete Time stabilizing gains
[K_dt_RE,rho_dt_RE,feas_dt_RE]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStruc_REBiStar);
rho_dt_RE
feas_dt_RE

eigenvaluesCT = eig(A+B*K_RE)
spectral_abscissa = round(max(real(eig(A))),10)

eigenvaluesDT = eig(F+G*K_dt_RE)
spectral_radius = round(max(abs(eig(F))),10)

%CT and DT H2 gains
[K_RE2,rho_RE2,feas_RE2]=LMI_CT_H2(A,Bdec,Cdec,N,ContStruc_REBiStar,Bw);
rho_RE2
feas_RE2
[K_dt_RE2,rho_dt_RE2,feas_dt_RE2]=LMI_DT_H2(F,Gdec,Hdec,N,ContStruc_REBiStar,Gw);
rho_dt_RE2
feas_dt_RE2

eigenvaluesCT = eig(A+B*K_RE2)
spectral_abscissa = round(max(real(eig(A))),10)

eigenvaluesDT = eig(F+G*K_dt_RE2)
spectral_radius = round(max(abs(eig(F))),10)

%% OPTIMIZED structure
% doesnt make much sense graphically speaking. Apart from the abscissa, the
% rest of the eig is in general worse than the Reinforced Bi Star
% Nevertheless this implementation gets the best rho!

%Reinforced Bidirectional Star topology
ContStruc_OPT= [ 1 1 0 0 1
                 0 1 1 0 0
                 1 1 1 1 1
                 0 0 0 1 0  %%PROBLEM: optimization outputs a matrix with no inputs on the 4th system! 
                 1 0 0 1 1];

CFM_OPT=di_fixed_modes(A,Bdec,Cdec,N,ContStruc_OPT, rounding_n)
%no fixed modes in Continuous
DFM_OPT=di_fixed_modes(F,Gdec,Hdec,N,ContStruc_OPT, rounding_n)
%no fixed modes in Discrete

%Continuous Time stabilizing gains
[K_OPT,rho_OPT,feas_OPT]=LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStruc_OPT);
rho_OPT 
feas_OPT
%Discrete Time stabilizing gains
[K_dt_OPT,rho_dt_OPT,feas_dt_OPT]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStruc_OPT);
rho_dt_OPT
feas_dt_OPT

eigenvaluesCT = eig(A+B*K_OPT)
spectral_abscissa = round(max(real(eig(A))),10)

eigenvaluesDT = eig(F+G*K_dt_OPT)
spectral_radius = round(max(abs(eig(F))),10)

%CT and DT H2 gains
[K_OPT2,rho_OPT2,feas_OPT2]=LMI_CT_H2(A,Bdec,Cdec,N,ContStruc_OPT,Bw);
rho_OPT2
feas_OPT2
[K_dt_OPT2,rho_dt_OPT2,feas_dt_OPT2]=LMI_DT_H2(F,Gdec,Hdec,N,ContStruc_OPT,Gw);
rho_dt_OPT2
feas_dt_OPT2

eigenvaluesCT = eig(A+B*K_OPT2)
spectral_abscissa = round(max(real(eig(A))),10)

eigenvaluesDT = eig(F+G*K_dt_OPT2)
spectral_radius = round(max(abs(eig(F))),10)

%% Bidirectional Cyclic structure
% for the total system this is the one that makes most sense in my opinion,
% as there is no clear hierarchical division between the subsystem. In case
% one wants to choose  a specific subsystem as the one with most
% production, the Reinforced Bi Star is more suitable
ContStruc_cyc= [ 1 1 0 0 1
                 1 1 1 0 0
                 0 1 1 1 0
                 0 0 1 1 1   
                 1 0 0 1 1];

CFM_cyc=di_fixed_modes(A,Bdec,Cdec,N,ContStruc_cyc, rounding_n)
%no fixed modes in Continuous
DFM_cyc=di_fixed_modes(F,Gdec,Hdec,N,ContStruc_cyc, rounding_n)
%no fixed modes in Discrete

%Continuous Time stabilizing gains
[K_cyc,rho_cyc,feas_cyc]=LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStruc_cyc);
rho_cyc
feas_cyc
%Discrete Time stabilizing gains
[K_dt_cyc,rho_dt_cyc,feas_dt_cyc]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStruc_cyc);
rho_dt_cyc
feas_dt_cyc

eigenvaluesCT = eig(A+B*K_cyc)
spectral_abscissa = round(max(real(eig(A))),10)

eigenvaluesDT = eig(F+G*K_dt_cyc)
spectral_radius = round(max(abs(eig(F))),10)

%CT and DT H2 gains
[K_cyc2,rho_cyc2,feas_cyc2]=LMI_CT_H2(A,Bdec,Cdec,N,ContStruc_cyc,Bw);
rho_cyc2
feas_cyc2
[K_dt_cyc2,rho_dt_cyc2,feas_dt_cyc2]=LMI_DT_H2(F,Gdec,Hdec,N,ContStruc_cyc,Gw);
rho_dt_cyc2
feas_dt_cyc2

eigenvaluesCT = eig(A+B*K_cyc2)
spectral_abscissa = round(max(real(eig(A))),10)

eigenvaluesDT = eig(F+G*K_dt_cyc2)
spectral_radius = round(max(abs(eig(F))),10)

