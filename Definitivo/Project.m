% Power Network system
% Marcello Farina, 19/12/2018
clc
clearvars
close all

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
D=[];
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

[F,G,H,L,h]=ssdata(c2d(ss(A,B,C,[]),h));


for i=1:N
    Bdec{i}=B(:,i);
    Cdec{i}=C(n_states*(i-1)+1:n_states*i,:);
    Gdec{i}=G(:,i);
    Hdec{i}=H(n_states*(i-1)+1:n_states*i,:);
end

%% spectral abscissa and radius
eigenvaluesCT = eig(A)
spectral_abscissa = round(max(real(eig(A))),10)
plot_eig_CT(A, 0, -1)

eigenvaluesDT = eig(F)
spectral_radius = round(max(abs(eig(F))),10)
plot_eig_DT(F ,0, -1)

%the system is simply stable as we have only one eig on the border (zero for CT , 1 for DT),
%in both cases
%% Open-loop Free Motion simulation CT

Tfinal=60;
T=[0:0.01:Tfinal];
% x0=[];
% for i=1:N
%     x0=[x0;randn(n_states,1)];     %random vector of initial states
% end
load x0.mat;
k=0;
for t=T
    k=k+1;
     x_OpenLoop_CT(:,k)=expm(A*t)*x0;
end
     for v=1:n_states
         switch v
            case 1
                figure
                for i=1:N
                    subplot(N,1,i)
                    hold on
                    grid on
                    title(['\Delta\theta_{',num2str(i),'}'])
                    plot(T,[x_OpenLoop_CT((i)*4-(4-v),:)],'k')
                    legend('Open-loop free motion')
                end
            case 2
                figure
                for i=1:N
                    subplot(N,1,i)
                    hold on
                    grid on
                    title(['\Delta\omega_{',num2str(i),'}'])
                    plot(T,[x_OpenLoop_CT((i)*4-(4-v),:)],'k')
                    legend('Open-loop free motion')
                end
            case 3
                figure
                for i=1:N
                    subplot(N,1,i)
                    hold on
                    grid on
                    title(['\DeltaP_{m,',num2str(i),'}'])
                    plot(T,[x_OpenLoop_CT((i)*4-(4-v),:)],'k')
                    legend('Open-loop free motion')
                end
            case 4
                figure
                for i=1:N
                    subplot(N,1,i)
                    hold on
                    grid on
                    title(['\DeltaP_{v,',num2str(i),'}'])
                    plot(T,[x_OpenLoop_CT((i)*4-(4-v),:)],'k')
                    legend('Open-loop free motion')
                end
        end
     end   
%% Open-loop Free Motion simulation DT

for k=1:Tfinal/h
     x_OpenLoop_DT(:,k)=(F^k)*x0;
end
     for v=1:n_states
         switch v
            case 1
                figure
                for i=1:N
                    subplot(N,1,i)
                    hold on
                    grid on
                    title(['\Delta\theta_{',num2str(i),'}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_OpenLoop_DT((i)*4-(4-v),:)],'k')
                    legend('Open-loop free motion')
                end
            case 2
                figure
                for i=1:N
                    subplot(N,1,i)
                    hold on
                    grid on
                    title(['\Delta\omega_{',num2str(i),'}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_OpenLoop_DT((i)*4-(4-v),:)],'k')
                    legend('Open-loop free motion')
                end
            case 3
                figure
                for i=1:N
                    subplot(N,1,i)
                    hold on
                    grid on
                    title(['\DeltaP_{m,',num2str(i),'}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_OpenLoop_DT((i)*4-(4-v),:)],'k')
                    legend('Open-loop free motion')
                end
            case 4
                figure
                for i=1:N
                    subplot(N,1,i)
                    hold on
                    grid on
                    title(['\DeltaP_{v,',num2str(i),'}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_OpenLoop_DT((i)*4-(4-v),:)],'k')
                    legend('Open-loop free motion')
                end
        end
     end   
%% Centralized Fixed Modes
rounding_n=3;
ContStrucC=ones(N,N);
%[ 1 1 1 1 1;
%  1 1 1 1 1;
%  1 1 1 1 1;
%  1 1 1 1 1;
%  1 1 1 1 1];  %centralized structure
[CFMC]=di_fixed_modes(A,Bdec,Cdec,N,ContStrucC, rounding_n) % no centralized fixed modes CT
[DFMC]=di_fixed_modes(F,Gdec,Hdec,N,ContStrucC, rounding_n) % no centralized fixed modes DT

%% Decentralized Fixed Modes
rounding_n=3;
ContStrucDe=diag(ones(N,1));
%[ 1 0 0 0 0;
%  0 1 0 0 0;
%  0 0 1 0 0;
%  0 0 0 1 0;
%  0 0 0 0 1];  %decentralized structure
[CFMDe]=di_fixed_modes(A,Bdec,Cdec,N,ContStrucDe, rounding_n) %No decentralized FM CT =NO FM in general
[DFMDe]=di_fixed_modes(F,Gdec,Hdec,N,ContStrucDe, rounding_n) %No dec FM DT

%% Distributed Fixed Modes
rounding_n=3;
ContStrucDi=   [ 1 1 0 0 1
                 0 1 1 0 0
                 1 1 1 1 1
                 0 0 0 1 0  %%PROBLEM: optimization outputs a matrix with no inputs on the 4th system! 
                 1 0 0 1 1];

[CFMDi]=di_fixed_modes(A,Bdec,Cdec,N,ContStrucDi, rounding_n) %as expected since no decentralized
[DFMDi]=di_fixed_modes(F,Gdec,Hdec,N,ContStrucDi, rounding_n) %as expected

%% #### Conituous time Control Gains ####
% Centralized
[K_C_CT,rho_C_CT,feas_C_CT]=LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStrucC); %control gains for stability only
[K_C2_CT,rho_C2_CT,feas_C2_CT]=LMI_CT_EIG_TRESH(A,Bdec,Cdec,N,ContStrucC); %control gains that put eigs of (A+B*K_c2) before -alpha
[K_C3_CT,rho_C3_CT,feas_C3_CT]=LMI_CT_EIG_CIRCLE(A,Bdec,Cdec,N,ContStrucC); %why unfeasible prob if it does what I want?
[K_C4_CT,rho_C4_CT,feas_C4_CT]=LMI_CT_REG_ALPHA_MINU(A,Bdec,Cdec,N,ContStrucC); %eig in region
[K_C5_CT,rho_C5_CT,feas_C5_CT]=LMI_CT_H2(A,Bdec,Cdec,N,ContStrucC); %minimize H2 norm[K_c,rho_c,feas_c]=LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStrucC); %control gains for stability only

% Decentralized
[K_De_CT,rho_De_CT,feas_De_CT]=LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStrucDe); %control gains for stability only
[K_De2_CT,rho_De2_CT,feas_De2_CT]=LMI_CT_EIG_TRESH(A,Bdec,Cdec,N,ContStrucDe); %control gains that put eigs of (A+B*K_c2) before -alpha
[K_De3_CT,rho_De3_CT,feas_De3_CT]=LMI_CT_EIG_CIRCLE(A,Bdec,Cdec,N,ContStrucDe); %why unfeasible prob if it does what I want?
[K_De4_CT,rho_De4_CT,feas_De4_CT]=LMI_CT_REG_ALPHA_MINU(A,Bdec,Cdec,N,ContStrucDe); %eig in region
[K_De5_CT,rho_De5_CT,feas_De5_CT]=LMI_CT_H2(A,Bdec,Cdec,N,ContStrucDe); %minimize H2 norm[K_c,rho_c,feas_c]=LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStrucC); %control gains for stability only

% Distributed
[K_Di_CT,rho_Di_CT,feas_Di_CT]=LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStrucDi); %control gains for stability only
[K_Di2_CT,rho_Di2_CT,feas_Di2_CT]=LMI_CT_EIG_TRESH(A,Bdec,Cdec,N,ContStrucDi); %control gains that put eigs of (A+B*K_c2) before -alpha
[K_Di3_CT,rho_Di3_CT,feas_Di3_CT]=LMI_CT_EIG_CIRCLE(A,Bdec,Cdec,N,ContStrucDi); %why unfeasible prob if it does what I want?
[K_Di4_CT,rho_Di4_CT,feas_Di4_CT]=LMI_CT_REG_ALPHA_MINU(A,Bdec,Cdec,N,ContStrucDi); %eig in region
[K_Di5_CT,rho_Di5_CT,feas_Di5_CT]=LMI_CT_H2(A,Bdec,Cdec,N,ContStrucDi); %minimize H2 norm[K_c,rho_c,feas_c]=LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStrucC); %control gains for stability only

%%  #### Discrete time Control Gains ####
% Centralized
[K_C_DT,rho_C_DT,feas_C_DT]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStrucC); %control gains for stability only
[K_C2_DT,rho_C2_DT,feas_C2_DT]=LMI_DT_EIG_TRESH(F,Gdec,Hdec,N,ContStrucC); %control gains that put eigs of (A+B*K_c2) before -alpha
[K_C3_DT,rho_C3_DT,feas_C3_DT]=LMI_DT_EIG_CIRCLE(F,Gdec,Hdec,N,ContStrucC); %why unfeasible prob if it does what I want?
[K_C5_DT,rho_C5_DT,feas_C5_DT]=LMI_DT_H2(F,Gdec,Hdec,N,ContStrucC); %minimize H2 norm[K_c,rho_c,feas_c]=LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStrucC); %control gains for stability only

% Decentralized
[K_De_DT,rho_De_DT,feas_De_DT]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStrucDe); %control gains for stability only
[K_De2_DT,rho_De2_DT,feas_De2_DT]=LMI_DT_EIG_TRESH(F,Gdec,Hdec,N,ContStrucDe); %control gains that put eigs of (A+B*K_c2) before -alpha
[K_De3_DT,rho_De3_DT,feas_De3_DT]=LMI_DT_EIG_CIRCLE(F,Gdec,Hdec,N,ContStrucDe); %why unfeasible prob if it does what I want?
[K_De5_DT,rho_De5_DT,feas_De5_DT]=LMI_DT_H2(F,Gdec,Hdec,N,ContStrucDe); %minimize H2 norm[K_c,rho_c,feas_c]=LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStrucC); %control gains for stability only

% Distributed
[K_Di_DT,rho_Di_DT,feas_Di_DT]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStrucDi); %control gains for stability only
[K_Di2_DT,rho_Di2_DT,feas_Di2_DT]=LMI_DT_EIG_TRESH(F,Gdec,Hdec,N,ContStrucDi); %control gains that put eigs of (A+B*K_c2) before -alpha
[K_Di3_DT,rho_Di3_DT,feas_Di3_DT]=LMI_DT_EIG_CIRCLE(F,Gdec,Hdec,N,ContStrucDi); %why unfeasible prob if it does what I want?
[K_Di5_DT,rho_Di5_DT,feas_Di5_DT]=LMI_DT_H2(F,Gdec,Hdec,N,ContStrucDi); %minimize H2 norm[K_c,rho_c,feas_c]=LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStrucC); %control gains for stability only

%% ----- Simulation Continuous Time -----
k=0;
for t=T
    k=k+1;
    % equations of centralized free movement
    x_C_STABILITY_CT(:,k)=expm((A+B*K_C_CT)*t)*x0;
    x_C_TRESH_CT(:,k)=expm((A+B*K_C2_CT)*t)*x0;
    x_C_DISK_CT(:,k)=expm((A+B*K_C3_CT)*t)*x0;
    x_C_REGION_CT(:,k)=expm((A+B*K_C4_CT)*t)*x0;
    x_C_H2_CT(:,k)=expm((A+B*K_C5_CT)*t)*x0;
    
    % equations of decentralized free movement
    x_De_STABILITY_CT(:,k)=expm((A+B*K_De_CT)*t)*x0;
    x_De_TRESH_CT(:,k)=expm((A+B*K_De2_CT)*t)*x0;
    x_De_DISK_CT(:,k)=expm((A+B*K_De3_CT)*t)*x0;
    x_De_REGION_CT(:,k)=expm((A+B*K_De4_CT)*t)*x0;
    x_De_H2_CT(:,k)=expm((A+B*K_De5_CT)*t)*x0;
    
    % equations of distributed free movement
    x_Di_STABILITY_CT(:,k)=expm((A+B*K_Di_CT)*t)*x0;
    x_Di_TRESH_CT(:,k)=expm((A+B*K_Di2_CT)*t)*x0;
    x_Di_DISK_CT(:,k)=expm((A+B*K_Di3_CT)*t)*x0;
    x_Di_REGION_CT(:,k)=expm((A+B*K_Di4_CT)*t)*x0;
    x_Di_H2_CT(:,k)=expm((A+B*K_Di5_CT)*t)*x0;
end

%now compute the control actions:
    u_C_STABILITY_CT = K_C_CT*x_C_STABILITY_CT;
    u_C_TRESH_CT = K_C2_CT*x_C_TRESH_CT;
    u_C_DISK_CT = K_C3_CT*x_C_DISK_CT;
    u_C_REGION_CT = K_C4_CT*x_C_REGION_CT;
    u_C_H2_CT = K_C5_CT*x_C_H2_CT;

    u_De_STABILITY_CT = K_De_CT*x_De_STABILITY_CT;
    u_De_TRESH_CT = K_De2_CT*x_De_TRESH_CT;
    u_De_DISK_CT = K_De3_CT*x_De_DISK_CT;
    u_De_REGION_CT = K_De4_CT*x_De_REGION_CT;
    u_De_H2_CT = K_De5_CT*x_De_H2_CT;

    u_Di_STABILITY_CT = K_Di_CT*x_Di_STABILITY_CT;
    u_Di_TRESH_CT = K_Di2_CT*x_Di_TRESH_CT;
    u_Di_DISK_CT = K_Di3_CT*x_Di_DISK_CT;
    u_Di_REGION_CT = K_Di4_CT*x_Di_REGION_CT;
    u_Di_H2_CT = K_Di5_CT*x_Di_H2_CT;

%% ----- PLOTS Continuous time -----

%% Stability CT
for v=1:n_states
    switch v
         case 1
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\theta_{',num2str(i),'}_{,C}'])
                    plot(T,[x_C_STABILITY_CT((i)*4-(4-v),:)],'k')
                    legend('Stability Centralized CT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\theta_{',num2str(i),'}_{,De}'])
                    plot(T,[x_De_STABILITY_CT((i)*4-(4-v),:)],'k')
                    legend('Stability Decentralized CT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\theta_{',num2str(i),'}_{,Di}'])
                    plot(T,[x_Di_STABILITY_CT((i)*4-(4-v),:)],'k')
                    legend('Stability Distributed CT')
             end
         case 2
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\omega_{',num2str(i),'}_{,C}'])
                    plot(T,[x_C_STABILITY_CT((i)*4-(4-v),:)],'k')
                    legend('Stability Centralized CT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\omega_{m,',num2str(i),'}_{,De}'])
                    plot(T,[x_De_STABILITY_CT((i)*4-(4-v),:)],'k')
                    legend('Stability Decentralized CT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\omega_{',num2str(i),'}_{,Di}'])
                    plot(T,[x_Di_STABILITY_CT((i)*4-(4-v),:)],'k')
                    legend('Stability Distributed CT')
             end
              
         case 3
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{m,',num2str(i),'}_{,C}'])
                    plot(T,[x_C_STABILITY_CT((i)*4-(4-v),:)],'k')
                    legend('Stability Centralized CT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{m,',num2str(i),'}_{,De}'])
                    plot(T,[x_De_STABILITY_CT((i)*4-(4-v),:)],'k')
                    legend('Stability Decentralized CT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{m,',num2str(i),'}_{,Di}'])
                    plot(T,[x_Di_STABILITY_CT((i)*4-(4-v),:)],'k')
                    legend('Stability Distributed CT')
             end
         case 4
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{v,',num2str(i),'}_{,C}'])
                    plot(T,[x_C_STABILITY_CT((i)*4-(4-v),:)],'k')
                    legend('Stability Centralized CT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{v,',num2str(i),'}_{,De}'])
                    plot(T,[x_De_STABILITY_CT((i)*4-(4-v),:)],'k')
                    legend('Stability Decentralized CT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{v,',num2str(i),'}_{,Di}'])
                    plot(T,[x_Di_STABILITY_CT((i)*4-(4-v),:)],'k')
                    legend('Stability Distributed CT')
             end
    end
end   


figure  %control actions 
for i=1:N
    subplot(N,3,1+(3*(i-1)))
    hold on
    grid on
    title(['u_{',num2str(i),'}_{,C}'])
    plot(T,[u_C_STABILITY_CT(i,:)],'k')
    legend('Stability Centralized CT')

    subplot(N,3,2+(3*(i-1)))
    hold on
    grid on
    title(['u_{',num2str(i),'}_{,De}'])
    plot(T,[u_De_STABILITY_CT(i,:)],'k')
    legend('Stability Decentralized CT')

    subplot(N,3,3+(3*(i-1)))
    hold on
    grid on
    title(['u_{',num2str(i),'}_{,Di}'])
    plot(T,[u_Di_STABILITY_CT(i,:)],'k')
    legend('Stability Distributed CT')
end

%% Tresh CT
for v=1:n_states
    switch v
         case 1
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\theta_{',num2str(i),'}_{,C}'])
                    plot(T,[x_C_TRESH_CT((i)*4-(4-v),:)],'k')
                    legend('Tresh Eig Centralized CT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\theta_{',num2str(i),'}_{,De}'])
                    plot(T,[x_De_TRESH_CT((i)*4-(4-v),:)],'k')
                    legend('Tresh Eig Decentralized CT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\theta_{',num2str(i),'}_{,Di}'])
                    plot(T,[x_Di_TRESH_CT((i)*4-(4-v),:)],'k')
                    legend('Tresh Eig Distributed CT')
             end
         case 2
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\omega_{',num2str(i),'}_{,C}'])
                    plot(T,[x_C_TRESH_CT((i)*4-(4-v),:)],'k')
                    legend('Tresh Eig Centralized CT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\omega_{m,',num2str(i),'}_{,De}'])
                    plot(T,[x_De_TRESH_CT((i)*4-(4-v),:)],'k')
                    legend('Tresh Eig Decentralized CT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\omega_{',num2str(i),'}_{,Di}'])
                    plot(T,[x_Di_TRESH_CT((i)*4-(4-v),:)],'k')
                    legend('Tresh Eig Distributed CT')
             end
         case 3
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{m,',num2str(i),'}_{,C}'])
                    plot(T,[x_C_TRESH_CT((i)*4-(4-v),:)],'k')
                    legend('Tresh Eig Centralized CT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{m,',num2str(i),'}_{,De}'])
                    plot(T,[x_De_TRESH_CT((i)*4-(4-v),:)],'k')
                    legend('Tresh Eig Decentralized CT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{m,',num2str(i),'}_{,Di}'])
                    plot(T,[x_Di_TRESH_CT((i)*4-(4-v),:)],'k')
                    legend('Tresh Eig Distributed CT')
             end
         case 4
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{v,',num2str(i),'}_{,C}'])
                    plot(T,[x_C_TRESH_CT((i)*4-(4-v),:)],'k')
                    legend('Tresh Eig Centralized CT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{v,',num2str(i),'}_{,De}'])
                    plot(T,[x_De_TRESH_CT((i)*4-(4-v),:)],'k')
                    legend('Tresh Eig Decentralized CT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{v,',num2str(i),'}_{,Di}'])
                    plot(T,[x_Di_TRESH_CT((i)*4-(4-v),:)],'k')
                    legend('Tresh Eig Distributed CT')
             end
    end
end   

figure
 for i=1:N
        subplot(N,3,1+(3*(i-1)))
        hold on
        grid on
        title(['u_{',num2str(i),'}_{,C}'])
        plot(T,[u_C_TRESH_CT(i,:)],'k')
        legend('Tresh Eig Centralized CT')
    
        subplot(N,3,2+(3*(i-1)))
        hold on
        grid on
        title(['u_{',num2str(i),'}_{,De}'])
        plot(T,[u_De_TRESH_CT(i,:)],'k')
        legend('Tresh Eig Decentralized CT')
    
        subplot(N,3,3+(3*(i-1)))
        hold on
        grid on
        title(['u_{',num2str(i),'}_{,Di}'])
        plot(T,[u_Di_TRESH_CT(i,:)],'k')
        legend('Tresh Eig Distributed CT')
 end

%% Disk CT
for v=1:n_states
    switch v
         case 1
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\theta_{',num2str(i),'}_{,C}'])
                    plot(T,[x_C_DISK_CT((i)*4-(4-v),:)],'k')
                    legend('Disk Eig Centralized CT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\theta_{',num2str(i),'}_{,De}'])
                    plot(T,[x_De_DISK_CT((i)*4-(4-v),:)],'k')
                    legend('Disk Eig Decentralized CT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\theta_{',num2str(i),'}_{,Di}'])
                    plot(T,[x_Di_DISK_CT((i)*4-(4-v),:)],'k')
                    legend('Disk Eig Distributed CT')
             end
         case 2
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\omega_{',num2str(i),'}_{,C}'])
                    plot(T,[x_C_DISK_CT((i)*4-(4-v),:)],'k')
                    legend('Disk Eig Centralized CT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                   title(['\Delta\omega_{m,',num2str(i),'}_{,De}'])
                    plot(T,[x_De_DISK_CT((i)*4-(4-v),:)],'k')
                    legend('Disk Eig Decentralized CT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\omega_{',num2str(i),'}_{,Di}'])
                    plot(T,[x_Di_DISK_CT((i)*4-(4-v),:)],'k')
                    legend('Disk Eig Distributed CT')
             end
         case 3
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{m,',num2str(i),'}_{,C}'])
                    plot(T,[x_C_DISK_CT((i)*4-(4-v),:)],'k')
                    legend('Disk Eig Centralized CT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{m,',num2str(i),'}_{,De}'])
                    plot(T,[x_De_DISK_CT((i)*4-(4-v),:)],'k')
                    legend('Disk Eig Decentralized CT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{m,',num2str(i),'}_{,Di}'])
                    plot(T,[x_Di_DISK_CT((i)*4-(4-v),:)],'k')
                    legend('Disk Eig Distributed CT')
             end
         case 4
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{v,',num2str(i),'}_{,C}'])
                    plot(T,[x_C_DISK_CT((i)*4-(4-v),:)],'k')
                    legend('Disk Eig Centralized CT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{v,',num2str(i),'}_{,De}'])
                    plot(T,[x_De_DISK_CT((i)*4-(4-v),:)],'k')
                    legend('Disk Eig Decentralized CT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{v,',num2str(i),'}_{,Di}'])
                    plot(T,[x_Di_DISK_CT((i)*4-(4-v),:)],'k')
                    legend('Disk Eig Distributed CT')
             end
    end
end

figure
 for i=1:N
        subplot(N,3,1+(3*(i-1)))
        hold on
        grid on
        title(['u_{',num2str(i),'}_{,C}'])
        plot(T,[u_C_DISK_CT(i,:)],'k')
        legend('Disk Eig Centralized CT')
    
        subplot(N,3,2+(3*(i-1)))
        hold on
        grid on
        title(['u_{',num2str(i),'}_{,De}'])
        plot(T,[u_De_DISK_CT(i,:)],'k')
        legend('Disk Eig Decentralized CT')
    
        subplot(N,3,3+(3*(i-1)))
        hold on
        grid on
        title(['u_{',num2str(i),'}_{,Di}'])
        plot(T,[u_Di_DISK_CT(i,:)],'k')
        legend('Disk Eig Distributed CT')
 end

%% Region CT
for v=1:n_states
    switch v
         case 1
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\theta_{',num2str(i),'}_{,C}'])
                    plot(T,[x_C_REGION_CT((i)*4-(4-v),:)],'k')
                    legend('Region Eig Centralized CT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\theta_{',num2str(i),'}_{,De}'])
                    plot(T,[x_De_REGION_CT((i)*4-(4-v),:)],'k')
                    legend('Region Eig Decentralized CT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\theta_{',num2str(i),'}_{,Di}'])
                    plot(T,[x_Di_REGION_CT((i)*4-(4-v),:)],'k')
                    legend('Region Eig Distributed CT')
             end
         case 2
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\omega_{',num2str(i),'}_{,C}'])
                    plot(T,[x_C_REGION_CT((i)*4-(4-v),:)],'k')
                    legend('Region Eig Centralized CT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\omega_{m,',num2str(i),'}_{,De}'])
                    plot(T,[x_De_REGION_CT((i)*4-(4-v),:)],'k')
                    legend('Region Eig Decentralized CT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\omega_{',num2str(i),'}_{,Di}'])
                    plot(T,[x_Di_REGION_CT((i)*4-(4-v),:)],'k')
                    legend('Region Eig Distributed CT')
             end
         case 3
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{m,',num2str(i),'}_{,C}'])
                    plot(T,[x_C_REGION_CT((i)*4-(4-v),:)],'k')
                    legend('Region Eig Centralized CT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{m,',num2str(i),'}_{,De}'])
                    plot(T,[x_De_REGION_CT((i)*4-(4-v),:)],'k')
                    legend('Region Eig Decentralized CT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{m,',num2str(i),'}_{,Di}'])
                    plot(T,[x_Di_REGION_CT((i)*4-(4-v),:)],'k')
                    legend('Region Eig Distributed CT')
             end
         case 4
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{v,',num2str(i),'}_{,C}'])
                    plot(T,[x_C_REGION_CT((i)*4-(4-v),:)],'k')
                    legend('Region Eig Centralized CT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{v,',num2str(i),'}_{,De}'])
                    plot(T,[x_De_REGION_CT((i)*4-(4-v),:)],'k')
                    legend('Region Eig Decentralized CT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{v,',num2str(i),'}_{,Di}'])
                    plot(T,[x_Di_REGION_CT((i)*4-(4-v),:)],'k')
                    legend('Region Eig Distributed CT')
             end
    end
end

figure
 for i=1:N
        subplot(N,3,1+(3*(i-1)))
        hold on
        grid on
        title(['u_{',num2str(i),'}_{,C}'])
        plot(T,[u_C_REGION_CT(i,:)],'k')
        legend('Region Eig Centralized CT')
    
        subplot(N,3,2+(3*(i-1)))
        hold on
        grid on
        title(['u_{',num2str(i),'}_{,De}'])
        plot(T,[u_De_REGION_CT(i,:)],'k')
        legend('Region Eig Decentralized CT')
    
        subplot(N,3,3+(3*(i-1)))
        hold on
        grid on
        title(['u_{u',num2str(i),'}_{,Di}'])
        plot(T,[u_Di_REGION_CT(i,:)],'k')
        legend('Region Eig Distributed CT')
 end

%% H2 CT
for v=1:n_states
    switch v
         case 1
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\theta_{',num2str(i),'}_{,C}'])
                    plot(T,[x_C_H2_CT((i)*4-(4-v),:)],'k')
                    legend('H2 Centralized CT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\theta_{',num2str(i),'}_{,De}'])
                    plot(T,[x_De_H2_CT((i)*4-(4-v),:)],'k')
                    legend('H2 Decentralized CT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\theta_{',num2str(i),'}_{,Di}'])
                    plot(T,[x_Di_H2_CT((i)*4-(4-v),:)],'k')
                    legend('H2 Distributed CT')
             end
         case 2
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\omega_{',num2str(i),'}_{,C}'])
                    plot(T,[x_C_H2_CT((i)*4-(4-v),:)],'k')
                    legend('H2 Centralized CT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\omega_{m,',num2str(i),'}_{,De}'])
                    plot(T,[x_De_H2_CT((i)*4-(4-v),:)],'k')
                    legend('H2 Decentralized CT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\omega_{',num2str(i),'}_{,Di}'])
                    plot(T,[x_Di_H2_CT((i)*4-(4-v),:)],'k')
                    legend('H2 Distributed CT')
             end
         case 3
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{m,',num2str(i),'}_{,C}'])
                    plot(T,[x_C_H2_CT((i)*4-(4-v),:)],'k')
                    legend('H2 Centralized CT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{m,',num2str(i),'}_{,De}'])
                    plot(T,[x_De_H2_CT((i)*4-(4-v),:)],'k')
                    legend('H2 Decentralized CT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{m,',num2str(i),'}_{,Di}'])
                    plot(T,[x_Di_H2_CT((i)*4-(4-v),:)],'k')
                    legend('H2 Distributed CT')
             end
         case 4
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{v,',num2str(i),'}_{,C}'])
                    plot(T,[x_C_H2_CT((i)*4-(4-v),:)],'k')
                    legend('H2 Centralized CT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{v,',num2str(i),'}_{,De}'])
                    plot(T,[x_De_H2_CT((i)*4-(4-v),:)],'k')
                    legend('H2 Decentralized CT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{v,',num2str(i),'}_{,Di}'])
                    plot(T,[x_Di_H2_CT((i)*4-(4-v),:)],'k')
                    legend('H2 Distributed CT')
             end
    end
end

figure
 for i=1:N
        subplot(N,3,1+(3*(i-1)))
        hold on
        grid on
        title(['u_{',num2str(i),'}_{,C}'])
        plot(T,[u_C_H2_CT(i,:)],'k')
        legend('H2 Centralized CT')
    
        subplot(N,3,2+(3*(i-1)))
        hold on
        grid on
        title(['u_{',num2str(i),'}_{,De}'])
        plot(T,[u_De_H2_CT(i,:)],'k')
        legend('H2 Decentralized CT')
    
        subplot(N,3,3+(3*(i-1)))
        hold on
        grid on
        title(['u_{',num2str(i),'}_{,Di}'])
        plot(T,[u_Di_H2_CT(i,:)],'k')
        legend('H2 Distributed CT')
 end

%% ----- Simulation Discrete Time -----
for k=1:Tfinal/h
    % equations of centralized free movement
    x_C_STABILITY_DT(:,k)=((F+G*K_C_DT)^k)*x0;
    x_C_TRESH_DT(:,k)=((F+G*K_C2_DT)^k)*x0;
    x_C_DISK_DT(:,k)=((F+G*K_C3_DT)^k)*x0;
    x_C_H2_DT(:,k)=((F+G*K_C5_DT)^k)*x0;
    
    % equations of decentralized free movement
    x_De_STABILITY_DT(:,k)=((F+G*K_De_DT)^k)*x0;
    x_De_TRESH_DT(:,k)=((F+G*K_De2_DT)^k)*x0;
    x_De_DISK_DT(:,k)=((F+G*K_De3_DT)^k)*x0;
    x_De_H2_DT(:,k)=((F+G*K_De5_DT)^k)*x0;
    
    % equations of distributed free movement
    x_Di_STABILITY_DT(:,k)=((F+G*K_Di_DT)^k)*x0;
    x_Di_TRESH_DT(:,k)=((F+G*K_Di2_DT)^k)*x0;
    x_Di_DISK_DT(:,k)=((F+G*K_Di3_DT)^k)*x0;
    x_Di_H2_DT(:,k)=((F+G*K_Di5_DT)^k)*x0;
end

% now compute the control actions:

    u_C_STABILITY_DT = K_C_DT*x_C_STABILITY_DT;
    u_C_TRESH_DT = K_C2_DT*x_C_TRESH_DT;
    u_C_DISK_DT = K_C3_DT*x_C_DISK_DT;
    u_C_H2_DT = K_C5_DT*x_C_H2_DT;

    u_De_STABILITY_DT = K_De_DT*x_De_STABILITY_DT;
    u_De_TRESH_DT = K_De2_DT*x_De_TRESH_DT;
    u_De_DISK_DT = K_De3_DT*x_De_DISK_DT;
    u_De_H2_DT = K_De5_DT*x_De_H2_DT;

    u_Di_STABILITY_DT = K_Di_DT*x_Di_STABILITY_DT;
    u_Di_TRESH_DT = K_Di2_DT*x_Di_TRESH_DT;
    u_Di_DISK_DT = K_Di3_DT*x_Di_DISK_DT;
    u_Di_H2_DT = K_Di5_DT*x_Di_H2_DT;



%% ----- PLOTS Discrete Time -----

%% Stability DT
for v=1:n_states
    switch v
         case 1
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\theta_{',num2str(i),'}_{,C}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_C_STABILITY_DT((i)*4-(4-v),:)],'k')
                    legend('Stability Centralized DT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\theta_{',num2str(i),'}_{,De}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_De_STABILITY_DT((i)*4-(4-v),:)],'k')
                    legend('Stability Decentralized DT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\theta_{',num2str(i),'}_{,Di}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_Di_STABILITY_DT((i)*4-(4-v),:)],'k')
                    legend('Stability Distributed DT')
             end
         case 2
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\omega_{',num2str(i),'}_{,C}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_C_STABILITY_DT((i)*4-(4-v),:)],'k')
                    legend('Stability Centralized DT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\omega_{',num2str(i),'}_{,De}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_De_STABILITY_DT((i)*4-(4-v),:)],'k')
                    legend('Stability Decentralized DT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\omega_{',num2str(i),'}_{,Di}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_Di_STABILITY_DT((i)*4-(4-v),:)],'k')
                    legend('Stability Distributed DT')
             end
         case 3
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{m,',num2str(i),'}_{,C}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_C_STABILITY_DT((i)*4-(4-v),:)],'k')
                    legend('Stability Centralized DT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{m,',num2str(i),'}_{,De}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_De_STABILITY_DT((i)*4-(4-v),:)],'k')
                    legend('Stability Decentralized DT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{m,',num2str(i),'}_{,Di}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_Di_STABILITY_DT((i)*4-(4-v),:)],'k')
                    legend('Stability Distributed DT')
             end
         case 4
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{v,',num2str(i),'}_{,C}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_C_STABILITY_DT((i)*4-(4-v),:)],'k')
                    legend('Stability Centralized DT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{v,',num2str(i),'}_{,De}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_De_STABILITY_DT((i)*4-(4-v),:)],'k')
                    legend('Stability Decentralized DT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{v,',num2str(i),'}_{,Di}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_Di_STABILITY_DT((i)*4-(4-v),:)],'k')
                    legend('Stability Distributed DT')
             end
    end
end   

figure %control starts acting in k=1, not in k=0 !!!
 for i=1:N
        subplot(N,3,1+(3*(i-1)))
        hold on
        grid on
        title(['u_{',num2str(i),'}_{,C}'])
        plot([1:h:Tfinal],u_C_STABILITY_DT(i,:),'k')
        legend('Stability Centralized DT')
    
        subplot(N,3,2+(3*(i-1)))
        hold on
        grid on
        title(['u_{',num2str(i),'}_{,De}'])
        plot([1:h:Tfinal],u_De_STABILITY_DT(i,:),'k')
        legend('Stability Decentralized DT')
    
        subplot(N,3,3+(3*(i-1)))
        hold on
        grid on
        title(['u_{',num2str(i),'}_{,Di}'])
        plot([1:h:Tfinal],u_Di_STABILITY_DT(i,:),'k')
        legend('Stability Distributed DT')
 end


%% Tresh DT
for v=1:n_states
    switch v
         case 1
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\theta_{',num2str(i),'}_{,C}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_C_TRESH_DT((i)*4-(4-v),:)],'k')
                    legend('Tresh Eig Centralized DT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\theta_{',num2str(i),'}_{,De}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_De_TRESH_DT((i)*4-(4-v),:)],'k')
                    legend('Tresh Eig Decentralized DT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\theta_{',num2str(i),'}_{,Di}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_Di_TRESH_DT((i)*4-(4-v),:)],'k')
                    legend('Tresh Eig Distributed DT')
             end
         case 2
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\omega_{',num2str(i),'}_{,C}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_C_TRESH_DT((i)*4-(4-v),:)],'k')
                    legend('Tresh Eig Centralized DT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\omega_{',num2str(i),'}_{,De}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_De_TRESH_DT((i)*4-(4-v),:)],'k')
                    legend('Tresh Eig Decentralized DT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\omega_{',num2str(i),'}_{,Di}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_Di_TRESH_DT((i)*4-(4-v),:)],'k')
                    legend('Tresh Eig Distributed DT')
             end
         case 3
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{m,',num2str(i),'}_{,C}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_C_TRESH_DT((i)*4-(4-v),:)],'k')
                    legend('Tresh Eig Centralized DT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{m,',num2str(i),'}_{,De}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_De_TRESH_DT((i)*4-(4-v),:)],'k')
                    legend('Tresh Eig Decentralized DT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{m,',num2str(i),'}_{,Di}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_Di_TRESH_DT((i)*4-(4-v),:)],'k')
                    legend('Tresh Eig Distributed DT')
             end
         case 4
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{v,',num2str(i),'}_{,C}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_C_TRESH_DT((i)*4-(4-v),:)],'k')
                    legend('Tresh Eig Centralized DT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{v,',num2str(i),'}_{,De}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_De_TRESH_DT((i)*4-(4-v),:)],'k')
                    legend('Tresh Eig Decentralized DT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{v,',num2str(i),'}_{,Di}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_Di_TRESH_DT((i)*4-(4-v),:)],'k')
                    legend('Tresh Eig Distributed DT')
             end
    end
end   

figure
 for i=1:N
        subplot(N,3,1+(3*(i-1)))
        hold on
        grid on
        title(['u_{',num2str(i),'}_{,C}'])
        plot([1:h:Tfinal],u_C_TRESH_DT(i,:),'k')
        legend('Tresh Eig Centralized DT')
    
        subplot(N,3,2+(3*(i-1)))
        hold on
        grid on
        title(['u_{',num2str(i),'}_{,De}'])
        plot([1:h:Tfinal],u_De_TRESH_DT(i,:),'k')
        legend('Tresh Eig Decentralized DT')
    
        subplot(N,3,3+(3*(i-1)))
        hold on
        grid on
        title(['u_{',num2str(i),'}_{,Di}'])
        plot([1:h:Tfinal],u_Di_TRESH_DT(i,:),'k')
        legend('Tresh Eig Distributed DT')
 end
%% Disk DT
for v=1:n_states
    switch v
         case 1
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\theta_{',num2str(i),'}_{,C}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_C_DISK_DT((i)*4-(4-v),:)],'k')
                    legend('Disk Eig Centralized DT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\theta_{',num2str(i),'}_{,De}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_De_DISK_DT((i)*4-(4-v),:)],'k')
                    legend('Disk Eig Decentralized DT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\theta_{',num2str(i),'}_{,Di}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_Di_DISK_DT((i)*4-(4-v),:)],'k')
                    legend('Disk Eig Distributed DT')
             end
         case 2
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\omega_{',num2str(i),'}_{,C}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_C_DISK_DT((i)*4-(4-v),:)],'k')
                    legend('Disk Eig Centralized DT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                   title(['\Delta\omega_{m,',num2str(i),'}_{,De}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_De_DISK_DT((i)*4-(4-v),:)],'k')
                    legend('Disk Eig Decentralized DT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\omega_{',num2str(i),'}_{,Di}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_Di_DISK_DT((i)*4-(4-v),:)],'k')
                    legend('Disk Eig Distributed DT')
             end
         case 3
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{m,',num2str(i),'}_{,C}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_C_DISK_DT((i)*4-(4-v),:)],'k')
                    legend('Disk Eig Centralized DT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{m,',num2str(i),'}_{,De}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_De_DISK_DT((i)*4-(4-v),:)],'k')
                    legend('Disk Eig Decentralized DT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{m,',num2str(i),'}_{,Di}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_Di_DISK_DT((i)*4-(4-v),:)],'k')
                    legend('Disk Eig Distributed DT')
             end
         case 4
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{v,',num2str(i),'}_{,C}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_C_DISK_DT((i)*4-(4-v),:)],'k')
                    legend('Disk Eig Centralized DT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{v,',num2str(i),'}_{,De}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_De_DISK_DT((i)*4-(4-v),:)],'k')
                    legend('Disk Eig Decentralized DT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{v,',num2str(i),'}_{,Di}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_Di_DISK_DT((i)*4-(4-v),:)],'k')
                    legend('Disk Eig Distributed DT')
             end
    end
end

figure
 for i=1:N
        subplot(N,3,1+(3*(i-1)))
        hold on
        grid on
        title(['u_{',num2str(i),'}_{,C}'])
        plot([1:h:Tfinal],u_C_DISK_DT(i,:),'k')
        legend('Disk Eig Centralized DT')
    
        subplot(N,3,2+(3*(i-1)))
        hold on
        grid on
        title(['u_{',num2str(i),'}_{,De}'])
        plot([1:h:Tfinal],u_De_DISK_DT(i,:),'k')
        legend('Disk Eig Decentralized DT')
    
        subplot(N,3,3+(3*(i-1)))
        hold on
        grid on
        title(['u_{',num2str(i),'}_{,Di}'])
        plot([1:h:Tfinal],u_Di_DISK_DT(i,:),'k')
        legend('Disk Eig Distributed DT')
 end
%% H2 DT
for v=1:n_states
    switch v
         case 1
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\theta_{',num2str(i),'}_{,C}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_C_H2_DT((i)*4-(4-v),:)],'k')
                    legend('H2 Centralized DT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\theta_{',num2str(i),'}_{,De}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_De_H2_DT((i)*4-(4-v),:)],'k')
                    legend('H2 Decentralized DT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\theta_{',num2str(i),'}_{,Di}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_Di_H2_DT((i)*4-(4-v),:)],'k')
                    legend('H2 Distributed DT')
             end
         case 2
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\omega_{',num2str(i),'}_{,C}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_C_H2_DT((i)*4-(4-v),:)],'k')
                    legend('H2 Centralized DT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\omega_{',num2str(i),'}_{,De}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_De_H2_DT((i)*4-(4-v),:)],'k')
                    legend('H2 Decentralized DT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\Delta\omega_{',num2str(i),'}_{,Di}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_Di_H2_DT((i)*4-(4-v),:)],'k')
                    legend('H2 Distributed DT')
             end
         case 3
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{m,',num2str(i),'}_{,C}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_C_H2_DT((i)*4-(4-v),:)],'k')
                    legend('H2 Centralized DT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{m,',num2str(i),'}_{,De}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_De_H2_DT((i)*4-(4-v),:)],'k')
                    legend('H2 Decentralized DT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{m,',num2str(i),'}_{,Di}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_Di_H2_DT((i)*4-(4-v),:)],'k')
                    legend('H2 Distributed DT')
             end
         case 4
             figure
             for i=1:N
                    subplot(N,3,1+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{v,',num2str(i),'}_{,C}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_C_H2_DT((i)*4-(4-v),:)],'k')
                    legend('H2 Centralized DT')
                
                    subplot(N,3,2+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{v,',num2str(i),'}_{,De}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_De_H2_DT((i)*4-(4-v),:)],'k')
                    legend('H2 Decentralized DT')
                
                    subplot(N,3,3+(3*(i-1)))
                    hold on
                    grid on
                    title(['\DeltaP_{v,',num2str(i),'}_{,Di}'])
                    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_Di_H2_DT((i)*4-(4-v),:)],'k')
                    legend('H2 Distributed DT')
             end
    end
end

figure
 for i=1:N
        subplot(N,3,1+(3*(i-1)))
        hold on
        grid on
        title(['u_{',num2str(i),'}_{,C}'])
        plot([1:h:Tfinal],u_C_H2_DT(i,:),'k')
        legend('H2 Centralized DT')
    
        subplot(N,3,2+(3*(i-1)))
        hold on
        grid on
        title(['u_{',num2str(i),'}_{,De}'])
        plot([1:h:Tfinal],u_De_H2_DT(i,:),'k')
        legend('H2 Decentralized DT')
    
        subplot(N,3,3+(3*(i-1)))
        hold on
        grid on
        title(['u_{',num2str(i),'}_{,Di}'])
        plot([1:h:Tfinal],u_Di_H2_DT(i,:),'k')
        legend('H2 Distributed DT')
 end

%% Eigenvalues comparion on different control laws
eigenvaluesDT = eig(F+G*K_C_DT);
    figure()
    grid on
    hold on
    plot(real(eigenvaluesDT), imag(eigenvaluesDT), '*')
    disk = circle(0,0,1);
    hold off
    title('Eigenvalues')
    xlabel('Re')
    ylabel('Im')

    eigenvaluesDT = eig(F+G*K_C5_DT);
    figure()
    grid on
    hold on
    plot(real(eigenvaluesDT), imag(eigenvaluesDT), '*')
    disk = circle(0,0,1);
    hold off
    title('Eigenvalues')
    xlabel('Re')
    ylabel('Im')
%% Grafici sovrapposti sample
         % case 1
         %     figure
         %     for i=1:N
         %            subplot(N,1,i)
         %            hold on
         %            grid on
         %            title(['\Delta\theta_{',num2str(i),'}'])
         %            plot(T,[x_C_STABILITY_CT((i)*4-(4-v),:)],'k')
         %            plot(T,[x_De_STABILITY_CT((i)*4-(4-v),:)],'k')
         %            plot(T,[x_Di_STABILITY_CT((i)*4-(4-v),:)],'k')
         %            legend('Stability Centralized CT','Stability Decentralized CT','Stability Distributed CT')
         %     end         %     end