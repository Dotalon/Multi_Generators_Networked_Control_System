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

rounding_n=3;

ContStruc_Pij=[ 1 1 0 0 0
                 1 1 1 0 1
                 0 1 1 1 0
                 0 0 1 1 1 
                 0 1 0 1 1];

ContStruc_REBiStar=[ 1 1 1 1 1
                       1 1 0 0 1
                       1 1 1 0 0
                       1 0 1 1 0 
                       1 0 0 1 1];
ContStrucC=ones(N,N);
ContStrucDe=diag(ones(N,1));

load('x0.mat')
%% Centralized
Tfinal=7;
T=0:0.01:Tfinal;

[K_C_CT_1,rho_C_CT_1,feas_C_CT_1]=LMI_CT_EIG_CIRCLE_VARIABLE(A,Bdec,Cdec,N,ContStrucC, 5, 1);
[K_C_CT_2,rho_C_CT_2,feas_C_CT_2]=LMI_CT_EIG_CIRCLE_VARIABLE(A,Bdec,Cdec,N,ContStrucC, 5, 2);
[K_C_CT_3,rho_C_CT_3,feas_C_CT_3]=LMI_CT_EIG_CIRCLE_VARIABLE(A,Bdec,Cdec,N,ContStrucC, 5, 3);
[K_C_CT_4,rho_C_CT_4,feas_C_CT_4]=LMI_CT_EIG_CIRCLE_VARIABLE(A,Bdec,Cdec,N,ContStrucC, 5, 4);

k=0;
for t=T
    k=k+1;
    x_C_DISK_CT_1(:,k)=expm((A+B*K_C_CT_1)*t)*x0;
    x_C_DISK_CT_2(:,k)=expm((A+B*K_C_CT_2)*t)*x0;
    x_C_DISK_CT_3(:,k)=expm((A+B*K_C_CT_3)*t)*x0;
    x_C_DISK_CT_4(:,k)=expm((A+B*K_C_CT_4)*t)*x0;
end

 figure
 for i=1:N
        subplot(N,1,1+((i-1)))
        hold on
        grid on
        title(['\Delta\omega_{',num2str(i),'}_{,Centralized}'])
        plot(T,[x_C_DISK_CT_1((i)*4-(4-2),:)],'m')
        plot(T,[x_C_DISK_CT_2((i)*4-(4-2),:)],'r')
        plot(T,[x_C_DISK_CT_3((i)*4-(4-2),:)],'g')
        plot(T,[x_C_DISK_CT_4((i)*4-(4-2),:)],'b')
        legend('radius = 1','radius = 2','radius = 3','radius = 4')
 end

%% Decentralized
[K_De_CT_1,rho_De_CT_1,feas_De_CT_1]=LMI_CT_EIG_CIRCLE_VARIABLE(A,Bdec,Cdec,N,ContStrucDe, 5, 1);
[K_De_CT_2,rho_De_CT_2,feas_De_CT_2]=LMI_CT_EIG_CIRCLE_VARIABLE(A,Bdec,Cdec,N,ContStrucDe, 5, 2);
[K_De_CT_3,rho_De_CT_3,feas_De_CT_3]=LMI_CT_EIG_CIRCLE_VARIABLE(A,Bdec,Cdec,N,ContStrucDe, 5, 3);
[K_De_CT_4,rho_De_CT_4,feas_De_CT_4]=LMI_CT_EIG_CIRCLE_VARIABLE(A,Bdec,Cdec,N,ContStrucDe, 5, 4);

k=0;
for t=T
    k=k+1;
    x_De_DISK_CT_1(:,k)=expm((A+B*K_De_CT_1)*t)*x0;
    x_De_DISK_CT_2(:,k)=expm((A+B*K_De_CT_2)*t)*x0;
    x_De_DISK_CT_3(:,k)=expm((A+B*K_De_CT_3)*t)*x0;
    x_De_DISK_CT_4(:,k)=expm((A+B*K_De_CT_4)*t)*x0;
end

 figure
 for i=1:N
        subplot(N,1,1+((i-1)))
        hold on
        grid on
        title(['\Delta\omega_{',num2str(i),'}_{,Decentralized}'])
        plot(T,[x_De_DISK_CT_1((i)*4-(4-2),:)],'m')
        plot(T,[x_De_DISK_CT_2((i)*4-(4-2),:)],'r')
        plot(T,[x_De_DISK_CT_3((i)*4-(4-2),:)],'g')
        plot(T,[x_De_DISK_CT_4((i)*4-(4-2),:)],'b')
        legend('radius = 1','radius = 2','radius = 3','radius = 4')
 end

%% Distributed 1
[K_Di1_CT_1,rho_Di1_CT_1,feas_Di1_CT_1]=LMI_CT_EIG_CIRCLE_VARIABLE(A,Bdec,Cdec,N,ContStruc_Pij, 5, 1);
[K_Di1_CT_2,rho_Di1_CT_2,feas_Di1_CT_2]=LMI_CT_EIG_CIRCLE_VARIABLE(A,Bdec,Cdec,N,ContStruc_Pij, 5, 2);
[K_Di1_CT_3,rho_Di1_CT_3,feas_Di1_CT_3]=LMI_CT_EIG_CIRCLE_VARIABLE(A,Bdec,Cdec,N,ContStruc_Pij, 5, 3);
[K_Di1_CT_4,rho_Di1_CT_4,feas_Di1_CT_4]=LMI_CT_EIG_CIRCLE_VARIABLE(A,Bdec,Cdec,N,ContStruc_Pij, 5, 4);

k=0;
for t=T
    k=k+1;
    x_Di1_DISK_CT_1(:,k)=expm((A+B*K_Di1_CT_1)*t)*x0;
    x_Di1_DISK_CT_2(:,k)=expm((A+B*K_Di1_CT_2)*t)*x0;
    x_Di1_DISK_CT_3(:,k)=expm((A+B*K_Di1_CT_3)*t)*x0;
    x_Di1_DISK_CT_4(:,k)=expm((A+B*K_Di1_CT_4)*t)*x0;
end

 figure
 for i=1:N
        subplot(N,1,1+((i-1)))
        hold on
        grid on
        title(['\Delta\omega_{',num2str(i),'}_{,Distributed 1}'])
        plot(T,[x_Di1_DISK_CT_1((i)*4-(4-2),:)],'m')
        plot(T,[x_Di1_DISK_CT_2((i)*4-(4-2),:)],'r')
        plot(T,[x_Di1_DISK_CT_3((i)*4-(4-2),:)],'g')
        plot(T,[x_Di1_DISK_CT_4((i)*4-(4-2),:)],'b')
        legend('radius = 1','radius = 2','radius = 3','radius = 4')
 end

%% Distributed 2
[K_Di2_CT_1,rho_Di2_CT_1,feas_Di2_CT_1]=LMI_CT_EIG_CIRCLE_VARIABLE(A,Bdec,Cdec,N,ContStruc_REBiStar, 5, 1);
[K_Di2_CT_2,rho_Di2_CT_2,feas_Di2_CT_2]=LMI_CT_EIG_CIRCLE_VARIABLE(A,Bdec,Cdec,N,ContStruc_REBiStar, 5, 2);
[K_Di2_CT_3,rho_Di2_CT_3,feas_Di2_CT_3]=LMI_CT_EIG_CIRCLE_VARIABLE(A,Bdec,Cdec,N,ContStruc_REBiStar, 5, 3);
[K_Di2_CT_4,rho_Di2_CT_4,feas_Di2_CT_4]=LMI_CT_EIG_CIRCLE_VARIABLE(A,Bdec,Cdec,N,ContStruc_REBiStar, 5, 4);

k=0;
for t=T
    k=k+1;
    x_Di2_DISK_CT_1(:,k)=expm((A+B*K_Di2_CT_1)*t)*x0;
    x_Di2_DISK_CT_2(:,k)=expm((A+B*K_Di2_CT_2)*t)*x0;
    x_Di2_DISK_CT_3(:,k)=expm((A+B*K_Di2_CT_3)*t)*x0;
    x_Di2_DISK_CT_4(:,k)=expm((A+B*K_Di2_CT_4)*t)*x0;
end

 figure
 for i=1:N
        subplot(N,1,1+((i-1)))
        hold on
        grid on
        title(['\Delta\omega_{',num2str(i),'}_{,Distributed 2}'])
        plot(T,[x_Di2_DISK_CT_1((i)*4-(4-2),:)],'m')
        plot(T,[x_Di2_DISK_CT_2((i)*4-(4-2),:)],'r')
        plot(T,[x_Di2_DISK_CT_3((i)*4-(4-2),:)],'g')
        plot(T,[x_Di2_DISK_CT_4((i)*4-(4-2),:)],'b')
        legend('radius = 1','radius = 2','radius = 3','radius = 4')
 end
%% Eigenvalues
 plot_eig_CT(A+B*K_C_CT_1, 5, 1)
 plot_eig_CT(A+B*K_C_CT_2, 5, 2)
 plot_eig_CT(A+B*K_C_CT_3, 5, 3)
 plot_eig_CT(A+B*K_C_CT_4, 5, 4)

 plot_eig_CT(A+B*K_De_CT_1, 5, 1)
 plot_eig_CT(A+B*K_De_CT_2, 5, 2)
 plot_eig_CT(A+B*K_De_CT_3, 5, 3)
 plot_eig_CT(A+B*K_De_CT_4, 5, 4)

 plot_eig_CT(A+B*K_Di1_CT_1, 5, 1)
 plot_eig_CT(A+B*K_Di1_CT_2, 5, 2)
 plot_eig_CT(A+B*K_Di1_CT_3, 5, 3)
 plot_eig_CT(A+B*K_Di1_CT_4, 5, 4)

 plot_eig_CT(A+B*K_Di2_CT_1, 5, 1)
 plot_eig_CT(A+B*K_Di2_CT_2, 5, 2)
 plot_eig_CT(A+B*K_Di2_CT_3, 5, 3)
 plot_eig_CT(A+B*K_Di2_CT_4, 5, 4)