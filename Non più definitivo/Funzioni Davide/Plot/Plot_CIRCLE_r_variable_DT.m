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
%CTot is an identity! all states are measurable!

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
Tfinal=30;
T=0:0.01:Tfinal;

[K_C_DT_1,rho_C_DT_1,feas_C_DT_1]=LMI_DT_EIG_CIRCLE_VARIABLE(F,Gdec,Hdec,N,ContStrucC, 0.5, 0.2);
[K_C_DT_2,rho_C_DT_2,feas_C_DT_2]=LMI_DT_EIG_CIRCLE_VARIABLE(F,Gdec,Hdec,N,ContStrucC, 0.5, 0.3);
[K_C_DT_3,rho_C_DT_3,feas_C_DT_3]=LMI_DT_EIG_CIRCLE_VARIABLE(F,Gdec,Hdec,N,ContStrucC, 0.5, 0.4);

k=0;
for k=1:Tfinal/h
    x_C_DISK_DT_1(:,k)= ((F+G*K_C_DT_1)^k)*x0;
    x_C_DISK_DT_2(:,k)= ((F+G*K_C_DT_2)^k)*x0;
    x_C_DISK_DT_3(:,k)= ((F+G*K_C_DT_3)^k)*x0;
end

 figure
 for i=1:N
        subplot(N,1,1+((i-1)))
        hold on
        grid on
        title(['\Delta\omega_{',num2str(i),'}_{,Centralized}'])
        plot([0:h:Tfinal],[x0(i*4-2),x_C_DISK_DT_1((i)*4-(4-2),:)],'m')
        plot([0:h:Tfinal],[x0(i*4-2),x_C_DISK_DT_2((i)*4-(4-2),:)],'r')
        plot([0:h:Tfinal],[x0(i*4-2),x_C_DISK_DT_3((i)*4-(4-2),:)],'g')
        legend('radius = 0.2','radius = 0.3','radius = 0.4')
 end

%% Decentralized
[K_De_DT_1,rho_De_DT_1,feas_De_DT_1]=LMI_DT_EIG_CIRCLE_VARIABLE(F,Gdec,Hdec,N,ContStrucDe, 0.5, 0.2);
[K_De_DT_2,rho_De_DT_2,feas_De_DT_2]=LMI_DT_EIG_CIRCLE_VARIABLE(F,Gdec,Hdec,N,ContStrucDe, 0.5, 0.3);
[K_De_DT_3,rho_De_DT_3,feas_De_DT_3]=LMI_DT_EIG_CIRCLE_VARIABLE(F,Gdec,Hdec,N,ContStrucDe, 0.5, 0.4);

k=0;
for k=1:Tfinal/h
    x_De_DISK_DT_1(:,k)= ((F+G*K_De_DT_1)^k)*x0;
    x_De_DISK_DT_2(:,k)= ((F+G*K_De_DT_2)^k)*x0;
    x_De_DISK_DT_3(:,k)= ((F+G*K_De_DT_3)^k)*x0;
end

 figure
 for i=1:N
        subplot(N,1,1+((i-1)))
        hold on
        grid on
        title(['\Delta\omega_{',num2str(i),'}_{,Decentralized}'])
        plot([0:h:Tfinal],[x0(i*4-2),x_De_DISK_DT_1((i)*4-(4-2),:)],'m')
        plot([0:h:Tfinal],[x0(i*4-2),x_De_DISK_DT_2((i)*4-(4-2),:)],'r')
        plot([0:h:Tfinal],[x0(i*4-2),x_De_DISK_DT_3((i)*4-(4-2),:)],'g')
        legend('radius = 0.2','radius = 0.3','radius = 0.4')
 end

%% Distributed 1
[K_Di1_DT_1,rho_Di1_DT_1,feas_Di1_DT_1]=LMI_DT_EIG_CIRCLE_VARIABLE(F,Gdec,Hdec,N,ContStruc_Pij, 0.5, 0.2);
[K_Di1_DT_2,rho_Di1_DT_2,feas_Di1_DT_2]=LMI_DT_EIG_CIRCLE_VARIABLE(F,Gdec,Hdec,N,ContStruc_Pij, 0.5, 0.3);
[K_Di1_DT_3,rho_Di1_DT_3,feas_Di1_DT_3]=LMI_DT_EIG_CIRCLE_VARIABLE(F,Gdec,Hdec,N,ContStruc_Pij, 0.5, 0.4);

k=0;
for k=1:Tfinal/h
    x_Di1_DISK_DT_1(:,k)= ((F+G*K_Di1_DT_1)^k)*x0;
    x_Di1_DISK_DT_2(:,k)= ((F+G*K_Di1_DT_2)^k)*x0;
    x_Di1_DISK_DT_3(:,k)= ((F+G*K_Di1_DT_3)^k)*x0;
end

 figure
 for i=1:N
        subplot(N,1,1+((i-1)))
        hold on
        grid on
        title(['\Delta\omega_{',num2str(i),'}_{,Distributed 1}'])
        plot([0:h:Tfinal],[x0(i*4-2),x_Di1_DISK_DT_1((i)*4-(4-2),:)],'m')
        plot([0:h:Tfinal],[x0(i*4-2),x_Di1_DISK_DT_2((i)*4-(4-2),:)],'r')
        plot([0:h:Tfinal],[x0(i*4-2),x_Di1_DISK_DT_3((i)*4-(4-2),:)],'g')
        legend('radius = 0.2','radius = 0.3','radius = 0.4')
 end

%% Distributed 2
[K_Di2_DT_1,rho_Di2_DT_1,feas_Di2_DT_1]=LMI_DT_EIG_CIRCLE_VARIABLE(F,Gdec,Hdec,N,ContStruc_REBiStar, 0.5, 0.2);
[K_Di2_DT_2,rho_Di2_DT_2,feas_Di2_DT_2]=LMI_DT_EIG_CIRCLE_VARIABLE(F,Gdec,Hdec,N,ContStruc_REBiStar, 0.5, 0.3);
[K_Di2_DT_3,rho_Di2_DT_3,feas_Di2_DT_3]=LMI_DT_EIG_CIRCLE_VARIABLE(F,Gdec,Hdec,N,ContStruc_REBiStar, 0.5, 0.4);

k=0;
for k=1:Tfinal/h
    x_Di2_DISK_DT_1(:,k)= ((F+G*K_Di2_DT_1)^k)*x0;
    x_Di2_DISK_DT_2(:,k)= ((F+G*K_Di2_DT_2)^k)*x0;
    x_Di2_DISK_DT_3(:,k)= ((F+G*K_Di2_DT_3)^k)*x0;
end

 figure
 for i=1:N
        subplot(N,1,1+((i-1)))
        hold on
        grid on
        title(['\Delta\omega_{',num2str(i),'}_{,Distributed 2}'])
        plot([0:h:Tfinal],[x0(i*4-2),x_Di2_DISK_DT_1((i)*4-(4-2),:)],'m')
        plot([0:h:Tfinal],[x0(i*4-2),x_Di2_DISK_DT_2((i)*4-(4-2),:)],'r')
        plot([0:h:Tfinal],[x0(i*4-2),x_Di2_DISK_DT_3((i)*4-(4-2),:)],'g')
        legend('radius = 0.2','radius = 0.3','radius = 0.4')
 end
%% Eigenvalues
 plot_eig_DT(F+G*K_C_DT_1, 0.5, 0.2)
 plot_eig_DT(F+G*K_C_DT_2, 0.5, 0.3)
 plot_eig_DT(F+G*K_C_DT_3, 0.5, 0.4)

 plot_eig_DT(F+G*K_De_DT_1, 0.5, 0.2)
 plot_eig_DT(F+G*K_De_DT_2, 0.5, 0.3)
 plot_eig_DT(F+G*K_De_DT_3, 0.5, 0.4)

 plot_eig_DT(F+G*K_Di1_DT_1, 0.5, 0.2)
 plot_eig_DT(F+G*K_Di1_DT_2, 0.5, 0.3)
 plot_eig_DT(F+G*K_Di1_DT_3, 0.5, 0.4)

 plot_eig_DT(F+G*K_Di2_DT_1, 0.5, 0.2)
 plot_eig_DT(F+G*K_Di2_DT_2, 0.5, 0.3)
 plot_eig_DT(F+G*K_Di2_DT_3, 0.5, 0.4)