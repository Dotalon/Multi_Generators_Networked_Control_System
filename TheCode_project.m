% Power Network system
% Marcello Farina, 19/12/2018

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

eigenvaluesDT = eig(F)
spectral_radius = round(max(abs(eig(F))),10)

%the system is simply stable as we have only one eig on the border (zero for CT , 1 for DT),
%in both cases

%% centralized fixed modes
rounding_n=3;
ContStrucC=ones(N,N);
%[ 1 1 1 1 1
%  1 1 1 1 1
%  1 1 1 1 1
%  1 1 1 1 1];  %centralized structure
[CFMC]=di_fixed_modes(A,Bdec,Cdec,N,ContStrucC, rounding_n) % no centralized fixed modes CT
[DFMC]=di_fixed_modes(F,Gdec,Hdec,N,ContStrucC, rounding_n) % no centralized fixed modes DT

[K_c,rho_c,feas_c]=LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStrucC); %control gains for stability only
[K_c2,rho_c2,feas_c2]=LMI_CT_EIG_TRESH(A,Bdec,Cdec,N,ContStrucC); %control gains that put eigs of (A+B*K_c2) before -alpha
[K_c3,rho_c3,feas_c3]=LMI_CT_CIRCLE_EIG(A,Bdec,Cdec,N,ContStrucC); %why unfeasible prob if it does what I want?
[K_c4,rho_c4,feas_c4]=LMI_CT_REGION(A,Bdec,Cdec,N,ContStrucC); %eig in region
[K_c5,rho_c5,feas_c5]=LMI_CT_H2(A,Bdec,Cdec,N,ContStrucC); %minimize H2 norm

%Can't make Hinf work

%% decentralized fixed modes
ContStrucDe=diag(ones(N,1));
[CFMDe]=di_fixed_modes(A,Bdec,Cdec,N,ContStrucDe, rounding_n) %No decentralized FM CT =NO FM in general
[DFMDe]=di_fixed_modes(F,Gdec,Hdec,N,ContStrucDe, rounding_n) %No dec FM DT

%Same LMIs as before
[K_c_DE,rho_c_DE,feas_c_DE]=LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStrucDe);
[K_c2_DE,rho_c2_DE,feas_c2_DE]=LMI_CT_EIG_TRESH(A,Bdec,Cdec,N,ContStrucDe);
[K_c3_DE,rho_c3_DE,feas_c3_DE]=LMI_CT_CIRCLE_EIG(A,Bdec,Cdec,N,ContStrucDe);
[K_c4_DE,rho_c4_DE,feas_c4_DE]=LMI_CT_REGION(A,Bdec,Cdec,N,ContStrucDe);
[K_c5_DE,rho_c5_DE,feas_c5_DE]=LMI_CT_H2(A,Bdec,Cdec,N,ContStrucDe);


%% distributed fixed modes
ContStrucDi=[ 1 1 0 0 0
              1 1 1 0 1
              0 1 1 1 0
              0 0 1 1 1
              0 1 0 1 1];
[CFMDi]=di_fixed_modes(A,Bdec,Cdec,N,ContStrucDi, rounding_n) %as expected since no decentralized
[DFMDi]=di_fixed_modes(F,Gdec,Hdec,N,ContStrucDi, rounding_n) %as expected

%Same LMIs as before
[K_c_Di,rho_c_Di,feas_c_Di]=LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStrucDi);
[K_c2_Di,rho_c2_Di,feas_c2_Di]=LMI_CT_EIG_TRESH(A,Bdec,Cdec,N,ContStrucDi);
[K_c3_Di,rho_c3_Di,feas_c3_Di]=LMI_CT_CIRCLE_EIG(A,Bdec,Cdec,N,ContStrucDi);
[K_c4_Di,rho_c4_Di,feas_c4_Di]=LMI_CT_REGION(A,Bdec,Cdec,N,ContStrucDi);
[K_c5_Di,rho_c5_Di,feas_c5_Di]=LMI_CT_H2(A,Bdec,Cdec,N,ContStrucDi);


% %% Apply control gains for stability only
% [K_c,rho_c,feas_c]=LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStrucC); %control gains for stability only
% [K_c_DT,rho_c_DT,feas_c_DT]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStrucC);
% eig(F+G*K_c_DT)
% 
% %% alpha-stability
% [K_c2,rho_c2,feas_c2]=LMI_CT_EIG_TRESH(A,Bdec,Cdec,N,ContStrucC); %control gains that put eigs of (A+B*K_c2) before -alpha
% 
% [K_c2_DT,rho_c2_DT,feas_c2_DT]=LMI_DT_EIG_CIRCLE(F,Gdec,Hdec,N,ContStrucC);
% 
% 
% %% Confine eigs in a circle 
% [K_c3,rho_c3,feas_c3]=LMI_CT_CIRCLE_EIG(A,Bdec,Cdec,N,ContStrucC); %why unfeasible prob if it does what I want?
% 
% %[eigenvalues_wrong] = eig(A+B*K_c3);
% %plot([eigenvalues_wrong])
% 
% %No discrete time cause it doesn't make sense
% 
% %% Confine eigs in region WORKS
% [K_c4,rho_c4,feas_c4]=LMI_CT_REGION(A,Bdec,Cdec,N,ContStrucC);
% 
% %No discrete time cause it doesn't make sense
% 
% %% H2 norm minimization
% [K_c5,rho_c5,feas_c5]=LMI_CT_H2(A,Bdec,Cdec,N,ContStrucC);
% 
% % [K_c5_DT,rho_c5_DT,feas_c5_DT]=LMI_DT_H2(F,Gdec,Hdec,N,ContStrucC);
% 
% 
% %% Hinf non funziona, serve guardare
% [K_c6,rho_c6,feas_c6]=LMI_CT_DeDicont_Hinf(A,Bdec,Cdec,N,ContStrucC);
% 


%% don't actually know why this part is here but I'm scared of deleting it
Gtot=[];
Htot=[];
Btot=[];
Ctot=[];
for i=1:N
    Btot=[Btot,Bdec{i}];
    Ctot=[Ctot
        Cdec{i}];
    Gtot=[Gtot,Gdec{i}];
    Htot=[Htot
        Hdec{i}];
end



%% simulation
%   Display on every figure 3 cols: centralized, decentralized and distributed
%   schemes, for now in continuous time (CT)

Tfinal=45;
T=[0:0.01:Tfinal];
x0=[];
for i=1:N
    x0=[x0;randn(n_states,1)];     %random col vector of initial states
end

k=0;
for t=T
    k=k+1;
    
    % equations of centralized free movement
    x_c(:,k)=expm((A+B*K_c)*t)*x0;
    x_alpha_stab(:,k)=expm((A+B*K_c2)*t)*x0;
    x_disc(:,k)=expm((A+B*K_c3)*t)*x0;
    x_region(:,k)=expm((A+B*K_c4)*t)*x0;
    x_H2(:,k)=expm((A+B*K_c5)*t)*x0;
    
    % equations of decentralized free movement
    x_simple_stab_DE(:,k)=expm((A+B*K_c_DE)*t)*x0;
    x_alpha_stab_DE(:,k)=expm((A+B*K_c2_DE)*t)*x0;
    x_disc_DE(:,k)=expm((A+B*K_c3_DE)*t)*x0;
    x_region_DE(:,k)=expm((A+B*K_c4_DE)*t)*x0;
    x_H2_DE(:,k)=expm((A+B*K_c5_DE)*t)*x0;
    
    % equations of distributed free movement
    x_simple_stab_Di(:,k)=expm((A+B*K_c_Di)*t)*x0;
    x_alpha_stab_Di(:,k)=expm((A+B*K_c2_Di)*t)*x0;
    x_disc_Di(:,k)=expm((A+B*K_c3_Di)*t)*x0;
    x_region_Di(:,k)=expm((A+B*K_c4_Di)*t)*x0;
    x_H2_Di(:,k)=expm((A+B*K_c5_Di)*t)*x0;
end

%%CODE TO USE FOR DT
% for k=1:Tfinal/h
% %     x_c_DT(:,k)=expm((F+G*K_c_DT)^k)*x0;
%     x_c_DT(:,k)=((F+G*K_c_DT)^k)*x0;
%     x_circle_DT(:,k)=((F+G*K_c2_DT)^k)*x0;
% 
% end

%     subplot(N,2,2*(i))            copy this code when doing CentralizedDT
%     hold on
%     grid on
%     title(['\speed_DT{',num2str(i),'}'])
%     plot([0:h:Tfinal-1],[x_circle_DT((i)*4-2,:)],'k')
%% Plot of the movements for every control law used (simple stability, alpha-stab,ecc) 

figure
for i=1:N
    subplot(N,3,1+(3*(i-1)))
    hold on
    grid on
    title(['\speed_{',num2str(i),'}'])
    plot(T,[x_c((i)*4-2,:)],'k')
    legend('Simple stab CENTR CT')

    subplot(N,3,2+(3*(i-1)))
    hold on
    grid on
    title(['\speed_DE{',num2str(i),'}'])
    plot(T,[x_simple_stab_DE((i)*4-2,:)],'k')
    legend('simple stab DE CT')

    subplot(N,3,3+(3*(i-1)))
    hold on
    grid on
    title(['\speed_Di{',num2str(i),'}'])
    plot(T,[x_simple_stab_Di((i)*4-2,:)],'k')
    legend('simple stab Di CT')
   
end

figure 
for i=1:N
    subplot(N,3,1+(3*(i-1)))
    hold on
    grid on
    title(['\speed_{',num2str(i),'}'])
    plot(T,[x_alpha_stab((i)*4-2,:)],'k')
    legend('Alpha stab CENTR CT')

    subplot(N,3,2+(3*(i-1)))
    hold on
    grid on
    title(['\speed_DE{',num2str(i),'}'])
    plot(T,[x_alpha_stab_DE((i)*4-2,:)],'k')
    legend('Alpha stab DE CT')
    
    subplot(N,3,3+(3*(i-1)))
    hold on
    grid on
    title(['\speed_Di{',num2str(i),'}'])
    plot(T,[x_alpha_stab_Di((i)*4-2,:)],'k')
    legend('alpha stab Di CT')

end

figure
for i=1:N
    subplot(N,3,1+3*(i-1))
    hold on
    grid on
    title(['\speed_{',num2str(i),'}'])
    plot(T,[x_disc((i)*4-2,:)],'k')
    legend('Eig in disc CENTR CT')

    subplot(N,3,2+(3*(i-1)))
    hold on
    grid on
    title(['\speed_DE{',num2str(i),'}'])
    plot(T,[x_disc_DE((i)*4-2,:)],'k')
    legend('Eig in disc DE CT')

    subplot(N,3,3+(3*(i-1)))
    hold on
    grid on
    title(['\speed_Di{',num2str(i),'}'])
    plot(T,[x_disc_Di((i)*4-2,:)],'k')
    legend('eig in disc Di CT')
end


figure
for i=1:N
    subplot(N,3,1+3*(i-1))
    hold on
    grid on
    title(['\speed_{',num2str(i),'}'])
    plot(T,[x_region((i)*4-2,:)],'k')
    legend('Eig in region CENTR CT')

    subplot(N,3,2+(3*(i-1)))
    hold on
    grid on
    title(['\speed_DE{',num2str(i),'}'])
    plot(T,[x_region_DE((i)*4-2,:)],'k')
    legend('eig in region DE CT')

    subplot(N,3,3+(3*(i-1)))
    hold on
    grid on
    title(['\speed_Di{',num2str(i),'}'])
    plot(T,[x_region_Di((i)*4-2,:)],'k')
    legend('eig in region Di CT')
end


figure
for i=1:N
    subplot(N,3,1+(3*(i-1)))
    hold on
    grid on
    title(['\speed_{',num2str(i),'}'])
    plot(T,[x_H2((i)*4-2,:)],'k')
    legend('H2 centr')


    subplot(N,3,2+(3*(i-1)))
    hold on
    grid on
    title(['\speed_DE{',num2str(i),'}'])
    plot(T,[x_H2_DE((i)*4-2,:)],'k')
    legend('H2 DE CT')
    
    subplot(N,3,3+(3*(i-1)))
    hold on
    grid on
    title(['\speed_Di{',num2str(i),'}'])
    plot(T,[x_H2_Di((i)*4-2,:)],'k')
    legend('H2 Di CT')
end






%% decentralized fixed modes
ContStrucDe=diag(ones(N,1));
[CFMDe]=di_fixed_modes(A,Bdec,Cdec,N,ContStrucDe, rounding_n)
[DFMDe]=di_fixed_modes(F,Gdec,Hdec,N,ContStrucDe, rounding_n)

[K_c_DE,rho_c_DE,feas_c_DE]=LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStrucDe);
[K_c2_DE,rho_c2_DE,feas_c2_DE]=LMI_CT_EIG_TRESH(A,Bdec,Cdec,N,ContStrucDe);
[K_c5_DE,rho_c5_DE,feas_c5_DE]=LMI_CT_H2(A,Bdec,Cdec,N,ContStrucDe);


%% distributed fixed modes
ContStrucDi=[ 1 1 0 0 0
              1 1 1 0 1
              0 1 1 1 0
              0 0 1 1 1
              0 1 0 1 1];
[CFMDi]=di_fixed_modes(A,Bdec,Cdec,N,ContStrucDi, rounding_n)
[DFMDi]=di_fixed_modes(F,Gdec,Hdec,N,ContStrucDi, rounding_n)

%% distributed with another structure
ContStrucDi_2=[ 1 0 0 0 0
                0 0 0 0 0
                0 0 0 0 0
                0 0 0 0 0
                0 0 0 0 1];
[CFMDi_2]=di_fixed_modes(A,Bdec,Cdec,N,ContStrucDi, rounding_n)
[DFMDi_2]=di_fixed_modes(F,Gdec,Hdec,N,ContStrucDi, rounding_n)