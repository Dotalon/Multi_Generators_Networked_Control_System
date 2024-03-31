clc
clearvars
close all
%%
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
D=zeros(20,5);
N_dt = zeros(20,5);
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

% plot_eig_CT(A, 0, -1)

eigenvaluesDT = eig(F)
spectral_radius = round(max(abs(eig(F))),10)

% plot_eig_DT(F ,0, -1)

%the system is simply stable as we have only one eig on the border [zero for CT (but one is really close to 0), 1 for DT],
%in both cases

%% centralized fixed modes
rounding_n=3;
ContStrucC=ones(N,N);
%[ 1 1 1 1 1
%  1 1 1 1 1
%  1 1 1 1 1
%  1 1 1 1 1];  %centralized structure
[CFMC]=di_fixed_modes(A,Bdec,Cdec,N,ContStrucC, rounding_n)
[DFMC]=di_fixed_modes(F,Gdec,Hdec,N,ContStrucC, rounding_n)

% no centralized fixed modes
%% decentralized fixed modes
ContStrucDe=diag(ones(N,1));
[CFMDe]=di_fixed_modes(A,Bdec,Cdec,N,ContStrucDe, rounding_n)
[DFMDe]=di_fixed_modes(F,Gdec,Hdec,N,ContStrucDe, rounding_n)

%% distributed fixed modes
ContStrucDi=[ 1 1 0 0 0
              1 1 1 0 1
              0 1 1 1 0
              0 1 0 1 1
              0 1 0 1 1];
[CFMDi]=di_fixed_modes(A,Bdec,Cdec,N,ContStrucDi, rounding_n)
[DFMDi]=di_fixed_modes(F,Gdec,Hdec,N,ContStrucDi, rounding_n)

[K_C_CT,rho_C_CT,feas_C_CT]=LMI_CT_DeDicont_Hinf(A,Bdec,Cdec,N,ContStrucC);
sys_ss_CT=ss(A+B*K_C_CT, B, C, D);
sys_CT=tf(sys_ss_CT);
bode(sys_CT)
[mag_CT,phase_CT,w_CT] = bode(sys_CT);
h_CT=hinfnorm(sys_CT)

[K_C_DT,rho_C_DT,feas_C_DT]=LMI_DT_DeDicont_Hinf(F,Gdec,Hdec,N,ContStrucC);
sys_ss_DT=ss(F+G*K_C_DT, G, H, N_dt);
sys_DT=tf(sys_ss_DT);
bode(sys_DT)
[mag_DT,phase_DT,w_DT] = bode(sys_DT);
h_DT=hinfnorm(sys_DT)

% 
% %% Decentralized Control
% % Continuous-time Stability
% [K_De_C,rho_De_C,feas_De_C]=LMI_DeDicont(A,Bdec,Cdec,N,ContStrucDe, 'CT');
% spectral_abscissa_alpha = round(max(real(eig(A+B*K_De_C))),10)
% plot_eig_CT(A+B*K_De_C)
% K_De_C_Stability = K_De_C;
% % Continuous-time alpha-Stability
% alpha = input("Value of alpha (positive):")
% [K_De_C,rho_De_C,feas_De_C]=LMI_DeDicont(A,Bdec,Cdec,N,ContStrucDe, 'CT', alpha);
% spectral_abscissa_alpha = round(max(real(eig(A+B*K_De_C))),10)
% plot_eig_CT(A+B*K_De_C, alpha)
% K_De_C_Alpha_Stability =K_De_C;
% % Continuous-time Eigenvalues in a disk (alpha=0 for limitations on spectral radius)
% alpha = input("Value of alpha (positive):")
% radius = input("Value of radius (positive):")
% [K_De_C,rho_De_C,feas_De_C]=LMI_DeDicont(A,Bdec,Cdec,N,ContStrucDe, 'CT', alpha, radius);
% spectral_abscissa_alpha = round(max(real(eig(A+B*K_De_C))),10)
% plot_eig_CT(A+B*K_De_C, alpha, radius)
% K_De_C_Disk = K_De_C;
% % Discrete-time stability
% [K_De_DT,rho_De_DT,feas_De_DT]=LMI_DeDicont(F,Gdec,Hdec,N,ContStrucDe, 'DT');
% spectral_radius_rho = round(max(abs(eig(F+G*K_De_DT))),10)
% plot_eig_DT(F+G*K_De_DT)
% K_De_DT_Stability = K_De_DT;
% % Discrete-time Eigenvalues in a disk (alpha=0 for limitations on spectral radius)
% alpha = input("Value of alpha (positive):")
% radius = input("Value of radius (positive):")
% [K_De_DT,rho_De_DT,feas_De_DT]=LMI_DeDicont(F,Gdec,Hdec,N,ContStrucDe, 'DT', alpha, radius);
% spectral_radius_rho = round(max(abs(eig(F+G*K_De_DT))),10)
% plot_eig_DT(F+G*K_De_DT, alpha, radius)
% K_De_DT_Disk = K_De_DT;
% %%
% %%%%%%%%%%%%%
% %   Plots   %
% %%%%%%%%%%%%%
% 
% j = input('Which K_De_C do you want to use? (1 for Stability, 2 for Alpha stability, 3 for Disk)')
% if j == 1
%     K_De_C = K_De_C_Stability;
% elseif j == 2
%     K_De_C = K_De_C_Alpha_Stability;
% elseif j == 3
%     K_De_C = K_De_C_Disk;
% end
% 
% j = input('Which K_De_DT do you want to use? (1 for Stability, 2 for Disk)')
% 
% if j == 1
%     K_De_DT = K_De_DT_Stability;
% elseif j == 2
%     K_De_DT = K_De_DT_Disk;
% end
% Gtot=[];
% Htot=[];
% Btot=[];
% Ctot=[];
% for i=1:N
%     Btot=[Btot,Bdec{i}];
%     Ctot=[Ctot
%         Cdec{i}];
%     Gtot=[Gtot,Gdec{i}];
%     Htot=[Htot
%         Hdec{i}];
% end
% 
% % simulation data
% Tfinal=6;
% T=[0:0.01:Tfinal];
% x0=[];
% for i=1:N
%     x0=[x0;randn(n_states,1)];
% end
% Atot = A;
% Ftot = F;
% k=0;
% for t=T
%     k=k+1;
%     x_De_C(:,k)=expm((Atot+Btot*K_De_C)*t)*x0;
% end
% for k=1:Tfinal/h
%     x_De_DT(:,k)=((Ftot+Gtot*K_De_DT)^k)*x0;
% end
% for v=1:n_states
%     switch v
%         case 1
%             state = '\Delta\theta_{'
%             plot_freemotion(N,v,state, T, x_De_C, x_De_DT, x0, h, Tfinal)
%         case 2
%             state = '\Delta\omega_{'
%             plot_freemotion(N,v,state, T, x_De_C, x_De_DT, x0, h, Tfinal)
%         case 3
%             state = '\DeltaP_{m,'
%             plot_freemotion(N,v,state, T, x_De_C, x_De_DT, x0, h, Tfinal)
%         case 4
%             state = '\DeltaP_{v,'
%             plot_freemotion(N,v,state, T, x_De_C, x_De_DT, x0, h, Tfinal)
%     end
% end
%% Functions
% 
% function [K,rho,feas]=LMI_DeDicont(A,B,C,N,ContStruc, Mode, alpha, radius)
% % Computes, using LMIs, the distributed "state feedback" control law for the continuous-time system, with reference to the control
% % information structure specified by 'ContStruc'.
% %
% % Inputs:
% % - A: system matrix.
% % - B: input matrices (i.e., B{1},..., B{N} are the input matrices of the decomposed system, one for each channel).
% % - C: output matrices  (i.e., C{1},..., C{N} are the output matrices of the decomposed system, one for each channel, where [Cdec{1}',...,
% % Cdec{N}']=I).
% % - N: number of subsystems.
% % - ContStruc: NxN matrix that specifies the information structure
% % constraints (ContStruc(i,j)=1 if communication is allowed between channel
% % j to channel i, ContStruc(i,j)=0 otherwise).
% % - Mode: CT for Continuos Time, DT for Discrete Time
% % - Alpha: Center of the disk region. If there is no input for alpha, it is set to 0. If there is no input for radius 
% %   alpha is the limit value for the spectral abscissa instead of the center of a disk
% % - Radius: Radius of of the disk region.
% % Output:
% % - K: structured control gain
% % - rho: spectral abscissa of matrix (A+B*K) - note that [C{1}',...,
% % C{N}']=I
% % - feas: feasibility of the LMI problem (=0 if yes)
% % - B_total: Matrix to calculate the eigenvalues of the closed loop (used for the plots)
% if Mode == 'CT'
%     if nargin < 7
%      alpha = 0;
%     end
% 
%     if nargin < 8
%      radius = -1;
%     end
% 
%     Btot=[];
%     for i=1:N
%     m(i)=size(B{i},2);
%     n(i)=size(C{i},1);
%     Btot=[Btot,B{i}];
% end
%     ntot=size(A,1);
%     mtot=sum(m);
% 
%     yalmip clear
% 
%     if ContStruc==ones(N,N)
%         % Centralized design
%         P=sdpvar(ntot);
%         L=sdpvar(mtot,ntot);
%     else
%         % Decentralized/distributed design
%         P=[];
%         L=sdpvar(mtot,ntot);
%         minc=0;
%         for i=1:N
%             P=blkdiag(P,sdpvar(n(i)));
%             ninc=0;
%             for j=1:N
%                 if ContStruc(i,j)==0
%                     L(minc+1:minc+m(i),ninc+1:ninc+n(j))=zeros(m(i),n(j));
%                 end
%                 ninc=ninc+n(j);
%             end
%             minc=minc+m(i);
%         end  
%     end
%    if radius==-1
%     LMIconstr=[P*A'+A*P+Btot*L+L'*Btot'+ 2*alpha*P<=-1e-2*eye(ntot)]+[P>=1e-2*eye(ntot)];
%    else
%     LMIconstr=[[(radius^2-alpha^2)*P-A*P*A'-A*L'*Btot'-Btot*L*A'-alpha*(P*A'+A*P+L'*Btot'+Btot*L), Btot*L;
%                  L'*Btot', P]>=1e-2*eye(ntot*2)];
%    end
%     options=sdpsettings('solver','sedumi');
%     J=optimize(LMIconstr,[],options);
%     feas=J.problem;
%     L=double(L);
%     P=double(P);
% 
%     K=L/P;
%     rho=max(real(eig(A+Btot*K)));
% elseif Mode == 'DT'
%     if nargin < 7
%      alpha = 0;
%      radius = 1;
%     end
% 
%     if nargin < 8
%      radius = 1;
%     end
% 
%     F = A;
%     G = B;
%     H = C;
% 
%     Gtot=[];
%     for i=1:N
%      m(i)=size(G{i},2);
%      n(i)=size(H{i},1);
%      Gtot=[Gtot,G{i}];
%     end
%     ntot=size(F,1);
%     mtot=sum(m);
% 
%     yalmip clear
% 
%     if ContStruc==ones(N,N)
%       % Centralized design
%       P=sdpvar(ntot);
%       L=sdpvar(mtot,ntot);
%     else
%       % Dentralized/distributed design
%       P=[];
%       L=sdpvar(mtot,ntot);
%       minc=0;
%       for i=1:N
%           P=blkdiag(P,sdpvar(n(i)));
%           ninc=0;
%          for j=1:N
%               if ContStruc(i,j)==0
%                   L(minc+1:minc+m(i),ninc+1:ninc+n(j))=zeros(m(i),n(j));
%                end
%                ninc=ninc+n(j);
%          end
%           minc=minc+m(i);
%      end
%     end
% 
% if radius==-1
%     LMIconstr=[[P-F*P*F'-F*L'*Gtot'-Gtot*L*F' Gtot*L;
%                L'*Gtot' P]>=1e-2*eye(ntot*2)];
%    else
%     LMIconstr=[[(radius^2-alpha^2)*P-F*P*F'-F*L'*Gtot'-Gtot*L*F'-alpha*(P*F'+F*P+L'*Gtot'+Gtot*L), Gtot*L;
%                  L'*Gtot', P]>=1e-2*eye(ntot*2)];
%    end
%     options=sdpsettings('solver','sedumi');
%     J=optimize(LMIconstr,[],options);
%     feas=J.problem;
%     L=double(L);
%     P=double(P);
% 
%     K=L/P;
%     rho=max(abs(eig(F+Gtot*K)));
% else print('Mode not vailid, please use CT for continuous time or DT for Discrete time')
% end
% end
% function [Difm]=di_fixed_modes(A,B,C,N,ContStruc,rounding_n)
% 
% Btot=[];
% Ctot=[];
% for i=1:N
%     m(i)=size(B{i},2);
%     p(i)=size(C{i},1);
%     Btot=[Btot,B{i}];
%     Ctot=[Ctot
%         C{i}];
% end
% 
% m_tot=size(Btot,2);
% p_tot=size(Ctot,1);
% 
% Difm=round(eig(A),rounding_n);
% nelD=length(Difm);
% 
% kend=1000;
% k=0;
% while (nelD~=0)&&(k<=kend)
%     k=k+1;
%     K=zeros(m_tot,p_tot);
%     m_inc=0;
%     for i=1:N
%         p_inc=0;
%         for j=1:N
%             if ContStruc(i,j)~=0
%                 K(m_inc+1:m_inc+m(i),p_inc+1:p_inc+p(j))=100*randn(m(i),p(j));
%             end
%             p_inc=p_inc+p(j);
%         end
%         m_inc=m_inc+m(i);
%     end
%     eF=round(eig(A+Btot*K*Ctot),rounding_n);
%     C=intersect(Difm,eF);
%     Difm=C;
%     nelD=length(Difm);
% end
% end
% function plot_eig_CT(A, alpha, radius)
%  if nargin < 2
%      alpha = 0;
%  end
% 
%  if nargin < 3
%      radius=-1;
%  end
% 
% eigenvaluesCT = eig(A);
% 
%  if radius == -1
%     figure()
%     grid on
%     hold on
%     plot(real(eigenvaluesCT), imag(eigenvaluesCT), '*')
%     xline(-alpha)
%     hold off
%     title('Eigenvalues')
%     xlabel('Re')
%     ylabel('Im')
%  else
%     figure()
%     grid on 
%     hold on
%     plot(real(eigenvaluesCT), imag(eigenvaluesCT), '*')
%     unit_c = circle(-1*alpha,0,radius);
%     hold off
%     title('Eigenvalues')
%     xlabel('Re')
%     ylabel('Im')
% 
%  end
% end
% function plot_eig_DT(F, alpha, radius)
%  if nargin < 2
%      alpha = 0;
%  end
% 
%  if nargin < 3
%      radius = -1;
%  end
% 
% eigenvaluesDT = eig(F);
% 
%  if radius == -1
%     figure()
%     grid on
%     hold on
%     plot(real(eigenvaluesDT), imag(eigenvaluesDT), '*')
%     disk = circle(0,0,1);
%     hold off
%     title('Eigenvalues')
%     xlabel('Re')
%     ylabel('Im')
%  else
%     figure()
%     grid on 
%     hold on
%     plot(real(eigenvaluesDT), imag(eigenvaluesDT), '*')
%     disk = circle(-1*alpha,0,radius);
%     hold off
%     title('Eigenvalues')
%     xlabel('Re')
%     ylabel('Im')
% 
%  end
% end
% function circle = circle(x,y,r)
% hold on
% th = 0:pi/50:2*pi;
% xunit = r * cos(th) + x;
% yunit = r * sin(th) + y;
% circle = plot(xunit, yunit);
% end
% function plot_freemotion(N, v, state, T, x_De_C, x_De_DT, x0, h, Tfinal)
% figure()
% for i=1:N
%     subplot(N,2,2*(i-1)+1)
%     hold on
%     grid on
%     title([state,num2str(i),'}'])
%     plot(T,[x_De_C((i)*4-(4-v),:)],'m')
%     axis([0 T(end) min(x0)-4 max(x0)+4])
% 
%     subplot(N,2,2*i)
%     hold on
%     grid on
%     title([state,num2str(i),'}'])
%     plot([0:h:Tfinal],[x0(i*4-(4-v)),x_De_DT((i)*4-(4-v),:)],'m.-')
%     axis([0 T(end) min(x0)-4 max(x0)+4])
% end
% legend('Decentralized')
% end