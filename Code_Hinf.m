% Centralized Hinf
[K_C_CT,rho_C_CT,feas_C_CT]=LMI_CT_DeDicont_Hinf(A,Bdec,Cdec,N,ContStrucC);
sys_ss_CT=ss(A+B*K_C_CT, B, C, []);
sys_CT=tf(sys_ss_CT);
h_C_CT=norm(sys_CT,Inf)

[K_C_DT,rho_C_DT,feas_C_DT]=LMI_DT_DeDicont_Hinf(F,Gdec,Hdec,N,ContStrucC);
sys_ss_DT=ss(F+G*K_C_DT, G, H, []);
sys_DT=tf(sys_ss_DT);
h_C_DT=hinfnorm(sys_DT)

% Decentralized Hinf
[K_De_CT,rho_De_CT,feas_De_CT]=LMI_CT_DeDicont_Hinf(A,Bdec,Cdec,N,ContStrucDe);
sys_ss_CT=ss(A+B*K_De_CT, B, C, []);
sys_CT=tf(sys_ss_CT);
h_De_CT=norm(sys_CT,Inf)

[K_De_DT,rho_De_DT,feas_De_DT]=LMI_DT_DeDicont_Hinf(F,Gdec,Hdec,N,ContStrucDe);
sys_ss_DT=ss(F+G*K_De_DT, G, H, []);
sys_DT=tf(sys_ss_DT);
h_De_DT=hinfnorm(sys_DT)

% Distributed Hinf
[K_Di_CT,rho_Di_CT,feas_Di_CT]=LMI_CT_DeDicont_Hinf(A,Bdec,Cdec,N,ContStrucDi);
sys_ss_CT=ss(A+B*K_Di_CT, B, C, []);
sys_CT=tf(sys_ss_CT);
h_Di_CT=norm(sys_CT,Inf)

[K_Di_DT,rho_Di_DT,feas_Di_DT]=LMI_DT_DeDicont_Hinf(F,Gdec,Hdec,N,ContStrucDi);
sys_ss_DT=ss(F+G*K_Di_DT, G, H, []);
sys_DT=tf(sys_ss_DT);
h_Di_DT=hinfnorm(sys_DT)

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
    x_Hinf_C(:,k)=expm((A+B*K_C_CT)*t)*x0;
    
    % equations of decentralized free movement
    x_Hinf_De(:,k)=expm((A+B*K_De_CT)*t)*x0;
    
    % equations of distributed free movement
    x_Hinf_Di(:,k)=expm((A+B*K_Di_CT)*t)*x0;
end

%% Plot of the movements for every control law used (simple stability, alpha-stab,ecc) 

figure
for i=1:N
    subplot(N,3,1+(3*(i-1)))
    hold on
    grid on
    title(['SpeedC_{',num2str(i),'}'])
    plot(T,[x_Hinf_C((i)*4-2,:)],'k')
    legend('H_{inf} CENTR CT')

    subplot(N,3,2+(3*(i-1)))
    hold on
    grid on
    title(['SpeedDe_{',num2str(i),'}'])
    plot(T,[x_Hinf_De((i)*4-2,:)],'k')
    legend('H_{inf} DECENTR CT')

    subplot(N,3,3+(3*(i-1)))
    hold on
    grid on
    title(['SpeedDi_{',num2str(i),'}'])
    plot(T,[x_Hinf_Di((i)*4-2,:)],'k')
    legend('H_{inf} DISTR CT')
   
end