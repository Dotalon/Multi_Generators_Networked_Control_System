function plot_freemotion(N, v, state, T, x_De_C, x_De_DT, x0, h, Tfinal)
% Plots the freemotions of the Decentralized Structure
% Inputs:
% - N: number of subsystems
% - v: number of the considered state
% - state: Name of the state (it needs to end with '_{' or '_{additional_subscript' )
% - x_De_C: Free motion (Continuous case)
% - x_De_DT: Free motion (Discrete case)
% - x0: starting conditions
% - h: sampling time
% - Tfinal: final time of the simulation 

figure()
for i=1:N
    subplot(N,2,2*(i-1)+1)
    hold on
    grid on
    title([state,num2str(i),'}'])
    plot(T,[x_De_C((i)*4-(4-v),:)],'m')
    axis([0 T(end) min(x0)-4 max(x0)+4])
    
    subplot(N,2,2*i)
    hold on
    grid on
    title([state,num2str(i),'}'])
    plot([0:h:Tfinal],[x0(i*4-(4-v)),x_De_DT((i)*4-(4-v),:)],'m.-')
    axis([0 T(end) min(x0)-4 max(x0)+4])
end
legend('Decentralized')
end