function plot_eig_DT(F, alpha, radius)
% Plots the eigenvalues of matrix A in the discrete time case.
% Inputs:
% - A: System Matrix
% - Alpha: Center of the disk region. If there is no input for alpha, it is set to 0. If there is no input for radius 
%   alpha is the limit value for the spectral abscissa instead of the center of a disk
% - Radius: Radius of of the disk region.

 if nargin < 2
     alpha = 0;
 end

 if nargin < 3
     radius = -1;
 end

eigenvaluesDT = eig(F);

 if radius == -1
    figure()
    grid on
    hold on
    plot(real(eigenvaluesDT), imag(eigenvaluesDT), '*')
    disk = circle(0,0,1);
    hold off
    title('Eigenvalues')
    xlabel('Re')
    ylabel('Im')
 else
    figure()
    grid on 
    hold on
    plot(real(eigenvaluesDT), imag(eigenvaluesDT), '*')
    disk = circle(-1*alpha,0,radius);
    hold off
    title('Eigenvalues')
    xlabel('Re')
    ylabel('Im')

 end
end