function plot_eig_CT(A, alpha, radius)
% Plots the eigenvalues of matrix A in the continuous time case.
% Inputs:
% - A: System Matrix
% - Alpha: Center of the disk region. If there is no input for alpha, it is set to 0. If there is no input for radius 
%   alpha is the limit value for the spectral abscissa instead of the center of a disk
% - Radius: Radius of of the disk region.

 if nargin < 2
     alpha = 0;
 end

 if nargin < 3
     radius=-1;
 end

eigenvaluesCT = eig(A);

 if radius == -1
    figure()
    grid on
    hold on
    plot(real(eigenvaluesCT), imag(eigenvaluesCT), '*')
    xline(-alpha)
    hold off
    title('Eigenvalues')
    xlabel('Re')
    ylabel('Im')
 else
    figure()
    grid on 
    hold on
    plot(real(eigenvaluesCT), imag(eigenvaluesCT), '*')
    unit_c = circle(-1*alpha,0,radius);
    hold off
    title('Eigenvalues')
    xlabel('Re')
    ylabel('Im')

 end
end