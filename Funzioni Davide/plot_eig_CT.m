function plot_eig_CT(A, alpha, radius)
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
    unit_c = circle(-alpha,0,radius);
    hold off
    title('Eigenvalues')
    xlabel('Re')
    ylabel('Im')

 end
end