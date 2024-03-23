function plot_eig_DT(F, alpha, radius)
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
    disk = circle(-alpha,0,radius);
    hold off
    title('Eigenvalues')
    xlabel('Re')
    ylabel('Im')

 end
end
function circle = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
circle = plot(xunit, yunit);
end