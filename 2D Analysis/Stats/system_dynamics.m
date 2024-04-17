%We pass the current x and the system parameters and get the corresponding x_dot.
%It can be used to model both symmetric and non symmetric systems (with the dynamics prescribed on the paper)
function x_dot=system_dynamics(t, x, a00, a01, a10, a11, b00, b01, b10, b11)
    x1=x(1);
    x2=x(2);
    
    x_dot1=(a11*x1*x2 + a10*x1*(1-x2) + a01*(1-x1)*x2 + a00*(1-x1)*(1-x2))*x1*(1-x1);
    x_dot2=(b11*x1*x2 + b10*x2*(1-x1) + b01*(1-x2)*x1 + b00*(1-x1)*(1-x2))*x2*(1-x2);
    
    x_dot=[x_dot1;x_dot2];
end