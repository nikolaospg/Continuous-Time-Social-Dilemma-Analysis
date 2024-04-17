%We pass the current x and the system parameters and get the corresponding x_dot.
%It can be used to model both symmetric and non symmetric systems (with the dynamics prescribed on the paper)
%The _vec parameters are 8x1 matrices specifying the parameters of the system
function x_dot=system_dynamics3D(t, x, a_vec, b_vec, c_vec)

    %"Extracting" easily handable scalar values from the vector
    x1=x(1);
    x2=x(2);
    x3=x(3);
    
    a000=a_vec(1);a001=a_vec(2);a010=a_vec(3);a011=a_vec(4);a100=a_vec(5);a101=a_vec(6);a110=a_vec(7);a111=a_vec(8);
    b000=b_vec(1);b001=b_vec(2);b010=b_vec(3);b011=b_vec(4);b100=b_vec(5);b101=b_vec(6);b110=b_vec(7);b111=b_vec(8);
    c000=c_vec(1);c001=c_vec(2);c010=c_vec(3);c011=c_vec(4);c100=c_vec(5);c101=c_vec(6);c110=c_vec(7);c111=c_vec(8);
    %Finished getting the scalar values
    
    %Calculating the derivatives :
    x_dot1= (a000*(1-x1)*(1-x2)*(1-x3) + a001*(1-x1)*(1-x2)*x3 + a010*(1-x1)*x2*(1-x3) + a011*(1-x1)*x2*x3 ...
          + a100*x1*(1-x2)*(1-x3) + a101*x1*(1-x2)*x3 + a110*x1*x2*(1-x3) + a111*x1*x2*x3)*x1*(1-x1);
      
    x_dot2= (b000*(1-x1)*(1-x2)*(1-x3) + b001*(1-x1)*(1-x2)*x3 + b010*(1-x1)*x2*(1-x3) + b011*(1-x1)*x2*x3 ...
          + b100*x1*(1-x2)*(1-x3) + b101*x1*(1-x2)*x3 + b110*x1*x2*(1-x3) + b111*x1*x2*x3)*x2*(1-x2);
      
    x_dot3= (c000*(1-x1)*(1-x2)*(1-x3) + c001*(1-x1)*(1-x2)*x3 + c010*(1-x1)*x2*(1-x3) + c011*(1-x1)*x2*x3 ...
          + c100*x1*(1-x2)*(1-x3) + c101*x1*(1-x2)*x3 + c110*x1*x2*(1-x3) + c111*x1*x2*x3)*x3*(1-x3);   
      
      
    x_dot=[x_dot1;x_dot2;x_dot3];
end