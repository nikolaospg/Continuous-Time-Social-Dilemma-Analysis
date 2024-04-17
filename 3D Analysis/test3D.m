clc
clear
close all

%% Params etc
% System parameters:
a_vec=[-1;-1;-1;1;1;1;-1;1];
b_vec=[1;-1;1;1;-1;1;-1;1];
c_vec=[-1;-1;-1;-1;-1;1;1;-1];
x_dot=@(t,state)system_dynamics3D(t, state, a_vec, b_vec, c_vec) ;     %Function handle for the dynamical system 

%Simulation and phase portait parameters
x_step=0.05;
y_step=0.05;
z_step=0.05;
make_phase_portait=1;             %Flag on whether to make the portaits or not
zero_threshold=0.008;             %This is a parameter in the create_phase_plots func

time_step=0.01;
time_lim=10;
x_init=[0.5; 0.5; 0.5];                  %The initial state 
make_simulation_plots=1;
%%Finished with the parameter initialisation


%%  Creating the x_range and y_range values and calling the create_phase_plots function
x_range=0:x_step:1;
y_range=0:y_step:1;
z_range=0:z_step:1;
if(make_phase_portait==1)
    [x1_dot_grid, x2_dot_grid, x3_dot_grid, x1_dot_grid_norm, x2_dot_grid_norm, x3_dot_grid_norm]=create_phase_plots3D(x_dot, x_range, y_range, z_range, zero_threshold);
end
%% Finished with the mesh grid 



%% Now simulating 
time_vec=0:time_step:time_lim;
[t, state]=ode45(x_dot, time_vec, x_init);

if(make_simulation_plots==1)
    simulate_plots3D(x_dot, x_init, t, state)
end
