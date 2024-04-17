clc
clear
close all

%% Params etc
% System parameters:
a00=1;
a01=-1;
a10=-1;
a11=-1;

b00=1;
b01=-1;
b10=-1;
b11=-1;

%Simulation and phase portait parameters
x_step=0.01;
y_step=0.01;
make_phase_portait=1;             %Flag on whether to make the portaits or not
zero_threshold=0.008;             %This is a parameter in the create_phase_plots func

time_step=0.01;
time_lim=10;
x_init=[0.5; 0.5];                  %The initial state 
make_simulation_plots=1;
%%Finished with the parameter initialisation

x_dot=@(t,state)system_dynamics(t, state, a00, a01, a10, a11, b00, b01, b10, b11) ;     %Function handle for the dynamical system 

%%  Creating the x_range and y_range values and calling the create_phase_plots function
x_range=0:x_step:1;
y_range=0:y_step:1;
if(make_phase_portait==1)
    [x1_dot_grid, x2_dot_grid, x1_dot_grid_norm, x2_dot_grid_norm]=create_phase_plots(x_dot, x_range, y_range, zero_threshold);
end
%% Finished with the mesh grid 


%% Now simulating 
time_vec=0:time_step:time_lim;
[t, state]=ode45(x_dot, time_vec, x_init);

if(make_simulation_plots==1)
    simulate_plots(x_init, t, state)
end
%% Finished simulating

%% Using the model I got from equilibria_study, I see whether the curve that is predicted actually has zero derivatives (i.e. the points predicted are actually equilibria)
norms=sqrt(x1_dot_grid.^2 + x2_dot_grid.^2);
x2=@(x1)(0.49952 - 0.47612*x1 -0.78643 * x1^2 + 0.76292 * x1^3 - 2.5558 * x1^4);

x_dots_norms=[];
for x1=0:0.01:0.5
    current_x1=x1;
    current_x2=x2(x1);
    x_dots_norms=[x_dots_norms; norm(x_dot(0,[current_x1; current_x2]))];
end

mean_absolute_xdot_vals=mean(abs(x_dots_norms))
mean_squared_xdot_vals=mean(x_dots_norms.^2);
root_mean_squared_xdot_vals=sqrt(mean(x_dots_norms.^2));

fprintf("The predicted curve (x1,x2) have xdot_norms with mean absolute value=%f mean squared value=%f root mean square value=%f\n",mean_absolute_xdot_vals, mean_squared_xdot_vals, root_mean_squared_xdot_vals)