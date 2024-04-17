%After the simulation is done I use this function to create some plots
function simulate_plots3D(x_dot, x_init, t, state)
  
    %First the 3D plot
    title=sprintf("3D representation, x0=(%.2f, %.2f,%.2f)", x_init(1), x_init(2), x_init(3));
    figure("name",title)
    plot3(state(:,1), state(:,2), state(:,3), "LineWidth", 3)
    hold on
    scatter3(x_init(1), x_init(2), x_init(3), 80, 'r', 'o')
    xlim([-0.01 1.01])
    ylim([-0.01 1.01])
    zlim([-0.01 1.01])
    [X_mesh,Y_mesh,Z_mesh] = meshgrid(0:0.05:1,0:0.05:1,0:0.05:1);
    scatter3(state(end,1), state(end,2), state(end,3), 80, 'r', 'x')
    scatter3(X_mesh(:), Y_mesh(:), Z_mesh(:), 2, 'k','filled')

    %Next the two 1D plots
    title=sprintf("x1 values through time, x0=(%.2f, %.2f, %.2f)", x_init(1), x_init(2), x_init(3));
    figure("name", title);
    plot(t, state(:, 1))
    ylabel("x1")
    xlabel("t: sec")

    title=sprintf("x2 values through time, x0=(%.2f, %.2f, %.2f)", x_init(1), x_init(2), x_init(3));
    figure("name", title);
    plot(t, state(:, 2))
    ylabel("x2")
    xlabel("t: sec")
  
    title=sprintf("x3 values through time, x0=(%.2f, %.2f, %.2f)", x_init(1), x_init(2), x_init(3));
    figure("name", title);
    plot(t, state(:, 3))
    ylabel("x3")
    xlabel("t: sec")
  
  
  
end
