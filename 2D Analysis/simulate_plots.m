%After the simulation is done I use this function to create some plots
function simulate_plots(x_init, t, state)
  
    %First the 2D plot
    title=sprintf("2D representation, x0=(%.2f, %.2f)", x_init(1), x_init(2));
    figure("name",title)
    plot(state(:,1), state(:,2), "LineWidth", 3)
    hold on
    scatter(x_init(1), x_init(2), 80, 'r', 'o')
    xlim([-0.01 1.01])
    ylim([-0.01 1.01])
    [X_mesh,Y_mesh] = meshgrid(0:0.05:1,0:0.05:1);
    scatter(state(end,1), state(end,2), 80, 'r', 'x')
    scatter(X_mesh(:), Y_mesh(:), 5, 'k','filled')

    %Next the two 1D plots
    title=sprintf("x1 values through time, x0=(%.2f, %.2f)", x_init(1), x_init(2));
    figure("name", title);
    plot(t, state(:, 1))
    ylabel("x1")
    xlabel("t: sec")

    title=sprintf("x2 values through time, x0=(%.2f, %.2f)", x_init(1), x_init(2));
    figure("name", title);
    plot(t, state(:, 2))
    ylabel("x2")
    xlabel("t: sec")
  
  
  
end
