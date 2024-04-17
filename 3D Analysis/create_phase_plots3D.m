%Function used to create the phase portait for the dyn system. It creates both a normalised (i.e. every vector of the plot is illustrated with the same length) and an unnormalised plot
%dynamics_handle-> The fun handle that gives me the x_dot 
%x_range, y_range, z_range-> Specifies the ranges for the 3 vars, in order to create the mesh 
%zero_threshold-> In the normalised plot It makes sense to set a threshold to judge on whether the norm of the x_dot vector should be regarded as zero or not.
% ----
%The function returns the values for the x1_dot and x2_dot (both normalised and unnormalised) on the grid points. Normalised means that the |x_dot|=1
function [U,Y, T, U_normalised, Y_normalised, T_normalised]=create_phase_plots3D(dynamics_handle, x_range, y_range, z_range, zero_threshold)
    [X_mesh,Y_mesh,Z_mesh] = meshgrid(x_range,y_range,z_range);
    U=zeros(size(X_mesh));                                  %The values of the x1_dot on the mesh
    Y=zeros(size(X_mesh));                                  %The values of the x2_dot on the mesh
    T=zeros(size(X_mesh));                                  %The values of the x3_dot on the mesh
    
    
    U_normalised=zeros(size(X_mesh));               %The values of the x1_dot on the mesh (normalised)
    Y_normalised=zeros(size(X_mesh));                 %The values of the x2_dot on the mesh (normalised)
    T_normalised=zeros(size(X_mesh));                                  %The values of the x3_dot on the mesh
    
    %In the following loop I calculate the x_dot values on the points specified by the mesh grid 
    for current_row_index=1:length(x_range)
        current_row=x_range(current_row_index);
        for current_col_index=1:length(y_range)
            current_col=y_range(current_col_index) ;
            for current_height_index=1:length(z_range)
                current_height=z_range(current_height_index) ;
                
                current_x_dot=dynamics_handle(1, [current_row; current_col;current_height] );
                U(current_row_index, current_col_index, current_height_index)=current_x_dot(1);
                Y(current_row_index, current_col_index, current_height_index)=current_x_dot(2);
                T(current_row_index, current_col_index, current_height_index)=current_x_dot(3);
                current_x_dot_norm=norm(current_x_dot);
                if(current_x_dot_norm<zero_threshold)                                      %Using the threshold so that the Visualisation actually makes sense in the normalised portait
                    U_normalised(current_row_index, current_col_index, current_height_index)=0;
                    Y_normalised(current_row_index, current_col_index, current_height_index)=0;
                    T_normalised(current_row_index, current_col_index, current_height_index)=0;
                else
                    U_normalised(current_row_index, current_col_index, current_height_index)=current_x_dot(1)/current_x_dot_norm;
                    Y_normalised(current_row_index, current_col_index, current_height_index)=current_x_dot(2)/current_x_dot_norm;
                    T_normalised(current_row_index, current_col_index, current_height_index)=current_x_dot(3)/current_x_dot_norm;
                end
            end
            
        end
    end
    
    %Creating the figures:
    figure("name", "Phase Portait")
    quiver3(X_mesh, Y_mesh, Z_mesh, U, Y, T)
    xlim([-0.01 1.01])
    ylim([-0.01 1.01])
    zlim([-0.01 1.01])
    
    
    figure("name", "Normalised Phase Portait")
    quiver3(X_mesh, Y_mesh,Z_mesh, U_normalised, Y_normalised, T_normalised, 1)
    
    xlim([-0.01 1.01])
    ylim([-0.01 1.01])
    zlim([-0.01 1.01])
  
end
