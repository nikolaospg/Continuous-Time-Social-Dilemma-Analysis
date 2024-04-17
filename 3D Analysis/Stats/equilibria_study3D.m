clc
clear
close all

%In this script we use the compute_equilibria function to numerically get the equilibria of our system.

%% Params etc
load_equilibria_flag=1;
equilibria_file_name="equilibria2.mat";
states_file_name="states2.mat";
classification_analysis=0;         %Path on whether to apply the classification analysis on the cluster labels or not 

a_vec=[-1;-1;-1;1;1;1;-1;1];
b_vec=[1;-1;1;1;-1;1;-1;1];
c_vec=[-1;-1;-1;-1;-1;1;1;-1];
x_dot=@(t,state)system_dynamics3D(t, state, a_vec, b_vec, c_vec) ;     %Function handle for the dynamical system 


space_step=0.01;                %The step size for x,y,z
num_interior_points=300000;
time_step=0.01;
time_lim=200;
%%Finished with the parameters


%Calling the compute_equilibria function (if the flag says so, else we just
%load the already calculated ones)
x_range=0:space_step:1;
y_range=0:space_step:1;
time_vec=0:time_step:time_lim;

if(load_equilibria_flag==1)
    equilibria=load(equilibria_file_name);
    init_states=load(states_file_name);
    equilibria=equilibria.equilibria;
    init_states=init_states.init_states;
else
    init_states=get_init_states(space_step, num_interior_points);
    equilibria=compute_equilibria(x_dot, time_vec, init_states);
    save(equilibria_file_name, 'equilibria'); 
    save(states_file_name, 'init_states'); 
end
%Finished getting the equilibria 




%Doing classification analysis - AFTER GETTING THE LABELS
if(classification_analysis)
    classify(init_states, "my_labels.csv", 0)
    classify(init_states, "my_labels.csv", 1)
end
%Finished with the classification analysis

%%Applying PC analysis
labels=readtable("my_labels.csv");
labels=labels(2:end, 2);
y=table2array(labels)+2;              %Adding 2 so that we have only positive labels   
max_labels=max(y);

%Each iteration of the following loop corresponds to one cluster
for label=2:max_labels
    indices=find(y==label);
    current_init_states=init_states(indices,:);
    current_equilibria=equilibria(indices,:);
    fprintf("For label=%d:\n\n",label-2);
    [~,~,latent,~,explained] = pca(current_equilibria);
    fprintf("Perc Variance Explained:\n")
    disp(explained)
    fprintf("Variances:\n")
    disp(latent)
end
fprintf("***CONCLUSIONS FROM THE PCA***\n")
fprintf("As we see, in all of the clusters, the variance is less than 1e-6 ...\nAll except one, as in the label==1 the variance is clearly larger than the others\n")
fprintf("This leads us to the conclusion that all of the clusters, except the one with label==1, must represent points\n")
fprintf("The ones with label==-1 are classified as noise by OPTICS. This means that the points are either equilibria,\nor points that due to the computational nature of the network could not reach their steady state.\n")
fprintf("The cluster with label==1 should be either a surface or a curve. Looking at the percentage of variances explained by each variable\non the PC analysis,we see that the first variable on the PCA analysis is clearly larger than the other two\n")
fprintf("This leads us to the final conclusion, that according to the PC analysis the cluster with label==1 might be a curve!\n")
fprintf("We therefore make the hypothesis that this cluster respesents an attractor that is a curve, and all the others are point equilibria.\n")
%%Finished with the PCA


%% Modeling the curve:
final_indices=find(y==3);
final_init_states=init_states(final_indices,:);
final_equilibria=equilibria(final_indices,:);
model_curve(final_equilibria)
%% Finished with the modeling


%%Final results:
figure("Name", "The calculated equilibria")
scatter3(equilibria(:,1), equilibria(:,2), equilibria(:,3))
xlabel("x1")
ylabel("x2")
zlabel("x3")

%Plotting the equilibria:
figure("Name", "The equilibria and the curves")
scatter3(final_equilibria(:,1), final_equilibria(:,2), final_equilibria(:,3))
xlabel("x1")
ylabel("x2")
zlabel("x3")
hold on
%Plotting the r1 curve
x2_range_plots=[0:0.01:1]';
R1=[ones(length(x2_range_plots),1), x2_range_plots, 0.5*ones(length(x2_range_plots),1)];
scatter3(R1(:,1), R1(:,2), R1(:,3), 'r')

%Plotting the r2 curve
final_equilibria_dims=size(final_equilibria);
count=1;
for index=1:final_equilibria_dims(1)
   current_equilibrium=final_equilibria(index,:);
   if(current_equilibrium(1)<0.99)
       R2(count,1)=current_equilibrium(1);
       R2(count,2)=1;
       R2(count,3)=0.0375 + 0.4626*current_equilibrium(1);
       count=count+1;
   end
end
scatter3(R2(:,1), R2(:,2), R2(:,3), 'k')
zlim([0,1])
legend("Observations", "First Curve", "Second Curve")
%FInished plotting the equilibria


% %---------------------------------FUNCTIONS---------------------------------------------------------%

%This function simulates the system with a variety of initial states (generated from the grid), and returns the final state. 
%This should be able to help us study the equilibria and the steady state behaviour.
function results_matrix=compute_equilibria(fun_handle, time_vec, init_states)
    
    init_states_dims=size(init_states);
    results_matrix=ones(init_states_dims(1),3);
    total_points=init_states_dims(1);
    
    parfor s0_index=1:total_points
        current_s0=init_states(s0_index, :);
        [~, state]=ode45(fun_handle, time_vec, current_s0);
        final_state=state(end,:);
        results_matrix(s0_index,:)=final_state;
        if(mod(s0_index, 100)==0)
            fprintf("point_count %d out of %d\n",s0_index, total_points)
        end
    end
  
end

%This function opens the labels we get by using the dbscan algorithm and
%applies classification analysis using Multinomial Logistic Regression.
%This way we try to find some way to be able to get information about the
%basins of attraction for the 
function B=classify(init_states, cluster_file_name, polynomial_flag)


    %loading the data
    labels=readtable(cluster_file_name);
    labels=labels(2:end, 2);
    y=table2array(labels)+2;              %Adding 2 so that we have only positive labels
    %Finished loading the data
    
    %Creating the X
    if(polynomial_flag==0)
        fprintf("Linear Classification Model\n")
        X=init_states;
    else
        fprintf("Quadratic Classification Model\n")
        X=zeros(length(init_states(:,1)), 5);
        X(:, [1,2,3])=init_states;
        X(:,4)=init_states(:,1).^2;                             %The x1 squared
        X(:,5)=init_states(:,2).^2;                             %The x2 squared
        X(:,6)=init_states(:,3).^2;                             %The x3 squared
        X(:,7)=init_states(:,2).*init_states(:,1);              %The x1x2
        X(:,8)=init_states(:,3).*init_states(:,1);              %The x2x3
        X(:,9)=init_states(:,2).*init_states(:,3);              %The x1x3
    end
    %Finished creating the X
    
    %Train test split
    partition_object=cvpartition(y, "Holdout", 0.3);            %Applying stratified sampling, Holdout with 0.7/0.3
    X_train=X(partition_object.training, :);
    X_test=X(partition_object.test, :);
    y_train=y(partition_object.training,:);
    y_test=y(partition_object.test,:);
    %Finished with the split
    
    %Training 
    fprintf("Training MLG model\n")
    B=mnrfit(X_train,y_train);
    fprintf("Finished training MLG model\n")
    train_pred_proba=mnrval(B, X_train);           %Probability estimates 
    [~,train_preds]=max(train_pred_proba');
    train_preds=train_preds(:);
    %Finished training
    
    %Making predictions on the test set
    test_pred_proba=mnrval(B, X_test);           %Probability estimates 
    [~,test_preds]=max(test_pred_proba');
    test_preds=test_preds(:);
    %Finished with the predictions on the test set
   
    %Printing confusion matrices
    train_acc=length(find(y_train==train_preds))/length(y_train);
    if(polynomial_flag==0)
        fig_title=sprintf("Train Conf Matrix, Linear, Acc=%f", train_acc);
    else
        fig_title=sprintf("Train Conf Matrix, Quadratic, Acc=%f", train_acc);
    end
    figure("Name", fig_title)
    conf_train=confusionchart(y_train-2, train_preds-2);
    
    test_acc=length(find(y_test==test_preds))/length(y_test);
    if(polynomial_flag==0)
        fig_title=sprintf("Test Conf Matrix, Linear, Acc=%f", test_acc);
    else
        fig_title=sprintf("Test Conf Matrix, Quadratic, Acc=%f", test_acc);
    end
    figure("Name", fig_title)
    conf_test=confusionchart(y_test-2, test_preds-2);
    %Finished printing the confusion matrices

end

%We pass the step size in space (x_step, y_step, z_step), and the number of
%interior points and get a bunch of initial states 
function init_states=get_init_states(step_size, num_interior_points)

    %Creating a matrix with a lot various initial states
    x_range=0:step_size:1;
    y_range=0:step_size:1;

    %First getting the data from the borders
    var_couples=zeros(length(x_range) * length(y_range), 2);
    count=1;
    for current_x_index=1:length(x_range)
        current_x=x_range(current_x_index);
        for current_y_index=1:length(y_range)
            current_y=y_range(current_y_index);
            var_couples(count,:)=[current_x, current_y];
            count=count+1;
        end
    end
    X_border_1=[ones(length(x_range) * length(y_range), 1), var_couples];
    X_border_2=[zeros(length(x_range) * length(y_range), 1), var_couples];
    Z_border_1=[var_couples, ones(length(x_range) * length(y_range), 1)];
    Z_border_2=[var_couples, zeros(length(x_range) * length(y_range), 1)];
    Y_border_1=[var_couples(:,1), zeros(length(x_range) * length(y_range), 1), var_couples(:,2)];
    Y_border_2=[var_couples(:,1), ones(length(x_range) * length(y_range), 1), var_couples(:,2)];
    %Finished with the data from the borders

    %Now getting data from the interior, using a uniform distribution
    interiors=rand([num_interior_points,3]);
    %Finished with the interior points

    %Merging all of the above in a final matrix
    init_states=[X_border_1; X_border_2; Y_border_1; Y_border_2; Z_border_1; Z_border_2; interiors];


end

%The code to be executed when we try to model the equilibria
function model_curve(equilibria)

    fprintf("\n***REGRESSION ANALYSIS***\n")
    X1=equilibria(:,1);
    X2=equilibria(:,2);
    X3=equilibria(:,3);
    
    %Creating all the 2D scatters
    figure("Name", "x2 vs x1")
    scatter(X1,X2)
    ylabel("x2")
    xlabel("x1")

    figure("Name", "x3 vs x1")
    scatter(X1,X3)
    ylabel("x3")
    xlabel("x1")

    figure("Name", "x3 vs x2")
    scatter(X2,X3)
    ylabel("x3")
    xlabel("x2")
    
    fprintf("Looking at the scatter plots, it seems that I could use a r(x1)=[x1, f2(x1), f3(x1)]\n")
    %Finished creating the 2D scatters 
    
    %Creating the augmented dataset, getting the limits
    max_X1=max(X1);
    min_X1=min(X1);
    limits=[min_X1, max_X1];
    fprintf("The limits of the parameter that models the curve, x1, are:")
    disp(limits)
    X1_augmented=[ ones(length(X1), 1), X1];
    %Finished with the augmented dataset.
    
    
    %%Beginning with x3=fun(x1)
    Q=X1_augmented'*X1_augmented;
    W_first=inv(Q)*X1_augmented'*X3;
    fprintf("The parameters of the first model, x3=linear(x1) are:")
    disp(W_first');
    
    %Now making a scatter plot with the observations to see if the
    %regression is succesful
    x3_pred=X1_augmented*W_first;
    figure("Name","Results of x3 vs x1 regression")
    xlabel("x1")
    ylabel("x3")
    scatter(X1,X3)
    hold on
    scatter(X1, x3_pred, 'r')
    legend("Observations", "Predictions")
    %%Finished with x3=fun(x1)
    
    
    %Now doing x2=fun(x1)
    %I will need a high degree polynomial. I try a 10 degree one
    degree=20;
    count=1;
    X1_augmented_pol=ones(length(X1),1);
    for current_degree=1:degree
        X1_augmented_pol=[X1_augmented_pol, X1.^current_degree];
    end
    Q_decond=X1_augmented_pol'*X1_augmented_pol;
    W_second=inv(Q_decond)*X1_augmented_pol'*X2;
    fprintf("The parameters of the second model, x2=pol(x1) deg=5 are:")
    disp(W_second');
    
    %Making the scatter plot to see if the regression is succesful
    x2_pred=X1_augmented_pol*W_second;
    figure("Name","Results of x2 vs x1 polynomial regression")
    xlabel("x1")
    ylabel("x2")
    scatter(X1,X2)
    hold on
    scatter(X1, x2_pred, 'r')
    legend("Observations", "Predictions")
    fprintf("As we see, the polynomial regression of x2=fun(x1) is unsuccesful... We are forced to use 2 different cases\n")
    

end

