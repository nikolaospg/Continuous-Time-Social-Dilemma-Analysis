clc
clear
close all

%In this script we use the compute_equilibria function to numerically get the equilibria of our system.

%% Params etc
load_equilibria_flag=1;
save_file_name="equilibria.mat";
classification_analysis=1;         %Flag on whether to apply the classification analysis on the cluster labels or not 

a00=1;
a01=-1;
a10=-1;
a11=-1;

b00=1;
b01=-1;
b10=-1;
b11=-1;


x_step=0.01;
y_step=0.01;
time_step=0.01;
time_lim=15;
%%Finished with the parameters

%Creating a matrix with a lot various initial states
x_range=0:x_step:1;
y_range=0:y_step:1;
init_states=zeros(length(x_range) * length(y_range), 2);
count=1;
for current_x_index=1:length(x_range)
    current_x=x_range(current_x_index);
    for current_y_index=1:length(y_range)
        current_y=y_range(current_y_index);
        init_states(count,:)=[current_x, current_y];
        count=count+1;
    end
end
%Finished creating the initial states 

%Calling the compute_equilibria function
x_dot=@(t,state)system_dynamics(t, state, a00, a01, a10, a11, b00, b01, b10, b11) ;     %Function handle for the dynamical system 
x_dot_norm=@(state)norm(x_dot(0,state));
x_range=0:x_step:1;
y_range=0:y_step:1;
time_vec=0:time_step:time_lim;

if(load_equilibria_flag==1)
    equilibria=load(save_file_name);
    equilibria=equilibria.equilibria;
else
    equilibria=compute_equilibria(x_dot, time_vec, init_states);          %TO BE CHANGED
    save(save_file_name, 'equilibria'); 
end
%Finished computing the equilibria



%In the following for loop I get rid of the equilibria that are on the 4
%edges. I also calculate some values for every equilibrium that I need
%later on.
norms=[];                               %The norm of each equilibrium
eq_coord_sums=[];                 %The sum of x1+x2 for each equilibrium
equilibria_reduced=[];            %The list of equilibria after getting rid of the bad ones 
accepted_index=[];                  %The indexes of the "accepted equilibria" on the original equilibria list
for i=1:length(equilibria)
    current_equilibrium=equilibria(i,:);

    if(norm(current_equilibrium-[1,1])<0.1)
        continue
    end
    
    if(norm(current_equilibrium-[0,0])<0.1)
        continue
    end
    
    if(norm(current_equilibrium-[1,0])<0.1)
        continue
    end
    
    if(norm(current_equilibrium-[0,1])<0.1)
        continue
    end
    
    if(current_equilibrium(1)>0.5)
        continue
    end
    
    if(current_equilibrium(2)>0.5)
        continue
    end
    
    current_eq_coord_sum=current_equilibrium(1)+current_equilibrium(2);
    current_norm=norm(current_equilibrium);
    
    equilibria_reduced=[equilibria_reduced;current_equilibrium];
    norms=[norms;current_norm];
    eq_coord_sums=[eq_coord_sums;current_eq_coord_sum];
end
%Finished finding the equilibria that are not on the 4 edges.

%Making some scatter plots to understand what is going on:
figure("Name", "The calculated equilibria")
scatter(equilibria(:,1), equilibria(:,2))


figure("Name", "The reduced equilibria")
scatter(equilibria_reduced(:,1), equilibria_reduced(:,2))
xlim([0,1])
ylim([0,1])
%From the scatter plot we see that there is a curve

%%Now trying to model the curve:

%First hypothesis: The curve is the arc of the circle with radius=0.5 and center=0
figure("Name", "Histogram of the norms")
hist(norms, 500)
%The first hypothesis is rejected, because looking at the norms figure the norm of the equilibria is not constant=0.5

%Second hypothesis: The equilibria are on the locus of the points where the sum of the coordinates is constant and equal to 0.5
figure("Name", "Histogram of the sum of the coordinates ")
hist(eq_coord_sums, 500)
%The hypothesis is rejected as the scatter plot does not seem to show constant values




%Regression analysis:
%Third hypothesis: y=a1+a2x+a3x^2 + ... a5x^4    where y=x2, x=x1 (for simplicity)
Y=zeros(length(norms),1);
X=zeros(length(norms),5);

%Creating the data matrices:
for i=1:length(Y)
    current_y=equilibria_reduced(i,2);
    current_x1=equilibria_reduced(i,1);
    current_X=[1, current_x1, current_x1^2, current_x1^3, current_x1^4];
    
    Y(i)=current_y;
    X(i,:)=current_X;
end


%Applying the Ordinary Least Squares Solution:
X_trans=X';
Q=inv(X_trans*X);
solution=Q*X_trans*Y;
preds=X*solution;
errors=equilibria_reduced(:,2)-preds;

%Regression results:
figure("Name", "Regression Residuals")
hist(errors, 500)

MSE=mean(errors.^2);
RMSE=sqrt(MSE);
fprintf("MSE=%f Rmse=%f\n", MSE, RMSE)

figure("Name", "Regression Results")
scatter(equilibria_reduced(:,1), equilibria_reduced(:,2), 0.1, 'r', 'x')
hold on
scatter(equilibria_reduced(:,1), preds, 0.1, 'k', 'o')
legend("Actual points", "Modeled Curve")
%Finished with the curve modeling


%Solving the classification problem:
if(classification_analysis)
    classify(init_states, 0)
    classify(init_states, 1)
end
%Solving with the classification




%---------------------------------FUNCTIONS---------------------------------------------------------%

%This function simulates the system with a variety of initial states (generated from the grid), and returns the final state. 
%This should be able to help us study the equilibria and the steady state behaviour.
function results_matrix=compute_equilibria(fun_handle, time_vec, init_states)
    
    init_states_dims=size(init_states);
    results_matrix=ones(init_states_dims(1),2);
    point_count=1;
    
    for s0_index=1:init_states_dims(1)
        current_s0=init_states(s0_index, :);
        [t, state]=ode45(fun_handle, time_vec, current_s0);
        final_state=state(end,:);
        results_matrix(point_count,:)=final_state;
        point_count=point_count+1;
        if(mod(s0_index, 100)==0)
            fprintf("point_count %d out of %d\n",point_count, init_states_dims(1))
        end
    end
  
end



%This function opens the labels we get by using the dbscan algorithm and
%applies classification analysis using Multinomial Logistic Regression.
%This way we try to find some way to be able to get information about the
%basins of attraction for the system
function B=classify(init_states, polynomial_flag)


    %loading the data
    labels=readtable("clusters.csv");
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
        X(:, [1,2])=init_states;
        X(:,3)=init_states(:,1).^2;                 %The x1 squared
        X(:,4)=init_states(:,2).^2;                 %The x2 squared
        X(:,5)=init_states(:,2).*init_states(:,1);
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
