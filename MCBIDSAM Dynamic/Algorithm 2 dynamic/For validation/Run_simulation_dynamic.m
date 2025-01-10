%% File for simulating model for any set of values of inputs
% Please ensure you have the following in the same directory as this file:
% 1) "simulate_model_dynamic" function file, 2) "BasisInterpretation"
% function file and 3) "solution_workspace" mat file.

% Step 1: declare the set of inputs (u) for which you desire to use the model
% to predict outputs. This can be for just a time instant or multiple
% instances e.g load('u.mat'), u=[2 5 5]'....
constr = true;
constr = false;
selection =  1550:1900;
data = dd(:,selection); u=uu(:,selection);

% Step 2: declare the rank (slnrank) of the model you desire to run from the set of
% hierarchically ranked models the algorithm has identified. e.g 
slnrank = 12;

% Step 3: declare initial values of output variables e.g
y_initial = data(:,1);

% Step 3: run simulation
y = simulate_model_dynamic2(u, y_initial, slnrank, constr);



for k=1:N
    figure
    plot(y(k,:),'Color',"#77AC30")
    hold on
    plot(data(k,:),'-.','Color',"#D95319")
    hold off
    legend('BIDSAM model','Data')
    title(['Training plot for variable y', num2str(k)])
    xlabel('time indices')
    ylabel(['y',num2str(k)])
end

er=data-y(:,1:end-1);
ersq=er.*er;
SSE_all = sum(ersq,'all')
SSE = sum(ersq,2);
MSE = SSE/size(data,2);
RMSE = sqrt(MSE)


load('holdup_true.mat')
holdsel = holdup_true(selection);
nci = [1 0 1 0];
ncio = [1 0 1 0];
for j=1:numel(selection)-1
    c_in(j) = nci*u(:,j);
    c_out(j) = ncio*y(1:N,j+1);
    c_acc(j) = nci*u(:,j)-ncio*data(1:N,j+1); %holdsel(1,j);
    c_bal(j) = c_in(j) - c_out(j); - c_acc(j);
end
c_bal_rel = 100*c_bal./c_in;

