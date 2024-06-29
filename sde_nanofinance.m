% /////////////////////////////////////////
% SDE simulation of the NF model
% author: Rolando Gonzales Martinez
% v. May 2024
% ////////////////////////////////////////
cd 'Z:\UpWork\UpWork2024\Muyu_2\Simulation'
close all; clear; clc;
nsimulations = 10000; % number of SDE simulations
% Specify the Excel file name
excelFileName = 'model_parameters_and_KLdistance.xlsx';
% ------------------------------------------------------------------------
figure_name = 'figure_5035.png';
startingRow = 2; % starting row
% Parameters of the model
    mu_e = 1.7;        % drift external
    mu_i = .9;        % drift internal
    sigma_e = .8;     % difusion external
    sigma_i = .9;     % difusion internal
    phi = 0.7894;      % weight of the drift process for external factors
    lambda = phi;      % weight of the difussion process for external factors 
% ------------------------------------------------------------------------    
% Defining the extended SDE: savings
    F_s = @(t,S)(phi*mu_e + (1 - phi)*mu_i)*S;
    G_s = @(t,S)(lambda*sigma_e + (1 - lambda)*sigma_i)*S;
    obj_S = sde(F_s,G_s);
% Defining the extended SDE: returns on savings
    F_r = @(t,S) phi*mu_e + (1 - phi)*mu_i - (lambda*sigma_e + (1 - lambda)*sigma_i)^2 /2;
    G_r = @(t,S) lambda*sigma_e + (1 - lambda)* sigma_i;
    obj_r = sde(F_r,G_r);
% ------------------------------------------------------------------------    
% Simulations
    nPeriods = 12;      % cycle of 12 months     
    T = 1:nPeriods+1;   
    dt = 1;             % 1 month
    simulations = nsimulations; %
           % Pre-allocating for speed
           S = zeros(nPeriods+1,simulations);
           r = zeros(nPeriods+1,simulations);
         ROS = zeros(nPeriods+1,simulations);
rng(666)           
    for s = 1:simulations
        S(:,s) = simulate(obj_S, nPeriods, 'DeltaTime', dt);
        r(:,s) = simulate(obj_r, nPeriods, 'DeltaTime', dt);
      ROS(:,s) = r(:,s)./S(:,s);
    end
% Figure
    S_mean   = mean(S,2);
    S_median = median(S,2);
    r_mean   = mean(r,2);
    r_median = median(r,2);
  ROS_mean   = mean(ROS,2);
  ROS_median = median(ROS,2);
% ------------------------------------------------------------------------
% Figure
set(gcf,'Position',[100 100 900 400])    
subplot(1,4,1)
     hold on
            plot(T, S)
            plot(T, S_median, 'k','LineWidth',2)
            xlabel('Time in the cycle (t)')
            ylabel('Savings accumulation (S_t)')
            box on
            grid on
            xlim([1 12])
            title('(i)','FontWeight','normal')
    hold off
 subplot(1,4,2)
    hold on
        plot(T, r)
        plot(T, r_median, 'k','LineWidth',2)
        xlabel('Time in the cycle (t)')
        ylabel('returns (r_t)')
        box on
        grid on
        xlim([1 12])
        title('(ii)','FontWeight','normal')
    hold off
% Returns on savings
subplot(1,4,3)
    ROS_simulated = reshape(ROS,[],1)*100; % Simulated values
    boxplot(ROS_simulated)
    title('(iii)','FontWeight','normal')
    xlabel('Simulated ROS')
    ylabel('percentage (%)')
    xticklabels({''}); % Change x-axis tick labels
    box on
    grid on
    % ylim([0, 140])
subplot(1,4,4)
load ROS_BDM
    boxplot(ROS_savix)
    title('(iv)','FontWeight','normal')
    xlabel('Observed ROS')
    ylabel('percentage (%)')
    xticklabels({''}); % Change x-axis tick labels
    box on
    grid on    
    ylim([0, 140])
% Save the figure
saveas(gcf,figure_name); % Replace 'your_figure_name' with the desired file name and extension    
% ------------------------------------------------------------------------
% KL distance
% Add a column with values 0 for ROS_savix and 1 for ROS_simulated
% I = [zeros(size(ROS_savix, 1), 1); ones(size(ROS_simulated, 1), 1)];
I = [false(size(ROS_savix, 1), 1); true(size(ROS_simulated, 1), 1)];
% Concatenate the new column to the existing data
R = [ROS_savix; ROS_simulated];
X = [R, I];
KL = relativeEntropy(X,I);
KL_distance = KL(1);
% ------------------------------------------------------------------------
% Store the results in Excel:
    % Create a matrix with parameter values as rows
    parametersMatrix = [mu_e, mu_i, sigma_e, sigma_i, phi, lambda, KL_distance];
    headers = {'mu_e', 'mu_i', 'sigma_e', 'sigma_i', 'phi', 'lambda', 'KL distance'};
    writecell(headers, excelFileName, 'Range', 'A1:G1');
    writematrix(parametersMatrix, excelFileName, 'Range', ['A', num2str(startingRow)]);    
disp('simulation concluded and results stored');
disp(['KL_distance: ', mat2str(KL_distance)]);    
% -----------------  e n d   o f   s c r i p t  ---------------------------
    
