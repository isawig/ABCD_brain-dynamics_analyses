%% PLSR model to explore relation between LEiDA and FCH
% Part A2, to run after with code_A1_plsr_leida6_10comp.m

% Calculates cross-validated MSE and find the number of PLS components that 
% minimizes it.

% Jetro J. Tuulari, jetro.tuulari@utu.fi

% Adapted by:
% Isabella L.C. Mariani Wigley, 09 / 2025; ilmawi@utu.fi
% Aurora Berto, 09 / 2025; aurber@utu.fi

close all

% Define the state names
stateNames = {'State 1', 'State 2', 'State 3', 'State 4', 'State 5', 'State 6'};

% Define the number of folds for cross-validation
numFolds = 10;

% Define the maximum number of PLS components to consider
maxNumComponents = 10;

% Preallocate a vector to hold the cross-validated MSE for each number of components
cvMSE = zeros(maxNumComponents, 1);
cvMSE_rand = zeros(maxNumComponents, 1);
for i_state = 1:6
    % Extract the predictors and responses for this state
    Y = V(:, i_state); % responses

    for i_comp = 1:maxNumComponents
        [~, ~, ~, ~, ~, ~, MSE_rand] = plsregress(H_rand, Y, i_comp, 'CV', numFolds);
        cvMSE_rand(i_comp) = mean(MSE_rand(2, :)); % take the mean of the second row of MSE
    end

    for i_comp = 1:maxNumComponents
        [~, ~, ~, ~, ~, ~, MSE] = plsregress(H, Y, i_comp, 'CV', numFolds);
        cvMSE(i_comp) = mean(MSE(2, :)); % take the mean of the second row of MSE
    end
% Plot the cross-validated MSE
    fig = figure;
    plot(1:maxNumComponents, cvMSE, '-o', 1:maxNumComponents, cvMSE_rand, 'r^','LineWidth',1.5);
%plot(1:maxNumComponents, cvMSE_rand, 'r^')
    xlabel('Number of PLS components','FontSize',16);
    ylabel('Cross-validated MSE','FontSize',16);
    title(['Number of PLS components vs. Cross-validated MSE for ' stateNames{i_state}],'FontSize',17);

% Find the number of components that gives the minimum MSE
    [~, optimalNumComponents] = min(cvMSE);
    disp(['The optimal number of PLS components for state ' num2str(i_state) ' is ' num2str(optimalNumComponents) ' (MSE =' num2str(min(cvMSE)) ')']);
% Save the figure
    saveas(fig, sprintf(['cross_validation%d' stateNames{i_state} '.png'], i_state));
    
    % Close the figure (optional, to save memory)
    close(fig);

end