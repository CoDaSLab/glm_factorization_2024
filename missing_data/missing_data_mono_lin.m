%% Simulation of missing data with a modelled interaction term.
%
% This script utilises a method for simulating data with significant terms
% as per Camacho et al. Camacho, José,
% et al. 'Least-squares approximation of a space distribution for a given
% covariance and latent sub-space.' Chemometrics and Intelligent Laboratory
% Systems 105.2 (2011): 171-180.
%
% DEPENDENCIES: MEDA Toolbox v1.4 from
% https://github.com/josecamachop/MEDA-Toolbox/tree/v1.4
%
% The data is simulated according to an n x m matrix of observations and
% features as:
% X = A + E
% where E is the error term, which is in this case randomly generated. The
% data is modelled according to:
% X = A + B + AB
% where B is a second linear term in an ANOVA-like model, and AB is the
% interacting term. The effect of missing values, and the different methods
% for missing value imputation relative to the true p values is measured in
% this script. 
%
% coded by: Michael Sorochan Armstrong (mdarmstr@go.ugr.es) with input from
% José Camacho
% last modification: 17-04-2023
%
% Copyright (C) 2023  Universidad de Granada
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


%% SYNTHETIC DATA GENERATION
close all; clear all; %#ok

% Useful parameters
colors = lines(3);
R = 10;
perm = 1e3;
labels = {'UMR','CMR','pCMR'};
reps = 8;
vars = 400;
levels = {[1,2,3,4],[1,2,3]}; % First factor (groups), second factor (time points)

% Create the design matrix
F = create_design(levels, reps);

% Initialize data matrix
X_original = zeros(size(F,1), vars);

% Simulate data for each group in the first factor
for i = 1:length(levels{1})
    idx = find(F(:,1) == levels{1}(i));
    X_original(idx,:) = simuleMV(length(idx), vars, 8) + repmat(randn(1, vars), length(idx), 1);
end

%Nominal p values
[~,parglmo_lin] = parglm(X_original, F, [], 2, perm);
[~,parglmo_int] = parglm(X_original, F, [1,2], 2, perm);

% Define missingness levels
missingness_levels = [0, 1, 5, 10, 15, 20];

num_entries = numel(X_original);

% Loop over missingness levels
for ii = 1:length(missingness_levels)
    % Copy the original data matrix for this level of missingness
    X = X_original;

    % Calculate the percentage of patients to drop out based on the missingness level
    num_patients = size(F, 1) / length(levels{2}); % Number of unique patients
    % Calculate number of patients that will drop out - at least one for 1% missingness
    num_dropout_patients = ceil((missingness_levels(ii) / 100) * num_patients); 

    % Select random patients for dropout
    dropout_patients = randperm(num_patients, num_dropout_patients);

    % Introduce monotonic missing data for selected dropout patients
    for p = dropout_patients
        % Find rows corresponding to the selected patient
        patient_rows = find(F(:,1) == ceil(p / reps)); % Find patient rows (based on design matrix)

        % Randomly select a time point for dropout - either 2 or 3. Patient
        % must be observed at least once for this to make sense.
        dropout_time = randi([2, length(levels{2})]);

        % Introduce missingness starting from the dropout time point onward
        for t = dropout_time:length(levels{2})
            time_idx = find(F(:,2) == levels{2}(t)); % Find rows for the time point
            patient_time_rows = intersect(patient_rows, time_idx); % Get specific rows for the patient at this time point
            
            % Calculate the subset of `patient_time_rows` to drop out
            % Drop out a percentage of rows within the patient_time_rows
            num_missing_rows = min(ceil(missingness_levels(ii) * length(patient_time_rows)), length(patient_time_rows)); % Cap to available rows
            if num_missing_rows > 0
                missing_rows = randsample(patient_time_rows, num_missing_rows); % Random subset
                % Set the selected rows to NaN
                X(missing_rows, :) = NaN;
            end
        end
    end
    for jj = 1:R

        X2 = X;

        [~, parglmo_grand_lin] = parglm_grand(X2, F, [],2,perm);
        [~, parglmo_cell_lin] = parglm_cell(X2, F, [],2,perm);
        [~, parglmo_cell2_lin] = parglm_cell2(X2, F, [],2,perm);

        [~, parglmo_grand_int] = parglm_grand(X2, F, [1,2],2,perm);
        [~, parglmo_cell_int] = parglm_cell(X2, F, [1,2],2,perm);
        [~, parglmo_cell2_int] = parglm_cell2(X2, F, [1,2],2,perm);

        results_lin(jj,ii,1,:) = (parglmo_grand_lin.p - parglmo_lin.p)./parglmo_lin.p;
        results_lin(jj,ii,2,:) = (parglmo_cell_lin.p - parglmo_lin.p)./parglmo_lin.p;
        results_lin(jj,ii,3,:) = (parglmo_cell2_lin.p - parglmo_lin.p)./parglmo_lin.p;

        results_int(jj,ii,1,:) = (parglmo_grand_int.p - parglmo_int.p)./parglmo_int.p;
        results_int(jj,ii,2,:) = (parglmo_cell_int.p - parglmo_int.p)./parglmo_int.p;
        results_int(jj,ii,3,:) = (parglmo_cell2_int.p - parglmo_int.p)./parglmo_int.p;

        lvls_actual(jj,ii) = (sum(isnan(X2),'all')/num_entries)*100;

        fprintf('Replicate %d of %d complete \n',jj, R)
    end

    fprintf('Replicates for lvl %d are complete \n', ii)

    % Save or analyze the data matrix `X` for this missingness level
    fprintf('Missingness level: %d%% completed.\n', missingness_levels(ii));

    % Optionally save the dataset for further analysis
    %save(sprintf('X_missingness_%d.mat', miss_level), 'X', 'F');
end


%% ITERATIVELY REMOVE ZEROS, CALCULATE P VALUES
% for ii = 1:length(lvls)
% 
%     for jj = 1:R
% 
%         X2 = X;
%         rng('shuffle');
%         idx2 = randperm(numel(idxl2),round(numel(idxl2)*lvls(ii)*0.5/100)); 
%         idx3 = randperm(numel(idxl3),round(numel(idxl3)*lvls(ii)/100));
% 
%         X2(idx2) = nan; X2(idx3) = nan;
% 
%         [~, parglmo_grand] = parglm_grand(X2, F, [],2,perm);
%         [~, parglmo_cell] = parglm_cell(X2, F, [],2,perm);
%         [~, parglmo_cell2] = parglm_cell2(X2, F, [],2,perm);
% 
%         results(jj,ii,1,:) = (parglmo_grand.p - parglmo.p)./parglmo.p;
%         results(jj,ii,2,:) = (parglmo_cell.p - parglmo.p)./parglmo.p; 
%         results(jj,ii,3,:) = (parglmo_cell2.p - parglmo.p)./parglmo.p; 
% 
%         lvls_actual(jj,ii) = (sum(isnan(X2),'all')/num_entries)*100;
% 
%     fprintf('Replicate %d of %d complete \n',jj, R)
%     end
% 
% fprintf('Replicates for lvl %d are complete \n', ii)
% 
% end

subplot(3,2,1)
hold on;
sc1 = scatter(lvls_actual,results_lin(:,:,1,1),50,'MarkerFaceColor', colors(1, :), ...
    'MarkerEdgeColor', colors(1, :), ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'DisplayName', 'Method 1');
sc2 = scatter(lvls_actual,results_lin(:,:,2,1),50,'MarkerFaceColor', colors(2, :), ...
    'MarkerEdgeColor', colors(2, :), ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'DisplayName', 'Method 1');
sc3 = scatter(lvls_actual,results_lin(:,:,3,1),50,'MarkerFaceColor', colors(3, :), ...
    'MarkerEdgeColor', colors(3, :), ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'DisplayName', 'Method 1');
hold off;

%h = [sc1(1);sc2(1);sc3(1)];
%legend(h,labels,'Location','southwest')
subtitle('2 Factor Model - A (significant)')
ylabel('ERROR');
xlabel('% Missing Data');

subplot(3,2,3)
hold on;
sc1 = scatter(lvls_actual,results_lin(:,:,1,2),50,'MarkerFaceColor', colors(1, :), ...
    'MarkerEdgeColor', colors(1, :), ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'DisplayName', 'Method 1');
sc2 = scatter(lvls_actual,results_lin(:,:,2,2),50,'MarkerFaceColor', colors(2, :), ...
    'MarkerEdgeColor', colors(2, :), ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'DisplayName', 'Method 1');
sc3 = scatter(lvls_actual,results_lin(:,:,3,2),50,'MarkerFaceColor', colors(3, :), ...
    'MarkerEdgeColor', colors(3, :), ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'DisplayName', 'Method 1');
hold off;

%h = [sc1(1);sc2(1);sc3(1)];
%legend(h,labels,'Location','southwest')
subtitle('2 Factor Model - B (not significant)')
ylabel('ERROR');
xlabel('% Missing Data');

sgtitle('Imputation Methods, No Interaction', 'FontSize', 14);

% Save figures
savefig('resultslin_mono.fig');
fig = gcf;
saveas(fig, 'resultslin_mono.png');


% %% WITH INTERACTION MODELLED
% [T,parglmo] = parglm(X, F, [1,2], 2, perm);
% 
% true_p = parglmo.p;
% 
% results = zeros(reps,length(lvls),3,length(true_p));
% 
% X_size = size(X);
% 
% row_indices_f2_lvl2 = find(all(F(:,2) == 2,2)); %only for the selection in f_2
% row_indices_f2_lvl3 = find(all(F(:,2) == 3,3));
% all_columns = 1:X_size(2);
% [rowsl2,colsl2] = ndgrid(row_indices_f2_lvl2,all_columns);
% [rowsl3,colsl3] = ndgrid(row_indices_f2_lvl3,all_columns);
% 
% idxl2 = sub2ind(X_size,rowsl2,colsl2);
% idxl3 = sub2ind(X_size,rowsl3,colsl3);
% 
% num_entries = numel(X);
% 
% %% ITERATIVELY REMOVE ZEROS, CALCULATE P VALUES
% for ii = 1:length(lvls)
% 
%     for jj = 1:R
% 
%         X2 = X;
%         rng('shuffle');
%         idx2 = randperm(numel(idxl2),round(numel(idxl2)*lvls(ii)*0.5/100)); 
%         idx3 = randperm(numel(idxl3),round(numel(idxl3)*lvls(ii)/100));
% 
%         X2(idx2) = nan; X2(idx3) = nan;
% 
%         [~, parglmo_grand] = parglm_grand(X2, F, [1,2],2,perm);
%         [~, parglmo_cell] = parglm_cell(X2, F, [1,2],2,perm);
%         [~, parglmo_cell2] = parglm_cell2(X2, F, [1,2],2,perm);
% 
%         results(jj,ii,1,:) = (parglmo_grand.p - parglmo.p)./parglmo.p;
%         results(jj,ii,2,:) = (parglmo_cell.p - parglmo.p)./parglmo.p; 
%         results(jj,ii,3,:) = (parglmo_cell2.p - parglmo.p)./parglmo.p; 
% 
%         lvls_actual(jj,ii) = (sum(isnan(X2),'all')/num_entries)*100;
% 
%     fprintf('Replicate %d of %d complete \n',jj, R)
%     end
% 
% fprintf('Replicates for lvl %d are complete \n', ii)
% 
% end
% 
subplot(3,2,2)
hold on;
sc1 = scatter(lvls_actual,results_int(:,:,1,1),50,'MarkerFaceColor', colors(1, :), ...
    'MarkerEdgeColor', colors(1, :), ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'DisplayName', 'Method 1');
sc2 = scatter(lvls_actual,results_int(:,:,2,1),50,'MarkerFaceColor', colors(2, :), ...
    'MarkerEdgeColor', colors(2, :), ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'DisplayName', 'Method 1');
sc3 = scatter(lvls_actual,results_int(:,:,3,1),50,'MarkerFaceColor', colors(3, :), ...
    'MarkerEdgeColor', colors(3, :), ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'DisplayName', 'Method 1');
hold off;

%h = [sc1(1);sc2(1);sc3(1)];
%legend(h,labels,'Location','southwest')
subtitle('2 Factor Model - A (significant)')
ylabel('ERROR');
xlabel('% Missing Data');

subplot(3,2,4)
hold on;
sc1 = scatter(lvls_actual,results_int(:,:,1,2),50,'MarkerFaceColor', colors(1, :), ...
    'MarkerEdgeColor', colors(1, :), ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'DisplayName', 'Method 1');
sc2 = scatter(lvls_actual,results_int(:,:,2,2),50,'MarkerFaceColor', colors(2, :), ...
    'MarkerEdgeColor', colors(2, :), ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'DisplayName', 'Method 1');
sc3 = scatter(lvls_actual,results_int(:,:,3,2),50,'MarkerFaceColor', colors(3, :), ...
    'MarkerEdgeColor', colors(3, :), ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'DisplayName', 'Method 1');
hold off;

% h = [sc1(1);sc2(1);sc3(1)];
% legend(h,labels,'Location','southwest')
subtitle('2 Factor Model - B (not significant)')
ylabel('ERROR');
xlabel('% Missing Data');

subplot(3,2,6)
hold on;
sc1 = scatter(lvls_actual,results_int(:,:,1,3),50,'MarkerFaceColor', colors(1, :), ...
    'MarkerEdgeColor', colors(1, :), ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'DisplayName', 'Method 1');
sc2 = scatter(lvls_actual,results_int(:,:,2,3),50,'MarkerFaceColor', colors(2, :), ...
    'MarkerEdgeColor', colors(2, :), ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'DisplayName', 'Method 1');
sc3 = scatter(lvls_actual,results_int(:,:,3,3),50,'MarkerFaceColor', colors(3, :), ...
    'MarkerEdgeColor', colors(3, :), ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'DisplayName', 'Method 1');
hold off;

h = [sc1(1);sc2(1);sc3(1)];

lgd = legend(h,labels);
set(lgd, 'Position', [0.2, 0.1, 0.2, 0.2]); % Adjust this for precise centering
subtitle('2 Factor Model - AB (not significant)')
ylabel('ERROR');
xlabel('% Missing Data');


% Save figures
savefig('resultsint_mono.fig');
fig = gcf;
saveas(fig, 'resultsint_mono.png');










