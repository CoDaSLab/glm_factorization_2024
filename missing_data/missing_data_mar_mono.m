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
% coded by: Michael Sorochan Armstrong (mdarmstr@go.ugr.es) 
%           José Camacho (josecamacho@ugr.es)
% last modification: 23-12-2024
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
colors = lines(4);
R = 10;
perm = 1e3;
labels = {'UMR','CMR','pCMR','pCMR_2'};
reps = 25;
vars = 400;
levels = {[1,2,3,4],[1,2,3],1:reps}; % First factor (groups), second factor (time points)

% Create the design matrix
F2 = create_design(levels, 1);
F = F2(:,1:end-1);

% Initialize data matrix
X_original = zeros(size(F,1), vars);

% Simulate data for each group in the first factor
for i = 1:length(levels{1})
    idx = find(F(:,1) == levels{1}(i));
    X_original(idx,:) = simuleMV(length(idx), vars, 8) + repmat(randn(1, vars), length(idx), 1);
end

% Initialize data matrix
X_original2 = zeros(size(F,1), vars);

% Simulate data for each group in the second factor
for i = 1:length(levels{2})
    idx = find(F(:,2) == levels{2}(i));
    X_original2(idx,:) = simuleMV(length(idx), vars, 8) + repmat(randn(1, vars), length(idx), 1);
end

%Nominal p values
[~,parglmo_lin] = parglm(X_original, F, [], 2, perm);
[~,parglmo_lin2] = parglm(X_original2, F, [], 2, perm);

% Define missingness levels
missingness_levels = [1, 5, 10];

num_entries = numel(X_original);
    
patient_index = (F(:,1)-1)*reps + F2(:,3);

%% Loop over missingness levels
results_lin = zeros(R,length(missingness_levels),3,2);
results_lin2 = zeros(R,length(missingness_levels),3,2);
for ii = 1:length(missingness_levels)

    num_patients = reps*length(levels{1}); 
    
    num_dropout_patients = round((missingness_levels(ii) / 100) * num_patients); 

    for jj = 1:R
        % Copy the original data matrix for this level of missingness
        X = X_original;
        X2 = X_original2;
        rng(jj);
        % Select random patients for dropout
        dropout_patients = randperm(num_patients, num_dropout_patients);

        % Introduce monotonic missing data for selected dropout patients
        for kk = dropout_patients
            % Find rows corresponding to the selected patient
            patient_rows = find(patient_index==kk); % Find patient rows (based on design matrix)

            dropout_time = 2; 

            % Introduce missingness starting from the dropout time point onward
            for ll = dropout_time:length(levels{2})
                time_idx = find(F(:,2) == levels{2}(ll)); % Find rows for the time point
                patient_time_rows = intersect(patient_rows, time_idx); % Get specific rows for the patient at this time point

                X(patient_time_rows,:) = NaN;
                X2(patient_time_rows,:) = NaN;
            end
        end

        [~, parglmo_grand_lin] = parglm_grand(X, F, [],2,perm);
        [~, parglmo_cell_lin] = parglm_cell(X, F, [],2,perm);
        [~, parglmo_cell2_lin] = parglm_cell2(X, F, [],2,perm);
        [~, parglmo_cell3_lin] = parglm_cell3(X, F, patient_index, [],2,perm);

        [~, parglmo_grand_lin2] = parglm_grand(X2, F, [],2,perm);
        [~, parglmo_cell_lin2] = parglm_cell(X2, F, [],2,perm);
        [~, parglmo_cell2_lin2] = parglm_cell2(X2, F, [],2,perm);
        [~, parglmo_cell3_lin2] = parglm_cell3(X2, F, patient_index, [],2,perm);

        results_lin(jj,ii,1,:) = (parglmo_grand_lin.p - parglmo_lin.p)./parglmo_lin.p;
        results_lin(jj,ii,2,:) = (parglmo_cell_lin.p - parglmo_lin.p)./parglmo_lin.p;
        results_lin(jj,ii,3,:) = (parglmo_cell2_lin.p - parglmo_lin.p)./parglmo_lin.p;
        results_lin(jj,ii,4,:) = (parglmo_cell3_lin.p - parglmo_lin.p)./parglmo_lin.p;

        results_lin2(jj,ii,1,:) = (parglmo_grand_lin2.p - parglmo_lin2.p)./parglmo_lin2.p;
        results_lin2(jj,ii,2,:) = (parglmo_cell_lin2.p - parglmo_lin2.p)./parglmo_lin2.p;
        results_lin2(jj,ii,3,:) = (parglmo_cell2_lin2.p - parglmo_lin2.p)./parglmo_lin2.p;
        results_lin2(jj,ii,4,:) = (parglmo_cell3_lin2.p - parglmo_lin2.p)./parglmo_lin2.p;

        lvls_actual(jj,ii) = (sum(isnan(X2(:)))/num_entries)*100;

        fprintf('Replicate %d of %d complete \n',jj, R)
    end

    fprintf('Replicates for lvl %d are complete \n', ii)

    % Save or analyze the data matrix `X` for this missingness level
    fprintf('Missingness level: %d%% completed.\n', missingness_levels(ii));

    % Optionally save the dataset for further analysis
    %save(sprintf('X_missingness_%d.mat', miss_level), 'X', 'F');
    
    save results_mar_mono

end

%% Calculate mean and standard deviation for error bars
mean_lin = mean(results_lin, 1); % Mean across replicates
std_lin = std(results_lin, [], 1); % Standard deviation across replicates
mean_lvls = mean(lvls_actual, 1); % Average missingness levels across replicates - little variantion

subplot(3,2,1)
hold on;
sc1 = errorbar(missingness_levels,squeeze(mean_lin(1,:,1,1)),squeeze(std_lin(1,:,1,1)),'-o','Color',[colors(1,:),0.5]);
sc2 = errorbar(missingness_levels,squeeze(mean_lin(1,:,2,1)),squeeze(std_lin(1,:,2,1)),'-o','Color',[colors(2,:),0.5]);
sc3 = errorbar(missingness_levels,squeeze(mean_lin(1,:,3,1)),squeeze(std_lin(1,:,3,1)),'-o','Color',[colors(3,:),0.5]);
sc4 = errorbar(missingness_levels,squeeze(mean_lin(1,:,4,1)),squeeze(std_lin(1,:,4,1)),'-o','Color',[colors(4,:),0.5]);
hold off;

title('A (significant)','FontSize', 10)
ylabel('ERROR','FontSize', 10);
xlabel('% Missing Data','FontSize', 10);

subplot(3,2,3)
hold on;
sc1 = errorbar(missingness_levels,squeeze(mean_lin(1,:,1,2)),squeeze(std_lin(1,:,1,2)),'-o','Color',[colors(1,:),0.5]);
sc2 = errorbar(missingness_levels,squeeze(mean_lin(1,:,2,2)),squeeze(std_lin(1,:,2,2)),'-o','Color',[colors(2,:),0.5]);
sc3 = errorbar(missingness_levels,squeeze(mean_lin(1,:,3,2)),squeeze(std_lin(1,:,3,2)),'-o','Color',[colors(3,:),0.5]);
sc4 = errorbar(missingness_levels,squeeze(mean_lin(1,:,4,2)),squeeze(std_lin(1,:,4,2)),'-o','Color',[colors(4,:),0.5]);
hold off;

title('B (not significant)','FontSize', 10)
ylabel('ERROR','FontSize', 10);
xlabel('% Missing Data','FontSize', 10);

mean_lin2 = mean(results_lin2, 1); % Mean across replicates
std_lin2 = std(results_lin2, [], 1); % Standard deviation across replicates
mean_lvls = mean(lvls_actual, 1); % Average missingness levels across replicates - little variantion

subplot(3,2,2)
hold on;
sc1 = errorbar(missingness_levels,squeeze(mean_lin2(1,:,1,1)),squeeze(std_lin2(1,:,1,1)),'-o','Color',[colors(1,:),0.5]);
sc2 = errorbar(missingness_levels,squeeze(mean_lin2(1,:,2,1)),squeeze(std_lin2(1,:,2,1)),'-o','Color',[colors(2,:),0.5]);
sc3 = errorbar(missingness_levels,squeeze(mean_lin2(1,:,3,1)),squeeze(std_lin2(1,:,3,1)),'-o','Color',[colors(3,:),0.5]);
sc4 = errorbar(missingness_levels,squeeze(mean_lin2(1,:,4,1)),squeeze(std_lin2(1,:,4,1)),'-o','Color',[colors(4,:),0.5]);
hold off;

title('A (not significant)','FontSize', 10)
ylabel('ERROR','FontSize', 10);
xlabel('% Missing Data','FontSize', 10);

subplot(3,2,4)
hold on;
sc1 = errorbar(missingness_levels,squeeze(mean_lin2(1,:,1,2)),squeeze(std_lin2(1,:,1,2)),'-o','Color',[colors(1,:),0.5]);
sc2 = errorbar(missingness_levels,squeeze(mean_lin2(1,:,2,2)),squeeze(std_lin2(1,:,2,2)),'-o','Color',[colors(2,:),0.5]);
sc3 = errorbar(missingness_levels,squeeze(mean_lin2(1,:,3,2)),squeeze(std_lin2(1,:,3,2)),'-o','Color',[colors(3,:),0.5]);
sc4 = errorbar(missingness_levels,squeeze(mean_lin2(1,:,4,2)),squeeze(std_lin2(1,:,4,2)),'-o','Color',[colors(4,:),0.5]);
hold off;

title('B (significant)','FontSize', 10)
ylabel('ERROR','FontSize', 10);
xlabel('% Missing Data','FontSize', 10);


%sgtitle('Imputation Methods, monotone MAR Simulation', 'FontSize', 14);

h = [sc1(1);sc2(1);sc3(1)]; 
lgd = legend(h,labels);
set(lgd, 'Position', [0.2, 0.1, 0.2, 0.2]);

% Save figures
savefig('results_mar_mono.fig');
fig = gcf;
saveas(fig, 'results_mar_mono.png');








