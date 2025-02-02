%% Simulation of missing data with a modelled interaction term.
%
% This script utilises a method for simulating data with significant terms
% as per Camacho et al. Camacho, Jos�,
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
% Jos� Camacho
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
reps = 4;
vars = 400;
levels = {[1,2,3,4],[1,2,3]}; % First factor (groups), second factor (time points)

% Create the design matrix
F = create_design(levels, reps);

% Initialize data matrix
X_original = zeros(size(F,1), vars);
X_size = size(X_original);

row_indices_f2_lvl2 = find(all(F(:,2) == 2,2)); %only for the selection in f_2
row_indices_f2_lvl3 = find(all(F(:,2) == 3,3));
all_columns = 1:size(X_original,2);
[rowsl2,colsl2] = ndgrid(row_indices_f2_lvl2,all_columns);
[rowsl3,colsl3] = ndgrid(row_indices_f2_lvl3,all_columns);

idxl2 = sub2ind(X_size,rowsl2,colsl2);
idxl3 = sub2ind(X_size,rowsl3,colsl3);

% Simulate data for each group in the first factor
for i = 1:length(levels{1})
    idx = find(F(:,1) == levels{1}(i));
    X_original(idx,:) = simuleMV(length(idx), vars, 8) + repmat(randn(1, vars), length(idx), 1);
end

%Nominal p values
[~,parglmo_lin] = parglm(X_original, F, [], 2, perm);
[~,parglmo_int] = parglm(X_original, F, [1,2], 2, perm);

% Define missingness levels
missingness_levels = [1, 5, 10, 15, 20];

num_entries = numel(X_original);


%% ITERATIVELY REMOVE ZEROS, CALCULATE P VALUES
results_lin = zeros(R,length(missingness_levels),3,2);
results_int = zeros(R,length(missingness_levels),3,3);
for ii = 1:length(missingness_levels)

    for jj = 1:R

        X2 = X_original;
        rng(jj);
        idx2 = randperm(numel(idxl2),round(numel(idxl2)*missingness_levels(ii)*0.5/100)); 
        idx3 = randperm(numel(idxl3),round(numel(idxl3)*missingness_levels(ii)/100));

        X2(idx2) = nan; X2(idx3) = nan;

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

        lvls_actual(jj,ii) = (sum(isnan(X2(:)))/num_entries)*100;

    fprintf('Replicate %d of %d complete \n',jj, R)
    end

fprintf('Replicates for lvl %d are complete \n', ii)

end

save results_mar_rand


%% Calculate mean and standard deviation for error bars

mean_lin = mean(results_lin, 1); % Mean across replicates
std_lin = std(results_lin, [], 1); % Standard deviation across replicates
mean_lvls = mean(lvls_actual, 1); % Average missingness levels across replicates - little variantion

subplot(3,2,1)
hold on;
sc1 = errorbar(mean_lvls,squeeze(mean_lin(1,:,1,1)),squeeze(std_lin(1,:,1,1)),'-o','Color',[colors(1,:),0.5]);
sc2 = errorbar(mean_lvls,squeeze(mean_lin(1,:,2,1)),squeeze(std_lin(1,:,2,1)),'-o','Color',[colors(2,:),0.5]);
sc3 = errorbar(mean_lvls,squeeze(mean_lin(1,:,3,1)),squeeze(std_lin(1,:,3,1)),'-o','Color',[colors(3,:),0.5]);
hold off;

title('A (significant)','FontSize', 10)
ylabel('ERROR','FontSize', 10);
xlabel('% Missing Data','FontSize', 10);

subplot(3,2,3)
hold on;
sc1 = errorbar(mean_lvls,squeeze(mean_lin(1,:,1,2)),squeeze(std_lin(1,:,1,2)),'-o','Color',[colors(1,:),0.5]);
sc2 = errorbar(mean_lvls,squeeze(mean_lin(1,:,2,2)),squeeze(std_lin(1,:,2,2)),'-o','Color',[colors(2,:),0.5]);
sc3 = errorbar(mean_lvls,squeeze(mean_lin(1,:,3,2)),squeeze(std_lin(1,:,3,2)),'-o','Color',[colors(3,:),0.5]);
hold off;

title('B (not significant)','FontSize', 10)
ylabel('ERROR','FontSize', 10);
xlabel('% Missing Data','FontSize', 10);

mean_int = mean(results_int, 1); % Mean across replicates
std_int = std(results_int, [], 1); % Standard deviation across replicates
mean_lvls = mean(lvls_actual, 1); % Average missingness levels across replicates - little variantion

subplot(3,2,2)
hold on;
sc1 = errorbar(mean_lvls,squeeze(mean_int(1,:,1,1)),squeeze(std_int(1,:,1,1)),'-o','Color',[colors(1,:),0.5]);
sc2 = errorbar(mean_lvls,squeeze(mean_int(1,:,2,1)),squeeze(std_int(1,:,2,1)),'-o','Color',[colors(2,:),0.5]);
sc3 = errorbar(mean_lvls,squeeze(mean_int(1,:,3,1)),squeeze(std_int(1,:,3,1)),'-o','Color',[colors(3,:),0.5]);
hold off;

title('A (significant)','FontSize', 10)
ylabel('ERROR','FontSize', 10);
xlabel('% Missing Data','FontSize', 10);

subplot(3,2,4)
hold on;
sc1 = errorbar(mean_lvls,squeeze(mean_int(1,:,1,2)),squeeze(std_int(1,:,1,2)),'-o','Color',[colors(1,:),0.5]);
sc2 = errorbar(mean_lvls,squeeze(mean_int(1,:,2,2)),squeeze(std_int(1,:,2,2)),'-o','Color',[colors(2,:),0.5]);
sc3 = errorbar(mean_lvls,squeeze(mean_int(1,:,3,2)),squeeze(std_int(1,:,3,2)),'-o','Color',[colors(3,:),0.5]);
hold off;

title('B (not significant)','FontSize', 10)
ylabel('ERROR','FontSize', 10);
xlabel('% Missing Data','FontSize', 10);

subplot(3,2,6)
hold on;
sc1 = errorbar(mean_lvls,squeeze(mean_int(1,:,1,3)),squeeze(std_int(1,:,1,3)),'-o','Color',[colors(1,:),0.5]);
sc2 = errorbar(mean_lvls,squeeze(mean_int(1,:,2,3)),squeeze(std_int(1,:,2,3)),'-o','Color',[colors(2,:),0.5]);
sc3 = errorbar(mean_lvls,squeeze(mean_int(1,:,3,3)),squeeze(std_int(1,:,3,3)),'-o','Color',[colors(3,:),0.5]);
hold off;

title('AB (not significant)','FontSize', 10)
ylabel('ERROR','FontSize', 10);
xlabel('% Missing Data','FontSize', 10);

%sgtitle('Imputation Methods, haphazard MAR Missingness', 'FontSize', 14);

h = [sc1(1);sc2(1);sc3(1)]; 
lgd = legend(h,labels);
set(lgd, 'Position', [0.2, 0.1, 0.2, 0.2]);

% Save figures
savefig('results_mar_rand.fig');
fig = gcf;
saveas(fig, 'results_mar_rand.png');










