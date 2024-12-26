%% Simulation of missing data relative to a General Linear Model
%
% This script utilises a method for simulating data with significant terms
% as per Camacho, José,
% et al. 'Least-squares approximation of a space distribution for a given
% covariance and latent sub-space.' Chemometrics and Intelligent Laboratory
% Systems 105.2 (2011): 171-180.
%
% DEPENDENCIES: MEDA Toolbox v1.4 from
% https://github.com/josecamachop/MEDA-Toolbox/tree/v1.4
%
% The data is simulated according to an n x m matrix of observations and
% features as: X = A + E 
% where E is the error term, which is in this case
% randomly generated. The data is modelled according to: 
% X = A + B 
% where B is a second linear term in an ANOVA-like model, The effect of missing
% values, and the different methods for missing value imputation relative
% to the true p values is measured in this script.
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

clear all; close all;

%% SYNTHETIC DATA GENERATION
colors = lines(3);
reps = 4;
vars = 400;
levels = {[1,2,3,4],[1,2,3]};
R = 10;
lvls = [1, 5, 10, 15, 20];
labels = {'UMR','CMR','pCMR'};

F = create_design(levels,reps);

X = zeros(size(F,1),vars);
for i = 1:length(levels{1})
    X(find(F(:,1) == levels{1}(i)),:) = simuleMV(length(find(F(:,1) == levels{1}(i))),vars,8) + repmat(randn(1,vars),length(find(F(:,1) == levels{1}(i))),1); %#ok
end

perm = 1e3;
[~,parglmo_lin] = parglm(X,F,[],2,perm);
[~,parglmo_int] = parglm(X,F,[1,2],2,perm);

results = zeros(reps,length(lvls),3);
num_entries = numel(X);


%% ITERATIVELY REMOVE ZEROS, CALCULATE P VALUES
results_lin = zeros(R,length(lvls),3,2);
results_int = zeros(R,length(lvls),3,3);
for ii = 1:length(lvls)

    for jj = 1:R

        X2 = X;
        rng(jj);
        idx = randperm(numel(X),round(numel(X)*lvls(ii)/100));
        X2(idx) = nan;

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

save results_mcar_rand

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

%sgtitle('Imputation Methods, haphazard MCAR Missingness', 'FontSize', 14);

h = [sc1(1);sc2(1);sc3(1)]; 
lgd = legend(h,labels);
set(lgd, 'Position', [0.2, 0.1, 0.2, 0.2]);

% Save figures
savefig('results_mcar_rand.fig');
fig = gcf;
saveas(fig, 'results_mcar_rand.png');







