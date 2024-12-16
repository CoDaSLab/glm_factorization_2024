%% Simulation of missing data
%
% This script utilises a method for simulating data with significant terms
% as per Camacho et al. Camacho, Jose,
% et al. 'Least-squares approximation of a space distribution for a given
% covariance and latent sub-space.' Chemometrics and Intelligent Laboratory
% Systems 105.2 (2011): 171-180.
%
% DEPENDENCIES: MEDA Toolbox v1.3 from
% https://github.com/josecamachop/MEDA-Toolbox/tree/v1.3
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
% Jose Camacho
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
R = 10;
lvls = [0,1,5,10,15,20,35];
labels = {'UMR','CMR','pCMR'};

reps = 4;
vars = 400;
levels = {[1,2,3,4],[1,2,3]};

F = create_design(levels,reps);

X = zeros(size(F,1),vars);
for i = 1:length(levels{1})
    X(find(F(:,1) == levels{1}(i)),:) = simuleMV(length(find(F(:,1) == levels{1}(i))),vars,8) + repmat(randn(1,vars),length(find(F(:,1) == levels{1}(i))),1); %#ok
end

perm = 1e3;
[T,parglmo] = parglm(X, F, [],2,perm);

true_p = parglmo.p;

results = zeros(reps,length(lvls),3);

X_size = size(X);

row_indices_f2_lvl2 = find(all(F(:,2) == 2,2)); %only for the selection in f_2
row_indices_f2_lvl3 = find(all(F(:,2) == 3,3));
all_columns = 1:X_size(2);
[rowsl2,colsl2] = ndgrid(row_indices_f2_lvl2,all_columns);
[rowsl3,colsl3] = ndgrid(row_indices_f2_lvl3,all_columns);

idxl2 = sub2ind(X_size,rowsl2,colsl2);
idxl3 = sub2ind(X_size,rowsl3,colsl3);

num_entries = numel(X);

%% ITERATIVELY REMOVE ZEROS, CALCULATE P VALUES
for ii = 1:length(lvls)

    for jj = 1:R

        X2 = X;
        rng('shuffle');
        idx2 = randperm(numel(idxl2),round(numel(idxl2)*lvls(ii)*0.5/100)); 
        idx3 = randperm(numel(idxl3),round(numel(idxl3)*lvls(ii)/100));
        
        X2(idx2) = nan; X2(idx3) = nan;

        [~, parglmo_grand] = parglm_grand(X2, F, [],2,perm);
        [~, parglmo_cell] = parglm_cell(X2, F, [],2,perm);
        [~, parglmo_cell2] = parglm_cell2(X2, F, [],2,perm);

        results(jj,ii,1) = sum((parglmo_grand.p - parglmo.p)./parglmo.p);
        results(jj,ii,2) = sum((parglmo_cell.p - parglmo.p)./parglmo.p); 
        results(jj,ii,3) = sum((parglmo_cell2.p - parglmo.p)./parglmo.p); 
        
        lvls_actual(jj,ii) = (sum(isnan(X2),'all')/num_entries)*100;

    fprintf('Replicate %d of %d complete \n',jj, R)
    end

fprintf('Replicates for lvl %d are complete \n', ii)

end

%% PLOT RESULTS
% Calculate mean and standard deviation along the second mode
mns = squeeze(mean(results, 2)); % Mean along replicates
sds = squeeze(std(results, 0, 2)); % Standard deviation along replicates

% Colors for scatter plots (adjust or add more as needed)
colors = lines(size(mns, 2)); % Generates a colormap for unique colors

% Scatter plot results for each method in the 3rd mode
figure;
hold on;
    % Scatter plot for mean values
sc1 =    scatter(lvls_actual, results(:,:,1), 100, 'MarkerFaceColor', colors(1, :), ...
        'MarkerEdgeColor', colors(1, :),...
        'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);

sc2 =    scatter(lvls_actual, results(:,:,2), 100, 'MarkerFaceColor', colors(2, :), ...
        'MarkerEdgeColor', colors(2, :), ...
        'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);

sc3 =    scatter(lvls_actual, results(:,:,3), 100, 'MarkerFaceColor', colors(3, :), ...
        'MarkerEdgeColor', colors(3, :),...
        'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);

h = [sc1(1);sc2(1);sc3(1)];


% Zero-line
plot(lvls_actual, zeros(size(lvls_actual)), 'k:', 'LineWidth', 1.5);

% Labels and title
title('Imputation Methods, no Interaction', 'FontSize', 12);
ylabel('ERROR');
xlabel('% Missing Data');
legend(h,labels,'Location','southwest')
% Save figures
savefig('resultsint_mar.fig');
fig = gcf;
saveas(fig, 'resultsint_mar.png');

hold off;






