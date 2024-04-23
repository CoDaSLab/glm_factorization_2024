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

%% REPLACE WITH PATH TO MEDA TOOLBOX
addpath('MEDA\')

%% SYNTHETIC DATA GENERATION
reps = 4;
vars = 400;
levels = {[1,2,3,4],[1,2,3]};
R = 10;
lvls = [0,1,5,10,15,20];
labels = {'UMR','CMR','pCMR'};

F = create_design(levels,reps);

X = zeros(size(F,1),vars);

for i = 1:length(levels{1})
    X(find(F(:,1) == levels{1}(i)),:) = simuleMV(length(find(F(:,1) == levels{1}(i))),vars,8) + repmat(randn(1,vars),length(find(F(:,1) == levels{1}(i))),1); %#ok
end

[~,parglmo] = parglm(X,F,{[1,2]});

true_p = parglmo.p;

results = zeros(reps,length(lvls),3);

%% ITERATIVELY REMOVE ZEROS, CALCULATE P VALUES
for ii = 1:length(lvls)

    for jj = 1:R

        X2 = X;
        rng('shuffle');
        idx = randperm(numel(X),round(numel(X)*lvls(ii)/100)); 
        X2(idx) = nan;

        [~, parglmo_grand] = parglm_grand(X2, F, [1,2],2,10000);
        [~, parglmo_cell] = parglm_cell(X2, F, [1,2],2,10000);
        [~, parglmo_cell2] = parglm_cell2(X2, F, [1,2],2,10000);

        results(jj,ii,1) = sum((parglmo_grand.p - parglmo.p)./parglmo.p);
        results(jj,ii,2) = sum((parglmo_cell.p - parglmo.p)./parglmo.p); 
        results(jj,ii,3) = sum((parglmo_cell2.p - parglmo.p)./parglmo.p); 

    fprintf('Replicate %d of %d complete \n',jj, R)
    end

fprintf('Replicates for lvl %d are complete \n', ii)

end

%% PLOT RESULTS
mns = squeeze(mean(results,1));
sds = squeeze(std(results,0,1));

plot(lvls,mns(:,1),':','LineWidth',2); hold on;
plot(lvls,mns(:,2),'--','LineWidth',2);
plot(lvls,mns(:,3),'.-',LineWidth=2);

e = errorbar(lvls,mns(:,1),sds(:,1),'LineStyle','none');
e.Color = [0 0.4470 0.7410];
e = errorbar(lvls,mns(:,2),sds(:,2),'LineStyle','none');
e.Color = [0.8500 0.3250 0.0980];
e = errorbar(lvls,mns(:,3),sds(:,3),'LineStyle','none');
e.Color = [0.9290 0.6940 0.1250];

yline(0, 'LineWidth',1.5, alpha=0.75); hold off;

title('Imputation methods, with interaction','FontSize',12)
ylabel('ERROR')
xlabel('% Missing Data')
legend(labels,'location','southwest');

savefig('resultsint.fig')
fig = gcf;
saveas(fig,'resultsint.png')






