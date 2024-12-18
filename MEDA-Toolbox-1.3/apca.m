function apcao = apca(paranovao)

% ANOVA-PCA is a data analysis algorithm for designed experiments. It does a 
% principal component analysis on the level averages of each experimental 
% factor in a designed experiment with balanced data. Interactions between 
% two factors can also be calculated.
%
% Related routines: parglm, paranova, apca, gasca, create_design
%
% apcao = apca(paranovao)   % complete call
%
%
% INPUTS:
%
% paranovao (structure): structure with the factor and interaction
% matrices, p-values and explained variance. Obtained with parallel anova
% o parallel general linear model
%
%
% OUTPUTS:
%
% apcao (structure): structure that contains scores, loadings, singular
% values and projections of the factors and interactions.
%
%
% EXAMPLE OF USE (copy and paste the code in the command line)
%   Random data, two factors, with 4 and 3 levels, but only the first one 
%   is significative, and 4 replicates:
%
% reps = 4;
% vars = 400;
% levels = {[1,2,3,4],[1,2,3]};
%
% F = create_design(levels,reps);
%
% X = zeros(size(F,1),vars);
% for i = 1:length(levels{1}),
%     X(find(F(:,1) == levels{1}(i)),:) = simuleMV(length(find(F(:,1) == levels{1}(i))),vars,8) + repmat(randn(1,vars),length(find(F(:,1) == levels{1}(i))),1);
% end
%
% [table, parglmo] = parglm(X, F);
% table
% 
% apcao = apca(parglmo);
%
% for i=1:2, % Note, the second factor is shown for the sake of illustration
%   scores(apcao.factors{i},[],[],sprintf('Factor %d',i),[],apcao.design(:,i));
% end
%
%
% EXAMPLE OF USE (copy and paste the code in the command line)
%   Random data, two significative factors, with 4 and 3 levels, and 4 replicates:
%
% reps = 4;
% vars = 400;
% levels = {[1,2,3,4],[1,2,3]};
%
% F = create_design(levels,reps);
%
% X = zeros(size(F,1),vars);
% for i = 1:length(levels{1}),
%     fi{i} = randn(1,vars);
% end
% for j = 1:length(levels{2}),
%     fj{j} = randn(1,vars);
% end
% for i = 1:length(levels{1}),
%     for j = 1:length(levels{2}),
%         X(find(F(:,1) == levels{1}(i) & F(:,2) == levels{2}(j)),:) = simuleMV(reps,vars,8) + repmat(fi{i} + fj{j},reps,1);
%     end
% end
%
% [table, parglmo] = parglm(X, F, [1 2]);
% table
% 
% apcao = apca(parglmo);
%
% for i=1:2,
%   scores(apcao.factors{i},[],[],sprintf('Factor %d',i),[],apcao.design(:,i));
% end
%
%
% EXAMPLE OF USE (copy and paste the code in the command line)
%   Random data, two factors with 4 and 3 levels, and 4 replicates, with 
%   significant interaction:
%
% reps = 4;
% vars = 400;
% levels = {[1,2,3,4],[1,2,3]};
%
% F = create_design(levels,reps);
%
% X = zeros(size(F,1),vars);
% for i = 1:length(levels{1}),
%     for j = 1:length(levels{2}),
%         X(find(F(:,1) == levels{1}(i) & F(:,2) == levels{2}(j)),:) = simuleMV(reps,vars,8) + repmat(randn(1,vars),reps,1);
%     end
% end
%
% [table, parglmo] = parglm(X, F, [1 2]);
% table
% 
% apcao = apca(parglmo);
%
% M = apcao.factors{1}.matrix + apcao.factors{2}.matrix + apcao.interactions{1}.matrix;
% code_levels = F(:,1)*10+F(:,2);
% scores_pca(M,1:2,X,0,101,[],code_levels);
% legend(num2str(unique(code_levels)))
%
%
% coded by: Jos� Camacho (josecamacho@ugr.es)
% last modification: 17/Oct/22
%
% Copyright (C) 2022  University of Granada, Granada
% Copyright (C) 2022  Jose Camacho Paez
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

%% Arguments checking

% Set default values
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);


%% Main code

apcao = paranovao;

%Do PCA on level averages for each factor
for factor = 1 : apcao.n_factors
    
    xf = apcao.factors{factor}.matrix+apcao.residuals;
    [p,t] = pca_pp(xf,1:rank(apcao.factors{factor}.matrix));
    
    apcao.factors{factor}.var = trace(xf'*xf);
    apcao.factors{factor}.lvs = 1:size(p,2);
    apcao.factors{factor}.loads = p;
    apcao.factors{factor}.scores = t;
end

%Do PCA on interactions
for interaction = 1 : apcao.n_interactions
    
    xf = apcao.interactions{interaction}.matrix+apcao.residuals;
    p = pca_pp(xf,1:rank(apcao.interactions{interaction}.matrix));
    
    apcao.factors{factor}.var = trace(xf'*xf);
    apcao.interactions{interaction}.lvs = 1:size(p,2);
    apcao.interactions{interaction}.loads = p;
    apcao.interactions{interaction}.scores = t;
end
