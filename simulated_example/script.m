%% Example 1 of Considerations for missing data, outliers and 
% transformations in permutation testing for ANOVA, ASCA and related 
% factorizations, by Polushkina et al. Submitted to Chemolab
% Simulation with two factors, one significant, and no significant 
% interaction. Missing data is artificially generated in the data and
% imputed with Unconditional Mean Replacement (UMR), Conditional Mean 
% Replacement (CMR) and Trimmed Score Regression (TSR). 
%
% coded by: José Camacho (josecamacho@ugr.es)
% last modification: 09/Apr/24
%
% Copyright (C) 2024  Universidad de Granada
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


for rep = 1:10
    
rep
    
rng(rep);

reps = 4;
vars = 400;
levels = {[1,2,3,4],[1,2,3]};

F = create_design(levels,'Replicates',reps);

X = zeros(size(F,1),vars);
for i = 1:length(levels{1})
    X(find(F(:,1) == levels{1}(i)),:) = simuleMV(length(find(F(:,1) == levels{1}(i))),vars,8) + repmat(randn(1,vars),length(find(F(:,1) == levels{1}(i))),1);
end

T = parglm(X, F, [], 0);


%% Create +5% missing elements and compute SS

A = randn(size(X));
[~,ind] = sort(A(:));

Xmis = X;
Xmis(ind(1:1000)) = nan;

m = nanmean(Xmis);

u1 = unique(F(:,1));
for f1 = 1:length(u1)
    ind = find(F(:,1)==u1(f1));
    mf1(f1,:) = nanmean(Xmis(ind,:)) - m;
    nf1(f1,:) = sum(~isnan(Xmis(ind,:)));  
end     
SSf1 = sum(sum(nf1.*(mf1.^2)));

u2 = unique(F(:,2));
for f2 = 1:length(u2)
    ind = find(F(:,2)==u2(f2));
    mf2(f2,:) = nanmean(Xmis(ind,:)) - m;
    nf2(f2,:) = sum(~isnan(Xmis(ind,:))); 
end     
SSf2 = sum(sum(nf2.*(mf2.^2)));

Xmisc = Xmis;
SStotal = sum(Xmisc(~isnan(Xmisc(:))).^2);
SSerr = SStotal - SSf1 - SSf2;

Tmiss = [SSf1 SSf2 SSerr SStotal]';

%% Impute with UMR

[row,col] = find(isnan(Xmis));

Xumr = Xmis;
cv = unique(col);
for c = 1:length(cv)
    ind = find(col==cv(c));
    Xumr(row(ind),cv(c)) = nanmean(Xmis(:,cv(c)));
end

Tumr = parglm(Xumr, F, [], 0);


%% Impute with CMR

[row,col] = find(isnan(Xmis));

Xcmr = Xmis;
fv = unique(F,'rows');
cv = unique(col);
for c = 1:length(cv)
    for f = 1:length(fv)
        indf = ismember(fv(f,:),F,'rows'); 
        indc = find(col(indf)==cv(c));
        Xcmr(row(indf(indc)),cv(c)) = nanmean(Xmis(indf,cv(c)));
    end
end

Tcmr = parglm(Xcmr, F, [], 0);


%% Impute with TSR

[row,col] = find(isnan(Xmis));

Xtsr = missTSR2D(Xmis,1:3,0);

Ttsr = parglm(Xtsr, F, [], 0);


%% Comparison of SS distribution

Trep(:,:,rep) = [T{2:end,2},Tmiss,Tumr{2:end,2},Tcmr{2:end,2},Ttsr{2:end,2}];

end

Ttot = array2table(mean(Trep,3),'RowNames', T{2:end,1},'VariableNames', {'Original','Missing','UMR','CMR','TSR'})

