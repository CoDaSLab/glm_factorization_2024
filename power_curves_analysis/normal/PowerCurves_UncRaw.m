%% Computiation of Power Curves for unconstrained permutations on raw 
% observations, residuals follow a uniform distribution.
% 
% Experimental Design: X = m + A + B + AB + C(A) 
% where: 
%   - A: Treatment 
%   - B: Time 
%   - C(A): Subject, nested in Treatment
%   - AB: Interaction A & B
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 17/Jul/2023
%
% Copyright (C) 2023  University of Granada, Granada
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

%% Inicialization

clear
close all
clc

tic

% Import your own F matrix
% load ../blood_hemo_vasca 
 load ../biomark_vasca 
% load ../bacteria_vasca  

obs_l = 1:size(F);

uc = unique(F(:,1));
ut = unique(F(:,2));
up = unique(F(:,3));

perm_tot = 100;

%% Repite permutations

alpha = 0:0.005:0.05; % This controls the compromise of true significance vs random
eD = zeros(5,length(alpha),4,perm_tot);
addpath '..'

for i2=1:perm_tot
    
    disp(i2)
    
    rng(i2);
    
    % Build data
    
    Xpac = randn(length(up),length(var_l));
    Xpac = Xpac/norm(Xpac);
    Xtime = randn(length(ut),length(var_l));
    Xtime = Xtime/norm(Xtime);
    for i = 1:length(obs_l)
        if F(i,1)==2
            Xstruct(i,:) = Xpac(F(i,3),:) + Xtime(F(i,2),:);
        else
            Xstruct(i,:) = Xpac(F(i,3),:);
        end
    end
    
    Xnoise = randn(length(obs_l),length(var_l)); % This is the only change in this branch
    Xnoise = Xnoise/norm(Xnoise);

    for a = 1:length(alpha)
        
        % Make a blend
         
        Xm = alpha(a)*Xstruct + (1-alpha(a))*Xnoise; 
        
        % raw data perm
        [T, paranovao] = parglm(Xm,F,[1 2],1,[],0);
        reo = [4 2 1 3];
        for o = 1:length(reo)
            if paranovao.p(reo(o))<0.05
                eD(1,a,o,i2) = 1;
            end
        end

        % raw data perm F
        [T, paranovao] = parglm(Xm,F,[1 2],1,[],1);
        reo = [4 2 1 3];
        for o = 1:length(reo)
            if paranovao.p(reo(o))<0.05
                eD(2,a,o,i2) = 1;
            end
        end
        
        % raw data perm F nested
        [T, paranovao] = parglm(Xm,F,[1 2],0,[],1,[],[],[],[1 3]);
        reo = [4 2 1 3];
        for o = 1:length(reo)
            if paranovao.p(reo(o))<0.05
                eD(3,a,o,i2) = 1;
            end
        end
        
        % raw data perm F after boxcox
        Xm2 = Xm;
        Xm2 = Xm2 - min(min(Xm)) + 0.01;
        for i=1:size(Xm,2), Xm2(:,i) = boxcox(Xm2(:,i)); end
        [T, paranovao] = parglm(Xm2,F,[1 2],1,[],1);
        reo = [4 2 1 3];
        for o = 1:length(reo)
            if paranovao.p(reo(o))<0.05
                eD(4,a,o,i2) = 1;
            end
        end
        
        % raw data perm F after rank transform
        [T, paranovao] = parglm(rank_transform(Xm),F,[1 2],1,[],1);
        reo = [4 2 1 3];
        for o = 1:length(reo)
            if paranovao.p(reo(o))<0.05
                eD(5,a,o,i2) = 1;
            end
        end
        
        % Combine 2 and 5
        Xr =  rank_transform(Xm);
        [T, paranovao] = parglm([Xm/norm(Xm) Xr/norm(Xr)],F,[1 2],1,[],1);
        reo = [4 2 1 3];
        for o = 1:length(reo)
            if paranovao.p(reo(o))<0.05
                eD(6,a,o,i2) = 1;
            end
        end
        
    end
end

toc

e = mean(eD,4);
eST = std(eD,0,4);
save random_e e eST eD alpha

%% Compare Interaction

i=1;
figure, plot(alpha,e(1,:,i),'b'); hold on, plot(alpha,e(2,:,i),'r'); plot(alpha,e(3,:,i),'g.-'); plot(alpha,e(4,:,i),'k--'); plot(alpha,e(5,:,i),'m:'); plot(alpha,e(6,:,i),'co-'); plot(alpha,e(6,:,i),'co-');
xlabel('Alpha'),ylabel('Power'),title('Interaction Time x Class')
legend('Obs TreeFM','Obs TreeFM F','Nested Obs TreeFM F','BoxCox Obs TreeFM F', 'Rank Obs TreeFM F', 'Raw + Rank Obs TreeFM F')

%% Compare time

i=2;
figure, plot(alpha,e(1,:,i),'b'); hold on, plot(alpha,e(2,:,i),'r'); plot(alpha,e(3,:,i),'g.-'); plot(alpha,e(4,:,i),'k--'); plot(alpha,e(5,:,i),'m:'); plot(alpha,e(6,:,i),'co-');  
xlabel('Alpha'),ylabel('Power'),title('Factor Time')
 legend('Obs TreeFM','Obs TreeFM F','Nested Obs TreeFM F','BoxCox Obs TreeFM F', 'Rank Obs TreeFM F', 'Raw + Rank Obs TreeFM F')

%% Compare class

i=3;
figure, plot(alpha,e(1,:,i),'b'); hold on, plot(alpha,e(2,:,i),'r'); plot(alpha,e(3,:,i),'g.-'); plot(alpha,e(4,:,i),'k--'); plot(alpha,e(5,:,i),'m:'); plot(alpha,e(6,:,i),'co-');  
xlabel('Alpha'),ylabel('Power'),title('Factor Class')
 legend('Obs TreeFM','Obs TreeFM F','Nested Obs TreeFM F','BoxCox Obs TreeFM F', 'Rank Obs TreeFM F', 'Raw + Rank Obs TreeFM F')

%% Compare individual

i=4;
figure, plot(alpha,e(1,:,i),'b'); hold on, plot(alpha,e(2,:,i),'r'); plot(alpha,e(3,:,i),'g.-'); plot(alpha,e(4,:,i),'k--'); plot(alpha,e(5,:,i),'m:'); plot(alpha,e(6,:,i),'co-'); 
xlabel('Alpha'),ylabel('Power'),title('Factor Individual')
 legend('Obs TreeFM','Obs TreeFM F','Nested Obs TreeFM F','BoxCox Obs TreeFM F', 'Rank Obs TreeFM F', 'Raw + Rank Obs TreeFM F')

