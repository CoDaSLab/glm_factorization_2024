%% Considerations for missing data, outliers and transformations in permutation testing for ANOVA  with multivariate responses
% Oliver Polushkina Merchanskaya, Michael D. Sorochan Armstrong, Carolina Gómez Llorente, Patricia Ferrer, Sergi Fernandez-Gonzalez, Miriam Perez-Cruz, María Dolores Gómez-Roig, José Camacho
%
% Power curves in Section 5.2. Figure 3.
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 21/Jul/2023
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

%% Normal

cd normal

load random_e e alpha

% Compare class

i=3;
figure, hold on, plot(alpha,e(2,:,i),'g'); plot(alpha,e(4,:,i),'k--'); plot(alpha,e(5,:,i),'b-.');  plot(alpha,e(6,:,i),'r:'); 
xlabel('Effect Size'),ylabel('Power'),title('Factor Class')
legend('Raw data','BoxCox', 'Rank', 'Raw + Rank')

cd ..

saveas(gcf,'./Figures/normal')
saveas(gcf,'./Figures/normal.eps','epsc')
saveas(gcf,'./Figures/normal.png','png')

%% Uniform

cd uniform

load random_e e alpha

% Compare class

i=3;
figure, hold on, plot(alpha,e(2,:,i),'g'); plot(alpha,e(4,:,i),'k--'); plot(alpha,e(5,:,i),'b-.');  plot(alpha,e(6,:,i),'r:'); 
xlabel('Effect Size'),ylabel('Power'),title('Factor Class')
legend('Raw data','BoxCox', 'Rank', 'Raw + Rank')

cd ..

saveas(gcf,'./Figures/uniform')
saveas(gcf,'./Figures/uniform.eps','epsc')
saveas(gcf,'./Figures/uniform.png','png')

%% Exp^3

cd exp3

load random_e e alpha

% Compare class

i=3;
figure, hold on, plot(alpha,e(2,:,i),'g'); plot(alpha,e(4,:,i),'k--'); plot(alpha,e(5,:,i),'b-.');  plot(alpha,e(6,:,i),'r:'); 
xlabel('Effect Size'),ylabel('Power'),title('Factor Class')
legend('Raw data','BoxCox', 'Rank', 'Raw + Rank')

cd ..

saveas(gcf,'./Figures/exp3')
saveas(gcf,'./Figures/exp3.eps','epsc')
saveas(gcf,'./Figures/exp3.png','png')

%% Normal + outlier

cd 'normal - outlier'

load random_e e alpha

% Compare class

i=3;
figure, hold on, plot(alpha,e(2,:,i),'g'); plot(alpha,e(4,:,i),'k--'); plot(alpha,e(5,:,i),'b-.');  plot(alpha,e(6,:,i),'r:'); 
xlabel('Effect Size'),ylabel('Power'),title('Factor Class')
legend('Raw data','BoxCox', 'Rank', 'Raw + Rank')

cd ..

saveas(gcf,'./Figures/normalOut')
saveas(gcf,'./Figures/normalOut.eps','epsc')
saveas(gcf,'./Figures/normalOut.png','png')

 

