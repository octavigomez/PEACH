%% Sensitivity analysis for to determine seed value for PEACH calculation

%Authors: Octavi Gómez, Bruno Pace, Francesco Visini, Joanna Faure Walker, Oona Scotti
%Year: 2023

%FUNDING
%This work has been supported by two consecutive postdoctoral grants awarded to 
% Octavi Gómez-Novell in 2022: “Borsa di studio n. 2341/2021” funded by the INGEO 
% Department (Università degli Studi “G. d%Annunzio” di Chieti e Pescara) and “Margarita Salas
% grant” awarded by the University of Barcelona and funded by the Spanish Ministry of Universities 
% with the EU “Next Generation” and “Plan de recuperación, transformación y resiliencia” programs.

%LICENSING
%The code and derivatives are published under the Creative Commons license CC-BY-NC 4.0. 
% For more info: https://creativecommons.org/licenses/by-nc/4.0/deed.es This means that you 
% are free to copy, share and edit the material as long as you give credit to the authors, publish 
% it disclosing the same license and always for non-comercial purposes. 

%Start

clear all
close all

%% Create PDFs from the numerical dates and sample random values from them.

pathin = 'Inputs';
fault = "faultR.csv";  % Change file to the desired fault paleoseismic dataset
data= readtable(fullfile(pathin, fault));
time = min(round(data.Event_date_old,0)-round(data.Error,0).*3):2022;
sz = 100:100:10000; % Seeds and intervals to explore.

col = 0;
for s=sz;
    col = col+1;
    for x = 1:size(data, 1)
        data_old = data.Event_date_old(x);
        sd_data_old = data.Error(x);
        data_young = data.Event_date_young(x);
        sd_data_young = data.Error_1(x);
        
        
        y1=pdf('norm',time, data_old, sd_data_old);
        y2=pdf('norm', time, data_young, sd_data_young);
        n1 = round(normrnd(data_old, sd_data_old, s, 1), 0);
        n2 =round(normrnd(data_young, sd_data_young, s, 1), 0);
        
        for j=1:s                                               
            range = min([n1(j),n2(j)]):max([n1(j),n2(j)]);           
            [val, pos] = intersect(time, range);                
            h1(pos(1):pos(end), j) = 1;
        end
        
        sumh = sum(h1, 2)./s; 
        positive = find(sumh>0);
        sensit(1:length(sumh), x, col) = sumh;
    end
end

%% Analyze differences in mean probabilities an plot the results

for len = 1:size(sensit, 2)
    means (len, :) = mean(sensit(:,len,:));
end
means_final = mean(means,1);

figure (1)
ay.FontSize = 12;
ax.FontSize = 12;
plot(sz, means_final, "-", "Color", "black", "LineWidth",1)
subtitle("Sensitivity analysis of the " + fault)
xlabel ("Seed (n)", "FontSize", 12)
ylabel("Average PDF probabilities", "FontSize", 12) %closer to 1 means more detailed sampling
xlim ([min(sz), max(sz)])
