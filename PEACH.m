%% PEACH. Earthquake chronology modelling from paleoseismic data

%Authors: Octavi Gómez, Bruno Pace, Francesco Visini
%Year: 2022

%Start

clear all
close all

%% Path and read input data from table

set(0, 'defaultFigureRenderer', 'painters')
set(groot,'DefaultFigureGraphicsSmoothing','off')

pathin = 'Inputs';
x_allsites = readtable(fullfile(pathin, "Example.xlsx"));
site_specs = readtable (fullfile(pathin,"site_specs.txt"));

%% Read site specifications from input file and create respective variables
trunc = site_specs.sigma_level;
if trunc == 1
    trunc_level = [15.9 84.1];
elseif trunc == 2
    trunc_level = [2.3 97.7];
elseif trunc == 3
    trunc_level = [0.1 99.9];
elseif trunc == 0
    trunc_level = [0 100];
end

oldest_faulted= site_specs.oldest_faulted;
sd_oldest_faulteold = site_specs.sd_faulted;
s = site_specs.seed;


%% Calculate event PDFs based on the paleoseismic site variable

%The iteration is done site by site
[sites, IA, IC] = unique(x_allsites.Site, "stable");
for nsite = 1:size(sites)
u = strcmp(x_allsites.Site, sites(nsite));
x = x_allsites(u,:);
youngest_age=2022;
time = min(round(x_allsites.Event_date_old,0)-round(x_allsites.Error,0).*3):youngest_age;
column = 0;

if isnan(site_specs.oldest_faulted)
     [oldest_faulted, ind] = min(x.Event_date_old);
     sd_oldest_faulted = x.Error(ind);
     if isnan(oldest_faulted)
        [oldest_faulted, ind2] = min(x_allsites.Event_date_old);
        sd_oldest_faulted = x_allsites.Error(ind2);
   end
end

%Sites formed by all NaNs are removed from the dataset

 remove = find(isnan(x.Event_date_old) & isnan(x.Event_Date_young));
 if length(remove) == size(x,1)
 x(remove,:) =[];
 end

%Loop to extract the age data from their respective variables referenced to each site.
for ii = 1:size(x)-1
    if isnan(x.Event_Date_young(ii)) & isnan(x.Event_Date_young(ii+1)) & x.Event_date_old(ii) ~= x.Event_date_old(ii+1)
        x.Event_Date_young(ii+1) = x.Event_date_old(ii);
        x.Error_1(ii+1) = x.Error(ii);
    elseif isnan(x.Event_date_old(ii)) & isnan(x.Event_date_old(ii+1)) & x.Event_Date_young(ii) ~= x.Event_Date_young(ii+1)
        x.Event_date_old(ii) = x.Event_Date_young(ii+1);
        x.Error(ii) = x.Error_1(ii+1);
    end
end

for i=1:size(x)
    h1=[];
    h1(1:length(time),1) = 0; 
    data_old = round(x.Event_date_old(i),0);
    sd_data_old = round(x.Error(i),0);
    data_young = round(x.Event_Date_young(i),0);
    sd_data_young = round(x.Error_1(i),0);
    oldest_unfaulted = site_specs.oldest_unfaulted ;
    sd_oldest_unfaulted  = site_specs.sd_unfaulted;

% Conditions depending on the constraints on the input data.
if  ~isnan(data_old) & ~isnan(sd_data_old) & ~isnan(data_young) & ~isnan(sd_data_young)
    run paleo
elseif ~isnan(data_old) & ~isnan(sd_data_old) & isnan(data_young) & isnan(sd_data_young)
    run paleo_2
elseif isnan(data_old) & isnan(sd_data_old) & ~isnan(data_young) & ~isnan(sd_data_young)
    run paleo_3
end

%Store all event PDFs in a 3D variable.
column = column+1;
sumh = sumh/max(sumh);
outputs(:, column, nsite) = sumh;

%% Plot all event PDFs together

figure (100)
hold on
for i2 = nsite(1):nsite(end)   
    plot (time, outputs(:,:,nsite), "Color",[0.8 0.8 0.8]);
end
title("Mean PDF for the whole fault trace", "Interpreter","none", "FontSize",12, "Position", ...
    [median(time), 1.03]);
xlabel("Years (BCE/CE)");
ylabel("Normalized probability");
box on
end
end

%% The event PDFs are truncated based on the sigma truncation levels set by the user.
for kr=1:size(outputs, 3);
    for ar = 1:size(outputs, 2);
        condi = length(find(outputs(:, ar, kr)>0))>1;
        if condi ==1;
            cum_out = (cumsum(outputs(:,ar,kr))./sum(outputs(:,ar,kr)))*100;
            low_limit = (find(cum_out<=trunc_level(1),1,"last"))-1;
            low_limit(low_limit==0) = 1;
            upper_limit=(find(cum_out>=trunc_level(2),1, "first"))+1;
            upper_limit(upper_limit>size(outputs,1)) = size(outputs,1);
            outputs(1:low_limit, ar, kr) = 0;
            outputs(upper_limit:end, ar, kr) = 0;
        else
            %Do nothing.
        end
    end
end

%% Compute the mean probability curve of all sites in a fault section

%Maximum probabilities for each time and site are extracted and mean probability curve is computed.
for g = 1:size(outputs, 3);
    temp(:,g) = max(outputs(:,:,g),[], 2);
end
temp2=mean(temp, 2);
tempnorm = (temp2-min(temp2))./(max(temp2)-min(temp2));
plot (time, tempnorm, "-", "LineWidth",1, "Color","black");
hold off

%% Find peaks in the mean probability curve

%Using the findpeaks function from MATLAB in both directions (left to right and right to left).
[peaks1, locs1, widths1, proms1] = findpeaks(tempnorm, time, "WidthReference","halfprom", ...
    "MinPeakProminence",0.03); 
pks = findpeaks(tempnorm, time, "Annotate","extents", "WidthReference","halfprom", ...
    "MinPeakProminence",0.03);
flipped_tempnrom =flipud(tempnorm);
flipped_time = fliplr(time);
[peaks2, locs2, widths2, proms2] = findpeaks(flipped_tempnrom, time, "WidthReference", ...
    "halfprom", "MinPeakProminence",0.03); 

for len=1:length(locs2)
    ind_locs2(len) = find(time==locs2(len));
    locs2(len) = flipped_time(ind_locs2(len));
end
locs =[locs1;fliplr(locs2)];
locs_plateau = mean(locs(:,:));

%Plot the peak detection in the mean curve
figure (900)
hold on
scatter(locs_plateau, peaks1,"black","filled","v")
plot(time, tempnorm, "Color", "black")
title("Peak identification in the mean PDF", "Interpreter","none", "FontSize",12,"position", ...
    [median(time), 1.03])
xlabel("Years (BCE/CE)")
ylabel("Normalized probability")
legend ("Peaks", "Location", "best")
box on
hold off

%% Compute final PDFs from the peaks and event PDFs

%Extraction of all event PDFs intersecting the peak positions.
check3=[];
for l =1:length(peaks1)
trial = round(locs_plateau(l));
pos1 = find(trial==time);
new_out = [];
    for k=1:size(outputs, 3);
        check = outputs(pos1, :, k);
        new_out = [new_out, outputs(:,find(check>0),k)];

        %Compute mean and standard deviation for each event PDF.
        for a = 1:size(outputs, 2);
            weights_indiv_pdf = outputs(:,a,k);
            mean_indiv_pdf(a,:) = (time*weights_indiv_pdf)./sum(weights_indiv_pdf);
            std_indiv_pdf(a,:) = std(time,weights_indiv_pdf);
            indiv_range(a,:,k) = [mean_indiv_pdf(a,:)-std_indiv_pdf(a,:), mean_indiv_pdf(a,:)+std_indiv_pdf(a,:)];
            mean_range(a,:,k)= [mean_indiv_pdf(a,:), indiv_range(a,:,k)];
            if isnan(mean_range(a,:,k));
                mean_range(a,:,k)=0;
            end
        end
        mean_range(mean_range==0) = NaN;
        mean_range(:,:,k) =sort(mean_range(:,:,k),1,"descend");
        row(:,:,k)= ~isnan(mean_range(:,1, k));
        nonans=find(row(:,1,k)==1);
        rownan(:,:,k)= isnan(mean_range(:,1, k));
        nans = find(rownan(:,1,k)==1);
        mean_range(:,:,k)= [mean_range(nonans,:,k);mean_range(nans,:,k)];

        %Analyze effective overlaps between event PDFs by using their mean and standard deviation bars.
        for a2 = 1:size(outputs,2)-1
            cond(a2,1,k) = mean_range(a2,1,k) >= mean_range(a2+1, 1, k) & mean_range(a2, 1, k)<= mean_range(a2+1, 3, k) & mean_range(a2,2,k)<=mean_range(a2+1,1,k);
            if cond(a2, 1, k) == 1
                  cond2(a2:a2+1, 1, k) =1;
            elseif cond (a2, 1, k) == 0
                  cond2(a2+1, 1, k)=0;
            end
        end

        %Distinguish different grups of overlaps
        idx_ovlp = find(cond2(:, :, k));
        check2= outputs(pos1, idx_ovlp, k);
        check2(check2~=0) = 1;
        no_zero_cols{l,k}= idx_ovlp'.*check2;
        no_zero_cols = cellfun(@(z) z(z~=0), no_zero_cols, "UniformOutput",false);

        %Find the number of PDFs in each site that have effective overlaps.
        tempnum_events(l,k)= sum(check2>0);
        if tempnum_events(l,k) == 0;
            tempnum_events(l,k) =1;
        end
    end

    % Compute product PDFs per event by multiplying each event PDF intersected at each peak.
    final(:, l) = prod(new_out, 2);
    finalnorm(:, l) = normalize(final(:, l), "range");
end


%Remove identical final PDFs coming from noisy peaks close in time.
[final2, indexes]= unique(final', "rows", "stable");
final2= final2';
finalnorm3 = normalize(final2, "range");

%Compute total number of PDFs that form each final PDF

[minnum_events, indx_max] = max(tempnum_events,[], 2);
Min_num_events = minnum_events(indexes);
Num_events_site = [Min_num_events,indx_max(indexes)];
[Num_events_single, indx_rep] = unique(Min_num_events, "stable");

% Effective_breaks: correspond to the breaks in sites with largest number of overlaps
breaks_str=[];
for br = 1:size(no_zero_cols,1)
breaks{br,1} = no_zero_cols{br, indx_max(br)};
breaks_str= char(breaks_str,num2str(breaks{br}));
end

breaks_str(1,:) =[];
breaks_str=breaks_str(indexes,:);
[id, IA1, IC2] = unique(breaks_str, "rows","stable");

%Separate groups of PDFs with the same number of PDFs contained
Min_num_events = [Min_num_events, IC2];
Min_num_events(end+1,:) = 0;
Min_num_events=[zeros(1, size(Min_num_events,2));Min_num_events];
Diff_consec = find(diff(Min_num_events(:,2))~=0);
[rx, b] = size(Diff_consec);
Start = Min_num_events(2,1);
Groups = cell(rx,1);

for rs = 1:rx
    End= Diff_consec(rs,1);
    Groups{rs}=Min_num_events(Start:End,1);
    Start =End+1;
end
Min_num_events(1,:) = [];
Min_num_events(end,:) =[];
Groups(1,:) = [];
n = cellfun(@(v)v(1),Groups);

%Create all combinations of events that add up to the number of events contained in a PDF
for group_ev= 1:length(n)
    count_events = n(group_ev);
    len=size(Groups{group_ev},1);
    if len>1
         if count_events>len
            Comb = cell(1,len);
            [Comb{1:len}] = ndgrid(1:len);
            Comb = allVL1(len, count_events);                                                                                                %Function created by Bruno Luong --> Bruno Luong (2022). All Permutations of integers with sum criteria (https://www.mathworks.com/matlabcentral/fileexchange/17818-all-permutations-of-integers-with-sum-criteria), MATLAB Central File Exchange. Retrieved November 25, 2022.
            Comb(any(Comb==0,2),:) =[];
         elseif count_events<=len
            Comb=ones(1,len);
         end
    else
       Comb=  count_events;
    end
    Peak_comb{group_ev}= Comb;
end

for opt = 1:size(Peak_comb, 2)
    Comb_idx= Peak_comb{opt};
    for opt2 = 1:size(Comb_idx, 1)
            Comb_final{opt,opt2}= Comb_idx(opt2,:);
    end
    Comb_sum_opt(opt,1)= sum(Comb_final{opt,1},2);
    Comb_final(cellfun(@isempty, Comb_final))= {NaN};
end

for opt3=1:size(Comb_final, 1);
    opt4=size(Comb_final,2);
     if isnan(Comb_final{opt3,opt4});
        Comb_final(opt3,2:end) = Comb_final(opt3,1);
     else
     end
end

%Total number of events in the fault
Comb_sum = sum(Comb_sum_opt,1);

%% Analyze whether the total number of events is higher than the most populated site.

most_populated = max(size(outputs, 2));
if most_populated >Comb_sum
    warning("There is a problem: the total number of final events is smaller than the " + ...
        "number of events in your most populated site. Check user manual for guidance.");
else
    %do nothing
end

%% Create the final output

finalnorm4=[time',finalnorm3];
for r=1:size(finalnorm4,1)
    if sum(finalnorm4(r,2:end))==0;
        finalnorm4(r,1)=0;
    end
end

% Correction of the lengths of the variables when the PDFs are truncated.
if trunc_level~=[0, 100]
    temp = sum(finalnorm4(:,2:end),2);
    index2 = find(temp(:,1)>0,1,"first");
    index3 = find(temp(:,1)>0,1,"last");
    time_range = [finalnorm4(index2,1)-1:finalnorm4(index3,1)+1];
    new_pdfs=finalnorm4(index2-1:index3+1, 2:end);
    final_pdfs=[time_range',new_pdfs]';
    Set = (1:length(final_pdfs(2:end,1)))';
elseif trunc_level  == [0, 100];
    temp = sum(finalnorm4(:,2:end),2);
    index2 = find(temp(:,1)>0,1,"first");
    index3 = find(temp(:,1)>0,1,"last");
    time_range = [finalnorm4(index2,1):finalnorm4(index3,1)];
    new_pdfs=finalnorm4(index2:index3, 2:end);
    final_pdfs=[time_range',new_pdfs]';
    Set = (1:length(final_pdfs(2:end,1)))';
end

%Generate the final output
for combo=1:size(Comb_final,2);
    Event_count_combos(combo,1:numel(cell2mat(Comb_final(:,combo)'))) = cell2mat(Comb_final(:,combo)');
end
Event_count_combos  = Event_count_combos';
final_outputs = [Set,Event_count_combos,new_pdfs'];

%% Plot final results (PDFs)

final_plot = final_outputs(:,size(Event_count_combos,2)+2:end);

% Plot final PDFs with boxplot indicating the statistics of each final PDF.
figure (4500);
box on
tiledlayout(2,1, "TileSpacing", "tight")
step3=0;
nexttile
hold on
non_zeros_stored =[];
for pdf_stats = 1:size(final_plot, 1);
    step3= step3+1;
    weights= final_plot(pdf_stats,:);
    means (pdf_stats,:)= (time_range*weights')./sum(weights);
    sigma_1(pdf_stats,:) = std(time_range, weights);
    output_stats(pdf_stats,:)= [means(pdf_stats,:), sigma_1(pdf_stats,:)];
    non_zeros = find(final_plot(pdf_stats,:)>0);
    non_zeros_stored{pdf_stats} = non_zeros;
    if length(time_range(non_zeros))>1;
        for_box = randpdf(final_plot(pdf_stats,non_zeros), time_range(non_zeros), ...
            [1, length(time_range(non_zeros))]);                                                                                            %Function created by Adam Niselony (Opole University of Technology, Poland) --> Adam Nieslony (2022). Random numbers from a user defined distribution (https://www.mathworks.com/matlabcentral/fileexchange/26003-random-numbers-from-a-user-defined-distribution), MATLAB Central File Exchange. Retrieved October, 2022.
        [fillhandle, msg] = jbfill(time_range(non_zeros), final_plot(pdf_stats, non_zeros)*0.7+ ...
        step3, final_plot(pdf_stats, non_zeros)*0+step3, [0.8, 0.8, 0.8], [0 0 0], 0,1);                      %Function created by John Bockstege. --> John Bockstege (2022). Shade area between two curves (https://www.mathworks.com/matlabcentral/fileexchange/13188-shade-area-between-two-curves), MATLAB Central File Exchange. Retrieved October, 2022.
        quantiles_final(pdf_stats,:) = quantile(for_box,[0.25, 0.5, 0.75]);
        mean_final(pdf_stats,:) = mean(for_box);
        devi_final(pdf_stats,:) = std(for_box);
        scatter (median(for_box), step3, "red", "filled", "diamond", SizeData=20);
    elseif length(time_range(non_zeros))==1;
        quantiles_final(pdf_stats,:) = quantile(time_range(non_zeros),3);
        mean_final(pdf_stats,:) = time_range(non_zeros);
        devi_final(pdf_stats,:) = 0;
        text(mean_final(pdf_stats,:), step3, "  Hist. ("+ mean_final(pdf_stats,:) + ")  ", ...
            "HorizontalAlignment","right")
    end
    scatter (means(pdf_stats,:), step3, "black", "filled", "o", SizeData=10);
end

yticks([1:step3]);
ytick = gca;
ytick.YGrid = "on";
yticklabels("PDF. "+ Set);
ylabel("Paleoearthquake chronology");
title("Final paleoarthquake chronology for the whole fault and comparison with site data", "FontSize",12, "Position", ...
    [median(time_range), step3+1+0.2]);
ylim ([0.5, step3+1]);
xlim([min(time_range), max(time_range)]);
xticklabels(" ")
box on
hold off
final_pdf_stats = [Set,mean_final, devi_final, quantiles_final, Event_count_combos];

%Plot the final PDFs together with the event PDFs in each site to validate the model.
nexttile
hold on
box on
for positive_pdfs=1:size(final_plot,1)
    positives1 = find(final_plot(positive_pdfs,:), 1, "first");
    positives2 = find(final_plot(positive_pdfs,:), 1, "last");
    if trunc_level ~= [0 100];
        positives1 = positives1-1;
        positives2 = positives2+1;
    end
    positives = positives1:positives2;
    hleg(positive_pdfs) =  plot(time_range(positives), ...
        (final_plot(positive_pdfs, positives)*length(IA)+1), "-",  "LineWidth",0.5, "Color",[0.7 0.1 0.2]);
end
rg=0;
for rr =1:length(IA)
    step = 0;
    ra = find(IC==rr);
    for xplot = length(ra):-1:1
        rg = rg+1;
        step = step+1;
        outputs_positive_1 = find(outputs(:,xplot,rr)>0,1, "first");
        outputs_positive_2 = find(outputs(:,xplot,rr)>0,1, "last");
        outputs_positive = outputs_positive_1:outputs_positive_2;
        if length(outputs_positive)>1
            if trunc_level ~= [0 100]
                outputs_positive_1 = outputs_positive_1-1;
                outputs_positive_2=outputs_positive_2+1;
            end 
           indiv_pdfs = plot([time(outputs_positive_1), time(outputs_positive_2)], ...
               [rr+step*(1/(length(ra)+1)),rr+step*(1/(length(ra)+1))], "Color", [0.6 0.6 0.6], ...
               "LineWidth",1, "LineStyle","-");
        elseif length(outputs_positive)==1
            indiv_pdfs2 = scatter(time(outputs_positive), rr+step*(1/(length(ra)+1)), ...
                "black", "filled", "diamond", "SizeData",10);
        end
    end
end
max_point =rr+1;
label_sites = string(sites)+"   " ;
set(gca,'TickLabelInterpreter', 'none');
yticks([1:rr]);
text(zeros(1,rr)+min(time_range), [1.6:rr+0.8], label_sites, "Interpreter","none", ...
     "HorizontalAlignment","right", "Color",[0 0 0], "Rotation",45)
yticklabels(" ")
ytickangle(0);
ytick = gca;
ytick.YRuler.TickLabelGapOffset = 5;
ytick.YGrid = "on";
xlabel ("Years (BCE/CE)")
xlim([min(time_range), max(time_range)]);
ylim([1, max_point]);
labels = num2str(final_outputs(:, 2:size(Event_count_combos,2)+1));
hold off
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.5 0.8]);

%% Ending: export outputs to folder

header_time = min(time_range):max(time_range);
contained_pdfs=1:size(Event_count_combos,2);
VarNames = ["PDF ID", "N. events hyp. "+ string(contained_pdfs), string(header_time)];
VarNames2 = ["PDF ID","Mean", "1σ", "Q25", "Q50", "Q75", "N. events hyp. "+contained_pdfs];
export_table=[VarNames; final_outputs];
stats_table = array2table(final_pdf_stats, 'VariableNames', VarNames2);

%Table.
writematrix(export_table, "Outputs/Final_PDFs.csv")
writetable(stats_table, "Outputs/Final_PDFs_stats.csv")

%Visualization figures.
saveas (figure(900), fullfile("Outputs/DetectedPeaks.eps"),"epsc")
saveas (figure(4500), fullfile("Outputs/FinalPDFs_statistics.eps"),"epsc")

%Summary
disp("Summary");
disp("- The total number of events identified in the fault is: "+ Comb_sum);
disp("- These events are contained in a total number of final PDFs of: " + size(final_plot, 1));
disp("For more information check out the outputs in the /Outputs folder.");

%End.
end_msg = msgbox("Your calculation is complete!", "Success", "help");