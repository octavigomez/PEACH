%% Individual event PDF computation when the older earthquake date is unknown

%Create PDFs from the numerical dates and sample random values from them.

% if historical_date<data_young
%     historical_date= data_young;
%     sd_historical = sd_data_young;
% end

y1=pdf('norm', time, oldest_faulted, sd_oldest_faulted);
y2=pdf('norm',time, data_young, sd_data_young);
n1 = round(normrnd(oldest_faulted, sd_oldest_faulted, s, 1), 0);
n2 = round(normrnd(data_young, sd_data_young, s, 1), 0);

%Establish ranges between samples at each numerical date limiting an event
for j=1:s
    range = min([n1(j), n2(j)]):max([n1(j), n2(j)]);
    [val, pos] = intersect(time, range);
    h1(pos(1):pos(end), j) = 1;
end

%Sum and normalize the probability of all ranges
sumh = sum(h1, 2)./s; 
positive = find(sumh>0);