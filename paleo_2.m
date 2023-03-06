%% Individual event PDF computation when the younger earthquake date is unknown

%Create PDFs from the numerical dates and sample random values from them.

if oldest_unfaulted  < data_old
    oldest_unfaulted  = data_old;
    sd_oldest_unfaulted  = sd_data_old;
elseif oldest_unfaulted <max(x_allsites.Event_date_young) & x_allsites.Error_1 ~=0
    [oldest_unfaulted , ind3] = max(x_allsites.Event_Date_young);
    sd_oldest_unfaulted  = x_allsites.Error_1(ind3);
end

y1=pdf('norm',time, data_old, sd_data_old);
y2=pdf('norm', time, oldest_unfaulted, sd_oldest_unfaulted);
n1 = round(normrnd(data_old, sd_data_old, s, 1), 0);
n2 =round(normrnd(oldest_unfaulted , sd_oldest_unfaulted , s, 1), 0);

%Establish ranges between samples at each numerical date limiting an event
for j=1:s
    range = min([n1(j),n2(j)]):max([n1(j),n2(j)]);
    [val, pos] = intersect(time, range);
    h1(pos(1):pos(end), j) = 1;
end

%Sum and normalize the probability of all ranges
sumh = sum(h1, 2)./s; 
positive = find(sumh>0);