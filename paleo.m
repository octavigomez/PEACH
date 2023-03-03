%% Individual event PDF computation when both age constraints are known

%Create PDFs from the numerical dates and sample random values from them.

y1=pdf('norm',time, data_old, sd_data_old);
y2=pdf('norm', time, data_young, sd_data_young);
n1 = round(normrnd(data_old, sd_data_old, s, 1), 0);
n2 =round(normrnd(data_young, sd_data_young, s, 1), 0);

%Establish ranges between samples at each numerical date limiting an event
for j=1:s                                               
    range = min([n1(j),n2(j)]):max([n1(j),n2(j)]);           
    [val, pos] = intersect(time, range);                
    h1(pos(1):pos(end), j) = 1;
end

%Sum and normalize the probability of all ranges
sumh = sum(h1, 2)./s; 
positive = find(sumh>0);