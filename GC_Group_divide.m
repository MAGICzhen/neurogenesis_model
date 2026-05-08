Act_GC=cell(1,num_ordors);
S_GC2MC = S(1:num_MC,1+num_MC:end);
S_MC2GC = S(1+num_MC:end,1:num_MC);
Granuels = num_MC+[1:num_GC];
for odor=1:num_ordors
    temp1 = abs(sum(S_GC2MC(Act_MC{odor},:),1));
    temp2 = sum(S_MC2GC(:,Act_MC{odor}'),2);
    gc1 = Granuels(temp1>prctile(temp1,85));%larger than 85% percentile
    gc2 = Granuels(temp2>prctile(temp2,85));%larger than 85% percentile
    Act_GC{odor} = intersect(gc1,gc2,'stable');
end