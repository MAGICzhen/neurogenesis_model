%Grouping of MCs according to glomerulus
num_aglo = num_aglo';
Act_glo=cell(1,num_ordors);
Act_MC=cell(1,num_ordors);
Act_MCa=cell(1,num_ordors);
Act_MCb=cell(1,num_ordors);
spn_act_MC = cell(1,num_ordors);
spn_Nonact_MC = cell(1,num_ordors);
Earliest_glo=[];
for odor=1: num_ordors
    Act_glo{odor} = find(Latency_glo(:,odor)<sniff_cycle);
    tempt = repmat((Act_glo{odor}'-1)*num_ob,num_ob,1)+repmat([1:num_ob]',1,num_aglo(odor));
    Act_MC{odor} = tempt(:);
    Earliest_glo(1,odor)=find(Latency_glo(:,odor) ==min(Latency_glo(:,odor) ));
    Act_MCa{odor}=1+(Earliest_glo(1,odor)-1)*num_ob:Earliest_glo(1,odor)*num_ob; %earliest MCs
    Act_MCb{odor}=setdiff(Act_MC{odor}, Act_MCa{odor},'stable');
    spn_act_MC{odor} = ismember([1:num_MC]',Act_MC{odor});
    spn_Nonact_MC{odor} = ismember([1:num_MC]',setdiff([1:num_MC]',Act_MC{odor},'stable'));
end


%Not grouping of PCc but simply summarizing the random feeddforward
%projections: which MCs connect to which PCs
Act_PC=cell(1,num_ordors);%among those PCs that have feedforward weights, Act_PC have weights larger than 70% percentile
Act_PCa=cell(1,num_ordors);%among Act_PC, Act_PCs are driven by Act_MCa
W_feedforward = W(1:num_pyramd,1:num_MC);

for odor=1:num_ordors
    temp = sum(W_feedforward(:,Act_MC{odor}'),2);
    thr = (max(temp)-mean(temp))*0.35 + mean(temp);
    Act_PC{odor} = Pyramcells(temp>thr);%larger than 90% percentile
    temp1 = sum(W_feedforward(:,Act_MCa{odor}'),2);
    thr1 = (max(temp1)-mean(temp1))*0.25 + mean(temp1);
    Act_PCa{odor} = intersect(Act_PC{odor}, Pyramcells(temp1>thr1),'stable');
end