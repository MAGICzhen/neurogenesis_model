clear;clc;
%% 
feedback=1;
dt=1;%[ms]
load('Weight_matrix_New\Neurogensis_newWeightMatrix_day0.mat');
load('Odor_denfinition\Neurogensis_newOdor_table_100new.mat');
save_folder = strcat('UnGroup_new_run10\');
num_ordors = 100;
Act_duration = 90;
t_bin=20;% window time bin:[ms]
trials=10;% number of trails for each odor
sniff_cycle=250;
time_odor=600;%odor onset time [ms]
num_sniff=1;
tspan=time_odor+num_sniff*sniff_cycle;%[ms]
Latency_glo(Latency_glo==0)=1000;
sniff_length = (tspan-time_odor-t_bin+2)/dt ;
S_base = S;
Sg_base = Sg;
mkdir(save_folder);
%% post-genesis connectivity
num_rep_day = round(num_GC*0.1);%number of GCs updated every day: 30%
num_days = 10;
GiantMatrixPC_mean =zeros(num_pyramd, num_ordors*(num_days+2)*sniff_length);%trial-averaged
GiantMatrixMC_mean =zeros(num_MC, num_ordors*(num_days+2)*sniff_length);%trial-averaged
MC_PC_Group_divide;
Chang_den_GC2MC;
para_name = strcat(save_folder,'parameters.mat');
save(para_name,'-v7.3');
%%
for day = -1:num_days
    if day>0
        GC_replace_New;
    end
%     Si = S(1:num_MC,1+num_MC:end);
%     sum(Si,2);
%     figure(1); hold on
%     plot(ans);
    GC_Group_divide;
    Simulate_network_responses;
    disp("Day: " + num2str(day));
    fname = strcat(save_folder,"Matrix_day_",num2str(day),'.mat');
    save(fname,'GiantMatrixMC_day','GiantMatrixPC_day','S','Sg','W','Counts_OB','Counts_Piri','-v7.3');
end

save(strcat(save_folder,'Day10_all_mean.mat'),'GiantMatrixMC_mean','GiantMatrixPC_mean','-v7.3');

