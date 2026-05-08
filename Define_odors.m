clear;clc;

num_glom = 50;
num_ordors = 300;
Latency_glo = zeros(num_glom,num_ordors);
Stim_glo = zeros(num_glom,num_ordors);
num_aglo = 4+randi(6,num_ordors,1);
odor_sti = 8.5 - 0.5+ 1*rand(num_ordors,1);

for odor = 1:num_ordors
    Glom1 = randperm(num_glom,num_aglo(odor));
    lat = [40:160/(num_aglo(odor)-1):200]+20*rand(1,num_aglo(odor)); 
    Latency_glo(Glom1,odor) = lat';
    largest_glom = Glom1(1);
    Stim_glo(largest_glom,odor) =  odor_sti(odor);
    Stim_glo(setdiff(Glom1,largest_glom),odor) = odor_sti(odor)*0.85+0.25*randn(num_aglo(odor)-1,1);
end
filename = strcat('Odor_denfinition\Neurogensis_newOdor_table_',num2str(num_ordors),'new.mat');
save(filename,'Latency_glo',"Stim_glo","num_ordors",'num_aglo','odor_sti');
%% heatmap plot
figure(1);
imagesc(Latency_glo);
colorbar
ice =cmocean('ice');
ice(1,:)=0;
ice=flip(ice);
colormap(tempo)
xlabel('Odors');
ylabel('Glomerulus');

figure(2);
imagesc(Stim_glo);
colorbar
tempo = cmocean('tempo');
tempo(end,:)=0;
colormap(tempo)
xlabel('Odors');
ylabel('Glomerulus');

%% Latency correlation
Latency_corr_plt;

%% Stimulus magnitude correlation
Magnitude_corr_plt;

%% Stimulus map
odor = randi(num_ordors)
Stimulus_map_plot;