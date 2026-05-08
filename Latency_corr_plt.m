Latency_glo1= Latency_glo;
Latency_glo1(Latency_glo1==0)=1000;
result1 = 1./Latency_glo1;
result1(result1 == 0.001) =0;
rho_lat = corr(result1);
rho_lat(end+1,: )=-1;
figure;
imagesc(rho_lat);
colormap(redblue)
colorbar;
ylim([0.5, 100.5]);
xlim([0.5, 100.5]);
xlabel('Odors');
ylabel('Odors');
title('Latency correlation');
%%
rho_lat(end,: )=[];
B1 = triu(rho_lat,1);
C1=[];
for i =1:size(B1,2)
    C1= [C1,B1(i,i+1:end)];
end
% C1 = B1(B1~=0);
figure(6);
histogram(C1,15,'Normalization','probability');
set(gca, 'YScale', 'log')
title('Latency histogram');
xlabel('Correlation coefficient');
ylabel('Probability');