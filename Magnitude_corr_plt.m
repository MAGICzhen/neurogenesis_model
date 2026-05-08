rho_amp = corr(Stim_glo);
rho_amp(end+1,1)=-1;
figure(3);
imagesc(rho_amp);
colormap(redblue)
colorbar;
ylim([0.5, 30.5])
xlim([0.5, 30.5])
xlabel('Odors');
ylabel('Odors');
title('Amplitude correlation');

rho_amp(end,:)=[];
B = triu(rho_amp,1);
% C = B(B~=0);
C=[];
for i =1:size(B,2)
    C= [C,B(i,i+1:end)];
end
figure(4);
% [f,xi] = ksdensity(C);
% plot(xi,f)
histogram(C,15,'Normalization','probability');
set(gca, 'YScale', 'log')
title('Amplitude histogram');
xlabel('Correlation coefficient');
ylabel('Probability');