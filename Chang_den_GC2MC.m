%GC inhibition changes over days
den_GCtoMC_days = zeros(length(den_GCtoMC),num_days+2);
den_GCtoMC_days(:,1) = den_GCtoMC';%day -1
tar_den = flip(den_GCtoMC);
for i=1:length(den_GCtoMC)
    den_GCtoMC_days(i,2:end) = linspace(den_GCtoMC(i),tar_den(i),num_days+1);
end
