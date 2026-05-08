GC_replacing = num_MC + randperm(num_GC,num_rep_day);

for i=1:num_glom %GC2MC
    temp = num_ob;
    j=1;
    while temp>0
        if temp >= num_sis_per_lv
            indx = (i-1)*num_ob + (j-1)*num_sis_per_lv + 1:(i-1)*num_ob + j*num_sis_per_lv;
            temp = temp-num_sis_per_lv;
        else
            indx = (i-1)*num_ob + (j-1)*num_sis_per_lv + 1:i*num_ob;
            temp=0;
            j=num_lvs;
        end        
        S(sis_MCs(indx), 1+num_MC:end) = str_GCtoMC*sprand11(length(indx),num_GC,den_GCtoMC_days(j,day+2));
        j=j+1;
    end
end

S(GC_replacing,1:num_MC) = str_MCtoGC*sprand11(num_rep_day,num_MC,den_MCtoGC);%MC2GC
S(GC_replacing,1+num_MC:end) = str_GCtoGC*sprand11(num_rep_day,num_GC,den_GCtoGC);%GC2GC: new GCs as receivers
S(1+num_MC:end,GC_replacing) = str_GCtoGC*sprand11(num_GC,num_rep_day,den_GCtoGC);%GC2GC: new GCs as senders
S=S-diag(diag(S));

Sg(GC_replacing,:) = str_PCtoGC*sprand11(num_rep_day,num_pyramd,den_PCtoGC);%PC2GC