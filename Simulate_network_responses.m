GiantMatrixMC_day = zeros(num_MC, num_ordors*sniff_length, trials);
GiantMatrixPC_day = zeros(num_pyramd, num_ordors*sniff_length, trials);
for odor=1:num_ordors
    firing_Pyram=[];
    firing_MC=[];
    firing_GC=[];
    Counts_OB=zeros(num_MC+num_GC,tspan/dt,trials);%% spike counts over all trails
    Counts_Piri=zeros(num_pyramd+num_FBIN+num_FFIN,tspan/dt,trials);
    for t=1:trials
        v_OB=zeros(num_MC+num_GC,tspan/dt);
        u_OB=zeros(num_MC+num_GC,tspan/dt);
        v_OB(:,1)=-65*ones(num_MC+num_GC,1);
        u_OB(:,1)=b1.*v_OB(:,1);
        
        v_PC=zeros(num_pyramd+num_FFIN+num_FBIN,tspan/dt);
        u_PC=zeros(num_pyramd+num_FFIN+num_FBIN,tspan/dt);
        v_PC(:,1)=-65*ones(num_pyramd+num_FFIN+num_FBIN,1);
        u_PC(:,1)=b2.*v_PC(:,1);
        
        OSNinput_times = [];
        for i=1:num_glom
            OSNinput_times=[OSNinput_times;time_odor+Latency_glo(i,odor)+round((-3+6*rand(num_ob,1)))];%the latency casued delay of stimulus for each MC
        end
        OSNinput_times=repmat(OSNinput_times,1,num_sniff);
        OSNinput_times=OSNinput_times+repmat(0:sniff_cycle:(num_sniff-1)*sniff_cycle,num_MC,1);
        %
        firing_OB=[];
        I_OSN_inputs=zeros(num_MC+num_GC,tspan/dt);
        I_MC_ex=zeros(num_MC+num_GC,tspan/dt);
        I_GC_in=zeros(num_MC+num_GC,tspan/dt);
        I_spn=zeros(num_MC+num_GC,tspan/dt);
        I_spnPCx=zeros(num_pyramd+num_FFIN+num_FBIN,tspan/dt);
        I_PC_ex=zeros(num_MC+num_GC,tspan/dt);
        
        firedPClen=zeros(1,tspan/dt);
        
        firing_piriform=[];
        I_ex_OB=zeros(num_pyramd+num_FFIN+num_FBIN,tspan/dt);
        I_piriform_ex=zeros(num_pyramd+num_FFIN+num_FBIN,tspan/dt);
        I_piriform_in=zeros(num_pyramd+num_FFIN+num_FBIN,tspan/dt);
        I_recurrent=zeros(num_pyramd,tspan/dt);
        I_FFI=zeros(num_pyramd,tspan/dt);
        I_FBI=zeros(num_pyramd,tspan/dt);
        
        % % % odor-evoked stimuli, as a glomerulus-correlated step current to each MC
        % % % this part decay with time
        I_dir = Stim_glo(:,odor);
        I_dir = repmat(I_dir,1,num_ob)';
        I_dir=I_dir(:)+0.5*randn(num_MC,1);
        I_corrsti=zeros(num_MC,(sniff_cycle)/dt);%within one sniff cycle
        seedcorr=0.5*randn(num_glom,(sniff_cycle)/dt);%within one sniff cycle
        alpha=0.8;%correlation degree of the stimulus
        sigma=1;
        for i=1:num_MC
            I_corrsti(i,:) = I_dir(i,1)+ alpha*sigma*seedcorr(floor((i-1)/num_ob)+1,:)+(1-alpha)*sigma*randn(1,(sniff_cycle)/dt);
        end
        I_corrsti=I_corrsti.*(repmat([0:dt:sniff_cycle-dt],num_MC,1)>=repmat(OSNinput_times(:,1)-time_odor,1,(sniff_cycle)/dt) & repmat([0:dt:sniff_cycle-dt],num_MC,1)<=repmat(OSNinput_times(:,1)-time_odor+Act_duration,1,(sniff_cycle)/dt)); % turned on only after the glomeruli activated

        I_corrsti= [zeros(num_MC,time_odor/dt),repmat(I_corrsti,1,num_sniff)];
        
        for i=2:tspan/dt
            fired1 = find(v_OB(:,i-1)>=v_thrshld_OB); % MCs and GCs that spike, no duplicate
            Counts_OB(fired1,i-1,t) = Counts_OB(fired1,i-1,t)+1;% spike counts at each time point over all trials
            firing_OB = [firing_OB; (i-2)*dt+0*fired1,fired1];
            v_OB(fired1,i-1) = v_reset_OB(fired1);% reset
            u_OB(fired1,i-1) = u_OB(fired1,i-1)+d1(fired1);% reset
            
            fired2 = find(v_PC(:,i-1)>=v_thrshld_PC); % piriform cells that spike, no duplicate
            Counts_Piri(fired2,i-1,t) = Counts_Piri(fired2,i-1,t)+1;% spike counts at each time point over all trials
            firing_piriform = [firing_piriform; (i-2)*dt+0*fired2,fired2];
            v_PC(fired2,i-1) =v_reset_PC(fired2);% reset
            u_PC(fired2,i-1) = u_PC(fired2,i-1)+d2(fired2);% reset
            %1.75, 0.8
%             if i*dt >= time_odor+min(Latency_glo(:,odor))-0 %after
%                 spn_col = 1.0*randn(num_MC,1).*spn_act_MC{odor} + 1.4*randn(num_MC,1).*spn_Nonact_MC{odor};
%                 I_spn(:,i)=I_spn(:,i-1)*exp(-dt/3)+[spn_col;0.8*randn(num_GC,1)];
%             else
%                 I_spn(:,i)=I_spn(:,i-1)*exp(-dt/3)+[1.5*randn(num_MC,1);0.8*randn(num_GC,1)];
%             end
            if feedback ==1
                if i*dt >= time_odor+min(Latency_glo(:,odor))-0 %after
                    spn_col = 1.0*randn(num_MC,1).*spn_act_MC{odor} + 1.35*randn(num_MC,1).*spn_Nonact_MC{odor};
                    I_spn(:,i)=I_spn(:,i-1)*exp(-dt/3)+[spn_col;0.8*randn(num_GC,1)];
                else
                    I_spn(:,i)=I_spn(:,i-1)*exp(-dt/10)+[0.65*randn(num_MC,1)+1.0*abs(sum(S(1:num_MC,fired1(fired1>num_MC)),2));0.8*randn(num_GC,1)];
                end
            else
                I_spn(:,i)=I_spn(:,i-1)*exp(-dt/3)+[1.5*randn(num_MC,1);0.8*randn(num_GC,1)];
            end
            %0.9,0.45
            if feedback ==1
                if i*dt <= time_odor+min(Latency_glo(:,odor))+5 %before
                    I_spnPCx(:,i) = I_spnPCx(:,i-1)*exp(-dt/3)+[1.0*randn(num_pyramd,1); 0.8*randn(num_FFIN,1);zeros(num_FBIN,1)];
                else
                    I_spnPCx(:,i) = I_spnPCx(:,i-1)*exp(-dt/3)+[1.0*randn(num_pyramd,1); 0.8*randn(num_FFIN,1);zeros(num_FBIN,1)];
                end
            else
                I_spnPCx(:,i) = I_spnPCx(:,i-1)*exp(-dt/3)+[1.0*randn(num_pyramd,1); 1.0*randn(num_FFIN,1);0*randn(num_FBIN,1)];
            end
            I_sti = I_corrsti(:,i).*exp(-((i-1)*dt-OSNinput_times(:,1+max(floor(((i-1)*dt-time_odor)/sniff_cycle),0)))/500);% exponential decay of step current stimulis from receptors
            I_OSN_inputs(:,i) = [I_sti;zeros(num_GC,1)];%OSN inputs
            
            fired_MC = fired1(fired1<=num_MC);
            fired_MCa = intersect(fired_MC,Act_MC{odor} );%firing MCs that are also activated by odor
            fired_MCb = setdiff(fired_MC,fired_MCa,'stable');
            fired_PC = fired2(fired2<=num_pyramd);
            fired_PCa = intersect(fired_PC,Act_PC{odor} );
            fired_PCb = setdiff(fired_PC,fired_PCa,'stable');
            
            if i*dt >= time_odor+min(Latency_glo(:,odor))-5
                I_ex_OB(1:num_pyramd,i)=I_ex_OB(1:num_pyramd,i-1)*exp(-dt/tau_ex2)+sum(W(1:num_pyramd,fired_MC),2);
                I_ex_OB(1+num_pyramd:end,i)=I_ex_OB(1+num_pyramd:end,i-1)*exp(-dt/tau_ex2)+sum(W(1+num_pyramd:end,fired_MC),2);
            else
                I_ex_OB(1:num_pyramd,i)=I_ex_OB(1:num_pyramd,i-1)*exp(-dt/tau_ex2)+0.8*sum(W(1:num_pyramd,fired_MC),2);
                I_ex_OB(1+num_pyramd:end,i)=I_ex_OB(1+num_pyramd:end,i-1)*exp(-dt/tau_ex2)+0.8*sum(W(1+num_pyramd:end,fired_MC),2);
            end

            I_piriform_ex(1:num_pyramd,i)=I_piriform_ex(1:num_pyramd,i-1)*exp(-dt/tau_ex2)+sum(W(1:num_pyramd,num_MC+fired2(fired2<=num_pyramd)),2);
            I_piriform_ex(num_pyramd+num_FFIN+1:end,i) = I_piriform_ex(num_pyramd+num_FFIN+1:end,i-1)*exp(-dt/tau_ex2)+sum(W(num_pyramd+num_FFIN+1:end,num_MC+fired2(fired2<=num_pyramd)),2);
             
            I_FFI(:,i)=I_FFI(:,i-1)*exp(-dt/tau_in2)+sum(W(1:num_pyramd,num_MC+fired2(fired2>num_pyramd & fired2<=num_pyramd+num_FFIN)),2);
            I_FBI(:,i)=I_FBI(:,i-1)*exp(-dt/tau_in2)+sum(W(1:num_pyramd,num_MC+fired2(fired2>num_pyramd+num_FFIN)),2);
            if feedback ==0 && i*dt>= time_odor+max(mink(Latency_glo(:,odor),2))+10
                ffin_temp = sum(50*W(:,num_MC+fired2(fired2>num_pyramd & fired2<=num_pyramd+num_FFIN)),2);
                fbin_temp = sum(W(:,num_MC+fired2(fired2>num_pyramd+num_FFIN)),2);
                I_piriform_in(:,i)=I_piriform_in(:,i-1)*exp(-dt/tau_in2)+ffin_temp+fbin_temp;
            else
                ffin_temp = sum(W(:,num_MC+fired2(fired2>num_pyramd & fired2<=num_pyramd+num_FFIN)),2);
                fbin_temp = sum(W(:,num_MC+fired2(fired2>num_pyramd+num_FFIN)),2);
                I_piriform_in(:,i)=I_piriform_in(:,i-1)*exp(-dt/tau_in2)+ffin_temp+fbin_temp;
%                 I_piriform_in(:,i)=I_piriform_in(:,i-1)*exp(-dt/tau_in2)+sum(W(:,num_MC+fired2(fired2>num_pyramd)),2);
            end
            v_PC(:,i)=max(v_min,v_PC(:,i-1)+dt*(0.04*v_PC(:,i-1).^2+5*v_PC(:,i-1)+140-u_PC(:,i-1)) + dt*(I_spnPCx(:,i)+I_ex_OB(:,i)+I_piriform_ex(:,i)+I_piriform_in(:,i)));
            u_PC(:,i)=u_PC(:,i-1)+dt*a2.*(b2.*v_PC(:,i-1)-u_PC(:,i-1));
            v_PC(fired2,i)=v_reset_PC(fired2);% refractoriness for 1ms
            
            if ~isempty(fired_MC)%% long latency of granule cells since they receive the excitatory inputs from MCs
                glom=floor((fired_MC-1)/num_ob)+1;%the glomeruli that fired_MCs belong to
                delay_time=i+round(Latency_GC_glom(:,unique(glom)));%onset time for these inhibitions of GCs
                EPSP=full(S(1+num_MC:num_MC+num_GC,fired_MC));% with reptitions of the same glomerulus
                if length(glom)>1
                    k = find([true;diff(glom(:))~=0;true]);
                    r = [k(1:end-1) diff(k)];
                    same_glom = r(r(:,2)>1,:);
                    del=[];
                    for n=1:size(same_glom,1)
                        EPSP(:,same_glom(n,1))=sum(EPSP(:,same_glom(n,1):same_glom(n,1)+same_glom(n,2)-1),2);
                        del=[del,same_glom(n,1)+1:same_glom(n,1)+same_glom(n,2)-1];
                    end
                    EPSP(:,del)=[];
                    [vals, col_idx] = sort(delay_time,2);
                    idx = bsxfun(@plus,(col_idx-1)*size(delay_time,1), (1:size(delay_time,1))');% idx is the one-dimensional index of elements of vals in A
                    is_duplicate=false(size(delay_time,1)*size(delay_time,2),1);%fix the size
                    is_duplicate(idx(:,2:end)) = vals(:,1:end-1) == vals(:,2:end);% duplicate latencies occure only from the second column
                    while ismember(1,is_duplicate)
                        is_duplicate = reshape(is_duplicate, size(delay_time));
                        delay_time(is_duplicate)=delay_time(is_duplicate)+1;
                        [vals, col_idx] = sort(delay_time,2);
                        idx = bsxfun(@plus,(col_idx-1)*size(delay_time,1), (1:size(delay_time,1))');
                        is_duplicate(idx(:,2:end)) = vals(:,1:end-1) == vals(:,2:end);
                    end
                end
                EPSP(delay_time>tspan-dt)=0;
                delay_time(delay_time>tspan-dt)=tspan-dt;
                delay_EPSP_MCtoGC=zeros(tspan/dt,num_GC);
                k=repmat(1:num_GC,1,size(delay_time,2));
                delay_EPSP_MCtoGC((k'-1)*tspan/dt+delay_time(:))=EPSP(:);
                delay_EPSP_MCtoGC=delay_EPSP_MCtoGC';
                I_MC_ex(num_MC+1:num_MC+num_GC,:)=I_MC_ex(num_MC+1:num_MC+num_GC,:)+delay_EPSP_MCtoGC;%% the decay of the previous instant fo time, plus the delay input of EPSP
            end
            I_MC_ex(:,i)=I_MC_ex(:,i-1)*exp(-dt/tau_ex1)+[sum(S(1:num_MC,fired_MC),2);I_MC_ex(num_MC+1:end,i)];
            
            if feedback==0
                I_GC_in(:,i)=I_GC_in(:,i-1)*exp(-dt/5)+sum(S(:,fired1(fired1>num_MC)),2);
                I_GC_in(fired_MCa,i)=I_GC_in(fired_MCa,i)+0.1*sum(S(fired_MCa,fired1(fired1>num_MC)),2);
            else
                if i*dt <= min(Latency_glo(Act_glo{odor},odor))+time_odor
                    I_GC_in(1:num_MC,i)=I_GC_in(1:num_MC,i-1)*exp(-dt/tau_in1)+sum(S(1:num_MC,fired1(fired1>num_MC)),2);
                else
                    I_GC_in(Act_MC{odor},i)=I_GC_in(Act_MC{odor},i-1)*exp(-dt/tau_in1)+sum(S(Act_MC{odor},fired1(fired1>num_MC)),2);
                    aa = fired1(fired1>num_MC);
                    if length(aa)>45
                        aa = aa(randperm(length(aa),45));
                    end
                    I_GC_in(setdiff(1:num_MC,Act_MC{odor}),i)=I_GC_in(setdiff(1:num_MC,Act_MC{odor}),i-1)*exp(-dt/tau_in1)+sum(S(setdiff(1:num_MC,Act_MC{odor}),aa),2);
                end
%                  I_GC_in(num_MC+1:num_MC+num_GC,i)=I_GC_in(num_MC+1:num_MC+num_GC,i)*exp(-dt/tau_in1)+sum(S(num_MC+1:num_MC+num_GC ,fired1(fired1>num_MC)),2);
                I_GC_in(Act_GC{odor} ,i)=I_GC_in(Act_GC{odor} ,i)*exp(-dt/tau_in1)+10*sum(S(Act_GC{odor} ,fired1(fired1>num_MC)),2);%linked to GC_divide
                I_GC_in(setdiff(num_MC+1:num_MC+num_GC,Act_GC{odor} ,'stable'),i)=I_GC_in(setdiff(num_MC+1:num_MC+num_GC,Act_GC{odor} ,'stable'),i)*exp(-dt/tau_in1)+sum(S(setdiff(num_MC+1:num_MC+num_GC,Act_GC{odor} ,'stable'),fired1(fired1>num_MC)),2);
            end
            
            firedPClen(i)=length(fired2(fired2<=num_pyramd));
            
            if feedback==1
                I_feedback = 0.65*randn(num_MC+num_GC,1); % stochastic feedback to GCs before the earliest PC fires
                I_feedback(1:num_MC) =0;                
                if firedPClen(i) > 2
                    temp=fired2(fired2<=num_pyramd);
                    temp=temp(randperm(firedPClen(i),2));
                    I_PC_ex(:,i)=I_PC_ex(:,i-1)*exp(-dt/10)+sum(Sg(:,temp),2) + I_feedback;
                else
                    I_PC_ex(:,i)=I_PC_ex(:,i-1)*exp(-dt/10)+sum(Sg(:,fired2(fired2<=num_pyramd)),2)+ I_feedback;
                end
            end
            v_OB(:,i)=max(v_min,v_OB(:,i-1)+dt*(0.04*v_OB(:,i-1).^2+5*v_OB(:,i-1)+140-u_OB(:,i-1)) + dt*(I_spn(:,i)+I_MC_ex(:,i)+I_GC_in(:,i)+I_OSN_inputs(:,i)+I_PC_ex(:,i)));
            u_OB(:,i)=u_OB(:,i-1)+dt*a1.*(b1.*v_OB(:,i-1)-u_OB(:,i-1));
        end
        % the last instant
        fired1=find(v_OB(:,i)>=v_thrshld_OB); % MCs and GCs that spike, no duplicate
        Counts_OB(fired1,i,t)=Counts_OB(fired1,i,t)+1;% spike counts at each time point over all trials
        firing_OB=[firing_OB; (i-1)*dt+0*fired1,fired1];
        v_OB(fired1,i)=v_reset_OB(fired1);% reset
        u_OB(fired1,i) = u_OB(fired1,i)+d1(fired1);% reset
        fired2=find(v_PC(:,i)>=v_thrshld_PC); % piriform cells that spike, no duplicate
        Counts_Piri(fired2,i,t)=Counts_Piri(fired2,i,t)+1;% spike counts at each time point over all trials
        firing_piriform=[firing_piriform; (i-1)*dt+0*fired2,fired2];
        v_PC(fired2,i)=v_reset_PC(fired2);% reset
        u_PC(fired2,i) = u_PC(fired2,i)+d2(fired2);% reset
        
    end
%     Counts_Piri = mean(Counts_Piri,3);
    FR_PC=zeros(num_pyramd, sniff_length,trials);% [Hz],FR of PC cells  
    for t=1:trials
        for j=time_odor:(tspan-t_bin+1)/dt
            FR_PC(:,j-time_odor+1,t)=sum(Counts_Piri(1:num_pyramd,j:j+t_bin-1,t),2)/(t_bin/1000);%
        end
    end

%     Counts_OB = mean(Counts_OB,3);
    FR_MC=zeros(num_MC, sniff_length,trials);% [Hz],FR of MC cells  
    for t=1:trials
        for j=time_odor:(tspan-t_bin+1)/dt
            FR_MC(:,j-time_odor+1,t)=sum(Counts_OB(1:num_MC,j:j+t_bin-1,t),2)/(t_bin/1000);%
        end
    end
%     Single_trials_TEST;                   
    GiantMatrixMC_day(:,1+(odor-1)*sniff_length : odor*sniff_length,:) =  FR_MC;
    GiantMatrixPC_day(:,1+(odor-1)*sniff_length : odor*sniff_length,:) =  FR_PC;
    GiantMatrixPC_mean(: , [1+(odor-1)*sniff_length : odor*sniff_length] +(day+2-1)* num_ordors*sniff_length) = mean(FR_PC,3);  
    GiantMatrixMC_mean(: , [1+(odor-1)*sniff_length : odor*sniff_length] +(day+2-1)* num_ordors*sniff_length) = mean(FR_MC,3);
end