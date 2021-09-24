clear all;
%FIND PATH
path0 = strsplit(pwd, '/');
path0 = strjoin(path0(1:end-1),'/');
%%ADD REQUIRED SUBFLODERS
path2=[path0 '/Functions'];
addpath(genpath(path2))
path3 = strsplit(path0, '/');
path3 = strjoin(path3(1:end-2),'/');
path3=[path3 '/FastICA'];
addpath(path3)
%%ADD REQUIRED SUBFLODERS
% path5=[ '~/Resting2/Serotonin/shared'];
% addpath(path5)
% path6=[ '~/Resting2/Serotonin/knn'];
% addpath(path6)
% path7=[ '~/Resting2/Serotonin/sqdistance'];
% addpath(path7)
% path8=[ '~/Resting2/Serotonin/Non_parametric_Entropy_Estimator'];
% addpath(path8)

co = HShannon_kNN_k_initialization(1);

load metadata.mat;

control=1:121;
control([ageMatchedCtlForSZ ageMatchedCtlForBD])=[];
ageMatchedCtlForAD=control;

load UCLA_all.mat;
N=68;
Isubdiag = find(tril(ones(N),-1));

SCHIZOS=122:171;
BIPOLAR=172:220;
ADHD=221:242; %260;

% %%%%%%%%%%%%%%

nfig=1;
for NA=1:4
    for NB=NA+1:4
        if NA==1 && NB==2
            Group1 = ageMatchedCtlForSZ;
            n_subjects_g1=length(Group1);
            Group2=SCHIZOS;
            n_subjects_g2=length(Group2);
        elseif NA==1 && NB==3
            Group1=ageMatchedCtlForBD;
            n_subjects_g1=length(Group1);
            Group2=BIPOLAR;
            n_subjects_g2=length(Group2);
        elseif NA==1 && NB==4
            Group1=ageMatchedCtlForAD;
            n_subjects_g1=length(Group1);
            Group2=ADHD;
            n_subjects_g2=length(Group2);
        elseif NA==2 && NB==3
            Group1=SCHIZOS;
            n_subjects_g1=length(Group1);
            Group2=BIPOLAR;
            n_subjects_g2=length(Group2);
        elseif NA==2 && NB==4
            Group1=SCHIZOS(1:length(ageMatchedCtlForAD));
            n_subjects_g1=length(Group1);
            Group2=ADHD;
            n_subjects_g2=length(Group2);
        elseif NA==3 && NB==4
            Group1=BIPOLAR(1:length(ageMatchedCtlForAD));
            n_subjects_g1=length(Group1);
            Group2=ADHD;
            n_subjects_g2=length(Group2);
        end
        
        nsub1=1;
        for nsub=Group1
            signaldata=(squeeze(time_series(:,:,nsub)));
            tc_aal{nsub1,1}=signaldata;
            nsub1=nsub1+1;
        end
        nsub1=1;
        for nsub=Group2
            signaldata=(squeeze(time_series(:,:,nsub)));
            tc_aal{nsub1,2}=signaldata;
            nsub1=nsub1+1;
        end
        
        n_Task=2;
        [N_areas, Tmax]=size(tc_aal{1,1});
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% 1 - Compute the ICs
        disp('Processing the ICs from BOLD data')
        
        Tmaxtotal=Tmax*(n_subjects_g1+n_subjects_g2);
        
        % Parameters of the data
        TR=2.;  % Repetition Time (seconds)
        
        % Preallocate variables to save FC patterns and associated information
        Time_all=zeros(2, Tmaxtotal); % vector with subject nr and task at each t
        t_all=0; % Index of time (starts at 0 and will be updated until n_Sub*Tmax)
        
        % Bandpass filter settings
        fnq=1/(2*TR);                 % Nyquist frequency
        flp = 0.04;                    % lowpass frequency of filter (Hz)
        fhi = 0.07;                    % highpass
        Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
        k=2;                          % 2nd order butterworth filter
        [bfilt,afilt]=butter(k,Wn);   % construct the filter
        clear fnq flp fhi Wn k
        
        n_Subjects(1)=n_subjects_g1;
        n_Subjects(2)=n_subjects_g2;
        
        timeseriedatatotal=zeros(Tmaxtotal,N);%length(Isubdiag)); % All leading eigenvectors
        
        kk3=1;
        timeseriebold=[];
        for task=1:n_Task
            for s=1:n_Subjects(task)
                % Get the BOLD signals from this subject in this task
                [N_areas, Tmax]=size(tc_aal{s,task});
                % Get the BOLD signals from this subject in this task
                BOLD = tc_aal{s,task};
                BOLD=BOLD(:,1:Tmax);
                Phase_BOLD=zeros(N_areas,Tmax);
                
                for seed=1:N_areas
                    BOLD(seed,:)=BOLD(seed,:)-mean(BOLD(seed,:));
                    signal_filt =filtfilt(bfilt,afilt,BOLD(seed,:));
                    timeseriedata(seed,:) = zscore(signal_filt);
                    Phase_BOLD(seed,:) = angle(hilbert(signal_filt));
                end
                
                for t=1:Tmax
                    iPH=zeros(N_areas,N_areas);
                    %Calculate the Instantaneous FC (BOLD Phase Synchrony)
                    for n=1:N_areas
                        for m=1:N_areas
                            iPH(n,m)=cos(Phase_BOLD(n,t)-Phase_BOLD(m,t));
                        end
                    end
                    t_all=t_all+1; % Update time
                    
                    [V1,~]=eigs(iPH,1);
                    timeseriedatatotal(t_all,:)=V1; %iPH(Isubdiag);
                    Time_all(:,t_all)=[s task]; % Information that at t_all, V1 corresponds to subject s in a given task
                end
                
                timeseriebold=horzcat(timeseriebold,timeseriedata);
                
            end
        end
        
        timeseriedatatotal=timeseriedatatotal';
        NumAssemblies_ph = NumberofIC(timeseriedatatotal);
        NumAssemblies = NumberofIC(timeseriebold);
        
        [Projection, ICcomp, W]=fastica25(timeseriebold,'numOfIC',NumAssemblies,'verbose','off');
        [Projection_ph, ICcomp_ph, Wph]=fastica25(timeseriedatatotal,'numOfIC',NumAssemblies_ph,'verbose','off');
        
        %%%%%%%
        
        for task=1:n_Task
            for ass=1:N_areas
                lista=[];
                for s=1:n_Subjects(task)
                    % Select the time points representing this subject and task
                    T=((Time_all(1,:)==s)+(Time_all(2,:)==task))>1;
                    lista=horzcat(lista,squeeze(timeseriebold(ass,T)));
                    entro(task,s,ass)=HShannon_kNN_k_estimation(squeeze(timeseriebold(ass,T)),co);
                    meanvalue(task,s,ass)=mean(squeeze(timeseriebold(ass,T)));
                end
                distribution{task,ass}=lista;
            end
        end
        
        for task=1:n_Task
            for ass=1:NumAssemblies
                lista=[];
                for s=1:n_Subjects(task)
                    % Select the time points representing this subject and task
                    T=((Time_all(1,:)==s)+(Time_all(2,:)==task))>1;
                    lista=horzcat(lista,squeeze(Projection(ass,T)));
                    entroIC(task,s,ass)=HShannon_kNN_k_estimation(squeeze(Projection(ass,T)),co);
                    meanvalueIC(task,s,ass)=mean(squeeze(Projection(ass,T)));
                end
                distributionIC{task,ass}=lista;
            end
        end
        
        for task=1:n_Task
            for ass=1:NumAssemblies_ph
                lista=[];
                for s=1:n_Subjects(task)
                    % Select the time points representing this subject and task
                    T=((Time_all(1,:)==s)+(Time_all(2,:)==task))>1;
                    lista=horzcat(lista,squeeze(Projection_ph(ass,T)));
                    entroICph(task,s,ass)=HShannon_kNN_k_estimation(squeeze(Projection_ph(ass,T)),co);
                    meanvalueICph(task,s,ass)=mean(squeeze(Projection_ph(ass,T)));
                end
                distributionICph{task,ass}=lista;
            end
        end
        
        %%%%%%%%%%%%
        sigvalue=0.01;
        
        figure(nfig);
        hold on;
        nfig=nfig+1;
        
        for ass=1:N_areas
            a=distribution{1,ass};
            b=distribution{2,ass};
            [aux pvalbold(ass)]=kstest2(a,b);
        end
        
        [h]=FDR_benjHoch(pvalbold,sigvalue);
        
        Boldsig=mean(timeseriebold(h,:),2);
        entroBoldsig=entro(:,:,h);
        meanvalueBoldsig=meanvalue(:,:,h);
        
        hh=zeros(N_areas,1);
        for ass=1:N_areas
            if ismember(ass,h)
                hh(ass)=0.9;
            else
                hh(ass)=0;
            end
        end
        
        subplot(3,1,1);
        plot(1:N_areas,pvalbold);
        hold
        for ass=1:N_areas
            if hh(ass)>0
                plot(ass,hh(ass),'*k');
            end
        end
        
        
        for ass=1:NumAssemblies
            a=distributionIC{1,ass};
            b=distributionIC{2,ass};
            [aux pvalIC(ass)]=kstest2(a,b);
        end
        
        [h]=FDR_benjHoch(pvalIC,sigvalue);
        
        ICcompsig=ICcomp(:,h);
        entroICcompsig=entroIC(:,:,h);
        meanvalueICcompsig=meanvalueIC(:,:,h);
        
        hh=zeros(NumAssemblies,1);
        for ass=1:NumAssemblies
            if ismember(ass,h)
                hh(ass)=0.9;
            else
                hh(ass)=0;
            end
        end
        
        subplot(3,1,2);
        plot(1:NumAssemblies,pvalIC);
        hold on;
        for ass=1:NumAssemblies
            if hh(ass)>0
                plot(ass,hh(ass),'*k');
            end
        end
        
        
        for ass=1:NumAssemblies_ph
            a=distributionICph{1,ass};
            b=distributionICph{2,ass};
            [aux pvalICph(ass)]=kstest2(a,b);
        end
        
        [h]=FDR_benjHoch(pvalICph,sigvalue);
        ICcompsig_ph=ICcomp_ph(:,h);
        entroICcompsig_ph=entroICph(:,:,h);
        meanvalueICcompsig_ph=meanvalueICph(:,:,h);
        
        numcomp=size(ICcompsig_ph,2);
        numind=size(ICcomp_ph,1);
        ICcompsig_ph_matrix=[];
        for ass=1:numcomp
            ICcompsig_ph_matrix(:,:,ass)=ICcompsig_ph(:,ass)*(ICcompsig_ph(:,ass))';
%             for i=1:numind
%                 [ii jj]=ind2sub([N_areas N_areas],Isubdiag(i));
%                 ICcompsig_ph_matrix(ii,jj,ass)=ICcompsig_ph(i,ass);
%                 ICcompsig_ph_matrix(jj,ii,ass)=ICcompsig_ph(i,ass);
%             end
        end
        
        hh=zeros(NumAssemblies_ph,1);
        for ass=1:NumAssemblies_ph
            if ismember(ass,h)
                hh(ass)=0.9;
            else
                hh_ph(ass)=0;
            end
        end
        
        subplot(3,1,3);
        plot(1:NumAssemblies_ph,pvalICph);
        hold on;
        for ass=1:NumAssemblies_ph
            if hh(ass)>0
                plot(ass,hh(ass),'*k');
            end
        end
        
        
        if NA==1 && NB==2
            save UCLAnets_Controls_Schizos.mat meanvalueBoldsig meanvalueICcompsig meanvalueICcompsig_ph entroBoldsig entroICcompsig entroICcompsig_ph Boldsig ICcompsig_ph_matrix ICcompsig ICcompsig_ph;
        elseif NA==1 && NB==3
            save UCLAnets_Controls_Bipolar.mat meanvalueBoldsig meanvalueICcompsig meanvalueICcompsig_ph entroBoldsig entroICcompsig entroICcompsig_ph Boldsig ICcompsig_ph_matrix ICcompsig ICcompsig_ph;
        elseif NA==1 && NB==4
            save UCLAnets_Controls_ADHD.mat meanvalueBoldsig meanvalueICcompsig meanvalueICcompsig_ph entroBoldsig entroICcompsig entroICcompsig_ph Boldsig ICcompsig_ph_matrix ICcompsig ICcompsig_ph;
        elseif NA==2 && NB==3
            save UCLAnets_Schizos_Bipolar.mat meanvalueBoldsig meanvalueICcompsig meanvalueICcompsig_ph entroBoldsig entroICcompsig entroICcompsig_ph Boldsig ICcompsig_ph_matrix ICcompsig ICcompsig_ph;
        elseif NA==2 && NB==4
            save UCLAnets_Schizos_ADHD.mat meanvalueBoldsig meanvalueICcompsig meanvalueICcompsig_ph entroBoldsig entroICcompsig entroICcompsig_ph Boldsig ICcompsig_ph_matrix ICcompsig ICcompsig_ph;
        elseif NA==3 && NB==4
            save UCLAnets_Bipolar_ADHD.mat meanvalueBoldsig meanvalueICcompsig meanvalueICcompsig_ph entroBoldsig entroICcompsig entroICcompsig_ph Boldsig ICcompsig_ph_matrix ICcompsig ICcompsig_ph;
        end
        
        clear distribution distributionIC distributionICph;
        clear pvalbold pvalIC pvalICph;
        clear NumAssemblies NumAssemblies_ph;
        clear Boldsig ICcompsig_ph_matrix ICcompsig ICcompsig_ph;
        clear h hh;
        
    end
end
