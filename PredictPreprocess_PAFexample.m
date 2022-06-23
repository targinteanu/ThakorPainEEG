%% Initial set up of variables
clear all 
%eeglab
home = '/Users/torenarginteanu/Documents/MATLAB/Thakor Pain EEG/Data_Chronic Pain';
cd(home)
addpath /data/andrew/ali_eeg/eeg_PAF_noICA/FieldTrip/fieldtrip-20121231/;
addpath /data/andrew/ali_eeg/eeg_PAF_noICA/FieldTrip/fieldtrip-20121231/utilities/;
subjects = dir('sub-*');

%% Initial Preprocessing in EEGLAB
for s = 1:size(subjects,1)
    cd(subjects(s,1).name);
    cd eeg
    clearvars datasets
    datasets = dir('*.vhdr');
    for d = 1:size(datasets,1);
        clearvars EEG rejected
        EEG = pop_fileio(datasets(d,1).name);
        EEG = pop_resample( EEG, 500);
        EEG=  pop_chanedit(EEG, 'lookup','/usr/local/eeglab13_6_5b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp');
        EEG = pop_select( EEG,'nochannel',{'GSR' 'HR' 'RESP'});
        EEG = pop_reref( EEG, []);
        EEG = pop_eegfiltnew(EEG, 2, 100, 826, 0, [], 1);
        EEG.event(1,end).type = 'End';
        EEG = pop_rmdat( EEG, {'Start/End'},[0 300] ,0);
        
         %inspect raw spectra and reject spectra that are odd looking
            
        figure; pop_spectopo(EEG, 1, [0  15000], 'EEG' , 'freqrange',[2 100],'electrodes','off');
        pop_eegplot( EEG, 1, 1, 1);
        rejected = input('Please identify sensor # to remove: ');
        if size(rejected,1) ~=0;
            EEG = pop_select( EEG,'nochannel',rejected);
            EEG = pop_reref( EEG, []); % KD - reref again after removing noisy channels
        end
        temp_file = [datasets(d,1).name(1:end-5) '_preprocessedAJF.mat'];
        save(temp_file,'EEG','rejected');
    end
    cd(home)
end

%% Now we go to FieldTrip
clearvars -except home subjects
for s = 1%:size(subjects,1)
    cd(subjects(s,1).name);
    cd eeg
    clearvars datasets
    datasets = dir('*preprocessedAJF.mat');
    for d = 1:size(datasets,1);
        clearvars data cfg X rejected data_pruned data_freq paf power
        load(datasets(d,1).name,'EEG');
        data=eeglab2fieldtrip(EEG,'preprocessing');
        data.label={EEG.chanlocs.labels}; %%%%% fix glitch that labels aren't imported
        cfg.length=5; %determine length of segment
        cfg.overlap= 0;
        data=ft_redefinetrial(cfg,data);
        %PCA to remove non-brain artifacts
        cfg.method       = 'runica';
        cfg.runica.pca=15;
        X=ft_componentanalysis(cfg,data);
        cfg.component=[1:15];
        cfg.layout='EEG1010.lay';
        ft_topoplotIC(cfg,X)
        ft_databrowser(cfg,X)
        prompt = 'Please Identify Component to reject in matrix form';
        rejected = input(prompt);
        prompt = 'Please Identify Component to reject in matrix form';
        rejected = input(prompt);
        cfg.component=rejected; % rejected components
        data_pruned= ft_rejectcomponent(cfg,X);
        %transform time
        cfg=[];
        cfg.method='mtmfft';
        cfg.taper   =  'hanning';
        cfg.foi =[2:.20:50];
        cfg.keeptrials='yes';
        data_freq=ft_freqanalysis(cfg,data_pruned);
        
        freq = 9:.2:11;
        paf = [];
        power = [];
        for z = 1:size(data_freq.powspctrm,1); % KD trial
            for y = 1:size(data_freq.powspctrm,2); % KD sensor
                temp = squeeze(data_freq.powspctrm(z,y,:));
                temp2 = zscore(temp);
                paf(y,z)= sum(freq'.*(((temp(36:46)))))/sum((((temp(36:46))))); %PAF
                power(y,z) = sum(temp(31:51));
            end
        end
        temp_file = [datasets(d,1).name(1:end-4) '_fftAJF.mat'];
        save(temp_file,'data','X','rejected','data_pruned','data_freq','paf','power');
    end
    cd(home)
end

