% https://github.com/BrainDynamicsUSYD/braintrak/wiki
% Estimation of multiscale neurophysiologic parameters by electroencephalographic means
% Physiologically based arousal state estimation and dynamics

clear all;close all;
addpath('D:\My Files\Work\BGU\scripts\nft\braintrak');
addpath(genpath('D:\My Files\Work\BGU\scripts\nft\corticothalamic-model'));

debugmode = false; %to test chain 

% fdata = load_br_data('demo');
% load('demo.mat');
% P = s_single;

[fn, fp] = uigetfile(['D:\My Files\Work\BGU\Datasets\DOC\mat\model_fit_figures\' '*.mat'], 'Select data file');
load([fp fn]);
% load("D:\My Files\Work\BGU\Datasets\DOC\mat\model_fit_figures\2_HEALTHY_0eyeo_ica_125Hz_for_modelling_RobinsonNFT_eegPSDsim_1ch_reg_minimizeYagradTNC_Abeysuriyav1test1.mat");
f = eegPatient.f(2:end); %exclude 0Hz
P = eegPatient.Psd(2:end)';

%init model
full_gab_model = bt.model.full_gab([]);

%get typical gab and nus
[f_typical,P_typical,idx,db_data,typical_chisq] = db_fit.quick_fit(full_gab_model,f,P);
typical_gab = db_data.gab(idx,:);
typical_nus = db_data.nus(idx,:);

% %use python fitted params to init braintrak model
% f_fit = eegSim.f(2:end); %exclude 0Hz
% P_fit = eegSim.Psd(2:end)';
% python_params = model.params;
% python_params.gab = [optVals.Model_Gee optVals.Model_Gei optVals.Model_Ges optVals.Model_Gse...
%     optVals.Model_Gsr typical_gab(6) optVals.Model_Gre optVals.Model_Grs]; 
% python_params.gabcd = [optVals.Model_Gee  optVals.Model_Gei  optVals.Model_Ges*optVals.Model_Gse...
%     optVals.Model_Ges*optVals.Model_Gsr*optVals.Model_Gre  optVals.Model_Gsr*optVals.Model_Grs];
% python_params.alpha(:) = optVals.Model_alpha;
% python_params.beta(:) = optVals.Model_beta;
% python_params.t0 = optVals.Model_t0;
% python_params.taues = python_params.t0/2;
% python_params.tause = python_params.t0/2;
% python_params.emg_a = eps;
% python_initial_vals = [python_params.gab(1:5) python_params.gab(7:8) python_params.alpha(1)...
%     python_params.beta(1) python_params.t0 python_params.emg_a];
% full_gab_model = bt.model.full_gab(python_initial_vals);

%fit corticothalamic model
% fit_result =  bt.core.fit_spectrum(bt.model.full,f,P,[],[],[],[],[],debugmode); %same as bt.fit
fit_result =  bt.core.fit_spectrum(full_gab_model,f,P,[],[],'30s',[],[],debugmode);
% fit_result = bt.fit_track(bt.model.full,'demo',1,1,'10s');
% fit_result = bt.core.fit_spectrum(bt.model.dummy,1,1,[],[],1e6);%dummy model
P_fit = fit_result.fit_data.fitted_P;
f_fit = fit_result.fit_data.target_f;

%plot
figure;
loglog(f,P, f_fit,P_fit, f_typical,P_typical);xlim([1 45]);xlabel('Hz');title('fitted spectrum');
legend('experimental','fitted', 'typical');

% fit_result.plot()
% fit_result.head_plot(false)
% fit_result.interactive_fit
% figure;
% loglog(f,P, f_fit,P_fit);xlim([1 45]);xlabel('Hz');title('fitted spectrum');
% legend('experimental','fitted');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = p_from_fitted_params(full_gab_model,model.params,fit_result.fit_data.fitted_params);

% params = model.params; % = fit_result.model.p.copy(); %params.m

% %abeysuria 2015 phia&nus
% % phivals = [p.phia(1) p.phia(1) p.phia(2) p.phia(3)];
% params.phia = [5.248361515 15.3960197 8.789733431];
% % id_map = [1 2 3 1 2 3 7 8 4 5 6];
% params.nus = [0.001525377176 -0.003022754434 0.0005674779589 0.003447358203 -0.001465128967 0.003593330094 0.0001695899041 5.070036187e-05];

% %braintrak fitted params
% Gee = fit_result.fit_data.fitted_params(strcmp(fit_result.model.param_names,'Gee')); 
% Gei = fit_result.fit_data.fitted_params(strcmp(fit_result.model.param_names,'Gei')); 
% Gese = fit_result.fit_data.fitted_params(strcmp(fit_result.model.param_names,'Gese')); 
% Gesre = fit_result.fit_data.fitted_params(strcmp(fit_result.model.param_names,'Gesre')); 
% Gsrs = fit_result.fit_data.fitted_params(strcmp(fit_result.model.param_names,'Gsrs'));

% %find typical params by setting params.gabcd from db_data.nus and comparing to GABCD
% GABCD = [Gee Gei Gese Gesre Gsrs];
% best_err = inf;
% typical_nus = [];
% typical_gab = [];
% for nus_idx = 1:size(db_data.nus,1)
%     params = model.params;
%     params.nus = db_data.nus(nus_idx,:);
%     for gabcd_idx = 1:size(params.gabcd ,1)
%         err = (params.gabcd(gabcd_idx,:) - GABCD)*(params.gabcd(gabcd_idx,:) - GABCD)';
%         if err < best_err
%             best_err = err;
%             typical_nus = params.nus;
%             typical_gab = params.gab(gabcd_idx,:);
%         end
%     end
% end

% %set gab using typical_gab
% % Gee
% % Gei
% Ges = typical_gab(3);
% Gse = Gese/Ges;
% Gsn = typical_gab(6);
% Gre = typical_gab(7);
% Gsr = Gesre/Ges/Gre;
% Grs = Gsrs/Gsr;

% %set gab naively
% Ges = Gese; Gse = 1; Gsr = Gsrs; Grs = 1; Gre = Gesre/Ges/Gsr; Gsn = 1;

% %set phia (Qe Qi Qs Qr),  nus (nuee nuei nues nuse nusr nusn nure nurs), gabcd (Gee Gei Gese Gesre Gsrs])
% params.gab = [Gee Gei Ges Gse Gsr Gsn Gre Grs];
% % params.nus = typical_nus; %sets gabcd, gab and phia

% %unnecessary (executed automatically when p.gab is set
% if ~isempty(params.phia)
%     params.complete_nuab(1);
% else
%     params.complete_nuab();
% end

% %ALSO TRY:
% parameter_sweep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%validate params
if chisq_fit>chisq_typical
    warning('worse fit that typical! Using typical gab')
    params.gab = typical_gab;
    P_fit = P_typical;
    f_fit = f_typical;
end
if size(params.nus,1) > 1 || size(params.phia,1) > 1
    disp('multiple roots');  
end
[params,is_realistic] = set_nus_with_realistic_phia(params,1);
if ~is_realistic
    error('Unrealistic phia!')
end
