% https://github.com/BrainDynamicsUSYD/braintrak/wiki/Spatial-fitting
%
clear all;close all;
addpath('D:\My Files\Work\BGU\scripts\nft\braintrak');
addpath('D:\My Files\Work\BGU\scripts\nft\corticothalamic-model');

debugmode = false; %to test chain 

%load data and compute spectra
[fn, fp] = uigetfile(['D:\My Files\Work\BGU\Datasets\DOC\mat\' '*.mat'], 'Select data file');
load([fp fn]);
EEGsect.data([EEGsect.chanlocs.radius] > 0.5,:) = [];
EEGsect.chanlocs([EEGsect.chanlocs.radius] > 0.5) = [];
[P, f] = compute_psd('fft', EEGsect.data, 3, EEGsect.srate, [1 45], false);
% spatial_data = load('demo_spatial.mat');
% f = spatial_data.f(2:end); P = spatial_data.P(2:end,:);
% figure;loglog(f,P)

% spatial_model = bt.model.spatial_express_t0();
% spatial_model.set_electrodes(spatial_data.electrodes);
spatial_model = bt.model.spatial_gab();
spatial_model.set_electrodes(EEGsect);

[init_default,prior_pp] = spatial_model.initialize_fit(f,P');
fit_result =  bt.core.fit_spectrum(spatial_model,f,P',prior_pp,init_default,'40s',[],[],debugmode);
f_fit = fit_result.fit_data.target_f;
P_fit = fit_result.fit_data.fitted_P';
    
% fit_result.plot()
figure;
loglog(f,mean(P,1), f_fit,mean(P_fit,1));xlabel('Hz');title('fitted spectrum');
legend('experimental','fitted');
fit_result.head_plot(true)
