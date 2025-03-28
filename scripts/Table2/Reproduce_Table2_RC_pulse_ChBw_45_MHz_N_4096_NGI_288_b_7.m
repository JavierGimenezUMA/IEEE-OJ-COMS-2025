clear,
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                 LOAD THE PARAMETERS                   % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
addpath('../../data/Table2')
load('Table2_Parameters_for_RC_pulse_ChBw_45_MHz_N_4096_NGI_288_b_7.mat');

Fsob = OFDM_param.Fsob; % Upsampling factor 
% Upsample in time-domain
N = OFDM_param.N*Fsob;
Ngi = OFDM_param.Ngi*Fsob;
b = (OFDM_param.b+1)*Fsob-1;
L = OFDM_param.L*Fsob;
Ms = OFDM_param.Ms;
g = OFDM_param.g;
us = OFDM_param.us;
ue = OFDM_param.ue;

ur = us + ue;

% Pd = GenPulse_param.Pd;
Pd_todas = GenPulse_param.Pd_todas;
% Paic = GenPulse_param.Paic;
% Ptk = GenPulse_param.Ptk;
% R_aic = GenPulse_param.R_aic;
% R_tk = GenPulse_param.R_tk;
% NexTemp = GenPulse_param.NexTemp;
% cell_coef_lim = GenPulse_param.cell_coef_lim;

adj_band_y = Adj_bands.adj_band_y;
left_adj_band_x = Adj_bands.left_adj_band_x;
right_adj_band_x = Adj_bands.right_adj_band_x;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            OBTAIN THE PULSES IN THE SCENARIO            % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

vn=(0:L+b-1)'-Ngi;  % Dicrete time axis vector
p_k_conv = exp(1j*(2*pi/N)*repmat(vn,1,length(Pd_todas)).*repmat(Pd_todas.',length(vn),1)).*repmat(g.',1,length(Pd_todas));




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     REPRESENTATION OF THE PULSES IN THE SCENARIO      % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Obtain the PSD

P_k = fft(p_k_conv,N*Ms);
psd_P = fftshift(sum(abs(P_k).^2,2));

x = (-N*Ms/2:N*Ms/2-1)/(Ms);
indices = find(ismember(x,[-N/2,N/2-1]));
indices = indices(1):1:indices(2);

% y_h = 10*log10(psd_H/max(psd_H));
y_p = 10*log10(psd_P/max(psd_P));

% Table 2
figure,
plot(x(indices),y_p(indices), 'DisplayName', 'PSD RC pulse-shaping', 'Color', 'r');hold on; grid on;
plot(left_adj_band_x, adj_band_y, 'DisplayName', 'Adjacent band limits', 'Color', "#A2142F", 'Linewidth', 1);
plot(right_adj_band_x, adj_band_y, 'HandleVisibility','off', 'Color', "#A2142F", 'Linewidth', 1);
xlabel('Carrier index (k)'); ylabel('Normalized PSD (dB)'); legend('show')
xlim([0, (right_adj_band_x(end))*1.05]); ylim([-100,0])
title(['Channel Bandwidth = 45 MHz; N = 4096; N_{GI} = 288; \beta = 7'])

