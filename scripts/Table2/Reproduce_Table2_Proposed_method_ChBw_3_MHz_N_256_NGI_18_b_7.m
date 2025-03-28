clear,
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                 LOAD THE PARAMETERS                   % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
addpath('../../data/Table2')
load('Table2_Parameters_for_Proposed_method_ChBw_3_MHz_N_256_NGI_18_b_7.mat');

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

Pd = GenPulse_param.Pd;
Pd_todas = GenPulse_param.Pd_todas;
Paic = GenPulse_param.Paic;
Ptk = GenPulse_param.Ptk;
R_aic = GenPulse_param.R_aic;
R_tk = GenPulse_param.R_tk;
NexTemp = GenPulse_param.NexTemp;
cell_coef_lim = GenPulse_param.cell_coef_lim;

adj_band_y = Adj_bands.adj_band_y;
left_adj_band_x = Adj_bands.left_adj_band_x;
right_adj_band_x = Adj_bands.right_adj_band_x;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            OBTAIN THE PULSES IN THE SCENARIO            % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

eta = (N+Ngi-b-1)/2; 
vn=(0:L+b-1)'-Ngi;  % Dicrete time axis vector
n0 = eta+b-Ngi;     % Phase delay due to Hermiticity


h_k = zeros((2*NexTemp+1)*L+b, N); % Matrix containing the proposed pulses
aic_k = zeros((2*NexTemp+1)*L+b, N);  % Matrix containing the AIC term
t_k = zeros((2*NexTemp+1)*L+b, N);    % Matrix containing the AST term

% Obtain h_k for every data carrier
for ind_k = 1:length(Pd);

    % % AIC terms
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Central term
    Pc = exp(1j*2*pi/N*repmat(vn,1,length(Paic{ind_k})).*repmat(Paic{ind_k}.',length(vn),1)).*repmat(g.',1,length(Paic{ind_k}));
    % Advanced and delayed terms
    Pc = repmat(Pc,1,2*NexTemp+1);

    Omega_c_n0 = diag(repmat(exp(-1j*2*pi*Paic{ind_k}*n0/N),1+2*NexTemp,1));
    aic_temp = Pc*Omega_c_n0*diag(R_aic{ind_k}).*exp(1j*2*pi*Pd(ind_k)*n0/N);
    aic_k((NexTemp*L+1):(NexTemp+1)*L+b,Pd(ind_k)+N/2+1) = sum(aic_temp(:,1:length(Paic{ind_k})),2);    % CENTRAL TERM
    
    m_temp = 0;
    for ind_m = 1:NexTemp;                                                                              % ADVANCED AND DELAYED TERMS
        aic_k(((NexTemp-ind_m)*L+1):((NexTemp-ind_m+1)*L+b),Pd(ind_k)+N/2+1) = aic_k(((NexTemp-ind_m)*L+1):((NexTemp-ind_m+1)*L+b),Pd(ind_k)+N/2+1) + sum(aic_temp(:,(length(Paic{ind_k})*(1+m_temp)) + (1:length(Paic{ind_k}))),2);
        aic_k(((NexTemp+ind_m)*L+1):((NexTemp+ind_m+1)*L+b),Pd(ind_k)+N/2+1) = aic_k(((NexTemp+ind_m)*L+1):((NexTemp+ind_m+1)*L+b),Pd(ind_k)+N/2+1) + sum(aic_temp(:,(length(Paic{ind_k})*(1+m_temp+1)) + (1:length(Paic{ind_k}))),2);
        m_temp = m_temp+2;
    end
    
    % % AST terms
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Central term
    Ts = exp(1j*2*pi/N*repmat(vn,1,length(Ptk{ind_k})).*repmat(Ptk{ind_k}.',length(vn),1)).*repmat(us.',1,length(Ptk{ind_k}));
    Te = exp(1j*2*pi/N*repmat(vn,1,length(Ptk{ind_k})).*repmat(Ptk{ind_k}.',length(vn),1)).*repmat(ue.',1,length(Ptk{ind_k}));
    T = [Ts, Te];
    
    % Advanced and delayed terms
    T = repmat(T,1,NexTemp+1);
    Omega_c_n0 = diag(repmat(exp(-1j*2*pi*Ptk{ind_k}*n0/N),2*(1+NexTemp),1));
    t_k_temp = T*Omega_c_n0*diag(R_tk{ind_k}).*exp(1j*2*pi*Pd(ind_k)*n0/N);
    t_k(((NexTemp)*L+1):((NexTemp+1))*L+b,Pd(ind_k)+N/2+1) = sum(t_k_temp(:,1:2*length(Ptk{ind_k})),2);

    m_temp = 0;
    for ind_m = 1:NexTemp;
            t_k(((NexTemp-ind_m)*L+1):(NexTemp-ind_m+1)*L+b,Pd(ind_k)+N/2+1) = t_k(((NexTemp-ind_m)*L+1):(NexTemp-ind_m+1)*L+b,Pd(ind_k)+N/2+1) + sum(t_k_temp(:,length(Ptk{ind_k})*(2+m_temp) + (1:length(Ptk{ind_k}))),2);
            t_k(((NexTemp+ind_m)*L+1):(NexTemp+ind_m+1)*L+b,Pd(ind_k)+N/2+1) = t_k(((NexTemp+ind_m)*L+1):(NexTemp+ind_m+1)*L+b,Pd(ind_k)+N/2+1) + sum(t_k_temp(:,length(Ptk{ind_k})*(2+m_temp+1) + (1:length(Ptk{ind_k}))),2);
            m_temp = m_temp + 2;
    end
end

h_k = aic_k + t_k;

p_k_conv = exp(1j*(2*pi/N)*repmat(vn,1,length(Pd_todas)).*repmat(Pd_todas.',length(vn),1)).*repmat(g.',1,length(Pd_todas));

h_k(NexTemp*L + (1:L+b),Pd_todas+N/2+1) = h_k(NexTemp*L + (1:L+b),Pd_todas+N/2+1) + p_k_conv;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     REPRESENTATION OF THE PULSES IN THE SCENARIO      % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Obtain the PSD
H_k = fft(h_k,N*Ms);

psd_H = fftshift(sum(abs(H_k).^2,2));

x = (-N*Ms/2:N*Ms/2-1)/(Ms);
indices = find(ismember(x,[-N/2,N/2-1]));
indices = indices(1):1:indices(2);

y_h = 10*log10(psd_H/max(psd_H));

% Table 2
figure,
plot(x(indices),y_h(indices), 'DisplayName', 'PSD proposed method', 'Color', "#4DBEEE"); hold on; grid on;
plot(left_adj_band_x, adj_band_y, 'DisplayName', 'Adjacent band limits', 'Color', "#A2142F", 'Linewidth', 1);
plot(right_adj_band_x, adj_band_y, 'HandleVisibility','off', 'Color', "#A2142F", 'Linewidth', 1);
xlabel('Carrier index (k)'); ylabel('Normalized PSD (dB)'); legend('show')
xlim([0, (right_adj_band_x(end))*1.05]); ylim([-100,0])
title(['Channel Bandwidth = 3 MHz; N = 256; N_{GI} = 18; \beta = 7'])

