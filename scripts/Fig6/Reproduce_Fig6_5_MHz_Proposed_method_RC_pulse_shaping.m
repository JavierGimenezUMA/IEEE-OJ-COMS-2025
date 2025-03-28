clear, close all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  PAPR ESTIMATION PARAMETER  % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

n_rounds = 20;          % Number of rounds of symbols simulated
n_sym_per_round = 5e3;  % Number of symbols per round

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                 LOAD THE PARAMETERS                   % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
addpath('../../data/Fig6')
load('Fig6_Parameters_for_5_MHz_Proposed_method_&_RC_pulse_shaping.mat');

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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %           ESTIMATE OF THE PAPR            % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

[ccdf_papr_p, papr_p, media_papr_p] = estima_PAPR(p_k_conv, L, b, 0, n_rounds, n_sym_per_round);
[ccdf_papr_h, papr_h, media_papr_h] = estima_PAPR(h_k, L, b, NexTemp, n_rounds, n_sym_per_round);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %          REPRESENTATION OF CCDF OF THE PAPR           % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%


x_p = 10*log10(papr_p);
y_p = ccdf_papr_p;
x_h = 10*log10(papr_h);
y_h = ccdf_papr_h;

figure,
semilogy(x_p,y_p, 'DisplayName', 'RC pulse shaping (5 MHz)', 'Color', 'r'); hold on;
semilogy(x_h,y_h, 'DisplayName', 'Proposed method (5 MHz)', 'Color', "#0072BD"); grid on;
xlabel('PAPR (dB)'); ylabel('CCDF'); legend('show')
xlim([6,12.5]); ylim([1e-4,1])
