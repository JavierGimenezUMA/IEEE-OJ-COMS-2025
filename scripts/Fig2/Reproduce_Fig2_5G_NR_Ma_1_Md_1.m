clear,
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                 LOAD THE PARAMETERS                   % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
addpath('../../data/Fig2')
load('Fig2_Parameters_for_5G_NR_Ma_1_Md_1.mat');

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


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     REPRESENTATION OF THE PULSES IN THE SCENARIO      % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Obtain the PSD
H_k = fft(h_k,N*Ms);
P_k = fft(p_k_conv,N*Ms);

psd_H = fftshift(sum(abs(H_k).^2,2));
psd_P = fftshift(sum(abs(P_k).^2,2));

x = (-N*Ms/2:N*Ms/2-1)/(Ms);
indices = find(ismember(x,[-N/2,N/2-1]));
indices = indices(1):1:indices(2);

x = x+N/2;
y_h = 10*log10(psd_H/max(psd_H));
y_p = 10*log10(psd_P/max(psd_P));

PSD_max = max(max(y_h(find(x>Pd_todas(end)+N/2+1.67))), max(y_h(find(x<Pd_todas(1)+N/2-1.67)))); %% PSD_MAX IS SOUGHT BEYOND THE 0.67 SCS GUARD BAND THAT IS LEFT BY EACH SIDE
x_PSD_max = x(find(y_h==PSD_max));

% Fig 2 (5G NR Ma=1, Md=1)
figure,
plot(x(indices),y_h(indices), 'DisplayName', 'h_k^{t-e}(M_a=M_d=1)', 'Color', "#4DBEEE"); hold on; grid on;
plot(x(indices),y_p(indices), 'DisplayName', 'RC pulse-shaping', 'Color', 'r');
plot(x_PSD_max, PSD_max, 'Marker', 'x', 'Color', 'k', 'LineStyle', 'none','LineWidth',1, 'DisplayName', ['PSD_{max} = ', num2str(PSD_max, '%.2f'), ' dB'])
xlabel('Carrier index (k)'); ylabel('Normalized PSD (dB)'); legend('show')
xlim(x_PSD_max+[-20,20]); ylim([-60,0])


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     CLARIFYING NOTE     % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

NOTE = ['NOTE: The maximum number of RB that can be used in a 5 MHz channel is 25, which is equivalent to 4.5 MHz (SCS = 15 kHz).',...
    ' That leaves a 500 kHz guard band, 250 kHz by each side of the transmission band.',...
    ' If this guard band is reduced by repurposing those carriers that were previously assigned to it for data transmission, it can be reduced by up to 16 subcarriers, leaving a residual guard band 10 kHz wide (2/3 of an SCS).',...
    ' Hence, the PSD peak can be sought beyond that guard band, which is 1.67 SCS beyond the outermost transmitting subcarrier by each end of the passband.'];
wrapped_NOTE = textwrap({NOTE}, 70);
disp('##################################################################################')
fprintf('\n')
disp(wrapped_NOTE)
disp('##################################################################################')
