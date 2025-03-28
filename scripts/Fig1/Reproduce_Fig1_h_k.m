clear,
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                 LOAD THE PARAMETERS                   % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
addpath('../../data/Fig1')
load('Fig1_Parameters_for_h_k.mat');

N = OFDM_param.N;
Ngi = OFDM_param.Ngi;
b = OFDM_param.b;
L = OFDM_param.L;
Ms = OFDM_param.Ms;
g = OFDM_param.g;
us = OFDM_param.us;
ue = OFDM_param.ue;

ur = us + ue;

Pd = GenPulse_param.Pd;
Paic = GenPulse_param.Paic;
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
        m_temp = m_temp+1;
    end
    
    % % AST terms
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Central term

    % Phase shift applied to optimal coefficients
    phase_shift_tk = exp(1j*2*pi*Pd(ind_k)*n0/N);
    R_tk{ind_k} = R_tk{ind_k}*phase_shift_tk;
    
    t_k(((NexTemp)*L+1):((NexTemp+1))*L+b,Pd(ind_k)+N/2+1) = ur.'.*[[eye(b),zeros(b)];zeros(L-b,2*b);[zeros(b),eye(b)]]*[flipud(R_tk{ind_k}(1:b));R_tk{ind_k}((b+1):2*b)];
    
    % Advanced and delayed terms
    m_temp = 0;
    for ind_m = 1:NexTemp
            t_k(((NexTemp-ind_m)*L+1):((NexTemp-ind_m))*L+b,Pd(ind_k)+N/2+1) = t_k(((NexTemp-ind_m)*L+1):((NexTemp-ind_m))*L+b,Pd(ind_k)+N/2+1) + fliplr(R_tk{ind_k}((2+m_temp)*b + (1:b))).*us(n_tk_s+1).';
            t_k(((NexTemp+ind_m)*L+1):((NexTemp+ind_m))*L+b,Pd(ind_k)+N/2+1) = t_k(((NexTemp+ind_m)*L+1):((NexTemp+ind_m))*L+b,Pd(ind_k)+N/2+1) + R_tk{ind_k}((2+m_temp+1)*b + (1:b)).*ue(n_tk_e+1).';
            m_temp = m_temp + 1;
    end
end

h_k = aic_k + t_k;

p_k_conv = exp(1j*(2*pi/N)*repmat(vn,1,length(Pd)).*repmat(Pd.',length(vn),1)).*repmat(g.',1,length(Pd));

h_k(NexTemp*L + (1:L+b),Pd+N/2+1) = h_k(NexTemp*L + (1:L+b),Pd+N/2+1) + p_k_conv;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     REPRESENTATION OF THE PULSES IN THE SCENARIO      % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Obtain the PSD
H_k = fft(h_k(:,Pd+N/2+1),N*Ms);
P_k = fft(p_k_conv,N*Ms);

psd_H = fftshift(sum(abs(H_k).^2,2));
psd_P = fftshift(sum(abs(P_k).^2,2));

x = (-N*Ms/2:N*Ms/2-1)/(Ms);
indices = find(ismember(x,[-N/2,N/2-1]));
indices = indices(1):1:indices(2);

x = x+N/2;
y_h = 10*log10(psd_H/max(psd_H));
y_p = 10*log10(psd_P/max(psd_P));

% Fig 1 (a)
figure,
plot(x(indices),y_h(indices), 'DisplayName', '|H_k(f)|','Color','#EDB120'); hold on; grid on;
plot(x(indices),y_p(indices), 'DisplayName', '|P_k(f)|','Color','r');
xlabel('Carrier k'); ylabel('Magnitude (dB)'); legend('show')
xlim([3009, 3029]); ylim([-125,0])

x = 0:(L+b-1);
y_h = abs(h_k(:,Pd+N/2+1));
y_p = abs(p_k_conv);

% Fig 1(b)
figure,
plot(x,y_h, 'DisplayName', '|h_k(n)|','Color','#EDB120'); grid on; hold on;
plot(x,y_p, 'DisplayName', '|p_k(n)|','Color','r');
xlabel('Sample index (n)'); ylabel('Amplitude'); legend('show')
xlim([-5000,11000]); ylim([-0.7,2.2])

