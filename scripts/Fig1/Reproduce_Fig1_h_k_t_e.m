clear,

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                 LOAD THE PARAMETERS                   % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
addpath('../../data/Fig1')
load('Fig1_Parameters_for_h_k_t-e.mat');

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
Ptk = GenPulse_param.Ptk;
R_aic = GenPulse_param.R_aic;
R_tk = GenPulse_param.R_tk;
NexTemp_a = GenPulse_param.NexTemp_a;
NexTemp_d = GenPulse_param.NexTemp_d;
cell_coef_lim = GenPulse_param.cell_coef_lim;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            OBTAIN THE PULSES IN THE SCENARIO            % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

vn=(0:L+b-1)'-Ngi;  % Dicrete time axis vector

h_k = zeros((NexTemp_a+NexTemp_d+1)*L+b, N); % Matrix containing the proposed pulses
aic_k = zeros((NexTemp_a+NexTemp_d+1)*L+b, N);  % Matrix containing the AIC term
t_k = zeros((NexTemp_a+NexTemp_d+1)*L+b, N);    % Matrix containing the AST term

% Obtain h_k for every data carrier
for ind_k = 1:length(Pd)

    % % AIC terms
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Central term
    Pc = exp(1j*2*pi/N*repmat(vn,1,length(Paic{ind_k})).*repmat(Paic{ind_k}.',length(vn),1)).*repmat(g.',1,length(Paic{ind_k}));
    % Advanced and delayed terms
    for ind_m = -NexTemp_a:NexTemp_d
        if ind_m ~= 0 % Skip the central term
            Pc = [Pc, exp(1j*2*pi/N*repmat(vn,1,length(Paic{ind_k})).*repmat(Paic{ind_k}.',length(vn),1)).*repmat(g.',1,length(Paic{ind_k}))];
        end
    end

    aic_temp = Pc*diag(R_aic{ind_k}); 
    aic_k((NexTemp_a*L+1):(NexTemp_a+1)*L+b,Pd(ind_k)+N/2+1) = sum(aic_temp(:,1:length(Paic{ind_k})),2);    % CENTRAL TERM

    m_temp = 0;                                                                                             
    for ind_m = -NexTemp_a:NexTemp_d                                                                      % aDVANCED AND DELAYED TERMS
        if ind_m ~= 0 % Skip the central term
            aic_k(((NexTemp_a+ind_m)*L+1):((NexTemp_a+ind_m+1)*L+b),Pd(ind_k)+N/2+1) = aic_k(((NexTemp_a+ind_m)*L+1):((NexTemp_a+ind_m+1)*L+b),Pd(ind_k)+N/2+1) + sum(aic_temp(:,(length(Paic{ind_k})*(1+m_temp)) + (1:length(Paic{ind_k}))),2);
            m_temp = m_temp+1;
        end
    end
    
    % % AST terms
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Central term
    Ts = exp(1j*2*pi/N*repmat(vn,1,length(Ptk{ind_k})).*repmat(Ptk{ind_k}.',length(vn),1)).*repmat(us.',1,length(Ptk{ind_k}));
    Te = exp(1j*2*pi/N*repmat(vn,1,length(Ptk{ind_k})).*repmat(Ptk{ind_k}.',length(vn),1)).*repmat(ue.',1,length(Ptk{ind_k}));
    T = [Ts, Te];
    
    % Advanced and delayed terms
    for ind_m = -NexTemp_a:-1
        T = [T, Ts];
    end
    for ind_m=1:NexTemp_d
        T = [T, Te];
    end
    t_k_temp = T*diag(R_tk{ind_k});
    t_k(((NexTemp_a)*L+1):((NexTemp_a+1))*L+b,Pd(ind_k)+N/2+1) = sum(t_k_temp(:,1:2*length(Ptk{ind_k})),2);
    
    m_temp = 0;
    for ind_m = -NexTemp_a:-1
            t_k(((NexTemp_a+ind_m)*L+1):(NexTemp_a+ind_m+1)*L+b,Pd(ind_k)+N/2+1) = t_k(((NexTemp_a+ind_m)*L+1):(NexTemp_a+ind_m+1)*L+b,Pd(ind_k)+N/2+1) + sum(t_k_temp(:,(2*length(Ptk{ind_k})+m_temp*length(Ptk{ind_k})) + (1:length(Ptk{ind_k}))),2);
            m_temp = m_temp + 1;
    end
    for ind_m = 1:NexTemp_d
            t_k(((NexTemp_a+ind_m)*L+1):(NexTemp_a+ind_m+1)*L+b,Pd(ind_k)+N/2+1) = t_k(((NexTemp_a+ind_m)*L+1):(NexTemp_a+ind_m+1)*L+b,Pd(ind_k)+N/2+1) + sum(t_k_temp(:,(2*length(Ptk{ind_k})+m_temp*length(Ptk{ind_k})) + (1:length(Ptk{ind_k}))),2);
            m_temp = m_temp + 1;
    end
end

h_k = aic_k + t_k;

p_k_conv = exp(1j*(2*pi/N)*repmat(vn,1,length(Pd)).*repmat(Pd.',length(vn),1)).*repmat(g.',1,length(Pd));

h_k(NexTemp_a*L + (1:L+b),Pd+N/2+1) = h_k(NexTemp_a*L + (1:L+b),Pd+N/2+1) + p_k_conv;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     REPRESENTATION OF THE PULSES IN THE SCENARIO      % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Obtain the PSD
H_k = fft(h_k(:,Pd+N/2+1),N*Ms);
P_k = fft(p_k_conv,N*Ms);

psd_H = fftshift(sum(abs(H_k).^2,2));
psd_P = fftshift(sum(abs(P_k).^2,2));

x = (-N*Ms/2:N*Ms/2-1)/(Ms);
indices = find(ismember(x,[-N/2, N/2-1]));
indices = indices(1):1:indices(2);

x = x+N/2;
y_h = 10*log10(psd_H/max(psd_H));
y_p = 10*log10(psd_P/max(psd_P));

% Fig 1 (a)
figure,
plot(x(indices),y_h(indices), 'DisplayName', '|H_k^{t-e}(f)|'); hold on; grid on;
plot(x(indices),y_p(indices), 'DisplayName', '|P_k(f)|','Color','r');
xlabel('Carrier k'); ylabel('Magnitude (dB)'); legend('show')
xlim([3009, 3029]); ylim([-125,0])

x_p = 0:(L+b-1);
x_h = -L:(2*L+b-1);
y_h = abs(h_k(:,Pd+N/2+1));
y_p = abs(p_k_conv);

% Fig 1(b)
figure,
plot(x_h,y_h, 'DisplayName', '|h_k^{t-e}(n)|'); grid on; hold on;
plot(x_p,y_p, 'DisplayName', '|p_k(n)|','Color','r');
xlabel('Sample index (n)'); ylabel('Amplitude'); legend('show')
xlim([-5000,11000]); ylim([-0.7,2.2])
