% Funcion que calcula la papr de una señal OFDM usando el pulso hNs

function [ccdf, papr, mean_papr] = estima_PAPR(h, L, b, NexTemp, num_tandas, num_sim_tanda)
    %   h:              Matriz que contiene los pulsos de todas las portadoras no nulas
    %                       usadas en el escenario
    %   L:              N + Ngi
    %   b:              Número de muestras de las transiciones entre símbolos
    %   NexTemp:        Número de extensiones temporales utilizadas
    %   num_tandas:     Número de tandas de símbolos a simular para la
    %                       estima
    %   num_sim_tanda:  Número de símbolos que se simulan por cada tanda


    Nu = size(h,2); % Número de portadoras no nulas

    % Inicializamos variables
    avg_power = 0;
    amp_max_2 = zeros(1,num_sim_tanda*num_tandas-num_tandas*2);


    textprogressbar('Rounds: ');
    for ind_tanda = 1:num_tandas
        textprogressbar(100*ind_tanda/num_tandas);
        x=zeros(1,num_sim_tanda*L+b+2*L*NexTemp);    %Ahora tendremos NexTemp pulsos adicionales antes y despues de los pulsos con datos
        d = ((randi(2,Nu,num_sim_tanda)-1.5)*2+1j*(randi(2,Nu,num_sim_tanda)-1.5)*2)/sqrt(2);

        for ind_sim=1:num_sim_tanda;
            x_i = h*d(:,ind_sim);
            x(1+(ind_sim-1)*L:ind_sim*L+b+2*L*NexTemp) = x(1+(ind_sim-1)*L:ind_sim*L+b+2*L*NexTemp) + x_i.';
        end;

        x_sims = reshape(x(L*(NexTemp+1)+1:end-b-L*(NexTemp+1)),L,num_sim_tanda-2); % Se descarta el primer y último símbolo de cada tanda. Se ponen los símbolos en columnas
        avg_power = avg_power + sum(sum(abs(x_sims).^2))/(L*num_sim_tanda); % Potencia de toda la señal generada
        amp_max_2(1+(ind_tanda-1)*(num_sim_tanda-2):ind_tanda*(num_sim_tanda-2))= max(abs(x_sims),[],1).^2; % Amplitud máxima de cada símbolo

    end;
    textprogressbar('done');

    %Estima la CCDF de la PAPR
    avg_power = avg_power/num_tandas;
    [cdf,papr]=ecdf(amp_max_2./avg_power);
    ccdf=1-cdf;

    mean_papr = mean(amp_max_2./avg_power);



end
