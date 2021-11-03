%% Calculation code for "fast and complete" energy transition


tic
clear all
close all
clc

%% simulation setup

t_L = 30;                               %a (lifetime of solar engine)
EROI = 20;                              %- (energy return on investment for solar engine)
EPBT = t_L/EROI;                        %a (energy payback time)
dt = .01;                               %a (time step for simulation)
t = 0:dt:100;                           %a (time frame for simulation)
n_t = size(t,2);

time_delay = 4;                         %a from 1.1.2018 until transition starts; during this time, emissions are at 2018 levels


alpha = [0 0.25 0.5 0.75 0.95 1];       % - (fossil replacement factor)
n_alpha = size(alpha,2);

tau(1,1,:) = t_L./((EROI).*(1-alpha));  % adjusted doubling time constant

gamma = [1 1.25 1.5 1.75 2];            % solar oversize factor
n_gamma = size(gamma,2);

beta_step = 0.001;
beta = (beta_step:beta_step:0.4)';      %- (fossil investment factor)
n_beta = size(beta,1);

P_demand = 6*10^(12);                   %W (final energy demand 2018, in electric energy equivalents, from non-renewable energy system)
dot_m_CO2_current = 35;                 %Gt/a (from non-renewable energy system)

%% calculation for remaining CO2 budget
%%% Statisitcal analysis of remaining carbon budget
%%% Finding a lognormal distribution that fits the IPCC values best


m_CO2_range = [0:1:8000];               % Gt (range for fitting log-normal distribution
n_range=size(m_CO2_range,2);

probability_data = xlsread('CO2 emissions and remaining budget.xlsx', 'remaining CO2 budget','D11:H11') ;
u = xlsread('CO2 emissions and remaining budget.xlsx', 'remaining CO2 budget','D13:H13') ;
m_CO2_1_5_data = xlsread('CO2 emissions and remaining budget.xlsx', 'remaining CO2 budget','D16:H16') ; %Gt
m_CO2_2_data = xlsread('CO2 emissions and remaining budget.xlsx', 'remaining CO2 budget','D23:H23') ; %Gt

m_CO2_2017=2313; %Gt
m_CO2_remaining_1_5_data = m_CO2_1_5_data - m_CO2_2017 ; %Gt
m_CO2_remaining_2_data = m_CO2_2_data - m_CO2_2017 ; %Gt



[mu_1_5_remaining, sigma_1_5_remaining,r2_1_5_remaining] = fitlog10normal(u,m_CO2_remaining_1_5_data,1000, 0.2);
m_1_5_mu_remaining = 10^(mu_1_5_remaining);

[mu_2_remaining, sigma_2_remaining,r2_2_remaining] = fitlog10normal(u,m_CO2_remaining_2_data,1000, 0.2);
m_2_mu_remaining = 10^(mu_2_remaining);

CDF_m_CO2_remaining_1_5 = zeros(1,n_range);
PDF_m_CO2_remaining_1_5 = zeros(1,n_range);
    for i=1:n_range
        CDF_m_CO2_remaining_1_5(1,i) = 1/2.*(1+ erf((log10(m_CO2_range(i))-log10(m_1_5_mu_remaining))./(sigma_1_5_remaining .* sqrt(2))));
        PDF_m_CO2_remaining_1_5(1,i) = 1/(m_CO2_range(i)*sigma_1_5_remaining*sqrt(2*pi)) * exp(-0.5*((log10(m_CO2_range(i))-log10(m_1_5_mu_remaining))./(sigma_1_5_remaining))^2); 
    end


CDF_m_CO2_remaining_2 = zeros(1,n_range);
PDF_m_CO2_remaining_2 = zeros(1,n_range);
    for i=1:n_range
        CDF_m_CO2_remaining_2(end,i) = 1/2.*(1+ erf((log10(m_CO2_range(i))-log10(m_2_mu_remaining))./(sigma_2_remaining .* sqrt(2))));
        PDF_m_CO2_remaining_2(end,i) = 1/(m_CO2_range(i)*sigma_2_remaining*sqrt(2*pi)) * exp(-0.5*((log10(m_CO2_range(i))-log10(m_2_mu_remaining))./(sigma_2_remaining))^2); 
    end


%% simulation

P_invest = beta.*P_demand;          %W (fossil power invested in solar growth)

P = zeros(n_beta,n_t,n_alpha);      %W (solar engine installed capacity)

for i=1:n_beta
    P(i,1,:) = 0;
    P(i,2,:) = P_invest(i,1) * t(1,2) / EPBT;
    for j=3:n_t
        for k=1:n_alpha
            P(i,j,k) = P(i,j-1,k) + P(i,2,k) .* 2.^(t(1,j-1)./tau(1,1,k));
        end
    end
end

% correction with EoL PV
P(:,(t_L/dt+1):n_t,:) = P(:,(t_L/dt+1):n_t,:) - P(:,1:(n_t - t_L/dt),:);

P_rep = zeros(n_beta,n_t,n_alpha);  %W (power to replace fossil engine)
for i=1:n_beta
    for j=1:n_t
        for k=1:n_alpha
            if alpha(1,k) .* P (i,j,k) >=  P_demand
                P_rep(i,j,k) = P_demand;   % cannot be larger than demand
            else
                P_rep(i,j,k) = alpha(1,k) .* P (i,j,k);
            end
        end
    end
end

% Calculating transition time
t_transition= zeros(n_beta,1,n_alpha,n_gamma);  %a (transition time)
dot_m_CO2 = zeros(n_beta,n_t,n_alpha,n_gamma);  %Gt/a (emissions during transition)
m_CO2 = zeros(n_beta,1,n_alpha,n_gamma);        %Gt (cumulative emissions from 2018 until t_transition)
for o=1:n_gamma
for k=1:n_alpha
for i=1:n_beta
    for j=1:n_t
        if P(i,j,k) > gamma(1,o) * P_demand
            t_transition(i,1,k,o) = t(1,j);
            break
        end
      
    end
    if P(i,end,k) <= gamma(1,o) * P_demand
         t_transition(i,1,k,o) = t(1,end);      % if transition is not possible within simulation timeframe
    end
end
end


% CO2 emissions
for i=1:n_beta
    for j=1:n_t
        for k=1:n_alpha
            if t(1,j) > t_transition(i,1,k,o)
                dot_m_CO2(i,j,k,o) = 0;         % no emissions after transition is completed
            else
                dot_m_CO2(i,j,k,o) = dot_m_CO2_current * (1 + beta(i,1) - P_rep (i,j,k)/ P_demand);
            end
        end
    end
end


m_CO2(:,:,:,o) = sum(dot_m_CO2(:,:,:,o) .* dt,2) + dot_m_CO2_current/0.8246 * time_delay; % also land use and other emissions stay constant during time delay
for i=1:n_beta
    for k=1:n_alpha
        if t_transition(i,1,k,o) == t(1,end)
            m_CO2(i,:,k,o) = 10000 + m_CO2(i,:,k,o); % if transition cannot be completed, cumulative emissions go to infinity for transition
        end
    end
end
for k=1:n_alpha
    [m_CO2_min_IR(1,1,k,o),n_i_min_IR(1,1,k,o)] = min(m_CO2(:,:,k,o),[],1); % finding minimum
end

end
[m_CO2_min,n_i_min] = min(m_CO2_min_IR,[],3); % finding minimum



%%    
toc