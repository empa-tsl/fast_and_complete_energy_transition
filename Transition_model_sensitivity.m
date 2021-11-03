%% Calculation code for sensitivity of "fast and complete" energy transition

tic
clear all
close all
clc

%% simulation setup

dt = .01;                   %a (time step)
t = 0:dt:100;               %a
n_t = size(t,2);

time_delay = 4;             %a from 1.1.2018 until transition starts; during this time, emissions are at 2018 levels

alpha = [0 1];              % - (fossil replacement factor)
n_alpha = size(alpha,2);

gamma = [2];                % solar oversize factor


P_demand = 6*10^(12);       %W (final energy demand 2018, in electric energy equivalents, from non-renewable energy system)
dot_m_CO2_current = 35; 	%Gt/a (from non-renewable energy system)
beta_step = 0.001;
beta = (beta_step:beta_step:0.4)';  %- (fossil investment factor)
n_beta = size(beta,1);



%% EROI calculations
n_omega = 3;    %sensitivity options

t = 0:dt:100; %a
n_t = size(t,2);

t_transition= zeros(n_beta,1,n_alpha,n_omega);
dot_m_CO2 = zeros(n_beta,n_t,n_alpha,n_omega);
m_CO2 = zeros(n_beta,1,n_alpha,n_omega);

for o=1:n_omega
p_ins = 164 * 365.25*24/1000; %kWh/m^2 a global average insolation per area
degradation = 0.005; %1/a annual loss of output oof PV system due to ageing
t_L(1,o) = 30; %a lifetime of PV system
if o==1
    eta(1,o) = 0.14; % global average of solar PV efficiency
    E_inv(1,o) = 350; %kWh/m^2 investment of energy to build PV system
    
elseif o==2
    eta(1,o) = 0.15; % global average of solar PV efficiency
    E_inv(1,o) = 300; %kWh/m^2 investment of energy to build PV system

elseif o==3
    eta(1,o) = 0.20; % global average of solar PV efficiency
    E_inv(1,o) = 280; %kWh/m^2 investment of energy to build PV system
end

EROI(1,o) = (p_ins * eta(1,o) / degradation)*(1-exp(-degradation*t_L(1,o)))/E_inv(1,o);
EPBT(1,o) = E_inv(1,o) / (p_ins * eta(1,o));

tau(1,1,:) = t_L(1,o)./((EROI(1,o)).*(1-alpha)); % adjusted doubling time constant



%% simulation

P_invest = beta.*P_demand; %W


P = zeros(n_beta,n_t,n_alpha);

for i=1:n_beta
    P(i,1,:) = 0;
    P(i,2,:) = P_invest(i,1) * t(1,2) / EPBT(1,o);
    for j=3:n_t
        for k=1:n_alpha
            P(i,j,k) = P(i,j-1,k) + P(i,2,k) .* 2.^(t(1,j-1)./tau(1,1,k));
        end
    end
end

% correction with EoL PV
P(:,(t_L(1,o)/dt+1):n_t,:) = P(:,(t_L(1,o)/dt+1):n_t,:) - P(:,1:(n_t - t_L(1,o)/dt),:);

P_rep = zeros(n_beta,n_t,n_alpha);
for i=1:n_beta
    for j=1:n_t
        for k=1:n_alpha
            if alpha(1,k) .* P (i,j,k) >=  P_demand
                P_rep(i,j,k) = P_demand;
            else
                P_rep(i,j,k) = alpha(1,k) .* P (i,j,k);
            end
        end
    end
end

% Calculating transition time


for k=1:n_alpha
for i=1:n_beta
    for j=1:n_t
        if P(i,j,k) > gamma * P_demand
            t_transition(i,1,k,o) = t(1,j);
            break
        end
      
    end
    if P(i,end,k) <= gamma * P_demand
         t_transition(i,1,k,o) = t(1,end);
    end
end
end


% CO2 emissions
for i=1:n_beta
    for j=1:n_t
        for k=1:n_alpha
            if t(1,j) > t_transition(i,1,k,o)
                dot_m_CO2(i,j,k,o) = 0;
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
    [m_CO2_min_IR(1,1,k,o),n_i_min_IR(1,1,k,o)] = min(m_CO2(:,:,k,o),[],1);
end

end
[m_CO2_min,n_i_min] = min(m_CO2_min_IR,[],3);



%%    
toc