% Maximum historic growth for PV
% Fit to logistic curve

clear variables
close all

tic

%% loading and preparing data
% time series for energy demand and CO2 emissions

EROI = 20; %- (energy return on investment)
t_L = 30; %a (lifetime)
dt = .01; %a (time step)
EPBT = t_L/EROI; % a (energy payback time)
alpha_step=0.01;
alpha = 0:alpha_step:1; % - (continuous replacement of fossil in the economy)
n_alpha = size(alpha,2);
tau(1,1,:) = t_L./((EROI).*(1-alpha)); % adjusted doubling time constant

P_demand = 6; %TW (in 2018 without renewable fraction, which doesn't need to be replaced)

beta_step = 0.001;
beta = (beta_step:beta_step:0.4)';
n_i = size(beta,1);


time_frame = xlsread('CO2 emissions and remaining budget.xlsx','RE','A6:A60'); %a
t_ext = 1980:1:2100;

P_PV = xlsread('CO2 emissions and remaining budget.xlsx','RE','K6:K60'); %TW
P_PV (isnan(P_PV))=0;

ATP_PV = 70; %TW

%% Define functions and perform non-linear regression
% logistic
logistic_curve = @(b_log,t) ATP_PV./(1+exp(-b_log(1).*(t-b_log(2))));

PV_fit = fitnlm(time_frame,P_PV, logistic_curve,[1 2020]);
b_log_estimates = PV_fit.Coefficients.Estimate;
R2_PV = PV_fit.Rsquared.Adjusted;

% exponential
exponential_curve = @(b,t) b(2) .* exp(b(3).*(t-b(1)));

PV_fit_exp = fitnlm(time_frame,P_PV, exponential_curve,[1980 10^(-7) 0.25]);
b_estimates_exp = PV_fit_exp.Coefficients.Estimate;
R2_PV_exp = PV_fit_exp.Rsquared.Adjusted;



%% Search for alpha and beta produce historically fitted growth rate

P_invest = beta.*P_demand; %W
t_PV = 2021:dt:2100; %a
n_t = size(t_PV,2);

P = zeros(n_i,n_t,n_alpha);

for i=1:n_i
    P(i,1,:) = 0;
    P(i,2,:) = P_invest(i,1) * (t_PV(1,2)-t_PV(1,1)) / EPBT;
    for j=3:n_t
        for k=1:n_alpha
            P(i,j,k) = P(i,j-1,k) + P(i,2,k) .* 2.^((t_PV(1,j-1)-t_PV(1,1))./tau(1,1,k));
            
        end
    end
end

% correction with EoL PV
P(:,(t_L/dt+1):n_t,:) = P(:,(t_L/dt+1):n_t,:) - P(:,1:(n_t - t_L/dt),:);

P_rep = zeros(n_i,n_t,n_alpha);
for i=1:n_i
    for j=1:n_t
        for k=1:n_alpha
            if P(i,j,k)>ATP_PV
                P(i,j,k)=ATP_PV;
            end
            if alpha(1,k) .* P (i,j,k) >=  P_demand
                P_rep(i,j,k) = P_demand;
            else
                P_rep(i,j,k) = alpha(1,k) .* P (i,j,k);
            end
        end
    end
end


exp_curve_cap = exponential_curve(b_estimates_exp,t_PV);
for j=1:n_t
    if exp_curve_cap(1,j)>ATP_PV
        exp_curve_cap(1,j)=ATP_PV;
    end
end

%% Finding best fit

t_fit = [2021:dt:2037];
n_t_fit = size(t_fit,2);

err_logistic_P = zeros(n_i, n_alpha);
err_exp_P = zeros(n_i, n_alpha);
for i=1:n_i
    for k=1:n_alpha
        err_logistic_P(i,k)  = immse((logistic_curve(b_log_estimates,t_fit)-logistic_curve(b_log_estimates,t_PV(1,1))),P(i,1:n_t_fit,k));
        err_exp_P(i,k)  = immse((exp_curve_cap(1,1:n_t_fit)-exp_curve_cap(1,1)),P(i,1:n_t_fit,k));
    end
end
[err_log_best_fit_beta,n_beta_fit_log] = min(err_logistic_P,[],1);
[err_log_best_fit,n_alpha_fit_log] = min(err_log_best_fit_beta,[],2);
beta_fit_log = beta(n_beta_fit_log(n_alpha_fit_log));
alpha_fit_log = alpha(n_alpha_fit_log);

[err_exp_best_fit_beta,n_beta_fit_exp] = min(err_exp_P,[],1);
[err_exp_best_fit,n_alpha_fit_exp] = min(err_exp_best_fit_beta,[],2);
beta_fit_exp = beta(n_beta_fit_exp(n_alpha_fit_exp));
alpha_fit_exp = alpha(n_alpha_fit_exp);


%% Calculating CO2 emissions for best fit (log)

dot_m_CO2_current = 35; %Gt/a



for k=1:n_alpha_fit_log
for i=1:n_beta_fit_log(n_alpha_fit_log)
    for j=1:n_t
        if P(i,j,k) > 2 * P_demand
            t_transition_fit = t_PV(1,j);
            break
        end
      
    end
    if P(i,end,k) <= 2 * P_demand
         t_transition_fit = t_PV(1,end);
    end
end
end


% CO2 emissions
for i=1:n_beta_fit_log(n_alpha_fit_log)
    for j=1:n_t
        for k=1:n_alpha_fit_log
            if t_PV(1,j) > t_transition_fit
                dot_m_CO2(1,j) = 0;
            else
                dot_m_CO2(1,j) = dot_m_CO2_current * (1 + beta(i,1) - P_rep (i,j,k)/ P_demand);
            end
        end
    end
end


m_CO2 = sum(dot_m_CO2 .* dt,2) + dot_m_CO2_current/0.8246 * 4; % also land use and other emissions stay constant during time delay


%%
toc