%% Historical and projected emissions

tic
clear variables

%% loading and preparing data
% time series for energy demand and CO2 emissions

time_frame = xlsread('CO2 emissions and remaining budget.xlsx','CO2 time series','A4:A172');

P_el_historic_by_source = xlsread('CO2 emissions and remaining budget.xlsx','CO2 time series','X4:AA172');
P_el_historic_by_source (isnan(P_el_historic_by_source))=0;

dot_m_CO2_historic_by_source = xlsread('CO2 emissions and remaining budget.xlsx','CO2 time series','B4:F172');
dot_m_CO2_historic_by_source (isnan(dot_m_CO2_historic_by_source))=0;

m_CO2_cumulative_by_source = cumsum(dot_m_CO2_historic_by_source,1);

IEA_energy_scenario_data = xlsread('CO2 emissions and remaining budget.xlsx','IEA','A21:I23');
IEA_emission_scenario_data = xlsread('CO2 emissions and remaining budget.xlsx','IEA','A27:D29');

IEA_range = 2018:1:2040;
for i=1:size(IEA_emission_scenario_data,2)-1
    IEA_emission_scenario_fit = fit(IEA_emission_scenario_data(:,1),IEA_emission_scenario_data(:,i+1),'linearinterp');
    IEA_emission_scenario(:,i) = IEA_emission_scenario_fit(IEA_range);
end

IEA_emission_scenario_cumulated = m_CO2_cumulative_by_source(end-1,1:3) + cumsum(IEA_emission_scenario,1);

IPCC_emission_scenario_data = xlsread('CO2 emissions and remaining budget.xlsx','IPCC','C7:F19')';
IPCC_emission_scenario_sigma = xlsread('CO2 emissions and remaining budget.xlsx','IPCC','K8:N13')';

IPCC_range = 2018:1:2100;
for i=1:12
    IPCC_emission_scenario_fit = fit(IPCC_emission_scenario_data(:,1),IPCC_emission_scenario_data(:,i+1),'linearinterp');
    IPCC_emission_scenario(:,i) = IPCC_emission_scenario_fit(IPCC_range);
end
for i=1:6
    IPCC_emission_scenario_FF(:,i) = IPCC_emission_scenario(:,(i*2-1));
    IPCC_emission_scenario_total(:,i) = IPCC_emission_scenario(:,(i*2));
end

for i=1:6
    IPCC_sigma_fit = fit(IPCC_emission_scenario_data(:,1),IPCC_emission_scenario_sigma(:,i),'linearinterp');
    IPCC_sigma(:,i) = IPCC_sigma_fit(IPCC_range);
end

IPCC_emission_scenario_cumulated_total = zeros(size(IPCC_range,2),6);
IPCC_emission_scenario_cumulated_FF = zeros(size(IPCC_range,2),6);
for i=1:6
    IPCC_emission_scenario_cumulated_FF(:,i) = sum(m_CO2_cumulative_by_source(end-1,1:4),2) + cumsum(IPCC_emission_scenario(:,(i*2-1)),1);
    IPCC_emission_scenario_cumulated_total(:,i) = sum(m_CO2_cumulative_by_source(end-1,:),2) + cumsum(IPCC_emission_scenario(:,(i*2)),1);
end


beta_disp = 0.1;
n_alpha_disp = 1;

%% Uncertainty of historic emissions

n_runs = 10^5;

mu_ff = sum(dot_m_CO2_historic_by_source(:,1:4),2); % fossil and industry emissions
sigma_ff = 0.05 .* mu_ff;   % sigma is proportional to mu (Friedlingstein et al. 2019)

mu_l = dot_m_CO2_historic_by_source(:,5);   % land use emissions
sigma_l = 2.56; % Gt/a (Friedlingstein et al. 2019)

dot_m_rnd = zeros(size(mu_ff,1),2,n_runs);
for i=1:size(mu_ff,1)
    dot_m_rnd(i,1,:) = normrnd(mu_ff(i,1),sigma_ff(i,1),n_runs,1);
    dot_m_rnd(i,2,:) = normrnd(mu_l(i,1),sigma_l,n_runs,1);
end

m_cum_rnd = cumsum(dot_m_rnd,1);


CI = 0.99;      % confidence intervall (can be chosen freely)
m_cum(:,2) = quantile(sum(m_cum_rnd,2),0.5,3);
m_cum(:,1) = quantile(sum(m_cum_rnd,2),(1-CI)/2,3);
m_cum(:,3) = quantile(sum(m_cum_rnd,2),1-(1-CI)/2,3);


mu_IPCC = IPCC_emission_scenario(:,[2 4 6 8 10 12]);    % total emissions
sigma_IPCC = IPCC_sigma;

dot_m_rnd_IPCC = zeros(size(IPCC_range,2),6,n_runs);
for i=1:size(IPCC_range,2)
    for j=1:6
        dot_m_rnd_IPCC(i,j,:) = normrnd(mu_IPCC(i,j),sigma_IPCC(i,j),n_runs,1);
    end
end

m_cum_rnd_IPCC = sum(m_cum_rnd(end-1,:,:),2) + cumsum(dot_m_rnd_IPCC,1);

m_cum_IPCC(:,:,2) = quantile(m_cum_rnd_IPCC,0.5,3);
m_cum_IPCC(:,:,1) = quantile(m_cum_rnd_IPCC,(1-CI)/2,3);
m_cum_IPCC(:,:,3) = quantile(m_cum_rnd_IPCC,1-(1-CI)/2,3);

%% Calculating negative emissions in IPCC pathways

CDR_IPCC = IPCC_emission_scenario_total - IPCC_emission_scenario_FF;
for i=1:size(IPCC_range,2)
    for j=1:6
        if CDR_IPCC(i,j) > 0
            CDR_IPCC(i,j) = 0;
        end
    end
end
CDR_cum = cumsum(CDR_IPCC,1);

Net_negative = max(IPCC_emission_scenario_cumulated_total,[],1) - IPCC_emission_scenario_cumulated_total(end,:);

%% Calculation of target violation

% uncertainty distribution of 1.5°C and 2°C targets

m_CO2_2017 = 2313;
mu_1_5 = 2.774;         %data from regression (see "Transition model")
sigma_1_5 = 0.3967;

mu_2 = 3.196;
sigma_2 = 0.3267;

m_target_1_5 = lognrnd(log(10^(mu_1_5)),log(10^(sigma_1_5)),n_runs,1) + m_CO2_2017;
m_target_2 = lognrnd(log(10^(mu_2)),log(10^(sigma_2)),n_runs,1) + m_CO2_2017;


[IPCC_max, i_t_max] = max(IPCC_emission_scenario_cumulated_total,[],1);
t_max = IPCC_range(1,i_t_max);

%% probability of violation
n_scenarios = 6;

for i=1:n_scenarios
    n_2100=0;
    m_2100=0;
    n_max=0;
    m_max=0;
    for j=1:n_runs
        if m_cum_rnd_IPCC(end,i,j)>= m_target_1_5(j,1);
            n_2100=n_2100+1;
        end
        if m_cum_rnd_IPCC(end,i,j)>= m_target_2(j,1);
            m_2100=m_2100+1;
        end
        if m_cum_rnd_IPCC(i_t_max(1,i),i,j)>= m_target_1_5(j,1);
            n_max=n_max+1;
        end
        if m_cum_rnd_IPCC(i_t_max(1,i),i,j)>= m_target_2(j,1);
            m_max=m_max+1;
        end
    end
    P_v_1_5_2100(i) = n_2100/n_runs;
    P_v_2_2100(i) = m_2100/n_runs;
    P_v_1_5_max(i) = n_max/n_runs;
    P_v_2_max(i) = m_max/n_runs;
end
      
    
%%
toc