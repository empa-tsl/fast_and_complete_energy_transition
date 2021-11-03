% Maximum historic growth for RE
% Fit to logistic curve

clear variables
close all

tic

%% loading and preparing data
% time series for energy demand and CO2 emissions

dt = 0.01; %a

[~,labels,~] = xlsread('CO2 emissions and remaining budget.xlsx','RE','B2:F2');

time_frame = xlsread('CO2 emissions and remaining budget.xlsx','RE','A6:A60'); %a
t_ext = 1980:1:2100;

P_RE = xlsread('CO2 emissions and remaining budget.xlsx','RE','K6:N60'); %TW
P_RE (isnan(P_RE))=0;
P_RE(:,end+1) = sum(P_RE,2);

ATP_RE = xlsread('CO2 emissions and remaining budget.xlsx','RE','K2:N2'); %TW
ATP_RE(:,end+1) = sum(ATP_RE,2);
n_r = size(ATP_RE,2);

% IEA scenarios
time_frame_scenario = xlsread('CO2 emissions and remaining budget.xlsx','RE','A66:A69'); %a
P_RE_SPS = xlsread('CO2 emissions and remaining budget.xlsx','RE','K66:N69'); %TW
P_RE_SDS = xlsread('CO2 emissions and remaining budget.xlsx','RE','O66:R69'); %TW

%% Define functions and perform non-linear regression

%logistic_curve = @(beta,t) ATP_PV./(1+exp(-beta(1).*(time_frame-beta(2))));
logistic_curve_PV = @(b_log,t) ATP_RE(1,1)./(1+exp(-b_log(1).*(t-b_log(2))));
logistic_curve_wind = @(b_log,t) ATP_RE(1,2)./(1+exp(-b_log(1).*(t-b_log(2))));
logistic_curve_others = @(b_log,t) ATP_RE(1,3)./(1+exp(-b_log(1).*(t-b_log(2))));
logistic_curve_hydro = @(b_log,t) ATP_RE(1,4)./(1+exp(-b_log(1).*(t-b_log(2))));
logistic_curve_RE = @(b_log,t) ATP_RE(1,5)./(1+exp(-b_log(1).*(t-b_log(2))));
exponential_curve = @(b,t) b(2) .* exp(b(3).*(t-b(1)));


PV_fit = fitnlm(time_frame,P_RE(:,1), logistic_curve_PV,[1 2020]);
wind_fit = fitnlm(time_frame,P_RE(:,2), logistic_curve_wind,[1 2020]);
others_fit = fitnlm(time_frame,P_RE(:,3), logistic_curve_others,[1 2020]);
hydro_fit = fitnlm(time_frame,P_RE(:,4), logistic_curve_hydro,[0.1 1960]);
RE_fit = fitnlm(time_frame,P_RE(:,5), logistic_curve_RE,[0.1 2020]);
b_log_estimates(1,:) = PV_fit.Coefficients.Estimate;
R2_PV(1,:) = PV_fit.Rsquared.Adjusted;
b_log_estimates(2,:) = wind_fit.Coefficients.Estimate;
R2_PV(2,:) = wind_fit.Rsquared.Adjusted;
b_log_estimates(3,:) = others_fit.Coefficients.Estimate;
R2_PV(3,:) = others_fit.Rsquared.Adjusted;
b_log_estimates(4,:) = hydro_fit.Coefficients.Estimate;
R2_PV(4,:) = hydro_fit.Rsquared.Adjusted;
b_log_estimates(5,:) = RE_fit.Coefficients.Estimate;
R2_PV(5,:) = RE_fit.Rsquared.Adjusted;

b_start = [1980 10^-7 0.25; 1980 10^-7 0.25; 1980 10^-7 0.25; 1965 0.1 0.1; 1965 0.1 0.25];
for r=1:n_r
PV_fit_exp = fitnlm(time_frame,P_RE(:,r), exponential_curve,b_start(r,:));
b_estimates_exp(r,:) = PV_fit_exp.Coefficients.Estimate;
R2_exp(r,:) = PV_fit_exp.Rsquared.Adjusted;
end

%% Preparing data for display

time_display = 1960:dt:2100;


logistic_curve_fits(:,1) = logistic_curve_PV(b_log_estimates(1,:),time_display);
logistic_curve_fits(:,2) = logistic_curve_wind(b_log_estimates(2,:),time_display);
logistic_curve_fits(:,3) = logistic_curve_others(b_log_estimates(3,:),time_display);
logistic_curve_fits(:,4) = logistic_curve_hydro(b_log_estimates(4,:),time_display);
logistic_curve_fits(:,5) = logistic_curve_RE(b_log_estimates(5,:),time_display);

for r=1:n_r
    exp_curve_fits(:,r) = exponential_curve(b_estimates_exp(r,:),time_display);
end

%%
toc
