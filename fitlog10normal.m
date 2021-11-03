function [mu, sigma,r2] = fitlog10normal(probability_data,data,test_runs, test_range)

n = size(probability_data,2);

mu_init = log10(data(3));
sigma_init = log10(data(end))-mu_init;

mu_test = ((mu_init*(1-test_range)):(2*test_range*mu_init/(test_runs-1)):(mu_init*(1+test_range)))';
sigma_test = (sigma_init*(1-test_range)):(2*test_range*sigma_init/(test_runs-1)):(sigma_init*(1+test_range));

for i=1:size(probability_data,2)
    for j=1:test_runs
        % error not normalized
        Squared_error(:,j,i) = ((10.^(mu_test + probability_data(i)*sigma_test(j))-data(i))).^2;
    end
end
Sum_squared_error = sum(Squared_error,3);
SS_tot = sum((data-mean(data,2)).^2,2);
r2_fit = 1-Sum_squared_error ./ SS_tot; 
[min_error,n_min_mu] = min(Sum_squared_error,[],1);
[~,n_min_sigma] = min(min_error,[],2);

mu = mu_test(n_min_mu(n_min_sigma));
sigma = sigma_test(n_min_sigma);
r2 = r2_fit(n_min_mu(n_min_sigma),n_min_sigma);
