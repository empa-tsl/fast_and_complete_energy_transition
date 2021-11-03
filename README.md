# fast_and_complete_energy_transition

Calculation code for
"Reducing climate risks with fast and complete energy transitions: applying the precautionary principle to the Paris agreement"
Harald Desing and Rolf Widmer

Matlab R2018a code, developed and tested under Windows 10, Version 1809.

Fill the empty green cells in the excel table "CO2 emissions and remaining budget.xlsx" (these data are protected by copyright, so they cannot be shared directly).

Each code file (.m) can be run independently, except "fitlog10normal.m", which is a customized function of fitting data points with known probability to a log-normal distribution and required in the other code files. It takes the data from the two excel files provided and saves the results in the workspace. Runtime ranges, depending on the code, from few seconds to minutes. Output are matrices containing the results for each simulation run and time step.
