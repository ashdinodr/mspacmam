% This is the file used to store the constants and fixed parameters
% for the MSPACMAM program.

% Constants
secInDay = 86400; % number of seconds in day
secInYr = 31536000; % number of seconds in year
g = 9.81; % gravitational acceleration (m/s^2)

% Fixed parameters (see Dinauer et al. [2022] for more details)
rho_om = 1.06; % density of particulate organic matter (g/cm^3)
rho_caco3 = 2.71; % density of CaCO3 (g/cm^3)
rho_opal = 2.10; % density of opal (g/cm^3)
MWc = 12; % molecular weight of carbon (g/mol)
MWcaco3 = 100; % molecular weight of CaCO3 (g/mol)
MWsio2 = 60; % molecular weight of opal (g/mol)
alpha = 2.7; % mass ratio of organic matter to carbon (g POM/g POC)
k_calcUp = 10^(-14.3); % calcite dissolution rate constant (mol/cm^2/s) for 0.8<omega<1
n_calcUp = 0.11; % calcite dissolution reaction order for 0.8<omega<1
k_calcLow = 10^(-10.8); % calcite dissolution rate constant (mol/cm^2/s) for omega<=0.8
n_calcLow = 4.7; % calcite dissolution reaction order for omega<=0.8
k_arag = 0.013; % aragonite dissolution rate constant (g/g/d)
n_arag = 1.37; % aragonite dissolution reaction order
SSA_calcSm = 1040; % SSA of small particle calcite (m^2/mol)
SSA_calcLg = 430; % SSA of large particle calcite (m^2/mol)
SSA_arag = 217; % SSA of aragonite (m^2/mol)