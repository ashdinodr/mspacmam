% This is the config file used to set the parameters and settings
% for the MSPACMAM program.

% Model experiment name, which is used as name of directory in which output is saved
exp_name = 'some_short_expname'; % experiment name, with words separated by underscores

% Horizontal grid points at which model is run
my_coords = 1; % flag specifying whether to read coordinates from 'my_ocean_gridpoints.txt' (1)
               % or use all available coordinates within the ocean model domain (0)

% Model configuration (see Dinauer et al. [2022] for more details)
protect_flag = 0; % flag specifying whether to protect fraction of POC (1) or not (0)
fragment_flag = 0; % flag specifying whether to fragment large particles (1) or not (0)

% Model input parameters (see Dinauer et al. [2022] for more details)
r_sm = 15.625; % small particle radius (micron)
r_lg = 250; % large particle radius (micron)
phi_sm = 0; % small particle porosity (fraction)
phi_lg = 0.97; % large particle porosity (fraction)
k_POC = 0.122; % POC remineralization rate constant (1/d)
Tref = 16.47; % reference temperature for k_POC (degC)
Q10 = 2.293; % remineralization rate increase per 10degC increase (unitless)
K_O2 = 8.25; % half-saturation constant for O2 uptake (micromol/L)
kappa_lg = 1; % proportionality constant for export partitioning (fraction)

% Numerical integration specs
dt = (137.5-112.5)/400; % time step size (d)
runlen = ((5500-100)/10)*2; % time span of integration (d)
drift_thresh = 0.1; % equilibrium drift criteria (percent)
integ_method = 2; % flag specifying whether to use first-order forward Euler method (1)
                  % or second-order Runge-Kutta (RK2) midpoint method (2)
                  % for numerically integrating (solving) the differential equations