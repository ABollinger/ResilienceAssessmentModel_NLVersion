%% LINKS

% http://www.energy.ca.gov/2012publications/CEC-500-2012-057/CEC-500-2012-057.pdf
% https://github.com/stevenblair/ieee738matlab
% IEEE738: IEEE Standard for Calculating the Current-Temperature of Bare Overhead Conductors
% http://www.pdc-cables.com/oh_limits_powerflow.pdf
% http://mitei.mit.edu/system/files/Electric_Grid_2_Enhancing_Transmission_Network_System_Operations.pdf
% http://www.wire1002.ch/fileadmin/user_upload/Major_events/WS_Nice_2011/Spec._presentations/hosek.pdf
% http://itee.uq.edu.au/~aupec/aupec06/htdocs/content/pdf/206.pdf
% http://www.waset.org/journals/waset/v61/v61-67.pdf


%% SETTINGS

R_T_high = 8.688e-5;  % conductor unit resistance at high temperature reference (ohm/m)
R_T_low = 7.283e-5;   % conductor unit resistance at low temperature reference (ohm/m)
T_high = 75.0;        % high temperature reference (degrees)
T_low = 25.0;         % low temperature reference (degrees)
Ta = 35.0;            % ambient temperature (degrees Celcius)
rho_f = 1.029;        % density of air (kg/m^3)
D = 22.8;             % conductor diameter (mm)
epsilon = 0.5;        % emissivity
alpha = 0.5;          % solar absorptivity

H_e = 25;             % elevation of conductor above sea level (m)
H_c = 72.5;           % altitude of sun (degrees)
Q_s = -42.2391 + 63.8044*H_c - 1.9220*H_c^2 + 3.46921e-2*H_c^3 - 3.61118e-4*H_c^4 + 1.94318e-6*H_c^5 - 4.07608e-9*H_c^6;
K_solar = 1 + 1.148e-4*H_e - 1.108e-8*H_e^2;
Q_se = K_solar * Q_s;
Z_c = 139;            % azimuth of sun (degrees)
Z_l = 90.0;           % azimuth of electrical line (degrees); 90 degrees (or 270 degrees) for east-west
theta = acos(cos(H_c) * cos(Z_c - Z_l));
area = D/1000;        % projected area of conducter per unit length (m^2/m)

Tc_test = 100;

IT = real(getI(Tc_test, R_T_high, R_T_low, T_high, T_low, Ta, rho_f, D, epsilon, alpha, Q_se, theta, area));
K = sqrt(3) * IT * 380 * 0.95

