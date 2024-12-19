function [O2, O2eq] = O2_SVU(calphase, temp, psal, pres, lon, lat, SVU_coef, S0)
%DEFINITION OF OXYGEN CALCULATION FUNCTION
% function [O2] = calculate_O2(P0, PC, Ksv, Scorr)
    % Calculate O2 concentration using the provided formula
    % Formula: [O2] = ((P0 / PC) - 1) / Ksv * Scorr
    %O2 = ((P0 / PC) - 1) / Ksv * Scorr;
% end

% Args:
%   P0: Phase offset (calibration coefficient)
%   PC: Calibrated phase factor
%   Ksv: Temperature-dependent calibration factor for O2 sensitivity
%   Scorr: Salinity correction factor
%   calphase: Calibrated phase
%   temp: Temperature (degC)
%   psal: Practical Salinity
%   pres: Pressure (dbar)
%   lon: Longitude (degE)
%   lat: Latitude (degN)
%   SVU_coef: Coefficients (SVUCoef0 ... SVUCoef6)
%   S0: Internally set salinity value (default: 0)

% nargin: Number of input arguments passed to the function
% Check if the optional parameter S0 is provided
  
 if nargin < 8
      S0 = 0; % Default value for S0 if not provided
 end

% Constants 
K0 = 273.15;
% insert manufacturer coefficients for each optode
o2coef_b = [-6.24097E-3, -6.93498E-3, -6.90358E-3, -4.29155E-3];
o2coef_c0 = -3.11680E-7;

% Physical variables using GSW toolbox (Gibbs SeaWater (GSW) Oceanographic Toolbox)
SA = gsw_SA_from_SP(psal, pres, lon, lat); % Calculates Absolute Salinity from Practical Salinity, pressure, longitude, and latitude.
CT = gsw_CT_from_t(SA, temp, pres); % Converts in-situ temperature to Conservative Temperature using Absolute Salinity and pressure.
O2eq = gsw_O2sol(SA, CT, pres, lon, lat); % Computes equilibrium oxygen solubility in seawater.rho0 = 1000 + gsw_sigma0(SA, CT); % Calculates in-situ density anomaly from Absolute Salinity and Conservative Temperature.

%SVU calibration
Ksv = SVU_coef(1) + SVU_coef(2)*temp + SVU_coef(3)*temp^2; % Temperature-dependent calibration factor for O2 sensitivity. 
P0 = SVU_coef(4) + SVU_coef(5)*temp; % Temperature-dependent phase offset.
PC = SVU_coef(6) + SVU_coef(7) * calphase; % Calibration phase factor dependent on input calibrated phase.
optode_uM = ((P0 / PC) - 1) / Ksv; % Calculated oxygen concentration in micromolar units based on calibration.

% Scaled temperature
temps = log((K0 + 25 - temp) / (K0 + temp)); % Computes a scaled logarithmic temperature used in salinity correction.

% Salinity correction
Scorr = exp((psal - S0) * polyval(o2coef_b, temps) + o2coef_c0 * (psal^2 - S0^2)); % Adjusts oxygen concentration for salinity differences.

% Oxygen concentration calculation
O2 = 1000 * Scorr * optode_uM / rho0; % Final O2 concentration in umol/kg, normalized by density.

end