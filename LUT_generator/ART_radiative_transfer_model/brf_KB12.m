%% Snow BRF function from Kokhanovsky et al. 
% Described in DOI:10.1109/LgrS.2012.2185775. 
% Calculate BRF value for given illumination / viewing angles, wavelengths and SSA.
% Code adapted from on REDRESS(python) - M. Lamare, M. Dumont, G. Picard 
%   
% Parameters
% theta_i: illumination zenith angle
% theta_v: viewing zenith angle
% phi: relative azimuth (illumination - viewing)
% SSA: Specific surface area of snow
% b: 13, as L = 13d in kokhanovsky's paper

function [BRF albedo wvl ]= brf_KB12(theta_i, theta_v, phi, opt_radius, M, n, ni, wvl_specim)

b = 13; % b: 13, as L = 13d in kokhanovsky's paper

%convert opt_radius from microns to meters
opt_radius =opt_radius*10^-6;

%convert inputs from degree to radians
theta_i=deg2rad(theta_i);
theta_v=deg2rad(theta_v);
phi=deg2rad(phi);

% r0 in kokhanovsky's paper
r = brf0(theta_i, theta_v, phi);

% k0 for theta_v and theta i
k0v = (3 / 7) * (1 + 2 * cos(theta_v));
k0i = (3 / 7) * (1 + 2 * cos(theta_i));

for i=1:length(n)
    gamma(i) = 4*pi * (ni(i) + M) / (wvl_specim(i));
    alpha(i) = sqrt(gamma(i) * b * 2* opt_radius);

    % r(theta_i, theta_v, phi)
    rr(i) = r * exp(-alpha(i) * k0i * k0v / r);


    %calculate albedo
    Q= (k0i * k0v) / r;
    A(i) = exp(-1/Q * log(r/rr(i))); % Kokhanovsky & Schreier 2009

end

BRF = rr;
albedo=A;