function [theta_i_eff, theta_v_eff, raa_eff] =  local_viewing_angle(theta_i, phi_i, theta_v, phi_v, slope, aspect)
%     function based on snowoptics from Ghislain Picard and Tuzet F
%     compute the effective viewing angle of incident and observer  as well as the relative azimuth
%     angle between incident and observer for a given slope according to dumont et al.2011
%     :param theta_i: Incident zenith angle (radians)
%     :param phi_i: Incident azimuth angle (radians)
%     :param theta_v: Observer zenith angle (radians)
%     :param phi_v: Observer azimuth angle (radians)
%     :param slope: Slope inclination (radians)
%     :param aspect: Slope aspect (radians)

%convert degrees to radians
theta_i = deg2rad(theta_i);
theta_v = deg2rad(theta_v);
phi_i = deg2rad(phi_i);
phi_v = deg2rad(phi_v);
slope = deg2rad(slope);
aspect = deg2rad(aspect);

% Local incident zenith angle
mu_i = cos(theta_i) * cos(slope) + sin(theta_i) *  sin(slope) * cos(phi_i - aspect);
    if mu_i < 0.000001  % Grazing angle, unstable
       mu_i = nan;
    end

% Local viewing zenith angle
    mu_v = cos(theta_v) * cos(slope) + sin(theta_v) * sin(slope) * cos(phi_v - aspect);

    theta_i_eff = acos(mu_i);
    theta_v_eff = acos(mu_v);
   
% % Remove part of the polar representation that correspond to an observer behind the slope
    if theta_v_eff > deg2rad(90)
        theta_v_eff = nan;
    end
% Local relative azimuth angle (dumont et al.2011)
    mu_az_numerator = (cos(theta_v) * cos(theta_i) + sin(theta_v) * sin(theta_i) * cos(phi_v-phi_i)- mu_i * mu_v);
    mu_az_denominator = sin(theta_i_eff) * sin(theta_v_eff);

% When illumination or observator is at nadir (in the new reference), set RAA to zero
if mu_az_denominator == 0
    mu_az = 0;
else
    mu_az =  mu_az_numerator /  mu_az_denominator; 
end

if mu_az < -1
    mu_az = -1;
elseif mu_az > 1
    mu_az = 1;
end

raa_eff = acos(mu_az);

%convert radians to degrees for output
theta_i_eff = rad2deg(theta_i_eff);
theta_v_eff = rad2deg(theta_v_eff);
raa_eff = rad2deg(raa_eff);
  