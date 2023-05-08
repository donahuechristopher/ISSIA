%%  """Calculate the r0 of the BRF.
%        See brf function for details.
%        theta_i: illumination zenith angle
%        :param theta_v: viewing zenith angle
%        :param phi: relative azimuth angle (illumination - viewing)
%        :return: r0

function r0 = brf0(theta_i, theta_v, phi)
new_phi=pi-phi; 
theta = rad2deg(acos(-cos(theta_i) * cos(theta_v) + sin(theta_i) * sin(theta_v) * cos(new_phi))  );

%   theta = acos(cos(theta_i) * cos(theta_v) + sin(theta_i) * sin(theta_v) * cos(phi)) * 180 / pi;
phase = 11.1 * exp(-0.087 * (theta)) + 1.1 * exp(-0.014 *(theta));
rr = 1.247 + 1.186 * (cos(theta_i) + cos(theta_v)) + 5.157 * (cos(theta_i) * cos(theta_v)) + phase;
rr = rr / (4 * (cos(theta_i) + cos(theta_v)));
r0=rr;

