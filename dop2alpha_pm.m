%%%%% Phase Matching for AR2018 %%%%%
%%%%% Eric Wang, Oct 2018       %%%%%
% INPUT:
% (i)  file_name: Output excess phase file from nret folder
%                /ags/projects/hiaper/2018.023_ar2018/nret/2018.026_rf01/sept
%                prn_gXXx_gXX_XXXXXX_phase_ar18_rf01_v2_shift_v2_cor.txt 
% (ii) nrec: in-situ refractivity reading (N-unit) at the aircraft
% (iii)amin: specify the impact parameter range (km) -> from amin to amin+20
%
% OUTPUT:
% (i)  Output_alpha_neg_%8s.txt: Bending angle profile
% (ii) Output_v_neg_%8s.txt: v profile -> only for debugging





% clear all; close all; clc;
% 
% addpath('C:\Work\Research\Jennifer_misc\ar_2018\data\');
% file_ind = 0;
% 
% switch(file_ind)
%     case 0
%         file_name = 'prn_g09s_g28_520208_phase_ar18_rf01_v2_shift_v2_cor.txt';
%         nrec = 49.7201;
%         amin = 6370;
%     case 1
%         file_name = 'prn_g11s_g23_500808_phase_ar18_rf01_v2_shift_v2_cor.txt';
%         nrec = 62.4017;
%         amin = 6370;
%     case 2
%         file_name = 'prn_g14s_g23_500808_phase_ar18_rf01_v2_shift_v2_cor.txt';
%         nrec = 61.6354;
%         amin = 6390;
%     case 3
%         file_name = 'prn_g22s_g09_505818_phase_ar18_rf01_v2_shift_v2_cor.txt';
%         nrec = 60.3737;
%         amin = 6370;
%     case 4
%         file_name = 'prn_g26s_g07_508553_phase_ar18_rf01_v2_shift_v2_cor.txt';
%         nrec = 60.2080;
%         amin = 6370;
%     case 5
%         file_name = 'prn_g31s_g23_503284_phase_ar18_rf01_v2_shift_v2_cor.txt';
%         nrec = 60.4768;
%         amin = 6380;
%     case 6
%         file_name = 'prn_g05r_g09_509408_phase_ar18_rf01_v2_shift_v2_cor.txt';
%         nrec = 60.2080;
%         amin = 6380;
%     case 7
%         file_name = 'prn_g30r_g23_507216_phase_ar18_rf01_v2_shift_v2_cor.txt';
%         nrec = 60.2080;
%         amin = 6380;
% end



function dop2alpha_pm(file_name, nrec, amin)

k = (1.57542e9 / 299792.458) * 2*pi;  % rad/km
% wavenumber = f/c

[time, x_rec, y_rec, z_rec, u_rec, v_rec, w_rec, x_gps, y_gps, z_gps, u_gps, v_gps, w_gps, ex_ph, loss, ex_dop] = read_ar_v1(file_name);
r_gps_rec = sqrt((x_gps-x_rec).^2 + (y_gps-y_rec).^2 + (z_gps-z_rec).^2); % vector
r_rec     = sqrt(x_rec.^2 + y_rec.^2 + z_rec.^2);
r_gps     = sqrt(x_gps.^2 + y_gps.^2 + z_gps.^2);  % km

theta = acos( (x_gps.*x_rec + y_gps.*y_rec + z_gps.*z_rec) ./ r_gps ./ r_rec );  % rad
dist  = sqrt( (x_gps-x_rec).^2 + (y_gps-y_rec).^2 + (z_gps-z_rec).^2); %vs r_gps_rec ??
phi_geo = k * dist;                        % rad  (geometric phase)
ex_ph = cumsum(ex_dop);                 % Bing didn't correct ex_ph, only ex_dop
% ex_ph = -1*np.cumsum(ex_dop+0.05)   % case 6

% figure();
% plot(time, ex_dop);

% Calculate the elevation angle
elev = zeros(length(time),1);
for i=1:1:length(elev) 
    elev(i) = pi/2-acos(dot([x_gps(i)-x_rec(i), y_gps(i)-y_rec(i), z_gps(i)-z_rec(i)],[x_rec(i), y_rec(i), z_rec(i)])/r_gps_rec(i)/r_rec(i));
end
% figure();
% plot(time,elev);


% Works for both setting and rising
i_0 = 1;
while elev(i_0)>0
    i_0 = i_0 + 1;
end



if file_name(8)=='s'  % Setting negative
    time  = time(i_0:end);
    theta = theta(i_0:end);
    r_gps = r_gps(i_0:end);
    r_rec = r_rec(i_0:end);
    phi   = phi_geo + ex_ph*1e-3*k;          % rad  (total phase)
    phi   = phi(i_0:end);
    
else    % Rising negative
    time  = flip(time(i_0:end));
    theta = flip(theta(i_0:end));
    r_gps = flip(r_gps(i_0:end));
    r_rec = flip(r_rec(i_0:end));
    phi   = phi_geo - ex_ph*1e-3*k;          % rad  (total phase)
    phi   = flip(phi(i_0:end));
    
end


% amin = 6370
amax = amin+20;
a = amin:1e-3:amax;
[a_out, alpha, v] = pm_v2_air(time, ones(length(elev)-i_0+1,1), a, theta, r_gps, r_rec, phi, k, 1+nrec*1e-6, 1);
% norp=1 => negative elevation,   norp=0 => positive elevation
% a_out, alpha, v = pm.pm_v2_air(time[1500:i_0], np.ones(i_0-1500), a, theta[1500:i_0], r_gps[1500:i_0], r_rec[1500:i_0], phi[1500:i_0], k, 1+nrec*1e-6, 0)
% a_out, alpha, v = pm.pm_v2_air(time[i_0:], np.ones(len(time)-i_0), a, theta[i_0:], r_gps[i_0:], r_rec[i_0:], phi[i_0:], k, 1+nrec*1e-6, 1)


% figure();
% plot(alpha, a_out(2:end));


output_alpha = sprintf('Output_alpha_neg_%8s.txt',file_name(5:12));
output_v     = sprintf('Output_v_neg_%8s.txt',file_name(5:12));
file_alpha = fopen(output_alpha,'w');
file_v     = fopen(output_v,'w');
fprintf(file_alpha,'%8.3f  %12.8f\n',[a_out(2:end); alpha']);
fprintf(file_v,'%8.3f  %12.8f  %12.8f\n',[a_out; real(v'); imag(v')]);
fclose(file_alpha);
fclose(file_v);
