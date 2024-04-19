%%%%% Phase Matching for AR2018 %%%%%
%%%%% Eric Wang, Oct 2018       %%%%%
%
% Plots the bending angle profile using phase matching
%
%
%% INPUT:
% (i)  file_name: Output excess phase file from nret folder
%                /ags/projects/hiaper/2018.023_ar2018/nret/2018.026_rf01/sept
%                prn_gXXx_gXX_XXXXXX_phase_ar18_rf01_v2_shift_v2_cor.txt 
% (ii) nrec: in-situ refractivity reading (N-unit) at the aircraft
% (iii)amin: specify the impact parameter range (km) -> from amin to amin+20
% (iv) norp: negative or positive elevation angle. positive: 0, negative: 1
% (v)  plt_bad: a switch. If true, plots all the alpha(a) data. If false,
% excludes vertical data at top and all points below the maximum alpha
%
%% OUTPUT:
% (i)  Output_alpha_neg_%8s.txt: Bending angle profile
% (ii) Output_v_neg_%8s.txt: v profile -> only for debugging
% 
%% Dependencies: 
%   pm_v4_air.m
%   read_ar_v1_test.m
%




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



function dop2alpha_pm_v4(file_name, nrec, amin, norp, plt_bad)

k = (1.57542e9 / 299792.458) * 2*pi;  % rad/km
std_cutoff = 1*10^-5;
% wavenumber = f/c1

%plt_bad is a switch. If true, plots all the alpha(a) data. If false,
%excludes vertical data at top and all points below the maximum alpha

%set default value for plt_bad
if nargin < 5
    plt_bad = true;
end

tic

%retrieve time, position, velocity, and phase info
[time, x_rec, y_rec, z_rec, u_rec, v_rec, w_rec, x_gps, y_gps, z_gps, u_gps, v_gps, w_gps, ex_ph, loss, ex_dop] = read_ar_v1_test(file_name);
r_gps_rec = sqrt((x_gps-x_rec).^2 + (y_gps-y_rec).^2 + (z_gps-z_rec).^2); % vector
r_rec     = sqrt(x_rec.^2 + y_rec.^2 + z_rec.^2);
r_gps     = sqrt(x_gps.^2 + y_gps.^2 + z_gps.^2);  % km

theta = acos( (x_gps.*x_rec + y_gps.*y_rec + z_gps.*z_rec) ./ r_gps ./ r_rec );  % rad
dist  = sqrt( (x_gps-x_rec).^2 + (y_gps-y_rec).^2 + (z_gps-z_rec).^2); %vs r_gps_rec ??
phi_geo = k * dist;                        % rad  (geometric phase)
ex_ph = cumsum(ex_dop);                 % Bing didn't correct ex_ph, only ex_dop
% ex_ph = -1*np.cumsum(ex_dop+0.05)   % case 6

% Calculate the elevation angle
elev = zeros(length(time),1);
for i=1:1:length(elev) 
    elev(i) = pi/2-acos(dot([x_gps(i)-x_rec(i), y_gps(i)-y_rec(i), z_gps(i)-z_rec(i)],[x_rec(i), y_rec(i), z_rec(i)])/r_gps_rec(i)/r_rec(i));
end

disp(["first elev", elev(1)])
disp(["first time", time(1)])

%% Modifications to Eric's code - Kate
rors = file_name(8); %rising or setting

if rors == 'r' %rising
    disp("RISING")
    i_0 = length(elev);
    disp(elev(i_0))
    disp("FIRST")
    while elev(i_0)>0
        i_0 = i_0 - 1;
    end
    % now i_0 and before is negative, after i_0 is positive


else           %setting
   
    i_0 = 1;
    while elev(i_0)>0
        i_0 = i_0 + 1;
    end
    % now i_0 and after is negative, before i_0 is positive

end

%% Eric's original code
% Works for both setting and rising
% i_0 = 1;
% while elev(i_0)>0
%     i_0 = i_0 + 1;
% end

%  
% 
% if file_name(8)=='s'  % Setting negative
%     time  = time(i_0:end);
%     theta = theta(i_0:end);
%     r_gps = r_gps(i_0:end);
%     r_rec = r_rec(i_0:end);
%     phi   = phi_geo + ex_ph*1e-3*k;          % rad  (total phase)
%     phi   = phi(i_0:end);
%     
% else    % Rising negative
%     %the original code
%     time  = flip(time(i_0:end));
%     theta = flip(theta(i_0:end));
%     r_gps = flip(r_gps(i_0:end));
%     r_rec = flip(r_rec(i_0:end));
%     phi   = phi_geo - ex_ph*1e-3*k;          % rad  (total phase)
%     phi   = flip(phi(i_0:end));
% 
% end
% 

%% Modifications to Eric's code - Kate
if norp == 1 %negative
    if rors == 'r' %rising
        %normal

%         ind_range = [1:i_0];
%         phi   = phi_geo + ex_ph*1e-3*k; % rad  (total phase)
        phi   = phi_geo - ex_ph*1e-3*k; %rising - ?
        
        time  = time(1:i_0);
        theta = theta(1:i_0);
        r_gps = r_gps(1:i_0);
        r_rec = r_rec(1:i_0);
        phi   = phi(1:i_0);

%         time  = time(i_0:end);
%         theta = theta(i_0:end);
%         r_gps = r_gps(i_0:end);
%         r_rec = r_rec(i_0:end);
%         phi   = phi(i_0:end);

%         disp("NEG RISING")
%         scatter(time,elev(1:i_0));
%         legend("all", "negative rising portion");
%         hold off;

%         save("kate.mat","time", "theta")

    else           %setting
        %flip
%         ind_range = [i_0:end]
        disp("NEGATIVE SETTING")
%         phi   = phi_geo - ex_ph*1e-3*k;          % rad  (total phase)
        
        %this works with eric/kate i_0
        phi   = phi_geo + ex_ph*1e-3*k; %setting + ?
        time  = time(i_0:end);
        theta = theta(i_0:end);
        r_gps = r_gps(i_0:end);
        r_rec = r_rec(i_0:end);
        phi   = phi(i_0:end);

%         time  = time(1:i_0);
%         theta = theta(1:i_0);
%         r_gps = r_gps(1:i_0);
%         r_rec = r_rec(1:i_0);
%         phi   = phi(1:i_0);
        
        disp(["i_0 here", i_0])
        disp(["length of phi", length(phi)])


    end


else         %positive
    if rors == 'r' %rising
        %flip
%         ind_range = [i_0+1:end];
        phi   = phi_geo - ex_ph*1e-3*k;          % rad  (total phase)
        %maybe works????? - if flipped in read_ar, kate i_0
        time  = time(i_0+1:end);
        theta = theta(i_0+1:end);
        r_gps = r_gps(i_0+1:end);
        r_rec = r_rec(i_0+1:end);
        phi   = phi(i_0+1:end);
        disp("POS RISING")

    else           %setting
        %normal
%         ind_range = [1:i_0-1];
        phi   = phi_geo + ex_ph*1e-3*k;          % rad  (total phase)

        %seems to work with eric/kate i_0
%         time  = flip(time(1:i_0-1));
%         theta = flip(theta(1:i_0-1));
%         r_gps = flip(r_gps(1:i_0-1));
%         r_rec = flip(r_rec(1:i_0-1));
%         phi   = flip(phi(1:i_0-1));

        %flip does not seem to matter
        time  = time(1:i_0-1);
        theta = theta(1:i_0-1);
        r_gps = r_gps(1:i_0-1);
        r_rec = r_rec(1:i_0-1);
        phi   = phi(1:i_0-1);

        disp("POS SETTING")
    end
end

% plot(time, theta) %this looks good
% title("Theta vs Time")
% % 
% % plot(r_gps, r_rec) %this looks good
% % title("r's")
% 
% disp(["i_0", i_0])
% disp(["elev length",length(elev)])
% disp(["time length", length(time)])

% plot(phi)

%% End modifications to Eric's code - Kate

% amin = 6370
amax = amin+20;
a = amin:1e-3:amax;

% set amplitude equal to one
% amp_1 = ones(length(elev)-i_0+1,1);
amp_1 = ones(length(time),1); %Kate

% convert refractivity to refractive index
n_i_rec = 1+nrec*1e-6;

% save("dop2alp_midfile_var");

% [a_out, alpha, v] = pm_v3_air(time, ones(length(elev)-i_0+1,1), a, theta, r_gps, r_rec, phi, k, 1+nrec*1e-6, 1);
[a_out, alpha, v] = pm_v4_air(time, amp_1, a, theta, r_gps, r_rec, phi, k, n_i_rec, norp);
% norp=1 => negative elevation,   norp=0 => positive elevation
% a_out, alpha, v = pm.pm_v2_air(time[1500:i_0], np.ones(i_0-1500), a, theta[1500:i_0], r_gps[1500:i_0], r_rec[1500:i_0], phi[1500:i_0], k, 1+nrec*1e-6, 0)
% a_out, alpha, v = pm.pm_v2_air(time[i_0:], np.ones(len(time)-i_0), a, theta[i_0:], r_gps[i_0:], r_rec[i_0:], phi[i_0:], k, 1+nrec*1e-6, 1)

split_name = split(file_name,'/');
name_segment = char(split_name(end));

if ~ plt_bad

    %negative section
    if norp == 1
        %exclude top portion that is roughly vertical
    
        % loop through array looking for three adjacent 0 values - will take
        % this as the start of the waste at the top
        found_start = false;
        zero_ind_array = find(alpha == 0);
        array_ind = 1;
    
        while ~ found_start
    
            if zero_ind_array(array_ind+1) == zero_ind_array(array_ind)+1
                if zero_ind_array(array_ind+2) == zero_ind_array(array_ind)+2
                    found_start = true;
                else
                    array_ind = array_ind + 3;
                end
            else
                array_ind = array_ind + 2;
            end

        end
        
        zero_ind = zero_ind_array(array_ind);
    
        %exlcude bottom waste %always 0?
    
        %find alpha maximum
        maximum = max(alpha);
        ind = find(alpha == maximum);
    
        figure();
%         plot(alpha, a_out(2:end));
%         hold on;
    
        %bottom at beginning, top at end (for setting and rising)
        alpha_new = alpha(ind:zero_ind-1);
        a_out_new = a_out(ind:zero_ind-1);
% 
%         disp(["alpha length",length(alpha)])
%         disp(["a length", length(a)])

        
        plot(alpha_new, a_out_new);
    
    %positive section    
    else 
        
        differences = diff(alpha);
        diff_std = std(differences);

        if diff_std < std_cutoff
        
        
                % exclude beginning and end where slope of a changes dramatically
                % maybe just the first and second???
        
                % loop forwards and backwards through array starting at middle 
                % point looking for three adjacent 0 values - will take
                % this as the start of the waste at the top
                found_start = false;
                found_end = false;
                zero_ind_array = find(alpha == 0);
                array_ind = length(zero_ind_array);
                halfway = fix(length(alpha)/2);
        
        
                %%% for "normal" < 10**5, end alpha at zero_ind_array(1)
        
        
                %start at beginning and go backwards towards beginning
        %         while ~ found_start
        %             disp(array_ind)
        %             if zero_ind_array(array_ind-1) == zero_ind_array(array_ind)-1
        %                 if zero_ind_array(array_ind-2) == zero_ind_array(array_ind)-2
        %                     found_start = true;
        %                 else
        %                     array_ind = array_ind - 3;
        %                 end
        %             else
        %                 array_ind = array_ind - 2;
        %             end
        %         end
        
                %change to just look for first <= 0 - alpha should not be 0
            
        %         while ~ found_start
        %             disp(array_ind)
        %             if zero_ind_array(array_ind-1) == zero_ind_array(array_ind)-1
        %                 if zero_ind_array(array_ind-2) == zero_ind_array(array_ind)-2
        %                     found_start = true;
        %                 else
        %                     array_ind = array_ind - 3;
        %                 end
        %             else
        %                 array_ind = array_ind - 2;
        %             end
        %         end
                
                zero_ind_start = zero_ind_array(array_ind);
                array_ind = fix(length(zero_ind_array)/2);
        
                %start at middle and move towards end
        %         while ~ found_end
        %     
        %             if zero_ind_array(array_ind+1) == zero_ind_array(array_ind)+1
        %                 if zero_ind_array(array_ind+2) == zero_ind_array(array_ind)+2
        %                     found_end = true;
        %                 else
        %                     array_ind = array_ind + 3;
        %                 end
        %             else
        %                 array_ind = array_ind + 2;
        %             end
        %         end
                
        %         zero_ind_end = zero_ind_array(array_ind);
        % 
        %         alpha_new = alpha(zero_ind_start:zero_ind_end-1);
        %         a_out_new = a_out(zero_ind_start:zero_ind_end-1);
        
                figure();
                plot(alpha, a_out(2:end))
                hold on
                scatter(alpha, a_out(2:end),20,linspace(10,length(alpha),length(alpha)));
                
%                 figure()
%                 plot(alpha)
%                 plot(alpha(1:zero_ind_array(1)))
%                 plot(diff(alpha))

%                 title(["Bending Angle vs Impact Parameter",name_segment(5:11), "STD ", diff_std]);
%                 xlabel("Alpha (rad)")
%                 ylabel("Impact Parameter (km)")

        end
    
    end

else
    figure();
    plot(alpha, a_out(2:end));


end


title(["Bending Angle vs Impact Parameter",name_segment(5:11)]);
xlabel("Alpha (rad)")
ylabel("Impact Parameter (km)")


if norp == 1
    output_alpha = sprintf('Output_alpha_neg_%7s.txt',name_segment(5:11));
    output_v     = sprintf('Output_v_neg_%7s.txt',name_segment(5:11));
    file_alpha = fopen(output_alpha,'w');
    file_v     = fopen(output_v,'w');
else
    output_alpha = sprintf('Output_alpha_pos_%7s.txt',name_segment(5:11));
    output_v     = sprintf('Output_v_pos_%7s.txt',name_segment(5:11));
    file_alpha = fopen(output_alpha,'w');
    file_v     = fopen(output_v,'w');
end


if ~ plt_bad
    fprintf(file_alpha,'%8.3f  %12.8f\n',[a_out_new; alpha_new']);
else
    fprintf(file_alpha,'%8.3f  %12.8f\n',[a_out(2:end); alpha']);
end
% fprintf(file_v,'%8.3f  %12.8f  %12.8f\n',[a_out; real(v'); imag(v')]);
fclose(file_alpha);
fclose(file_v);

toc
