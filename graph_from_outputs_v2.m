%function graph_from_outputs_v2(occ_code, norp) 
%
% Plots bending angle vs impact height for results from phase matching,
% geometric optics, and forward model on one graph
%  
%% Inputs: 
%   occ_code: a string, the name of the folder containing the
%   phase-matching files. Takes the form g04s_30
%   norp: means negative or positive elevation angle. 0 for positive, 1 for
%   negative
%
%% Instructions: 
%   User must change the hard-coded section to the appropriate directories
%   for phase matching, geometric optics, and forward model data
%
%% Dependencies
%  No direct dependencies, but must already have generated BA profiles
%  using phase matching, geometric optics, and forward modelling. Phase
%  matching can be done with dop2alpha_pm_v4.m for one occultation, or
%  loop_dop2alpha.m for a directory of occultations
%
%% Example:
%
% graph_from_outputs_v2('g04s_30',1)
%


function graph_from_outputs_v2(occ_code, norp)

%occ_name of form g04s_30
%a_ropp is singles, not doubles - will this cause problems?

tic
%% Hard-coded

%phase matching
pm_path_start = 'output_no_waste/Output_alpha_';
pm_path_end = '.txt';
info_path = '~/ags/projects/hiaper/2022.305_ar2023/nret/2023.015_iop16/';
% info_path = '~/ags/projects/hiaper/2021.019_ar2021/nret_pppar/2021.023_iop04/';

% will call the no-filt as pm
% pm_path_start = '~/ags/projects/hiaper/2023.304_ar2024/aro/6_nret_nofilt/2024.015_iop_17_n_49_t_no_filt/';
% pm_path_end = '.txt';
% info_path = '~/ags/projects/hiaper/2023.304_ar2024/aro/6_nret_nofilt/2024.015_iop_17_n_49_t_no_filt/';

occ_start = 1;
occ_end = 4;
sec_in_day = 86400;
start_day = 15; %2023.015
delta_height_pm = .4; %km

go_path = '~/ags/projects/hiaper/2022.305_ar2023/atmPrf/2023.015_iop16/';
% % go_path = '~/ags/projects/hiaper/2021.019_ar2021/atmPrf_pppar/2021.023_iop04/';
go_name_start = 'atmPrf_N49T.2023.';
% % go_name_start = 'atmPrf_N49T.2021.';

% call original filter go
% go_path = '~/ags/projects/hiaper/2023.304_ar2024/aro/6_nret_nofilt/2024.015_iop17_n49t';
% go_name_start = 'atmPrf_N49T.2024.';

ropp_path_start = '~/ags/projects/hiaper/2022.305_ar2023/ropp_bangle_2d/2023.015_iop16/';
% ropp_path_start = '~/ags/projects/hiaper/2021.019_ar2021/ropp_bangle_2d/2021.023_iop04/';
ropp_name_start = 'fm2d_N49T.2023.';
% ropp_name_start = 'fm2d_N49T.2021.';
ropp_file_type = 'nc';
ropp_addition = 5; % number points below pm and go to plot ropp

%plot params
% lw = 4;
% fs = 20;
% tick_fs = 20;
% tick_fs_small = 16;
% inset_length = 200;
% ropp_color = 'b';
% go_color = '#f29633'; %'#f29633'
% pm_color = '#3fbf24'; %'#3fbf24'

lw = 2;
fs = 20;
tick_fs = 20;
tick_fs_small = 16;
inset_length = 200;
ropp_color = 'b';
go_color = '#f29633'; %'#f29633'
pm_color = '#3fbf24'; %'#3fbf24'

%% Get data from files

if norp == 1 %neg
    pm_insert = 'neg';
elseif norp == 0 %pos
    pm_insert = 'pos';
else %pos and neg
%     come back to this
end

% get the directory for the info file

occ_name = occ_code(occ_start:occ_end);
info_file_path = [info_path,occ_name,'*'];
info_dir = dir(info_file_path);
folder_name = info_dir.name;

%calc date from sow and sec_in_week
%gps week starts on sunday
%2023.01.15 is a sunday


% info file

yaml_data = yaml.loadFile([info_path,folder_name,'/info.yaml']);
sow = yaml_data.sow_low;
Rc = yaml_data.Rc;
% sow = 609706;

% get hour and day to differentiate GO files

dow = fix(sow/sec_in_day);
day = start_day+dow; %should be 15 or 16

sod = rem(sow,sec_in_day); %seconds of the day
hour = fix(sod/3600); %hours of the day
% hod_rem = rem(sod,3600); %how many seconds leftover after full hours
% min = fix(hod_rem/60); %how many minutes leftover after full hours
% sec = rem(min,60); %how many seconds leftover after full hours

day_padded = sprintf('%03d', day);
% day_padded = '024';
hour_padded = sprintf('%02d', hour);

sat_name = occ_name(1:3);

%geometric optics

go_file_path = strcat(go_path,go_name_start,day_padded,'.',hour_padded,'*',upper(sat_name),'*');
go_dir = dir(go_file_path);
go_file_name = go_dir.name;

alpha_go = ncread([go_path,go_file_name],"Bend_ang");
a_go = ncread([go_path,go_file_name],"Impact_parm");


%ropp

ropp_path = strcat(ropp_path_start,ropp_name_start,day_padded,'.',hour_padded,'*',upper(sat_name),'*',ropp_file_type);
ropp_dir = dir(ropp_path);
ropp_file_name = ropp_dir.name;

% alpha_ropp = ncread([ropp_path_start,ropp_file_name],"bangle");
% a_ropp = ncread([ropp_path_start,ropp_file_name],"impact")/10^3;
alpha_ropp = ncread([ropp_path_start,ropp_file_name],"bangle_bg"); %bg is simulations
a_ropp = ncread([ropp_path_start,ropp_file_name],"impact_bg")/10^3;


%phase matching             

pm_data = importdata([pm_path_start,pm_insert,'_',occ_code,pm_path_end]);
% pm_data = importdata('output_no_waste/Output_alpha_neg_r02r_r0.txt');
a_pm = pm_data(:,1);
alpha_pm = pm_data(:,2);

%% Plotting modifications

% cut off phase matching profile at top
if occ_name == 'g14s'
    closest_ht_ind = 8;
else
    closest_ht = interp1(a_pm,a_pm,max(a_pm)-delta_height_pm,'nearest');
    closest_ht_ind = find(a_pm == closest_ht);
end

a_pm = a_pm(1:closest_ht_ind);
alpha_pm = alpha_pm(1:closest_ht_ind);

%cut off ropp profile where alpha <= 0
ropp_zero = find(alpha_ropp == 0);

if isempty(ropp_zero)
    ropp_start = 1;
else
    ropp_start = ropp_zero(end)+1;
end

alpha_ropp = alpha_ropp(ropp_start:end); 
a_ropp = a_ropp(ropp_start:end);

% modified from eric's code? here

alpha_go = flip(alpha_go);
a_go = flip(a_go);

%low is first
%a_ropp is usually much lower
if min([a_ropp(1),a_go(1),a_pm(1)]) == a_ropp(1)
    min_a = min(a_go(1),a_pm(1));
    closest = interp1(a_ropp,a_ropp,min_a,'nearest');
    low_ind = find(a_ropp == closest);

    if low_ind > ropp_addition
        a_ropp = a_ropp(low_ind-ropp_addition:end);
        alpha_ropp = alpha_ropp(low_ind-ropp_addition:end);
    end

    
end

% end mod here


% Rc = 6361.7111;

figure
% plot(alpha_ropp*10^3,a_ropp-Rc,LineWidth=lw,Color=ropp_color)
% hold on
% plot(alpha_go*10^3,a_go-Rc,LineWidth=lw,Color=go_color) %LineStyle="--"
% plot(alpha_pm*10^3,a_pm-Rc,LineWidth=lw,Color=pm_color)

plot(alpha_ropp*10^3,a_ropp-Rc,LineWidth=lw)
hold on
plot(alpha_go*10^3,a_go-Rc,LineWidth=lw) %LineStyle="--"
plot(alpha_pm*10^3,a_pm-Rc,LineWidth=lw, color = 'black')

ax = gca;
ax.FontSize = tick_fs;

% plot(alpha_ropp*10^3,a_ropp,LineWidth=lw)
% hold on
% plot(alpha_go*10^3,a_go,LineWidth=lw,LineStyle="--")
% plot(alpha_pm*10^3,a_pm,LineWidth=lw)

lgd = legend("Forward Model","Geometric Optics","Phase Matching");
xl = xlabel("Bending Angle (millirad)");
yl = ylabel("Impact Parameter (km)");
ttl = title(["Bending Angle vs Impact Parameter,",occ_name]);
fontsize(lgd,fs,'points')
fontsize(xl,fs,'points')
fontsize(yl,fs,'points')
fontsize(ttl,fs,'points')


%upsample so can subtract - PM is longest

% up_alpha_ropp = interp1(a_ropp,alpha_ropp,a_pm);
% up_alpha_go = interp1(a_go,alpha_go,a_pm);

% start = fix(length(a_ropp)/2);

%% Plots an zoomed in inset

% start = fix(length(a_pm)/2);
% stop = start+inset_length;
% 
% axes('Position',[.2 .2 .3 .3])
% box on
% plot(up_alpha_ropp(start:stop)*10^3,a_pm(start:stop),LineWidth=lw+2,Color=ropp_color)
% hold on
% plot(up_alpha_go(start:stop)*10^3,a_pm(start:stop),LineWidth=lw+2,Color=go_color) %LineStyle="--"
% plot(alpha_pm(start:stop)*10^3,a_pm(start:stop),LineWidth=lw+2,Color=pm_color)
% yticks(linspace(round(a_pm(start),1),round(a_pm(stop),1),3));
% 
% ax = gca;
% ax.FontSize = tick_fs_small;


%% Plots some light differences

% 
% plot(up_alpha_ropp(start:stop)*10^3,a_pm(start:stop),LineWidth=lw+2)
% hold on
% plot(up_alpha_go(start:stop)*10^3,a_pm(start:stop),LineWidth=lw+2,LineStyle="--")
% plot(alpha_pm(start:stop)*10^3,a_pm(start:stop),LineWidth=lw+2)

% figure
% plot(up_alpha_ropp-up_alpha_go,a_pm)
% hold on
% plot(up_alpha_ropp-alpha_pm, a_pm)
% legend("ROPP - GO", "ROPP - PM")

toc
