%% plot forward model, GO, PM
%% PM not cut off at top

function graph_from_outputs(occ_code, norp)

%occ_name of form g04s_30
%a_ropp is singles, not doubles - will this cause problems?
tic
%% hard-coded
%phase matching
pm_path_start = 'output_no_waste/Output_alpha_';
pm_path_end = '.txt';
info_path = '~/ags/projects/hiaper/2022.305_ar2023/nret/2023.015_iop16/';
occ_start = 1;
occ_end = 4;
sec_in_day = 86400;
start_day = 15; %2023.015

go_path = '~/ags/projects/hiaper/2022.305_ar2023/atmPrf/2023.015_iop16/';
go_name_start = 'atmPrf_N49T.2023.';

ropp_path_start = '~/ags/projects/hiaper/2022.305_ar2023/ropp_bangle_2d/2023.015_iop16/';
ropp_name_start = 'fm2d_N49T.2023.';
ropp_file_type = 'nc';
ropp_addition = 5; % number points below pm and go to plot ropp

%plot params
% lw = 2;
% fs = 14;
lw = 4;
fs = 20;
tick_fs = 20;
tick_fs_small = 16;
inset_length = 200;
ropp_color = 'b';
go_color = '#f29633'; %'#f29633'
pm_color = '#3fbf24'; %'#3fbf24'



if norp == 1 %neg
    pm_insert = 'neg';
elseif norp == 0 %pos
    pm_insert = 'pos';
else %pos and neg
%     come back to this
end

occ_name = occ_code(occ_start:occ_end);
info_file_path = [info_path,occ_name,'*'];
info_dir = dir(info_file_path);
folder_name = info_dir.name;

%calc date from sow and sec_in_week
%gps week starts on sunday
%2023.01.15 is a sunday

yaml_data = yaml.loadFile([info_path,folder_name,'/info.yaml']);
sow = yaml_data.sow_low;
dow = fix(sow/sec_in_day);
day = start_day+dow; %should be 15 or 16

sod = rem(sow,sec_in_day); %seconds of the day
hour = fix(sod/3600); %hours of the day
% hod_rem = rem(sod,3600); %how many seconds leftover after full hours
% min = fix(hod_rem/60); %how many minutes leftover after full hours
% sec = rem(min,60); %how many seconds leftover after full hours

%geometric optics

day_padded = sprintf('%03d', day);
hour_padded = sprintf('%02d', hour);

sat_name = occ_name(1:3);

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
a_pm = pm_data(:,1);
alpha_pm = pm_data(:,2);

%plotting

ropp_zero = find(alpha_ropp == 0);

if isempty(ropp_zero)
    ropp_start = 1;
else
    ropp_start = ropp_zero(end)+1;
end

alpha_ropp = alpha_ropp(ropp_start:end); 
a_ropp = a_ropp(ropp_start:end);

%% mod here

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


%% end mod here

figure
plot(alpha_ropp*10^3,a_ropp,LineWidth=lw,Color=ropp_color)
hold on
plot(alpha_go*10^3,a_go,LineWidth=lw,Color=go_color) %LineStyle="--"
plot(alpha_pm*10^3,a_pm,LineWidth=lw,Color=pm_color)

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

up_alpha_ropp = interp1(a_ropp,alpha_ropp,a_pm);
up_alpha_go = interp1(a_go,alpha_go,a_pm);

% start = fix(length(a_ropp)/2);
start = fix(length(a_pm)/2);
stop = start+inset_length;

axes('Position',[.2 .2 .3 .3])
box on
plot(up_alpha_ropp(start:stop)*10^3,a_pm(start:stop),LineWidth=lw+2,Color=ropp_color)
hold on
plot(up_alpha_go(start:stop)*10^3,a_pm(start:stop),LineWidth=lw+2,Color=go_color) %LineStyle="--"
plot(alpha_pm(start:stop)*10^3,a_pm(start:stop),LineWidth=lw+2,Color=pm_color)
yticks(linspace(round(a_pm(start),1),round(a_pm(stop),1),3));

ax = gca;
ax.FontSize = tick_fs_small;
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