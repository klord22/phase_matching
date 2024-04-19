%% plots relative error for all occs in flight

%function plot_all(path, norp, plt_go) 
%
% Plots bending angle ERROR vs impact height for results from two of: phase matching,
% geometric optics, and forward model
%  
%% Inputs: 
%   
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
% plot_all('output/',1,false)
%

function plot_all_v0(path, norp, plt_go)

% plt_go is true if comparing pm to go, false if comparing pm to ropp
%cut off at height of aircraft

% need to cut off the ropp op

%occ_name of form g04s_30
tic

% hard-coded
%phase matching
% pm_path_start = 'output_no_waste/Output_alpha_';
% pm_path_end = '.txt';
pm_path_start = [path,'Output_alpha_'];
pm_path_end = '.txt';
info_path = '~/ags/projects/hiaper/2022.305_ar2023/nret/2023.015_iop16/';
occ_start = 1;
occ_end = 4;
sec_in_day = 86400;
start_day = 15; %2023.015
delta_height_pm = .4; %km

go_path = '~/ags/projects/hiaper/2022.305_ar2023/atmPrf/2023.015_iop16/';
go_name_start = 'atmPrf_N49T.2023.';
delta_go = .1;

%should I load the mat files instead? Did Pawel say something?
ropp_path_start = '~/ags/projects/hiaper/2022.305_ar2023/ropp_bangle_2d/2023.015_iop16/';
ropp_name_start = 'fm2d_N49T.2023.';
ropp_file_type = 'nc';
ropp_addition = 5; % number points below pm and go to plot ropp
delta_ropp = .15; %.1 plots but has gap
bin_width_ropp = .2; %km

ropp_levels_path = '~/ags/projects/hiaper/code_ropp/ropp_11.0/ropp-11.0/ropp_io/data/ropp_thin_eum-247.dat';

if plt_go
    round_num = 1; %number points to right of decimal to round impact heights
else
    round_num = 2;
    table = readtable(ropp_levels_path);
    ropp_heights = table.Var2;
end

tolerance = 0.0001; %say doubles are equal if difference is less than this

%plot params
% lw = 2;
% fs = 14;
lw = 4;
fs = 20;
tick_fs = 20;
% tick_fs_small = 16;
% inset_length = 200;
clr = '#8ab5e3';
wide_clr = '#607e9e';



if norp == 1 %neg
    pm_insert = 'neg';
elseif norp == 0 %pos
    pm_insert = 'pos';
else %pos and neg
%     come back to this
end

% same above %%

%hard-coded
file_start = 'Output_alpha';
occ_index_start = 18;
occ_index_end = 24;

len = length(file_start);
directory = dir(path);

count = 0;

figure
hold on


% get info for every occultation
for file = directory'
 
    tic
    
    if length(file.name) > 2

        if file.name(1:len) == file_start

            count = count + 1;
%             disp(count)

            occ_code = file.name(occ_index_start:occ_index_end);

            %below is graph_from_outputs
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

            Rc = yaml_data.Rc;
            aircraft_height = yaml_data.hgt_hi;
            aircraft_height_low = yaml_data.hgt_lo;
            
            sod = rem(sow,sec_in_day); %seconds of the day
            hour = fix(sod/3600); %hours of the day
            % hod_rem = rem(sod,3600); %how many seconds leftover after full hours
            % min = fix(hod_rem/60); %how many minutes leftover after full hours
            % sec = rem(min,60); %how many seconds leftover after full hours

            %phase matching             
            
            pm_data = importdata([pm_path_start,pm_insert,'_',occ_code,pm_path_end]);
            a_pm = pm_data(:,1);
            alpha_pm = pm_data(:,2);
            s(count).pm_alpha = alpha_pm;

            %%% can delete
            s(count).pm_imp_param = a_pm;
            
            %geometric optics
            
            day_padded = sprintf('%03d', day);
            hour_padded = sprintf('%02d', hour);
            
            sat_name = occ_name(1:3);

            s(count).occ = occ_name;
%             s(count).impact_param = a_pm;

            if plt_go

                % retrieve geometric optics information
                go_file_path = strcat(go_path,go_name_start,day_padded,'.',hour_padded,'*',upper(sat_name),'*');
                go_dir = dir(go_file_path);
                go_file_name = go_dir.name;
                
                alpha_go = ncread([go_path,go_file_name],"Bend_ang");
                a_go = ncread([go_path,go_file_name],"Impact_parm");

                % geometric optics output has arrays in opposite order
                alpha_go = flip(alpha_go);
                a_go = flip(a_go);


%                 delta_a = linspace(0,a_go(1)-a_go(end),length(a_go));

%                 s(count).impact_param = delta_a;

                s(count).impact_param = a_go;
                s(count).base_alpha = alpha_go;

%                 up_alpha_go = interp1(a_go,alpha_go,a_pm);

                %downsample phase matching so can compare to geometric
                %optics
                down_alpha_pm = interp1(a_pm,alpha_pm,a_go);

%                 disp(["a_pm:",length(a_pm),"a_go:", length(a_go),"down_alpha_pm", length(down_alpha_pm)])
%                 diff = down_alpha_pm-alpha_go;

            
            else

                % retrieve info from ropp model
                ropp_path = strcat(ropp_path_start,ropp_name_start,day_padded,'.',hour_padded,'*',upper(sat_name),'*',ropp_file_type);
                ropp_dir = dir(ropp_path);
                ropp_file_name = ropp_dir.name;
                
                alpha_ropp = ncread([ropp_path_start,ropp_file_name],"bangle_bg");
                a_ropp = ncread([ropp_path_start,ropp_file_name],"impact_bg")/10^3;

                %start after any nans in the bending angle
                ropp_nan = find(isnan(alpha_ropp) == 1);
                if ~isempty(ropp_nan)
                    % disp("nans")
                    % display(ropp_nan)
                    no_nan_ind = ropp_nan(end)+1;
                else
                    no_nan_ind = 1;
                end
            
                alpha_ropp = alpha_ropp(no_nan_ind:end);
                a_ropp = a_ropp(no_nan_ind:end); 

                %check if there are zeros in bending angle, and start after
                %there are zeros
                %this doesn't seem to be the case for any???
                ropp_zero = find(alpha_ropp == 0);
                if isempty(ropp_zero)
                    ropp_start = 1;
                else
                    ropp_start = ropp_zero(end)+1;
                    disp("zeros")
                end
                disp(occ_name)
                % disp(length(alpha_ropp))
                alpha_ropp = alpha_ropp(ropp_start:end);
                a_ropp = a_ropp(ropp_start:end); 

                disp(alpha_ropp(1))

                %maybe can ignore this
                % disp("is a_ropp = ropp_heights ?")
                % if isempty(find(abs(ropp_heights - a_ropp(1)) < tolerance))
                %     disp("no :(")
                % else
                %     disp("yes")
                % end
               
                
                
                %low is first
                %a_ropp is usually much lower
                %find where to start a_ropp (not too far from start of
                %a_pm)
                if min([a_ropp(1),a_pm(1)]) == a_ropp(1)

                    % disp([a_ropp(1),a_pm(1)])
                    min_a = a_pm(1);

                    %find the closest value in a from the model to the 
                    %minimum in a from PM
                    closest = interp1(a_ropp,a_ropp,min_a,'nearest');
                    low_ind = find(a_ropp == closest);
                
                    if low_ind > ropp_addition
                        a_ropp = a_ropp(low_ind-ropp_addition:end);
                        alpha_ropp = alpha_ropp(low_ind-ropp_addition:end);
                   
               
                    else
                        disp("OTHER")
                        disp(occ_name)
                        disp(count)
                    end
                end


%                 delta_a = linspace(0,a_ropp(1)-a_ropp(end),length(a_ropp));

                s(count).impact_param = a_ropp;
%                 s(count).impact_param = delta_a;
                s(count).base_alpha = alpha_ropp;

                %upsample so can subtract - PM is longest
%                 up_alpha_ropp = interp1(a_ropp,alpha_ropp,a_pm);

                %downsample so can compare
                down_alpha_pm = interp1(a_pm,alpha_pm,a_ropp);
%                 s(count).alpha()

%                 diff = down_alpha_pm-alpha_ropp;


            end
        
            % disp(occ_name)
                
            %put this somewhere else for ropp?
            %remove nans from PM alpha (and also from impact param

            % disp("DOWN ALPHA PM")
            % disp(find(isnan(down_alpha_pm) == 1))

            % disp("original lengths")
            % disp([length(down_alpha_pm),length(s(count).impact_param),length(s(count).base_alpha)])

            down_alpha_pm_nonan = down_alpha_pm(~isnan(down_alpha_pm));
            impact_param_nonan = s(count).impact_param(~isnan(down_alpha_pm));
            base_alpha_nonan = s(count).base_alpha(~isnan(down_alpha_pm));
%             closest_ht = interp1(s(count).impact_height,s(count).impact_height,max(s(count).impact_height)-delta_height_pm,'nearest');

            % disp("updated lengths")
            % disp([length(down_alpha_pm_nonan),length(impact_param_nonan),length(base_alpha_nonan)])
            
%             disp(["a_go no nan",length(impact_param_nonan)])

            if occ_name == 'g14s'
                closest_ht_ind = 8;
            else
                closest_ht = interp1(impact_param_nonan,impact_param_nonan,max(impact_param_nonan)-delta_height_pm,'nearest');
                closest_ht_ind = find(impact_param_nonan == closest_ht);
            end


            final_alpha_pm = down_alpha_pm_nonan(1:closest_ht_ind);
            final_a_other_unrounded = impact_param_nonan(1:closest_ht_ind);
            final_alpha_other = base_alpha_nonan(1:closest_ht_ind);

%             s(count).a_rounded = round(final_a_other_unrounded,1);

%             disp(["final_a_other",length(final_a_other)])

            del_alpha = final_alpha_pm-final_alpha_other; %is this appropriate here for ropp?
            
            %             s(count).alpha_diff = diff;
            s(count).relative_err = del_alpha./final_alpha_other*100; %convert to percent
            s(count).impact_height = final_a_other_unrounded - Rc;
            s(count).impact_height_round = round(s(count).impact_height,round_num);
            
            %find the differences between adjacent impact heights
            test_diff = diff(s(count).impact_height_round);
            
            if ~isempty(find(test_diff > .2))
                disp("BIG DIFF")
                disp(nanmean(test_diff))
            end


%             disp(["final_a_other:",length(final_a_other),"s(count).relative_err",length(s(count).relative_err)])

%             hold on
%             scatter(linspace(1,length(final_a_pm),length(final_a_pm)),final_a_pm,10,linspace(1,length(final_a_pm),length(final_a_pm)))
%             scatter(linspace(1,length(s(count).impact_param),length(s(count).impact_param)),s(count).impact_param,10,linspace(1,length(s(count).impact_param),length(s(count).impact_param)))
%             plot(final_a_other,Color = 'r')
%             plot(s(count).impact_param,Color = clr)
%             hold off

%             legend("PM","Other")
%             disp(["PM:",length(final_alpha_pm),"Other:",length(s(count).impact_param)])

%             plot(s(count).relative_err,s(count).impact_height,Color=clr)


%             plot(linspace(1,length(s(count).impact_height),length(s(count).impact_height)),s(count).impact_height,Color=clr)
            
%             disp(closest_ht)
%             plot(down_alpha_pm,s(count).impact_height,Color=clr)
%             hold on
%             plot(dap_nonan(1:7),imph_nonan(1:7),Color='#f29633')
%            
%             plot(dap_nonan(1:closest_ht),imph_nonan(1:closest_ht),Color='#f29633')
%             plot(dap_nonan(closest_ht_ind:end),imph_nonan(closest_ht_ind:end),Color=clr)
%             hold off

%             if closest_ht_ind ~= 5
%                 disp("NOT FIVE")
%                 disp(closest_ht_ind)
%             end





%             figure
%             scatter(s(count).base_alpha, s(count).impact_param,15,linspace(1,length(s(count).base_alpha),length(s(count).base_alpha)))
%             scatter(down_alpha_pm, s(count).impact_param,10,linspace(1,length(s(count).base_alpha),length(s(count).base_alpha)))

            %%% uncomment this
            plot(s(count).relative_err,s(count).impact_height,Color=clr)

            % disp(["Ropp min",min(s(count).impact_param)])
            % disp(["Rc", Rc])
            % disp(["Ropp max",max(s(count).impact_param)])
            % disp(["high height", aircraft_height])

            % hold off
            % scatter(linspace(1,length(s(count).pm_imp_param),length(s(count).pm_imp_param)),s(count).pm_imp_param,20,'red')
            % hold on
            % scatter(linspace(1,length(s(count).impact_height)-1,length(s(count).impact_height)-1),diff(s(count).impact_height))
            % hold off

            % hold off
            % % scatter(s(count).pm_alpha,s(count).pm_imp_param,20,'blue')
            % hold on
            % scatter(s(count).base_alpha,s(count).impact_param,20,'red')
            % scatter(down_alpha_pm,s(count).impact_param,20,'blue')
            % hold off

            % scatter(count,min(s(count).pm_imp_param),20,'red')
            % scatter(count,max(s(count).pm_imp_param),20,'red')
            % 
            % scatter(count,min(s(count).impact_param),20,'blue')
            % scatter(count,max(s(count).impact_param),20,'blue')

%             plot(diff*10^3,a_pm,LineWidth=lw,Color=clr)
%             plot(diff*10^3,a_pm)
%             plot(s(count).comp_impact*10^3,a_pm)
%             hold on 

        end

    end

    toc

end
          
%find min and max of impact params

all_errs = cat(1,s.relative_err);
all_heights = cat(1,s.impact_height_round);
max_height = max(all_heights);
min_height = min(all_heights);


% go impact height differences: 0.1


ax = gca;
ax.FontSize = tick_fs;

% lgd = legend("Forward Model","Geometric Optics","Phase Matching");


if plt_go

    delta = delta_go;

    xl = xlabel("(PM-GO)/GO * 100%");
    ttl = title("Bending Angle Error vs Impact Height, Geometric Optics");

else

    delta = delta_ropp;
    % min_bin = min_height-bin_width_ropp/2;

    xl = xlabel("(PM-M)/M * 100%");
    ttl = title("Bending Angle Error vs Impact Height, Model");

end
%% 


% figure
hold on

if plt_go

    %make rounded list of impact params
    a_list = round(min_height,1):delta:round(max_height,1);
    len_a_list = length(a_list);
    % scatter(linspace(1,len_a_list,len_a_list),a_list, Color='#4ddb73')
    % disp(a_list(1))


    for c = 1:count
        disp(["c is",c])
    
        %match the first first term in each impact_height_round to a_list
    %     first_ind = find(a_list == s(c).impact_height_round(1));
        first = s(c).impact_height_round(1);
    
        if ~ plt_go
            first_a = a_list(1);
            firsts_diff = first - first_a;
            flr = delta*floor(firsts_diff/delta);
            cel = delta*ceil(firsts_diff/delta);
    
            offset_flr = flr+firsts_diff;
            offset_cel = cel+firsts_diff;
    
            %now see which is below the tolerance
    
            first_ind = find(abs(a_list - offset_flr) < tolerance);
            if isempty(first_ind)
                first_ind = find(abs(a_list - offset_cel) < tolerance);
            end
        
        else
       
            first_ind = find(abs(a_list - first) < tolerance);
    %     first_ind = find(abs(a_list - s(c).impact_height_round) < tolerance);
        end
        
        disp("a_list(first_ind)")
        disp(a_list(first_ind))
        disp("s(c).impact_height_round(1)")
        disp(s(c).impact_height_round(1))
    %     disp("differences")
    %     disp(a_list - s(c).impact_height_round)
        %pad with nans before and after
    %     disp(NaN(first_ind-1,1))
    %     disp(s(count).impact_height_round)
    %     disp(NaN(len_a_list-length(s(c).impact_height_round),1))
    
    %     disp("first_ind-1")
    %     disp(first_ind-1)
    %     disp("len_a_list-length(s(c).impact_height_round)-(first_ind-1)")
    %     disp(len_a_list-length(s(c).impact_height_round)-(first_ind-1))
    %     disp("len_a_list")
    %     disp(len_a_list)
    %     disp("length(s(c).impact_height_round)")
    %     disp(length(s(c).impact_height_round))
    
    
    
    
    %%% I think error is here
    %%% DONT NEED TO ROUND HERE
    
    
    
        s(c).padded_height = cat(1,NaN(low_ind-1,1),s(c).impact_height_round,NaN(len_a_list-length(s(c).impact_height_round)-(low_ind-1),1));
        s(c).padded_error = cat(1,NaN(low_ind-1,1),s(c).relative_err,NaN(len_a_list-length(s(c).relative_err)-(low_ind-1),1));
    
        % s(c).padded_height = cat(1,NaN(first_ind-1,1),s(c).impact_height_round,NaN(len_a_list-length(s(c).impact_height_round)-(first_ind-1),1));
        % s(c).padded_error = cat(1,NaN(first_ind-1,1),s(c).relative_err,NaN(len_a_list-length(s(c).relative_err)-(first_ind-1),1));
    
        %     disp([first_ind-1,length(s(c).impact_height_round),len_a_list-length(s(c).impact_height_round),"=",first_ind-1+length(s(c).impact_height_round)+len_a_list-length(s(c).impact_height_round)])
    %     disp([len_a_list,c,length(s(c).padded_height)])
    
    %     scatter(linspace(1,length(s(c).impact_height_round),length(s(c).impact_height_round)),s(c).impact_height_round)
    
        if length(s(c).padded_height) ~= len_a_list
            disp(length(s(c).padded_height))
            disp("PANIC")
            % disp(["first_ind",first_ind])
            % disp(["s(c).impact_height_round(1)",s(c).impact_height_round(1)])
            % disp(["length imp height",length(s(c).impact_height_round)])
            % disp(["s(c).impact_height_round(end)",s(c).impact_height_round(end)])
            % disp(["a_list end", a_list(end)])
        end
    
        if isempty(s(c).padded_height)
            disp(c)
            disp("EMPTY HEIGHT")
        end
        if isempty(s(c).padded_error)
            disp(c)
            disp("EMPTY ERROR")
        end
        disp(["size height",size(s(c).padded_height)])
        disp(["size error",size(s(c).padded_error)])

    end

    error_matrix = cat(2,s.padded_error);
    % error_mean = nanmean(error_matrix,2);
    %weighting is standard(0), dimension is 2, ignore nans
    [std_error, mean_error] = std(error_matrix,0,2,"omitnan");

   



    %make sure mean_error and a_list are in the same direction!!!
    plot(mean_error,a_list, LineWidth = lw, Color=wide_clr);
    plot(mean_error+std_error,a_list, LineWidth = lw-2, Color=wide_clr);
    plot(mean_error-std_error,a_list, LineWidth = lw-2, Color=wide_clr);


else 
    %make a matrix of all the bending angles and impact parameters padded
    %with NaNs AT THE END - just to make them all the same length 

    % take the mean and std of each bin and put into array
    min_bin = min_height-bin_width_ropp/3;
    num_bins = fix((max_height - min_bin)/bin_width_ropp);

    % would be nice to be able to pass in s.impact_height and
    % s.relative_err
    % need to be able to take means and stds of everything
    [mean_error, std_error, edges] = binned_std(all_heights,all_errs,num_bins);
    % for c = 1:count
    %     [m, s] = binned_std(s(c).impact_height,s(c).relative_err,num_bins);
    % 
    % end
    % a_bins = discretize(s.impact_height,num_bins);

    midpoints = edges(1:end-1)+diff(edges)/2;

    plot(mean_error,midpoints, LineWidth = lw, Color=wide_clr);
    plot(mean_error+std_error,midpoints, LineWidth = lw-2, Color=wide_clr);
    plot(mean_error-std_error,midpoints, LineWidth = lw-2, Color=wide_clr);

    % scatter(mean_error,midpoints, 15,'filled',Color=wide_clr);
    % scatter(mean_error+std_error,midpoints, 10,'filled', Color=wide_clr);
    % scatter(mean_error-std_error,midpoints, 10,'filled', Color=wide_clr);
    
    % for i = 1:length(edges)
    %     plot([-20,20],[edges(i),edges(i)],LineWidth = 1, Color = 'cyan')
    % end

% then use results to bin alpha
% treat each point like its own thing, devoid of occultation
end




% save("model_plot_all_vars")


yl = ylabel("Impact Height (km)");
% fontsize(lgd,fs,'points')
fontsize(xl,fs,'points')
fontsize(yl,fs,'points')
fontsize(ttl,fs,'points')






% start = fix(length(a_ropp)/2);
% start = fix(length(a_pm)/2);
% stop = start+inset_length;
% 
% axes('Position',[.2 .2 .3 .3])
% box on
% plot(up_alpha_ropp(start:stop)*10^3,a_pm(start:stop),LineWidth=lw+2,Color=ropp_color)
% hold on
% plot(up_alpha_go(start:stop)*10^3,a_pm(start:stop),LineWidth=lw+2,Color=go_color,LineStyle="--")
% plot(alpha_pm(start:stop)*10^3,a_pm(start:stop),LineWidth=lw+2,Color=pm_color)
% yticks(linspace(round(a_pm(start),1),round(a_pm(stop),1),3));
% 
% ax = gca;
% ax.FontSize = tick_fs_small;
% 
% plot(up_alpha_ropp(start:stop)*10^3,a_pm(start:stop),LineWidth=lw+2)
% hold on
% plot(up_alpha_go(start:stop)*10^3,a_pm(start:stop),LineWidth=lw+2,LineStyle="--")
% plot(alpha_pm(start:stop)*10^3,a_pm(start:stop),LineWidth=lw+2)


toc
