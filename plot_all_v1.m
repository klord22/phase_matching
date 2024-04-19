%% plots relative error for all occs in flight

%function plot_all(path, norp, plt_which) 
%
% Plots bending angle ERROR vs impact height for results from two of: phase matching,
% geometric optics, and forward model
%  
%% Inputs: 
%   path: string, path to directory containing phase-matching output files.
%   norp: means negative or positive elevation angle. 0 for positive, 1 for
%   negative
%   plt_which: 0 for pm vs go, 1 for pm vs ropp, 2 for go vs ropp
%
%% Algorithm:
%   Phase-matching and geometric optics have output equally spaced in impact
%   parameter. ERA5 model output does not.
%       * read in output files
%       * interpolate phase-matching to other array
%       For geometric optics:
%       * round heights and create grid of impact heights
%       * match each array of impact heights to subsection of grid
%       For model:
%       * bin points by impact height
%       For all:
%       * find relative error (%) for two methods being compared, and
%       calculate mean and standard deviation at every height bin
%
%% Instructions: 
%   User must change the hard-coded section to the appropriate directories
%   for phase matching, geometric optics, and forward model data
%
%% Dependencies
%  binned_std.m
%  Must already have generated BA profiles using phase matching, geometric 
%  optics, and forward modelling. Phase matching can be done with 
%  dop2alpha_pm_v4.m for one occultation, or loop_dop2alpha.m for a 
%  directory of occultations
%
%% Example:
%
% plot_all_v1('output/',1,1)
%

function plot_all_v1(path, norp, plt_which)

% plt_go is true if comparing pm to go, false if comparing pm to ropp
%cut off at height of aircraft

% need to cut off the ropp op

%occ_name of form g04s_30
tic

%%% hard-coded

%phase matching

pm_path_start = [path,'Output_alpha_'];
pm_path_end = '.txt';
info_path = '~/ags/projects/hiaper/2022.305_ar2023/nret/2023.015_iop16/';
occ_start = 1; %indices for occultation name
occ_end = 4;
sec_in_day = 86400;
start_day = 15; %2023.015
delta_height_pm = .4; %km

%geometric optics

go_path = '~/ags/projects/hiaper/2022.305_ar2023/atmPrf/2023.015_iop16/';
go_name_start = 'atmPrf_N49T.2023.';
delta_go = .1; %difference between adjacent heights used for grid
round_num = 1; %number points to right of decimal to round impact heights

%model
%should I load the mat files instead? Did Pawel say something?
ropp_path_start = '~/ags/projects/hiaper/2022.305_ar2023/ropp_bangle_2d/2023.015_iop16/';
% ropp_path_start = '~/ags/projects/hiaper/arivers/kate/code_mat_ropp_aro/2023.015_iop16/';
ropp_name_start = 'fm2d_N49T.2023.';
ropp_file_type = 'nc';
ropp_addition = 5; % number points below pm and go to plot ropp
bin_width_ropp = .2; % in km, width of height bins for taking means and stds

% ropp_levels_path = '~/ags/projects/hiaper/code_ropp/ropp_11.0/ropp-11.0/ropp_io/data/ropp_thin_eum-247.dat';

% if plt_go
%     round_num = 1; %number points to right of decimal to round impact heights
% else
%     round_num = 2;
%     table = readtable(ropp_levels_path);
%     ropp_heights = table.Var2;
% end

tolerance = 0.0001; %say doubles are equal if difference is less than this

%plot params

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

file_start = 'Output_alpha';
occ_index_start = 18;
occ_index_end = 24;

%%% end of hard-coded section

len = length(file_start);
directory = dir(path);

% read_pm = false;
% read_go = false;
% read_ropp = false;
% 
% if plt_which == 0 %pm vs go
%     read_pm = true;
%     read_gp = false;


count = 0;

figure
hold on

num_files = 0;
% count to preallocate
disp("time just for counting")
% tic
for file = directory'

    if length(file.name) > 2

        if file.name(1:len) == file_start

            num_files = num_files + 1;
        end
    end
end
% toc

s=struct('occ',cell(num_files,1), 'pm_alpha',cell(num_files,1), 'impact_param',cell(num_files,1),'base_alpha',cell(num_files,1),'relative_err',cell(num_files,1),'impact_height',cell(num_files,1),'impact_height_round',cell(num_files,1),'padded_height',cell(num_files,1),'padded_error',cell(num_files,1));    
    % occ
    % pm_alpha
    % impact_param
    % base_alpha
    % relative_err
    % impact_height
    % impact_height_round
    % padded_height
    % padded_error


% get info for every occultation
for file = directory'
 
    % tic
    
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
            
            day_padded = sprintf('%03d', day);
            hour_padded = sprintf('%02d', hour);
            
            sat_name = occ_name(1:3);

            s(count).occ = occ_name;
%             s(count).impact_param = a_pm;

            disp(occ_name)

            if plt_which == 0 || plt_which == 1

                %phase matching

                pm_data = importdata([pm_path_start,pm_insert,'_',occ_code,pm_path_end]);
                a_pm = pm_data(:,1);
                alpha_pm = pm_data(:,2);
                s(count).pm_alpha = alpha_pm;

            end
            % %%% can delete
            % s(count).pm_imp_param = a_pm;
            
            %geometric optics
            


            if plt_which == 0 || plt_which == 2

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
                
                %save impact parameter and bending angle to struct
                s(count).impact_param = a_go;
                s(count).base_alpha = alpha_go;

%                 up_alpha_go = interp1(a_go,alpha_go,a_pm);

                if plt_which == 0
                    %downsample phase matching so can compare to geometric
                    %optics
                    down_alpha = interp1(a_pm,alpha_pm,a_go);

                end

%                 disp(["a_pm:",length(a_pm),"a_go:", length(a_go),"down_alpha_pm", length(down_alpha_pm)])
%                 diff = down_alpha_pm-alpha_go;

            end

            if plt_which == 1 || plt_which == 2

                %forward model

                % retrieve info from ropp model results
                ropp_path = strcat(ropp_path_start,ropp_name_start,day_padded,'.',hour_padded,'*',upper(sat_name),'*',ropp_file_type);
                ropp_dir = dir(ropp_path);

                %%% PROBLEM 
                if isempty(ropp_dir)
                    continue

                else
                    disp("NOT EMPTY")
                end

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
                %would this be needed for go data?
                ropp_zero = find(alpha_ropp == 0);
                if isempty(ropp_zero)
                    ropp_start = 1;
                else
                    ropp_start = ropp_zero(end)+1;
                    disp("zeros")
                end
                
                alpha_ropp = alpha_ropp(ropp_start:end);
                a_ropp = a_ropp(ropp_start:end); 

                if plt_which == 1
                    min_a = a_pm(1);
                    a_2b_down = a_pm;
                    alpha_2b_down = alpha_pm;
                elseif plt_which == 2
                    min_a = a_go(1);
                    a_2b_down = a_go;
                    alpha_2b_down = alpha_go;
                end

                %find where to start a_ropp (not too far from start of obs)
                %a_ropp is usually much lower 
                if min([a_ropp(1),min_a]) == a_ropp(1)
    
                    %find the closest value in a from the model to the 
                    %minimum in a from obs
                    closest = interp1(a_ropp,a_ropp,min_a,'nearest');
                    low_ind = find(a_ropp == closest);
                
                    % shorten the ropp data if the closest value to the
                    % start of the obs data is further than 5 points from
                    % the start of the ropp data
                    if low_ind > ropp_addition
                        a_ropp = a_ropp(low_ind-ropp_addition:end);
                        alpha_ropp = alpha_ropp(low_ind-ropp_addition:end);                 
               
                    else
                        disp("OTHER")
                        disp(occ_name)
                    end

                    disp(["low_ind top", low_ind])
                end
    
                %downsample so can compare
                down_alpha = interp1(a_2b_down,alpha_2b_down,a_ropp);

                %%% to be changed
    %             if plt_which == 1
    %                 %find where to start a_ropp (not too far from start of
    %                 %a_pm)
    %                 %a_ropp is usually much lower                
    %                 if min([a_ropp(1),a_pm(1)]) == a_ropp(1)
    % 
    %                     min_a = a_pm(1);
    % 
    %                     %find the closest value in a from the model to the 
    %                     %minimum in a from PM
    %                     closest = interp1(a_ropp,a_ropp,min_a,'nearest');
    %                     low_ind = find(a_ropp == closest);
    % 
    %                     % shorten the ropp data if the closest value to the
    %                     % start of the pm data is further than 5 points from
    %                     % the start of the ropp data
    %                     if low_ind > ropp_addition
    %                         a_ropp = a_ropp(low_ind-ropp_addition:end);
    %                         alpha_ropp = alpha_ropp(low_ind-ropp_addition:end);                 
    % 
    %                     % else
    %                     %     disp("OTHER")
    %                     %     disp(occ_name)
    %                     end
    % 
    %                     disp(["low_ind top", low_ind])
    %                 end
    % 
    % 
    % %                 delta_a = linspace(0,a_ropp(1)-a_ropp(end),length(a_ropp));
    % 
    %                 %upsample so can subtract - PM is longest
    % %                 up_alpha_ropp = interp1(a_ropp,alpha_ropp,a_pm);
    % 
    %                 %downsample so can compare
    %                 down_alpha_pm = interp1(a_pm,alpha_pm,a_ropp);
    % 
    %             % where this should be
    %             elseif plt_which == 2
    %                 % figure 
    %                 % hold on 
    %                 % scatter(alpha_ropp, a_ropp)
    %                 % scatter(alpha_go, a_go)
    %                 % legend("ropp","go")
    % 
    %                 %find where to start a_ropp (not too far from start of
    %                 %a_go)
    %                 %a_ropp is usually much lower                
    %                 if min([a_ropp(1),a_go(1)]) == a_ropp(1)
    % 
    %                     min_a = a_pm(1);
    % 
    %                     %find the closest value in a from the model to the 
    %                     %minimum in a from PM
    %                     closest = interp1(a_ropp,a_ropp,min_a,'nearest');
    %                     low_ind = find(a_ropp == closest);
    % 
    %                     % shorten the ropp data if the closest value to the
    %                     % start of the pm data is further than 5 points from
    %                     % the start of the ropp data
    %                     if low_ind > ropp_addition
    %                         a_ropp = a_ropp(low_ind-ropp_addition:end);
    %                         alpha_ropp = alpha_ropp(low_ind-ropp_addition:end);                 
    % 
    %                     else
    %                         disp("OTHER")
    %                         disp(occ_name)
    %                     end
    % 
    %                     disp(["low_ind top", low_ind])
    %                 end
    % 
    %                 %downsample so can compare
    %                 down_alpha_pm = interp1(a_pm,alpha_pm,a_ropp);
                % 
                % end
                %%% to be changed

                % end
                
                %bring out, unindent
                s(count).impact_param = a_ropp;
%                 s(count).impact_param = delta_a;
                s(count).base_alpha = alpha_ropp;
            end    
                
                %%% start change here
                % 
                % 
                % if plt_which == 1
                %     min_a = a_pm(1);
                %     a_down = a_pm;
                %     alpha_down = alpha_pm;
                % elseif plt_which == 2
                %     min_a = a_go(1);
                %     a_down = a_go;
                %     alpha_down = alpha_go;
                % end
                
                %put this somewhere else for ropp?
                %remove nans from PM alpha (and also the same points from 
                %impact param and bending angle

            down_alpha_nonan = down_alpha(~isnan(down_alpha)); 
            impact_param_nonan = s(count).impact_param(~isnan(down_alpha));
            base_alpha_nonan = s(count).base_alpha(~isnan(down_alpha));

            closest_ht = interp1(impact_param_nonan,impact_param_nonan,max(impact_param_nonan)-delta_height_pm,'nearest');
            closest_ht_ind = find(impact_param_nonan == closest_ht);
    
            final_down_alpha = down_alpha_nonan(1:closest_ht_ind);
            final_base_a_unrounded = impact_param_nonan(1:closest_ht_ind);
            final_base_alpha = base_alpha_nonan(1:closest_ht_ind);
    
    %             s(count).a_rounded = round(final_a_other_unrounded,1);
    
    %             disp(["final_a_other",length(final_a_other)])
    
            del_alpha = final_down_alpha-final_base_alpha; %is this appropriate here for ropp?


    %             down_alpha_pm_nonan = down_alpha_pm(~isnan(down_alpha_pm));
    %             impact_param_nonan = s(count).impact_param(~isnan(down_alpha_pm));
    %             base_alpha_nonan = s(count).base_alpha(~isnan(down_alpha_pm));
    % 
    %             closest_ht = interp1(impact_param_nonan,impact_param_nonan,max(impact_param_nonan)-delta_height_pm,'nearest');
    %             closest_ht_ind = find(impact_param_nonan == closest_ht);
    % 
    %             final_alpha_pm = down_alpha_pm_nonan(1:closest_ht_ind);
    %             final_a_other_unrounded = impact_param_nonan(1:closest_ht_ind);
    %             final_alpha_other = base_alpha_nonan(1:closest_ht_ind);
    % 
    % %             s(count).a_rounded = round(final_a_other_unrounded,1);
    % 
    % %             disp(["final_a_other",length(final_a_other)])
    % 
    %             del_alpha = final_alpha_pm-final_alpha_other; %is this appropriate here for ropp?
    %             %%% end change here



                %             s(count).alpha_diff = diff;
            s(count).relative_err = del_alpha./final_base_alpha*100; %convert to percent
            s(count).impact_height = final_base_a_unrounded - Rc;
            s(count).impact_height_round = round(s(count).impact_height,round_num);
          
                %wrong place, but good for debugging
                % elseif plt_which == 2
                %     figure 
                %     hold on 
                %     scatter(alpha_ropp, a_ropp)
                %     scatter(alpha_go, a_go)
                %     legend("ropp","go")
                %     title(occ_name)
                % 
                %     %go is shorter (does not extend to bottom)
                %     %go has higher freq in impact parameter space - should
                %     %downsample go to ropp

                % end
            % end % delete this? this is where it used to end
            plot(s(count).relative_err,s(count).impact_height,Color=clr)
            
            % figure
            % hold on
            % scatter(alpha_pm, a_pm)
            % scatter(alpha_ropp,a_ropp, 'filled')
            % % scatter(alpha_go, a_go, 'filled')
            % title(["PM and Model,", occ_name])
            % legend("PM", "Model")
            % hold off

        end

    end

    % toc

end
          
%find min and max of impact params

all_errs = cat(1,s.relative_err);
all_heights = cat(1,s.impact_height_round);
max_height = max(all_heights);
min_height = min(all_heights);


% go impact height differences: 0.1


ax = gca;
ax.FontSize = tick_fs;


% figure
hold on

if plt_which == 0 %pm vs go

    delta = delta_go;

    xl = xlabel("(PM-GO)/GO * 100%");
    ttl = title("Bending Angle Error vs Impact Height, PM vs GO");

    %make rounded list of impact params
    a_list = round(min_height,1):delta:round(max_height,1);
    len_a_list = length(a_list);

    %we have made a list of impact params. Because the impact params have
    %all been rounded to the same place, a_list is a "ladder" of heights,
    %and every impact_height_round is a subset of contiguous "rungs". We
    %want to slide each impact_height_round so it matches a subset of
    %a_list, and then pad the beginning and end with NaNs so it is the same
    %length.
    for c = 1:count
        disp(["c is",c])
    
        %match the first first term in each impact_height_round to a_list
        first = s(c).impact_height_round(1);
        first_ind = find(abs(a_list - first) < tolerance);
        
        %add NaNs to beginning and end of impact height and error arrays
        s(c).padded_height = cat(1,NaN(first_ind-1,1),s(c).impact_height_round,NaN(len_a_list-length(s(c).impact_height_round)-(first_ind-1),1));
        s(c).padded_error = cat(1,NaN(first_ind-1,1),s(c).relative_err,NaN(len_a_list-length(s(c).relative_err)-(first_ind-1),1));
     
        %check for issues with padded arrays
        if length(s(c).padded_height) ~= len_a_list
            disp(length(s(c).padded_height))
            disp("PANIC")
        end
    
        if isempty(s(c).padded_height)
            disp(c)
            disp("EMPTY HEIGHT")
        end

        if isempty(s(c).padded_error)
            disp(c)
            disp("EMPTY ERROR")
        end

        % disp(["size height",size(s(c).padded_height)])
        % disp(["size error",size(s(c).padded_error)])

    end

    error_matrix = cat(2,s.padded_error);

    %weighting is standard(0), dimension is 2, ignore nans
    [std_error, mean_error] = std(error_matrix,0,2,"omitnan");

    plot(mean_error,a_list, LineWidth = lw, Color=wide_clr);
    plot(mean_error+std_error,a_list, LineWidth = lw-2, Color=wide_clr);
    plot(mean_error-std_error,a_list, LineWidth = lw-2, Color=wide_clr);


elseif plt_which == 1 || plt_which == 2  
    
    if plt_which == 1 %pm vs ropp
        xl = xlabel("(PM-M)/M * 100%");
        ttl = title("Bending Angle Error vs Impact Height, PM vs Model");

    elseif plt_which == 2 %go vs ropp
        xl = xlabel("(GO-M)/M * 100%");
        ttl = title("Bending Angle Error vs Impact Height, GO vs Model");
    end

    % because the ERA5 model output does not have consistently-spaced
    % impact parameters, we cannot use the same "ladder" approach as for
    % geometric optics. We must bin the points by impact height, and then
    % take the means and stds of the relative error in these height bins.
    % To do this, we must treat each point like its own thing, devoid of
    % occultation information.

    %make a matrix of all the bending angles and impact parameters padded
    %with NaNs AT THE END - just to make them all the same length 
    %actually don't need to do this??????

    % take the mean and std of each bin and put into array

    %set the lowest edge of bins to be 1/3 of the bin width below the
    %lowest impact parameter
    min_bin = min_height-bin_width_ropp/3;
    num_bins = fix((max_height - min_bin)/bin_width_ropp);

    % would be nice to be able to pass in s.impact_height and
    % s.relative_err
    % need to be able to take means and stds of everything
    % all_heights and all_errors are both long 1D arrays, but each height
    % and error still have the same index, so each point is preserved.
    [mean_error, std_error, edges] = binned_std(all_heights,all_errs,num_bins);

    %find the midpoints of the bins, because makes more sense to plot
    midpoints = edges(1:end-1)+diff(edges)/2;

    plot(mean_error,midpoints, LineWidth = lw, Color=wide_clr);
    plot(mean_error+std_error,midpoints, LineWidth = lw-2, Color=wide_clr);
    plot(mean_error-std_error,midpoints, LineWidth = lw-2, Color=wide_clr);


end

% save("model_plot_all_vars")

yl = ylabel("Impact Height (km)");
fontsize(xl,fs,'points')
fontsize(yl,fs,'points')
fontsize(ttl,fs,'points')

% disp(["COUNT IS",count])
% 
% disp("not pre-allocated")
toc
