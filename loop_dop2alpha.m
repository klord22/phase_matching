function loop_dop2alpha(path)

%hard-coded:
occ_start = 5;
occ_end = 11;
segment_start = 13;
segment_end = 27;
name_segment = 'phase_shift_cor';
info_path = '~/ags/projects/hiaper/2022.305_ar2023/nret/2023.015_iop16/';
% info_path = '/ags/projects/hiaper/2022.305_ar2023/nret/2023.015_iop16/';
n_surf = 1.0003;

% NOTE: This is written assuming the relevant prn files are in a folder
% called input. This shoud be changed to grab those files from their
% original directory, but need to weed out the "bad" ones!

% Currently not working because the file name needs to have "input/" in
% front of it, but that freaks out dop2alpha


%general coded

directory = dir(path);

tic

for file = directory'

    % get name of excess phase file
    name = file.name;
    
    % weed out random short names
    if length(name) >= segment_end

        % make sure it is an excess phase file
        if name(segment_start:segment_end) == name_segment

            % get name of occultation folder
            occ_name = name(occ_start:occ_end);
            info_file_path = [info_path,occ_name,'*'];
            info_dir = dir(info_file_path);
            disp("info_file_path")
	        disp(info_file_path)
	        disp("info_dir")
	        disp(info_dir)
	        disp("name")
	        disp(info_dir(1).name)%check this if problems
	        folder_name = info_dir.name;
              
            data = yaml.loadFile([info_path,folder_name,'/info.yaml']);

            % get and calculate needed quantities for pm function
            Nrec = data.Nrec;
            Rc = data.Rc;

            amin = Rc*n_surf;

            disp([Nrec, amin]) %for debugging dop2alp purposes

            % run phase matching for positive and negative bending angles

            % function dop2alpha_pm_v4(file_name, nrec, amin,norp)

%             dop2alpha_pm_v4(['input/',name], Nrec, amin, 0, false)
            dop2alpha_pm_v4(['input/',name], Nrec, amin, 1, false)

        end

    end
end

toc


