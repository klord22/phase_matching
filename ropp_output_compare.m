%compare Pawel's and my ropp output 

address_pawel = '~/ags/projects/hiaper/2022.305_ar2023/ropp_bangle_2d/2023.015_iop16/';
% address_kate = '~/ags/projects/hiaper/arivers/kate/code_mat_ropp_aro/2023.015_iop16/';
address_kate = '~/ags/projects/hiaper/arivers/kate/code_mat_ropp_aro/2d_no_drift/';
% dir_pawel = dir(strcat(address_pawel,'*.mat'));
% dir_kate = dir(strcat(address_kate,'*.mat'));
dir_pawel = dir(strcat(address_pawel,'*.nc'));
dir_kate = dir(strcat(address_kate,'*.nc'));

info_path = '~/ags/projects/hiaper/2022.305_ar2023/nret/2023.015_iop16/';

for file = dir_kate' %bc kate has fewer files
    occ_name = lower(file.name(26:28));
    disp(occ_name) 
    info_file_path = [info_path,occ_name,'*'];
    info_dir = dir(info_file_path);
    folder_name = info_dir.name;

    yaml_data = yaml.loadFile([info_path,folder_name,'/info.yaml']);
    Rc = yaml_data.Rc;

    % load(strcat(string(file.folder),'/',string(file.name)),"ro_data")
    % bangle_kate = ro_data.bangle; %make sure these are the correct ones (ie not obs or whatever)
    % impact_kate = ro_data.impact;
    % clear ro_data

    bangle_kate = ncread(strcat(address_kate,string(file.name)), 'bangle');
    impact_kate = ncread(strcat(address_kate,string(file.name)), 'impact');

    bangle_pawel = ncread(strcat(address_pawel,string(file.name)), 'bangle');
    impact_pawel = ncread(strcat(address_pawel,string(file.name)), 'impact');

    % load(strcat(address_pawel,string(file.name)),"ro_data")

    % bangle_diff = bangle_kate - ro_data.bangle(1:length(bangle_kate));

    figure
    hold on
    % plot(bangle_kate,impact_kate)
    % plot(ro_data.bangle, ro_data.impact)
    % scatter(bangle_kate,impact_kate, "filled")
    % scatter(ro_data.bangle, ro_data.impact, "filled")
    % plot(bangle_diff, ro_data.impact)

    impact_kate = impact_kate*10^(-3) - Rc;
    impact_pawel = impact_pawel*10^(-3) - Rc;
    
    scatter(bangle_pawel*10^3,impact_pawel, "filled")
    scatter(bangle_kate*10^3,impact_kate, "filled")

    title("BA Profile from different runs of ropp2d, ",occ_name)
    xlabel("Bending Angle (mrad)")
    ylabel("Impact Height (km)")
    legend("Kate","Pawel")
    clear ro_data

end
