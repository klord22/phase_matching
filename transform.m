% Transform the olchanlook.mat to geom file format
% So that we can use new version phase matching to process it (AR2015 and AR2018)
% Eric Wang 2018 Oct

function transform(rf, ch, prn_occ, prn_h, rors_occ)

exph_path = sprintf('/ags/projects/hiaper/2010.215_predict/exphs_ol');
olpd_path = sprintf('%s/%s/olpredict/CH%d_cli', exph_path, rf, ch);   % olpredict Path
olcl_path = sprintf('%s/%s/olchanlook/CH%d_cli', exph_path, rf, ch);  % olchanlook Path
olhe_path = sprintf('%s/%s/olchanlook/CH1', exph_path, rf);  % High elevation olchanlook Path

%if rf < 10 && rf > 0
%    rf_str = sprintf('0%d',rf);
%else
%    rf_str = sprintf('%d',rf);
%end
if prn_occ < 10 && prn_occ > 0
    prn_occ_str = sprintf('0%d',prn_occ);
else
    prn_occ_str = sprintf('%d',prn_occ);
end
if prn_h < 10 && prn_h > 0
    prn_h_str = sprintf('0%d',prn_h);
else
    prn_h_str = sprintf('%d',prn_h);
end


geofile_air       = sprintf('%s/geom_apx_gps_shift_prn%d%s', olpd_path, prn_occ, rors_occ);
geofile_air_ol    = sprintf('%s/geom_apx_gps_shift_ol_prn%d%s', olcl_path, prn_occ, rors_occ);
geofile_air_ol    = sprintf('%s/prn_g%s%s_g%s_geom.txt', olcl_path, prn_occ_str, rors_occ, prn_h_str);
load(sprintf('%s/OLchanlookSV%2s%s', olcl_path, prn_occ_str, rors_occ));
gps_sow_vecR_o    = gps_sow_vecR;
CorrCorrectAmpR_o = ex_CorrCorrectAmpR;
CorrCorrectPhaseR_o = ex_CorrCorrectPhaseR;
load(sprintf('%s/OLchanlookSV%2sh', olhe_path, prn_h_str));


fid = fopen(geofile_air);
if (fid == -1)
    error('Data file "%s" not found.  Check current directory or path.', datafilename);
end
i     =1;
line    = fgets(fid);
while ~feof(fid)
    line    = fgets(fid);
    time(i) = str2double(line(1:9));
    i = i + 1;
end


gps_sow_min       = min(time);
gps_sow_max       = max(time);

if min(gps_sow_vecR_o)>ceil(gps_sow_min)
    str_ind_o = 26;
    gps_sow_min = ceil(min(gps_sow_vecR_o));
else
    str_ind_o = find(gps_sow_vecR_o==ceil(gps_sow_min));
end
if min(gps_sow_vecR)>ceil(gps_sow_min)
    str_ind_h = 26;
else
    str_ind_h = find(gps_sow_vecR==ceil(gps_sow_min));
end
if max(gps_sow_vecR_o)<floor(gps_sow_max)+1
    end_ind_o = length(gps_sow_vecR_o)-24;
    gps_sow_max = floor(max(gps_sow_vecR_o))-1;
else
    end_ind_o = find(gps_sow_vecR_o==floor(gps_sow_max)-1);
end
if max(gps_sow_vecR)<floor(gps_sow_max)
    end_ind_h = length(gps_sow_vecR)-24;
else
    end_ind_h = find(gps_sow_vecR  ==floor(gps_sow_max));
end

%str_ind_o         = find(gps_sow_vecR_o==ceil(gps_sow_min));
%end_ind_o         = find(gps_sow_vecR_o==floor(gps_sow_max));
%str_ind_h         = find(gps_sow_vecR  ==ceil(gps_sow_min));
%end_ind_h         = find(gps_sow_vecR  ==floor(gps_sow_max));
for i=1:1:gps_sow_max-gps_sow_min+1
    amp(i)        = mean(CorrCorrectAmpR_o(str_ind_o+(i-1)*50-25:str_ind_o+(i-1)*50+24));
    phase_ex(i)   = -1*(mean(CorrCorrectPhaseR_o(str_ind_o+(i-1)*50-25:str_ind_o+(i-1)*50+24))...
                       -mean(CorrCorrectPhaseR(str_ind_h+(i-1)*50-25:str_ind_h+(i-1)*50+24)));
end
dp_ex = gradient(phase_ex);

k = (1.57542e+9/299792.458)*2*pi;   % rad/km
phase_ex = phase_ex/k*1000;   % m
dp_ex = dp_ex/k*1000;         % m/s

fid  = fopen(geofile_air);
fid2 = fopen(geofile_air_ol,'w');
i     =1;
line    = fgets(fid);
fprintf(fid2, line);
while ~feof(fid)
    line    = fgets(fid);
    if str2double(line(1:9))>=gps_sow_min && str2double(line(1:9))<=gps_sow_max
        fprintf(fid2, '%s  %15.6f  %15.6f  %11.6f\n', line(1:184), phase_ex(i), amp(i), dp_ex(i));
        
        i = i + 1;
    end

end

    

%fprintf(i);











