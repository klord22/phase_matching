function [time, x_rec, y_rec, z_rec, u_rec, v_rec, w_rec, x_gps, y_gps, z_gps, u_gps, v_gps, w_gps, ex_ph, loss, ex_dop] = read_ar_v2(file_name)

fid = fopen(file_name);
if (fid == -1)
    error('Data file "%s" not found.  Check current directory or path.', file_name);
end

fprintf('Reading data from %s...\n', file_name);

i     =1;
line    = fgets(fid);
while ~feof(fid)
    
    line          = fgets(fid);
    C             = strsplit(line);
    time(i)       = str2double(C(2));
    x_rec(i)      = str2double(C(3)); 
    y_rec(i)      = str2double(C(4)); 
    z_rec(i)      = str2double(C(5)); 
    u_rec(i)      = str2double(C(6)); 
    v_rec(i)      = str2double(C(7)); 
    w_rec(i)      = str2double(C(8));
    x_gps(i)      = str2double(C(9)); 
    y_gps(i)      = str2double(C(10)); 
    z_gps(i)      = str2double(C(11)); 
    u_gps(i)      = str2double(C(12)); 
    v_gps(i)      = str2double(C(13)); 
    w_gps(i)      = str2double(C(14)); 
    ex_ph(i)      = str2double(C(15)); 
    loss(i)       = str2double(C(16)); 
    ex_dop(i)     = str2double(C(17)); 
    
    i=i+1;
end

if file_name(8)=='r' && time(2)-time(1)>0
    time   = fliplr(time);
    x_rec  = fliplr(x_rec);
    y_rec  = fliplr(y_rec);
    z_rec  = fliplr(z_rec);
    u_rec  = fliplr(u_rec);
    v_rec  = fliplr(v_rec);
    w_rec  = fliplr(w_rec);
    x_gps  = fliplr(x_gps);
    y_gps  = fliplr(y_gps);
    z_gps  = fliplr(z_gps);
    u_gps  = fliplr(u_gps);
    v_gps  = fliplr(v_gps);
    w_gps  = fliplr(w_gps);
    ex_ph  = fliplr(ex_ph);
    loss   = fliplr(loss);
    ex_dop = fliplr(ex_dop);
    
end

    


    