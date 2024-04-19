function [a_out, alpha, v] = pm_v3_air(t, snr, a, theta, r_1, r_2, phi, k, n, norp)


% t = time;
% snr = amp_1;
% r_1 = r_gps;
% r_2 = r_reec;
% n = n_i_rec;
% norp = 1;

% figure;
% plot(t, snr);
% title("SNR vs time");
% 
% figure;
% plot(t, theta);
% title("Theta vs a")
% 
% figure;
% plot(t, phi);
% title("phi vs t");

disp(["t length", length(t)])

v = zeros(length(a),1);
rate_inp = 50;
% rate_inp = 1;
t_inp = t(1):(t(2)-t(1))/rate_inp:t(length(t));

% interp1(t,snr,t_inp); this doesn't work - need t and snr to be the same
% length
% interp1(t,phi,t_inp); this works

sig_inp = interp1(t,snr,t_inp).*exp(1j*interp1(t,phi,t_inp)); %amplitude * e^(i *phase)

% figure;
% scat1 = scatter(linspace(1, length(sig_inp),length(sig_inp)), sig_inp,10,linspace(1, length(sig_inp),length(sig_inp)));
% scat1.Marker = '.';
% title("Sig inp");

% sig_inp = snr.*exp(1j*phi');
% 
% re = real(sig_inp);
% im = imag(sig_inp);

% figure;
% plot(t, re);
% title("real");
% 
% figure;
% plot(t_inp, im);
% title("imaginary sig inp vs t inp");
% 
% figure;
% plot(t, phi);
% title("phi vs t");
% 
% figure;
% plot(t_inp, interp1(t,phi,t_inp));
% title("interpolated phi");

% sig_inp_re = snr.*

% figure;
% plot(t_inp, phi);
% title("phase vs t_inp");
% 
% figure;
% plot(t,snr);
% title("SNR");
% 
% figure;
% plot(t,sig_inp);
% title("sig inp");


% a_ind = 9000;
for a_ind=1:1:length(a)
%     disp(a(a_ind));

    if norp == 1  % negative elevation
        s = a(a_ind)*theta + sqrt(r_1.^2-a(a_ind).^2) + sqrt(n.^2*r_2.^2-a(a_ind).^2)...
            - a(a_ind)*(acos(a(a_ind)./r_1)+acos(a(a_ind)./(n.*r_2)));
    else          % positive elevation
        s = a(a_ind)*theta + sqrt(r_1.^2-a(a_ind).^2) - sqrt(n.^2*r_2.^2-a(a_ind).^2)...
            - a(a_ind)*(acos(a(a_ind)./r_1)-acos(a(a_ind)./(n.*r_2)));     
    end
    

% 
%    figure;
%    plot(t, a(a_ind)*theta);
%    plot(t, sqrt(r_1.^2-a(a_ind).^2))
%    plot(t, sqrt(n.^2*r_2.^2-a(a_ind).^2))
%    plot(t, -a(a_ind)*(acos(a(a_ind)./r_1)))
%    plot(t, acos(a(a_ind)./(n.*r_2)))
%    legend("a(a_ind)*theta","sqrt(r_1.^2-a(a_ind).^2)","sqrt(n.^2*r_2.^2-a(a_ind).^2)","-a(a_ind)*(acos(a(a_ind)./r_1))","acos(a(a_ind)./(n.*r_2))")

%     figure;
%     plot(t, imag(s));
%     title("imag s");


    % Make sure s is real
    if sum(imag(s))== 0
        s_inp = exp(-1j*interp1(t,k*s,t_inp));
        v(a_ind) = sum(sig_inp.*s_inp);
        
    else
        v(a_ind) = 0;
    end
        
    
end

% figure;
% plot(a);
% title("a");
% 
% figure();
% plot(t, theta);
% title("theta vs t");
% 
% figure();
% plot(t, r_1);
% title("r_1 vs t");

% figure();
% plot(t, r_2);
% title("r_2 vs t");

disp(["size of phi", size(phi)])
disp(["size of s", size(s)])

% figure;
% plot(t,s);
% title("s vs t");
% plot(t, phi);
% legend("s", "phi");

% figure;
% plot(t,phi);
% title("phi vs t")
% 
% figure;
% plot(t,s);
% title("s vs t");

% figure; %this looks good
% scatter(linspace(1, length(s_inp), length(s_inp)), s_inp,10,linspace(1, length(s_inp), length(s_inp)),'.');
% title("s_inp");

% figure;
% scatter(linspace(1, length(v), length(v)), v,10,linspace(1, length(v), length(v)),'.');
% title("v");

%     def running_mean(x, N):
%         cumsum = np.cumsum(np.insert(x, 0, 0))
%         return (cumsum[N:] - cumsum[:-N]) / N

    %ag_v  = sv(np.unwrap(np.angle(v)),199,3)
ag_v  = movmean(unwrap(angle(v)), 200);
alpha = (-1) * 1/k * diff(ag_v) / 1e-3;
a_out = a;

% figure
