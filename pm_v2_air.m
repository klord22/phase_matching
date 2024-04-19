function [a_out, alpha, v] = pm_v2_air(t, snr, a, theta, r_1, r_2, phi, k, n, norp)



v = zeros(length(a),1);
rate_inp = 50;
t_inp = t(1):(t(2)-t(1))/rate_inp:t(length(t));
sig_inp = interp1(t,snr,t_inp).*exp(1j*interp1(t,phi,t_inp));

a_ind = 9000;
% for a_ind=1:1:length(a)

    disp(a(a_ind));
    
    if norp == 1  % negative elevation
        s = a(a_ind)*theta + sqrt(r_1.^2-a(a_ind).^2) + sqrt(n.^2*r_2.^2-a(a_ind).^2)...
            - a(a_ind)*(acos(a(a_ind)./r_1)+acos(a(a_ind)./(n.*r_2)));
    else          % positive elevation
        s = a(a_ind)*theta + sqrt(r_1.^2-a(a_ind).^2) - sqrt(n.^2*r_2.^2-a(a_ind).^2)...
            - a(a_ind)*(acos(a(a_ind)/r_1)-acos(a(a_ind)/(n.*r_2)));     
    end
    
    % Make sure s is real
    if sum(imag(s))== 0
        s_inp = exp(-1j*interp1(t,k*s,t_inp));
        v(a_ind) = sum(sig_inp.*s_inp);
        
    else
        v(a_ind) = 0;
    end
        
    
% end

%     def running_mean(x, N):
%         cumsum = np.cumsum(np.insert(x, 0, 0))
%         return (cumsum[N:] - cumsum[:-N]) / N

    %ag_v  = sv(np.unwrap(np.angle(v)),199,3)
    ag_v  = movmean(unwrap(angle(v)), 200);
    alpha = (-1) * 1/k * diff(ag_v) / 1e-3;
    a_out = a;
