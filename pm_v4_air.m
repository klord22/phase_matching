%function [a_out, alpha, v] = pm_v4_air(t, snr, a, theta, r_1, r_2, phi, k, n, norp) 
%
% Calculate bending angle profile using phase matching equations. Called in
% dop2alpha_pm_v4.m
%  
%% Inputs: 
% t, snr, a, theta, r_1, r_2, phi, k, n, norp 
% 
%   norp: means negative or positive elevation angle. 0 for positive, 1 for
%   negative
%
%% Outputs
%   a_out: impact parameter array
%   alpha: bending angle array (as function of impact parameter array)
%   v: from equation NUMBER ????
%

function [a_out, alpha, v] = pm_v4_air(t, snr, a, theta, r_1, r_2, phi, k, n, norp)

v = zeros(length(a),1);
rate_inp = 50; %v3
% rate_inp = 10; %works with rate >=4
t_inp = t(1):(t(2)-t(1))/rate_inp:t(length(t));

% interp1(t,snr,t_inp); this doesn't work - need t and snr to be the same
% length
% interp1(t,phi,t_inp); this works

sig_inp = interp1(t,snr,t_inp).*exp(1j*interp1(t,phi,t_inp)); %amplitude * e^(i *phase)

% fake_sig_inp = snr.*exp(1j*phi);

% disp(fake_sig_inp == sig_inp);
% 
% disp(["sig inp length", length(sig_inp)])

% sig_inp = snr.*exp(1j*phi');


for a_ind=1:1:length(a)
%     disp(a(a_ind));

    if norp == 1  % negative elevation
        s = a(a_ind)*theta + sqrt(r_1.^2-a(a_ind).^2) + sqrt(n.^2*r_2.^2-a(a_ind).^2)...
            - a(a_ind)*(acos(a(a_ind)./r_1)+acos(a(a_ind)./(n.*r_2)));
    else          % positive elevation
        s = a(a_ind)*theta + sqrt(r_1.^2-a(a_ind).^2) - sqrt(n.^2*r_2.^2-a(a_ind).^2)...
            - a(a_ind)*(acos(a(a_ind)./r_1)-acos(a(a_ind)./(n.*r_2)));     
    end
    

    % Make sure s is real
    if sum(imag(s))== 0
        s_inp = exp(-1j*interp1(t,k*s,t_inp));
        v(a_ind) = sum(sig_inp.*s_inp); % represents eq 24 or 25 (s_inp is neg, so sub)
        % but where is beginning term? sqrt? maybe just eq 24
%         fake_s_inp = exp(-1j*(k*s));
        
    else
        v(a_ind) = 0;
    end
        
    
end


%     def running_mean(x, N):
%         cumsum = np.cumsum(np.insert(x, 0, 0))
%         return (cumsum[N:] - cumsum[:-N]) / N

    %ag_v  = sv(np.unwrap(np.angle(v)),199,3)

%change this?
% adjust param, set at top

% ag_v  = movmean(unwrap(angle(v)), 2000);
ag_v  = movmean(unwrap(angle(v)), 200); %v3
alpha = (-1) * 1/k * diff(ag_v) / 1e-3;
a_out = a;

