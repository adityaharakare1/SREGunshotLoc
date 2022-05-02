clc;
clear;

% Muzzle
A=1; 	    				% peak amplitude
Tr_m = 0.2e-3;				% rise time
Td_m = 1e-3;				% decay time
Tdur_m = 10*Td_m;			% signal duration
Fs = 80000;
Ts = 1/Fs;
t_m = (0:Ts:Tdur_m-Ts)';	% time vector
xt1_muzzle = A*(t_m(t_m<=Tr_m)/Tr_m);
%plot(xt1_muzzle)
xt2_muzzle = A*(1 - (t_m(t_m>Tr_m)-Tr_m)/Td_m).*exp((Tr_m-t_m(t_m>Tr_m))/Td_m);
%plot(xt2_muzzle)
xt_muzzle = [xt1_muzzle;xt2_muzzle];

Nf=1024;
fk = (0:Nf-1)'*pi/Nf;
Xf_muzzle = freqz(xt_muzzle,1,fk);

% shock
B = 1;	    			    	% peak amplitude
Tr_s = 0.15e-3;			    	% rise time
Td_s = 0.15e-3;					% decay time
%Tdur_s = 2*Tr_s + Td_s+10*Ts;	% signal duration
Tdur_s = 0.01;
t_s = (0:Ts:Tdur_s-Ts)';	% time vector
xt1_shock = B*(t_s(t_s<=Tr_s)/Tr_s);
xt2_shock = B*(1 -  2*(t_s(t_s>Tr_s & t_s<=Tr_s+Td_s)-Tr_s)/Td_s);
xt3_shock = B*((t_s(t_s>Tr_s+Td_s & t_s<=2*Tr_s+Td_s) - Tr_s - Td_s)/Tr_s - 1);
xt4_shock = 0*t_s(t_s>2*Tr_s+Td_s);
xt_shock = [xt1_shock;xt2_shock;xt3_shock;xt4_shock];

Xf_shock = freqz(xt_shock,1,fk);

req_len = 1500;
k = zeros(1, req_len - length(xt_muzzle))';
xt_muzzle = [k; xt_muzzle];

k = zeros(1, req_len - length(xt_shock))';
xt_shock = [k; xt_shock];

% Plots
%plot(xt_muzzle)
t_m = (0:Ts:req_len/Fs-Ts)';
t_s = (0:Ts:req_len/Fs-Ts)';
f1 = figure;
plot(t_m/1e-3,xt_muzzle,'b'); grid on;
xlabel('Time (ms)'); ylabel('Amplitude');
title('Muzzle Blast');
axis([-inf inf -0.2 1.2])
f2 = figure;
plot(t_s/1e-3,xt_shock,'r'); grid on;
xlabel('Time (ms)'); ylabel('Amplitude');
title('Shock Wave');
axis([-inf inf -1.2 1.2])


% Noise Addition
muz_npower = -40;
swave_npower = -40;
muzzle_noise = wgn(req_len,1,muz_npower);
shock_noise = wgn(req_len,1,swave_npower);
xt_muzzle_n = xt_muzzle + muzzle_noise;
xt_shock_n = xt_shock + shock_noise;

f3 = figure;
plot(t_m/1e-3,xt_muzzle_n,'b'); grid on;
xlabel('Time (ms)'); ylabel('Amplitude');
title('Muzzle blast with Noise');
axis([-inf inf -0.2 1.2])

f4 = figure;
plot(t_s/1e-3,xt_shock_n,'b'); grid on;
xlabel('Time (ms)'); ylabel('Amplitude');
title('Shock Wave with Noise');
axis([-inf inf -1.2 1.2])

% Adding delay to the 4 microphones
f5 = figure;
plot(xt_muzzle_n(500:1199),'b');
hold on;
plot(xt_muzzle_n(440:1139),'r');
plot(xt_muzzle_n(540:1239),'g');
plot(xt_muzzle_n(600:1299),'y');
hold off;

% Saving data to mat file
data_muzzle = [xt_muzzle_n(500:1199),xt_muzzle_n(440:1139),xt_muzzle_n(540:1239),xt_muzzle_n(600:1299)];
data_muzzle = data_muzzle';
save('data_muzzle.mat','data_muzzle')
