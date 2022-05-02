% this script shows a different setup then previous one. here,
% we make use of target range as well, hence the magnitude of the
% received signal is dependent on range, and the delay experienced by
% the signal at each sensor is due to both the distance and angle. moreover,
% we see only muzzle-blast signature here.

clc;
clear;
pkg load communications;

% muzzle
A=1; 					% peak amplitude
Tr_m=10e-6;				% rise time
Td_m=500e-6;					% decay time
Tdur_m=4*Td_m;				% signal duration
Ts=5e-6;				% sampling interval
Fs = 1/Ts;
t_m = (0:Ts:Tdur_m-Ts)';	% time vector
xt1_muzzle = A*(t_m(t_m<=Tr_m)/Tr_m);
xt2_muzzle = A*(1 - t_m(t_m>Tr_m)/Td_m).*exp(-t_m(t_m>Tr_m)/Td_m);
xt_muzzle = [xt1_muzzle;xt2_muzzle];

Nf=1024;
fk = (0:Nf-1)'*pi/Nf;
Xf_muzzle = freqz(xt_muzzle,1,fk);

% define matched-filter {h(t) = x*(T-t)}
h_matched_muzzle = flipud(conj(xt_muzzle));

% array parameters
fc=100e3;				% maximum freq component in the signal
Pmax=20;			% defines the area in which sensors are deployed
c=343;				% speed of sound
tht_target=0.2*pi;			% we are interested in only planar angle location and not elevation
R_target=50;					% target range (range is considered from reference)

% we are taking into consideration that the source can be within the max
% distance range. hence, the angle is different for all sensors.

steer = [cos(tht_target);sin(tht_target)];	% target direction vector
t = (0:Ts:0.3-Ts)';			% 20 ms capturing time
Nt = numel(t);
nvar=-50;						% SNR (db)
N=10;
%P = [[0;0],[0;1],[1;0],[1;1]]*Pmax;
P = [[0;0],rand(2,N-1)*Pmax];

% signals captured using sensors
dR = sqrt(sum((P - R_target*steer).^2,1));			% target distance from each sensor
tau_muzzle_R = dR/c;
tau_muzzle = tau_muzzle_R;

for(k1=1:N)
	muzzle_noise = wgn(1,Nt,nvar);
	h_muzzle = A/sqrt(4*pi*dR(k1)^2)*sinc(fc*(t - tau_muzzle(k1)));
	x_muzzle(k1,:) = filter(xt_muzzle,1,h_muzzle)' + muzzle_noise;
end
x_received = x_muzzle;

% find peaks of matched-filter output
for(k2=1:N)
	y_matchedOut_muzzle(k2,:) = filter(h_matched_muzzle,1,x_received(k2,:));
	[~,max_muzzle_ind(k2)] = max(abs(y_matchedOut_muzzle(k2,:)));
end

% estimated delay
tau_estim_muzzle = max_muzzle_ind*Ts;

% range-sphere estimates
R_estim_muzzle = c*tau_estim_muzzle;

% localise the target (solve Vx=y for x)
V = [-ones(N,2),2*P'];
y = (sum(P.^2) - R_estim_muzzle.^2).';
x = V\y;
target_loc_estim = x(3:4);

% plots
subplot(2,2,1);
plot(t/1e-3,x_received); grid on;
xlabel('Time (ms)'); ylabel('Received signal amplitude at each sensor');

subplot(2,2,2);
plot(t/1e-3,y_matchedOut_muzzle); grid on;
xlabel('Time (ms)'); ylabel('Matched-filter output for muzzle blast');

subplot(2,2,3:4);
plot(P(1,:),P(2,:),'bx',R_target*steer(1),R_target*steer(2),'ro',target_loc_estim(1),target_loc_estim(2),'kd'); grid on;
xlabel('X'); ylabel('Y');
%xlim([-2,Pmax*1.1]); ylim([-2,Pmax*1.1]);
legend('Element locations','True target location','Estimated location','Location','Southwest');



















