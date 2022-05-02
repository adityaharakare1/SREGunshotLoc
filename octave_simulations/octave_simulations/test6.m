% for previous test, let's find the estimation error.

clc;
clear;
pkg load communications;

% muzzle
A=14; 					% peak amplitude
Tr_m=10e-6;				% rise time
Td_m=0.4e-3;					% decay time
Tdur_m=10*Td_m;				% signal duration
Ts=5e-6;				% sampling interval
Fs = 1/Ts;
t_m = (0:Ts:Tdur_m-Ts)';	% time vector
xt1_muzzle = A*(t_m(t_m<=Tr_m)/Tr_m);
xt2_muzzle = A*(1 - t_m(t_m>Tr_m)/Td_m).*exp(-t_m(t_m>Tr_m)/Td_m);
xt_muzzle = [xt1_muzzle;xt2_muzzle];

% FT characteristics
Nf=1024;
fk = (0:Nf-1)'*pi/Nf;
Xf_muzzle = freqz(xt_muzzle,1,fk);

% anti-talk high-pass filter characteristics
Nord=4;									% filter order
Wc = 4.4e3/Fs*2*pi;						% cut-off
[b_talk,a_talk] = butter(Nord,Wc,'high');			% design highpass filter
xt_muzzle_waveform = filter(b_talk,a_talk,xt_muzzle);
Xf_muzzle_waveform = freqz(xt_muzzle_waveform,1,fk);

% define matched-filter {h(t) = x*(T-t)}
% now we have anti-talk filter, so the matched-filter waveform changes
h_matched_muzzle = flipud(conj(xt_muzzle_waveform))/norm(xt_muzzle_waveform);

% reverb characteristics
D=1;				% reflection interval (in samples)
decay=1;			% decay
Necho=1;			% only 1 reflection considered
b_rev=upsample((1-decay).^(0:Necho)',D+1);		% reverb filter

% array parameters
fc=100e3;				% maximum freq component in the signal
Pmax=20;			% defines the area in which sensors are deployed
c=343;				% speed of sound
tht_target=0.25*pi;			% we are interested in only planar angle location and not elevation
R_target=20;					% target range (range is considered from reference)

% we are taking into consideration that the source can be within the max
% distance range. hence, the angle is different for all sensors.
steer = [cos(tht_target);sin(tht_target)];	% target direction vector
target_loc = R_target*steer;

t = (0:Ts:0.3-Ts)';			% 300 ms capturing time = max range: ~100 m
Nt = numel(t);
nvar=-10;						% noise power (db)

Nvec=4:10;
err_loc_estimate = zeros(numel(Nvec),1);
for(k3=1:numel(Nvec))
	Nelements=Nvec(k3);
	%P = [[0;0],[0;1],[1;0],[1;1]]*Pmax;
	P = [[0;0],rand(2,Nelements-1)*Pmax];

	% signals captured using sensors
	dR = sqrt(sum((P - R_target*steer).^2,1));			% target distance from each sensor
	tau_muzzle_R = dR/c;
	tau_muzzle = tau_muzzle_R;
	for(k1=1:Nelements)
		muzzle_noise = wgn(1,Nt,nvar);
		h_muzzle = A/sqrt(4*pi*dR(k1)^2)*sinc(fc*(t - tau_muzzle(k1)));
		x_delayed = filter(xt_muzzle,1,h_muzzle)';
		x_reverb = filter(b_rev,1,x_delayed);
		x_muzzle(k1,:) = filter(b_talk,a_talk,x_reverb) + muzzle_noise;
	end
	x_received = x_muzzle;

	% find peaks of matched-filter output
	for(k2=1:Nelements)
		y_matchedOut_muzzle(k2,:) = filter(h_matched_muzzle,1,x_received(k2,:));
		[~,max_muzzle_ind(k2)] = max(abs(y_matchedOut_muzzle(k2,:)));
	end

	% estimated delay
	% practically, the delays are measured compared to the minimum-delay
	% element, so we modify all the delays accordingly. we know which sensor
	% received the signal first.
	[min_ind,ref_sensor] = min(max_muzzle_ind);
	muzzle_ind_shifted = max_muzzle_ind - min_ind;

	tau_estim_muzzle = muzzle_ind_shifted*Ts;

	% range-sphere estimates
	R_estim_muzzle = c*tau_estim_muzzle;

	% localise the target (solve Vx=y for x)
	V = 2*[P(1,:) - P(1,ref_sensor); P(2,:) - P(2,ref_sensor); R_estim_muzzle]';
	y = (sum(P.^2) - sum(P(:,ref_sensor).^2) - R_estim_muzzle.^2).';
	V(ref_sensor,:)=[]; y(ref_sensor)=[];
	x = V\y;
	target_loc_estim = x(1:2);
	err_loc_estimate(k3) = norm(target_loc_estim - target_loc);
end


% plots
%subplot(2,2,1);
%plot(t/1e-3,x_received); grid on;
%xlabel('Time (ms)'); ylabel('Received signal amplitude at each sensor');
%
%subplot(2,2,2);
%plot(t/1e-3,y_matchedOut_muzzle); grid on;
%xlabel('Time (ms)'); ylabel('Matched-filter output for muzzle blast');
%
%subplot(2,2,3:4);
%plot(P(1,:),P(2,:),'bx',R_target*steer(1),R_target*steer(2),'ro',target_loc_estim(1),target_loc_estim(2),'kd'); grid on;
%xlabel('X'); ylabel('Y');
%xlim([-2,Pmax*1.1]); ylim([-2,Pmax*1.1]);
%legend('Element locations','True target location','Estimated location','Location','Southwest');
plot(Nvec,err_loc_estimate/1e-3); grid on;
xlabel('Number of sensors'); ylabel('Range error (mm)');








