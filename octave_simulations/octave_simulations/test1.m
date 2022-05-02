% this script shows the simulated gunshot signals for muzzle-blast 
% and balistic shock waves.
% muzzle-blast is simulated using 'Friedlander wave' equation and shock-wave 
% is generated using 'N-wave' equation. 
% then we design an arbitrary array (sensors at random locations) and
% simulate the signal received by all the sensors.
% the captured signal is then passed through matched-filters and then the AoA
% is estimated.

clc;
clear;
pkg load communications;

% muzzle
A=1; 					% peak amplitude
Tr_m=10e-6;				% rise time
Td_m=50e-6;					% decay time
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

% shock
B=1;					% peak amplitude
Tr_s=10e-6;				% rise time
Td_s=50e-6;					% decay time
Tdur_s=2*Tr_s + Td_s+10*Ts;				% signal duration
t_s = (0:Ts:Tdur_s-Ts)';	% time vector
xt1_shock = B*(t_s(t_s<=Tr_s)/Tr_s);
xt2_shock = B*(1 -  2*t_s(t_s>Tr_s & t_s<=Td_s)/Td_s);
xt3_shock = B*((t_s(t_s>Td_s & t_s<=2*Tr_s+Td_s) - Tr_s - Td_s)/Tr_s - 1);
xt4_shock = 0*t_s(t_s>2*Tr_s+Td_s);
xt_shock = [xt1_shock;xt2_shock;xt3_shock;xt4_shock];

Xf_shock = freqz(xt_shock,1,fk);

% define matched-filter {h(t) = x*(T-t)}
h_matched_muzzle = flipud(conj(xt_muzzle));
h_matched_shock = flipud(conj(xt_shock));

% array parameters
fc=100e3;				% maximum freq component in the signal
Pmax=5;			% defines the area in which sensors are deployed
c = 343;				% speed of sound
c_muzzle = c;
Mach=3;					% mach number
c_shock = Mach*c;
tht=0.25*pi;			% we are interested in only planar angle location and not elevation
steer = [cos(tht);sin(tht)];
t = (0:Ts:0.02-Ts)';			% 20 ms capturing time
Nt = numel(t);
nvar=-50;

Nvec=4:10;
for(k3=1:numel(Nvec))
	N=Nvec(k3);
	% sensor locations (first is always reference)
	% the 'referece' is considered a central node farthest from the source
	% and at origin. all other nodes measure delay w.r.t. the central node.
	%P = [[0;0],[0;1],[1;0],[1;1]]*Pmax;
	P = [[0;0],rand(2,N-1)*Pmax];

	% signals captured using sensors
	tau_muzzle = steer'*P/c_muzzle;
	tau_shock = steer'*P/c_shock;

	for(k1=1:N)
		muzzle_noise = wgn(1,Nt,nvar);
		h_muzzle = sinc(fc*(t - tau_muzzle(k1)));
		x_muzzle(k1,:) = filter(xt_muzzle,1,h_muzzle)' + muzzle_noise;
		
		shock_noise = wgn(1,Nt,nvar);
		h_shock = sinc(fc*(t - tau_shock(k1)));
		x_shock(k1,:) = filter(xt_shock,1,h_shock)' + shock_noise;
	end
	x_received = x_muzzle + x_shock;

	% processing
	% find peaks of the matched-filter output
	max_muzzle_ind=zeros(N,1);
	max_shock_ind=zeros(N,1);
	for(k2=1:N)
		y_matchedOut_muzzle(k2,:) = filter(h_matched_muzzle,1,x_received(k2,:));
		y_matchedOut_shock(k2,:) = filter(h_matched_shock,1,x_received(k2,:));
		[~,max_muzzle_ind(k2)] = max(abs(y_matchedOut_muzzle(k2,:)));
		[~,max_shock_ind(k2)] = max(abs(y_matchedOut_shock(k2,:)));
	end

	% estimated delay
	tau_estim_muzzle = max_muzzle_ind*Ts;
	tau_estim_shock = max_shock_ind*Ts;

	% least-square solution of steering vector
	steer_estim_muzzle = (c_muzzle*tau_estim_muzzle)\(P');
	steer_estim_shock = (c_shock*tau_estim_shock)\(P');
	tht_estim_muzzle = arg(complex(steer_estim_muzzle(1),steer_estim_muzzle(2)));
	tht_estim_shock = arg(complex(steer_estim_shock(1),steer_estim_shock(2)));

	% location error
	err_tht_estimate(k3) = abs(tht_estim_muzzle - tht);
end

% plots
subplot(2,3,1);
plot(t_m/1e-3,xt_muzzle,'b',t_s/1e-3,xt_shock,'r--'); grid on;
xlabel('Time (ms)'); ylabel('Amplitude');
legend('muzzle blast','shock wave');

subplot(2,3,2);
plot(t/1e-3,x_received); grid on;
xlabel('Time (ms)'); ylabel('Received signal amplitude at each sensor');

subplot(2,3,3);
plot(P(1,:),P(2,:),'x'); grid on;
hold on;
quiver(0,0,steer(1),steer(2),1);
quiver(0,0,steer_estim_muzzle(1),steer_estim_muzzle(2),1);
quiver(0,0,steer_estim_shock(1),steer_estim_shock(2),1);
hold off;
xlabel('X'); ylabel('Y');
xlim([-2,Pmax*1.1]); ylim([-2,Pmax*1.1]);
legend('Element locations','True target direction','Estimated direction (muzzle)','Estimated direction (shock)','Location','Southwest');

subplot(2,3,4);
plot(t/1e-3,y_matchedOut_muzzle); grid on;
xlabel('Time (ms)'); ylabel('Matched-filter output for muzzle blast');

subplot(2,3,5);
plot(t/1e-3,y_matchedOut_shock); grid on;
xlabel('Time (ms)'); ylabel('Matched-filter output for shockwave');

subplot(2,3,6);
plot(Nvec,err_tht_estimate*180/pi,'b'); grid on;
xlabel('Number of elements'); ylabel('Angle estimation error (\circ)');







