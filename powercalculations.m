
clear
close all
out = load("rigidbodyforces_withecc.mat");
out_no_exc = load("rigidbodyforces_noecc.mat");

%extra gegevens definiëren
R_0 = (out.bcr + out.rof)*10^(-3);
N = out.normalforce_tot;
alpha = out.pressure_angle;
exc = out.exc*10^(-3);
S = out.S * 10^(-3);
R = sqrt(out.xpitch.^2 + out.ypitch.^2)*10^(-3);


% zonder eccentriciteit
R_0_no_exc = (out_no_exc.bcr + out_no_exc.rof)*10^(-3);
N_no_exc = out_no_exc.normalforce_tot;
alpha_no_exc = out_no_exc.pressure_angle;
S_no_exc = out_no_exc.S * 10^(-3);
R_no_exc = sqrt(out_no_exc.xpitch.^2 + out_no_exc.ypitch.^2)*10^(-3);


%test
% R_test = sqrt(out.xpitch.^2 + out.ypitch.^2)*10^(-3);
% R_err = R_test - R_0 - S;
% %plot(out.theta, R_err);
% DR = diff(R);

% vermogen met eccentriciteit = 0
P_no_exc2 = N_no_exc.*sin(alpha_no_exc).*(R_no_exc).*out_no_exc.w;

%vermogen met gekozen eccentriciteit
gamma = atan(exc ./ R);
P_exc2 = N.*cos(alpha).*((sqrt(R_0^2 - exc^2) + S).*tan(alpha) + exc)*out.w;

%fout op vermogen
P_err2 = P_no_exc2 - P_exc2;
figure
subplot(311)
plot(out.thetadegree, P_no_exc2)
subplot(312)
plot(out.thetadegree, P_exc2)
subplot(313)
plot(out.thetadegree, P_err2)

subplot(3,1,1)
xlim([3 403])
ylim([-700 700])
xlabel('theta (°)')
ylabel('Power no ecc (W)')
subplot(3,1,2)
xlim([3 403])
ylim([-700 700])
xlabel('theta (°)')
ylabel('Power with ecc(W)')
subplot(3,1,3)
xlabel('theta (°)')
ylabel('Error (\)')



%gemiddeld vermogen
P_av2 = mean(P_exc2)

out.w % zijn dus allebei cte
out.rpm
out.S
Rb = out.bcr


K = 0.1
M = P_exc2./ out.w
M_av = P_av2 ./ out.w
M_delta = M_av - M
% dus altijd wnnr M > M_av vertraagt de cam en wnnr M < M_av versnelt de
% cam. M is hier het lastmoment dat verandert in de tijd en M_av stelt het
% veronderstelde geleverde vermogen door de motor voor. 

figure
plot(out.theta, M_delta, out.theta,ones(size(out.theta)) * M_av)

figure
plot(out.thetadegree, M, out.thetadegree,ones(size(out.theta)) * M_av)
% 
xlabel('cam angle (degree)')
ylabel('torque (Nm)')
legend({'load torque','average'})
%speed variation berekenen met formule voor A 
A = cumtrapz(out.theta, M_delta)


figure
plot(out.thetadegree, A)

xlabel('theta (degree)')
ylabel('A (J)')
%maximum work surplus berekenen: 
[A_max,I_Amax] = max(A)
[A_min,I_Amin] = min(A)
A_maxvar = A_max - A_min
B2 = out.thetadegree(I_Amin)
B1 = out.thetadegree(I_Amax)



A_MAX = trapz(out.theta(I_Amin:size(A,2)), M_delta(I_Amin:size(A,2))) + trapz(out.theta(1:I_Amax), M_delta(1:I_Amax));


%massatraagheidsmoment berekenen
I = A_maxvar / (K*out.w^2)

%vliegwiel dimensioneren (schijftype)
rho =   7700 % kg/m^3 voor stainless steel
r = 110*10^(-3) % ietsje kleiner dan basisstraal
t = 2*I / (pi*rho*r^4)

%I opnieuw berekenen
I2 = pi*r^4*t*rho/2
