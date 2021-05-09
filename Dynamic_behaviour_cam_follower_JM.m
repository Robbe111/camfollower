%%%%%%%%%%%%%%%    Dynamic behaviour cam-follower     %%%%%%%%%%%%%%%%%

%Matthijs Coninx
%Jarne De Nies
close all
clear all
out = load('rigidbodyforces_withecc.mat');

%%%%%%%%%%%%%%%     Single Rise Analysis      %%%%%%%%%%%%%%%%%%%%

%%% Parameters
zeta = 0.058;            %[-], given
lambda = 0.75/zeta;           %[-], condition Single Rise requires lambda > 8.33333, dimensionless resonant frequency
beta = 20;              %[degrees], critical rise
w = out.w;              %[rad/s]
t1 = 2*pi*beta/(360*w); %[s]

time = out.theta(20000:22000)/out.w;
t =  time-time(1);



tau_sr = t/t1;

theta_sr = ((out.S(20000:22000)*0.001)-0.0125)/0.0025;%[-], dimensionless follower input motion

%%% numerical analysis
num = (2*pi*lambda)^2;
den = [1, 2*zeta*(2*pi*lambda), (2*pi*lambda)^2];

theta0 = 1; % vul hier zelf de initiele dimensieloze heffing in
theta_dot0 = (-0.25*1e-3*180/pi*out.w)/0.0025; % vul hier zelf de initiele dimensieloze snelheid in

[A,B,C,D] = tf2ss(num,den);
X0 = [1/C(2)*theta_dot0; 1/C(2)*theta0];
gamma_sr = lsim(A,B,C,D, theta_sr, tau_sr, X0);

%calculation of amplitude exponential envelope
x0 = gamma_sr(2000)-1;                             %[-], gamma(tau = 1) -1
xd0 = (gamma_sr(2001)-gamma_sr(1999))/(2/2000);       %[-], numerical derivative, 2*step = 2/3500
lambdad = lambda * sqrt(1-zeta^2);              %[-] 
A1 = sqrt(((x0 * (2*pi*lambdad))^2 + (xd0+zeta*(2*pi*lambda)*x0)^2)/(2*pi*lambdad)^2) %formula 6b slide 13 H9

figure('Name','Single Rise: Dimensionless follower input and output motion')
plot(tau_sr,theta_sr);
hold on
plot(tau_sr,gamma_sr);
hold off
title('Single Rise: Dimensionless follower input and output motion')
xlabel('\tau [-]')
ylabel('\theta(\tau) and \gamma_{num}(\tau) [-]')
legend('Input','Output')

figure('Name', 'Single Rise: Deformation of the follower')
subplot(1,2,1)
plot(tau_sr, theta_sr - gamma_sr)
title('Single Rise: Dimensionless deformation of the follower')
xlabel('\tau [-]')
ylabel('\theta(\tau) - \gamma_{num}(\tau) [-]')

subplot(1,2,2)
plot(tau_sr, 10*(theta_sr - gamma_sr))
title('Single Rise: Deformation of the follower in mm')
xlabel('\tau [-]')
ylabel('\theta(\tau) - \gamma_{num}(\tau) [mm]')

%%% Approximate analysis
Q = (2*pi)^2;                                       % [-], for cycloidal cams
N = 3;                                              % [-], for cycloidal cams
A1_tilde = Q/(2*pi*lambda)^N * sqrt(1/(1-zeta^2))   % [-], amplitude for exponential envelope 
rel_error = abs((A1-A1_tilde)/A1)                   % [-], relative error between numerical and approximate solution





%%%%%%%%%%%%%%%%%   Multi Rise Analysis    %%%%%%%%%%%%%%%%%%%%%%%%


%%% parameters
T = 2*pi/w;                         %[s], period of cam rotation
step = 1/length(out.thetadegree);   % length vectors
tn = t1/lambda;                     %[s], same as for single rise analysis 
lambda_tilde = T/tn;                %[-], number of periods of the undamped resonance frequency during one rotation of the cam
wn = 2*pi /tn;                      %[rad/s], undamped resonace frequency
m = 21;                             %[kg] equivalent mass
kf = m* wn^2 / 1000                %[N/mm], equivalent stiffness follower

M = 25;                              % number of rotations
tau = [0:step:M-step];               %[-], non dimensional time, M rotations
theta = repmat(out.S/25,1,M);        % [-], follower input motion

%%% Numerical solution
num = (2*pi*lambda_tilde)^2;
den = [1, 2*zeta*(2*pi*lambda_tilde), (2*pi*lambda_tilde)^2];
sys = tf(num, den);

gamma_num = lsim(sys,theta,tau)';


%%% Analytical solution
K = 100;
%calculation coefficients Fourier Analysis
a0 = 2/T * sum(theta(1:36000)*step);
hhh = [1:1:K];
ak = ones(size(hhh));
bk = ones(size(hhh));
ck = ones(size(hhh));
dk = ones(size(hhh));


for i = 1:1:K
    a = 2/T * sum(theta(1:36000).*cos(2*pi*i*tau(1:36000))*step);
    b = 2/T * sum(theta(1:36000).*sin(2*pi*i*tau(1:36000))*step);
    ak(i) = a;
    bk(i) = b;
end 


for i = 1:1:K
    d = (2*zeta*i/lambda_tilde * ak(i) + (1-(i/lambda_tilde)^2)*bk(i))/((2*zeta*i/lambda_tilde)^2 + (1-(i/lambda_tilde)^2)^2);
    c = (-2*zeta*i/lambda_tilde * bk(i) + (1-(i/lambda_tilde)^2)*ak(i))/((2*zeta*i/lambda_tilde)^2 + (1-(i/lambda_tilde)^2)^2);
    ck(i) = c;
    dk(i) = d;
end

g =a0/2;
for i = 1:1:K
    g = g + ck(i)*cos(2*pi*i*tau) + dk(i)*sin(2*pi*i*tau);
end

gamma_anal = g;

figure
plot(tau(24*36000 +1:25*36000), gamma_anal(24*36000+1:25*36000))
hold on
plot(tau(24*36000+1:25*36000),gamma_num(24*36000+1:25*36000))
hold off
title('Dimensionless follower output motion (analytic and numeric)')
xlabel('\tau [-]')
ylabel('\gamma_{anal}(\tau) en \gamma_{num}(\tau)[-]')
legend('Analytical', 'Numeric')


figure
plot(tau(24*36000 +1:25*36000), gamma_anal(24*36000 +1:25*36000) - gamma_num(24*36000 + 1:25*36000))
title('Error analytic and numeric follower output')
xlabel('\tau [-]')
ylabel('\gamma_{anal}(\tau) - \gamma_{num}(\tau) [-]')

%%%deformation
figure
subplot(1,2,1)
plot(tau(24*36000+1:25*36000),(theta(24*36000+1:25*36000) - gamma_num(24*36000+1:25*36000)))
title('Dimensionless deformation of the flexible follower')
xlabel('\tau [-]')
ylabel('\theta(\tau) - \gamma_{num}(\tau)')
subplot(1,2,2)
plot(tau(24*36000+1:25*36000),(25*theta(24*36000+1:25*36000) - 25*gamma_num(24*36000+1:25*36000)))
title('Deformation of the flexible follower in mm')
xlabel('\tau [-]')
ylabel('\theta(\tau) - \gamma_{num}(\tau) [mm]')

figure('Name','Transient response (difference 1st and 25th period)')
plot(tau(1:36000),gamma_num(24*36000+1:25*36000)-gamma_num(1:36000))
title('Transient response (difference 1st and 25th period)')
xlabel('\tau [-]')
ylabel('\gamma(\tau)_{25} - \gamma_(\tau)_1')

extra_normal_force = (theta(24*36000:25*36000-1) - gamma_num(24*36000 +1:25*36000))*25*kf./cos(out.pressure_angle);
total_normal_force = extra_normal_force + out.normalforce_tot;

figure('Name','Normal force')
subplot(1,2,1)
plot(tau(24*36000+1:25*36000),extra_normal_force)
title('Extra normal force due to the flexible follower [N]')
xlabel('\tau [-]')
ylabel('Extra normal force [N]')

subplot(1,2,2)
plot(tau(24*36000+1:25*36000),out.normalforce_tot)
hold on
plot(tau(24*36000+1:25*36000),total_normal_force)
hold off
title('Normal force with and without vibrations[N]')
xlabel('\tau [-]')
ylabel('Normal force [N]')
legend('No vibrations', 'Vibrations')

%%% plots
figure('Name','Single Rise vs. Multi Rise')
plot(tau_sr,10*(theta_sr - gamma_sr))
hold on
plot(tau_sr,25*(theta(4500:9000)-gamma_num(4500:9000)))
hold off

figure('Name','Single Rise vs. Multi Rise')
plot(tau_sr,10*gamma_sr - 25*gamma_num((4500+24*36000):(9000 + 24*36000)))
title('Single rise vs. Multi rise')
xlabel('\tau [-]')
ylabel('10 * \gamma_{sr}(\tau) - 25 *\gamma_{num}(\tau) [mm]')


%%%% Kv and Fvo with vibrations
Fv0_min = 100;
Fv0_step = 1;
Fv0_max = 1200;

Kv_min = 13;
Kv_step = 0.1;
Kv_max = 30;

R = (out.ypitch.^2 + out.xpitch.^2).^(1/2);

Power_min = 10000000;
for Kv = Kv_min:Kv_step:Kv_max
    for Fv0 = Fv0_min:Fv0_step:Fv0_max
        if (total_normal_force + (Fv0 + Kv * out.S) ./ cos(out.pressure_angle)) > 0
            N = total_normal_force + (Fv0 + Kv * out.S)./cos(out.pressure_angle);
            Power = max(abs(out.w.* (N .* sin(out.pressure_angle) .* R)));
            if Power < Power_min
                Power_min = Power;
                Kv_set = Kv;
                Fv0_set = Fv0;
            end 
        
        end
    end
end

Power_min
Kv_set
Fv0_set
