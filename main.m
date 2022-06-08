clc
clear
close all

tStart = tic;
addpath('..')
fprintf('%s\n', 'Calculating...');

%% Specify the fiber parameters

% Note: diameter is specified in micrometers, wavelength in nanometers

% Fibre materials (core, cladding)
materials = {'silica'; 'silicaclad'};

% Fibre structure description
fibre = struct(...
    'materials', {materials});

%% Create the task for dispersion curves calculation

% Argument for dispersion curves calculation
argument = struct(...
    'type', 'wvl',... % calculate vs. wavelength
    'harmonic', 1,... % required
    'min', 1000,... % calculate from
    'max', 2000); % calculate to
% Specify which modes to search for

modeTask = struct(...
    'nu', [0 1 2],... % first modal index
    'type', {{'hybrid', 'te', 'tm'}},... % mode types
    'maxmode', 3,... % how many modes of each type and NU to calculate
    'diameter', 9);%,... % parameter, structure diameter, if argument is wavelength
    %'region', 'cladding');

%% Find guided modes and calculate the dispersion curves

infomode = false;
modes = buildModes(argument, fibre, modeTask, infomode);

% Display calculated dispersion curves
showModes(modes, 'n_{eff} in a fiber');

% Show which modes have been calculated
fprintf('Modes found:\n');
for i = 1:numel(modes)
    fprintf('%s\n', modeDescription(modes(i), infomode));
end

%% Field distribution

poiType = 'wvl';
poi = 1550;
par = 9;
task = struct('modetype', 'HYBRID', 'modeindex', [1 1]);
infomode = false;
getWholeMode = false;
n = neff(poiType, poi, par, fibre, task, infomode, getWholeMode);

%
d = par;
window = 5 * d;
lambda = poi;
Nr = [1000, 2000]; % inn, out
Nphi = [64, 64]; % inn, out
dr_inn = d/2 / Nr(1);
dr_out = (window/2 - d/2) / Nr(2);

dphi_inn = 2 * pi / Nphi(1);
dphi_out = 2 * pi / Nphi(2);

F = struct('dr', [dr_inn dr_out], 'dphi', [dphi_inn dphi_out], ...
    'diam', [d window], 'E1', [], 'H1', [], 'E2', [], 'H2', []);
FG = fieldGrid(F);
[F.E1, F.H1] = modeField(lambda, d, n, fibre, task, FG.R1, FG.PHI1);
[F.E2, F.H2] = modeField(lambda, d, n, fibre, task, FG.R2, FG.PHI2);
F.FG = FG;

displayField2(F, d)

%% Mode dispersion in a nanofibre

% Calculate mode dispersion diagram for a fixed wavelength and varying fibre diameter
% Specify the fibre parameters

materials = {'silica'; 'silicaclad'};

nanofibre = struct(...
    'materials', {materials});

argument = struct(...
    'type', 'dia',... % calculate vs. wavelength
    'min', 5,... % calculate from
    'max', 15); % calculate to

modeTask = struct(...
    'nu', [0 1 2],... % first modal index
    'type', {{'hybrid', 'te', 'tm'}},... % mode types
    'maxmode', 3,... % how many modes of each type and NU to calculate
    'lambda', 1550);%,... % parameter, structure diameter, if argument is wavelength
    %'region', 'cladding');

modes = buildModes(argument, nanofibre, modeTask, false);

%% Display calculated dispersion curves

figure;
showModes(modes, 'Modal dispersion in a nanofibre');

%% Show which modes have been calculated

fprintf('Modes found:\n');
for i = 1:numel(modes)
    fprintf('%s\n', modeDescription(modes(i), false));
end

% Unlike in the weakly guidance example above, this nanofibre has a high
% refractive index step (about 0.5) and is therefore strongly guiding.
% All four modes HE11, HE21, TE01 and TM01 are therefore clearly distinct.
% The LP approximation is not valid in this case any more. The Optical
% Fibre Toolbox correctly simulates this situation using the full vector
% solution of the Maxwell equations.

toc(tStart)

%% Nonlinear Schrodinger Equation with Split-Step Method2

%idu/dz-sgn(beta2)/2 d^2u/d(tau)^2 + N^2*¦u¦^2*u = 0

clear all;

%Specify input parameters
distance = 10;%input('Enter fiber length (in units of L_D) ='); %
beta2 = -1;%input('dispersion: 1 for normal -1 for anomalous');%
N = 10;%input('Nonlinear parameter N ='); % Soliton order
mshape = 0;%input('m=0 for sech, m>0 for super-Gaussian =');
chirp0 = 0; % input pulse chirp (default value)

%Simulation Parameters
nt=2048; Tmax=32; %FFT ponts and window size
step_num = round(20*distance*N^2); %No. of z steps to
deltaz = distance/step_num; %Step size in z
dtau = 2*Tmax/nt; % Step size in tau

%tau and omega arrays
tau = (-nt/2:nt/2-1)*dtau; % temporal grid
omega = (pi/Tmax)*[(0:nt/2-1) (-nt/2:-1)]; %frequency grid

if mshape==0
    uu = sech(tau).*exp(-0.5i*chirp0*tau.^2);
else
    uu = exp(-0.5*(1+1i*chirp0).*tau.^2);
end

temp= fftshift(ifft(uu)).*(nt*dtau)/sqrt(2*pi); %spectrum
figure; subplot(2,1,1);
plot(tau, abs(uu).^2,'--k'); hold on;
axis([-5 5 0 inf]);
xlabel('Normalized Time');
ylabel('Normalized Power');
subplot(2,1,2);
plot(fftshift(omega)/(2*pi), abs(temp).^2, '--k'); hold on;
axis([-5 5 0 inf]);
xlabel('Normalized Freq');
ylabel('Spectral Power');

%---storde dispersive shifts to speed up the code
dispersion = exp(i*0.5*beta2*omega.^2*deltaz); %phase factor
hhz = 1i*N^2*deltaz; %nonlonear phase factor

% Main Loop
%Scheme: 1/2N-> D -> 1/2N first half step nonlinear
temp = uu.*exp(abs(uu).^2.*hhz/2);

for n=1:step_num
    f_temp=ifft(temp).*dispersion;
    uu=fft(f_temp);
    temp = uu.*exp(abs(uu).^2.*hhz);
end

uu = temp.*exp(-abs(uu).^2.*hhz/2); %final field
temp = fftshift(ifft(uu)).*(nt*dtau)/sqrt(2*pi); %final spectrum
% End of main loop

% Plot output pulse shape and spectrum
hold on
subplot(2,1,1)
plot(tau,abs(uu).^2,'--r');
subplot(2,1,2)
plot(fftshift(omega)/(2*pi),abs(temp).^2,'--r');
