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
