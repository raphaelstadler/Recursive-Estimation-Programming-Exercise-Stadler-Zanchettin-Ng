function const = KnownConstants()
% const = KnownConstants()
% 
% Define the physical constants that are available to the estimator.
%
%
% Class:
% Recursive Estimation
% Spring 2013
% Programming Exercise 1
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Michael Muehlebach
% michaemu@ethz.ch
%
% --
% Revision history
% [15.04.11, ST]    first version by Sebastian Trimpe
% [30.04.12, PR]    adapted version 2012, added unknown wheel radius
% [06.05.13, MH]    adapted version 2013
% [24.04.15, MM]    version 2015



%% Robot kinematic constants

% The nominal right and left wheel radius (W_0), in meter.
const.NominalWheelRadius = 0.1;

% The variability in the wheel radius is captured by a uniform distribution, 
% with the following bound (\bar{\gamma}), in meters.
const.WheelRadiusError = 0.05; % const.WheelRadiusRelativeError = \bar{gamma}

% The wheel base (B), in meters.
const.WheelBase = 0.5;

%% Noise properties

% The compass sensor noise (w_r), normally distributed with zero mean 
% and variance \sigma_r^2, units rad^2.
const.CompassNoise = 0.05; % const.CompassNoise = \sigma_r^2

% The distance sensor noise (w_d), triangular distribution in 
% [-\bar{w}_d,\bar{w}_d], units m.
const.DistNoise = 0.1; % const.DistNoise = \bar{w}_d

% Power spectral density of noise on wheel angular velocity commands (Q_v); 
% multiplicative noise, unit (rad/s)^2/Hz
% This is only relevant for estimator design part 2.
const.VelocityInputPSD = 0.5; % const.VelocityInputNoise = Q_v

% Power spectral density of noise on driving wheel angle commands (Q_r); 
% units rad^2/Hz.
% This is only relevant for estimator design part 2.
const.AngleInputPSD = 0.1; % const.AngleInputNoise = Q_r

%% Starting point

% The robot nominally starts at the origin, uniformly distributed with the
% following bound (\bar{p}), in meters.
const.TranslationStartBound = 1.0; % const.TranslationStartBound = \bar{p}

% The nominal orientation is also 0, and has a uniform distribution with
% the following bound (\bar{r}), in rad.
const.RotationStartBound = pi/2; % const.RotationStartBound = \bar{r}