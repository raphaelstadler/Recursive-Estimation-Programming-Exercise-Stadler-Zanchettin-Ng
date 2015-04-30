% Class UKC: UnKnown Constants (only contains constant properties)
%
% Sets the values of some simulation parameters that are unknown
% to the estimator.
%
% Class:
% Recursive Estimation
% Spring 2015
% Programming Exercise 2
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Mark Mueller
% mwm@ethz.ch
%

classdef UKC
    properties (Constant)
        % Number of discrete-time simulation steps:
        simSteps = 200;
        % Random control input segment length range:
        roboVelSegTimeBounds = [5,30]; % discrete-time steps
        % Range that random robot velocities are uniformly distributed on:
        roboVelRange = [1,2]; % (m/s)
        % Probability of a measurement by a sensor in a time step:
        pMeas = 0.2; % (-)
    end
end
