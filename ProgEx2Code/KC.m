% Class KC: Known Constants (only contains constant properties)
% 
% Sets the values of the constants that are available to the estimator.
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

classdef KC
   properties (Constant)
      % Side-length of square room:
      L = 5; % (meters)
      % Constant time interval in which the estimator is called:
      ts = 0.1; % (seconds)
      % Process Noise Parameter \bar{v} of Quadratic Distribution
      vbar = 1; % (-)
      % Sensor Noise Parameter \bar{w} of Triangular Distribution
      wbar = 1; % (meters)
      % Probability of a sensor reporting distance to wrong robot
      sbar = 0.1; % (-)
   end
end
