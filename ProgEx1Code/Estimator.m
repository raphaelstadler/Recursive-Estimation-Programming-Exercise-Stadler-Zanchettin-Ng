function [posEst,oriEst,radiusEst, posVar,oriVar,radiusVar,estState] = Estimator(estState,actuate,sense,tm,knownConst,designPart)
% [posEst,oriEst,posVar,oriVar,baseEst,baseVar,estState] =
% 	Estimator(estState,actuate,sense,tm,knownConst,designPart)
%
% The estimator.
%
% The Estimator function shall be used for both estimator design parts; the
% input argument designPart is used to distinguish the two:
%   designPart==1  -> Part 1
%   designPart==2  -> Part 2
%
% The function will be called in two different modes:
% If tm==0, the estimator is initialized; otherwise the estimator does an
% iteration step (compute estimates for the time step k).
%
% Inputs:
%   estState        previous estimator state (time step k-1)
%                   May be defined by the user (for example as a struct).
%   actuate         control input u(k-1), [1x2]-vector
%                   actuate(1): u_v(k-1), drive wheel angular velocity
%                   actuate(2): u_r(k-1), drive wheel angle
%   sense           sensor measurements z(k), [1x2]-vector, INF if no
%                   measurement
%                   sense(1): z_d(k), distance measurement
%                   sense(2): z_r(k), orientation measurement
%   tm              time, scalar
%                   If tm==0 initialization, otherwise estimator
%                   iteration step.
%   knownConst      known constants (from KnownConstants.m)
%   designPart      variable to distinguish the estimator design part
%                       designPart==1  -> Part 1
%                       designPart==2  -> Part 2
%
% Outputs:
%   posEst          position estimate (time step k), [1x2]-vector
%                   posEst(1): x position estimate
%                   posEst(2): y position estimate
%   oriEst          orientation estimate (time step k), scalar
%   radiusEst       estimate of wheel radius W (time step k), scalar
%   posVar          variance of position estimate (time step k), [1x2]-vector
%                   posVar(1): x position variance
%                   posVar(2): y position variance
%   oriVar          variance of orientation estimate (time step k), scalar
%   radiusVar       variance of wheel radius estimate (time step k), scalar
%   estState        current estimator state (time step k)
%                   Will be input to this function at the next call.
%
%
% Class:
% Recursive Estimation
% Spring 2015
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
% [19.04.11, ST]    first version by Sebastian Trimpe
% [30.04.12, PR]    adapted version for spring 2012, added unknown wheel
%                   radius
% [06.05.13, MH]    2013 version
% [24.04.15, MM]    2015 version


%% Mode 1: Initialization
if (tm == 0)
    % posEst / posVar
    % =====================
    % Initial position estimation posEst = [x(0), y(0)] = [x0, y0]
    % (x0, y0) uniformly distributed between [-p_bar,p_bar]
    p_bar = knownConst.TranslationStartBound;
    x0 = 0;     %   E[ unifrnd(-p_bar,p_bar) ]
    y0 = 0;     %   E[ unifrnd(-p_bar,p_bar) ]
    posEst = [x0, y0];
    
    posVar = [(p_bar^2)/3, (p_bar^2)/3];    % Var(unif(-p_bar,p_bar)) can be calculated like that
    
    % oriEst
    % =====================
    r_bar = knownConst.RotationStartBound;
    oriEst = 0; %   E[ unifrnd(-r_bar, r_bar) ]
    
    oriVar = (r_bar^2)/3;
    
    % radiusEst / radiusVar
    % =====================
    gamma = knownConst.WheelRadiusError;
    radiusEst = knownConst.NominalWheelRadius;  % E[ W0 + unifrnd(-gamma,gamma) ]
    
    radiusVar = (gamma^2)/3;
    
    field1 = 'Est';
    value1 = {[posEst';oriEst;radiusEst]};
    field2 = 'Var';
    value2 = {[posVar(1),         0,      0,         0;
                       0, posVar(2),      0,         0;
                       0,         0, oriVar,         0;
                       0,         0,      0, radiusVar] };
    field3 = 'Time';
    value3 = {tm};    
    estState = struct(field1,value1,field2,value2,field3,value3);
    return;
end


%% Mode 2: Estimator iteration.
% If we get this far tm is not equal to zero, and we are no longer
% initializing.  Run the estimator.

B = knownConst.WheelBase;

% TODO v = [gamma];
% X = [x,y,r,W]
% q(t,x) = X_dot    (Note that, actuate u is constant during time period, noise v set to zero)
q = @(t,x) [ x(4)*actuate(1)*cos(actuate(2))*cos(x(3)); % x_dot = s_v*cos(u_r)*cos(r) = W*u_v*cos(u_r)*cos(r)
             x(4)*actuate(1)*cos(actuate(2))*sin(x(3)); % y_dot = s_t*sin(r) = W*u_v*cos(u_r)*sin(r)
                    -x(4)*actuate(1)*sin(actuate(2))/B; % r_dot = s_r = -s_v*sin(u_r)/B = -W*u_v*sin(u_r)/B
                                                     0];% W_dot = 0
                                
% A = partial_der(q,X)
A = @(x) [ 0, 0, -x(4)*actuate(1)*cos(actuate(2))*sin(x(3)), actuate(1)*cos(actuate(2))*cos(x(3));  % partial_der(x,X)
           0, 0,  x(4)*actuate(1)*cos(actuate(2))*cos(x(3)), actuate(1)*cos(actuate(2))*sin(x(3));  % partial_der(y,X)
           0, 0,                                          0,        -actuate(1)*sin(actuate(2))/B;  % partial_der(r,X)
           0, 0,                                          0,                                    0]; % partial_der(W,X)

if designPart==1
    % The motion of the robot is corrupted by process noise, the properties
    % of which are unknown.
    % .: Do not take any process noise into account in the model
    %L = @(x) zeros(4);
    %Q = zeros(4);
    L = @(x) eye(4);
    Q = zeros(4,4);
    Q(1,1) = .1;
    Q(2,2) = .1;
    Q(3,3) = .01;
    
elseif designPart==2
    % A model of the process noise is available. This model takes
    % non-idealities in the actuation mechanism into account
    
    % L = partial_der(q,v) | v = 0, where v = [v_v, v_r]'
      L = @(x) x(4)*actuate(1)*...
          [ cos(actuate(2))*cos(x(3)),  -sin(actuate(2))*cos(x(3)); % partial_der(q(1),v)
            cos(actuate(2))*sin(x(3)),  -sin(actuate(2))*sin(x(3)); % partial_der(q(2),v)
            -sin(actuate(2))/B,         -cos(actuate(2))/B;         % partial_der(q(3),v)
            0,                          0                        ]; % partial_der(q(4),v)
    
    % Additionally to the information of the EKF in design part 1 above,
    % this EKF has acces to the constants Q_v and Q_r
    Q_v = knownConst.VelocityInputPSD;
    Q_r = knownConst.AngleInputPSD;
    Q = [Q_v,   0;
           0,     Q_r];
else
    error('Invalid designPart chosen. Choose designPart==1 or designPart==2.');
end


%TODO This is my idea to resolve between (k-1)T and tm=kT????
tspan = [estState.Time,tm];

%% Step 1 (S1): Prior update/Prediction step
%
xm = estState.Est;
[~,Yx] = ode45(q,tspan,xm'); % Solve x_dot = q(t,x) for t in tspan
xp = Yx(end,:)';

Pm = estState.Var;
% Solve Riccati Matrix Equation:
% P_dot = A*P + P*A' + L*Q*L'
[~,Yp] = ode45(@(t,X)mRiccati(t, X, A, L, Q,Yx,tspan), tspan, Pm(:));
Pp = Yp(end,:)';
Pp = reshape(Pp, size(A(xp)));

h = @(x) [ sqrt(x(1)^2+x(2)^2);
                         x(3)];

%% Step 2 (S2): A posteriori update/Measurement update step
if(not(sense(1) == Inf))
    % Measurement for sensor 1 is available
    if(not(sense(2) == Inf))
        % Measurement for sensor 2 is available
        % Both sensors provide useable measurements.
        % h_k = h(x(kT),w=0) = z(kT)
        H = @(x) [ x(1)/(sqrt(x(1)^2+x(2)^2)), x(2)/(sqrt(x(1)^2+x(2)^2)), 0, 0
                                            0,                          0, 1, 0];
    else
        % Measurement for sensor 2 is not available
        % Only sensor 1 provides useable measurement
        % H_k = partial_der(h_k, x)
        H = @(x) [x(1)/(sqrt(x(1)^2+x(2)^2)), x(2)/(sqrt(x(1)^2+x(2)^2)), 0, 0;
                                           0,                          0, 0, 0];
        hh = h(xp);
        sense(2) = hh(2);
    end
else
    % Measurement for sensor 1 is not available
    if(not(sense == Inf))
        % Measurement for sensor 2 is available
        % Only sensor 2 provides useable measurement
        % H_k = partial_der(h_k, x)
        H = @(x) [0, 0, 0, 0;
                  0, 0, 1, 0];
        hh = h(xp);
        sense(1) = hh(1);
    else
    % Measurement for sensor 2 is not available
    % Both sensors 1 and 2 do not provide useable measurements
    % No new information available!
    H = @(x) [0, 0, 0, 0;
              0, 0, 0, 0];
    sense = h(xp)';
    end
end
       
% M_k = partial_der(h_k, w)
M = eye(2);

% R is the co-variance matrix of the noise w
% Calculated from triangular pdf of CRV d
sigma_d_sq = (knownConst.DistNoise^2)/6;
sigma_r_sq = knownConst.CompassNoise;
R = [ sigma_d_sq,          0;
               0, sigma_r_sq];
           
% K: Kalman Gain matrix
K = Pp*H(xp).'/(H(xp)*Pp*H(xp).' + M*R*M.');

% Measurement update
xm = xp + K*(sense' - h(xp));
Pm = (eye(4) - K*H(xp))*Pp;
        
posEst = [xm(1) xm(2)];
posVar = [Pm(1,1),Pm(2,2)];
oriEst = xm(3);
oriVar = Pm(3,3);
radiusEst = xm(4);
radiusVar = Pm(4,4);

%%TODO ??
field1 = 'Est';
value1 = {xm};
field2 = 'Var';
value2 = {Pm};
field3 = 'Time';
value3 = {tm};
estState = struct(field1,value1,field2,value2,field3,value3);
end

function dP = mRiccati(t, P, A, L, Q, Y, tspan)
    detat = (tspan(2)-tspan(1))/(size(Y,1)-1); 
    xx = interp1(tspan(1):detat:tspan(2),Y,t);
    A = A(xx');
    L = L(xx');
    P = reshape(P, size(A)); %Convert from "n^2"-by-1 to "n"-by-"n"
    dP = A*P + P*A.' + L*Q*L.'; %Determine derivative
    dP = dP(:); %Convert from "n"-by-"n" to "n^2"-by-1
end