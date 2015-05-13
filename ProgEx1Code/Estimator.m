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
%                   uv: u_v(k-1), drive wheel angular velocity
%                   ur: u_r(k-1), drive wheel angle
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
    % posEst
    % =====================
    p_bar = knownConst.TranslationStartBound;
    x0 = 0;     %   E[ unifrnd(-p_bar,p_bar) ]
    y0 = 0;     %   E[ unifrnd(-p_bar,p_bar) ]
    posEst = [x0, y0];
    posVar = [(p_bar^2)/3, (p_bar^2)/3];  
    
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

    % SAVE initial state
    field1 = 'X';
    value1 = {[posEst'; oriEst; radiusEst]};
    field2 = 'P';
    value2 = {diag([posVar(1), posVar(2), oriVar, radiusVar]) };
    field3 = 'T';
    value3 = {tm};
    estState = struct(field1,value1,field2,value2,field3,value3);
    return;
end


%% Mode 2: Estimator iteration.
% If we get this far tm is not equal to zero, and we are no longer
% initializing.  Run the estimator.
B = knownConst.WheelBase;
uv = actuate(1);
ur = actuate(2);

%% DYNAMIC OF THE SYSTEM
% X = [x,y,r,W]
% q(t,x) = X_dot    (Note that, actuate u is constant during time period, noise v set to zero)
q = @(t,x) x(4)*uv*...
           [ cos(ur)*cos(x(3));  % x_dot = s_v*cos(u_r)*cos(r) = W*u_v*cos(u_r)*cos(r)
             cos(ur)*sin(x(3));  % y_dot = s_t*sin(r) = W*u_v*cos(u_r)*sin(r)
                    -sin(ur)/B;  % r_dot = s_r = -s_v*sin(u_r)/B = -W*u_v*sin(u_r)/B
                             0]; % W_dot = 0
                     
%% A = partial_der(q,X)
A = @(x) uv*...
         [ 0, 0, -x(4)*cos(ur)*sin(x(3)), cos(ur)*cos(x(3));  % partial_der(x,X)
           0, 0,  x(4)*cos(ur)*cos(x(3)), cos(ur)*sin(x(3));  % partial_der(y,X)
           0, 0,                       0,        -sin(ur)/B;  % partial_der(r,X)
           0, 0,                       0,                 0]; % partial_der(W,X)
    
%% SYSTEM NOISE
if (designPart==1)
    % The motion of the robot is corrupted by process noise, the properties
    % of which are unknown.
    % L = partial_der(q,v) | v = 0, where v = [v_1, v_2, v_3, v_4]'
    L = @(x) eye(4);
    % Founded with genetic algorithm to minimize final error
    Q = diag([0.07, 0.07, 0.01, 0]);
    
elseif (designPart==2)
    % A model of the process noise is available. This model takes
    % non-idealities in the actuation mechanism into account
    % L = partial_der(q,v) | v = 0, where v = [v_v, v_r]'
    L = @(x) x(4)*uv*...
          [ cos(ur)*cos(x(3)), -sin(ur)*cos(x(3));  % partial_der(q(1),v)
            cos(ur)*sin(x(3)), -sin(ur)*sin(x(3));  % partial_der(q(2),v)
                   -sin(ur)/B,         -cos(ur)/B;  % partial_der(q(3),v)
                            0,                  0]; % partial_der(q(4),v)
    
    % Additionally to the information of the EKF in design part 1 above,
    % this EKF has acces to the constants Q_v and Q_r
    Q_v = knownConst.VelocityInputPSD;
    Q_r = knownConst.AngleInputPSD;
    Q = diag([Q_v, Q_r]);
else
     error('Invalid designPart chosen. Choose designPart==1 or designPart==2.');
end

%% Resolve the SYSTEM DYNAMIC between [(k-1)T,kT]
tspan = [estState.T,tm];

%% Step 1 (S1): Prior update/Prediction step
% Solve the equations of motion
x0 = estState.X;
[Xx,Yx] = ode45(q,tspan,x0'); % Solve x_dot = q(t,x) for t in tspan
xp = Yx(end,:)';

P0 = estState.P;
% Solve Riccati Matrix Equation:
% P_dot = A*P + P*A' + L*Q*L'
[~,Yp] = ode45(@(t,X)mRiccati(t, X, A, L, Q, Yx, Xx), tspan, P0(:));
Pp = reshape(Yp(end,:)',size(A(xp)));


%% DYNAMIC OF THE MEASUREMENT;
% h_k = h(x(kT),w=0) = z(kT)
h = [ sqrt(xp(1)^2+xp(2)^2);
                     xp(3)];

if (h(1)== 0)
  sense(1) = inf;  
end

%% Step 2 (S2): A posteriori update/Measurement update step
if (not(sense(1) == inf))
    % Measurement for sensor 1 is available
    if (not(sense(2) == inf))
        % Measurement for sensor 2 is available
        % Both sensors provide useable measurements.
        % H_k = partial_der(h_k, x)
        H = [ xp(1)/h(1), xp(2)/h(1), 0, 0;
                       0,          0, 1, 0];
    else
        % Measurement for sensor 2 is not available
        % Only sensor 1 provides useable measurement
        % H_k = partial_der(h_k, x)
        H = [ xp(1)/h(1), xp(2)/h(1), 0, 0;
                       0,          0, 0, 0];
        sense(2) = h(2);
    end
else
    % Measurement for sensor 1 is not available
    if (not(sense(2) == inf))
        % Measurement for sensor 2 is available
        % Only sensor 2 provides useable measurement
        % H_k = partial_der(h_k, x)
        H = [0, 0, 0, 0;
             0, 0, 1, 0];
        sense(1) = h(1);
    else
    % Measurement for sensor 2 is not available
    % Both sensors 1 and 2 do not provide useable measurements
    % No new information available!
    H = [0, 0, 0, 0;
         0, 0, 0, 0];
    sense = h';
    end
end
       
% M_k = partial_der(h_k, w)
M = eye(2);

% R is the co-variance matrix of the noise w
% Calculated from triangular pdf of CRV d
sigma_d_sq = (knownConst.DistNoise^2)/6;
sigma_r_sq = knownConst.CompassNoise;
R = diag([ sigma_d_sq , sigma_r_sq]);
           
%% K: Kalman Gain matrix 
K = Pp*H'/(H*Pp*H' + M*R*M');

%% Measurement update
xm = xp + K*(sense' - h);
Pm = (eye(4) - K*H)*Pp;

%% SAVE new state
posEst = [xm(1), xm(2)];
posVar = [Pm(1,1), Pm(2,2)];
oriEst = xm(3);
oriVar = Pm(3,3);
radiusEst = xm(4);
radiusVar = Pm(4,4);

field1 = 'X';
value1 = {xm};
field2 = 'P';
value2 = {Pm};
field3 = 'T';
value3 = {tm};
estState = struct(field1,value1,field2,value2,field3,value3);

end

function dP = mRiccati(t, P, A, L, Q, Y, X)
    %A is time dependent, has to be calculate along 
    %the trajectory of xp(t)
    xx = interp1(X,Y,t);
    A = A(xx');
    L = L(xx');

    P = reshape(P,size(A)); %Convert from "n^2"-by-1 to "n"-by-"n"
    dP = A*P + P*A' + L*Q*L'; %Determine derivative
    dP = dP(:); %Convert from "n"-by-"n" to "n^2"-by-1
end