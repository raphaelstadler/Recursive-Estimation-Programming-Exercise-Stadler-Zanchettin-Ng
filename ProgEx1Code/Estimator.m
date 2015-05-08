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
q = @(t,x) [ x(4)*actuate(1)*cos(actuate(2))*cos(x(3)); % x_dot = s_v*cos(u_r*cos(r) = W*u_v*cos(u_r)*cos(r)
             x(4)*actuate(1)*cos(actuate(2))*sin(x(3)); % y_dot = s_t*sin(r) = W*u_v*cos(u_r)*cos(r)
                    -x(4)*actuate(1)*sin(actuate(2))/B; % r_dot = s_r = -s_v*sin(u_r)/B = -W*u_v*sin(u_r)/B
                                                     0];% W_dot = 0
                                
% A = partial_der(q,X)
A = @(x) [ 0, 0, -x(4)*actuate(1)*cos(actuate(2))*sin(x(3)), actuate(1)*cos(actuate(2))*cos(x(3));  % partial_der(x,X)
           0, 0,  x(4)*actuate(1)*cos(actuate(2))*cos(x(3)), actuate(1)*cos(actuate(2))*sin(x(3));  % partial_der(y,X)
           0, 0,                                          0,        -actuate(1)*sin(actuate(2))/B;  % partial_der(r,X)
           0, 0,                                          0,                                    0]; % partial_der(W,X)

% TODO it is right that gamma is a noise??
% TODO: Probably not, it should be considered as a constant bias
% L = partial_der(q,v)
% L = @(x) [actuate(1)*cos(actuate(2))*cos(x(3)); 
%           actuate(1)*cos(actuate(2))*sin(x(3)); 
%                  -actuate(1)*sin(actuate(2))/B; 
%                                              0];
L = 0;

% Q = variance gamma                                         
% Q = (knownConst.WheelRadiusError^2)/3;
% Q = E[ v(t) ]
Q = 0;

%TODO This is my idea to resolve between (k-1)T and tm=kT????
tspan = [estState.Time,tm];

%% Step 1 (S1): Prior update/Prediction step
%
xm = estState.Est;
[~,Yx] = ode45(q,tspan,xm); % Solve x_dot = q(t,x) for t in tspan
xp = Yx(end,:)';

Pm = estState.Var;
% Solve Riccati Matrix Equation:
% P_dot = A*P + P*A' + L*Q*L'
[~,Yp] = ode45(@(t,X)mRiccati(t, X, A(xp), L, Q), tspan, Pm(:));
Pp = Yp(end,:)';
Pp = reshape(Pp, size(A(xp)));


%% Step 2 (S2): A posteriori update/Measurement update step
%
if(not(sense(1) == Inf))
% Measurement for sensor 1 is available
    if(not(sense(2) == Inf))
    % Measurement for sensor 2 available:
    % Both sensors provide useable measurements.
        % h_k = h(x(kT),w(kT)) = z(kT)
        h = @(x) [ x(3);
           sqrt(x(1)^2+x(2)^2)];
       
        % H_k = partial_der(h_k, x)
        H = @(x) [  0,                          0,                          1, 0;
                    x(1)/(sqrt(x(1)^2+x(2)^2)), x(2)/(sqrt(x(1)^2+x(2)^2)), 0, 0];
        
        % M_k = partial_der(h_k, w)
        M = eye(2);
        
        % R is the co-variance matrix of the noise w
        sigma_r_sq = knownConst.CompassNoise;
        sigma_d_sq = (knownConst.DistNoise^2)/6;    % Calculated from triangular pdf of CRV d
        R = [   sigma_r_sq, 0;
                0,          sigma_d_sq];
        
        % K: Kalman Gain matrix
        K = Pp*H(xp).'/(H(xp)*Pp*H(xp).' + M*R*M.');

        % Measurement update
        xm = xp + K*(sense' - h(xp));
        Pm = (eye(4) - K*H(xp))*Pp;
    else
    % NO measurement for sensor 2:
    % Only sensor 1 provides useable measurement
        % h_k = h(x(kT),w(kT)) = z(kT)
        h = @(x) x(3);
        
        % H_k = partial_der(h_k, x)
        H = @(x) [0, 0, 1, 0];
        
        % M_k = partial_der(h_k, w)
        M = 1;
        
        % R is the co-variance matrix of the noise w
        R = knownConst.CompassNoise;
        
        % K: Kalman Gain matrix
        K = Pp*H(xp).'/(H(xp)*Pp*H(xp).' + M*R*M.');
        
        % Measurement update
        xm = xp + K*(sense(1) - h(xp));
        Pm = (1 - K*H(xp))*Pp;  
    end
else
% NO measurement for sensor 1
    if(not(sense(2) == Inf))
    % Measurement for sensor 2 available:
    % Only sensor 2 provides useable measurement
        % h_k = h(x(kT),w(kT)) = z(kT)
        h = @(x) sqrt(x(1)^2+x(2)^2);
        
        % H_k = partial_der(h_k, x)
        H = @(x) [x(1)/(sqrt(x(1)^2+x(2)^2)), x(2)/(sqrt(x(1)^2+x(2)^2)), 0, 0];
        
        % M_k = partial_der(h_k, w)
        M = 1;
        
        % R is the co-variance matrix of the noise w
        R = (knownConst.DistNoise^2)/6;
        
        % K: Kalman Gain matrix
        K = Pp*H(xp).'/(H(xp)*Pp*H(xp).' + M*R*M.');
       
        % Measurement update        
        xm = xp + K*(sense(2) - h(xp));
        Pm = (1 - K*H(xp))*Pp;
    else
    % NO measurement for sensor 2 available:
    % Both sensors 1 and 2 do not provide useable measurements
        % Measurement update
        % No new information available!
        xm = xp;
        Pm = Pp;
    end 
end

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

function dP = mRiccati(t, P, A, L, Q)
    P = reshape(P, size(A)); %Convert from "n^2"-by-1 to "n"-by-"n"
    dP = A*P + P*A.' + L*Q*L.'; %Determine derivative
    dP = dP(:); %Convert from "n"-by-"n" to "n^2"-by-1
end