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
    % Do the initialization of your estimator here!
    %uniform distribution between [-b1,b1]
    %posEst = [unifrnd(-b1,b1) , unifrnd(-b1,b1)];
    b1 = knownConst.TranslationStartBound;
    posEst = [0,0];
    posVar = [(b1^2)/3 , (b1^2)/3];
    %uniform distribution between [-b2,b2]
    %oriEst = unifrnd(0,b2) + unifrnd(0,b2) - b2;
    b2 = knownConst.RotationStartBound;
    oriEst = 0;
    oriVar = (b2^2)/3;
    %radius W0 plups uniform distribution between [-b3,b3]
    %radiusEst = knownConst.NominalWheelRadius + unifrnd(-b3,b3);
    b3 = knownConst.WheelRadiusError;
    radiusEst = knownConst.NominalWheelRadius;
    radiusVar = (b3^2)/3;
    
    %%TODO ??
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
% x = [x,y,r,W]
% TODO v = [gamma];

B = knownConst.WheelBase;

q = @(t,x) [ x(4)*actuate(1)*cos(actuate(2))*cos(x(3));
             x(4)*actuate(1)*cos(actuate(2))*sin(x(3));
                    -x(4)*actuate(1)*sin(actuate(2))/B;
                                                     0];
                                                 
A = @(x) [ 0, 0, -x(4)*actuate(1)*cos(actuate(2))*sin(x(3)), actuate(1)*cos(actuate(2))*cos(x(3));
           0, 0,  x(4)*actuate(1)*cos(actuate(2))*cos(x(3)), actuate(1)*cos(actuate(2))*sin(x(3));           
           0, 0,                                          0,        -actuate(1)*sin(actuate(2))/B;
           0, 0,                                          0,                                    0];

% TODO it is right that gamma is a noise??
L = @(x) [actuate(1)*cos(actuate(2))*cos(x(3)); 
          actuate(1)*cos(actuate(2))*sin(x(3)); 
                 -actuate(1)*sin(actuate(2))/B; 
                                             0];
                                         
Q = 0;%(knownConst.WheelRadiusError^2)/3;


%TODO This is my idea to resolve between (k-1)T and tm=kT????
tspan = [estState.Time,tm];

xm = estState.Est;
[~,Y] = ode45(q,tspan,xm);
xp = Y(end,:)';

Pm = estState.Var;
[~,Y] = ode45(@(t,X)mRiccati(t, X, A(xp), L(xp), Q), tspan, Pm(:));
Pp = Y(end,:)';
Pp = reshape(Pp, size(A(xp)));


if(not(sense(1) == Inf))
    if(not(sense(2) == Inf))
        h = @(x) [ x(3);
           sqrt(x(1)^2+x(2)^2)];
        H = @(x) [                          0,                          0, 1, 0;
           x(1)/(sqrt(x(1)^2+x(2)^2)), x(2)/(sqrt(x(1)^2+x(2)^2)), 0, 0];
        M = eye(2);
        R = [ knownConst.CompassNoise,                          0; 
                            0, (knownConst.DistNoise^2)/6];
        K = Pp*H(xp).'/(H(xp)*Pp*H(xp).' + M*R*M.');

        xm = xp + K*(sense' - h(xp));
        Pm = (eye(4) - K*H(xp))*Pp;
    else
        h = @(x) x(3);
        H = @(x) [0, 0, 1, 0];
        M = 1;
        R = knownConst.CompassNoise;
        K = Pp*H(xp).'/(H(xp)*Pp*H(xp).' + M*R*M.');
        xm = xp + K*(sense(1) - h(xp));
        Pm = (1 - K*H(xp))*Pp;  
    end
else
    if(not(sense == Inf))
        h = @(x) sqrt(x(1)^2+x(2)^2);
        H = @(x) [x(1)/(sqrt(x(1)^2+x(2)^2)), x(2)/(sqrt(x(1)^2+x(2)^2)), 0, 0];
        M = 1;
        R = (knownConst.DistNoise^2)/6;
        K = Pp*H(xp).'/(H(xp)*Pp*H(xp).' + M*R*M.');
        xm = xp + K*(sense(2) - h(xp));
        Pm = (1 - K*H(xp))*Pp;
    else
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
    %Convert from "n^2"-by-1 to "n"-by-"n"
    P = reshape(P, size(A));
    %Determine derivative
    dP = A*P + P*A.' + L*Q*L.'; 
    %Convert from "n"-by-"n" to "n^2"-by-1
    dP = dP(:); 
end