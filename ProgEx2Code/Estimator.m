function [postParticles] = Estimator(prevPostParticles, sens, act, init)
% [postParticles] = Estimator(prevPostParticles, sens, act, init)
%
% The estimator function. The function will be called in two different
% modes: If init==1, the estimator is initialized. If init == 0, the
% estimator does an iteration for a single sample time interval Ts (KC.ts)
% using the previous posterior particles passed to the estimator in
% prevPostParticles and the sensor measurements and control inputs.
%
% You must edit this function.
%
% Inputs:
%   prevPostParticles   previous posterior particles at discrete time k-1,
%                       which corresponds to continuous time t = (k-1)*Ts
%                       The variable is a struct whose fields are arrays
%                       that correspond to the posterior particle states.
%                       The fields are: (N is number of particles)
%                       .x = (2xN) array with the x-locations (metres)
%                       .y = (2xN) array with the y-locations (metres)
%                       .h = (2xN) array with the headings (radians)
%                       The first row in the arrays corresponds to robot A.
%                       The second row corresponds to robot B.
%
%   sens                Sensor measurements at discrete time k (t = k*Ts),
%                       [4x1]-array, an Inf entry indicates no measurement
%                       of the corresponding sensor.
%                       sens(1): distance reported by sensor 1 (metres)
%                       sens(2): distance reported by sensor 2 (metres)
%                       sens(3): distance reported by sensor 3 (metres)
%                       sens(4): distance reported by sensor 4 (metres)
%
%   act                 Control inputs u at discrete time k-1, which are
%                       constant during a time interval Ts:
%                       u(t) = u(k-1) for (k-1)*Ts <= t < k*Ts
%                       [2x1]-array:
%                       act(1): velocity of robot A, u_A(k-1) (metres/second)
%                       act(2): velocity of robot B, u_B(k-1) (metres/second)
%
%   init                Boolean variable indicating wheter the estimator
%                       should be initialized (init = 1) or if a regular
%                       estimator update should be performed (init = 0).
%                       OPTIONAL ARGUMENT. By default, init = 0.
%
% Outputs:
%   postParticles       Posterior particles at discrete time k, which
%                       corresponds to the continuous time t = k*Ts.
%                       The variable is a struct whose fields are arrays
%                       that correspond to the posterior particle states.
%                       The fields are: (N is number of particles)
%                       .x = (2xN) array with the x-locations (metres)
%                       .y = (2xN) array with the y-locations (metres)
%                       .h = (2xN) array with the headings (radians)
%                       The first row in the arrays corresponds to robot A.
%                       The second row corresponds to robot B.
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

% Check if init argument was passed to estimator:
if(nargin < 4)
    % if not, set to default value:
    init = 0;
end

% Side-length of square room
L = KC.L;
% Constant time interval in which the estimator is called:
dt = KC.ts;

%% Mode 1: Initialization
% Set number of particles:
N = 10; % obviously, you will need more particles than 10.
if (init)
    % Do the initialization of your estimator here!
    % These particles are the posterior particles at discrete time k = 0
    % which will be fed into your estimator again at k = 1
    % Replace the following:
    
    % At time t = 0:
    % - robot A starts out in one of the corners where the sensors S_1 = (L,0) and S_2 = (L,L) are located
    % - robot B starts out in one of the corners where the sensors S_3 = (0,L) and S_4 = (0,0) are located
    
    % A at sensors S_1 and B at sensor S_3
    N_half = floor(N/2);
    postParticles.x(:,1:N_half)  = repmat([L; 0],1,N_half);
    postParticles.y(:,1:N_half)  = repmat([0; L],1,N_half);
    postParticles.h(:,1:N_half)  = repmat([3*pi/4; -pi/4],1,N_half);
    
    % A at sensors S_2 and B at sensor S_4
    postParticles.x(:,(N_half+1):N)  = repmat([L;0],1,N-N_half);
    postParticles.y(:,(N_half+1):N)  = repmat([L;0],1,N-N_half);
    postParticles.h(:,(N_half+1):N)  = repmat([-3*pi/4; pi/4],1,N-N_half); % headings: theta_A, theta_B
    
    % and leave the function
    return;
end % end init

%% Mode 2: Estimator iteration.
% If init = 0, we perform a regular update of the estimator.

%% DYNAMIC OF THE SYSTEM

%% Step 1 (S1): Prior update/Prediction step

% Calculate resulting approximation for f_xp
BinSize = 0.05; % Corresponds to 
NumOfBins = ceil(L/BinSize);
f_xP = zeros(6,NumOfBins);

dX = BinSize*L;    % x and y in 0..L
dH = BinSize*2*pi; % theta in -pi..pi

% Apply the process equation to the particles
for n = 1:N
    % Define synonyms for easier understanding
    xA = prevPostParticles.x(1,n);
    xB = prevPostParticles.x(2,n);
    yA = prevPostParticles.y(1,n);
    yB = prevPostParticles.y(2,n);
    hA = prevPostParticles.h(1,n);
    hB = prevPostParticles.h(2,n);
    uA = act(1);
    uB = act(2);

    % Process Equation
    %
    hA = newHeading(hA,xA,yA,uA); % new hA needs old xA and old yA, so update hA first.
    xA = xA + dt*(uA*cos(hA));
    yA = yA + dt*(uA*sin(hA));
    
    hB = newHeading(hB,xB,yB,uB); % new hB needs old xB and old yB, so update hB first.
    xB = xB + dt*(uB*cos(hB));
    yB = yB + dt*(uB*sin(hB));
    
    for binId = 1:NumOfBins
       if (xA > (binId-1)*dX) && xA <= binId*dX)
         f_xP(1,binId) = f_xP(1,binId) + 1;
       end
    end
    
    % Assign to new variables
    % Theses variables are considered as the
    % x_p (prior variables)
    postParticles.x(1,n) = xA;
    postParticles.x(2,n) = xB;
    postParticles.y(1,n) = yA;
    postParticles.y(2,n) = yB;
    postParticles.h(1,n) = hA;
    postParticles.h(2,n) = hB;
end

%% Step 2 (S2): A posteriori update/Measurement update step

% Noise-free measurement

z_correctRobot = zeros(N,4); % N x 4: For each particle there are 4 measurements

for k = 1:N
% For all particles

    w = drawTriangularRVSample(4);  % parameter defines the vector length of the RV
    z_correctRobot(k,1) = sqrt((xA - L)^2 + yA^2) + w(1);
    z_correctRobot(k,2) = sqrt((xA - L)^2 + yA^2) + w(2);
    z_correctRobot(k,3) = sqrt((xA - L)^2 + yA^2) + w(3);
    z_correctRobot(k,4) = sqrt((xA - L)^2 + yA^2) + w(4);
end

% f_xM(xi) = sum(beta_n*delta(xi-x_p)
z_correctRobot = 0;

% Actual measurement
z_bar = sens;

% Likelihood of measurement given position: f(z[k]=z_bar|x_p[k])
%
% s_bar,    wrong robot is measured:    z_bar[k] consistent with z_correctRobot[k]
% 1-s_bar,  correct robot is measured:  z_bar[k] not consistent with z_correctRobot[k]
likelihoodZBar = zeros(N,1);

for i=1:4
    if z_correctRobot[i] == z_bar[i]
        likelihoodZBar[i] = s_bar;
    else
        likelihoodZBar[i] = 1-s_bar;
    end
end


% Measurement Likelihood beta_n

% Resampling:
% ---------------------------------
% Draw N uniform samples r_n
% Pick particle n_bar, such that:
% sum(beta_n,n,1,n_bar) >= r_n and
% sum(beta_n,n,1,n_bar-1) < r_n


% Sample Impoverishment: Roughening
% ----------------------------------
% Perturb the particles after resampling

function newHeading = newHeading(oldHeading, oldX, oldY, oldU)
    
    alpha = -1;
    newHeading = oldHeading;
    % TODO: Check if orientation changes because of a bouncing and add noise to heading
    
    % Upper wall:
    if (oldY == L) && (oldU*sin(oldHeading) > 0)
        if oldU*cos(oldHeading) > 0
            alpha = oldHeading;
        else
            alpha = pi - oldHeading;
        end
    end
    
    % Right wall:
    if (oldX == L) && (oldU*cos(oldHeading) > 0)
        if oldU*sin(oldHeading) > 0
            alpha = 0.5*pi - oldHeading;
        else
            alpha = -(-0.5*pi - oldHeading);
        end
    end
    
    % Lower wall:
    if (oldY == 0) && (oldU*sin(oldHeading) < 0)
        if oldU*cos(oldHeading) > 0
            alpha = -oldHeading;
        else
            alpha = -(-pi - oldHeading);
        end
    end
    
    % Left wall:
    if (oldX == 0) && (oldU*cos(oldHeading) < 0)
        if oldU*sin(oldHeading) > 0
            alpha = oldHeading - 0.5*pi;
        else
            alpha = -(-pi - oldHeading);
        end
    end
    
end

    function qRV = drawQuadraticRVSample(size)
        u = rand(size);
        c = 3/(2*(KC.vbar^3));
        qRV = zeros(size);
        for ind = 1:size
            if (3/c)*u(ind) - KC.vbar^3 > 0
                qRV(ind) = abs(((3/c)*u(ind) - KC.vbar^3)^(1/3));
            else
                qRV(ind) = -abs(((3/c)*u(ind) - KC.vbar^3)^(1/3));
            end
        end
    end

    function tRV = drawTriangularRVSample(size)
        u = rand(size);
        a = -KC.wbar;
        b = KC.wbar;
        c = 0;
        tRV = zeros(size);
        for ind = 1:size
            if u(ind) < (c-a)/(b-a)
                tRV(ind) = a + sqrt(u(ind)*(b-a)*(c-a));
            else
                tRV(ind) = b - sqrt((1-u(ind))*(b-a)*(b-c));        
            end
        end
    end

end % end estimator

