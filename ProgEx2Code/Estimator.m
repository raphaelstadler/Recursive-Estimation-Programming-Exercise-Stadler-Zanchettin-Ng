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
    % TODO: Draw a uniform RV to get initial headings
    postParticles.h(:,1:N_half)  = repmat([3*pi/4; -pi/4],1,N_half);
    
    % A at sensors S_2 and B at sensor S_4
    postParticles.x(:,(N_half+1):N)  = repmat([L;0],1,N-N_half);
    postParticles.y(:,(N_half+1):N)  = repmat([L;0],1,N-N_half);
    % TODO: Draw a uniform RV to get initial headings
    postParticles.h(:,(N_half+1):N)  = repmat([-3*pi/4; pi/4],1,N-N_half); % headings: theta_A, theta_B
    
    % and leave the function
    return;
end % end init

%% Mode 2: Estimator iteration.
% If init = 0, we perform a regular update of the estimator.

%% Step 1 (S1): Prior update/Prediction step

% Define synonyms for easier understanding
xA = prevPostParticles.x(1,:);
yA = prevPostParticles.y(1,:);
hA = prevPostParticles.h(1,:);

xB = prevPostParticles.x(2,:);
yB = prevPostParticles.y(2,:);
hB = prevPostParticles.h(2,:);

uA = act(1); uB = act(2);

vA = drawQuadraticRVSample([N,1]); % draw noise from quadratic pdf for post-bounce angles
vB = drawQuadraticRVSample([N,1]);

% Initialize variables which will be assigned after prior update
xA_P = zeros(1,N); yA_P = zeros(1,N); hA_P = zeros(1,N);
xB_P = zeros(1,N); yB_P = zeros(1,N); hB_P = zeros(1,N);

% TODO: Optimize newHeading(.) function to be able to treat a vector as a
% whole

% Apply the process equation to the particles
for n = 1:N
    % DYNAMIC OF THE SYSTEM: Process Equation
    %
    % Propagate N particles x_m through process dynamics, to get new
    % particles x_p
    hA_P(n) = newHeading(hA(n),xA(n),yA(n),uA,vA(n)); % new hA needs old xA and old yA, so update hA first
    xA_P(n) = xA(n) + dt*(uA*cos(hA(n)));
    yA_P(n) = yA(n) + dt*(uA*sin(hA(n)));
    
    hB_P(n) = newHeading(hB(n),xB(n),yB(n),uB,vB(n)); % new hB needs old xB and old yB, so update hB first
    xB_P(n) = xB(n) + dt*(uB*cos(hB(n)));
    yB_P(n) = yB(n) + dt*(uB*sin(hB(n)));
end

%% Step 2 (S2): A posteriori update/Measurement update step

% Noise-free measurement: 
% Apply measurement equation to particles and calculate which measurement you would expect.
% This represents the probability of the measurement (e.g. your current measurement) given the prior states.

% Size of noisefree measurement variables: 4 x N - For each particle there are 4 measurements
z_noiseFree_correctRobot = zeros(4,N); z_noiseFree_wrongRobot = zeros(4,N);
f_zm_xp = zeros(4,N);

% TODO: Get rid of for...
for n = 1:N
    % For all particles
    z_noiseFree_correctRobot(1,n)   = sqrt((xA_P(n) - L)^2 + yA_P(n)^2);        % s(1,n) is 1: measure robot A
    z_noiseFree_wrongRobot(1,n)     = sqrt((xB_P(n) - L)^2 + yB_P(n)^2);        % s(1,n) is 0: measure robot B
    
    z_noiseFree_correctRobot(2,n)   = sqrt((xA_P(n) - L)^2 + (yA_P(n) - L)^2);  % s(2,n) is 1: measure robot A
    z_noiseFree_wrongRobot(2,n)     = sqrt((xB_P(n) - L)^2 + (yB_P(n) - L)^2);  % s(2,n) is 0: measure robot B
    
    z_noiseFree_correctRobot(3,n)   = sqrt(xB_P(n)^2 + (yB_P(n) - L)^2);        % s(3,n) is 1: measure robot B
    z_noiseFree_wrongRobot(3,n)     = sqrt(xA_P(n)^2 + (yA_P(n) - L)^2);        % s(3,n) is 0: measure robot A
    
    z_noiseFree_correctRobot(4,n)   = sqrt(xB_P(n)^2 + yB_P(n)^2);              % s(1,n) is 1: measure robot B
    z_noiseFree_wrongRobot(4,n)     = sqrt(xA_P(n)^2 + yA_P(n)^2);              % s(1,n) is 0: measure robot A   

    for sensId = 1:4
        % Measurment likelihood: 4 x N
        % Using actual measurements sens(:)
        
        % TODO: Consider cases where we don't have all the sensor
        % measurements (INF)
        f_zm_xp(sensId,n) = GetProbabilityOutOfTriangularPDF(z_noiseFree_correctRobot(sensId,n),sens(sensId)).*(1-KC.sbar) + ...
                            GetProbabilityOutOfTriangularPDF(z_noiseFree_wrongRobot(sensId,n),sens(sensId)).*KC.sbar;
    end
end

% Calculate normalization constants alphas - size(alpha): 4 x 1
%
% For every sensor, you have a normalization constant.
% Calculate measurement likelihood beta which will be used to represent f_xM lateron
% beta should have size 1 x N
alpha = 1./sum(f_zm_xp,2);
beta = prod(diag(alpha)*f_zm_xp,1);% .* prod(alpha,1);

% Resampling
% ---------------------------------
% Initialize variables which will be assigned after measurement update
xA_M = zeros(N,1); yA_M = zeros(N,1); hA_M = zeros(N,1);
xB_M = zeros(N,1); yB_M = zeros(N,1); hB_M = zeros(N,1);

% Draw N uniform samples r_n
n_bar = ones(6);
accumSum = zeros(6,1);

% TODO: Get rid of for...
% TODO: Use "cumsum" to calculate accumulated sum and "find" to access cumsum
for n = 1:N
    r = rand(6,1);
    
    % Pick particle n_bar, such that:
    % sum(beta_n,n,1,n_bar) >= r_n and
    % sum(beta_n,n,1,n_bar-1) < r_n
    for k = 1:6
        for ind = 1:N
            if (accumSum(k) < r(k)) && ((accumSum(k) + beta(ind)) >= r(k))
                n_bar(k) = ind;
                break;
            else
                accumSum(k) = accumSum(k) + beta(ind);
            end
        end
    end
    
    xA_M(n) = xA_P(n_bar(1));
    yA_M(n) = yA_P(n_bar(2));
    hA_M(n) = hA_P(n_bar(3));
    
    xB_M(n) = xB_P(n_bar(4));
    yB_M(n) = yB_P(n_bar(5));
    hB_M(n) = hB_P(n_bar(6));
end % for...n

% Sample Impoverishment: Roughening
% ----------------------------------
% Perturb the particles after resampling

% TODO: Use formula provided in slides
% TODO: Normal distribution with zero-mean
% and std-dev K*E_i*N^(-1/d),
% K: tuning parameter; d: 6; E_i: maximim inter-sample variability (maximal
% difference of the particles); N^(-1/d): spacing between nodes of a
% corresponding uniform, square grid.
perturb = 0.01;

xA_M = xA_M + perturb*(rand(N,1)-0.5*ones(N,1));
yA_M(:) = yA_M + perturb*(rand(N,1)-0.5*ones(N,1));
hA_M(:) = hA_M + perturb*(rand(N,1)-0.5*ones(N,1));

xB_M(:) = xB_M + perturb*(rand(N,1)-0.5*ones(N,1));
yB_M(:) = yB_M + perturb*(rand(N,1)-0.5*ones(N,1));
hB_M(:) = hB_M + perturb*(rand(N,1)-0.5*ones(N,1));

% Assign to new variables
postParticles.x(1,n) = xA_M(n);
postParticles.y(1,n) = yA_M(n);
postParticles.h(1,n) = hA_M(n);

postParticles.x(2,n) = xB_M(n);
postParticles.y(2,n) = yB_M(n);
postParticles.h(2,n) = hB_M(n);

function newHeading = newHeading(oldHeading, oldX, oldY, oldU, noiseV)
    newHeading = oldHeading;
    
    % Upper wall:
    if (oldY == L) && (oldU*sin(oldHeading) > 0)
        if oldU*cos(oldHeading) > 0
            bounce = oldHeading; % alpha angle
            bounce = bounce*(1 + noiseV);
            newHeading = -bounce;
        else
            bounce = pi - oldHeading;
            bounce = bounce*(1 + noiseV);
            newHeading = -pi + bounce;
        end
    end
    
    % Right wall:
    if (oldX == L) && (oldU*cos(oldHeading) > 0)
        if oldU*sin(oldHeading) > 0
            bounce = 0.5*pi - oldHeading;
            bounce = bounce*(1 + noiseV);
            newHeading = 0.5*pi + bounce;
        else
            bounce = -(-0.5*pi - oldHeading);
            bounce = bounce*(1 + noiseV);
            newHeading = -0.5*pi - bounce;
        end
    end
    
    % Lower wall:
    if (oldY == 0) && (oldU*sin(oldHeading) < 0)
        if oldU*cos(oldHeading) > 0
            bounce = -oldHeading;
            bounce = bounce*(1 + noiseV);
            newHeading = bounce;
        else
            bounce = -(-pi - oldHeading);
            bounce = bounce*(1 + noiseV);
            newHeading = pi - bounce;
        end
    end
    
    % Left wall:
    if (oldX == 0) && (oldU*cos(oldHeading) < 0)
        if oldU*sin(oldHeading) > 0
            bounce = oldHeading - 0.5*pi;
            bounce = bounce*(1 + noiseV);
            newHeading = 0.5*pi - bounce;
        else
            bounce = -(-pi - oldHeading);
            bounce = bounce*(1 + noiseV);
            newHeading = -0.5*pi + bounce;
        end
    end    
end

function probab = GetProbabilityOutOfTriangularPDF(centerOfPDF,evaluationPoint)
    probab = 0;
    if evaluationPoint >= centerOfPDF - KC.wbar && ...
       evaluationPoint < centerOfPDF
        % linear INCREASING
        probab = 1/(KC.wbar^2)*(evaluationPoint+centerOfPDF)+1/KC.wbar;
    
    elseif evaluationPoint >= centerOfPDF && ...
       evaluationPoint <= centerOfPDF + KC.wbar
        % linear DECREASING
        probab = -1/(KC.wbar^2)*(evaluationPoint+centerOfPDF)+1/KC.wbar;
    end
end

% Noise used for bouncing angle
function qRV = drawQuadraticRVSample(size)
    u = rand(size);
    c = 3/(2*(KC.vbar^3));
    qRV = zeros(size);
    for index = 1:size
        if (3/c)*u(index) - KC.vbar^3 > 0
            qRV(index) = abs(((3/c)*u(index) - KC.vbar^3)^(1/3));
        else
            qRV(index) = -abs(((3/c)*u(index) - KC.vbar^3)^(1/3));
        end
    end
end

end % end estimator

