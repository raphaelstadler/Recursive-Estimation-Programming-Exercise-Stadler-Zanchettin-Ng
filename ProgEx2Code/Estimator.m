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
N = 100; % obviously, you will need more particles than 10.    
N_half = floor(N/2);
if (init)
    % Initialization of estimator:
    % These particles are the posterior particles at discrete time k = 0
    % which will be fed into your estimator again at k = 1
    % Replace the following:
    
    % At time t = 0:
    % - robot A starts out in one of the corners where the sensors S_1 = (L,0) and S_2 = (L,L) are located
    % - robot B starts out in one of the corners where the sensors S_3 = (0,L) and S_4 = (0,0) are located
    
    % RV in -pi/4..pi/4 and function handle used for initial heading assignment
    angle_RV = rand(2,N)*0.5*pi-0.25*pi; 
    
    getRandomHeading = @(rowNum, startIndCol, endIndCol, offset) angle_RV(rowNum,startIndCol:endIndCol) + offset*ones(1,(endIndCol-startIndCol)+1);
    
    % A at sensors S_1 and B at sensor S_3
    postParticles.x(:,1:N_half)  = repmat([L; 0],1,N_half);
    postParticles.y(:,1:N_half)  = repmat([0; L],1,N_half);
    % Draw a uniform RV to get initial headings theta_A, theta_B
    postParticles.h(:,1:N_half)  = [getRandomHeading(1,1,N_half,3*pi/4); getRandomHeading(2,1,N_half,-pi/4)];
    
    % A at sensors S_2 and B at sensor S_4
    postParticles.x(:,(N_half+1):N)  = repmat([L;0],1,N-N_half);
    postParticles.y(:,(N_half+1):N)  = repmat([L;0],1,N-N_half);
    % Draw a uniform RV to get initial headings theta_A, theta_B
    postParticles.h(:,(N_half+1):N)  = [getRandomHeading(1,N_half+1,N,-3*pi/4); getRandomHeading(2,N_half+1,N,pi/4)];
    
    % Leave the function
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
xA_P = xA; yA_P = yA; hA_P = hA;
xB_P = xB; yB_P = yB; hB_P = hB;

% TODO: Optimize newHeading(.) function to be able to treat a vector as a
% whole

% Apply the process equation to the particles
for n = 1:N
    % DYNAMIC OF THE SYSTEM: Process Equation
    %
    % Propagate N particles x_m through process dynamics, to get new
    % particles x_p

    % TODO: Adapt newHeading function: Currently even at the beginning
    % there is a detected collision!
    hA_P(n) = newHeading(hA(n),xA(n),yA(n),uA,vA(n)); % new hA needs old xA and old yA, so update hA first
    hB_P(n) = newHeading(hB(n),xB(n),yB(n),uB,vB(n)); % new hB needs old xB and old yB, so update hB first
    
    [xA(n), yA(n)] = shiftParticlesToValidBounds(xA(n),yA(n));
    [xB(n), yB(n)] = shiftParticlesToValidBounds(xB(n),yB(n));
    
    xA_P(n) = xA(n) + dt*(uA*cos(hA(n)));
    yA_P(n) = yA(n) + dt*(uA*sin(hA(n)));
    
    xB_P(n) = xB(n) + dt*(uB*cos(hB(n)));
    yB_P(n) = yB(n) + dt*(uB*sin(hB(n)));
    
    [xA(n), yA(n)] = shiftParticlesToValidBounds(xA(n),yA(n));
    [xB(n), yB(n)] = shiftParticlesToValidBounds(xB(n),yB(n));
end

%% Step 2 (S2): A posteriori update/Measurement update step

% if true % Uncomment this, to skip measurement update completely
if sens == Inf*ones(size(sens))
    % No sensor measurements available:
    % Completely skip measurement update step   
    postParticles.x(1,:) = xA_P(:);
    postParticles.y(1,:) = yA_P(:);
    postParticles.h(1,:) = hA_P(:);

    postParticles.x(2,:) = xB_P(:);
    postParticles.y(2,:) = yB_P(:);
    postParticles.h(2,:) = hB_P(:);
    return;
end

xA_P_mean = mean(xA_P); yA_P_mean = mean(yA_P);
xB_P_mean = mean(xB_P); yB_P_mean = mean(yB_P);

% Substitute missing measurements for robot A
if sens(1) == Inf && sens(2) == Inf
    sens(1) = sqrt((xA_P_mean - L)^2 + yA_P_mean^2);
    sens(2) = sqrt((xA_P_mean - L)^2 + (yA_P_mean - L)^2);
elseif sens(1) == Inf && sens(2) ~= Inf
    sens(1) = sqrt((xA_P_mean - L)^2 + yA_P_mean^2);
elseif sens(1) ~= Inf && sens(2) == Inf
    sens(2) = sqrt((xA_P_mean - L)^2 + (yA_P_mean - L)^2);
end

% Substitute missing measurements for robot B
if sens(3) == Inf && sens(4) == Inf
    sens(3) = sqrt(xB_P_mean^2 + (yB_P_mean - L)^2);
    sens(4) = sqrt(xB_P_mean^2 + yB_P_mean^2);
elseif sens(3) ~= Inf && sens(4) == Inf
    sens(3) = sqrt(xB_P_mean^2 + (yB_P_mean - L)^2);
elseif sens(3) == Inf && sens(4) ~= Inf
    sens(4) = sqrt(xB_P_mean^2 + yB_P_mean^2);
end

% Noise-free measurement: 
% Apply measurement equation to particles and calculate which measurement you would expect.
% This represents the probability of the measurement (e.g. the current measurement) given the prior states.

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
        % Note, that in the case where not all sensor measurements were
        % available, sens(.) contains another meaningful value.
        f_zm_xp(sensId,n) = GetProbabilityOutOfTriangularPDF(z_noiseFree_correctRobot(sensId,n),sens(sensId)).*(1-KC.sbar) + ...
                            GetProbabilityOutOfTriangularPDF(z_noiseFree_wrongRobot(sensId,n),sens(sensId)).*KC.sbar;
    end
end

% Calculate normalization constants alphas - size(alpha): 4 x 1
%
% For every sensor, you have a normalization constant.
% Calculate measurement likelihood beta which will be used to represent f_xM lateron
% beta should have size 1 x N

% TODO: row sum of f_zm_xp can be zero, which results in alpha = NaN

alpha_test = 1./sum(f_zm_xp,2);
validRows = find(alpha_test ~= Inf & ~isnan(alpha_test));

beta = zeros(1,N);
if ~isempty(validRows)
    % At least one of the rows of f_zm_xp is non-zero
    alpha = 1./sum(f_zm_xp(validRows,:),2);
    beta = prod(diag(alpha)*f_zm_xp(validRows,:),1);
end

if(sum(beta) > 0)
    beta = beta/sum(beta);
else
    % All rows of f_zm_xp are zero
    
    % Use the sensor measurements to reassamble particles in regions
    % meaningful for the sensor values
    unif = rand(4,N);
    
    hA_P = unif(2,:).*2*pi - pi*ones(1,N);
    hB_P = unif(4,:).*2*pi - pi*ones(1,N);
    if (sens(1) >= L/2) && (sens(2) >= L/2)
        % The "measurement circles" of the two sensors intercept 
        yA = (sens(1)^2-sens(2)^2+L^2)/(2*L);
        
        xA1 = (2*L+sqrt(4*L^2-4*(L^2+yA^2-sens(1)^2)))/2;
        xA2 = (2*L-sqrt(4*L^2-4*(L^2+yA^2-sens(1)^2)))/2;
        
        if (xA1 >= 0)
            xA = xA1;
        elseif (xA2 >= 0)
            xA = xA2;
        else
            error('The triangulation failed. Invalid use of formula.');
        end
        
        xA_P = xA*ones(1,N);
        yA_P = yA*ones(1,N);        
    else
        % The "measurement circles" do not intercept
        % Assemble half of the particles around quarter-circle of sensor 3
        % and the others around quarter-circle of sensor 4
        
        % First half of particles around sensor 1
        angle = pi/2*unif+pi/2;
        xA_P(1:N_half) = sens(1).*cos(angle(1,1:N_half)) + L*ones(1,N_half);
        yA_P(1:N_half) = sens(1).*sin(angle(1,1:N_half));
        
        % Second half of particles around sensor 2
        angle = pi/2*unif+pi;
        xA_P(N_half+1:N) = sens(2).*cos(angle(1,N_half+1:N)) + L*ones(1,N_half);
        yA_P(N_half+1:N) = sens(2).*sin(angle(1,N_half+1:N)) + L*ones(1,N_half);
    end
    
    if (sens(3) >= L/2) && (sens(4) >= L/2)
        % The "measurement circles" of the two sensors intercept
        yB = (sens(4)^2-sens(3)^2+L^2)/(2*L);
        xB = sqrt(sens(4)^2-yB);

        xB_P = xB*ones(1,N);
        yB_P = yB*ones(1,N); 
    else
        % The "measurement circles" do not intercept
        % Assemble half of the particles around quarter-circle of sensor 3
        % and the others around quarter-circle of sensor 4
        
        % First half of particles around sensor 3
        angle = pi/2*unif+3*pi/2;
        xA_P(1:N_half) = sens(3).*cos(angle(2,1:N_half));
        yA_P(1:N_half) = sens(3).*sin(angle(2,1:N_half)) + L*ones(1,N_half);
        
        % Second half of particles around sensor 4
        angle = pi/2*unif;
        xA_P(N_half+1:N) = sens(4).*cos(angle(2,N_half+1:N));
        yA_P(N_half+1:N) = sens(4).*sin(angle(2,N_half+1:N));
    end

    % Put KC.sbar percent of the particles xA to xB and
    % KC.sbar percent of the particles from xB to xA
    permutedInd = randperm(N);
    deltaInd = ceil(KC.sbar*N);

    xA_P(permutedInd(1:deltaInd))   = xB_P(permutedInd(1:deltaInd));
    yA_P(permutedInd(1:deltaInd))   = yB_P(permutedInd(1:deltaInd));

    xB_P(permutedInd(N-deltaInd:N)) = xA_P(permutedInd(N-deltaInd:N));
    yB_P(permutedInd(N-deltaInd:N)) = yA_P(permutedInd(N-deltaInd:N));
    
    beta = 1/N*ones(1,N);
end

% Resampling
% ---------------------------------
% Initialize variables which will be assigned after measurement update
xA_M = zeros(N,1); yA_M = zeros(N,1); hA_M = zeros(N,1);
xB_M = zeros(N,1); yB_M = zeros(N,1); hB_M = zeros(N,1);

% build cumulative sum of particle measurement likelihood
cumulativeSum = cumsum(beta);

% Draw N uniform samples r and choose subset of xA_P according to
% cumulative pdf of beta
for i = 1:N
    r = rand;
    n_bar = find(cumulativeSum >= r,1,'first');

    xA_M(i,1) = xA_P(n_bar);
    yA_M(i,1) = yA_P(n_bar);
    hA_M(i,1) = hA_P(n_bar);

    xB_M(i,1) = xB_P(n_bar);
    yB_M(i,1) = yB_P(n_bar);
    hB_M(i,1) = hB_P(n_bar); 
end

% Sample Impoverishment: Roughening
% ----------------------------------
% Perturb the particles after resampling and assign to new variables
[postParticles.x(1,:), postParticles.y(1,:), postParticles.h(1,:), ...
 postParticles.x(2,:), postParticles.y(2,:), postParticles.h(2,:)] = performRoughening(xA_M, yA_M, hA_M, xB_M, yB_M, hB_M);

function newHeading = newHeading(oldHeading, oldX, oldY, oldU, noiseV)
    newHeading = oldHeading;
    
    % Upper wall:
    if (oldY >= L) && (oldU*sin(oldHeading) > 0)
        bounce = oldHeading; % alpha angle
        bounce = bounce*(1 + noiseV);
        newHeading = -bounce;
    end
    
    % Right wall:
    if (oldX >= L) && (oldU*cos(oldHeading) > 0)
        bounce = 0.5*pi - oldHeading;
        bounce = bounce*(1 + noiseV);
        newHeading = 0.5*pi + bounce;
    end
    
    % Lower wall:
    if (oldY <= 0) && (oldU*sin(oldHeading) < 0)
        bounce = -oldHeading;
        bounce = bounce*(1 + noiseV);
        newHeading = bounce;
    end
    
    % Left wall:
    if (oldX <= 0) && (oldU*cos(oldHeading) < 0)
        bounce = oldHeading - 0.5*pi;
        bounce = bounce*(1 + noiseV);
        newHeading = 0.5*pi - bounce;
    end    
end

function probab = GetProbabilityOutOfTriangularPDF(centerOfPDF,evaluationPoint)
    probab = 0;
    if evaluationPoint >= (centerOfPDF - KC.wbar) && ...
       evaluationPoint < centerOfPDF
        % linear INCREASING
        probab = 1/(KC.wbar^2)*(evaluationPoint-centerOfPDF)+1/KC.wbar;
    elseif evaluationPoint >= centerOfPDF && ...
       evaluationPoint <= (centerOfPDF + KC.wbar)
        % linear DECREASING
        probab = -1/(KC.wbar^2)*(evaluationPoint-centerOfPDF)+1/KC.wbar;
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

function [xValid, yValid] = shiftParticlesToValidBounds(xTest, yTest)
    xValid = xTest;
    yValid = yTest;
    if xTest < 0
        xValid = -xTest;
    elseif xTest > L
        xValid = L - (xTest-L);
    end
    if yTest < 0
        yValid = -yTest;
    elseif yTest > L
        yValid = L - (yTest-L);
    end
end

function [xA_r, yA_r, hA_r, xB_r, yB_r, hB_r] = performRoughening(xA, yA, hA, xB, yB, hB)    
    % Normal distribution with zero-mean
    % and std-dev K*E_i*N^(-1/d),
    completeMatrix = [xA';yA';hA';xB';yB';hB'];

    % K: tuning parameter
    K = 0.3;
    % D: Dimension of state space
    D = 6;
    % E_i: maximim inter-sample variability
    Ei = max(completeMatrix,[],2) - min(completeMatrix,[],2);
    % N^(-1/d): spacing between nodes of a corresponding uniform, square grid.

    sigma = K*diag(Ei)*N^(-1/D)*ones(6,N); % std dev
    mu = zeros(6,N);             % expectation

    perturb = normrnd(mu,sigma,6,N);
    
    completeMatrix = completeMatrix + perturb;
    
    xA_r = completeMatrix(1,:); yA_r = completeMatrix(2,:); hA_r = completeMatrix(3,:);
    xB_r = completeMatrix(4,:); yB_r = completeMatrix(5,:); hB_r = completeMatrix(6,:);    
end
end % end estimator

