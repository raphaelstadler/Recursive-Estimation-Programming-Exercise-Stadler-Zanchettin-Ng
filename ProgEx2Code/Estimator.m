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
N = 5000; % obviously, you will need more particles than 10.    
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

vA = drawQuadraticRVSample([1, N]); % draw noise from quadratic pdf for post-bounce angles
vB = drawQuadraticRVSample([1, N]);

% DYNAMIC OF THE SYSTEM: Process Equation
%
% Propagate N particles x_m through process dynamics, to get new
% particles x_p
hA_P = newHeading(hA,xA,yA,uA,vA); % new hA needs old xA and old yA, so update hA first
hB_P = newHeading(hB,xB,yB,uB,vB); % new hB needs old xB and old yB, so update hB first

[xA, yA] = shiftParticlesToValidBounds(xA,yA);
[xB, yB] = shiftParticlesToValidBounds(xB,yB);

xA_P = xA + dt*(uA*cos(hA));
yA_P = yA + dt*(uA*sin(hA));

xB_P = xB + dt*(uB*cos(hB));
yB_P = yB + dt*(uB*sin(hB));

%% Step 2 (S2): A posteriori update/Measurement update step

%if true % Uncomment this, to skip measurement update completely
if sens == Inf*ones(size(sens))
    % No sensor measurements available:
    % Completely skip measurement update step   
    postParticles.x(1,:) = xA_P;
    postParticles.y(1,:) = yA_P;
    postParticles.h(1,:) = hA_P;

    postParticles.x(2,:) = xB_P;
    postParticles.y(2,:) = yB_P;
    postParticles.h(2,:) = hB_P;
    return;
end

% Store which sensor values are meaningful
validRowsSensors = find(sens ~= Inf);
validRowsSensorA = intersect([1;2],validRowsSensors);
validRowsSensorB = intersect([3;4],validRowsSensors);

% Noise-free measurement: 
% Apply measurement equation to particles and calculate which measurement you would expect.
% This represents the probability of the measurement (e.g. the current measurement) given the prior states.

% Size of noisefree measurement variables: 4 x N - For each particle there are 4 measurements
z_noiseFree_correctRobot = zeros(4,N); z_noiseFree_wrongRobot = zeros(4,N);
f_zm_xp = zeros(4,N);

% For all particles
z_noiseFree_correctRobot(1,:)   = sqrt((xA_P - L).^2 + yA_P.^2);        % s(1,:) is 1: measure robot A
z_noiseFree_wrongRobot(1,:)     = sqrt((xB_P - L).^2 + yB_P.^2);        % s(1,:) is 0: measure robot B

z_noiseFree_correctRobot(2,:)   = sqrt((xA_P - L).^2 + (yA_P - L).^2);  % s(2,:) is 1: measure robot A
z_noiseFree_wrongRobot(2,:)     = sqrt((xB_P - L).^2 + (yB_P - L).^2);  % s(2,:) is 0: measure robot B

z_noiseFree_correctRobot(3,:)   = sqrt(xB_P.^2 + (yB_P - L).^2);        % s(3,:) is 1: measure robot B
z_noiseFree_wrongRobot(3,:)     = sqrt(xA_P.^2 + (yA_P - L).^2);        % s(3,:) is 0: measure robot A

z_noiseFree_correctRobot(4,:)   = sqrt(xB_P.^2 + yB_P.^2);              % s(1,:) is 1: measure robot B
z_noiseFree_wrongRobot(4,:)     = sqrt(xA_P.^2 + yA_P.^2);              % s(1,:) is 0: measure robot A   

for sensId = 1:4
    % Measurment likelihood: 4 x N
    % Using actual measurements sens
    % Note, that in the case where not all sensor measurements were
    % available, sens(.) contains another meaningful value.
    f_zm_xp(sensId,:) = GetProbabilityOutOfTriangularPDF(z_noiseFree_correctRobot(sensId,:), sens(sensId)).*(1-KC.sbar) + ...
                        GetProbabilityOutOfTriangularPDF(z_noiseFree_wrongRobot(sensId,:), sens(sensId)).*KC.sbar;
end

% Calculate normalization constants alphas - size(alpha): 4 x 1
%
% For every sensor, you have a normalization constant.
% Calculate measurement likelihood beta which will be used to represent f_xM lateron
% beta should have size 1 x N

% Note: Row sum of f_zm_xp can be zero, which results in alpha = NaN
alpha_test = 1./sum(f_zm_xp,2);
validRowsProbab = find(alpha_test ~= Inf & ~isnan(alpha_test));

validRowsProbabA = intersect([1;2],validRowsProbab);
validRowsProbabB = intersect([3;4],validRowsProbab);

alpha1 = 1./sum(f_zm_xp(validRowsProbabA,:),2);
beta(1,:) = prod(diag(alpha1)*f_zm_xp(validRowsProbabA,:),1);
alpha2 = 1./sum(f_zm_xp(validRowsProbabB,:),2);
beta(2,:) = prod(diag(alpha2)*f_zm_xp(validRowsProbabB,:),1);

sumCloseToOne = sum(beta, 2) > 0.9;

if (isempty(alpha1) || (sumCloseToOne(1) ~= 1)) % beta(1,:) is not > 0
    validRowsProbab = intersect(validRowsProbab, [3;4]);
    beta(1,:) = zeros(1,N);
end
if (isempty(alpha2) || (sumCloseToOne(2) ~= 1)) % beta(2,:) is not > 0
    validRowsProbab = intersect(validRowsProbab, [1;2]);
    beta(2,:) = zeros(1,N);
end

validRowsProbabA = intersect([1;2],validRowsProbab);
validRowsProbabB = intersect([3;4],validRowsProbab);

doMeasurementUpdate = ones(2,1);

%% ------------------------------------------------------------------------
% Robot A
%
if isempty(validRowsSensorA)
    % Skip measurement update
    % No sensor measurements available for robot A
    % Skip measurement update step for robot A  
    postParticles.x(1,:) = xA_P;
    postParticles.y(1,:) = yA_P;
    postParticles.h(1,:) = hA_P;
    
    doMeasurementUpdate(1) = 0;
else
    % Measurements available for robot A
    if ~isempty(validRowsProbabA)
        % beta for robot A is meaningful an can be taken for resampling
        % afterwards
        
        % TODO: It is also possible that only for 1 of the sensors the
        % likelihood (beta) is zero
    else
        % beta for robot A is meaningless
        % Redistribute particles according to measurement
        
        % Assign equal probabililites to all the new distr particles
        beta(1,:) = 1/N;
        
        unif = rand(2,N);
        %hA_P = unif(1,:).*2*pi - pi;
        
        if length(validRowsSensorA) == 1
            % Only 1 sensor measurement for robot A available
            angle = zeros(2,N);
            angle(1,:) = pi/2*unif(1,:)+pi/2;
            angle(2,:) = pi/2*unif(2,:)-pi;
            
            rR = drawTriangularRVSample([2,N]);
            
            if intersect(1, validRowsSensorA)
                % Distribute particles around sensor 1
                xA_P = (sens(1)+rR(1,:)).*cos(angle(1,:)) + L;
                yA_P = (sens(1)+rR(1,:)).*sin(angle(1,:));
                
                % TODO: Perform KC.sbar part of the particles should be sampled from (xB_P, yB_P)
            elseif intersect(2, validRowsSensorA)
                % Distribute particles around sensor 2
                xA_P = (sens(2)+rR(2,:)).*cos(angle(2,:)) + L;
                yA_P = (sens(2)+rR(2,:)).*sin(angle(2,:)) + L;
                
                % TODO: Perform KC.sbar part of the particles should be sampled from (xB_P, yB_P)
            else
                error('Invalid state of particles.');
            end            
        elseif length(validRowsSensorA) == 2
             % 2 sensor measurements for robot A available
  
             if (sens(1) + sens(2)) >= L  
                  
                postParticles.x(1,:) = xA_P;
                postParticles.y(1,:) = yA_P;
                postParticles.h(1,:) = hA_P;

                doMeasurementUpdate(1) = 0;
%                 % If measurement circles intersect: Distribute particles around circle intersection
%                 yA = (sens(1)^2-sens(2)^2+L^2)/(2*L);
% 
%                 xA1 = (2*L+sqrt(4*L^2-4*(L^2+yA^2-sens(1)^2)))/2; % +
%                 xA2 = (2*L-sqrt(4*L^2-4*(L^2+yA^2-sens(1)^2)))/2; % -
% 
%                 if (xA1 >= 0 && xA1 <= L)
%                     xA = xA1;
%                 elseif (xA2 >= 0 && xA2 <= L)
%                     xA = xA2;
%                 else
%                     error('The triangulation failed. Invalid use of formula.');
%                 end
% 
%                 xA_P = xA*ones(1,N);
%                 yA_P = yA*ones(1,N);
% 
%                 % Add noise / blur
%                 rA = sqrt(xA_P.^2 + yA_P.^2);
%                 gamma = asin(yA_P ./ rA);
%                 rR = drawTriangularRVSample([1, N]);
%                 
%                 xA_P = (rA + rR).*cos(gamma);
%                 yA_P = (rA + rR).*sin(gamma);
            else
                % Distribute particles around the 2 measurment quarter-circles
                unif = rand(2,N);
                angle = zeros(2,N);
                
                angle(1,:) = pi/2*unif(1,:)+pi/2;
                angle(2,:) = pi/2*unif(2,:)-pi/2;
                rR = drawTriangularRVSample([2,N]);

                % Distribute half of particles around sensor 1
                xA_P(1:N_half) = (sens(1)+rR(1,1:N_half)).*cos(angle(1,1:N_half)) + L;
                yA_P(1:N_half) = (sens(1)+rR(1,1:N_half)).*sin(angle(1,1:N_half));
                % And half of particles around sensor 2
                xA_P(N_half+1:N) = (sens(2)+rR(2,N_half+1:N)).*cos(angle(2,N_half+1:N)) + L;
                yA_P(N_half+1:N) = (sens(2)+rR(2,N_half+1:N)).*sin(angle(2,N_half+1:N)) + L; 
            end
        else
            error('Invalid State of Sensor A');
        end
    end
end

%% ------------------------------------------------------------------------
% Robot B
%
if isempty(validRowsSensorB)
    % Skip measurement update
    % No sensor measurements available for robot B
    % Skip measurement update step for robot B  
    postParticles.x(2,:) = xB_P;
    postParticles.y(2,:) = yB_P;
    postParticles.h(2,:) = hB_P;
    
    doMeasurementUpdate(2) = 0;
else
    % Measurements available for robot B
    if ~isempty(validRowsProbabB)
        % beta for robot A is meaningful an can be taken for resampling
        % afterwards
    else
        % beta for robot B is meaningless
        % Redistribute particles according to measurement
        
        % Assign equal probabililites to all the new distr particles
        beta(2,:) = 1/N;
        
        unif = rand(2,N);
        %hB_P = unif(1,:).*2*pi - pi;
        
        if length(validRowsSensorB) == 1
            % Only 1 sensor measurement for robot B available
            %angle = pi/2*unif(2,:)+pi/2;
            angle = zeros(2,N);
            angle(1,:) = pi/2*unif(1,:);
            angle(2,:) = pi/2*unif(2,:)-pi/2;
            
            rR = drawTriangularRVSample([2,N]);
            
            if intersect(3, validRowsSensorB)
                % Distribute particles around sensor 3
                xB_P = (sens(3)+rR(1,:)).*cos(angle(1,:));
                yB_P = (sens(3)+rR(1,:)).*sin(angle(1,:)) + L;
            elseif intersect(4, validRowsSensorB)
                % Distribute particles around sensor 4
                xB_P = (sens(4)+rR(2,:)).*cos(angle(2,:));
                yB_P = (sens(4)+rR(2,:)).*sin(angle(2,:));
            else
                error('Invalid state of particles.');
            end
        elseif length(validRowsSensorB) == 2     
            % 2 sensor measurements for robot B available
            if (sens(3) + sens(4)) >= L
                postParticles.x(2,:) = xB_P;
                postParticles.y(2,:) = yB_P;
                postParticles.h(2,:) = hB_P;

                doMeasurementUpdate(2) = 0;
                
%                 % If measurement circles intersect: Distribute particles around circle intersection
%                 yB = (sens(4)^2-sens(3)^2+L^2)/(2*L);
%                 xB = sqrt(sens(4)^2-yB);
% 
%                 xB_P = xB*ones(1,N);
%                 yB_P = yB*ones(1,N);
%                 
%                 % Add noise / blur
%                 rB = sqrt(xB_P.^2 + yB_P.^2);
%                 gamma = asin(yB_P ./ rB);
%                 rR = drawTriangularRVSample([1, N]);
%                 
%                 xB_P = (rB + rR).*cos(gamma);
%                 yB_P = (rB + rR).*sin(gamma);
            else
                % Distribute particles around the 2 measurment quarter-circles
                unif = rand(2,N);
                angle = zeros(2,N);
                
                angle(1,:) = pi/2*unif(1,:)-pi/2;
                angle(2,:) = pi/2*unif(2,:);
                
                rR = drawTriangularRVSample([2,N]);

                % Distribute half of particles around sensor 3
                xB_P(1:N_half) = (sens(3)+rR(1,1:N_half)).*cos(angle(1,1:N_half));
                yB_P(1:N_half) = (sens(3)+rR(1,1:N_half)).*sin(angle(1,1:N_half)) + L;
                % And half of particles around sensor 4
                xB_P(N_half+1:N) = (sens(4)+rR(2,N_half+1:N)).*cos(angle(2,N_half+1:N));
                yB_P(N_half+1:N) = (sens(4)+rR(2,N_half+1:N)).*sin(angle(2,N_half+1:N));
            end
        else
            error('Invalid State of Sensor B');
        end
    end
end

%     % Put KC.sbar percent of the particles xA to xB and
%     % KC.sbar percent of the particles from xB to xA
%     permutedInd = randperm(N);
%     deltaInd = ceil(KC.sbar*N);
% 
%     xA_P(permutedInd(1:deltaInd))   = xB_P(permutedInd(1:deltaInd));
%     yA_P(permutedInd(1:deltaInd))   = yB_P(permutedInd(1:deltaInd));
% 
%     xB_P(permutedInd(N-deltaInd:N)) = xA_P(permutedInd(N-deltaInd:N));
%     yB_P(permutedInd(N-deltaInd:N)) = yA_P(permutedInd(N-deltaInd:N));


%% ------------------------------------------------------------------------
% Resampling
% ---------------------------------
% Initialize variables which will be assigned after measurement update
xA_M = zeros(1,N); yA_M = zeros(1,N); hA_M = zeros(1,N);
xB_M = zeros(1,N); yB_M = zeros(1,N); hB_M = zeros(1,N);

% build cumulative sum of particle measurement likelihood
cumulativeSum = cumsum(beta,2); % 2 x N

n_bar = zeros(2,1);

% Draw N uniform samples r and choose subset of xA_P according to
% cumulative pdf of beta
for i = 1:N
    r = rand(2,1);
    
    if doMeasurementUpdate(1) == 1        
        n_bar(1) = find(cumulativeSum(1,:) >= r(1),1,'first');
        xA_M(i) = xA_P(n_bar(1));
        yA_M(i) = yA_P(n_bar(1));
        hA_M(i) = hA_P(n_bar(1));    
    end
    
    if doMeasurementUpdate(2) == 1
        n_bar(2) = find(cumulativeSum(2,:) >= r(2),1,'first');
        xB_M(i) = xB_P(n_bar(2));
        yB_M(i) = yB_P(n_bar(2));
        hB_M(i) = hB_P(n_bar(2)); 
    end
end

%% ------------------------------------------------------------------------
% Sample Impoverishment: Roughening
% ----------------------------------
% Perturb the particles after resampling and assign to new variables
[xA_M, yA_M, hA_M, ...
 xB_M, yB_M, hB_M] = performRoughening(xA_M, yA_M, hA_M, xB_M, yB_M, hB_M);

if doMeasurementUpdate(1) == 1
    postParticles.x(1,:) = xA_M;
    postParticles.y(1,:) = yA_M;
    postParticles.h(1,:) = hA_M;
end
if doMeasurementUpdate(2) == 1
    postParticles.x(2,:) = xB_M;
    postParticles.y(2,:) = yB_M;
    postParticles.h(2,:) = hB_M;
end

function newHeading = newHeading(oldHeading, oldX, oldY, oldU, noiseV)
    newHeading = oldHeading;
    
    % Upper wall:
    upperWallCollInd = intersect(find(oldY >= L), find(oldU*sin(oldHeading) > 0));
    if ~isempty(upperWallCollInd)
        specificCollisionType = find(oldU*cos(oldHeading) > 0);
        if intersect(upperWallCollInd, specificCollisionType)
            bounce = oldHeading(upperWallCollInd); % alpha angle
            bounce = bounce.*(1 + noiseV(upperWallCollInd));
            newHeading(upperWallCollInd) = -bounce;
        else
            bounce = oldHeading(upperWallCollInd); % alpha angle
            bounce = bounce.*(1 + noiseV(upperWallCollInd));
            newHeading(upperWallCollInd) = -pi + bounce;
        end
            
    end
        
    % Right wall:
    rightWallCollInd = intersect(find(oldX >= L), find(oldU*cos(oldHeading) > 0));
    if ~isempty(rightWallCollInd)
        specificCollisionType = find(oldU*sin(oldHeading) > 0);
        if intersect(rightWallCollInd, specificCollisionType)
            bounce = 0.5*pi - oldHeading(rightWallCollInd);
            bounce = bounce.*(1 + noiseV(rightWallCollInd));
            newHeading(rightWallCollInd) = 0.5*pi + bounce;
        else
            bounce = -(-0.5*pi*ones(size(oldHeading(rightWallCollInd))) - oldHeading(rightWallCollInd));
            bounce = bounce.*(1 + noiseV(rightWallCollInd));
            newHeading(rightWallCollInd) = -0.5*pi - bounce;
        end
    end
    
    % Lower wall:
    lowerWallCollInd = intersect(find(oldY <= 0), find(oldU*sin(oldHeading) < 0));
    if ~isempty(lowerWallCollInd)
        specificCollisionType = find(oldU*cos(oldHeading) > 0);
        if intersect(lowerWallCollInd, specificCollisionType)
            bounce = -oldHeading(lowerWallCollInd);
            bounce = bounce.*(1 + noiseV(lowerWallCollInd));
            newHeading(lowerWallCollInd) = bounce;
        else
            bounce = -(-pi -oldHeading(lowerWallCollInd));
            bounce = bounce.*(1 + noiseV(lowerWallCollInd));
            newHeading(lowerWallCollInd) = pi - bounce;
        end
    end
    
    % Left wall:
    leftWallCollInd = intersect(find(oldX <= 0), find(oldU*cos(oldHeading) < 0));
    if ~isempty(leftWallCollInd)
        specificCollisionType = find(oldU*sin(oldHeading) > 0);
        if intersect(lowerWallCollInd, specificCollisionType)
            bounce = oldHeading(leftWallCollInd) - 0.5*pi;
            bounce = bounce.*(1 + noiseV(leftWallCollInd));
            newHeading(leftWallCollInd) = 0.5*pi - bounce;
        else
            bounce = -(-pi - oldHeading(leftWallCollInd) - 0.5*pi);
            bounce = bounce.*(ones(size(bounce)) + noiseV(leftWallCollInd));
            newHeading(leftWallCollInd) = -0.5*pi + bounce;
        end
    end    
end

function probab = GetProbabilityOutOfTriangularPDF(centerOfPDF,evaluationPoint)
    % Note: Function is designed to handle vector valued centerOfPDF, and scalar valued evaluationPoint
    probab = zeros(size(centerOfPDF));
    
    linearIncrInd = find(evaluationPoint >= (centerOfPDF - KC.wbar) & ...
                         evaluationPoint < centerOfPDF);
    if ~isempty(linearIncrInd)
        % linear INCREASING
        probab(linearIncrInd) = 1/(KC.wbar^2)*(evaluationPoint-centerOfPDF(linearIncrInd))+1./KC.wbar;
    end
    
    linearDecrInd = find(evaluationPoint >= centerOfPDF & ...
                         evaluationPoint <= (centerOfPDF + KC.wbar));
    if ~isempty(linearDecrInd)
        % linear DECREASING
        probab(linearDecrInd) = -1/(KC.wbar^2)*(evaluationPoint-centerOfPDF(linearDecrInd))+1./KC.wbar;
    end
end

% Noise used for bouncing angle
function qRV = drawQuadraticRVSample(size)
    u = rand(size);
    c = 3/(2*(KC.vbar^3));
    qRV = zeros(size);
    for index = 1:size(2)
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
    
    xValid(xTest < 0) = -xTest(xTest < 0);
    xValid(xTest > L) = L *ones(size(xTest(xTest > L))) - ( xTest(xTest > L) - L*ones(size(xTest(xTest > L))) );
    
    yValid(yTest < 0) = -yTest(yTest < 0);
    yValid(yTest > L) = L *ones(size(yTest(yTest > L))) - ( yTest(yTest > L) - L*ones(size(yTest(yTest > L))) );
end

function [xA_r, yA_r, hA_r, xB_r, yB_r, hB_r] = performRoughening(xA, yA, hA, xB, yB, hB)    
    % Normal distribution with zero-mean
    % and std-dev K*E_i*N^(-1/d),
    completeMatrix = [xA;yA;hA;xB;yB;hB];

    % K: tuning parameter
    K = 0.001;
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

function tRV = drawTriangularRVSample(size)
    u = rand(size);
    a = -KC.wbar;
    b = KC.wbar;
    c = 0;
    tRV = zeros(size);
    for ind = 1:size(2)
        if u(ind) < (c-a)/(b-a)
            tRV(ind) = a + sqrt(u(ind)*(b-a)*(c-a));
        else
            tRV(ind) = b - sqrt((1-u(ind))*(b-a)*(b-c));        
        end
    end
end

end % end estimator

