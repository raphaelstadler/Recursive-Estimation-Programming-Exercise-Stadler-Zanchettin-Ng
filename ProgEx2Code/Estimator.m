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

%% Step 1 (S1): Prior update/Prediction step

% Define synonyms for easier understanding
xA = prevPostParticles.x(1,:);
yA = prevPostParticles.y(1,:);
hA = prevPostParticles.h(1,:);

xB = prevPostParticles.x(2,:);
yB = prevPostParticles.y(2,:);
hB = prevPostParticles.h(2,:);

uA = act(1);
uB = act(2);

vA = drawQuadraticRVSample(N); % draw noise from quadratic pdf for post-bounce angles
vB = drawQuadraticRVSample(N);

% Initialize variables which will be assigned after prior update
xA_P = zeros(1,N); yA_P = zeros(1,N); hA_P = zeros(1,N);
xB_P = zeros(1,N); yB_P = zeros(1,N); hB_P = zeros(1,N);

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

% Represent prior PDF with the Monte-Carlo sampling (as done above) of the
% the measurement variables x_m

% Please note, that f_xP is in fact never used !

[f_xP,~] = ApproximatePDF(...
            [0;0;-pi;0;0;-pi], ...                  % min       of xA_P, yA_P, hA_P, xB_P, yB_P, hB_P
            [L;L;pi;L;L;pi], ...                    % max       of xA_P, yA_P, hA_P, xB_P, yB_P, hB_P
            100, ...                                % numOfBins of xA_P, yA_P, hA_P, xB_P, yB_P, hB_P
            [xA_P; yA_P; hA_P; xB_P; yB_P; hB_P] ); % sample vector


%% Step 2 (S2): A posteriori update/Measurement update step

% Noise-free measurement: 
% Apply measurement equation to particles and calculate which measurement you would expect.
% This represents the probability of any measurement (e.g. your current
% measurement) given the prior state variable.

z_correctRobot = zeros(4,N); % 4 x N: For each particle there are 4 measurements

for n = 1:N
    % For all particles
    w = drawTriangularRVSample(4);  % parameter defines the vector length of the RV
    s = drawBooleanRVSample(4);     % parameter defines whether the sensors detected the correct robot.
    
    z_correctRobot(1,n) = s(1)*sqrt((xA_P(n) - L)^2 + yA_P(n)^2) + (1 - s(1))*sqrt((xB_P(n) - L)^2 + yB_P(n)^2) + w(1);
    z_correctRobot(2,n) = s(2)*sqrt((xA_P(n) - L)^2 + (yA_P(n) - L)^2) + (1 - s(2))*sqrt((xB_P(n) - L)^2 + (yB_P(n) - L)^2) + w(2);
    z_correctRobot(3,n) = s(3)*sqrt(xB_P(n)^2 + (yB_P(n) - L)^2) + (1 - s(3))*sqrt(xA_P(n)^2 + (yA_P(n) - L)^2) + w(3);
    z_correctRobot(4,n) = s(4)*sqrt(xB_P(n)^2 + yB_P(n)^2) + (1 - s(4))*sqrt(xA_P(n)^2 + yA_P(n)^2) + w(4);
end

% Approximate PDF of z_correctRobot (which was calculated given x_p)
% size(f_zm_xp): 4 x numberOfBins
[f_zm_xp, z_grid] = ApproximatePDF(...
            [0;0;0;0], ...              % min       of sensors S1, S2, S3, S4
            sqrt(2)*L*[1;1;1;1], ...    % max       of sensors S1, S2, S3, S4
            100, ...                    % numOfBins of sensors S1, S2, S3, S4
            z_correctRobot );           % sample vector

% Calculate normalization constants alphas
% size(alpha): 4 x 1
% For every state variable (in this case sensor), you have a normalization constant.
alpha = 1./sum(f_zm_xp,2);

beta = zeros(6,N);

% Calculate measurement likelihood beta_n which will be used to represent f_xm lateron
% Beta should have size 4 x N
% Using actual measurements sens(1:4)
for n = 1:N
    % How likely is measurement sense(1:4) given xP
    probab_sens_xp = GetProbabilityFromPDF(f_zm_xp, z_grid, sens); % 4 x 1
    
    % TODO: How to reflect dependency of beta with n
    beta(1,n) = probab_sens_xp(1) * alpha(1); % sens(1)
    beta(2,n) = probab_sens_xp(2) * alpha(2); % sens(2)
    beta(3,n) = probab_sens_xp(3) * alpha(3); % sens(3)
    beta(4,n) = probab_sens_xp(4) * alpha(4); % sens(4)
end

% Perform scaling of original x_p to finally represent PDF of x_m
% Probability of the prior variables given the measurement

% TODO: scale with beta (matrix dimensions are currently wrong and all
% information of beta should be used.

% xA_m = xA .* beta(1,:);
% yA_m = yA .* beta(1,:);
% hA_m = hA .* beta(1,:);
% 
% xB_m = xB .* beta(2,:);
% yB_m = yB .* beta(2,:);
% hB_m = hB .* beta(2,:);
xA_m = xA;
yA_m = yA;
hA_m = hA;

xB_m = xB;
yB_m = yB;
hB_m = hB;

% Now represent PDF of x_m
[f_xm,~] = ApproximatePDF(...
            [0;0;-pi;0;0;-pi], ...                  % min       of xA_m, yA_m, hA_m, xB_m, yB_m, hB_m
            [L;L;pi;L;L;pi], ...                    % max       of xA_m, yA_m, hA_m, xB_m, yB_m, hB_m
            100, ...                                % numOfBins  of xA_m, yA_m, hA_m, xB_m, yB_m, hB_m
            [xA_m; yA_m; hA_m; xB_m; yB_m; hB_m] ); % sample vector

% Actual measurement
z_bar = sens;

% Resampling
% ---------------------------------
xA_M = zeros(N,1);
yA_M = zeros(N,1);
hA_M = zeros(N,1);

xB_M = zeros(N,1);
yB_M = zeros(N,1);
hB_M = zeros(N,1);

% Draw N uniform samples r_n
n_bar = ones(6);
accumSum = zeros(6,1);
for n = 1:N
    r = rand(6,1);
    
    % Pick particle n_bar, such that:
    % sum(beta_n,n,1,n_bar) >= r_n and
    % sum(beta_n,n,1,n_bar-1) < r_n

    for k = 1:6
        for ind = 1:N
            % TODO: Fix indexing error of beta
            % Should be something like this:
            % if (accumSum(k) < r(k)) && ((accumSum(k) + beta(k,1)) >= r(k))
            % AND
            % accumSum(k) = accumSum(k) + beta(k,ind);
            if (accumSum(k) < r(k)) && ((accumSum(k) + beta(k,1)) >= r(k))
            % if (accumSum(k) < r(k)) && ((accumSum(k) + beta(k,ind)) >= r(k))
                n_bar(k) = ind;
                break;
            else
                accumSum(k) = accumSum(k) + beta(k,1);
                % accumSum(k) = accumSum(k) + beta(k,ind);
                
            end
        end
    end
    
    xA_M(n) = xA_m(n_bar(1));
    yA_M(n) = yA_m(n_bar(2));
    hA_M(n) = hA_m(n_bar(3));
    
    xB_M(n) = xB_m(n_bar(4));
    yB_M(n) = yB_m(n_bar(5));
    hB_M(n) = hB_m(n_bar(6));
end % for...n

% Assign to new variables
% Theses variables are considered as the
% x_p (prior variables)
postParticles.x(1,n) = xA_M(n);
postParticles.y(1,n) = yA_M(n);
postParticles.h(1,n) = hA_M(n);

postParticles.x(2,n) = xB_M(n);
postParticles.y(2,n) = yB_M(n);
postParticles.h(2,n) = hB_M(n);

% Sample Impoverishment: Roughening
% ----------------------------------
% Perturb the particles after resampling

function newHeading = newHeading(oldHeading, oldX, oldY, oldU, noiseV)
    
    alpha = -1;
    newHeading = oldHeading;
    % TODO: Check if orientation changes because of a bouncing and add noise to heading
    
    % Upper wall:
    if (oldY == L) && (oldU*sin(oldHeading) > 0)
        if oldU*cos(oldHeading) > 0
            alpha = oldHeading;
            alpha = alpha*(1 + noiseV);
            newHeading = -alpha;
        else
            alpha = pi - oldHeading;
            alpha = alpha*(1 + noiseV);
            newHeading = -pi + alpha;
        end
    end
    
    % Right wall:
    if (oldX == L) && (oldU*cos(oldHeading) > 0)
        if oldU*sin(oldHeading) > 0
            alpha = 0.5*pi - oldHeading;
            alpha = alpha*(1 + noiseV);
            newHeading = 0.5*pi + alpha;
        else
            alpha = -(-0.5*pi - oldHeading);
            alpha = alpha*(1 + noiseV);
            newHeading = -0.5*pi - alpha;
        end
    end
    
    % Lower wall:
    if (oldY == 0) && (oldU*sin(oldHeading) < 0)
        if oldU*cos(oldHeading) > 0
            alpha = -oldHeading;
            alpha = alpha*(1 + noiseV);
            newHeading = alpha;
        else
            alpha = -(-pi - oldHeading);
            alpha = alpha*(1 + noiseV);
            newHeading = pi - alpha;
        end
    end
    
    % Left wall:
    if (oldX == 0) && (oldU*cos(oldHeading) < 0)
        if oldU*sin(oldHeading) > 0
            alpha = oldHeading - 0.5*pi;
            alpha = alpha*(1 + noiseV);
            newHeading = 0.5*pi - alpha;
        else
            alpha = -(-pi - oldHeading);
            alpha = alpha*(1 + noiseV);
            newHeading = -0.5*pi + alpha;
        end
    end    
end

function [f_vec, gridVec] = ApproximatePDF(minVec, maxVec, numOfBins, sampleVec)
    % Please note, that with the current implementation, the numOfBins is a
    % scalar (it has to be the same for all sample types (e.g. xA,yA,hA,...))
    
    numOfSamples = size(sampleVec,2);
    
    f_vec = zeros(length(minVec),numOfBins);
    
    deltaVec = (maxVec-minVec)./numOfBins; %dX, dH
    gridVec = zeros(length(deltaVec), numOfBins);
    
    for i = 1:length(deltaVec)
        gridVec(i,:) = minVec(i)+0.5*deltaVec(i):deltaVec(i):maxVec(i)-0.5*deltaVec(i);
    end
    
    for sampleId = 1:numOfSamples
        for vecId = 1:length(minVec) % for xA, yA, hA, xB, yB, hB
            for binId = 1:numOfBins-1
               if   (sampleVec(vecId,sampleId) >= (binId-1)*deltaVec(vecId)) && ...
                    (sampleVec(vecId,sampleId) < binId*deltaVec(vecId))
               
                    f_vec(vecId,binId) = f_vec(vecId,binId) + 1/numOfSamples;
                
               end
            end
            
            % Specially treat last bin to also include right border of bin
            if   (sampleVec(vecId,sampleId) >= (numOfBins-1)*deltaVec(vecId)) && ...
                (sampleVec(vecId,sampleId) <= numOfBins*deltaVec(vecId))

                f_vec(vecId,numOfBins) = f_vec(vecId,numOfBins) + 1/numOfSamples;
            end 
        end % for...vecId
    end % for...sampleId     
end % end of function

function probabValues = GetProbabilityFromPDF(pdf_func, gridVec, evaluatePoints)
    if size(pdf_func,1) ~= length(evaluatePoints)
        error('Dimensions of PDF and vectors of points where PDF is evaluated have to comfirm');
        return;
    end

    % This implementation assumes constant deltaVec over whole grid
    deltaVec = gridVec(:,2) - gridVec(:,1);
    probabValues = zeros(length(deltaVec),1);
    probabInd = ones(length(deltaVec),1);
    
    for vecId = 1:length(evaluatePoints)
        % Walk through gridVec to see which value should be taken to
        % extract probability
        for binId = 1:size(pdf_func,2)-1
            if evaluatePoints(vecId) >= gridVec(binId) && ...
               evaluatePoints(vecId) < gridVec(binId+1)
           
                dist_ID = norm(evaluatePoints(vecId) - gridVec(binId),2);
                dist_IDplusOne = norm(evaluatePoints(vecId) - gridVec(binId+1),2);
                
                if (dist_ID < dist_IDplusOne)
                    probabInd(vecId) = binId;
                else
                    probabInd(vecId) = binId+1;
                end
                break;
            end % if 
        end % for...binId
        
        % Assign actual probability values to output variable
        probabValues(vecId) = pdf_func(vecId, probabInd(vecId)) / deltaVec(vecId);
        
    end % for...vecId
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

function bRV = drawBooleanRVSample(size)
    u = rand(size);
    bRV = (u > KC.sbar);
end

end % end estimator

