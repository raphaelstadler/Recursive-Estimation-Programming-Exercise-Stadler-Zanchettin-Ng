function [meanError, tEstAvg] = run(setSeed, showPlots)
% [meanError, tEstAvg] = run(setSeed, showPlots)
%
% Main function for the Particle Filter exercise, which is used to simulate
% the actual robot, call the estimator, and plot and evaluate the results.
%
% You must not edit this file.
%
%   Inputs (optional):
%   setSeed:            Boolean variable to enable resetting of random
%                       seed. If setSeed == 1, the seed is reset.
%                       Otherwise, the seed is not reset. If the variable
%                       is not provided, the default value is setSeed = 0.
%
%   showPlots:          Boolean variable to enable plotting of results.
%                       If showPlots == 0, no plots are shown. This may be
%                       useful when running several experiments that
%                       produce the performance and timing outputs. By
%                       default, showPlots = 1, and the plots are shown.
%
%   Outputs:
%   meanError:          Performance measure in meters. The mean error is
%                       computed as follows: At each time step, the
%                       Euclidean distance of all robot position particles
%                       to their respective actual robot is computed. Then,
%                       the distances are sorted in ascending order, and
%                       the mean of the 50% of the smallest distances is
%                       computed. The meanError is then the mean over these
%                       means at all time steps. The heading error is not
%                       evaluated.
%
%   tEstAvg:            Average execution time in seconds of a single
%                       estimator update, averaged over all time steps.
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

% use default inputs
if(nargin < 2)
    % by default, show plots
    showPlots = 1;
    if(nargin < 1)
        % by default, do not reset the random number generator seed
        setSeed = 0;
    end
end

% reset seed, if requested
if(setSeed)
    stream = RandStream.getGlobalStream;
    reset(stream);
end

%% Simulation: Generate data
% Simulation of the system to obtain measurements and the true
% states.  The measurements will be visible to the estimator.
[roboA, roboB, meas, act] = simulateRobots();


%% Estimator
% Run the estimator.

% Initialize the estimator.
init = 1;
[posterior] = Estimator([],[],[],init);

% Get number of particles:
Nparticles = size(posterior.x,2);

% initialize storage arrays for x, y, and heading h:
robotsX = zeros(2,Nparticles,UKC.simSteps + 1);
robotsX(:,:,1) = posterior.x;
robotsY = zeros(2,Nparticles,UKC.simSteps + 1);
robotsY(:,:,1) = posterior.y;
robotsH = zeros(2,Nparticles,UKC.simSteps + 1);
robotsH(:,:,1) = posterior.h;

% Call the estimator for each time step and measure the total execution
% time (including the updates of the storage arrays):
if showPlots
    disp('Generating data...')
end
nextDispFracComplete = 0.1;
tstart = tic;
for k = 2:(UKC.simSteps + 1)
    if showPlots
        if (k/(UKC.simSteps + 1)) > nextDispFracComplete
            nextDispFracComplete = nextDispFracComplete + 0.1;
            disp([num2str(round((k/(UKC.simSteps + 1))*100)),'% complete'])
        end
    end
    [posterior] = Estimator(posterior, meas(:,k), act(:,k-1));
    % update storage arrays:
    robotsX(:,:,k) = posterior.x;
    robotsY(:,:,k) = posterior.y;
    robotsH(:,:,k) = posterior.h;
end
tEstAvg = toc(tstart)/UKC.simSteps;
if showPlots
    disp('Done. Making plots.')
end

%% The results
% Plot particles and evaluate mean error of 50% particles closest to actual
% robot (just positions).

% mean error array:
err = zeros(1, UKC.simSteps + 1);

% Relevant User Parameters:
delay = 0.1;        % user settable delay between frames, slow down or
                    % or speed up visualization with this parameter
SHOWSENSORS = 1;    % show the sensor measurements (use for debugging)
                    % the green circle is the sensor measurement, and the
                    % red circles are +/- the width of the triangular
                    % distribution

% Plotting Parameters:
colorRobotParticlesA = 'm'; % color of particles for robot A
colorActualRobotA = 'm';    % color of actual robot A
colorRobotParticlesB = 'b'; % color of particles for robot B
colorActualRobotB = 'b';    % color of actual robot B
axisScaling = 0.1;          % scaling for axes, to show some extra space
                            % about the bounding box defined by KC.L
nbins = 40;                 % number of bins in the histograms

% the colours of the sensor's current measurement
sensorColours = [colorActualRobotA, colorActualRobotA, colorActualRobotB, colorActualRobotB];     

if(SHOWSENSORS) % set up circle arrays for sensor visualization
    theta = linspace(0,2*pi,100);
    xcircle = cos(theta);
    ycircle = sin(theta);
end

% some preparation of figure axes
if(showPlots)
    try
        close(1)
    catch
    end
    figure(1)
    axis([-axisScaling*KC.L, (1+ axisScaling)*KC.L, -axisScaling*KC.L, (1+ axisScaling)*KC.L])
    ax = axis;
    set(gca,'DataAspectRatio',[1 1 1])
end
for k = 1:(UKC.simSteps + 1)
    if(showPlots)
        figure(1)
        subplot(2,2,[1,3])
        hold off
        hpa = plot(robotsX(1,:,k),robotsY(1,:,k),'Color',colorRobotParticlesA,'Marker','.','LineStyle','none');
        hold on
        hpb = plot(robotsX(2,:,k),robotsY(2,:,k),'Color',colorRobotParticlesB,'Marker','.','LineStyle','none');
        ha = plot(roboA(1,k),roboA(2,k),'Color', 'k', 'MarkerFaceColor',colorActualRobotA,'Marker','o','LineWidth',2,'LineStyle','none');
        hb = plot(roboB(1,k),roboB(2,k),'Color', 'k', 'MarkerFaceColor',colorActualRobotB,'Marker','o','LineWidth',2,'LineStyle','none');
        if(any(~isinf(meas(:,k))) && SHOWSENSORS)
            if(~isinf(meas(1,k)))
                plot(meas(1,k)*xcircle+KC.L,meas(1,k)*ycircle,sensorColours(1))
                plot((KC.wbar + meas(1,k))*xcircle+KC.L,(KC.wbar + meas(1,k))*ycircle, [sensorColours(1), ':'])
                if(-KC.wbar + meas(1,k) > 0)
                    plot((-KC.wbar + meas(1,k))*xcircle+KC.L,(-KC.wbar + meas(1,k))*ycircle,[sensorColours(1), ':'])
                end
            end
            if(~isinf(meas(2,k)))
                plot(meas(2,k)*xcircle+KC.L,meas(2,k)*ycircle+KC.L,sensorColours(2))
                plot((KC.wbar + meas(2,k))*xcircle+KC.L,(KC.wbar + meas(2,k))*ycircle+KC.L,[sensorColours(2), ':'])
                if(-KC.wbar + meas(2,k) > 0)
                    plot((-KC.wbar + meas(2,k))*xcircle+KC.L,(-KC.wbar + meas(2,k))*ycircle+KC.L,[sensorColours(2), ':'])
                end
            end
            if(~isinf(meas(3,k)))
                plot(meas(3,k)*xcircle,meas(3,k)*ycircle+KC.L,sensorColours(3))
                plot((KC.wbar + meas(3,k))*xcircle,(KC.wbar + meas(3,k))*ycircle+KC.L,[sensorColours(3), ':'])
                if(-KC.wbar + meas(3,k) > 0)
                    plot((-KC.wbar + meas(3,k))*xcircle,(-KC.wbar + meas(3,k))*ycircle+KC.L,[sensorColours(3), ':'])
                end
            end
            if(~isinf(meas(4,k)))
                plot(meas(4,k)*xcircle,meas(4,k)*ycircle,sensorColours(4))
                plot((KC.wbar + meas(4,k))*xcircle,(KC.wbar + meas(4,k))*ycircle,[sensorColours(4), ':'])
                if(-KC.wbar + meas(4,k) > 0)
                    plot((-KC.wbar + meas(4,k))*xcircle,(-KC.wbar + meas(4,k))*ycircle,[sensorColours(4), ':'])
                end
            end
        end
        % plot the bounding box
        line([0 KC.L KC.L 0; KC.L KC.L 0 0], [0 0 KC.L KC.L; 0 KC.L KC.L 0],'Color','k'), hold on
        % set axes:
        axis(ax)
        set(gca,'DataAspectRatio',[1 1 1])
        title(['Robots at Discrete Time Step ', int2str(k-1)])
        xlabel('X (meters)')
        ylabel('Y (meters)')
        legend([ha,hpa,hb,hpb],{'Robot A', 'Particles A', 'Robot B', 'Particles B'}, 'Location', 'southoutside')
        
        % plot heading histograms
        subplot(2,2,2)
        [n,xout] = hist(robotsH(1,:,k), linspace(-pi, pi, nbins));
        bar(xout,n,colorRobotParticlesA)
        hold on;
        plot(roboA(3,k), 0, 'Color', colorActualRobotA, 'Marker','o', 'MarkerSize', 25,'LineWidth',2);
        xlim([-pi, pi]);
        title(['Heading Histogram Robot A, Time Step ', int2str(k-1)])
        xlabel('Angle (rad)')
        hold off;
        subplot(2,2,4)
        [n,xout] = hist(robotsH(2,:,k), linspace(-pi, pi, nbins));
        bar(xout,n,colorRobotParticlesB);
        hold on;
        plot(roboB(3,k), 0, 'Color', colorActualRobotB, 'Marker','o', 'MarkerSize', 25,'LineWidth',2);
        xlabel('Angle (rad)')
        xlim([-pi, pi]);
        title('Heading Histogram Robot B');
        hold off;
        drawnow
        pause(delay);
    end
    % Evaluate average distance error of 50% closest particles:
    % distances of robot A to particles A, and analogous for B
    dA = sqrt(sum([robotsX(1,:,k) - roboA(1,k); robotsY(1,:,k) - roboA(2,k)].^2));
    dB = sqrt(sum([robotsX(2,:,k) - roboB(1,k); robotsY(2,:,k) - roboB(2,k)].^2));
    dAsorted = sort(dA,'ascend'); % sort distances
    dBsorted = sort(dB,'ascend'); % sort distances
    err(k) = mean([dAsorted(1:floor(length(dAsorted)/2)),dBsorted(1:floor(length(dBsorted)/2))]);
end
meanError = mean(err);
end
