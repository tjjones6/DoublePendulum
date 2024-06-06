%% Double Pendulum RK45 Fractal Evolution Map
% Author:      Tyler Jones
% Contact:     tjjones6@wisc.edu
% Institution: University of Wisconsin-Madison
% Date:          05.27.2024
%
% Description: 
% This script iteratively calculates the time it takes for the
% second pendulum to 'flip' based off of a matrix of initial angular positions. 
% The goal is to visualize and create order from chaos in its purest form.

clear all; close all; clc; 

%% User Input
g = 9.8035;  % Acceleration due to gravity (m/s^2)
L1 = 1.0;  % Length of the first pendulum arm (m)
L2 = 1.0;  % Length of the second pendulum arm (m)
m1 = 1.0;  % Mass of the first pendulum bob (kg)
m2 = 1.0;  % Mass of the second pendulum bob (kg)

% Initial Conditions
fidelity = 750;
theta1 = linspace(-pi,pi,fidelity);  % Initial angle of the first pendulum (radians)
theta2 = linspace(-pi,pi,fidelity);  % Initial angle of the second pendulum (radians)
omega1_0 = 0;      % Initial angular velocity of the first pendulum (rad/s)
omega2_0 = 0;      % Initial angular velocity of the second pendulum (rad/s)

% Time stepping increment 
time_step = 0.25; % (s)
t_final = 25; % (s)
num_iter = t_final/time_step;

figure('units','normalized','Position',[0.1 0.1 .8 .8])
% myWriter = VideoWriter('DoublePendulumFractalv2.mp4', 'MPEG-4');
% myWriter.FrameRate = 10;
% open(myWriter);

for k = 1:num_iter
    clf

    % Time Settings
    tspan = [0 time_step*k];  % Extend the time span for simulation (seconds)

    % Initialize a matrix to store the flip times
    flipTimes = NaN(fidelity, fidelity);

    %% Numerical Integration using ODE45
    for i = 1:fidelity
        % Clear command window and update user on runtime
        clc
        fprintf('Current Time Step Completed: %.2f%%\n', i/fidelity*100)
        fprintf('Total Simulation Completed: %.2f%%\n', k/num_iter*100-1)
        
        for j = 1:fidelity
            % Set initial conditions for the current pair of angles
            initial_conditions = [theta1(i), omega1_0, theta2(j), omega2_0];

            % Use ODE45 to integrate until the event (flip) occurs
            options = odeset('RelTol', 1e-4, 'AbsTol', 1e-4, 'Events', @stopEvent);
            [~, ~, te, ~, ~] = ode45(@(t, y) doublePendulumODE(t, y, g, L1, L2, m1, m2), tspan, initial_conditions, options);

            % Store the time of the event (flip) if it occurred
            if ~isempty(te)
                flipTimes(i, j) = te(1);
            end
        end
    end

    %% Plotting Results
    %Plotting the fractal map using a heatmap
    % figure('units','normalized','Position',[0.1 0.1 .8 .8])
    hold on; axis tight;
    myfigpref
    imagesc(theta1, theta2, flipTimes)
    colormap jet; colorbar;
    fig_xytit('$\theta_1$ [rad]','$\theta_2$ [rad]','Double Pendulum Fractal')
    hold off

    drawnow

%     frame = getframe(gcf);
%     writeVideo(myWriter,frame);
end
% close(myWriter)

%% Define the function for the double pendulum ODEs
function dydt = doublePendulumODE(~, y, g, L1, L2, m1, m2)
    % Unpack the state variables
    theta1 = y(1);
    omega1 = y(2);
    theta2 = y(3);
    omega2 = y(4);

    % Store equations of motion
    dydt = zeros(4, 1);
    dydt(1) = omega1;
    dydt(2) = (-g * (2 * m1 + m2) * sin(theta1) - m2 * g * sin(theta1 - 2 * theta2) - 2 * sin(theta1 - theta2) * m2 * (omega2^2 * L2 + omega1^2 * L1 * cos(theta1 - theta2))) / (L1 * (2 * m1 + m2 - m2 * cos(2 * theta1 - 2 * theta2)));
    dydt(3) = omega2;
    dydt(4) = (2 * sin(theta1 - theta2) * (omega1^2 * L1 * (m1 + m2) + g * (m1 + m2) * cos(theta1) + omega2^2 * L2 * m2 * cos(theta1 - theta2))) / (L2 * (2 * m1 + m2 - m2 * cos(2 * theta1 - 2 * theta2)));
end

%% Event function to stop the integration
function [value, isterminal, direction] = stopEvent(~, y)
    theta2 = y(3);  % Angle of the second pendulum
    value = abs(theta2) - 2 * pi;  % Condition for the event
    isterminal = 1;  % Stop the integration
    direction = 0;  % Detect both increasing and decreasing directions
end

%% Reference Functions for Plotting
function myfigpref
    set(0, 'DefaultAxesFontSize', 20);
    set(0, 'DefaultAxesLineWidth', 2);
    set(0, 'DefaultLineLineWidth', 2);
    set(0, 'DefaultPatchLineWidth', .7);
    set(0, 'DefaultLineMarkerSize', 6);
    grid on;
    box on;
    h = gca;
    h.TickLabelInterpreter='latex';
    h.MinorGridAlpha=0.05;
    h.GridAlpha=0.05;
    h.FontSize=25;
    h.LineWidth=2;
    h = gcf;
    h.Color = [1,1,1];
end

function fig_xytit(xlab, ylab, tit)
    if nargin < 3
        tit = '';
    end
    xlabel(xlab, 'interpreter', 'latex');
    ylabel(ylab, 'interpreter', 'latex');
    title(tit, 'interpreter', 'latex');
end
