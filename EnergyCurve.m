close all; clear all; clc

% Define the range for theta_1
theta_1 = linspace(-pi, pi, 1000);

% Calculate the possible values for the argument of acos
cos_theta_2_arg = 1 - 2*cos(theta_1);

% Only keep values within the valid range for acos
valid_indices = cos_theta_2_arg >= -1 & cos_theta_2_arg <= 1.01;

% Filter the values of theta_1 and cos_theta_2_arg
theta_1_valid = theta_1(valid_indices);
cos_theta_2_arg_valid = cos_theta_2_arg(valid_indices);

% Calculate theta_2 for valid values of theta_1
theta_2_valid = acos(cos_theta_2_arg_valid);

% Create the plot
figure('units','normalized','Position',[0.1 0.1 .8 .8])
myfigpref;

% Plot and fill the area between the curves
fill([theta_1_valid, fliplr(theta_1_valid)], [theta_2_valid, -fliplr(theta_2_valid)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on; grid on;
plot(theta_1_valid, theta_2_valid, 'b', 'LineWidth', 2);
plot(theta_1_valid, -theta_2_valid, 'b', 'LineWidth', 2);

% Set the axis limits and labels
xlim([-pi pi])
ylim([-pi pi])
fig_xytit('$\theta_1$ [rad]','$\theta_2$ [rad]','Energetics Curve: $\theta_2 = \arccos(1-2\cos\theta_1)$')

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
