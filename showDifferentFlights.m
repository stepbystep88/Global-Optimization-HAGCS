%% This script is used to compare the movement of different flights
% Author: Bin She
% Affiliation: University of Electronic Science and Technology of China
% Email: bin.stepbystep@gmail.com
% Date: September, 2019
close all;
clear;

nMove = 500;
gaussianFlight = 2 * randn(2, nMove);
levyFlight = bsLevy(2, nMove, 1.5);

gf_points = zeros(2, nMove);
lf_points = zeros(2, nMove);

for i = 2 : nMove
    gf_points(:, i) = gf_points(:, i-1) + gaussianFlight(:, i);
    lf_points(:, i) = lf_points(:, i-1) + levyFlight(:, i);
end

%% show the results
bsNewFigure;
p1 = plot(gf_points(1, :), gf_points(2, :), 'b.-', 'linewidth', 2); hold on;
p2 = plot(lf_points(1, :), lf_points(2, :), 'r.-', 'linewidth', 2);



plot(gf_points(1, :), gf_points(2, :), 'k.', 'linewidth', 2);
plot(lf_points(1, :), lf_points(2, :), 'k.', 'linewidth', 2);
legend([p1, p2], 'Gaussian Flight', "Levy Flight");

bsPlotSetDefault(bsGetDefaultPlotSet());