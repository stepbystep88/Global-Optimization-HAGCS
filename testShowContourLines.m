%% This script is used to show the contour map of all test functions 
% Author: Bin She
% Affiliation: University of Electronic Science and Technology of China
% Email: bin.stepbystep@gmail.com
% Date: September, 2019

clear;
clc;
dimensions = 2;
plotSet = bsGetDefaultPlotSet();
plotSet.fontsize = 11;
plotSet.fontweight = 'normal';

benchmarks = {    
    @bsSphere, 'Sphere', @(dim)(zeros(dim, 1)), @(dim)(0), [-10 10], dimensions;
    @bsSchwefel2_22, "Schwefel's P2.22", @(dim)(zeros(dim, 1)), @(dim)(0), [-10 10], dimensions;
    @bsRosenbrock, 'Rosenbrock', @(dim)(ones(dim, 1)), @(dim)(0), [-10 10], dimensions;
    {@bsCFSingle, @(dim)([0.1 1 0.5])}, 'CF1', @(dim)(ones(dim, 1)), @(dim)(0), [-10 10], dimensions;
    {@bsStochasticRosenbrock, @(dim)(abs(rand(dim, 1)) + 0.5)}, 'Stochastic Rosenbrock', @(dim)(ones(dim, 1)), @(dim)(0), [-10 10], dimensions;
    @bsXinSheYang, 'Xin-She Yang', @(dim)(zeros(dim, 1)), @(dim)(0), [-10 10], dimensions;
    @bsAckley, 'Ackley', @(dim)(zeros(dim, 1)), @(dim)(0), [-10 10], dimensions;
    @bsGriewank, 'Griewank', @(dim)(zeros(dim, 1)), @(dim)(0), [-10 10], dimensions;
    @bsRastrigin, 'Rastrigin', @(dim)(zeros(dim, 1)), @(dim)(0), [-10 10], dimensions; 
    {@bsShiftedRastrigin, @(dim)([5; 4])}, 'Shifted Rastrigin', @(dim)(zeros(dim, 1)), @(dim)(0), [-10 10], dimensions;
    @bsSchwefel, 'Schwefel', @(dim)(420.9687*ones(dim, 1)), @(dim)(0), [-500 500], dimensions;
    {@bsCFComplex, @(dim)([abs(rand(5, 1)+0.1)])}, 'CF2', @(dim)(ones(dim, 1)), @(dim)(0), [-10 10], dimensions;
};

nBenchmarks = size(benchmarks, 1);
% data = cell(1, nBenchmarks);
nDim = 2;

figure;
% set(gcf, 'position', [408         242        1043         718]);
set(gcf, 'position', [208         242        1243         718]);

for i = 1 : nBenchmarks
    
    bsSubPlotFit(3, 4, i, 0.95, 0.95, 0.03, 0.07, 0.05, 0.03);
    
    range = benchmarks{i, 5};
    
    x = linspace(range(1), range(2), 200);
    y = x;
    nx = length(x);
    ny = nx;
    
    [X, Y] = meshgrid(x, y);
    
    if isa( benchmarks{i, 1}, 'function_handle')
        objFunc = benchmarks{i, 1};
    else
        fcnpkgs = benchmarks{i, 1};
        fcn1 = fcnpkgs{1};
        fcn2 = fcnpkgs{2};
        o = fcn2(nDim);

        objFunc = @(x,y)(fcn1(x, y, o));
    end
                
    input = [reshape(X, 1, nx*ny); reshape(Y, 1, nx*ny)];
    output = objFunc(input, 0);
    Z = reshape(output, nx, ny);
    
    if i >= 3 && i <=5
        % we plot log scale for the third to fifth test functions 
        contour(X,Y,log(Z)); hold on;
    else
        contour(X,Y,Z); hold on;
    end
    
    if i == 3 || i == 5
        % the global optima of the third and 5-th test function is (1, 2),
        % i.e., x = 1, y = 1
        plot(1, 1, 'rp', 'linewidth', 1);
    elseif i == 11
        % the global optima of the 11-th test function is (420.9687, 420.9687)
        plot(420.9687, 420.9687, 'rp', 'linewidth', 1);   
    elseif i == 10
        % the global optima of the 10-th test function is (5, 4), it may be
        % different in different simulation, here, we just fixed it as (5, 4)
        plot(5, 4, 'rp', 'linewidth', 1);   
    else
        plot(0, 0, 'rp', 'linewidth', 1);
    end
    
    bsPlotSetDefault(plotSet);
    
    xlabel(sprintf('$\\textbf{(%s)}$ - $f_{%d}~x$', 'a'+i-1, i), 'fontweight', 'bold', 'interpreter', 'latex');

    if mod(i, 4) == 1
        ylabel(sprintf('$y$ '), 'fontweight', 'bold', 'interpreter', 'latex');
    end
end