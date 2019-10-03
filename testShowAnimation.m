%% This script is used to show the animations of the optimization process of different methods
% Author: Bin She
% Affiliation: University of Electronic Science and Technology of China
% Email: bin.stepbystep@gmail.com
% Date: September, 2019
close all;
clear;

%% Some basica settings such as stop criteria
% the maximum number of iteration, we generally don't use this item for the stop criterion 
maxIter = 200;              
% when f - f* < optimalFunctionTolerance, iteration will be stopped. f* is
% the global minimum function value
optimalFunctionTolerance = 1e-10;
% the number of simulations (for each problem) you'd like to run
nSimulation = 1;
% the dimension of the problems you'd like to test. 
dimensions = [2];
% whether to save the middle result of in each simulation
isSaveMiddleRes = false;
isSaveDetailUpdates = true;


%% Test function sets, choose one to show the animation of optimization process
benchmarks = {    
    %1. the objective function handle. It looks like [f, g] = obj(x, isGradient)
    %   see bsSphere.m for example. 
    %   Note that in some cases, for instance, the CF2 function, it needs
    %   to randomly generate a objective function, so the first item is
    %   also a cell data {objective function, extra data}
    %2. the name of this objective function
    %3. the global optimal solution. It is actually a funciton handle which
    %   returns a vector by inputting the dimension of the problem
    %4. the global minimum objecteve function value of the problem
    %5. the range of the parameters
    %6. the dimensions (that will be tested) of the test function

%     @bsSphere, 'Sphere', @(dim)(zeros(dim, 1)), @(dim)(0), [-100 100], dimensions, false;
%     @bsSchwefel2_22, "Schwefel's P2.22", @(dim)(zeros(dim, 1)), @(dim)(0), [-100 100], dimensions, false;
%     @bsRosenbrock, 'Rosenbrock', @(dim)(ones(dim, 1)), @(dim)(0), [-10 10], dimensions, true;
%     {@bsStochasticRosenbrock, @(dim)(abs(rand(dim, 1)) + 0.5)}, 'Stochastic Rosenbrock', @(dim)(ones(dim, 1)), @(dim)(0), [-100 100], dimensions, true;
%     {@bsXinSheYang, @(dim)(abs(rand(dim, 1)))}, 'Xin-She Yang', @(dim)(zeros(dim, 1)), @(dim)(0), [-100 100], dimensions, false;
%     @bsAckley, 'Ackley', @(dim)(zeros(dim, 1)), @(dim)(0), [-20 20], dimensions, false;
%     @bsGriewank, 'Griewank', @(dim)(zeros(dim, 1)), @(dim)(0), [-50 50], dimensions, false;
%     @bsRastrigin, 'Rastrigin', @(dim)(zeros(dim, 1)), @(dim)(0), [-10 10], dimensions, false; 
%     {@bsShiftedRastrigin, @(dim)(rand(dim, 1)-0.5)*50}, 'Shifted Rastrigin', @(dim)(zeros(dim, 1)), @(dim)(0), [-100 100], dimensions, false;
    @bsSchwefel, 'Schwefel', @(dim)(420.9687*ones(dim, 1)), @(dim)(0), [-500 500], dimensions, false;
%     {@bsCFRosenbrock, @(dim)([abs(rand(2, 1)+0.1); abs(rand(dim, 1)) + 0.5])}, 'CF1', @(dim)(ones(dim, 1)), @(dim)(0), [-100 100], dimensions, false;
%     {@bsCFComplex, @(dim)([abs(rand(5, 1)+0.1); abs(rand(dim, 1))])}, 'CF1', @(dim)(ones(dim, 1)), @(dim)(0), [-100 100], dimensions, false;

%     @bsPeak, "MATLAB Peak", @(dim)([0.2283; -1.6255]), @(dim)(-6.551133332835842), [-3 3], dimensions, false;
%     @bsDropWave, "Drop Wave", @(dim)([0; 0]), @(dim)(-1), [-20 20], dimensions, false;
%     @bsMichalewicz, "Michalewicz", @(dim)([2.20290552014618;1.57079632677565]), @(dim)(-1.80130341009855321), [0 4], dimensions, false;
%     @bsSquare, "Quadratic Curve", @(dim)([0; 0]), @(dim)(0), [-2 2], dimensions, false;
};

initNest = 25;

% test methods, you can only compare two methods at one time
testMethods = {
    %1. the name of the method
    %2. note of the method
    %3. special parameters of the corresponding method
    
%     'PSO', 'MATLAB', {100};     % 100 refers to the initial population size of the PSO 
%     'GA', 'MATLAB', {100};      % 100 refers to the initial population size of the PSO 
    % the special parameters of the following 2 methods are 
    % NInit
    % p_a (switching parameters)
    % alpha
    % population initial function
    'CS', 'Yang and Deb', {initNest 0.75 0.05, @bsGenerateInitialPopulationByRandom}; 
%     'GBCS', 'Fateen', {initNest, 0.75, 0.05, @bsGenerateInitialPopulationByRandom};
    % the special parameters of the following 3 methods are 
    % NInit
    % p_a (switching parameters)
    % alpha, population initial function
    % the number of calls of gradient-based local optimization (GBLO) during the optimization process
    % the maximum number iterations of the inner GBLO subroutine 
%     'AGBCS', 'She et al.', {initNest, 0.75, 0.05, @bsGenerateInitialPopulationByRandom, 50, 200};
%     'AHSACS', 'She et al.', {initNest, initNest, 500, @bsGenerateInitialPopulationByRandom};
    'HAGCS', 'She et al.', {initNest, initNest, 500, @bsGenerateInitialPopulationByRandom, 50, 200};
};


%% test
nTestMethods = size(testMethods, 1);
nBenchmarks = size(benchmarks, 1);
nProblem = 0;
for iBenchmark = 1 : nBenchmarks
    nProblem = nProblem + length(benchmarks{iBenchmark, 6});
end
nCases = nProblem * nTestMethods;
detailResults = cell(nBenchmarks, nTestMethods, nSimulation);
for iBenchmark = 1 : nBenchmarks
    % for each test function
    
    objFuncName = benchmarks{iBenchmark, 2};    % name of the test function
    xRange = benchmarks{iBenchmark, 5};         % range of the test function
    dimensions = benchmarks{iBenchmark, 6};     % dimensions of the test function
    
    for iDim = 1 : length(dimensions)
        % for each dimension case of the current test function
        
        nDim = dimensions(iDim);                % the dimension of the current problem
        
        % set the maximum number of funciton evaluations as 10^4*nDim 
        maxFevs = 1e4*nDim;      
        
        % the the lower and upper boundary
        if length(xRange(:)) == 2
            % in this case xRange loos like [lower, upper], lower and upper are the
            % same for all dimension of the parameters
            lower = ones(nDim, 1) * xRange(1);
            upper = ones(nDim, 1) * xRange(2);
        else
            % in this case xRange loos like [lower, upper], lower and upper
            % are two vectors
            lower = xRange(:, 1);
            upper = xRange(:, 2);
        end
        
        % obtain the global optimum solution of the current problem
        if isa(benchmarks{iBenchmark, 3},'function_handle')
            bestX = benchmarks{iBenchmark, 3}(nDim);
        else
            bestX = benchmarks{iBenchmark, 3};
        end
        
        % obtain the global minimum of the current problem
        if isa(benchmarks{iBenchmark, 4},'function_handle')
            minFVal = benchmarks{iBenchmark, 4}(nDim);
        else
            minFVal = benchmarks{iBenchmark, 4};
        end
       
        caseNames = cell(1, nTestMethods);
        
        for iTestMethod = 1 : nTestMethods
            % for each method to solve the curretn problem 
            
            methodName = testMethods{iTestMethod, 1};   % method name
            methodNote = testMethods{iTestMethod, 2};   % method note
            parameters = testMethods{iTestMethod, 3};   % special parameters of the current method
            methodInfo = sprintf('%s_%s', methodName, methodNote);
            
            % obtain the number of initial population size
            if isa(parameters{1},'function_handle')
                nInitNest = parameters{1}(nDim);
            else
                nInitNest = parameters{1};
            end
            
            
            for i = 1 : nSimulation
                % for each simulation of the problem

                % prepare the objective function
                if isa( benchmarks{iBenchmark, 1}, 'function_handle')
                    objFunc = benchmarks{iBenchmark, 1};
                else
                    % in this case, the objective funciton needs extra
                    % data, for exameple, see the CF2 test function
                    fcnpkgs = benchmarks{iBenchmark, 1};
                    fcn1 = fcnpkgs{1};
                    fcn2 = fcnpkgs{2};
                    data = fcn2(nDim);
                    
                    objFunc = @(x,y)(fcn1(x, y, data));
                end
                
                [xOut, funVal, exitFlag, OUTPUT] = bsRunOneMethodOneSimulation( ...
                    objFunc, ...
                    methodName, ... 
                    parameters, ...
                    lower, ...
                    upper, ...
                    optimalFunctionTolerance, ...
                    maxIter, ...
                    minFVal, ...
                    maxFevs, ...
                    isSaveMiddleRes, ...
                    isSaveDetailUpdates);       
                
%                 nIter = OUTPUT.iterations;
                nfev = OUTPUT.funcCount;
%                 midResults = OUTPUT.midResults;
                
                fprintf('Test function:%s, dimension:%d, using method:%s, simulation:%d nfev:%d...\n', objFuncName, nDim, caseNames{iTestMethod}, i, nfev);
                
                detailResults{iBenchmark, iTestMethod, i}.output = OUTPUT;
                detailResults{iBenchmark, iTestMethod, i}.xout = xOut;
                detailResults{iBenchmark, iTestMethod, i}.funval = funVal;
                
                
            end
            
        end
        
    end
end

% save detail_updates;
if(size(testMethods, 1) > 1)   
    showOneAnimation;
else
    show3DAnimationProcess;
end
