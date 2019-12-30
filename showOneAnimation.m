hl = figure;
lw = 0.96;
lh = 0.94;
mdw = 0.05;
mdh = 0.1;
left = 0.05;
below = 0.03;
set(gcf, 'position', [207   136   1054   700]);

%% prepare data
infos = cell(1, nTestMethods);
for iMethod = 1 : nTestMethods
    
    OUTPUT = detailResults{iMethod}.output;
    detailUpdate = OUTPUT.detailUpdates;
    nNest = length(detailUpdate);
    iters = zeros(1, nNest);
    for i = 1 : nNest
        iters(i) = size(detailUpdate{i}, 2);
    end
    maxIter = max(iters);
    infos{iMethod}.detailUpdate = detailUpdate;
    infos{iMethod}.maxIter = maxIter;
    infos{iMethod}.nNest = nNest;
    
    frames = moviein(maxIter);

    % plot the base contour map
    range = benchmarks{iBenchmark, 5};

    x = linspace(range(1), range(2), 500);
    y = x;
    nx = length(x);
    ny = nx;

    
    [X, Y] = meshgrid(x, y);

    if isa( benchmarks{iBenchmark, 1}, 'function_handle')
        objFunc = benchmarks{iBenchmark, 1};
    else
        fcnpkgs = benchmarks{iBenchmark, 1};
        fcn1 = fcnpkgs{1};
        fcn2 = fcnpkgs{2};
        o = fcn2(nDim);

        objFunc = @(x,y)(fcn1(x, y, o));
    end

    
    input = [reshape(X, 1, nx*ny); reshape(Y, 1, nx*ny)];
    output = objFunc(input, 0);
    Z = reshape(output, nx, ny);

    % plot the meshgrid
    bsSubPlotFit( nTestMethods, 3, 1+(iMethod-1)*3, lw, lh, mdw, mdh, left, below)
    mesh(X, Y, -Z, 'FaceAlpha',0.5);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title(sprintf('%s function', objFuncName));
    bsSetDefaultPlotSet(bsGetDefaultPlotSet());
    set(gca, 'view', [-22, 40]);
%     set(gca, 'view', [185, 50]);
    infos{iMethod}.view = get(gca, 'view');
    infos{iMethod}.reverse = false;
%     set(gca, 'zdir', 'reverse');
    hold on;
    s = scatter3(bestX(1), bestX(2), -minFVal, 'rp', 'filled','SizeData',200);
    alpha(s,.7)
%     plot3(bestX(1), bestX(2), -minFVal, 'rp', 'linewidth', 4, '');

    
    % plot the contour map
    bsSubPlotFit( nTestMethods, 3, 2+(iMethod-1)*3, lw, lh, mdw, mdh, left, below)
    if benchmarks{iBenchmark, 7}
        contour(X,Y,log(Z)); hold on;
    else
        contour(X,Y,Z); hold on;
    end
    xlabel('x');
    ylabel('y');
    
    title(sprintf('%s with population=%d', testMethods{iMethod, 1}, nNest));
    % plot global minimum;
    plot(bestX(1), bestX(2), 'rp', 'linewidth', 4);
    bsSetDefaultPlotSet(bsGetDefaultPlotSet());
    pls = [];

    %% obtain the range of the datas at each iteration
    maxVals = zeros(2, maxIter);
    minVals = zeros(2, maxIter);
    minFVals = zeros(1, maxIter);

    for i = 1 : maxIter
        maxX = -inf;
        maxY = -inf;
        minX = inf;
        minY = inf;
        minVal = inf;

        for j = 1 : nNest
            if i > size(detailUpdate{j}, 2)
                continue;
            end
            if maxX < detailUpdate{j}(3, i)
                maxX = detailUpdate{j}(3, i);
            end

            if maxY < detailUpdate{j}(4, i)
                maxY = detailUpdate{j}(4, i);
            end

            if minX > detailUpdate{j}(3, i)
                minX = detailUpdate{j}(3, i);
            end

            if minY > detailUpdate{j}(4, i)
                minY = detailUpdate{j}(4, i);
            end

            if minVal > detailUpdate{j}(2, i)
                minVal = detailUpdate{j}(2, i);
            end
        end
        maxVals(:, i) = [maxX; maxY];
        minVals(:, i) = [minX; minY];
        minFVals(i) = minVal;
    end
    
    infos{iMethod}.maxVals = maxVals;
    infos{iMethod}.minVals = minVals;
    infos{iMethod}.minFVals = minFVals;
end


%% show the animation of all methods
isExit = zeros(1, nTestMethods);
pls = cell(1, nTestMethods);
plotIds = zeros(1, nTestMethods);
maxIter = max(infos{1}.maxIter, infos{2}.maxIter);
for i = 1 : maxIter
    
    
        
    for iM = 1 : nTestMethods
        if isExit(iM)
            continue;
        end
        
        % clear the previous lines
        if exist('pls', 'var')
            delete(pls{iM});
        end
        pls{iM} = [];
    
        % plot the meshgrid
        bsSubPlotFit( nTestMethods, 3, 1+(iM-1)*3, lw, lh, mdw, mdh, left, below);
        points = zeros(infos{iM}.nNest, 3);
        
        for j = 1 : infos{iM}.nNest
            jniter = size(infos{iM}.detailUpdate{j}, 2);

            points(j, :) = infos{iM}.detailUpdate{j}(2:4, i);
        end
        
        tmp0 = plot3(points(:, 2), points(:, 3), -points(:, 1), 'b*', 'linewidth', 2); 
        
        pls{iM} = [pls{iM}, tmp0];
        
        
        l = min(points(:, 2));
        r = max(points(:, 2));
        dis = max(bestX(1) - l, r - bestX(1)) + 0.5;
%         set(gca, 'xlim', [bestX(1) - dis, bestX(1) + dis]);

        l = min(points(:, 3));
        r = max(points(:, 3));
        dis = max(bestX(2) - l, r - bestX(2)) + 0.5;
%         set(gca, 'ylim', [bestX(2) - dis, bestX(2) + dis]);

        l = -max(points(:, 1));
        r = -min(points(:, 1));
        dis = 1;
        set(gca, 'zlim', [l-0.01, -minFVal+0.01]);
        
        % set view
        set(gca, 'view', infos{iM}.view);
        
        
        if infos{iM}.view(2) < 3 
            infos{iM}.reverse = true;
        elseif  infos{iM}.view(2) > 85
            infos{iM}.reverse = false;
        end
            
       	if infos{iM}.reverse
            infos{iM}.view = infos{iM}.view - [1, -1];
        else
            infos{iM}.view = infos{iM}.view - [1, 1];
        end
        
        
        
        bsSubPlotFit( nTestMethods, 3, 2 + (iM-1)*3, lw, lh, mdw, mdh, left, below);
        
        
        % plot the contour map related convergence paths
        % set the title
        % title(sprintf('%s | Iter:%d | MinVal:%.4e', funcName, i, minFVals(i)));
    
        %     get the range of x
        if i >= 11
            l = min(infos{iM}.minVals(1, i-10:i));
            r = max(infos{iM}.maxVals(1, i-10:i));
            dis = max(bestX(1) - l, r - bestX(1)) + 0.001;
            set(gca, 'xlim', [bestX(1) - dis, bestX(1) + dis]);

            l = min(infos{iM}.minVals(2, i-10:i));
            r = max(infos{iM}.maxVals(2, i-10:i));
            dis = max(bestX(2) - l, r - bestX(2)) + 0.001;
            set(gca, 'ylim', [bestX(2) - dis, bestX(2) + dis]);
        end
        
        if infos{iM}.minFVals(i) - minFVal <= optimalFunctionTolerance
            isExit(iM) = true;
        end
            
        for j = 1 : infos{iM}.nNest
            jniter = size(infos{iM}.detailUpdate{j}, 2);

            s = i - 10;
            if s <= 0 
                s = 1;
            end

            e = i;
            if e > jniter
                e = size(infos{iM}.detailUpdate{j}, 2);
            end

            if s < e
                if e - s > 3
                    e;
                end
                nests = infos{iM}.detailUpdate{j}(:, s:e);
                
                if nests(2, end) == infos{iM}.minFVals(i) && isExit(iM)
                    tmp1 = plot(nests(3,:), nests(4,:), 'r-', 'linewidth', 2); hold on;
                else
                    tmp1 = plot(nests(3,:), nests(4,:), 'k-', 'linewidth', 1); hold on;
                end
                pls{iM} = [pls{iM}, tmp1];
            end


            if i > jniter
                continue;
            end

            nest = infos{iM}.detailUpdate{j}(:, i);

            if nest(1) == 1
                tmp2 = plot(nest(3), nest(4), 'rp', 'linewidth', 2); hold on;
                pause(0.1)
            else
                tmp2 = plot(nest(3), nest(4), 'b+', 'linewidth', 2); hold on;
            end
            pls{iM} = [pls{iM}, tmp2];

            
        end
        
        %% plot the convergence curves
        bsSubPlotFit( nTestMethods, 3, 3 + (iM-1)*3, lw, lh, mdw, mdh, left, below);
        if i == 1
            plotIds(iM) = plot(1, log10(infos{iM}.minFVals(1) - minFVal+1e-10), 'r-', 'linewidth', 2);
            xlabel('Iterations');
            ylabel('|F-F^*| (log scale)');
            title(sprintf('MinVal Plot | Iter:%d | MinVal:%.4e', i, infos{iM}.minFVals(i)));
            bsSetDefaultPlotSet(bsGetDefaultPlotSet());
            set(gca, 'ylim', [-10, 2]);
            set(gca, 'xlim', [1, maxIter]);
        else
            set(plotIds(iM), 'XData', 1:i);
            set(plotIds(iM), 'YData', log10(infos{iM}.minFVals(1:i) - minFVal + 1e-10));

            title(sprintf('MinVal Plot | Iter:%d | MinVal:%.4e', i, infos{iM}.minFVals(i)));
        end
    end
    
    
    %% draw and save the frame
    drawnow
    frames(:, i) = getframe(gcf);
%     pause(0.1)
end

bsSaveFramesToGif(frames, sprintf('./frames/Cmp_%s_%s_%s_%d.gif', testMethods{1, 1}, testMethods{2, 1}, objFuncName, nInitNest));

% plot(iter, fval, 'ks', 'linewidth', 2); hold on;
% drawnow
% 
% if iter == 1
%     xlabel('Iteration number');
%     ylabel('Objective function value');
%     bsSetDefaultPlotSet(bsGetDefaultPlotSet());
% end

