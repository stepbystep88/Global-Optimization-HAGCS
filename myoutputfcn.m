function stop = myoutputfcn(optimValues,state)
     stop = false;
     i = get(getCurrentTask(),'ID');
     
%      disp(i)
     
     switch state
         case 'init'
             midResults = [0; 0; optimValues.bestfval];
         case 'iter'
             midResults = evalin('base', sprintf('midResults_%d', i));
             midResults = [midResults, [optimValues.iteration; optimValues.funccount; optimValues.bestfval]];
         case 'done'
             midResults = evalin('base', sprintf('midResults_%d', i));
             if optimValues.iteration > midResults(1, end)
                midResults = [midResults, [optimValues.iteration; optimValues.funccount; optimValues.bestfval]];
             end
         otherwise
     end
     
     assignin('base', 'midResults', midResults);
     assignin('base', sprintf('midResults_%d', i), midResults);
 end