function [state,options,optchanged] = gamyoutputfcn(options,state,flag)
     optchanged = false;
%      global allMidResults;
     i = get(getCurrentTask(),'ID');

     switch flag
         case 'init'
             midResults = [0; 0; min(state.Score)];
         case 'iter'
             midResults = evalin('base', sprintf('midResults_%d', i));
             midResults = [midResults, [state.Generation; state.FunEval; min(state.Score)]];
         case 'done'
             midResults = evalin('base', sprintf('midResults_%d', i));
             if state.Generation > midResults(1, end)
                midResults = [midResults, [state.Generation; state.FunEval; min(state.Score)]];
             end
         otherwise
     end
% 
%     midResults
     assignin('base', sprintf('midResults_%d', i), midResults);
end