function this_lft = addPerformance(this_lft, performance)
%% ADDPERFORMANCE method for appending new performances to an LFT.
%
%     this_lft = addPerformance(this_lft, {perf1, perf2, perf3})
%
%     Variables:
%     ---------
%       Input:
%         this_lft : Ulft object :: the lft to add the performances to
%         performance : cell array of performance objects :: the performance objects to be added
%       Output:
%         this_lft : Ulft object :: the lft with the added performances
%
%     See also Ulft, addDisturbance.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

parser = inputParser;
isCellOfPerformances = @(arg) all(cellfun(@(a) isa(a, 'Performance'), arg));
parser.addRequired('performance',...
                   @(arg) iscell(arg) && isCellOfPerformances(arg));
parser.parse(performance);
new_performance = [this_lft.performance.performances, performance];

this_lft = Ulft(this_lft.a, this_lft.b, this_lft.c, this_lft.d, this_lft.delta,...
                'horizon_period', this_lft.horizon_period,...
                'disturbance', this_lft.disturbance,...
                'performance', new_performance);
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Jason Nezvadovitz (jason.nezvadovitz@ll.mit.edu)