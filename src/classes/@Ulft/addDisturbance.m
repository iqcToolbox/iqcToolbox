function this_lft = addDisturbance(this_lft, disturbance)
%% ADDDISTURBANCE method for appending new disturbances to an LFT.
%
%     this_lft = addDisturbance(this_lft, {dist1, dist2, dist3})
%
%     Variables:
%     ---------
%       Input:
%         this_lft : Ulft object :: the lft to add the disturbances to
%         disturbance : cell array of disturbance objects :: the disturbance objects to be added
%       Output:
%         this_lft : Ulft object :: the lft with the added disturbances
%
%     See also Ulft, addPerformance.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

parser = inputParser;
isCellOfDisturbances = @(arg) all(cellfun(@(a) isa(a, 'Disturbance'), arg));
parser.addRequired('disturbance',...
                   @(arg) iscell(arg) && isCellOfDisturbances(arg));
parser.parse(disturbance);
new_disturbance = [this_lft.disturbance.disturbances, disturbance];

this_lft = Ulft(this_lft.a, this_lft.b, this_lft.c, this_lft.d, this_lft.delta,...
                'horizon_period', this_lft.horizon_period,...
                'disturbance', new_disturbance,...
                'performance', this_lft.performance);
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Jason Nezvadovitz (jason.nezvadovitz@ll.mit.edu)