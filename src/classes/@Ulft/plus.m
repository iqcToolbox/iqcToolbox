function lft_out = plus(left_lft, right_lft)
%% + PLUS overloaded function for adding two Ulft objects.
%
%     lft_out = left_lft + right_lft
%     lft_out = plus(left_lft, right_lft)
%
%     Variables:
%     ---------
%       Input:
%         left_lft : Ulft object :: the lft on the left side of the "+" symbol
%         right_lft : Ulft object :: the lft on the right side of the "+" symbol
%       Output:
%         lft_out : Ulft object :: the resultant lft
%
%     See also Ulft, uplus, minus, uminus, mtimes.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Check correctness and consistency of inputs
if ~isa(left_lft, 'Ulft')
    left_lft = toLft(left_lft); 
end
if ~isa(right_lft, 'Ulft')
    right_lft = toLft(right_lft); 
end
validateattributes(left_lft, {'Ulft'}, {'nonempty'}, mfilename)
validateattributes(right_lft, {'Ulft'}, {'nonempty'}, mfilename)

if ~isequal(left_lft.horizon_period, right_lft.horizon_period)
    horizon_period = commonHorizonPeriod([left_lft.horizon_period;
                                          right_lft.horizon_period]);
    left_lft = matchHorizonPeriod(left_lft, horizon_period);
    right_lft = matchHorizonPeriod(right_lft, horizon_period);
end
assert(all(left_lft.horizon_period == right_lft.horizon_period),...
       'Ulft:plus',...
       'Both Ulft objects must have the same horizon_period')
   
assert(all(size(left_lft, 1) == size(right_lft, 1)),...
       'Ulft:plus',...
       'Both Ulft objects must have the same number of output channels')
assert(all(size(left_lft, 2) == size(right_lft, 2)),...
       'Ulft:plus',...
       'Both Ulft objects must have the same number of input channels')
%assert(isequaln(left_lft.performance, right_lft.performance),...
assert(all(strcmp(left_lft.performance.names, right_lft.performance.names)),...
       'Ulft:plus',...
       'Both Ulft objects must have the same performances')
%assert(isequaln(left_lft.disturbance, right_lft.disturbance),...
assert(all(strcmp(left_lft.disturbance.names, right_lft.disturbance.names)),...
       'Ulft:plus',...
       'Both Ulft objects must have the same disturbances')

%% Form Delta, a, b, c, and d matrices
delta = SequenceDelta(horzcat(left_lft.delta.deltas, right_lft.delta.deltas));

total_time = sum(left_lft.horizon_period);
d = cell(1, total_time);
c = cell(1, total_time);
b = cell(1, total_time);
a = cell(1, total_time);
for i = 1:total_time
    d{i} = left_lft.d{i} + right_lft.d{i};
    c{i} = [left_lft.c{i}, right_lft.c{i}];
    b{i} = [left_lft.b{i}; right_lft.b{i}];
    a{i} = blkdiag(left_lft.a{i}, right_lft.a{i});
end

%% Construct new lft
lft_out = Ulft(a, b, c, d, delta, ...
           'horizon_period', left_lft.horizon_period,...
           'performance', left_lft.performance,...
           'disturbance', left_lft.disturbance);

end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)