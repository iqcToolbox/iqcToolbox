function lft_out = minus(left_lft, right_lft)
%% - MINUS overloaded function for subtracting an Ulft object with another.
%
%     lft_out = left_lft - right_lft
%     lft_out = minus(left_lft, right_lft)
%
%     Variables:
%     ---------
%       Input:
%         left_lft : Ulft object :: the lft on the left side of the "-" symbol
%         right_lft : Ulft object :: the lft on the right side of the "-" symbol
%       Output:
%         lft_out : Ulft object :: the resultant lft
%
%     See also Ulft, uminus, uplus, plus, mtimes.

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
validateattributes(left_lft, 'Ulft', {'nonempty'}, mfilename)
validateattributes(right_lft, 'Ulft', {'nonempty'}, mfilename)

if ~isequal(left_lft.horizon_period, right_lft.horizon_period)
    horizon_period = commonHorizonPeriod([left_lft.horizon_period;
                                          right_lft.horizon_period]);
    left_lft = matchHorizonPeriod(left_lft, horizon_period);
    right_lft = matchHorizonPeriod(right_lft, horizon_period);
end
assert(all(left_lft.horizon_period == right_lft.horizon_period),...
       'Ulft:minus',...
       'Both Ulft objects must have the same horizon_period')
   
assert(all(size(left_lft, 1) == size(right_lft, 1)),...
       'Ulft:minus',...
       'Both Ulft objects must have the same number of output channels')
assert(all(size(left_lft, 2) == size(right_lft, 2)),...
       'Ulft:minus',...
       'Both Ulft objects must have the same number of input channels')
assert(isequaln(left_lft.performance, right_lft.performance),...
...assert(all(strcmp(left_lft.performance.names, right_lft.performance.names)),...
       'Ulft:minus',...
       'Both Ulft objects must have the same performances')
assert(isequaln(left_lft.disturbance, right_lft.disturbance),...
...assert(all(strcmp(left_lft.disturbance.names, right_lft.disturbance.names)),...
       'Ulft:minus',...
       'Both Ulft objects must have the same disturbances')

%% Conduct subtraction
lft_out = left_lft + (-right_lft);
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)