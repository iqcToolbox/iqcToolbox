function lft_out = mrdivide(left_lft, right_lft)
%% / MRDIVIDE overloaded function for dividing the left LFT with the right LFT.
%     lft_out = left_lft / right_lft
%     lft_out = mrdivide(left_lft, right_lft)
%
%     Variables:
%     ---------
%       Input:
%         left_lft : Ulft object :: the lft on the left side of the "/" symbol
%         right_lft : Ulft object :: the lft on the right side of the "/" symbol
%       Output:
%         lft_out : Ulft object :: the resultant lft
%
%  This method is no different than left_lft * inv(right_lft)
%
%     See also Ulft, rdivide, times, inv.

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
       'Ulft:mrdivide',...
       'Both Ulft objects must have the same horizon_period')
% All other necessary properties are checked in mtimes and inv

%% Form LFT from existing methods
lft_out = left_lft * inv(right_lft);                                            %#ok<MINV>
end

%%  CHANGELOG
% Dec. 02, 2021: Added after v0.6.0 - Micah Fry (micah.fry@ll.mit.edu)