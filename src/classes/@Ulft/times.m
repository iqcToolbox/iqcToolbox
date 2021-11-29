function lft_out = times(left_lft, right_lft)
%% .* TIMES overloaded function for multiplying two Ulft objects.
%     lft_out = left_lft .* right_lft
%     lft_out = times(left_lft, right_lft)
%
%     Variables:
%     ---------
%       Input:
%         left_lft : Ulft object :: the lft on the left side of the ".*" symbol
%         right_lft : Ulft object :: the lft on the right side of the ".*" symbol
%       Output:
%         lft_out : Ulft object :: the resultant lft
%
%     See also Ulft, mtimes.

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
       'Ulft:times',...
       'Both Ulft objects must have the same horizon_period')
   
defaultPerformance = @(lft) length(lft.performance.names) == 1 ...
                            && strcmp(lft.performance.names{1}, 'default_l2');
assert(defaultPerformance(right_lft) && defaultPerformance(left_lft),...
       'Ulft:times',...
       'Both Ulft objects must have a default performance measure')
defaultDisturbance = @(lft) length(lft.disturbance.names) == 1 ...
                            && strcmp(lft.disturbance.names{1}, 'default_l2');
assert(defaultDisturbance(left_lft),...
       'Ulft:times',...
       'The left Ulft object cannot have a non-default disturbance spec')
   
isScalar = @(lft) all(size(lft, 1) == 1) && all(size(lft, 2) == 1);
areConformable = @(llft, rlft) all(size(llft, 2) == size(rlft, 1));
assert(areConformable(left_lft, right_lft) || ...
       isScalar(left_lft) || ...
       isScalar(right_lft), ...
       'Ulft:times',...
       'Ulft objects do not have conformable dimensions, nor are either scalar') 

%% Determine if doing typical mtimes operation, or multiplying one LFT with a scalar
if areConformable(left_lft, right_lft)
    warning('Ulft:times',...
            ['.* is called, but it appears that * (or mtimes) is necessary.',...
             ' Will use mtimes instead'])
    lft_out = mtimes(left_lft, right_lft);
elseif isScalar(left_lft)
    dim_out_right = size(right_lft, 1);
    assert(all(dim_out_right == dim_out_right(1)),...
           'Ulft:times',...
           'Cannot perform .* multiplication if one LFT has time-varying dims')
    lft_array = repmat({left_lft}, 1, dim_out_right(1));
    left_lft = blkdiag(lft_array{:});
    lft_out = mtimes(left_lft, right_lft);
elseif isScalar(right_lft)
    dim_in_left = size(left_lft, 2);
    assert(all(dim_in_left == dim_in_left(1)),...
           'Ulft:times',...
           'Cannot perform .* multiplication if one LFT has time-varying dims')
    lft_array = repmat({right_lft}, 1, dim_in_left(1));
    right_lft = blkdiag(lft_array{:});
    lft_out = mtimes(left_lft, right_lft);
end
end

%%  CHANGELOG
% Nov. 29, 2021: Added after v0.6.0 - Micah Fry (micah.fry@ll.mit.edu)