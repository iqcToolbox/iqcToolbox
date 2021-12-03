function lft_out = rdivide(left_lft, right_lft)
%% ./ RDIVIDE overloaded function for dividing the left LFT with the right LFT.
%     lft_out = left_lft ./ right_lft
%     lft_out = rdivide(left_lft, right_lft)
%
%     Variables:
%     ---------
%       Input:
%         left_lft : Ulft object :: the lft on the left side of the "./" symbol
%         right_lft : Ulft object :: the lft on the right side of the "./" symbol
%       Output:
%         lft_out : Ulft object :: the resultant lft
%
%  This method differs from mrdivide only if right_lft is a scalar
%
%     See also Ulft, mrdivide, times, inv.

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
   
isScalar = @(lft) all(size(lft, 1) == 1) && all(size(lft, 2) == 1);
areConformable = @(llft, rlft) all(size(llft, 2) == size(rlft, 1));
assert(areConformable(left_lft, right_lft) || isScalar(right_lft),...
       'Ulft:rdivide',...
       ['Ulft objects do not have conformable dimensions, nor is the right lft',...
        ' scalar'])
% All other necessary properties are checking in mtimes and inv

%% Determine if doing typical mrdivide operation, or dividing one LFT with a scalar
if areConformable(left_lft, right_lft)
    warning('Ulft:rdivide',...
            ['./ is called, but it appears that / (or mrdivide) is necessary.',...
             ' Will use mrdivide instead'])
    lft_out = mrdivide(left_lft, right_lft);
elseif isScalar(right_lft)
    dim_in_left = size(left_lft, 2);
    assert(all(dim_in_left == dim_in_left(1)),...
           'Ulft:rdivide',...
           'Cannot perform ./ multiplication if one LFT has time-varying dims')
    lft_array = repmat({inv(right_lft)}, 1, dim_in_left(1));
    inv_right_lft = blkdiag(lft_array{:});
    lft_out = mtimes(left_lft, inv_right_lft);
end
end

%%  CHANGELOG
% Dec. 02, 2021: Added after v0.6.0 - Micah Fry (micah.fry@ll.mit.edu)