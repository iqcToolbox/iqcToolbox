function this_lft = normalizeLft(this_lft)
%% NORMALIZELFT method for normalizing LFTs that have Deltas which are not centered
%     about zero or are not bounded by one. The output LFT has modified a,
%     b, c, and d matrices, as well as normalized Deltas, such that the
%     output LFT is a normalized representation of the unnormalized LFT.
%
%     normalized_lft = normalizeLft(this_lft)
%
%     Variables:
%     ---------
%       Input:
%         this_lft : Ulft object :: the lft whose Deltas are to be normalized
%       Output:
%         this_lft : Ulft object :: the lft with normalized Deltas
%
%     See also Ulft.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%
if isempty(this_lft.delta.deltas)
    % No Deltas to normalize, return input
    return;
end
num_delta = length(this_lft.delta.deltas);
total_time = sum(this_lft.horizon_period);

% Generate normalized delta and modifier operators to a, b, c, d matrices
del_diff  = cell(num_delta, total_time);
del_ave   = cell(num_delta, total_time);
del_scale = cell(num_delta, total_time);
del_norm  = cell(1, num_delta);
for i = 1:num_delta
    [del_d, del_a, del_s, del_n] = normalizeDelta(this_lft.delta.deltas{i});
    if isempty(del_n)
        warning('Ulft:normalizeLft',...
                ['Delta "%s" is not normalized, because it''s class does ',...
                 'not define a normalization method'],...
                this_lft.delta.names{i})
        del_n = this_lft.delta.deltas{i};
    end
    del_diff(i, :)  = del_d;
    del_ave(i, :)   = del_a;
    del_scale(i, :) = del_s;
    del_norm{i}     = del_n;    
end

% Modify a, b, c, d matrices
a_norm = cell(1, total_time);
b_norm = cell(1, total_time);
c_norm = cell(1, total_time);
d_norm = cell(1, total_time);
for j = 1:total_time
    del_diff_mat = blkdiag(del_diff{:, j});
    del_ave_mat = blkdiag(del_ave{:, j});
    del_scale_mat = blkdiag(del_scale{:, j});
    a  = this_lft.a{j};
    b  = this_lft.b{j};
    c  = this_lft.c{j};
    d  = this_lft.d{j};
    mat_inv = inv(eye(size(a, 1)) - a * del_ave_mat);
    a_norm{j} = del_scale_mat * mat_inv * a * del_diff_mat;
    b_norm{j} = del_scale_mat * mat_inv * b;
    c_norm{j} = c * (del_diff_mat + del_ave_mat * mat_inv * a * del_diff_mat);
    d_norm{j} = d + c * del_ave_mat * mat_inv * b;
end

this_lft = Ulft(a_norm, b_norm, c_norm, d_norm, del_norm,...
                'horizon_period', this_lft.horizon_period,...
                'performance', this_lft.performance,...
                'disturbance', this_lft.disturbance);
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)