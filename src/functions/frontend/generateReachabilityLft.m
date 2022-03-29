function lft_out = generateReachabilityLft(lft_in, final_time)
%% GENERATEREACHABILITYLFT function which creates a time-varying LFT such
%  that IQC analysis on the resulting LFT will provide an l2-to-euclidean
%  norm of the input LFT at the given "final_time" index
%
%  lft_out = generateReachabilityLft(lft_in, final_time)
%
%  Variables:
%  ---------
%     Input:
%       lft_in : Ulft object :: the original LFT who's l2-to-euclidean norm is sought
%       final_time : nonnegative integer :: time instant [0, inf) for l2-to-euclidean norm
%     Output:
%       lft_out : Ulft object :: an LFT of [final_time + 1, 1] horizon_period, representing
%                                a mapping from disturbances to the performance at the final_time
%
%  See also Ulft

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

validateattributes(lft_in, {'Ulft'}, {'nonempty'})
validateattributes(final_time, {'numeric'}, {'nonnegative', 'nonempty'})

% Find the right horizon_period
common_hp = commonHorizonPeriod([lft_in.horizon_period; [final_time + 1, 1]]);
[ind_reach, horizon_period] = makeNewIndices(lft_in.horizon_period,...
                                             common_hp);

% Match horizon_period of delta, disturbance, and performance                                     
new_delta = matchHorizonPeriod(lft_in.delta, horizon_period);
new_disturbance = matchHorizonPeriod(lft_in.disturbance, horizon_period);
new_performance = matchHorizonPeriod(lft_in.performance, horizon_period);

% Create state-space matrices of reachability LFT
total_time = sum(horizon_period);
a = lft_in.a(1, ind_reach);
b = lft_in.b(1, ind_reach);
c = lft_in.c(1, ind_reach);
d = lft_in.d(1, ind_reach);
for i = 1:total_time
    if i < final_time + 1
        c{i} = zeros(size(c{i}));
        d{i} = zeros(size(d{i}));
    elseif i > final_time + 1
        a{i} = zeros(size(a{i}));
        b{i} = zeros(size(b{i}));
        c{i} = zeros(size(c{i}));
        d{i} = zeros(size(d{i}));
    end
end

% Construct LFT
lft_out = Ulft(a, b, c, d, new_delta,...
               'horizon_period', horizon_period,...
               'disturbance', new_disturbance,...
               'performance', new_performance);       
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)