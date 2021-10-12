classdef MultiplierPerformanceCombined < MultiplierPerformance
%% MULTIPLIERPERFORMANCECOMBINED for combining numerous multipliers
%  Extends the base class MultiplierPerformance
%
%    extended methods:
%       MultiplierPerformanceCombined
%
%  See also MultiplierPerformanceCombined.MultiplierPerformanceCombined

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

methods
function this_mult = MultiplierPerformanceCombined(mults_perf)
%% MULTIPLIERPERFORMANCECOMBINED constructor
%
%  this_mult = MultiplierPerformanceCombined(mults_perf)
%
%  Variables:
%  ---------
%    Input:
%      mults : array or singleton of MultiplierPerformance objects
%    Output:
%       this_mult : MultiplierPerformanceCombined object :: a single multiplier which combines
%                                                     all of the input multipliers
%
%  See also MultiplierPerformanceCombined

    validateattributes(mults_perf, 'MultiplierPerformance', {'nonempty'})

    horizon_periods = vertcat(mults_perf.horizon_period);
    matchingHorizonPeriods = @(hp) all(ismember(hp, hp(end,:), 'rows'));
    if ~isempty(horizon_periods)
        assert(matchingHorizonPeriods(horizon_periods),...
               ['MultiplierPerformanceCombined:',...
                'MultiplierPerformanceCombined'],...
               'All input multipliers must have the same horizon_period')
    end

    % Combine filters
    filters = horzcat(mults_perf.filter);
    a   = vertcat(filters.a);
    b1  = vertcat(filters.b1);
    b2  = vertcat(filters.b2);
    c1  = vertcat(filters.c1);
    c2  = vertcat(filters.c2);
    d11 = vertcat(filters.d11);
    d12 = vertcat(filters.d12);
    d21 = vertcat(filters.d21);
    d22 = vertcat(filters.d22);

    total_time = length(mults_perf(end).filter.a);
    filter.a   = cell(1, total_time);
    filter.b1  = cell(1, total_time);
    filter.b2  = cell(1, total_time);
    filter.c1  = cell(1, total_time);
    filter.c2  = cell(1, total_time);
    filter.d11 = cell(1, total_time);
    filter.d12 = cell(1, total_time);
    filter.d21 = cell(1, total_time);
    filter.d22 = cell(1, total_time);
    for i = 1:total_time
        filter.a{1,i}   = blkdiag(a{:,i});
        filter.b1{1,i}  = vertcat(b1{:,i});
        filter.b2{1,i}  = vertcat(b2{:,i});
        filter.c1{1,i}  = blkdiag(c1{:,i});
        filter.c2{1,i}  = blkdiag(c2{:,i});
        filter.d11{1,i} = vertcat(d11{:,i});
        filter.d12{1,i} = vertcat(d12{:,i});
        filter.d21{1,i} = vertcat(d21{:,i});
        filter.d22{1,i} = vertcat(d22{:,i});
    end

    % Combine quadratic terms
    quads = horzcat(mults_perf.quad);
    q11 = vertcat(quads.q11);
    q12 = vertcat(quads.q12);
    q21 = vertcat(quads.q21);
    q22 = vertcat(quads.q22);
    quad.q11 = cell(1, total_time);
    quad.q12 = cell(1, total_time);
    quad.q21 = cell(1, total_time);
    quad.q22 = cell(1, total_time);
    for i = 1:total_time
        quad.q11{1,i} = blkdiag(q11{:,i});
        quad.q12{1,i} = blkdiag(q12{:,i});
        quad.q21{1,i} = blkdiag(q21{:,i});
        quad.q22{1,i} = blkdiag(q22{:,i});
    end

    % Combine objective values
    objective = 0;
    for i = 1:length(mults_perf)
        objective = objective + ...
                    mults_perf(i).objective * mults_perf(i).objective_scaling;
    end

    name = {mults_perf.name};
    decision_vars = horzcat(mults_perf.decision_vars);
    constraints = horzcat(mults_perf.constraints);
    horizon_period = mults_perf(end).horizon_period;
    discrete = horzcat(mults_perf.discrete);
    if ~isempty(discrete)
        assert(all(discrete(1) == discrete),...
               'MultiplierPerformanceCombined:MultiplierPerformanceCombined',...
               'Cannot combine Multipliers that are discrete and continuous')
        discrete = discrete(1);
    end

    this_mult.name           = name;
    this_mult.filter         = filter;
    this_mult.decision_vars  = decision_vars;
    this_mult.constraints    = constraints;
    this_mult.quad           = quad;
    this_mult.horizon_period = horizon_period;
    this_mult.objective      = objective;
    this_mult.discrete       = discrete;
end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)