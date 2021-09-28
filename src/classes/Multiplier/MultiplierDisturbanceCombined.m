classdef MultiplierDisturbanceCombined < MultiplierDisturbance
%% MULTIPLIERDISTURBANCECOMBINED for combining numerous disturbance multipliers
%  Extends the base class MultiplierDisturbance
%
%    extended methods:
%       MultiplierDisturbanceCombined(mults_dis) :: Constructor for combining deltas
%
%  See also MultiplierDisturbanceCombined.MultiplierDisturbanceCombined

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

methods
    
function this_mult = MultiplierDisturbanceCombined(mults_dis)
%% MULTIPLIERDISTURBANCECOMBINED constructor
%
%  this_mult = MultiplierDisturbanceCombined(mults_dis)
%
%  Variables:
%  ---------
%    Input:
%       mults : array or singleton of MultiplierDisturbance objects
%    Output:
%       this_mult : MultiplierDisturbanceCombined object :: a single multiplier which combines
%                                                     all of the input multipliers
%
%  See also MultiplierDisturbanceCombined

    validateattributes(mults_dis, 'MultiplierDisturbance', {'nonempty'})

    horizon_periods = vertcat(mults_dis.horizon_period);
    matchingHorizonPeriods = @(hp) all(ismember(hp, hp(end,:), 'rows'));
    assert(matchingHorizonPeriods(horizon_periods),...
           'MultiplierDisturbanceCombined:MultiplierDisturbanceCombined',...
           'All input multipliers must have the same horizon_period')
                                     
    % Combine filters
    filters = horzcat(mults_dis.filter);
    a = vertcat(filters.a);
    b = vertcat(filters.b);
    c = vertcat(filters.c);
    d = vertcat(filters.d);

    total_time = length(mults_dis(end).filter.a);
    filter.a = cell(1, total_time);
    filter.b = cell(1, total_time);
    filter.c = cell(1, total_time);
    filter.d = cell(1, total_time);
    for i = 1:total_time
        filter.a{1, i} = blkdiag(a{:, i});
        filter.b{1, i} = vertcat(b{:, i});
        filter.c{1, i} = blkdiag(c{:, i});
        filter.d{1, i} = vertcat(d{:, i});
    end

    % Combine quadratic terms
    quads = horzcat(mults_dis.quad);
    q = vertcat(quads.q);
    quad.q = cell(1, total_time);
    for i = 1:total_time
        quad.q{1, i} = blkdiag(q{:, i});
    end

    name = {mults_dis.name};
    decision_vars = horzcat(mults_dis.decision_vars);
    constraints = horzcat(mults_dis.constraints);
    horizon_period = mults_dis(end).horizon_period;

    this_mult.name           = name;
    this_mult.filter         = filter;
    this_mult.decision_vars  = decision_vars;
    this_mult.constraints    = constraints;
    this_mult.quad           = quad;
    this_mult.horizon_period = horizon_period;
end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)