classdef MultiplierDeltaCombined < MultiplierDelta
%% MULTIPLIERDELTACOMBINED for combining numerous multipliers prior to iqcAnalysis
%  Extends the base class MultiplierDelta
%
%    extended methods:
%       MultiplierDeltaCombined(mults_del) :: Constructor for combining deltas
%
%    extended properties:
%       objective : double or sdpvar :: objective 
%
%  See also MultiplierDeltaCombined.MultiplierDeltaCombined

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    objective
end

methods
function this_mult = MultiplierDeltaCombined(mults_del)
%% MULTIPLIERDELTACOMBINED constructor
%
%  this_mult = MultiplierDeltaCombined(mults)
%
%  Variables:
%  ---------
%    Input:
%       mults : array or singleton of MultiplierDelta objects
%    Output:
%       this_mult : MultiplierDeltaCombined object :: a single multiplier which combines
%                                                     all of the input multipliers
%
%  See also MultiplierDeltaCombined

    validateattributes(mults_del, 'MultiplierDelta', {'nonempty'})

    horizon_periods = vertcat(mults_del.horizon_period);
    matchingHorizonPeriods = @(hp) all(ismember(hp, hp(end,:), 'rows'));
    if ~isempty(horizon_periods)
        assert(matchingHorizonPeriods(horizon_periods),...
               'MultiplierDeltaCombined:MultiplierDeltaCombined',...
               'All input multipliers must have the same horizon_period')
    end

    % Combine filters
    filters = horzcat(mults_del.filter);
    a   = vertcat(filters.a);
    b1  = vertcat(filters.b1);
    b2  = vertcat(filters.b2);
    c1  = vertcat(filters.c1);
    c2  = vertcat(filters.c2);
    d11 = vertcat(filters.d11);
    d12 = vertcat(filters.d12);
    d21 = vertcat(filters.d21);
    d22 = vertcat(filters.d22);

    total_time = length(mults_del(end).filter.a);
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
        filter.b1{1,i}  = blkdiag(b1{:,i});
        filter.b2{1,i}  = blkdiag(b2{:,i});
        filter.c1{1,i}  = blkdiag(c1{:,i});
        filter.c2{1,i}  = blkdiag(c2{:,i});
        filter.d11{1,i} = blkdiag(d11{:,i});
        filter.d12{1,i} = blkdiag(d12{:,i});
        filter.d21{1,i} = blkdiag(d21{:,i});
        filter.d22{1,i} = blkdiag(d22{:,i});
    end

    % Combine quadratic terms
    quads = horzcat(mults_del.quad);
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

    name = {mults_del.name};
    decision_vars = horzcat(mults_del.decision_vars);
    constraints = horzcat(mults_del.constraints);
    horizon_period = mults_del(end).horizon_period;

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