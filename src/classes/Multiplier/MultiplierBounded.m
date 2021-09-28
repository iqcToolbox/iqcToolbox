classdef MultiplierBounded < MultiplierDelta
%% MULTIPLIERBOUNDED class for bounded uncertainties (DeltaBounded).
% Extends the base class MultiplierDelta.
%
%   extended methods:
%     MultiplierBounded(delta) :: Constructor
%
%   extended properties:
%     upper_bound : double :: upper bound of uncertainty
%
%  See also MultiplierBounded.MultiplierBounded

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    upper_bound
end

methods
function this_mult = MultiplierBounded(delta)
%% MULTIPLIERBOUNDED constructor
%
%  this_mult = MultiplierBounded(delta)
%
%  Variables:
%  ---------
%    Input:
%       delta : DeltaBounded object :: uncertainty defining this multiplier
%    Output:
%       this_mult : MultiplierBounded object
%
%  See also MultiplierBounded

    validateattributes(delta, 'DeltaBounded', {'nonempty'}, mfilename)

    this_mult.name              = delta.name;
    this_mult.horizon_period    = delta.horizon_period;
    this_mult.upper_bound       = delta.upper_bound;        

    % Define filter
    total_time = sum(this_mult.horizon_period);
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
        filter.a{i}   = [];
        filter.b1{i}  = zeros(0, delta.dim_in(i));
        filter.b2{i}  = zeros(0, delta.dim_out(i));
        filter.c1{i}  = zeros(delta.dim_in(i), 0);
        filter.c2{i}  = zeros(delta.dim_out(i), 0);
        filter.d11{i} = this_mult.upper_bound(1) * eye(delta.dim_in(i));
        filter.d12{i} = zeros(delta.dim_in(i), delta.dim_out(i));
        filter.d21{i} = zeros(delta.dim_out(i), delta.dim_in(i));
        filter.d22{i} = eye(delta.dim_out(i));
    end
    this_mult.filter = filter;

    p = sdpvar(1);
    this_mult.decision_vars = {p};
    this_mult.constraints = (p >= 0):['Bounded Multiplier, ',...
                                      this_mult.name,...
                                      ' #1, p >= 0'];                           %#ok<BDSCA>

    % Define quadratic
    quad.q11 = cell(1, total_time);
    quad.q12 = cell(1, total_time);
    quad.q21 = cell(1, total_time);
    quad.q22 = cell(1, total_time);
    for i = 1:total_time
        quad.q11{i} = p * eye(delta.dim_in(i));
        quad.q12{i} = zeros(delta.dim_in(i), delta.dim_out(i));
        quad.q21{i} = quad.q12{i}';
        quad.q22{i} = -p * eye(delta.dim_out(i));
    end
    this_mult.quad = quad;
end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)