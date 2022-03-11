classdef MultiplierSltv < MultiplierDelta
%% MULTIPLIERSLTV class for static, linear, time-varying uncertainties 
%  (DeltaSltv).
%  Extends the base class MultiplierDelta.
%
%  extended methods:
%     MultiplierSltv(delta, varargin) :: Constructor
%
%  extended properties:
%     upper_bound : double :: upper bound of uncertainty
%     lower_bound : double :: lower bound of uncertainty
%     dim_outin : double :: input/output dimensions of uncertainty
%     quad_time_varying : logical :: true if desired to define a
%                         time-varying quadratic decision variable
%
%  See also MultiplierSltv.MultiplierSltv

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    upper_bound double
    lower_bound double
    dim_outin double
end

properties (SetAccess = immutable)
    quad_time_varying logical
end

methods
function this_mult = MultiplierSltv(delta, varargin)
%% MULTIPLIERSLTV constructor
%
%  this_mult = MultiplierSltv(delta, 'quad_time_varying', true)
%  this_mult = MutliplierSltv(delta) assumes the values provided above
%
%  Variables:
%  ---------
%    Input:
%       delta : The delta object characterized by this multiplier
%       quad_time_varying : logical :: Whether or not the decision variable
%                            quad is time-varying (if true, SDP will possibly
%                            be less conservative, but more computationally
%                            challenging)
%    Output:
%       this_mult : MultiplierSltv object
%
%  See also MultiplierSltv

% Inputs
input_parser = inputParser;
addRequired(input_parser,...
            'delta',...
            @(del) validateattributes(del, {'DeltaSltv'}, {'nonempty'}))
addParameter(input_parser,...
             'quad_time_varying',...
             true,...
             @(quad) validateattributes(quad, {'logical'}, {'nonempty'}))
addParameter(input_parser,... % This parameter is not used. Only defined for compatibility with other Multiplier constructor calls
             'discrete',...
             true,...
             @(disc) validateattributes(disc, {'logical'}, {'nonempty'}))

% Parsing inputs
parse(input_parser, delta, varargin{:})
delta               = input_parser.Results.delta;
quad_time_varying   = input_parser.Results.quad_time_varying;

assert(all(abs(delta.upper_bound + delta.lower_bound) < 1e-8, 'all'),...
       'MultiplierSltv:MultiplierSltv',...
       'MultiplierSltv currently does not support asymmetric bounds',...
       ' on DeltaSltv (lower_bound must equal -upper_bound). Normalize')

this_mult.name               = delta.name;
this_mult.horizon_period     = delta.horizon_period;
this_mult.upper_bound        = delta.upper_bound;  
this_mult.lower_bound        = delta.lower_bound;  
this_mult.dim_outin          = delta.dim_out;
this_mult.quad_time_varying  = quad_time_varying;

% Define filter
total_time      = sum(this_mult.horizon_period);
filter.a        = cell(1, total_time);
filter.b1       = cell(1, total_time);
filter.b2       = cell(1, total_time);
filter.c1       = cell(1, total_time);
filter.c2       = cell(1, total_time);
filter.d11      = cell(1, total_time);
filter.d12      = cell(1, total_time);
filter.d21      = cell(1, total_time);
filter.d22      = cell(1, total_time);
for i = 1:total_time
    filter.a{i}   = [];
    filter.b1{i}  = zeros(0, delta.dim_in(i));
    filter.b2{i}  = zeros(0, delta.dim_out(i));
    filter.c1{i}  = zeros(delta.dim_in(i), 0);
    filter.c2{i}  = zeros(delta.dim_out(i), 0);
    filter.d11{i} = this_mult.upper_bound(i) * eye(delta.dim_in(i));
    filter.d12{i} = zeros(delta.dim_in(i), delta.dim_out(i));
    filter.d21{i} = zeros(delta.dim_out(i), delta.dim_in(i));
    filter.d22{i} = eye(delta.dim_out(i));
end
this_mult.filter = filter;

% Define quadratic, constraints, and decision_vars
quad.q11 = cell(1, total_time);
quad.q12 = cell(1, total_time);
quad.q21 = cell(1, total_time);
quad.q22 = cell(1, total_time);
if this_mult.quad_time_varying
    q11_length = total_time;
else
    q11_length = 1;
end
q11_var = cell(1, q11_length);
q12_var = cell(1, q11_length);
ct = [];
for i = 1:total_time
    if i == 1 || this_mult.quad_time_varying
        q11 = sdpvar(this_mult.dim_outin(i));
        q12 = sdpvar(size(q11, 1), size(q11, 2), 'skew');
        ct = ct + ((q11 >= 0):['SLTV Multiplier, ',...
                               this_mult.name,...
                               ', q11 >= 0, #',...
                               num2str(i)]);                                    %#ok<BDSCA>
        q11_var{i} = q11;
        q12_var{i} = q12;        
    end
    quad.q11{i} = q11;
    quad.q12{i} = q12;
    quad.q21{i} = q12';
    quad.q22{i} = -q11;
end
this_mult.quad = quad;
this_mult.constraints = ct;
this_mult.decision_vars = reshape([q11_var; q12_var], 1, 2 * q11_length);
end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)