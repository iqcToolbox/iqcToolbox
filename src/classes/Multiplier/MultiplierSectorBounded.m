classdef MultiplierSectorBounded < MultiplierDelta
%% MUTLIPLIERSECTORBOUNDED class for sector-bounded (possibly nonlinear and/or
%  time-varying) uncertainties (DeltaSectorBounded). 
%  Extends the base class MultiplierDelta
%
%  extended methods:
%    MultiplierSectorBounded(delta, varargin) :: Constructor
%
%  extended properties:
%     upper_bound : double :: upper bound of uncertainty
%     lower_bound : double :: lower bound of uncertainty
%     dim_outin : double :: input/output dimensions of uncertainty
%     quad_time_varying : logical :: true if desired to define a
%                         time-varying quadratic decision variable
%
%  See also MultiplierSectorBounded.MultiplierSectorBounded

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
function this_mult = MultiplierSectorBounded(delta, varargin)
%% MULTIPLIERSECTORBOUNDED constructor
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
%      this_mult : MultiplierSectorBounded object

% Inputs
input_parser = inputParser;
addRequired(input_parser,...
            'delta',...
            @(del) validateattributes(del, {'DeltaSectorBounded'}, {'nonempty'}))
addParameter(input_parser,...
             'quad_time_varying',...
             true,...
             @(quad) validateattributes(quad, {'logical'}, {'nonempty'}))
addParameter(input_parser,... % This parameter is not used. Only defined for compatibility with other Multiplier constructor calls
             'discrete',...
             true,...
             @(disc) validateattributes(disc, {'logical'}, {'nonempty'}))
addParameter(input_parser,... % This parameter is not used. Only defined for compatibility with other Multiplier constructor calls
             'exponential',...
             [],...
             @(disc) validateattributes(disc, {'double'}, {'finite'})) 

% Parsing inputs
parse(input_parser, delta, varargin{:})
delta               = input_parser.Results.delta;
quad_time_varying   = input_parser.Results.quad_time_varying;

this_mult.name               = delta.name;
this_mult.horizon_period     = delta.horizon_period;
this_mult.dim_outin          = delta.dim_out;
this_mult.lower_bound        = delta.lower_bound;
this_mult.upper_bound        = delta.upper_bound;
this_mult.quad_time_varying  = quad_time_varying;

% Define filter
total_time = sum(this_mult.horizon_period);
filter.a     = cell(1, total_time);
filter.b1    = cell(1, total_time);
filter.b2    = cell(1, total_time);
filter.c1    = cell(1, total_time);
filter.c2    = cell(1, total_time);
filter.d11   = cell(1, total_time);
filter.d12   = cell(1, total_time);
filter.d21   = cell(1, total_time);
filter.d22   = cell(1, total_time);
quad.q11     = cell(1, total_time);
quad.q12     = cell(1, total_time);
quad.q21     = cell(1, total_time);
quad.q22     = cell(1, total_time);
if this_mult.quad_time_varying
    decision_vars = cell(1, total_time);
else
    decision_vars = cell(1);
end

constraints = [];
for i = 1:total_time
    % Define filter matrices
    filter.a{i}   = [];
    filter.b1{i}  = zeros(0, this_mult.dim_outin(i));
    filter.b2{i}  = zeros(0, this_mult.dim_outin(i));
    filter.c1{i}  = zeros(this_mult.dim_outin(i), 0);
    filter.c2{i}  = zeros(this_mult.dim_outin(i), 0);
    filter.d11{i} = -this_mult.lower_bound(i) * eye(this_mult.dim_outin(i));
    filter.d12{i} = eye(this_mult.dim_outin(i));
    filter.d21{i} = this_mult.upper_bound(i) * eye(this_mult.dim_outin(i));
    filter.d22{i} = -eye(this_mult.dim_outin(i));

    % Define quad
    if i == 1 || this_mult.quad_time_varying
        p = sdpvar(this_mult.dim_outin(i), 1);
        constraints = constraints + ((p >= 0):['Sector-bounded Multiplier, ',...
                                               this_mult.name,...
                                               ', p >= 0, #',...
                                               num2str(i)]);                    %#ok<BDSCA>
        decision_vars{i} = p;
    end
    quad.q12{i} = diag(p);
    quad.q21{i} = quad.q12{i}';
    quad.q11{i} = zeros(this_mult.dim_outin(i));
    quad.q22{i} = zeros(this_mult.dim_outin(i));    
end       
this_mult.filter        = filter;
this_mult.quad          = quad;
this_mult.constraints   = constraints;
this_mult.decision_vars = decision_vars;
end    
end

end

%%  CHANGELOG
% Dec. 1, 2021: Added after v0.6.0 - Micah Fry (micah.fry@ll.mit.edu)