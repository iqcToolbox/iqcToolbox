classdef MultiplierL2Induced < MultiplierPerformance
%% MUTLIPLIERL2INDUCED class. Used as a placeholder for l2-induced performance metrics (PerformanceL2Induced)
%  Extends the base class MultiplierPerformance
%
%  extended methods:
%    MultiplierL2Induced(performance, dim_out_lft, dim_in_lft, varargin) :: Constructor
%
%  extended properties:
%    chan_in : cell array of column vectors of naturals :: The input channels of the LFT pertaining 
%                                                          to this performance
%    chan_out : cell array of column vectors of naturals :: The output channels of the LFT pertaining 
%                                                          to this performance
%    dim_out : row of naturals :: The output dimension of the LFT pertaining to this performance
%    dim_in : row of naturals :: The input dimension of the LFT pertaining to this performance
%    gain : sdpvar or double :: The l2-induced norm to be optimized or checked for feasibility
%    
%  See also MultiplierL2Induced.MultiplierL2Induced

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    chan_in (1, :) cell
    chan_out (1, :) cell
    dim_out (1, :) double
    dim_in (1, :) double
    gain 
end

methods
    
function this_mult = MultiplierL2Induced(performance, dim_out_lft, dim_in_lft, varargin)
%% MULTIPLIERL2INDUCED constructor
%
%  this_mult = MultiplierL2Induced(performance, dim_out_lft, dim_in_lft, 'objective_scaling', 1)
%  this_mult = MultiplierL2Induced(performance, dim_out_lft, dim_in_lft) assumes the values above 
%
%  Variables:
%  ---------
%    Input:
%      performance : PerformanceL2Induced object :: performance defining this multiplier
%      dim_out_lft : row of naturals :: output dimension of LFT associated with performance
%      dim_in_lft : row of naturals :: input dimension of LFT associated with performance
%      objective_scaling : double scalar :: weight on how much this_mult.gain should contribute
%                                           to the combined performances objective
%    Output:
%      this_mult : MultiplierL2Induced object

% Inputs
input_parser = inputParser;
addRequired(input_parser,...
            'performance',...
            @(perf) validateattributes(perf,...
                                       'PerformanceL2Induced',...
                                       {'nonempty'},...
                                       mfilename));
addRequired(input_parser,...
            'dim_out_lft',...
            @(dim) validateattributes(dim,...
                                      'numeric',...
                                      {'integer', 'row', 'positive'},...
                                      mfilename));
addRequired(input_parser,...
            'dim_in_lft',...
            @(dim) validateattributes(dim,...
                                      'numeric',...
                                      {'integer', 'row', 'positive'},...
                                      mfilename));
addParameter(input_parser,...
             'objective_scaling',...
             1,...
             @(os) validateattributes(os,...
                                      'numeric',...
                                      {'nonnegative', 'real', 'nonempty'}))
parse(input_parser, performance, dim_out_lft, dim_in_lft, varargin{:})
performance       = input_parser.Results.performance;
dim_out_lft       = input_parser.Results.dim_out_lft;
dim_in_lft        = input_parser.Results.dim_in_lft;
objective_scaling = input_parser.Results.objective_scaling;

total_time   = sum(performance.horizon_period);
for i = 1:total_time
    if ~isempty(performance.chan_out{i})
    assert(dim_out_lft(i) >= max(performance.chan_out{i}),...
           'MultiplierL2Induced:MultiplierL2Induced',...
           'Perfomance output channel exceeds the output dimensions of lft')
    end
    if ~isempty(performance.chan_in{i})
    assert(dim_in_lft(i) >= max(performance.chan_in{i}),...
           'MultiplierL2Induced:MultiplierL2Induced',...
           'Perfomance input channel exceeds the input dimensions of lft')
    end
end

this_mult.name              = performance.name;
this_mult.chan_in           = performance.chan_in;
this_mult.chan_out          = performance.chan_out;
this_mult.gain              = performance.gain;
this_mult.horizon_period    = performance.horizon_period;
this_mult.dim_out           = dim_out_lft;
this_mult.dim_in            = dim_in_lft;
this_mult.objective_scaling = objective_scaling;

% Define filter
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
if ~isempty(this_mult.gain)
    gain_squared = this_mult.gain^2;
    decision_vars = [];
else
    gain_squared = sdpvar;
    decision_vars{1} = gain_squared;
    this_mult.gain = sqrt(gain_squared);
end
for i = 1:total_time
    % Define filter matrices
    filter.a{i}   = [];
    filter.b1{i}  = zeros(0, this_mult.dim_out(i));
    filter.b2{i}  = zeros(0, this_mult.dim_in(i));
    filter.c1{i}  = zeros(this_mult.dim_out(i), 0);
    filter.c2{i}  = zeros(this_mult.dim_in(i), 0);
    filter.d11{i} = eye(this_mult.dim_out(i));
    filter.d12{i} = zeros(this_mult.dim_out(i), this_mult.dim_in(i));
    filter.d21{i} = zeros(this_mult.dim_in(i), this_mult.dim_out(i));
    filter.d22{i} = eye(this_mult.dim_in(i));

    % Define quad
    if ~isempty(this_mult.chan_out{i})
        quad.q11{i} = double(diag(ismember(1:this_mult.dim_out(i),...
                                           this_mult.chan_out{i}')));
    else
        quad.q11{i} = eye(this_mult.dim_out(i));
    end
    quad.q12{i}   = zeros(this_mult.dim_out(i), this_mult.dim_in(i));
    quad.q21{i}   = zeros(this_mult.dim_in(i), this_mult.dim_out(i));
    if ~isempty(this_mult.chan_in{i})
        quad.q22{i} = double(diag(ismember(1:this_mult.dim_in(i),...
                                           this_mult.chan_in{i}')));
    else
        quad.q22{i} = eye(this_mult.dim_in(i));
    end
    quad.q22{i}   = - gain_squared * quad.q22{i};
end       

this_mult.objective     = gain_squared;
this_mult.filter        = filter;
this_mult.quad          = quad;
this_mult.decision_vars = decision_vars;

end    
end

end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)