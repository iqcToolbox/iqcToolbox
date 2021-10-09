classdef MultiplierPerformancePassive < MultiplierPerformance
%% MUTLIPLIERPASSIVE class. Used for the passive performance metric (PerformancePassive)
%  Extends the base class MultiplierPerformance
%
%  extended methods:
%    MultiplierPerformancePassive(performance, dim_out_lft, dim_in_lft) :: Constructor
%
%  extended properties:
%    chan_in : cell array of column vectors of naturals :: The input channels of the LFT pertaining 
%                                                          to this performance
%    chan_out : cell array of column vectors of naturals :: The output channels of the LFT pertaining 
%                                                          to this performance
%    dim_out : row of naturals :: The output dimension of the LFT pertaining to this performance
%    dim_in : row of naturals :: The input dimension of the LFT pertaining to this performance
%    
%  See also MultiplierPerformancePassive.MultiplierPerformancePassive

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    chan_in (1, :) cell
    chan_out (1, :) cell
    dim_out (1, :) double
    dim_in (1, :) double
end

methods
    
function this_mult = MultiplierPerformancePassive(performance, dim_out_lft, dim_in_lft)
%% MULTIPLIERPERFORMANCEPASSIVE constructor
%
%  this_mult = MultiplierPerformancePassive(performance, dim_out_lft, dim_in_lft)
%
%  Variables:
%  ---------
%    Input:
%      performance : PerformancePassive object :: performance defining this multiplier
%      dim_out_lft : row of naturals :: output dimension of LFT associated with performance
%      dim_in_lft : row of naturals :: input dimension of LFT associated with performance
%                                           
%    Output:
%      this_mult : MultiplierPerformancePassive object

% Inputs
input_parser = inputParser;
addRequired(input_parser,...
            'performance',...
            @(perf) validateattributes(perf,...
                                       'PerformancePassive',...
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
parse(input_parser, performance, dim_out_lft, dim_in_lft)
performance       = input_parser.Results.performance;
dim_out_lft       = input_parser.Results.dim_out_lft;
dim_in_lft        = input_parser.Results.dim_in_lft;

total_time   = sum(performance.horizon_period);
for i = 1:total_time
    if ~isempty(performance.chan_out{i})
    assert(dim_out_lft(i) >= max(performance.chan_out{i}),...
           'MultiplierPerformancePassive:MultiplierPerformancePassive',...
           'Perfomance output channel exceeds the output dimensions of lft')
    end
    if ~isempty(performance.chan_in{i})
    assert(dim_in_lft(i) >= max(performance.chan_in{i}),...
           'MultiplierPerformancePassive:MultiplierPerformancePassive',...
           'Perfomance input channel exceeds the input dimensions of lft')
    end
end

this_mult.name              = performance.name;
this_mult.chan_in           = performance.chan_in;
this_mult.chan_out          = performance.chan_out;
this_mult.horizon_period    = performance.horizon_period;
this_mult.dim_out           = dim_out_lft;
this_mult.dim_in            = dim_in_lft;

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
        quad.q12{i} = -double(diag(ismember(1:this_mult.dim_out(i),...
                                            this_mult.chan_out{i}')));
    else
        quad.q12{i} = -eye(this_mult.dim_out(i));
    end
    quad.q21{i}   = quad.q12{i}';
    quad.q11{i}   = zeros(this_mult.dim_out(i), this_mult.dim_in(i));
    quad.q22{i}   = zeros(this_mult.dim_in(i), this_mult.dim_out(i));   
end       

this_mult.objective         = 0;
this_mult.filter            = filter;
this_mult.quad              = quad;
this_mult.objective_scaling = 1;
end    
end

end

%%  CHANGELOG
% Oct. 7, 2021: Added after v0.6.0 - Micah Fry (micah.fry@ll.mit.edu)