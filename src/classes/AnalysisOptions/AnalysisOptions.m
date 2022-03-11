classdef AnalysisOptions
%% ANALYSISOPTIONS class for setting IQC SDP options.
%
%   methods:
%     options = AnalysisOptions('solver', 'sdpt3', 'debug', true, ...)
%
%   properties:
%     yalmip_settings : struct :: struct for defining yalmip sdpsettings()
%     performance : double/empty :: the LFT's worst-case performance.  If
%                                   double, solve a feasibility problem. If
%                                   empty, solve an optimization problem.
%     lmi_shift : double :: shift parameter to express (>,<) constraints from (>=,=<)
%     scale_state_obj : double :: scale factor for initial state uncertainty metric,
%                                   scalarizing a vector optimization
%                                   problem when considering uncertain
%                                   initial conditions
%     init_cond_ellipse : double mat / empty / nan
%                               :: parameter for defining initial condition
%                                  ellipse.
%     init_cond_states  : logical row vector / empty
%                               :: row vector to express which states are
%                                  have non-zero initial conditions
%     p0 : double mat / empty :: provides initial lyapunov matrix to solve
%                                IQC KYP LMIs
%
%   See also AnalysisOptions.AnalysisOptions

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%
    
properties
    yalmip_settings
    lmi_shift
    scale_state_obj
    init_cond_ellipse
    init_cond_states
    p0
end

methods
    function this_options = AnalysisOptions(varargin)
    %% ANALYSISOPTIONS constructor
    %
    %       this_options = AnalysisOptions('solver', 'sdpt3', 'debug', true, ..., keyN, valueN)
    %
    %       Variables:
    %       ---------
    %         Keys: 
    %          'solver' : May be 'sdpt3', 'mosek', etc. Determines SDP solver.
    %          'verbose' : May be true/false. Determines if solver progress is displayed.
    %          'debug' : May be true/false. Determines YALMIP debug setting.
    %          'lmi_shift' : May be positive scalar. This determines how much a semidefinite constraint (>= or <=) 
    %                is shifted to express a definite constraint (> or <)
    %          'init_cond_ellipse' : May be empty, nan, or positive definite matrix. Expresses the ellipse 
    %               in which the initial condition may reside: {x in R^n | x' * E * x <= 1}
    %          'init_cond_states' : May be logical array. Indicates which indices of the full state are 
    %               uncertain and constrained by the aforementioned ellipse
    %          'scale_state_obj' : May be non-negative scalar. Indicates the amount the SDP objective is weighed to 
    %               reduce the impact of uncertainties in the initial condition (as opposed to reducing the value
    %               associated with a performance multiplier)
    %          'p0' : May be empty or symmetric matrix. Specifies the 0th indexed lyapunov matrix which satisfies
    %               the IQC analysis KYP LMIs. Typically used for debugging.
    
    input_parser = inputParser;

    addParameter(input_parser,...
                 'solver',...
                 'sdpt3',...
                 @(in) validateattributes(in, {'char'}, {'nonempty'}))
    addParameter(input_parser,...
                 'verbose',...
                 true,...
                 @(in) validateattributes(in, {'logical'}, {'nonempty'}))
    addParameter(input_parser,...
                 'debug',...
                 false,...
                 @(in) validateattributes(in, {'logical'}, {'nonempty'}))
    addParameter(input_parser,...
                 'lmi_shift',...
                 1e-8,...
                 @(in) validateattributes(in,...
                                         {'numeric'},...
                                         {'nonnegative', 'nonempty'}))
    addParameter(input_parser,...
                 'init_cond_ellipse',...
                 nan,...  % [nan] corresponds to x(0) = 0 (i.e., zero initial condition)
                 @(in) validateattributes(in, {'numeric'}, {'nonempty'}))
    % Other potential inputs:  [inf] corresponds to an undefined ellipse (to appear in optimization objective)
    %                          [matrix] corresponds to constraining ellipse for x(0) (must be positive definite)

    addParameter(input_parser,...
                 'init_cond_states',...
                 logical.empty(),... % [] corresponds to x(0) = 0 (i.e., zero initial condition)
                 @(in) validateattributes(in, {'logical'}, {'binary'}))
    % Other potential inputs: boolean row indicating which states are non-zero

    addParameter(input_parser,...
                 'scale_state_obj',...
                 1,...
                 @(in) validateattributes(in,...
                                          {'numeric'},...
                                          {'nonnegative', 'nonempty'}))
    addParameter(input_parser,...
                 'p0',...
                 [],...  % [] corresponds to allowing the solver find any suitable solution.
                 @(in) validateattributes(in, {'numeric'}, {'finite'}))                                 
    % Other potential inputs: [matrix] which is used as the first lyapunov matrix in an attempt to satisfy the IQC KYP LMIs

    parse(input_parser, varargin{:})
    
    yalmip_settings         = sdpsettings();
    yalmip_settings.solver  = input_parser.Results.solver;
    yalmip_settings.verbose = input_parser.Results.verbose;
    yalmip_settings.debug   = input_parser.Results.debug;

    this_options.yalmip_settings   = yalmip_settings;
    this_options.lmi_shift         = input_parser.Results.lmi_shift;
    this_options.scale_state_obj   = input_parser.Results.scale_state_obj;
    this_options.init_cond_ellipse = input_parser.Results.init_cond_ellipse;
    this_options.init_cond_states  = input_parser.Results.init_cond_states;
    this_options.p0                = input_parser.Results.p0;
    
    end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)