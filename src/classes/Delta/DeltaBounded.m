classdef DeltaBounded < Delta
%% DELTABOUNDED class for uncertainties with a known l2-induced norm upper 
%  bound (may be nonlinear, dynamic, and/or time-varying). Extends the base 
%  class Delta.
%
%   extended methods:
%     DeltaBounded(name, dim_out, dim_in, upper_bound, horizon_period) :: Constructor
%     disp(this_delta) :: Display method
%     matchHorizonPeriod(this_delta, new_horizon_period) 
%                        :: Matches delta properties to new horizon_period
%     deltaToMultiplier(this_delta, varargin)
%                        :: Method for constructing a multiplier from a delta
%     normalizeDelta(this_delta) :: Changes unnormalized LFT to normalized one
%
%   extended properties:
%     upper_bound : double :: upper bound of uncertainty
%
%   See also Delta

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    upper_bound
end

methods
    function this_delta = DeltaBounded(name,...
                                       dim_out,...
                                       dim_in,...
                                       upper_bound,...
                                       horizon_period)
    %% DELTABOUNDED constructor
    %
    %  d = DeltaBounded(name, dim_out, dim_in, upper_bound, horizon_period)
    %  d = DeltaBounded(name, dim_out, dim_in, upper_bound) assumes horizon_period == [0, 1]
    %  d = DeltaBounded(name, dim_out, dim_in) also assumes upper_bound == 1
    %  d = DeltaBounded(name) also assumes dim_out == dim_in == 1
    %
    %  Variables:
    %  ---------
    %    Input:
    %       name : string :: unique ID of the uncertainty (ex. 'dforce')
    %       dim_out : natural :: output dimensions of uncertainty
    %       dim_in : natural :: input dimensions of uncertainty
    %       upper_bound : double :: upper_bound of uncertainty
    %       horizonperiod : 1 x 2 array of naturals :: [horizon, period] (in timesteps) 
    %                                                  of Delta properties
    %    Output:
    %       this_delta : DeltaBounded object :: the constructed DeltaBounded object
    %
    %  See also Delta, DeltaBounded, Delta.Delta

        % Defining defaults for missing arguments
        switch nargin
            case 1
                dim_out = 1;
                dim_in = 1;
                upper_bound = 1.0;
                horizon_period = [0, 1];
            case 3
                upper_bound = 1.0;
                horizon_period = [0, 1];
            case 4
                horizon_period = [0, 1];
        end
        % Calling Delta constructor
        this_delta@Delta(name, dim_out, dim_in, horizon_period);                    

        % Checking inputs for specialized properties of DeltaBounded
        validateattributes(upper_bound,...
                           'numeric',...
                           {'nonnan', 'finite', 'nonnegative'},...
                           mfilename);
        assert(all(upper_bound(1) == upper_bound),...
               'DeltaBounded:DeltaBounded',...
               'upper_bound must be constant for DeltaBounded objects')

        this_delta.upper_bound = upper_bound;
        this_delta = matchHorizonPeriod(this_delta);
    end

    function disp(this_delta)
    %% DISP function for DeltaBounded object
    %
    %  disp(delta_bounded_obj) (e.g., disp(DeltaBounded('d')) )
    %
    %  Variables:
    %  ---------
    %     Input:
    %        this_delta : DeltaBounded object
    %
    %  See also DeltaBounded, Delta.disp, Ulft.disp
    
        disp@Delta(this_delta, 'l2 norm-bounded uncertainty')  
        fprintf('%13s with an l2-induced norm less than: %3.1f \n',...
                '   ',...
                this_delta.upper_bound(1))            
    end

    function this_del = matchHorizonPeriod(this_del, new_horizon_period)
    %% MATCHHORIZONPERIOD function to ensure properties of DeltaBounded 
    %  object match its own horizon_period, or a new_horizon_period
    %
    %  this_del = matchHorizonPeriod(this_del, new_horizon_period) will change the horizon_period and associated properties of this_del
    %  this_del = matchHorizonPeriod(this_del) only ensures that each pertinent property of this_del matches this_del.horizon_period
    %
    %  Variables:
    %  ---------
    %     Input:
    %       this_del : DeltaBounded object
    %       new_horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of new sequence
    %    Output:
    %       this_delta : DeltaBounded object
    %
    %  See also DeltaBounded.
    
    if nargin == 1
    % Ensuring that this_del.horizon_period matches with other properties 
    % of this_del
    
        % Assumes properties are a mix of sequences of length 1 or
        % horizon_period
        total_time = sum(this_del.horizon_period);
        if length(this_del.dim_out) ~= total_time
            assert(length(this_del.dim_out) == 1,...
                   'DeltaBounded:matchHorizonPeriod',...
                   'dim_out of %s is not compatible w/ horizon_period',...
                   this_del.name);
            this_del.dim_out = this_del.dim_out * ones(1, total_time);
        end
        if length(this_del.dim_in) ~= total_time
            assert(length(this_del.dim_in) == 1,...
                   'DeltaBounded:matchHorizonPeriod',...
                   'dim_in of %s is not compatible w/ horizon_period',...
                   this_del.name)
            this_del.dim_in = this_del.dim_in * ones(1, total_time);
        end
        if length(this_del.upper_bound) ~= total_time
            assert(length(this_del.upper_bound) == 1,...
                   'DeltaBounded:matchHorizonPeriod',...
                   'upper_bound of %s is not compatible w/ horizon_period',...
                   this_del.name)
            this_del.upper_bound = ...
                this_del.upper_bound * ones(1, total_time);
        end
    else
    % Changing this_del.horizon_period and other properties of this_del to
    % a new horizon_period
        [indices, new_horizon_period] = ...
            makeNewIndices(this_del.horizon_period, new_horizon_period);
        
        % Set properties according to indices
        this_del.dim_out     = this_del.dim_out(indices);
        this_del.dim_in      = this_del.dim_in(indices);
        this_del.upper_bound = this_del.upper_bound(1) ...
                               * ones(1, sum(new_horizon_period));
        this_del.horizon_period = new_horizon_period;
        this_del = matchHorizonPeriod(this_del);
    end
    end

    function multiplier = deltaToMultiplier(this_del, varargin)                 %#ok<VANUS>
    %% DELTATOMULTIPLIER function to generate a multiplier from this object
    %
    %  multiplier = deltaToMultiplier(this_del)
    %
    %  Variables:
    %  ---------
    %     Input:
    %        this_del : DeltaBounded object
    %     Output:
    %        multiplier : MultiplierBounded object
    %
    %  See also DeltaBounded.
        multiplier = MultiplierBounded(this_del);
    end
    
    function [del_diff, del_ave, del_scale, del_norm] = normalizeDelta(this_delta)
    %% NORMALIZEDELTA function for DeltaBounded object. This extends the default operation
    %  defined by the Delta superclass, such that DeltaBounded uncertainties
    %  are normalized
    %
    %    [del_diff, del_ave, del_scale, del_norm] = normalizeDelta(this_delta)
    %
    %    Variables:
    %    ---------
    %      Input:
    %         this_delta : Delta object
    %      Output:
    %         del_diff : double matrix :: Matrix for modifying portions of LFT a, b, c, d
    %                                       matrices pertaining to this Delta
    %         del_ave : double matrix :: Matrix for modifying portions of LFT a, b, c, d
    %                                       matrices pertaining to this Delta
    %         del_scale : double matrix :: Matrix for modifying portions of LFT a, b, c, d
    %                                       matrices pertaining to this Delta
    %         del_norm : Delta object (or empty) :: Normalized Delta object. If empty, indicates
    %                                      that no normalization took place and a warning is thrown.
    %
    %    See also Ulft.normalizeLft, Delta.normalizeDelta
    
    % Define operations to modify a, b, c, d matrices of LFT
    total_time = sum(this_delta.horizon_period);
    del_diff   = cell(1, total_time);
    del_ave    = cell(1, total_time);
    del_scale  = cell(1, total_time);
    for i = 1:total_time
        del_diff{i}  = eye(this_delta.dim_out(i));
        del_ave{i}   = zeros(this_delta.dim_out(i), this_delta.dim_in(i));
        del_scale{i} = this_delta.upper_bound(i) * eye(this_delta.dim_in(i));
    end
    
    % Define normalized Delta
    del_norm = this_delta;
    del_norm.upper_bound = 1;
    del_norm = matchHorizonPeriod(del_norm, this_delta.horizon_period);
    end

    function value = sample(this_del, timestep)
        %% SAMPLE function for DeltaBounded.
        % The space of possible DeltaBounded is huge so we'll settle for these other deltas' subspaces of it
        if isempty(timestep)
            choice = randi([2,4]);
        else
            choice = randi([1,4]);
            if ~timestep
            assert(isequal(this_del.horizon_period, [0, 1]),...
                   'DeltaBounded:sample',...
                   ['Must have horizon_period = [0, 1] in order to ',...
                    'sample a continuous-time DeltaBounded'])
            end
        end
        switch choice
            case 1
                value = DeltaDlti(this_del.name,...
                                  this_del.dim_out(1),...
                                  this_del.dim_in(1),...
                                  this_del.upper_bound(1)).sample(timestep);
            case 2
                value = DeltaSlti(this_del.name,...
                                  this_del.dim_out(1),...
                                  -this_del.upper_bound(1),...
                                  this_del.upper_bound(1)).sample(timestep);
            case 3
                value = DeltaSltv(this_del.name,...
                                  this_del.dim_out,...
                                  -this_del.upper_bound,...
                                  this_del.upper_bound,...
                                  this_del.horizon_period).sample(timestep);
            case 4
                value = DeltaSltvRateBnd(this_del.name,...
                                         this_del.dim_out,...
                                         -this_del.upper_bound,...
                                         this_del.upper_bound,...
                                         -2*this_del.upper_bound*rand(),...
                                         2*this_del.upper_bound*rand(),...
                                         this_del.horizon_period).sample(timestep);
            otherwise
                error('DeltaBounded:DeltaBounded',...
                      'Switch-case mismatch (this is a backend bug)');
        end
        % Even if an LTI was sampled, the horizon_period should be matched back to the original
        % (unless continuous-time which currently only supports time-invariant)
        if isempty(timestep) || timestep(1)
            value = value.matchHorizonPeriod(this_del.horizon_period);
        end
        % Some of the other delta subspaces restrict to dim_out == dim_in, but we can generalize that
        if size(value, 2) ~= this_del.dim_in
            dim_fixer = rand(this_del.dim_out(1), this_del.dim_in(1)) - 0.5;
            dim_fixer = dim_fixer / (norm(dim_fixer, 2) + 1e-6);
            value = value * dim_fixer;
        end
    end

    function validateSample(this_delta, value, ~)
        %% VALIDATESAMPLE function for DeltaBounded.
        % Validate base attributes
        validateSample@Delta(this_delta, value, [], true);
        % Memoryless LFT has D-matrix tested for boundedness over time
        if isempty(value.delta.deltas)
            assert(...
                all(cellfun(@(d) norm(d,2) <= this_delta.upper_bound(1), value.d)),...
                'DeltaBounded:validateSample',...
                ['Specified value for delta "',this_delta.name,'" must have L2-induced norm within the bound.']);
        % Dynamic LFT gets analyzed for induced norm
        elseif isa(value.delta.deltas{1}, 'DeltaDelayZ') ||...
               isa(value.delta.deltas{1}, 'DeltaIntegrator')
            analysis = iqcAnalysis(value, 'analysis_options', AnalysisOptions('verbose', false));
            assert(...
                analysis.performance <= this_delta.upper_bound(1),...
                'DeltaBounded:validateSample',...
                ['Specified value for delta "',this_delta.name,'" must have L2-induced norm within the bound.']);
        end
    end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0): Added sample and validateSample methods - Jason Nezvadovitz (jason.nezvadovitz@ll.mit.edu)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)