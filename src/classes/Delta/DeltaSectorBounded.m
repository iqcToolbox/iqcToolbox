classdef DeltaSectorBounded < Delta
%% DELTASECTORBOUNDED class for operators (possibly nonlinear and/or 
%  time-varying) which are sector-bounded between upper and lower bounds. 
%  Extends the base class Delta.
%
%   extended methods:
%     DeltaSectorBounded(name, dim_outin, lower_bound, upper_bound, horizon_period) :: Constructor
%     disp(this_delta) :: Display method
%     matchHorizonPeriod(this_delta, new_horizon_period) 
%                        :: Matches delta properties to new horizon_period
%     deltaToMultiplier(this_delta, varargin)
%                        :: Method for constructing a multiplier from a delta
%
%   Note: Passive nonlinearities (see DeltaPassive) are sector-bounded within
%     the sector [0, Inf].
%
%   See also Delta, DeltaSectorBounded.DeltaSectorBounded, DeltaPassive

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    upper_bound (1, :) double {mustBeNonnegative, mustBeFinite}
    lower_bound (1, :) double {mustBeNonpositive, mustBeFinite}
end

methods
    function this_delta = DeltaSectorBounded(name, dim_outin, lower_bound, upper_bound, horizon_period)
    %% DELTASECTORBOUNDED constructor
    %
    %  d = DeltaSectorBounded(name, dim_outin, horizon_period)
    %  d = DeltaSectorBounded(name, dim_outin) assumes horizon_period == [0, 1]
    %  d = DeltaSectorBounded(name) also assumes dim_outin == 1
    %
    %  Variables:
    %  ---------
    %    Input:
    %       name : string :: unique ID of the uncertainty (ex. 'dforce')
    %       dim_outin : natural :: input/output dimensions of uncertainty
    %       horizonperiod : 1 x 2 array of naturals :: [horizon, period] (in timesteps) 
    %                                                  of Delta properties
    %    Output:
    %       this_delta : DeltaSectorBounded object :: the constructed DeltaSectorBounded object
    %
    %  See also Delta, DeltaSectorBounded, Delta.Delta

        % Defining defaults for missing arguments
        switch nargin
            case 1
                dim_outin = 1;
                lower_bound = -1;
                upper_bound = 1;
                horizon_period = [0, 1];
            case 2
                lower_bound = -1;
                upper_bound = 1;
                horizon_period = [0, 1];
            case 4
                horizon_period = [0, 1];
            case 5
            otherwise
                error('DeltaSectorBounded:DeltaSectorBounded',...
                      ['Must provide 1, 2, 4, or 5 arguments to construct',...
                       'DeltaSectorBounded objects'])
        end
        % Calling Delta constructor
        this_delta@Delta(name, dim_outin, dim_outin, horizon_period);
        
%         % Checking inputs for specialized properties of DeltaSectorBounded
%         validateattributes(lower_bound,...
%                            {'numeric'},...
%                            {'nonnan', 'finite', 'nonempty', 'nonpositive'})
%         validateattributes(upper_bound,...
%                            {'numeric'},...
%                            {'nonnan', 'finite', 'nonempty', 'nonnegative'})
        
        this_delta.lower_bound = lower_bound;
        this_delta.upper_bound = upper_bound;
                       
        this_delta = matchHorizonPeriod(this_delta);
    end

    function disp(this_delta)
    %% DISP function for DeltaSectorBounded object
    %
    %  disp(delta_sb_obj) (e.g., disp(DeltaSectorBounded('d')) )
    %
    %  Variables:
    %  ---------
    %     Input:
    %        this_delta : DeltaSectorBounded object
    %
    %  See also DeltaSectorBounded, Delta.disp, Ulft.disp
    
        disp@Delta(this_delta, 'sector bounded uncertainty')
        fprintf('%13s wherein the sector is given by: [%3.1f, %3.1f] \n',...
                '',...
                this_delta.lower_bound(1),...
                this_delta.upper_bound(1))  
    end

    function this_del = matchHorizonPeriod(this_del, new_horizon_period)
    %% MATCHHORIZONPERIOD function to ensure properties of DeltaSectorBounded 
    %  object match its own horizon_period, or a new_horizon_period
    %
    %  this_del = matchHorizonPeriod(this_del, new_horizon_period) will change the horizon_period and associated properties of this_del
    %  this_del = matchHorizonPeriod(this_del) only ensures that each pertinent property of this_del matches this_del.horizon_period
    %
    %  Variables:
    %  ---------
    %     Input:
    %       this_del : DeltaSectorBounded object
    %       new_horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of new sequence
    %    Output:
    %       this_delta : DeltaSectorBounded object
    %
    %  See also DeltaSectorBounded.
    
    if nargin == 1
    % Ensuring that this_del.horizon_period matches with other properties 
    % of this_del
    
        % Assumes properties are a mix of sequences of length 1 or
        % horizon_period
        total_time = sum(this_del.horizon_period);
        if length(this_del.dim_out) ~= total_time
            assert(length(this_del.dim_out) == 1,...
                   'DeltaSectorBounded:matchHorizonPeriod',...
                   'dim_out of %s is not compatible w/ horizon_period',...
                   this_del.name);
            this_del.dim_out = this_del.dim_out * ones(1, total_time);
        end
        if length(this_del.dim_in) ~= total_time
            assert(length(this_del.dim_in) == 1,...
                   'DeltaSectorBounded:matchHorizonPeriod',...
                   'dim_in of %s is not compatible w/ horizon_period',...
                   this_del.name)
            this_del.dim_in = this_del.dim_in * ones(1, total_time);
        end
        if length(this_del.lower_bound) ~= total_time
            assert(length(this_del.lower_bound) == 1,...
                   'DeltaSectorBounded:matchHorizonPeriod',...
                   'lower_bound of %s is not compatible w/ horizon_period',...
                   this_del.name)
            this_del.lower_bound = ...
                this_del.lower_bound * ones(1, total_time);
        end
        if length(this_del.upper_bound) ~= total_time
            assert(length(this_del.upper_bound) == 1,...
                   'DeltaSectorBounded:matchHorizonPeriod',...
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
        this_del.dim_out        = this_del.dim_out(indices);
        this_del.dim_in         = this_del.dim_in(indices);
        this_del.upper_bound    = this_del.upper_bound(indices);
        this_del.lower_bound    = this_del.lower_bound(indices);
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
    %        this_del : DeltaSectorBounded object
    %     Output:
    %        multiplier : MultiplierSectorBounded object
    %
    %  See also DeltaSectorBounded.
        multiplier = MultiplierSectorBounded(this_del, varargin{:});
    end
    
    function value = sample(this_del, ~)
    %% SAMPLE function for DeltaSectorBounded.
        total_time = sum(this_del.horizon_period);
        a = cell(1, total_time);
        b = cell(1, total_time);
        c = cell(1, total_time);
        d = cell(1, total_time);
        for t = 1 : total_time
            a{t} = zeros(0,0);
            b{t} = zeros(0,this_del.dim_in(t));
            c{t} = zeros(this_del.dim_out(t),0);
            d_mag = (this_del.upper_bound(t) - this_del.lower_bound(t))*rand()...
                    + this_del.lower_bound(t);
            d{t} = d_mag * eye(this_del.dim_out(t));
        end
        value = Ulft(a, b, c, d, {}, 'horizon_period', this_del.horizon_period);
    end

    function validateSample(this_delta, value, ~)
        %% VALIDATESAMPLE function for DeltaSectorBounded.
    % Validate base attributes
    validateSample@Delta(this_delta, value, [], true);
    % Memoryless
    assert(isempty(value.delta.deltas),...
           'DeltaSectorBounded:validateSample',...
           ['Specified value for delta "',this_delta.name,...
            '" must be memoryless.']);
    % Isotropic
    assert(all(cellfun(@(d) isdiag(d) && all(diag(d) == d(1,1)), value.d)),...
           'DeltaSectorBounded:validateSample',...
           ['Specified value for delta "',this_delta.name,...
            '" must be a scaling of the identity at all times.']);
    % Bounded
    lb = num2cell(this_delta.lower_bound);
    ub = num2cell(this_delta.upper_bound);
    assert(all(cellfun(@(d, lb, ub) (lb <= d(1,1)) && (d(1,1) <= ub),...
                       value.d, lb, ub)),...
           'DeltaSectorBounded:validateSample',...
           ['Specified value for delta "',this_delta.name,...
            '" must be within the bounds at all times.']);
    end
end
end

%%  CHANGELOG
% Dec. 1, 2021: Added after v0.6.0 - Micah Fry (micah.fry@ll.mit.edu)