classdef DeltaConstantDelay < Delta
%% DELTACONSTANTDELAY class for constant delay operators whose delay is within
%  a range: delay \in [0, delay_max]. For discrete-time systems, this delay
%  is measured in timesteps (delay must be a natural number). For continuous-time
%  systems, delay is measured in seconds. 
%
%  This class differs from DeltaConstantDelay2 in that the operator Delta is 
%  defined as u_delay = Delta (u).
%  This formulation means that Delta's nominal value of 0 implies u_delay(t) = 0
%  for all t.
%
%   extended methods:
%     DeltaConstantDelay(name, dim_outin, delay_max, horizon_period) :: Constructor
%     disp(this_delta) :: Display method
%     matchHorizonPeriod(this_delta, new_horizon_period) 
%                        :: Matches delta properties to new horizon_period
%     deltaToMultiplier(this_delta, varargin)
%                        :: Method for constructing a multiplier from a delta
%     sample(this_delta, timestep) 
%                        :: Produces an LFT object sampled from the admissible set of operators
%     validateSample(this_delta, value, timestep) 
%                        :: Checks if the provided value is a valid element of the uncertainty set
%
%   See also Delta, DeltaConstantDelay.DeltaConstantDelay, DeltaConstantDelay2

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    delay_max (1, :) double {mustBeNonnegative, mustBeFinite}    
end

methods
    function this_delta = DeltaConstantDelay(name, dim_outin, delay_max, horizon_period)
    %% DELTACONSTANTDELAY constructor
    %
    %  d = DeltaConstantDelay(name, dim_outin, delay_max, horizon_period)
    %  d = DeltaConstantDelay(name, dim_outin, delay_max) assumes horizon_period == [0, 1]
    %  d = DeltaConstantDelay(name, dim_outin) also assumes delay_max == 1
    %  d = DeltaConstantDelay(name) also assumes dim_outin == 1
    %
    %  Variables:
    %  ---------
    %    Input:
    %       name : string :: unique ID of the uncertainty (ex. 'dforce')
    %       delay_max : nonnegative double :: maximum delay of Delta
    %       dim_outin : natural :: input/output dimensions of uncertainty
    %       horizonperiod : 1 x 2 array of naturals :: [horizon, period] (in timesteps) 
    %                                                  of Delta properties
    %    Output:
    %       this_delta : DeltaConstantDelay object :: the constructed DeltaConstantDelay object
    %
    %  See also Delta, DeltaConstantDelay, Delta.Delta

        % Defining defaults for missing arguments
        switch nargin
            case 1
                dim_outin = 1;
                delay_max = 1;
                horizon_period = [0, 1];
            case 2
                delay_max = 1;
                horizon_period = [0, 1];
            case 3
                horizon_period = [0, 1];
            case 4
            otherwise
                error('DeltaConstantDelay:DeltaConstantDelay',...
                      ['Must provide 1, 2, 3, or 4 arguments to construct',...
                       'DeltaConstantDelay objects'])
        end
        % Calling Delta constructor
        this_delta@Delta(name, dim_outin, dim_outin, horizon_period);
        
        % Checking inputs for specialized properties of DeltaConstantDelay
        validateattributes(delay_max,...
                           {'numeric'},...
                           {'nonnan', 'finite', 'nonempty', 'nonnegative'})
        
        this_delta.delay_max = delay_max;
                       
        this_delta = matchHorizonPeriod(this_delta);
    end

    function disp(this_delta)
    %% DISP function for DeltaConstantDelay object
    %
    %  disp(delta_sb_obj) (e.g., disp(DeltaConstantDelay('d')) )
    %
    %  Variables:
    %  ---------
    %     Input:
    %        this_delta : DeltaConstantDelay object
    %
    %  See also DeltaConstantDelay, Delta.disp, Ulft.disp
    
        disp@Delta(this_delta, 'constant delay uncertainty')
        fprintf('%13s formulated as: u_delay = delay(u) \n', '')
        fprintf('%13s wherein the delay may be within: [0, %3.1f] \n',...
                '',...
                this_delta.delay_max(1))  
    end

    function this_del = matchHorizonPeriod(this_del, new_horizon_period)
    %% MATCHHORIZONPERIOD function to ensure properties of DeltaConstantDelay
    %  object match its own horizon_period, or a new_horizon_period
    %
    %  this_del = matchHorizonPeriod(this_del, new_horizon_period) will change the horizon_period and associated properties of this_del
    %  this_del = matchHorizonPeriod(this_del) only ensures that each pertinent property of this_del matches this_del.horizon_period
    %
    %  Variables:
    %  ---------
    %     Input:
    %       this_del : DeltaConstantDelay object
    %       new_horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of new sequence
    %    Output:
    %       this_delta : DeltaConstantDelay object
    %
    %  See also DeltaConstantDelay.
    
    if nargin == 1
    % Ensuring that this_del.horizon_period matches with other properties 
    % of this_del
    
        % Assumes properties are a mix of sequences of length 1 or
        % horizon_period
        total_time = sum(this_del.horizon_period);
        if length(this_del.dim_out) ~= total_time
            assert(length(this_del.dim_out) == 1,...
                   'DeltaConstantDelay:matchHorizonPeriod',...
                   'dim_out of %s is not compatible w/ horizon_period',...
                   this_del.name);
            this_del.dim_out = this_del.dim_out * ones(1, total_time);
        end
        if length(this_del.dim_in) ~= total_time
            assert(length(this_del.dim_in) == 1,...
                   'DeltaConstantDelay:matchHorizonPeriod',...
                   'dim_in of %s is not compatible w/ horizon_period',...
                   this_del.name)
            this_del.dim_in = this_del.dim_in * ones(1, total_time);
        end
        if length(this_del.delay_max) ~= total_time
            assert(length(this_del.delay_max) == 1,...
                   'DeltaConstantDelay:matchHorizonPeriod',...
                   'delay_max of %s is not compatible w/ horizon_period',...
                   this_del.name)
            this_del.delay_max = ...
                this_del.delay_max * ones(1, total_time);
        end
    else
    % Changing this_del.horizon_period and other properties of this_del to
    % a new horizon_period
        [indices, new_horizon_period] = ...
            makeNewIndices(this_del.horizon_period, new_horizon_period);
        
        % Set properties according to indices
        this_del.dim_out        = this_del.dim_out(indices);
        this_del.dim_in         = this_del.dim_in(indices);
        this_del.delay_max      = this_del.delay_max(indices);
        this_del.horizon_period = new_horizon_period;
        
        this_del = matchHorizonPeriod(this_del);
    end
    end

    function multiplier = deltaToMultiplier(this_del, varargin)                 
    %% DELTATOMULTIPLIER function to generate a multiplier from this object
    %
    %  multiplier = deltaToMultiplier(this_del)
    %
    %  Variables:
    %  ---------
    %     Input:
    %        this_del : DeltaConstantDelay object
    %     Output:
    %        multiplier : MultiplierConstantDelay object
    %
    %  See also DeltaConstantDelay.
        multiplier = MultiplierConstantDelay(this_del, varargin{:});
    end
    
    function value = sample(this_delta, timestep)
        %% SAMPLE function for DeltaConstantDelay.
        % Check if dynamic and discrete-time
        assert(~isempty(timestep),...
               'DeltaConstantDelay:sample',...
               'LFT must be dynamic to sample a DeltaConstantDelay');
        assert(logical(timestep),...
               'DeltaConstantDelay:sample',...
               ['LFT must be discrete-time to create a state-space sample of',...
                ' DeltaConstantDelay'])
        % Generate discrete-time delay within bounds
        delay_steps = randi([1, this_delta.delay_max]);
        total_time = sum(this_delta.horizon_period);
        one_step = DeltaDelayZ(this_delta.dim_in(1) * ones(1, total_time),...
                               timestep,...
                               this_delta.horizon_period);
        value = one_step ^ delay_steps;
    end

    function validateSample(this_delta, value, timestep)
        %% VALIDATESAMPLE function for DeltaConstantDelay.
        % Validate base attributes
        validateSample@Delta(this_delta, value, timestep, true);
        % Time invariant
        assert(all(cellfun(@(a) isequal(a, value.a{1}), value.a)) &&...
               all(cellfun(@(b) isequal(b, value.b{1}), value.b)) &&...
               all(cellfun(@(c) isequal(c, value.c{1}), value.c)) &&...
               all(cellfun(@(d) isequal(d, value.d{1}), value.d)),...
               'DeltaConstantDelay:validateSample',...
               ['Specified value for delta "',this_delta.name,...
                '" must be time-invariant']);
        % Can only validate sample for discrete-time 
        if ~isempty(timestep)
            assert(logical(timestep),...
                   'DeltaConstantDelay:sample',...
                   ['Cannot express a continuous-time state-space sample of',...
                    ' DeltaConstantDelay'])
        end
        dim_out = this_delta.dim_out(1);
        input = ones(dim_out, 1) * (1:100);
        output = value.simulate(input);
        first_match = all(abs(output(1, :) - ones(dim_out, 1))...
                           < 1e-8 * ones(dim_out, 1), 1); 
        delay = find(first_match, 1, 'first') - 1;
        assert(~isempty(delay),...
               'DeltaConstantDelay:validateSample',...
               'Sample of uncertainty is not a delay operator')
        assert(delay <= this_delta.delay_max,...
               'DeltaConstantDelay:validateSample',...
               'Sample of uncertainty produces a delay greater than delay_max')
        expected_output = [zeros(dim_out, delay),...
                           input(:, 1 : end - delay)];
        assert(all(all(abs(output - expected_output) < 1e-8)),...
               'DeltaConstantDelay:validateSample',...
               'Sample of uncertainty is not a constant delay operator')
    end
end
end

%%  CHANGELOG
% Mar. 29, 2022: Added after v0.9.0 - Micah Fry (micah.fry@ll.mit.edu)