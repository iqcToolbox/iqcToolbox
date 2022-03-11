classdef DeltaDelayZ < Delta
%% DELTADELAYZ class for the discrete-time delay operator Z, extends
% the base class Delta.
%
%   extended methods:
%     DeltaDelayZ(dim_state, timestep, horizon_period) :: Constructor
%     disp(this_delta) :: Display method
%     matchHorizonPeriod(this_delta, new_horizon_period) 
%                        :: Matches delta properties to new horizon_period
%     deltaToMultiplier(this_delta, varargin)
%                        :: Method for constructing a multiplier from a delta
%     normalizeDelta(this_delta) :: Changes unnormalized LFT to normalized one
%     toRct(this_delta) :: Converts Delta to object needed in lftToRct()
%
%   extended properties:
%     timestep : double :: timestep of time delay
%
%   See also Delta, DeltaIntegrator

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    timestep
end

methods
    function this_delta = DeltaDelayZ(dim_state, timestep, horizon_period)
    %% DELTADELAYZ constructor
    %
    %  d = DeltaDelayZ(dim_state, timestep, horizon_period)
    %  d = DeltaDelayZ(dim_state, timestep) assumes horizon_period == [0,1]
    %  d = DeltaDelayZ(dim_state) additionally leaves timestep undefined
    %  d = DeltaDelayZ() assumes dim_state == 1
    %
    %  Variables:
    %  ---------
    %    Input:
    %       dim_state : natural :: dimension of the state signal (i.e., x) may
    %                           be time-varying (i.e., an array)
    %       horizon_period : 1 x 2 array of int :: horizon and period for Z
    %       timestep : positive double :: timestep of delay in seconds
    %    Output:
    %       this_delta : DeltaDelayZ object
    %
    %  See also Delta, Delta.Delta, DeltaDelayZ

        % Defining defaults for missing arguments
        switch nargin
            case 0
                dim_state = 1;
                timestep = -1;
                horizon_period = [0, 1];
            case 1
                timestep = -1;
                horizon_period = [0, 1];
            case 2
                horizon_period = [0, 1];
        end
        validateattributes(dim_state,...
                           {'numeric'},...
                           {'integer', 'nonnegative', 'nonempty'})
        validateattributes(timestep,...
                           {'numeric'},...
                           {'nonempty'})
        
        % Checking inputs for specialized properties of DeltaDelayZ
        if length(dim_state) == 1
        % If state dimensions are time-invariant
            dim_out        = dim_state;
            dim_in         = dim_state;
        else
        % Map dim_state to dim_out & dim_in for time-varying dimensions
            dim_out = dim_state;
            dim_in  = circshift(dim_out, -1);
            dim_in(end) = dim_out(end - (horizon_period(2) - 1));
        end

        % Calling Delta constructor
        this_delta@Delta('DelayZ', dim_out, dim_in, horizon_period);
        % Setting specialized properties of DeltaDelayZ
        assert(all(timestep == -1) || all(timestep > 0),...
               'DeltaDelayZ:DeltaDelayZ',...
               'timestep must be -1 (unspecified) or positive')
        this_delta.timestep = timestep;
        
        this_delta = matchHorizonPeriod(this_delta);
    end

    function disp(this_delta)
    %% DISP function for DeltaDelayZ object
    %
    %  disp(delta_delay_z_obj) (e.g., disp(DeltaDelayZ()) )
    %
    %  Variables:
    %  ---------
    %     Input:
    %        this_delta : DeltaDelayZ object
    %
    %  See also DeltaDelayZ, Delta.disp, Ulft.disp

     
        disp@Delta(this_delta, 'delay operator')
        if this_delta.timestep == -1
            fprintf('%13s with an undefined timestep \n', '')
        else
            fprintf('%13s with a timestep of %3.1e seconds \n',...
                    '',...
                    this_delta.timestep(1))
        end
    end

    function this_del = matchHorizonPeriod(this_del, new_horizon_period)
    %% MATCHHORIZONPERIOD function to ensure properties of DeltaDelayZ 
    %  object match its own horizon_period, or a new_horizon_period
    %
    %  matchHorizonPeriod(this_del, new_horizon_period) will change the horizon_period and associated properties of this_del
    %  matchHorizonPeriod(this_del) only ensures that each pertinent property of this_del matches this_del.horizon_period
    %
    %  Variables:
    %  ---------
    %     Input:
    %       this_del : DeltaDelayZ object
    %       new_horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of new sequence
    %    Output:
    %       this_delta : DeltaDelayZ object
    %
    %  See also DeltaDelayZ.

    if nargin == 1
    % Ensuring that this_del.horizon_period matches with other properties 
    % of this_del
    
        % Assumes properties are a mix of sequences of length 1 or
        % horizon_period
        total_time = sum(this_del.horizon_period);
        if length(this_del.dim_out) ~= total_time
            assert(length(this_del.dim_out) == 1,...
                   'DeltaDelayZ:matchHorizonPeriod',...
                   'dim_out of %s is not compatible w/ horizon_period',...
                   this_del.name);
            this_del.dim_out = this_del.dim_out * ones(1, total_time);
        end
        if length(this_del.dim_in) ~= total_time
            assert(length(this_del.dim_in) == 1,...
                   'DeltaDelayZ:matchHorizonPeriod',...
                   'dim_in of %s is not compatible w/ horizon_period',...
                   this_del.name);
            this_del.dim_in = this_del.dim_in * ones(1, total_time);
        end
        if length(this_del.timestep) ~= total_time
            assert(length(this_del.timestep) == 1,...
               'DeltaDelayZ:matchHorizonPeriod',...
               'timestep of %s is not compatible w/ horizon_period',...
               this_del.name);
            this_del.timestep = this_del.timestep * ones(1, total_time);
        end
    else
    % Changing this_del.horizon_period and other properties of this_del to
    % a new horizon_period
    
        [indices, new_horizon_period] = ...
            makeNewIndices(this_del.horizon_period, new_horizon_period);

        % Set properties according to indices
        this_del.dim_out     = this_del.dim_out(indices);
        this_del.dim_in      = this_del.dim_in(indices);
        this_del.timestep    = this_del.timestep(indices);
        this_del.horizon_period = new_horizon_period;
        this_del = matchHorizonPeriod(this_del);
    end
    end
    
    function multiplier = deltaToMultiplier(this_del, varargin)                 %#ok<VANUS>
    %% DELTATOMULTIPLIER function to generate a multiplier from this object
    %  This function only exists to concretize the abstract class Delta.
    %  There is no multiplier for DeltaDelayZ, therefore this returns the
    %  default (or empty) Multiplier.
    %
    %  multiplier = deltaToMultiplier(this_del)
    %
    %  Variables:
    %  ---------
    %     Input:
    %       this_del : DeltaDelayZ object
    %     Output:
    %       multiplier : MultiplierDeltaDefault object
    %
    %  See also DeltaDelayZ.
        multiplier = MultiplierDeltaDefault(this_del);
    end
    
    function [del_diff, del_ave, del_scale, del_norm] = normalizeDelta(this_delta)
    %% NORMALIZEDELTA function for DeltaDelayZ object. This method extendeds
    %  the default method, but because DeltaDelayZ objects are not
    %  normalized, the extended behavior essentially mimics that of the default
    %  behavior without throwing a warning.
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
    
    [del_diff, del_ave, del_scale] = normalizeDelta@Delta(this_delta);
    del_norm = this_delta;
    end
    
    function output = toRct(this_delta)                                         %#ok<MANU>
    %% TORCT function for DeltaDelayZ object. This extends the base method by 
    %  not returning any RCT object, because the RCT does not have separate
    %  uncertainties to express the delay operation.
    %
    %    rct_obj = toRct(this_delta)
    % 
    %    Variables:
    %    ---------
    %      Input:
    %        this_delta : Delta object to be converted
    %      Output:
    %        rct_obj : empty array
    %
    %    See also Delta.toRct, lftToRct, rctToLft
        output = [];
    end

    function value = sample(~, ~)                                               %#ok<STOUT>
        %% SAMPLE function for DeltaDelayZ, which does not represent an uncertainty.
        error('DeltaDelayZ:sample', 'DeltaDelayZ cannot be sampled.');
    end

    function validateSample(~, ~, ~)
        %% VALIDATESAMPLE function for DeltaDelayZ, which does not represent an uncertainty.
        error('DeltaDelayZ:validateSample', 'DeltaDelayZ cannot be sampled.');
    end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0): Added sample and validateSample methods - Jason Nezvadovitz (jason.nezvadovitz@ll.mit.edu)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)