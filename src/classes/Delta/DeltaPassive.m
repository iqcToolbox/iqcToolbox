classdef DeltaPassive < Delta
%% DELTAPASSIVE class for uncertainties which are passive (may be dynamic and/or time-varying).
%  Extends the base class Delta.
%
%   extended methods:
%     DeltaPassive(name, dim_outin, horizon_period) :: Constructor
%     disp(this_delta) :: Display method
%     matchHorizonPeriod(this_delta, new_horizon_period) 
%                        :: Matches delta properties to new horizon_period
%     deltaToMultiplier(this_delta, varargin)
%                        :: Method for constructing a multiplier from a delta
%     normalizeDelta(this_delta) :: Changes unnormalized LFT to normalized one
%
%   See also Delta, DeltaPassive.DeltaPassive

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

methods
    function this_delta = DeltaPassive(name, dim_outin, horizon_period)
    %% DELTAPASSIVE constructor
    %
    %  d = DeltaPassive(name, dim_outin, horizon_period)
    %  d = DeltaPassive(name, dim_outin) assumes horizon_period == [0, 1]
    %  d = DeltaPassive(name) also assumes dim_outin == 1
    %
    %  Variables:
    %  ---------
    %    Input:
    %       name : string :: unique ID of the uncertainty (ex. 'dforce')
    %       dim_outin : natural :: input/output dimensions of uncertainty
    %       horizonperiod : 1 x 2 array of naturals :: [horizon, period] (in timesteps) 
    %                                                  of Delta properties
    %    Output:
    %       this_delta : DeltaPassive object :: the constructed DeltaPassive object
    %
    %  See also Delta, DeltaPassive, Delta.Delta

        % Defining defaults for missing arguments
        switch nargin
            case 1
                dim_outin = 1;
                horizon_period = [0, 1];
            case 2
                horizon_period = [0, 1];
        end
        % Calling Delta constructor
        this_delta@Delta(name, dim_outin, dim_outin, horizon_period);                    

        this_delta = matchHorizonPeriod(this_delta);
    end

    function disp(this_delta)
    %% DISP function for DeltaPassive object
    %
    %  disp(delta_passive_obj) (e.g., disp(DeltaPassive('d')) )
    %
    %  Variables:
    %  ---------
    %     Input:
    %        this_delta : DeltaPassive object
    %
    %  See also DeltaPassive, Delta.disp, Ulft.disp
    
        disp@Delta(this_delta, 'passive uncertainty')          
    end

    function this_del = matchHorizonPeriod(this_del, new_horizon_period)
    %% MATCHHORIZONPERIOD function to ensure properties of DeltaPassive 
    %  object match its own horizon_period, or a new_horizon_period
    %
    %  this_del = matchHorizonPeriod(this_del, new_horizon_period) will change the horizon_period and associated properties of this_del
    %  this_del = matchHorizonPeriod(this_del) only ensures that each pertinent property of this_del matches this_del.horizon_period
    %
    %  Variables:
    %  ---------
    %     Input:
    %       this_del : DeltaPassive object
    %       new_horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of new sequence
    %    Output:
    %       this_delta : DeltaPassive object
    %
    %  See also DeltaPassive.
    
    if nargin == 1
    % Ensuring that this_del.horizon_period matches with other properties 
    % of this_del
    
        % Assumes properties are a mix of sequences of length 1 or
        % horizon_period
        total_time = sum(this_del.horizon_period);
        if length(this_del.dim_out) ~= total_time
            assert(length(this_del.dim_out) == 1,...
                   'DeltaPassive:matchHorizonPeriod',...
                   'dim_out of %s is not compatible w/ horizon_period',...
                   this_del.name);
            this_del.dim_out = this_del.dim_out * ones(1, total_time);
        end
        if length(this_del.dim_in) ~= total_time
            assert(length(this_del.dim_in) == 1,...
                   'DeltaPassive:matchHorizonPeriod',...
                   'dim_in of %s is not compatible w/ horizon_period',...
                   this_del.name)
            this_del.dim_in = this_del.dim_in * ones(1, total_time);
        end
    else
    % Changing this_del.horizon_period and other properties of this_del to
    % a new horizon_period
        [indices, new_horizon_period] = ...
            makeNewIndices(this_del.horizon_period, new_horizon_period);
        
        % Set properties according to indices
        this_del.dim_out     = this_del.dim_out(indices);
        this_del.dim_in      = this_del.dim_in(indices);
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
    %        this_del : DeltaPassive object
    %     Output:
    %        multiplier : MultiplierPassive object
    %
    %  See also DeltaPassive.
        multiplier = MultiplierPassive(this_del);
    end
    
    function [del_diff, del_ave, del_scale, del_norm] = normalizeDelta(this_delta)
    %% NORMALIZEDELTA function for DeltaPassive object. This extends the default operation
    %  defined by the Delta superclass, such that DeltaPassive uncertainties
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
    
    [del_diff, del_ave, del_scale] = normalizeDelta@Delta(this_delta);
    del_norm = this_delta;
    end

    function value = sample(this_del, timestep)
        %% SAMPLE function for DeltaPassive.
        dim_in = num2cell(this_del.dim_in);
        if isempty(timestep)
        % Memoryless
            factor = cellfun(@(dim) randn(dim), dim_in, 'UniformOutput', false);
            d = cellfun(@(f) f' * f, factor, 'UniformOutput', false);
            value = toLft(d, this_del.horizon_period);
        else
        % Dynamic
            if timestep
                % Disrete-time
                system = drss(randi([1, 5]), 1, 1);
                system.a = system.a * 0.99;
            else
                % Continuous-time
                assert(isequal(this_del.horizon_period, [0, 1]),...
                       'DeltaPassive:sample',...
                       'Cannot create an LTV continuous-time system')
                system = drss(randi([1, 5]), 1, 1);
                system.a = system.a - 0.01 * eye(size(system.a));
            end
            % Force it positive
            real_part = nyquist(system + system');
            if min(real_part) < 0
                system = system - min(real_part / 2) * 1.01;
            end            
            system = system * eye(dim_in{1});
            value = toLft(system);
            value = matchHorizonPeriod(value, this_del.horizon_period);
        end
    end

    function validateSample(this_delta, value, ~)
        %% VALIDATESAMPLE function for DeltaPassive.
        % Validate base attributes
        validateSample@Delta(this_delta, value, [], true);
        % Square
        assert(all(this_delta.dim_in == this_delta.dim_out),...
               'DeltaPassive:validateSample',...
               ['Specified value for delta "', this_delta.name,...
                '"must have square dimensions.']);
        % Passivity    
        if isempty(value.timestep)
            dim_in = size(value, 2);
            total_time = sum(value.horizon_period);
            eigs = [];
            for i = 1:total_time
                d_eye = [value.d{i}; eye(dim_in(i))];
                mult = [zeros(dim_in(i)), eye(dim_in(i));
                        eye(dim_in(i)),   zeros(dim_in(i))];
                eigs = [eigs; eig(d_eye' * mult * d_eye)];
            end
            valid = min(eigs) >= 0;
        else
            value = value.addPerformance({PerformancePassive('passive')});
            options = AnalysisOptions('lmi_shift', 0, 'verbose', false);
            result = iqcAnalysis(value, 'analysis_options', options);
            valid = min(check(result.debug.constraints)) > -1e-10;
        end
        assert(valid,...
               'DeltaPassive:validateSample',...
               ['Specified value for delta "', this_delta.name,...
                '"must be passive.']);
    end
end
end

%%  CHANGELOG
% Oct. 10, 2021: Added after v0.6.0 - Micah Fry (micah.fry@ll.mit.edu)