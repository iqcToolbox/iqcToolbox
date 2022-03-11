classdef DeltaDlti < Delta
%% DELTADLTI class for dynamic, linear, time-invariant uncertainties, 
% extends the base class Delta.
%
%   extended methods:
%     DeltaDlti(name, dim_out, dim_in, upper_bound) :: Constructor
%     disp(this_delta) :: Display method
%     matchHorizonPeriod(this_delta, new_horizon_period) 
%                        :: Matches delta properties to new horizon_period
%     deltaToMultiplier(this_delta, varargin)
%                        :: Method for constructing a multiplier from a delta
%     normalizeDelta(this_delta) :: Changes unnormalized LFT to normalized one
%     toRct(this_delta) :: Converts Delta to object needed in lftToRct()
%
%   extended properties:
%     upper_bound : double :: upper bound of hinf norm of uncertainty
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
    function this_delta = DeltaDlti(name, dim_out, dim_in, upper_bound)
    %% DELTADLTI constructor
    %
    %  d = DeltaDlti(name, dim_out, dim_in, upper_bound)
    %  d = DeltaDlti(name, dim_out, dim_in) assumes upper_bound == 1
    %  d = DeltaDlti(name, dim_outin) also assumes dim_in == dim_out
    %  d = DeltaDlti(name) also assumes dim_out == dim_in == 1
    %
    %  Variables:
    %  ---------
    %    Input:
    %       name : string :: unique ID of the uncertainty (ex. 'dmass')
    %       dim_out : natural :: output dimensions of uncertainty
    %       dim_in : natural :: input dimensions of uncertainty
    %       upper_bound : double :: upper_bound of uncertainty
    %    Output:
    %       this_delta : DeltaDlti object
    %
    %  See also Delta, Delta.Delta, DeltaDlti

        % Defining defaults for missing arguments
        switch nargin
            case 0
                error('DeltaDlti:DeltaDlti',...
                      ['Must provide at least one input to specify the',...
                       ' name of the uncertainty'])
            case 1
                dim_out     = 1;
                dim_in      = 1;
                upper_bound = 1.0;                
            case 2
                dim_in      = dim_out;
                upper_bound = 1.0;
            case 3
                upper_bound = 1.0;                
        end
        horizon_period = [0, 1];
        validateattributes(dim_out, {'numeric'}, {'positive'});
        validateattributes(dim_in, {'numeric'}, {'positive'});
        assert(all(dim_out(1) == dim_out) && all(dim_in(1) == dim_in),...
               'DeltaDlti:DeltaDlti',...
               'in/out dimensions must be constant for DeltaDlti objects')
       

       % Calling Delta constructor
        this_delta@Delta(name, dim_out, dim_in, horizon_period);                    

        % Checking inputs for specialized properties of DeltaDlti
        validateattributes(upper_bound,...
                           {'numeric'},...
                           {'nonempty', 'nonnan', 'finite', 'positive'},...
                           mfilename);
        assert(all(upper_bound(1) == upper_bound),...
               'DeltaDlti:DeltaDlti',...
               'upper_bound must be constant for DeltaDlti objects')

        this_delta.upper_bound = upper_bound;
        this_delta = matchHorizonPeriod(this_delta);

    end

    function disp(this_delta)
    %% DISP function for DeltaDlti object
    %
    %  disp(delta_dlti_obj) (e.g., disp(DeltaDlti('d')) )
    %
    %  Variables:
    %  ---------
    %     Input:
    %        this_delta : DeltaDlti object
    %
    %  See also DeltaDlti, Delta.disp, SequenceDelta.disp
        disp@Delta(this_delta, 'DLTI uncertainty')  
        fprintf('%13s whose hinf norm is less than: %3.1f \n',...
                '   ',...
                this_delta.upper_bound(1))            
    end

    function this_del = matchHorizonPeriod(this_del, new_horizon_period)
    %% MATCHHORIZONPERIOD function to ensure properties of DeltaDlti 
    %  object match a total_time length
    %
    %  matchHorizonPeriod(this_del, new_horizon_period) will change the horizon_period and associated properties of this_del
    %  matchHorizonPeriod(this_del) only ensures that each pertinent property of this_del matches this_del.horizon_period
    %
    %  Variables:
    %  ---------
    %     Input:
    %       this_del : DeltaDlti object
    %       new_horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of new sequence
    %    Output:
    %       this_delta : DeltaDlti object
    %
    %  See also DeltaDlti

    if nargin == 1
    % Ensuring that this_del.horizon_period matches with other properties 
    % of this_del
    
        % Assumes properties are a mix of sequences of length 1 or
        % horizon_period
        total_time = sum(this_del.horizon_period);
        if length(this_del.dim_out) ~= total_time
            assert(length(this_del.dim_out) == 1,...
                   'DeltaDlti:matchHorizonPeriod',...
                   'dim_out of %s is not compatible w/ horizon_period',...
                   this_del.name);
            this_del.dim_out = this_del.dim_out * ones(1, total_time);
        end
        if length(this_del.dim_in) ~= total_time
            assert(length(this_del.dim_in) == 1,...
                   'DeltaDlti:matchHorizonPeriod',...
                   'dim_in of %s is not compatible w/ horizon_period',...
                   this_del.name)
            this_del.dim_in = this_del.dim_in * ones(1, total_time);
        end
        if length(this_del.upper_bound) ~= total_time
            assert(length(this_del.upper_bound) == 1,...
                   'DeltaDlti:matchHorizonPeriod',...
                   'upper_bound of %s is not compatible w/ horizon_period',...
                   this_del.name)
            this_del.upper_bound = ...
                this_del.upper_bound * ones(1, total_time);
        end
    else
    % Changing this_del.horizon_period and other properties of this_del to
    % a new horizon_period
        
        % Set properties according to indices
        total_time = sum(new_horizon_period);
        this_del.dim_out     = this_del.dim_out(1) * ones(1, total_time);
        this_del.dim_in      = this_del.dim_in(1) * ones(1, total_time);
        this_del.upper_bound = this_del.upper_bound(1) * ones(1, total_time);
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
    %       this_del : DeltaDlti object
    %     Output:
    %       multiplier : MultiplierDlti object
    %
    %  See also DeltaDlti.
        multiplier = MultiplierDlti(this_del, varargin{:});
    end
    
    
    function rct_obj = toRct(this_delta)
    %% TORCT function for DeltaDlti object. This extends the base method by 
    %  returning ultidyn objects having a specified Bound.
    %
    %    rct_obj = toRct(this_delta)
    % 
    %    Variables:
    %    ---------
    %      Input:
    %        this_delta : Delta object to be converted
    %      Output:
    %        rct_obj : ultidyn object (from the Robust Control Toolbox)
    %
    %    See also Delta.toRct, lftToRct, rctToLft
        rct_obj = ultidyn(this_delta.name,...
                         [this_delta.dim_out, this_delta.dim_in],...
                         'Bound', this_delta.upper_bound);
    end
    
    function [del_diff, del_ave, del_scale, del_norm] = normalizeDelta(this_delta)
    %% NORMALIZEDELTA function for DeltaDlti object. This extends the default operation
    %  defined by the Delta superclass, such that DeltaDlti uncertainties
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
    
    function value = sample(this_delta, timestep)
        %% SAMPLE function for DeltaDlti.
        % Check if dynamic
        assert(~isempty(timestep),...
               'DeltaDlti:sample',...
               'LFT must be dynamic to sample a DeltaDlti');
        % Generate unbounded
        if timestep
            memory = DeltaDelayZ(randi([1,10]), timestep(1));
        else
            assert(isequal(this_delta.horizon_period, [0, 1]),...
                   'DeltaDlti:sample',...
                   ['Must have horizon_period = [0, 1] in order to ',...
                    'sample a continuous-time DeltaDlti'])
            memory = DeltaIntegrator(randi([1,10]));
        end
        value = Ulft.random('dim_in', this_delta.dim_in(1),...
                            'dim_out', this_delta.dim_out(1),...
                            'req_deltas', {memory},...
                            'num_deltas', 1,...
                            'horizon_period', [0,1]);
        % Scale norm
        scaling = this_delta.upper_bound(1)*rand() / norm(lftToSs(value), 'inf');
        value = value * (scaling*eye(this_delta.dim_in(1)));
        value = value.matchHorizonPeriod(this_delta.horizon_period);
    end

    function validateSample(this_delta, value, timestep)
        %% VALIDATESAMPLE function for DeltaDlti.
        % Validate base attributes
        validateSample@Delta(this_delta, value, timestep, true);
        % Time invariant
        assert(all(cellfun(@(a) isequal(a, value.a{1}), value.a)) &&...
               all(cellfun(@(b) isequal(b, value.b{1}), value.b)) &&...
               all(cellfun(@(c) isequal(c, value.c{1}), value.c)) &&...
               all(cellfun(@(d) isequal(d, value.d{1}), value.d)),...
               'DeltaDlti:validateSample',...
               ['Specified value for delta "',this_delta.name,...
                '" must be time-invariant']);
        % Memory type
        if ~isempty(timestep)
            if timestep(1)
                assert(isa(value.delta.deltas{1}, 'DeltaDelayZ'),...
                       'DeltaDlti:validateSample',...
                       ['Specified value for delta "',this_delta.name,...
                        '" must be discrete-time dynamic.']);
                % Reduce to horizon_period [0,1] for bounded check
                value.delta.deltas{1} = ...
                    DeltaDelayZ(value.delta.deltas{1}.dim_out(1),...
                                value.delta.deltas{1}.timestep(1));
            else
                assert(isa(value.delta.deltas{1}, 'DeltaIntegrator'),...
                       'DeltaDlti:validateSample',...
                       ['Specified value for delta "',this_delta.name,...
                        '" must be continuous-time dynamic.']);
            end
        end
        % Bounded (check requires horizon_period [0,1])
        value = Ulft(value.a{1}, value.b{1}, value.c{1}, value.d{1}, value.delta);
        assert(norm(lftToSs(value), 'inf') <= this_delta.upper_bound(1),...
               'DeltaDlti:validateSample',...
               ['Specified value for delta "',this_delta.name,...
                '" must have Hinf-norm within the bound.']);
    end
end

end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0): Added sample and validateSample methods - Jason Nezvadovitz (jason.nezvadovitz@ll.mit.edu)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)