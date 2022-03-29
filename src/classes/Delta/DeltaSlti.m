classdef DeltaSlti < Delta
%% DELTASLTI class for static, linear, time-invariant uncertainties, extends
% the base class Delta.
%
%   extended methods:
%     DeltaSlti(name, dim_outin, lower_bound, upper_bound) :: Constructor
%     disp(this_delta) :: Display method
%     matchHorizonPeriod(this_delta, new_horizon_period) 
%                        :: Matches delta properties to new horizon_period
%     deltaToMultiplier(this_delta, varargin)
%                        :: Method for constructing a multiplier from a delta
%     normalizeDelta(this_delta) :: Changes unnormalized LFT to normalized one
%     toRct(this_delta) :: Converts Delta to object needed in lftToRct()
%
%   extended properties:
%     upper_bound : double :: upper bound of uncertainty
%     lower_bound : double :: lower bound of uncertainty
%
%   See also Delta.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    upper_bound
    lower_bound
end

methods
    function this_delta = DeltaSlti(name, dim_outin, lower_bound, upper_bound)
    %% DELTASLTI constructor
    %
    %  d = DeltaSlti(name, dim_outin, lower_bound, upper_bound)
    %  d = DeltaSlti(name, dim_outin) assumes lower_bound == -1,
    %                                         upper_bound == 1
    %  d = DeltaSlti(name) additionally assumes dim_outin == 1
    %
    %  Variables:
    %  ---------
    %    Input:
    %       name : string :: unique ID of the uncertainty (ex. 'dmass')
    %       dim_outin : natural :: output/input dimensions of uncertainty
    %       lower_bound : double :: lower_bound of uncertainty
    %       upper_bound : double :: upper_bound of uncertainty
    %    Output:
    %       this_delta : DeltaSlti object
    %
    %  See also Delta, Delta.Delta, DeltaSlti

        % Defining defaults for missing arguments
        switch nargin
            case 0
                error('DeltaSlti:DeltaSlti',...
                      ['Must provide at least one input to specify the',...
                       ' name of the uncertainty'])
            case 1
                dim_outin = 1;
                lower_bound = -1.0;
                upper_bound = 1.0;                
            case 2
                lower_bound = -1.0;
                upper_bound = 1.0;
        end
        horizon_period = [0, 1];
        assert(all(dim_outin(1) == dim_outin),...
               'DeltaSlti:DeltaSlti',...
               'dimensions must be constant for DeltaSlti objects')

           % Calling Delta constructor
        this_delta@Delta(name, dim_outin, dim_outin, horizon_period);                    

        % Checking inputs for specialized properties of DeltaSlti
        validateattributes(lower_bound,...
                           {'numeric'},...
                           {'nonnan', 'finite', 'nonempty'},...
                           mfilename);
        validateattributes(upper_bound,...
                           {'numeric'},...
                           {'nonnan', 'finite', 'nonempty'},...
                           mfilename);
        assert(lower_bound <= upper_bound,...
               'DeltaSlti:DeltaSlti',...
               'lower_bound is greater than upper_bound');
        assert(all(lower_bound(1) == lower_bound),...
               'DeltaSlti:DeltaSlti',...
               'lower_bound must be constant for DeltaSlti objects')
        assert(all(upper_bound(1) == upper_bound),...
               'DeltaSlti:DeltaSlti',...
               'upper_bound must be constant for DeltaSlti objects')

        this_delta.lower_bound = lower_bound;
        this_delta.upper_bound = upper_bound;
        this_delta = matchHorizonPeriod(this_delta);

    end

    function disp(this_delta)
    %% DISP function for DeltaSlti object
    %
    %  disp(delta_slti_obj) (e.g., disp(DeltaSlti('d')) )
    %
    %  Variables:
    %  ---------
    %     Input:
    %       this_delta : DeltaSlti object 
    %
    %  See also DeltaSlti, Delta.disp, SequenceDelta.disp

        disp@Delta(this_delta, 'SLTI uncertainty')  
        fprintf('%13s within the bounds: [%3.1f, %3.1f] \n',...
                '   ',...
                this_delta.lower_bound(1),...
                this_delta.upper_bound(1))            
    end

    function this_del = matchHorizonPeriod(this_del, new_horizon_period)
    %% MATCHHORIZONPERIOD function to ensure properties of DeltaSlti 
    %  object match a total_time length
    %
    %  matchHorizonPeriod(this_del, new_horizon_period) will change the horizon_period and associated properties of this_del
    %  matchHorizonPeriod(this_del) only ensures that each pertinent property of this_del matches this_del.horizon_period
    %
    %  Variables:
    %  ---------
    %     Input:
    %        this_del : DeltaSlti object
    %         new_horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of new sequence
    %    Output:
    %       this_delta : DeltaSlti object
    %
    %  See also DeltaSlti
    
    if nargin == 1
    % Ensuring that this_del.horizon_period matches with other properties 
    % of this_del
    
        % Assumes properties are a mix of sequences of length 1 or
        % horizon_period
        total_time = sum(this_del.horizon_period);
        if length(this_del.dim_out) ~= total_time
            assert(length(this_del.dim_out) == 1,...
                   'DeltaSlti:matchHorizonPeriod',...
                   'dim_out of %s is not compatible w/ horizon_period',...
                   this_del.name);
            this_del.dim_out = this_del.dim_out * ones(1, total_time);
        end
        if length(this_del.dim_in) ~= total_time
            assert(length(this_del.dim_in) == 1,...
                   'DeltaSlti:matchHorizonPeriod',...
                   'dim_in of %s is not compatible w/ horizon_period',...
                   this_del.name)
            this_del.dim_in = this_del.dim_in * ones(1, total_time);
        end
        if length(this_del.lower_bound) ~= total_time
            assert(length(this_del.lower_bound) == 1,...
                   'DeltaSlti:matchHorizonPeriod',...
                   'lower_bound of %s is not compatible w/ horizon_period',...
                   this_del.name)
            this_del.lower_bound = ...
                this_del.lower_bound * ones(1, total_time);
        end
        if length(this_del.upper_bound) ~= total_time
            assert(length(this_del.upper_bound) == 1,...
                   'DeltaSlti:matchHorizonPeriod',...
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
        this_del.lower_bound = this_del.lower_bound(1) * ones(1, total_time);
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
    %       this_del : DeltaSlti object
    %     Output:
    %       multiplier : MultiplierSlti object
    %
    %  See also DeltaSlti.
        multiplier = MultiplierSlti(this_del, varargin{:});
    end
    
    function [del_diff, del_ave, del_scale, del_norm] = normalizeDelta(this_delta)
    %% NORMALIZEDELTA function for DeltaSlti object. This extends the default operation
    %  defined by the Delta superclass, such that DeltaSlti uncertainties
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
        diff = (this_delta.upper_bound(i) - this_delta.lower_bound(i)) / 2;
        ave  = (this_delta.upper_bound(i) + this_delta.lower_bound(i)) / 2;
        del_diff{i}  = diff * eye(this_delta.dim_out(i));
        del_ave{i}   = ave * eye(this_delta.dim_out(i));
        del_scale{i} = eye(this_delta.dim_in(i));
    end
    
    % Define normalized Delta
    del_norm = this_delta;
    del_norm.lower_bound = -1;
    del_norm.upper_bound = 1;
    del_norm = matchHorizonPeriod(del_norm, this_delta.horizon_period);
    end
    
    function rct_obj = toRct(this_delta)
    %% TORCT function for DeltaSlti object. This extends the base method by 
    %  returning ureal objects.
    %
    %    rct_obj = toRct(this_delta)
    % 
    %    Variables:
    %    ---------
    %      Input:
    %        this_delta : DeltaSlti object to be converted
    %      Output:
    %        rct_obj : ureal or umat object (from the Robust Control Toolbox)
    %    See also Delta.toRct, lftToRct, rctToLft
        rct_obj = cell(this_delta.dim_in, 1);
        % DeltaSlti assumes the nominal value is zero
        nominal_value = 0;
        for i = 1:this_delta.dim_in
            rct_obj{i} = ureal(this_delta.name,...
                             nominal_value,...
                             'Range', [this_delta.lower_bound,...
                                       this_delta.upper_bound]);
        end
        rct_obj = blkdiag(rct_obj{:});    
    end

    function value = sample(this_delta, ~)
        %% SAMPLE function for DeltaSlti.
        mag = (this_delta.upper_bound(1) - this_delta.lower_bound(1))*rand()...
              + this_delta.lower_bound(1);
        value = toLft(mag*eye(this_delta.dim_out(1)));
        value = value.matchHorizonPeriod(this_delta.horizon_period);
    end

    function validateSample(this_delta, value, ~)
        %% VALIDATESAMPLE function for DeltaSlti.
        % Validate base attributes
        validateSample@Delta(this_delta, value, [], true);
        % Time invariant
        assert( all(cellfun(@(d) isequal(d, value.d{1}), value.d)),...
               'DeltaSlti:validateSample',...
               ['Specified value for delta "',this_delta.name,...
                '" must be time-invariant.']);
        % Memoryless
        assert(isempty(value.delta.deltas),...
               'DeltaSlti:validateSample',...
               ['Specified value for delta "',this_delta.name,...
               '" must be memoryless.']);
        % Isotropic
        assert(isdiag(value.d{1}) && all(diag(value.d{1}) == value.d{1}(1,1)),...
               'DeltaSlti:validateSample',...
               ['Specified value for delta "',this_delta.name,...
               '" must be a scaling of the identity.']);
        % Bounded
        assert((this_delta.lower_bound(1) <= value.d{1}(1,1)) && ...
               (value.d{1}(1,1) <= this_delta.upper_bound(1)),...
               'DeltaSlti:validateSample',...
               ['Specified value for delta "',this_delta.name,...
               '" must be within the bounds.']);
    end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0): Added sample and validateSample methods - Jason Nezvadovitz (jason.nezvadovitz@ll.mit.edu)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)