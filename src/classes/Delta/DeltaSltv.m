classdef DeltaSltv < Delta
%% DELTASLTV class for static, linear, time-varying uncertainties, extends
% the base class Delta.
%
%   extended methods:
%     DeltaSltv(name, dim_outin, lower_bound, upper_bound, horizon_period) :: Constructor
%     disp(this_delta) :: Display method
%     matchHorizonPeriod(this_delta, new_horizon_period) 
%                        :: Matches delta properties to new horizon_period
%     deltaToMultiplier(this_delta, varargin)
%                        :: Method for constructing a multiplier from a delta
%     normalizeDelta(this_delta) :: Changes unnormalized LFT to normalized one
%
%   extended properties:
%     upper_bound : double :: upper bound of uncertainty
%     lower_bound : double :: lower bound of uncertainty
%
%   See also Delta, DeltaSltvRepeated, DeltaSltvRateBnd

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    upper_bound
    lower_bound
end

methods
function this_delta = DeltaSltv(name,...
                                dim_outin,...
                                lower_bound,...
                                upper_bound,...
                                horizon_period)
%% DELTASLTV constructor
%
%  d = DeltaSltv(name, dim_outin, lower_bound, upper_bound, horizon_period)
%  d = DeltaSltv(name, dim_outin, lower_bound, upper_bound) assumes horizon_period == [0, 1]
%  d = DeltaSltv(name, dim_outin) also assumes lower_bound == -1,
%                                              upper_bound == 1
%  d = DeltaSltv(name) also assumes dim_outin == 1
%
%  Variables:
%  ---------
%    Input:
%       name : string :: unique ID of the uncertainty (ex. 'dmass')
%       dim_outin : natural :: output/input dimensions of uncertainty
%       lower_bound : double :: lower_bound of uncertainty
%       upper_bound : double :: upper_bound of uncertainty
%       horizon_period : 1 x 2 array of integers :: horizon and period of uncertainty
%    Output:
%       this_delta : DeltaSltv object
%
%  See also Delta, Delta.Delta, DeltaSltv

    % Defining defaults for missing arguments
    switch nargin
        case 0
            error('DeltaSltv:DeltaSltv',...
                  ['Must provide at least one input to specify the',...
                   ' name of the uncertainty'])
        case 1
            dim_outin = 1;
            lower_bound = -1.0;
            upper_bound = 1.0;
            horizon_period = [0, 1];
        case 2
            lower_bound = -1.0;
            upper_bound = 1.0;
            horizon_period = [0, 1];
        case 4
            horizon_period = [0, 1];
    end

    % Calling Delta constructor
    this_delta@Delta(name, dim_outin, dim_outin, horizon_period);                    

    % Checking inputs for specialized properties of DeltaSltv
    validateattributes(lower_bound,...
                       'numeric',...
                       {'nonnan', 'finite', 'nonempty'},...
                       mfilename);
    validateattributes(upper_bound,...
                       'numeric',...
                       {'nonnan', 'finite', 'nonempty'},...
                       mfilename);
    assert(all(lower_bound <= upper_bound, 'all'),...
           'DeltaSltv:DeltaSltv',...
           'lower_bound is greater than upper_bound');

    this_delta.lower_bound = lower_bound;
    this_delta.upper_bound = upper_bound;
    this_delta = matchHorizonPeriod(this_delta);

end

function disp(this_delta)
%% DISP function for DeltaSltv object
%
%  disp(delta_sltv_obj) (e.g., disp(DeltaSltv('d')) )
%
%  Variables:
%  ---------
%     Input:
%       this_delta : DeltaSltv object 
%
%  See also DeltaSltv, Delta.disp, SequenceDelta.disp    
    disp@Delta(this_delta, 'SLTV uncertainty')  
    fprintf('%13s within the bounds: [%3.1f, %3.1f] \n',...
            '   ',...
            this_delta.lower_bound(1),...
            this_delta.upper_bound(1))            
end

function this_del = matchHorizonPeriod(this_del, new_horizon_period)
%% MATCHHORIZONPERIOD function to ensure properties of DeltaSltv 
%  object match a total_time length
%
%  matchHorizonPeriod(this_del, new_horizon_period) will change the horizon_period and associated properties of this_del
%  matchHorizonPeriod(this_del) only ensures that each pertinent property of this_del matches this_del.horizon_period
%
%  Variables:
%  ---------
%     Input:
%       this_del : DeltaSltv object
%       new_horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of new sequence
%    Output:
%       this_delta : DeltaSltv object
%
%  See also DeltaSltv

if nargin == 1
% Ensuring that this_del.horizon_period matches with other properties 
% of this_del

    % Assumes properties are a mix of sequences of length 1 or
    % horizon_period
    total_time = sum(this_del.horizon_period);
    if length(this_del.dim_out) ~= total_time
        assert(length(this_del.dim_out) == 1,...
               'DeltaSltv:matchHorizonPeriod',...
               'dim_out of %s is not compatible w/ horizon_period',...
               this_del.name);
        this_del.dim_out = this_del.dim_out * ones(1, total_time);
    end
    if length(this_del.dim_in) ~= total_time
        assert(length(this_del.dim_in) == 1,...
               'DeltaSltv:matchHorizonPeriod',...
               'dim_in of %s is not compatible w/ horizon_period',...
               this_del.name)
        this_del.dim_in = this_del.dim_in * ones(1, total_time);
    end
    if length(this_del.lower_bound) ~= total_time
        assert(length(this_del.lower_bound) == 1,...
               'DeltaSltv:matchHorizonPeriod',...
               'lower_bound of %s is not compatible w/ horizon_period',...
               this_del.name)
        this_del.lower_bound = ...
            this_del.lower_bound * ones(1, total_time);
    end
    if length(this_del.upper_bound) ~= total_time
        assert(length(this_del.upper_bound) == 1,...
               'DeltaSltv:matchHorizonPeriod',...
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
    this_del.upper_bound = this_del.upper_bound(indices);
    this_del.lower_bound = this_del.lower_bound(indices);
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
%        multiplier : MultiplierSltv object
%
%  See also DeltaSltv.
    multiplier = MultiplierSltv(this_del, varargin{:});
end

function [del_diff, del_ave, del_scale, del_norm] = normalizeDelta(this_delta)
%% NORMALIZEDELTA function for DeltaSltv object. This extends the default operation
%  defined by the Delta superclass, such that DeltaSltv uncertainties
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
horizon_period = this_delta.horizon_period;
name           = this_delta.name;
dim_outin      = this_delta.dim_in;
lower_bound    = -ones(1, sum(horizon_period));
upper_bound    = ones(1, sum(horizon_period));

del_norm = DeltaSltv(name, dim_outin, lower_bound, upper_bound, horizon_period);
end

function value = sample(this_delta, ~)
    %% SAMPLE function for DeltaSltv.
    total_time = sum(this_delta.horizon_period);
    a = cell(1, total_time);
    b = cell(1, total_time);
    c = cell(1, total_time);
    d = cell(1, total_time);
    for t = 1 : total_time
        a{t} = zeros(0,0);
        b{t} = zeros(0,this_delta.dim_in(t));
        c{t} = zeros(this_delta.dim_out(t),0);
        d_mag = (this_delta.upper_bound(t) - this_delta.lower_bound(t))*rand()...
                + this_delta.lower_bound(t);
        d{t} = d_mag * eye(this_delta.dim_out(t));
    end
    value = Ulft(a, b, c, d, {}, 'horizon_period', this_delta.horizon_period);
end

function validateSample(this_delta, value, ~)
    %% VALIDATESAMPLE function for DeltaSltv.
    % Validate base attributes
    validateSample@Delta(this_delta, value, [], true);
    % Memoryless
    assert(isempty(value.delta.deltas),...
           'DeltaSltv:validateSample',...
           ['Specified value for delta "',this_delta.name,...
            '" must be memoryless.']);
    % Isotropic
    assert(all(cellfun(@(d) isdiag(d) && all(diag(d) == d(1,1)), value.d)),...
           'DeltaSltv:validateSample',...
           ['Specified value for delta "',this_delta.name,...
            '" must be a scaling of the identity at all times.']);
    % Bounded
    lb = num2cell(this_delta.lower_bound);
    ub = num2cell(this_delta.upper_bound);
    assert(all(cellfun(@(d, lb, ub) (lb <= d(1,1)) && (d(1,1) <= ub),...
                       value.d, lb, ub)),...
           'DeltaSltv:validateSample',...
           ['Specified value for delta "',this_delta.name,...
            '" must be within the bounds at all times.']);
end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0): Added sample and validateSample methods - Jason Nezvadovitz (jason.nezvadovitz@ll.mit.edu)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)