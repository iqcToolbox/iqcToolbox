classdef DeltaSltvRateBnd < Delta
%% DELTASLTVRATEBND class for static, linear, time-varying rate-bounded 
% uncertainties, extends the base class Delta.
%
%   extended methods:
%     DeltaSltvRateBnd(name, dim_outin, lower_bound, upper_bound, lower_rate, upper_rate, horizon_period, basis_length) :: Constructor
%     disp(this_delta) :: Display method
%     matchHorizonPeriod(this_delta, new_horizon_period) 
%                        :: Matches delta properties to new horizon_period
%     deltaToMultiplier(this_delta, varargin)
%                        :: Method for constructing a multiplier from a delta
%     normalizeDelta(this_delta) :: Changes unnormalized LFT to normalized one
%     recastMatricesAndDelta(this_delta) :: Method for changing initial LFT to the analyzed one
%
%   extended properties:
%     upper_bound : double :: upper bound of uncertainty
%     lower_bound : double :: lower bound of uncertainty
%     upper_rate  : double :: upper bound of uncertainty rate of change
%     lower_rate  : double :: lower bound of uncertainty rate of change
%     basis_length : double :: length of basis_function for multiplier
%
%   See also Delta, DeltaSltv, DeltaSltvRateBndImpl

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    upper_bound
    lower_bound
    upper_rate
    lower_rate
    basis_length
end

methods
function this_delta = DeltaSltvRateBnd(name,...
                                       dim_outin,...
                                       lower_bound,...
                                       upper_bound,...
                                       lower_rate,...
                                       upper_rate,...
                                       horizon_period,...
                                       basis_length)
%% DELTASLTVRATEBND constructor
%
%  d = DeltaSltvRateBnd(name, dim_outin, lower_bound, upper_bound, lower_rate, upper_rate, horizon_period, basis_length)
%  d = DeltaSltvRateBnd(name, dim_outin, lower_bound, upper_bound, lower_rate, upper_rate, horizon_period) assumes basis_length == 2
%  d = DeltaSltvRateBnd(name, dim_outin, lower_bound, upper_bound, lower_rate, upper_rate) also assumes horizon_period == [0, 1]
%  d = DeltaSltvRateBnd(name, dim_outin, lower_bound, upper_bound) also assumes lower_rate == -abs(ub - lb), upper_rate == abs(ub - lb)
%  d = DeltaSltvRateBnd(name, dim_outin) also assumes -lower_bound == upper_bound == 1
%  d = DeltaSltvRateBnd(name) also assumes dim_outin == 1
%
%  Variables:
%  ---------
%    Input:
%       name : string :: unique ID of the uncertainty (ex. 'dmass')
%       dim_outin : natural :: output/input dimensions of uncertainty
%       lower_bound : double :: lower_bound of uncertainty
%       upper_bound : double :: upper_bound of uncertainty
%       lower_rate : double :: lower_bound of uncertainty rate of change
%       upper_rate : double :: upper_bound of uncertainty rate of change
%       horizon_period : 1 x 2 array of integers :: horizon and period of uncertainty
%       basis_length : natural :: length of multiplier basis
%    Output:
%       this_delta : DeltaSltvRateBnd object
%
%  See also Delta, Delta.Delta, DeltaSltvRateBnd

    % Defining defaults for missing arguments
    switch nargin
        case 0
            error('DeltaSltvRateBnd:DeltaSltvRateBnd',...
                  ['Must provide at least the',...
                   ' name of the uncertainty'])
        case 1
            dim_outin = 1;
            lower_bound = -1.0;
            upper_bound = 1.0;
            lower_rate  = -2.0;
            upper_rate  = 2.0;
            horizon_period = [0, 1];
            basis_length = 2;
        case 2
            lower_bound = -1.0;
            upper_bound =  1.0;
            lower_rate  = -2.0;
            upper_rate  =  2.0;
            horizon_period = [0, 1];
            basis_length = 2;
        case 4
            lower_rate = -abs(upper_bound - lower_bound);
            upper_rate = abs(upper_bound - lower_bound);
            horizon_period = [0, 1];
            basis_length = 2;
        case 6
            horizon_period = [0, 1];
            basis_length = 2;
        case 7
            basis_length = 2;
    end

    % Calling Delta constructor
    this_delta@Delta(name, dim_outin, dim_outin, horizon_period);                    

    % Checking inputs for specialized properties of DeltaSltvRateBnd
    validateattributes(lower_bound,...
                       'numeric',...
                       {'nonnan', 'finite', 'nonempty'},...
                       mfilename);
    validateattributes(upper_bound,...
                       'numeric',...
                       {'nonnan', 'finite', 'nonempty'},...
                       mfilename);
    assert(all(lower_bound <= upper_bound, 'all'),...
           'DeltaSltvRateBnd:DeltaSltvRateBnd',...
           'lower_bound is greater than upper_bound');
    validateattributes(lower_rate,...
                       'numeric',...
                       {'nonnan', 'finite', 'nonempty'},...
                       mfilename);
    validateattributes(upper_rate,...
                       'numeric',...
                       {'nonnan', 'finite', 'nonempty'},...
                       mfilename);
    assert(all(lower_rate <= upper_rate, 'all'),...
           'DeltaSltvRateBnd:DeltaSltvRateBnd',...
           'lower_rate is greater than upper_rate');
    validateattributes(basis_length,...
                       'numeric',...
                       {'integer', 'nonempty', 'scalar', '>', 1})
                       
    this_delta.lower_bound  = lower_bound;
    this_delta.upper_bound  = upper_bound;
    this_delta.lower_rate   = lower_rate;
    this_delta.upper_rate   = upper_rate;
    this_delta.basis_length = basis_length;
    this_delta = matchHorizonPeriod(this_delta);

end

function disp(this_delta)
%% DISP function for DeltaSltvRateBnd object
%
%  disp(delta_sltv_obj) (e.g., disp(DeltaSltvRateBnd('d')) )
%
%  Variables:
%  ---------
%     Input:
%       this_delta : DeltaSltvRateBnd object  
%
%  See also DeltaSltvRateBnd, Delta.disp, SequenceDelta.disp    
    disp@Delta(this_delta, 'SLTV, rate-bounded uncertainty')  
    fprintf(['%13s within the bounds: [%3.1f, %3.1f] \n',...
             '%13s whose rate-of-change is within: [%3.1f, %3.1f] \n'],...
            '',...
            this_delta.lower_bound(1),...
            this_delta.upper_bound(1),...
            '',...
            this_delta.lower_rate(1),...
            this_delta.upper_rate(1))            
end

function this_del = matchHorizonPeriod(this_del, new_horizon_period)
%% MATCHHORIZONPERIOD function to ensure properties of DeltaSltvRateBnd 
%  object match a total_time length
%
%  matchHorizonPeriod(this_del, new_horizon_period) will change the horizon_period and associated properties of this_del
%  matchHorizonPeriod(this_del) only ensures that each pertinent property of this_del matches this_del.horizon_period
%
%  Variables:
%  ---------
%     Input:
%       this_del : DeltaSltvRateBnd object
%       new_horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of new sequence
%     Output:
%        this_delta : DeltaSltvRateBnd object
%
%  See also DeltaSltvRateBnd

if nargin == 1
% Ensuring that this_del.horizon_period matches with other properties 
% of this_del

    % Assumes properties are a mix of sequences of length 1 or
    % horizon_period
    total_time = sum(this_del.horizon_period);
    if length(this_del.dim_out) ~= total_time
        assert(length(this_del.dim_out) == 1,...
               'DeltaSltvRateBnd:matchHorizonPeriod',...
               'dim_out of %s is not compatible w/ horizon_period',...
               this_del.name);
        this_del.dim_out = this_del.dim_out * ones(1, total_time);
    end
    if length(this_del.dim_in) ~= total_time
        assert(length(this_del.dim_in) == 1,...
               'DeltaSltvRateBnd:matchHorizonPeriod',...
               'dim_in of %s is not compatible w/ horizon_period',...
               this_del.name)
        this_del.dim_in = this_del.dim_in * ones(1, total_time);
    end
    if length(this_del.lower_bound) ~= total_time
        assert(length(this_del.lower_bound) == 1,...
               'DeltaSltvRateBnd:matchHorizonPeriod',...
               'lower_bound of %s is not compatible w/ horizon_period',...
               this_del.name)
        this_del.lower_bound = ...
            this_del.lower_bound * ones(1, total_time);
    end
    if length(this_del.upper_bound) ~= total_time
        assert(length(this_del.upper_bound) == 1,...
               'DeltaSltvRateBnd:matchHorizonPeriod',...
               'upper_bound of %s is not compatible w/ horizon_period',...
               this_del.name)
        this_del.upper_bound = ...
            this_del.upper_bound * ones(1, total_time);
    end
    if length(this_del.lower_rate) ~= total_time
        assert(length(this_del.lower_rate) == 1,...
               'DeltaSltvRateBnd:matchHorizonPeriod',...
               'lower_rate of %s is not compatible w/ horizon_period',...
               this_del.name)
        this_del.lower_rate = ...
            this_del.lower_rate * ones(1, total_time);
    end
    if length(this_del.upper_rate) ~= total_time
        assert(length(this_del.upper_rate) == 1,...
               'DeltaSltvRateBnd:matchHorizonPeriod',...
               'upper_rate of %s is not compatible w/ horizon_period',...
               this_del.name)
        this_del.upper_rate = ...
            this_del.upper_rate * ones(1, total_time);
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
    this_del.upper_rate  = this_del.upper_rate(indices);
    this_del.lower_rate  = this_del.lower_rate(indices);
    this_del.horizon_period = new_horizon_period;
    this_del = matchHorizonPeriod(this_del);
end
end

function multiplier = deltaToMultiplier(this_del, varargin)                     %#ok<VANUS,INUSD,STOUT>
%% DELTATOMULTIPLIER function to generate a multiplier from this object. There is no
%  direct conversion from a DeltaSltvRateBnd object to a MultiplierSltvRateBnd object.
%  Instead, one must transform this Delta and the containing LFT to a 
%  DeltaSltvRateBndImpl + modified LFT pair.  A multiplier can then be generated
%  from the DeltaSltvRateBndImpl object
%
%  multiplier = deltaToMultiplier(this_del)
%
%  Variables:
%  ---------
%     Input:
%        this_del : DeltaSltvRateBnd object
%
%  See also DeltaSltvRateBnd, DeltaSltvRateBndImpl.deltaToMultiplier.
    error('DeltaSltvRateBnd:deltaToMultiplier',...
          ['DeltaSltvRateBnd cannot be directly expressed with a',...
           ' multiplier, the LFT with DeltaSltvRateBnd must first be',...
           ' converted to an extended LFT with a',...
           ' DeltaSltvRateBndImpl uncertainty']);
end

function [recastA, recastB, recastC, newDelta] = recastMatricesAndDelta(this_del)
%% RECASTMATRICESANDDELTA function to recast the LFT and DeltaSltvRateBnd such that
%  analysis results of former LFT is equivalent to the recast LFT. This
%  involves replacing DeltaSltvRateBnd with DeltaSltvRateBndImpl, and
%  modifying the LFT to have extra zero inputs
%
%    [newA, newB, newC, newD, newDelta] = recastMatricesAndDelta(this_delta)
%
%    Variables:
%    ---------
%      Input:
%         this_delta : DeltaSltvRateBnd object
%      Output:
%         recastA : function_handle :: function to transform a matrices of LFT
%         recastB : function_handle :: function to transform b matrices of LFT
%         recastC : function_handle :: function to transform c matrices of LFT
%         recastD : function_handle :: function to transform d matrices of LFT
%         recastDelta : DeltaSltvRateBndImpl object :: new Delta object for modified LFT
%
%    See also iqcAnalysis.modifyLft, DeltaSltvRateBndImpl.
    dim_out = num2cell(this_del.dim_out);
    dim_in  = num2cell(this_del.dim_in * (this_del.basis_length - 1));
    recastA = @(a_cell)...
        cellfun(@(a_mat, dim_out, dim_in) [a_mat, zeros(size(a_mat, 1), dim_in)],...
                a_cell,...
                dim_out,...
                dim_in,...
                'UniformOutput', false);
    recastB = [];
    recastC = @(c_cell)...
        cellfun(@(c_mat, dim_out, dim_in) [c_mat, zeros(size(c_mat, 1), dim_in)],...
                c_cell,...
                dim_out,...
                dim_in,...
                'UniformOutput', false);
    total_time = sum(this_del.horizon_period);
    box = cell(1, total_time);
    for i = 1 : total_time
        box{1, i} = [this_del.lower_bound(i), this_del.upper_bound(i);
                     this_del.lower_rate(i),  this_del.upper_rate(i)];
    end
    region_type = 'box';
    newDelta = DeltaSltvRateBndImpl(this_del.name,...
                                              this_del.dim_in,...
                                              region_type,...
                                              box,...
                                              this_del.horizon_period,...
                                              this_del.basis_length);
end

function [del_diff, del_ave, del_scale, del_norm] = normalizeDelta(this_delta)
%% NORMALIZEDELTA function for DeltaSltv object. This extends the default operation
%  defined by the Delta superclass, such that DeltaSltvRateBnd uncertainties
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
rate_scale = nan(1, total_time);
for i = 1:total_time
    diff = (this_delta.upper_bound(i) - this_delta.lower_bound(i)) / 2;
    ave  = (this_delta.upper_bound(i) + this_delta.lower_bound(i)) / 2;
    del_diff{i}   = diff * eye(this_delta.dim_out(i));
    del_ave{i}    = ave * eye(this_delta.dim_out(i));
    del_scale{i}  = eye(this_delta.dim_in(i));
    rate_scale(i) = 1 / diff;
end

% Define normalized Delta
horizon_period = this_delta.horizon_period;
name           = this_delta.name;
dim_outin      = this_delta.dim_in;
lower_bound    = -ones(1, sum(horizon_period));
upper_bound    = ones(1, sum(horizon_period));
lower_rate     = rate_scale .* this_delta.lower_rate;
upper_rate     = rate_scale .* this_delta.upper_rate;
basis_length   = this_delta.basis_length;

del_norm = DeltaSltvRateBnd(name,...
                            dim_outin,...
                            lower_bound,...
                            upper_bound,...
                            lower_rate,...
                            upper_rate,...
                            horizon_period,...
                            basis_length);
end

function value = sample(this_delta, ~)
    %% SAMPLE function for DeltaSltvRateBnd.
    % Currently works for discete-time LTV and continuous-time LTI, but will need modification for continuous-time LTV
    total_time = sum(this_delta.horizon_period);
    a = cell(1, total_time);
    b = cell(1, total_time);
    c = cell(1, total_time);
    d = cell(1, total_time);
    d_mag = (this_delta.upper_bound(1) - this_delta.lower_bound(1))*rand() ...
            + this_delta.lower_bound(1);
    for t = 1 : total_time
        a{t} = zeros(0,0);
        b{t} = zeros(0,this_delta.dim_in(t));
        c{t} = zeros(this_delta.dim_out(t),0);
        change = (this_delta.upper_rate(t) - this_delta.lower_rate(t))*rand() ...
                 + this_delta.lower_rate(t);
        d_mag = d_mag + change;
        d_mag = min(d_mag, this_delta.upper_bound(t));
        d_mag = max(d_mag, this_delta.lower_bound(t));
        d{t} = d_mag * eye(this_delta.dim_out(t));
    end
    value = Ulft(a, b, c, d, {}, 'horizon_period', this_delta.horizon_period);
end

function validateSample(this_delta, value, timestep)
    %% VALIDATESAMPLE function for DeltaSltvRateBnd.
    % Validate base attributes
    validateSample@Delta(this_delta, value, timestep, true);
    % Memoryless
    assert(isempty(value.delta.deltas),...
           'DeltaSltvRateBnd:validateSample',...
           ['Specified value for delta "',this_delta.name,...
           '" must be memoryless.']);
    % Isotropic
    assert(all(cellfun(@(d) isdiag(d) && all(diag(d) == d(1,1)), value.d)),...
           'DeltaSltvRateBnd:validateSample',...
           ['Specified value for delta "',this_delta.name,...
           '" must be a scaling of the identity at all times.']);
    % Bounded
    lb = num2cell(this_delta.lower_bound);
    ub = num2cell(this_delta.upper_bound);
    assert(all(cellfun(@(d, lb, ub) (lb <= d(1,1)) && (d(1,1) <= ub),...
                       value.d, lb, ub)),...
           'DeltaSltvRateBnd:validateSample',...
           ['Specified value for delta "',this_delta.name,...
           '" must be within the bounds at all times.']);
    % Rate bounded (currently continuous-time is always time-invariant)
    if timestep
        if sum(this_delta.horizon_period) > 1
            lr = num2cell(this_delta.lower_rate(1,1:end-1));
            ur = num2cell(this_delta.upper_rate(1,1:end-1));
            assert(all(cellfun(@(d1, d2, lr, ur) (lr <= (d2(1,1)-d1(1,1))) &&...
                                                 ((d2(1,1)-d1(1,1)) <= ur),...
                               value.d(1,1:end-1), value.d(1,2:end), lr, ur)),...
                   'DeltaSltvRateBnd:validateSample',...
                   ['Specified value for delta "',this_delta.name,...
                    '" must never change faster than the rate bounds.']);
            wrap_diff = value.d{end}(1,1) - value.d{value.horizon_period(1)+1}(1,1);
            assert((this_delta.lower_rate(1,end) <= wrap_diff) &&...
                   (wrap_diff <= this_delta.upper_rate(1,end)),...
                   'DeltaSltvRateBnd:validateSample',...
                   ['Specified value for delta "',this_delta.name,...
                   '" exceeds rate bound at period wrap.']);
        end
    end
end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0): Added sample and validateSample methods - Jason Nezvadovitz (jason.nezvadovitz@ll.mit.edu)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)