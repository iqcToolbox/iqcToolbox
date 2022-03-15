classdef DeltaSltvRateBndImpl < Delta
%% DELTASLTVRATEBNDIMPLEMENTATION class for representing the modified Delta
% which is used in IQC analysis for static, linear, time-varying rate-bounded 
% uncertainties, extends the base class Delta. This class should mainly be
% constructed by calling the recastMatricesAndDelta function handle from
% DeltaSltvRateBnd. As in: [~, ~, ~, recastDelta] = recastMatricesAndDelta(d_sltv_rb);
% This class draws significantly from the structure of DeltaSltvRepeated,
% since the rate-of-change can simply be considered as another SLTV parameter.
%
%   extended methods:
%     DeltaSltvRateBndImpl(name, dim_in, region_type, region_data, horizon_period, basis_length) :: Constructor
%     disp(this_delta) :: Display method
%     matchHorizonPeriod(this_delta, new_horizon_period) 
%                        :: Matches delta properties to new horizon_period
%     deltaToMultiplier(this_delta, varargin)
%                        :: Method for constructing a multiplier from a delta
%
%   extended properties:
%     region_type : string :: type of uncertainty region.  May be 'box', 
%                            'ellipse', or 'polytope'
%     region_data : cell of doubles :: specification of region, which will
%                                      differ by the region_type
%     basis_length : double :: length of basis_function for multiplier
%
%   extended dependent properties:
%     upper_bound : double :: upper bound of parameter
%     lower_bound : double :: lower bound of parameter
%     upper_rate : double :: upper bound of parameter's rate-of-change
%     lower_rate : double :: lower bound of parameter's rate-of-change
%     ellipses : double :: half-length of axes for centered ellipse
%     vertices : double :: set of vertices describing polytope
%
%   See also Delta, DeltaSltvRateBnd, DeltaSltvRepeated

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    region_type
    region_data
    basis_length
end

properties (Dependent)
    upper_bound
    lower_bound
    upper_rate
    lower_rate
    vertices
    ellipses
end

methods
function this_delta = DeltaSltvRateBndImpl(name,...
                                           dim_in,...
                                           region_type,...
                                           region_data,...
                                           horizon_period,...
                                           basis_length)
%% DELTASLTVRATEBNDIMPLEMENTATION constructor
%
%  d = DeltaSltvRateBndImpl(name, dim_in, region_type, region_data, horizon_period, basis_length)
%  d = DeltaSltvRateBndImpl(name, dim_in, region_type, region_data, horizon_period) assumes basis_length == 2
%  d = DeltaSltvRateBndImpl(name, dim_in, region_type, region_data) also assumes horizon_period == [0, 1]
%  d = DeltaSltvRateBndImpl(names, dim_outins) also assumes region_type == 'box'
%                                         and region_data == {[-ones(2, 1), ones(2, 1)]}
%                                  (lower_bound and lower_rate of -1, upper_bound and upper_rate of 1)
%  d = DeltaSltvRateBndImpl(names) also assumes dim_outins == ones(2, 1)
%
%  Variables:
%  ---------
%     Input:
%        name : string :: unique ID of the uncertainty (ex. 'dmass')
%        dim_in : natural :: input dimensions of uncertainty
%        region_type : string :: either 'box', 'polytope', or 'ellipse'
%        region_data : cell array of double arrays :: specification of the region.
%                     if region_type == 'box', region_data is a (1 x total_time) cell array of 
%                                              (2 x 2) double arrays, representing at each time 
%                                              the upper_bound {:}(1, 2), the lower_bound {:}(:, 1)
%                                              the upper_rate {:}(2, 2), and the lower_rate {:}(2, 1)
%                                              of the uncertainty throughout time.
%                     if region_type == 'polytope', region_data is a (1 x total_time) cell array
%                                              of (2 x N) arrays, representing at each time
%                                              the N vertices of the convex polytope
%                                              throughout time. The nth component of the ith 
%                                              vertex at the jth time instant is 
%                                              region_type{1, j}(n, i). For example, 
%                                              region_type{1, 1} = [[1, 1], [-1; -1]], 
%                                              region_type{1, 2} = [[0; 0], [0; 0]]
%                                              specifies that (delta1, delta2) is on the line 
%                                              segment from (-1, -1) to (1, 1) at the first time
%                                              instant, and then on the point (0, 0) at the 
%                                              second time instant
%                     if region_type == 'ellipse', region_data is a (1 x total_time) cell array 
%                                              of (2 x 1) arrays, representing at each time
%                                              step the length of the semi-minor and semi-major axes
%                                              of an axis-symmetric, origin-centered ellipse
%        horizon_period : 1 x 2 array of integers :: horizon and period of uncertainty
%        basis_length : natural :: length of multiplier basis
%     Output:
%        this_delta : DeltaSltvRateBndImpl object
%  
%  See also DeltaSltvRateBndImpl, DeltaSltvRepeated.DeltaSltvRepeated

    % Defining defaults for missing arguments
    switch nargin
        case 1
            dim_in = 1;
            region_type = 'box';
            region_data = {[-ones(2, 1), ones(2, 1)]};
            horizon_period = [0, 1];
            basis_length = 2;
        case 2
            region_type = 'box';
            region_data = {[-ones(2, 1), ones(2, 1)]};
            horizon_period = [0, 1];
            basis_length = 2;
        case 4
            horizon_period = [0, 1];
            basis_length = 2;
        case 5
            basis_length = 2;
        case 6
        otherwise
            error(['DeltaSltvRateBndImpl:',...
                   'DeltaSltvRateBndImpl'],...
                  'Incorrect number of arguments given')
    end
    dim_out = dim_in * (basis_length); % dim_in + dim_in * (basis_length - 1)

    % Calling Delta constructor
    this_delta@Delta(name, dim_out, dim_in, horizon_period);                    

    % Checking inputs for specialized properties of
    % DeltaSltvRateBndImpl    
    isValidRegionType = @(str) ischar(str)...
                               && (strcmp(str, 'ellipse')...
                                   || strcmp(str, 'box')...
                                   || strcmp(str, 'polytope'));
    assert(isValidRegionType(region_type),...
           ['DeltaSltvRateBndImpl:',...
            'DeltaSltvRateBndImpl'],...
           'region_type must be "ellipse", "box", or "polytope"')
    switch region_type
        case 'ellipse'
            validateattributes(region_data, {'cell'}, {'size', [1, NaN]})
            for i = 1:length(region_data)
                ellipse = region_data{1, i};
                validateattributes(ellipse, {'numeric'}, {'positive',...
                                                        'size', [2, 1],...
                                                        'nonnan',...
                                                        'finite'})
            end
        case 'box'
            validateattributes(region_data, {'cell'}, {'size', [1, NaN]})
            for i = 1:length(region_data)
                box = region_data{1, i};
                validateattributes(box,...
                                   {'numeric'},...
                                   {'size', [2, 2], 'nonnan', 'finite'})
                assert(all(all(box(:, 1) <= box(:, 2))),...
                   ['DeltaSltvRateBndImpl:'...
                    'DeltaSltvRateBndImpl'],...
                   ['region_data for "box" must have lower bounds (:, 1)',...
                    ' less than upper bounds (:, 2)'])
            end            
        case 'polytope'
            validateattributes(region_data, {'cell'}, {'size', [1, NaN]})
            for i = 1:length(region_data)
                polytope = region_data{1, i};
                validateattributes(polytope, {'numeric'}, {'size', [2, NaN],...
                                                         'nonnan',...
                                                         'finite'})
                [valid_polytope, points] = processConvexHullPoints(polytope);
                assert(valid_polytope,...
                       ['DeltaSltvRateBndImpl:',...
                        'DeltaSltvRateBndImpl'],...
                       ['region_data for "polytope" must have points whose',...
                        ' convex hull includes the origin'])
                region_data{1,i} = points;
            end
    end                       
    this_delta.region_type = region_type;
    this_delta.region_data = region_data;
    this_delta.basis_length = basis_length;
    this_delta = matchHorizonPeriod(this_delta);
end

function disp(this_delta)
%% DISP function for DeltaSltvRateBndImpl object
%
%  disp(delta_sltv_obj) (e.g., disp(DeltaSltvRateBndImpl('d')) )
%
%  Variables:
%  ---------
%     Input:
%        this_delta : DeltaSltvRateBndImpl object    
%
%  See also DeltaSltvRateBndImpl, Delta.disp, SequenceDelta.disp
    disp@Delta(this_delta, 'SLTV, rate-bounded uncertainty (implementation)')
    switch this_delta.region_type
        case 'ellipse'
            fprintf('%13s within a %8s', '   ', this_delta.region_type)
            str = sprintf('%3.1f, ', this_delta.region_data{1, 1});
            fprintf([' having axes: [',str(1 : end - 2),'] \n'])
        case 'box'
            fprintf(['%13s within the bounds: [%3.1f, %3.1f] \n',...
                     '%13s whose rate-of-change is within: [%3.1f, %3.1f] \n',...
                     '%13s with a multiplier basis length of: %3d \n'],...
                    '',...
                    this_delta.lower_bound(1),...
                    this_delta.upper_bound(1),...
                    '',...
                    this_delta.lower_rate(1),...
                    this_delta.upper_rate(1))
        case 'polytope'
            fprintf('%13s within a %8s', '   ', this_delta.region_type)
            fprintf(' described by %2d points having a max 2-norm of %3.1f \n',...
                    size(this_delta.region_data{1, 1}, 2),...
                    max(vecnorm(this_delta.region_data{1, 1}, 2)))
    end 
    fprintf('%13s with a multiplier basis length of: %3d \n',...
            '',...
            this_delta.basis_length)
end

function this_del = matchHorizonPeriod(this_del, new_horizon_period)
%% MATCHHORIZONPERIOD function to ensure properties of DeltaSltvRateBndImpl 
%  object match a total_time length
%
%  matchHorizonPeriod(this_del, new_horizon_period) will change the horizon_period and associated properties of this_del
%  matchHorizonPeriod(this_del) only ensures that each pertinent property of this_del matches this_del.horizon_period
%
%  Variables:
%  ---------
%     Input:
%       this_del : DeltaSltvRateBndImpl object
%       new_horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of new sequence
%     Output:
%       this_delta : DeltaSltvRateBndImpl object
%
%  See also DeltaSltvRateBndImpl

if nargin == 1
% Ensuring that this_del.horizon_period matches with other properties 
% of this_del

    % Assumes properties are a mix of sequences of length 1 or
    % horizon_period
    total_time = sum(this_del.horizon_period);
    if length(this_del.dim_out) ~= total_time
        assert(length(this_del.dim_out) == 1,...
               'DeltaSltvRateBndImpl:matchHorizonPeriod',...
               'dim_out of %s is not compatible w/ horizon_period',...
               this_del.name);
        this_del.dim_out = this_del.dim_out * ones(1, total_time);
    end
    if length(this_del.dim_in) ~= total_time
        assert(length(this_del.dim_in) == 1,...
               'DeltaSltvRateBndImpl:matchHorizonPeriod',...
               'dim_in of %s is not compatible w/ horizon_period',...
               this_del.name)
        this_del.dim_in = this_del.dim_in * ones(1, total_time);
    end
    if length(this_del.region_data) ~= total_time
        assert(length(this_del.region_data) == 1,...
               'DeltaSltvRateBndImpl:matchHorizonPeriod',...
               'region_data of %s is not compatible w/ horizon_period',...
               this_del.name)
        this_del.region_data = repmat(this_del.region_data, 1, total_time);
    end
else
% Changing this_del.horizon_period and other properties of this_del to
% a new horizon_period
    [indices, new_horizon_period] = ...
        makeNewIndices(this_del.horizon_period, new_horizon_period);

    % Set properties according to indices
    this_del.dim_out     = this_del.dim_out(indices);
    this_del.dim_in      = this_del.dim_in(indices);
    this_del.region_data    = this_del.region_data(indices);
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
%        this_del : DeltaSltvRateBndImpl object
%     Output:
%        multiplier : MultiplierSltvRateBndImpl object
%
%  See also DeltaSltvRateBndImpl, DeltaSltvRateBnd.deltaToMultiplier
      multiplier = MultiplierSltvRateBndImpl(this_del, varargin{:});
end

%% Getter methods for dependent properties
function lower_bound = get.lower_bound(this_delta)
    lower_bound = NaN;
    if strcmp(this_delta.region_type, 'box')
        lower_bound = cellfun(@(array) array(1, 1),...
                               this_delta.region_data);
    end
end

function upper_bound = get.upper_bound(this_delta)
    upper_bound = NaN;
    if strcmp(this_delta.region_type, 'box')
        upper_bound = cellfun(@(array) array(1, 2),...
                               this_delta.region_data);
    end
end

function lower_rate = get.lower_rate(this_delta)
    lower_rate = NaN;
    if strcmp(this_delta.region_type, 'box')
        lower_rate = cellfun(@(array) array(2, 1),...
                               this_delta.region_data);
    end
end

function upper_rate = get.upper_rate(this_delta)
    upper_rate = NaN;
    if strcmp(this_delta.region_type, 'box')
        upper_rate = cellfun(@(array) array(2, 2),...
                               this_delta.region_data);
    end
end

function ellipses = get.ellipses(this_delta)
    ellipses = NaN;
    if strcmp(this_delta.region_type, 'ellipse')
        ellipses = this_delta.region_data;
    end
end

function vertices = get.vertices(this_delta)
    vertices = NaN;
    if strcmp(this_delta.region_type, 'polytope')
        vertices = this_delta.region_data;
    end
end

end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)