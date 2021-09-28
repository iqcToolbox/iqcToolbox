classdef DeltaSltvRepeated < Delta
%% DELTASLTVREPEATED class for multiple static, linear, time-varying 
% uncertainties/parameters, extends the base class Delta.
%
%   extended methods:
%     DeltaSltvRepeated(names, dim_outins, region_type, region_data, horizon_period) :: Constructor
%     disp(this_delta) :: Display method
%
%   extended properties:
%     names : cell of char arrays :: list of each uncertainty group's name
%     dim_outins : double ::  dimensions of each uncertainty group (square)
%     region_type : string :: type of uncertainty region.  May be 'box', 
%                            'ellipse', or 'polytope'
%     region_data : cell of doubles :: specification of region, which will
%                                      differ by the region_type
%
%   extended dependent properties:
%     upper_bounds : double :: upper bounds of different parameters
%     lower_bounds : double :: lower bounds of different parameters
%     ellipses : double :: half-length of axes for centered ellipse
%     vertices : double :: set of vertices describing polytope
%
%   See also Delta, DeltaSltv, DeltaSltvRateBndImpl

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    names
    dim_outins
    region_type
    region_data
end

properties (Dependent)
    lower_bounds
    upper_bounds
    vertices
    ellipses    
end

methods
function this_delta = DeltaSltvRepeated(names,...
                                        dim_outins,...
                                        region_type,...
                                        region_data,...
                                        horizon_period)
%% DELTASLTVREPEATED constructor
%
%  d = DeltaSltvRepeated(names, dim_outins, region_type, region_data, horizon_period)
%  d = DeltaSltvRepeated(names, dim_outins, region_type, region_data) assumes horizon_period == [0, 1]
%  d = DeltaSltvRepeated(names, dim_outins) also assumes region_type == 'box'
%                                         and region_data == {[-ones(n, 1), ones(n, 1)]}
%                                  (lower_bounds of -1, upper_bounds of 1)
%  d = DeltaSltvRepeated(names) also assumes dim_outins == ones(n, 1)
%  n in the previous calls refers to the number of repeated Deltas (i.e., the length of the 'names' argument)
%
%  Variables:
%  ---------
%     Input:
%        names : cell array of strings :: unique IDs of the uncertainties (ex. 'dmass')
%        dim_outins : array of naturals :: output/input dimensions of uncertainties
%        region_type : string :: either 'box', 'polytope', or 'ellipse'
%        region_data : cell array of double arrays :: specification of the region.
%                     if region_type == 'box', region_data is a (1 x total_time) cell array of 
%                                              (n x 2) double arrays, representing at each time 
%                                              the upper_bounds {:}(:, 2) of the n 
%                                              uncertainties through time, and the lower_bounds 
%                                              {:}(:, 1) of the n uncertainties throughout time.
%                     if region_type == 'polytope', region_data is a (1 x total_time) cell array
%                                              of (n x N) arrays, representing at each time
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
%                                              of (n x 1) arrays, representing at each time
%                                              step the axis length of
%                                              the nth component of an axis-symmetric, 
%                                              origin-centered ellipse
%        horizon_period : 1 x 2 array of integers :: horizon and period of uncertainty
%     Output:
%        this_delta : DeltaSltvRepeated object
%  
%  See also DeltaSltvRepeated, DeltaSltvRateBndImpl.DeltaSltvRateBndImpl

    % Defining defaults for missing arguments
    if nargin == 0
        error('DeltaSltvRepeated:DeltaSltvRepeated',...
              ['Must provide at least one input to specify the',...
               ' name/s of the uncertainty/uncertainties'])
    else
    % Parsing 'names'
        isValidStr = @(str) ischar(str) && ~isempty(str) && ~contains(str, ',');
        isValidCellOfStrs = @(strs) iscell(strs) ...
                                    ...&& size(strs, 2) == 1 ...
                                    && all(cellfun(isValidStr, strs));
        if isValidStr(names)
            name = names;
            names = {name};
        elseif isValidCellOfStrs(names)
            name = strjoin(names, ', ');
        else
            error('DeltaSltvRepeated:DeltaSltvRepeated',...
                  ['the argument "name" must be a string for a single',...
                   ' uncertainty, or a (n x 1) cell of strings for multiple'])
        end
        n = length(names); % Number of uncertainties/parameters
    end
    switch nargin
        case 1
            dim_outins = ones(n, 1);
            region_type = 'box';
            region_data = {[-ones(n, 1), ones(n, 1)]};
            horizon_period = [0, 1];
        case 2
            region_type = 'box';
            region_data = {[-ones(n, 1), ones(n, 1)]};
            horizon_period = [0, 1];
        case 4
            horizon_period = [0, 1];
    end

    % Calling Delta constructor
    this_delta@Delta(name,...
                     sum(dim_outins, 1),...
                     sum(dim_outins, 1),...
                     horizon_period);                    

    % Checking inputs for specialized properties of DeltaSltv
    validateattributes(dim_outins, 'numeric', {'nonnegative',...
                                               'integer',...
                                               'size', [n, NaN]})
    assert(all(any(dim_outins, 2)),...
           'DeltaSltvRepeated:DeltaSltvRepeated',...
           'Dimensions of each uncertainty must be non-zero at least once')
    isValidRegionType = @(str) ischar(str)...
                               && (strcmp(str, 'ellipse')...
                                   || strcmp(str, 'box')...
                                   || strcmp(str, 'polytope'));
    assert(isValidRegionType(region_type),...
           'DeltaSltvRepeated:DeltaSltvRepeated',...
           'region_type must be "ellipse", "box", or "polytope"')
    switch region_type
        case 'ellipse'
            validateattributes(region_data, 'cell', {'size', [1, NaN]})
            for i = 1:length(region_data)
                ellipse = region_data{1, i};
                validateattributes(ellipse, 'numeric', {'positive',...
                                                        'size', [n, 1],...
                                                        'nonnan',...
                                                        'finite'})
            end
        case 'box'
            validateattributes(region_data, 'cell', {'size', [1, NaN]})
            for i = 1:length(region_data)
                box = region_data{1, i};
                validateattributes(box,...
                                   'numeric',...
                                   {'size', [n, 2], 'nonnan', 'finite'})
                assert(all(box(:, 1) <= box(:, 2), 'all'),...
                   'DeltaSltvRepeated:DeltaSltvRepeated',...
                   ['region_data for "box" must have lower bounds (:, 1)',...
                    ' less than upper bounds (:, 2)'])
            end
            
        case 'polytope'
            validateattributes(region_data, 'cell', {'size', [1, NaN]})
            for i = 1:length(region_data)
                polytope = region_data{1, i};
                validateattributes(polytope, 'numeric', {'size', [n,NaN],...
                                                         'nonnan',...
                                                         'finite'})
                [valid_polytope, points] = processConvexHullPoints(polytope);
                assert(valid_polytope,...
                       'DeltaSltvRepeated:DeltaSltvRepeated',...
                       ['region_data for "polytope" must have points whose',...
                        ' convex hull includes the origin'])
                region_data{1,i} = points;
            end
    end
    
    this_delta.names = names;
    this_delta.dim_outins = dim_outins;
    this_delta.region_type = region_type;
    this_delta.region_data = region_data;
    this_delta = matchHorizonPeriod(this_delta);
end

function disp(this_delta)
%% DISP function for DeltaSltvRepeated object
%
%  disp(delta_sltv_obj) (e.g., disp(DeltaSltv('d')) )
%
%  Variables:
%  ---------
%     Input:
%        this_delta : DeltaSltvRepeated object        
%
%  See also DeltaSltvRepeated, Delta.disp, SequenceDelta.disp
    disp@Delta(this_delta, 'repeated SLTV uncertainty')
    str = sprintf(' %3d, ', this_delta.dim_outins(:, 1));
    fprintf(['%13s with repetitions: [', str(1 : end - 2), '] \n'], '')
    fprintf('%13s within a %8s', '   ', this_delta.region_type)
    switch this_delta.region_type
        case 'ellipse'
            str = sprintf('%3.1f, ', this_delta.region_data{1, 1});
            fprintf([' having axes: [',str(1 : end - 2),'] \n'])
        case 'box'
            str = sprintf('[ %3.1f, %3.1f ], ', this_delta.region_data{1, 1}');
            fprintf([' having bounds: ', str(1 : end - 2), '\n'])
        case 'polytope'
            fprintf(' described by %2d points having a max 2-norm of %3.1f \n',...
                    size(this_delta.region_data{1, 1}, 2),...
                    max(vecnorm(this_delta.region_data{1, 1}, 2)))
    end            
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
%       this_del : DeltaSltvRepeated object
%       new_horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of new sequence
%     Output:
%       this_delta : DeltaSltvRepeated object
%
%  See also DeltaSltvRepeated

if nargin == 1
% Ensuring that this_del.horizon_period matches with other properties 
% of this_del

    % Assumes properties are a mix of sequences of length 1 or
    % horizon_period
    total_time = sum(this_del.horizon_period);
    if length(this_del.dim_out) ~= total_time
        assert(length(this_del.dim_out) == 1,...
               'DeltaSltvRepeated:matchHorizonPeriod',...
               'dim_out of %s is not compatible w/ horizon_period',...
               this_del.name);
        this_del.dim_out = this_del.dim_out * ones(1, total_time);
    end
    if length(this_del.dim_in) ~= total_time
        assert(length(this_del.dim_in) == 1,...
               'DeltaSltvRepeated:matchHorizonPeriod',...
               'dim_in of %s is not compatible w/ horizon_period',...
               this_del.name)
        this_del.dim_in = this_del.dim_in * ones(1, total_time);
    end
    if size(this_del.dim_outins, 2) ~= total_time
        assert(size(this_del.dim_outins, 2) == 1,...
               'DeltaSltvRepeated:matchHorizonPeriod',...
               'dim_outins of %s is not compatible w/ horizon_period',...
               this_del.name)
        this_del.dim_outins = ...
            this_del.dim_outins * ones(1, total_time);
    end
    if length(this_del.region_data) ~= total_time
        assert(length(this_del.region_data) == 1,...
               'DeltaSltvRepeated:matchHorizonPeriod',...
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
    this_del.dim_out        = this_del.dim_out(indices);
    this_del.dim_in         = this_del.dim_in(indices);
    this_del.dim_outins     = this_del.dim_outins(:, indices);
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
%       this_del : DeltaSltvRepeated object
%     Output:
%       multiplier : MultiplierSltvRepeated object
%
%  See also DeltaSltvRepeated
    multiplier = MultiplierSltvRepeated(this_del, varargin{:});
end

function [del_diff, del_ave, del_scale, del_norm] = normalizeDelta(this_del)
%% NORMALIZEDELTA function for DeltaSltvRepeated object. This extends the 
%  default operation defined by the Delta superclass, such that DeltaSltv 
%  uncertainties are normalized
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

if strcmp(this_del.region_type, 'box')
    % Define operations to modify a, b, c, d matrices of LFT
    total_time = sum(this_del.horizon_period);
    del_diff   = cell(1, total_time);
    del_ave    = cell(1, total_time);
    del_scale  = cell(1, total_time);
    num_params = length(this_del.names);
    for i = 1:total_time
        % Loop over different repeated parameters
        vec_diff = [];
        vec_ave  = [];
        for j = 1:num_params
            diff = (this_del.upper_bounds{i}(j)-this_del.lower_bounds{i}(j))/2;
            ave  = (this_del.upper_bounds{i}(j)+this_del.lower_bounds{i}(j))/2;
            vec_diff = [vec_diff; diff * ones(this_del.dim_outins(j, i), 1)];
            vec_ave  = [vec_ave; ave * ones(this_del.dim_outins(j, i), 1)];
        end
        del_diff{i}  = diag(vec_diff);
        del_ave{i}   = diag(vec_ave);
        del_scale{i} = eye(this_del.dim_in(i));
    end

    % Define normalized Delta
    horizon_period = this_del.horizon_period;
    names          = this_del.names;
    dim_outins     = this_del.dim_outins;
    type           = 'box';
    lower_bnd      = -ones(num_params, 1);
    upper_bnd      = ones(num_params, 1);
    data           = repmat({[lower_bnd, upper_bnd]}, 1, total_time);
    
    del_norm = DeltaSltvRepeated(names, dim_outins, type, data, horizon_period);
elseif strcmp(this_del.region_type, 'polytope')
    [del_diff, del_ave, del_scale] = normalizeDelta@Delta(this_del);
    del_norm = this_del;
    warning('DeltaSltvRepeated:normalizeDelta',...
            ['Delta "%s" is not normalized, because it has a "polytope"',...
             ' region type, which does not support normalization'],...
            this_del.name)    
elseif strcmp(this_del.region_type, 'ellipse')
    [del_diff, del_ave, del_scale] = normalizeDelta@Delta(this_del);
    del_norm = this_del;
    warning('DeltaSltvRepeated:normalizeDelta',...
            ['Delta "%s" is not normalized, because it has an "ellipse"',...
             ' region type, which does not support normalization'],...
            this_del.name) 
end
end

%% Getter methods for dependent properties
function lower_bounds = get.lower_bounds(this_delta)
    lower_bounds = NaN;
    if strcmp(this_delta.region_type, 'box')
        lower_bounds = cellfun(@(array) array(:, 1),...
                               this_delta.region_data,...
                               'UniformOutput', false);
    end
end

function upper_bounds = get.upper_bounds(this_delta)
    upper_bounds = NaN;
    if strcmp(this_delta.region_type, 'box')
        upper_bounds = cellfun(@(array) array(:, 2),...
                               this_delta.region_data,...
                               'UniformOutput', false);
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

function value = sample(this_del, ~)
    %% SAMPLE function for DeltaSltvRepeated.
    assert(...
        ~isnan(this_del.upper_bounds{1}(1)),...
        'DeltaSltvRepeated:sample',...
        'Sampling currently only supports box bounds.')
    components = cell(1, size(this_del.dim_outins, 1));
    total_time = sum(this_del.horizon_period);
    for i = 1 : length(components)
        a = cell(1, total_time);
        b = cell(1, total_time);
        c = cell(1, total_time);
        d = cell(1, total_time);
        for t = 1 : total_time
            a{t} = zeros(0,0);
            b{t} = zeros(0,this_del.dim_outins(i));
            c{t} = zeros(this_del.dim_outins(i),0);
            range = this_del.upper_bounds{t}(i) - this_del.lower_bounds{t}(i);
            d_mag = range * rand() + this_del.lower_bounds{t}(i);
            d{t} = d_mag * eye(this_del.dim_outins(i));
        end
        components{i} = Ulft(a, b, c, d, {}, 'horizon_period', this_del.horizon_period);
    end
    value = blkdiag(components{:});
end

function validateSample(this_delta, value, ~)
    %% VALIDATESAMPLE function for DeltaSltvRepeated.
    % Validate base attributes
    validateSample@Delta(this_delta, value, [], true);
    % Validate each DeltaSltv component block
    I = eye(this_delta.dim_out(1));
    p = 1;
    i = 1;
    for dim = this_delta.dim_outins'
        assert(...
            all(this_delta.dim_outins(i,:) == dim(1)),...
            'DeltaSltvRepeated:validateSample',...
            'DeltaSltvRepeated does not validate samples with time varying dims.');
        sel = I(:, p:p+dim-1);
        selft = toLft(sel).matchHorizonPeriod(this_delta.horizon_period);
        selft_tr = toLft(sel').matchHorizonPeriod(this_delta.horizon_period);
        sltv = DeltaSltv(this_delta.names{i},...
                         dim',...
                         cellfun(@(lb) lb(i), this_delta.lower_bounds),...
                         cellfun(@(ub) ub(i), this_delta.upper_bounds),...
                         this_delta.horizon_period);
        sltv.validateSample(selft_tr*value*selft);
        p = p + dim;
        i = i + 1;
    end
end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0): Added sample and validateSample methods - Jason Nezvadovitz (jason.nezvadovitz@ll.mit.edu)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)