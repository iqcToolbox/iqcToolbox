classdef MultiplierSltvRepeated < MultiplierDelta
%% MULTIPLIERSLTVREPEATED class for multiple static, linear, time-varying
%  uncertainties (DeltaSltvRepeated). Extends the base class MultiplierDelta.
%
%  extended methods:
%     MultiplierSltvRepeated(delta, varargin) :: Constructor method
%
%  extended properties:
%     names : cell of char arrays :: list of each uncertainty group's name
%     dim_outins : double ::  dimensions of each uncertainty group (square)
%     region_type : string :: type of uncertainty region.  May be 'box', 
%                            'ellipse', or 'polytope'
%     region_data : cell of doubles :: specification of region, which will
%                                      differ by the region_type
%     upper_bounds : double :: upper bounds of different parameters
%     lower_bounds : double :: lower bounds of different parameters
%     ellipses : double :: half-length of axes for centered ellipse
%     vertices : double :: set of vertices describing polytope
%     quad_time_varying : logical :: true if desired to define a
%                         time-varying quadratic decision variable
%
%  See also MultiplierSltvRepeated.MultiplierSltvRepeated

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    names cell
    dim_outins double
    region_type char
    region_data cell
end

properties (Dependent)
    upper_bounds double
    lower_bounds double
    vertices double
    ellipses double
end

properties (SetAccess = immutable)
    quad_time_varying logical
end

methods
function this_mult = MultiplierSltvRepeated(delta, varargin)
%% MULTIPLIERSLTVREPEATED constructor
%
%  this_mult = MultiplierSltvRepeated(delta, 'quad_time_varying', true)
%  this_mult = MutliplierSltvRepeated(delta) assumes the values provided above
%
%  Variables:
%  ---------
%    Input:
%       delta : The delta object characterized by this multiplier
%       quad_time_varying : logical :: Whether or not the decision variable
%                            quad is time-varying (if true, SDP will possibly
%                            be less conservative, but more computationally
%                            challenging)
%    Output:
%       this_mult : MultiplierSltvRepeated object
%
%  See also MultiplierSltvRepeated

% Inputs
input_parser = inputParser;
addRequired(input_parser,...
            'delta',...
            @(del) validateattributes(del, 'DeltaSltvRepeated', {'nonempty'}))
addParameter(input_parser,...
             'quad_time_varying',...
             true,...
             @(quad) validateattributes(quad, 'logical', {'nonempty'}))
addParameter(input_parser,... % This parameter is not used. Only defined for compatibility with other Multiplier constructor calls
             'discrete',...
             true,...
             @(disc) validateattributes(disc, 'logical', {'nonempty'}))

% Parsing inputs
parse(input_parser, delta, varargin{:})
delta               = input_parser.Results.delta;
quad_time_varying   = input_parser.Results.quad_time_varying;

for i = 1 : sum(delta.horizon_period)
    if strcmp(delta.region_type, 'box')
assert(all(abs(delta.upper_bounds{i} + delta.lower_bounds{i}) < 1e-8,'all'),...
       'MultiplierSltvRepeated:MultiplierSltvRepeated',...
       'MultiplierSltvRepeated currently does not support asymetric bounds',...
       ' for DeltaSltvRepeated with "region_type"=="box". (i.e., ',...
       'lower_bound must equal -upper_bound). Normalize')
    end
end

this_mult.name               = delta.name;
this_mult.horizon_period     = delta.horizon_period;
this_mult.names              = delta.names;
this_mult.dim_outins         = delta.dim_outins;
this_mult.region_type        = delta.region_type;
this_mult.region_data        = delta.region_data;
this_mult.quad_time_varying  = quad_time_varying;

% Define filter
total_time      = sum(this_mult.horizon_period);
filter.a        = cell(1, total_time);
filter.b1       = cell(1, total_time);
filter.b2       = cell(1, total_time);
filter.c1       = cell(1, total_time);
filter.c2       = cell(1, total_time);
filter.d11      = cell(1, total_time);
filter.d12      = cell(1, total_time);
filter.d21      = cell(1, total_time);
filter.d22      = cell(1, total_time);
n_dels          = length(this_mult.names);
for i = 1:total_time
    filter.a{i}   = [];
    filter.b1{i}  = zeros(0, delta.dim_in(i));
    filter.b2{i}  = zeros(0, delta.dim_out(i));
    filter.c1{i}  = zeros(delta.dim_in(i), 0);
    filter.c2{i}  = zeros(delta.dim_out(i), 0);
    if strcmp(this_mult.region_type, 'box')
        d11_blocks = cell(1, n_dels);
        for j = 1:n_dels
            d11_blocks{j} = this_mult.upper_bounds{i}(j, 1) *...
                            eye(this_mult.dim_outins(j, i));
        end
        filter.d11{i} = blkdiag(d11_blocks{:});
    else
        filter.d11{i} = eye(delta.dim_in(i));
    end
    filter.d12{i} = zeros(delta.dim_in(i), delta.dim_out(i));
    filter.d21{i} = zeros(delta.dim_out(i), delta.dim_in(i));
    filter.d22{i} = eye(delta.dim_out(i));
end
this_mult.filter = filter;

% Define quadratic, constraints, and decision_vars
quad.q11 = cell(1, total_time);
quad.q12 = cell(1, total_time);
quad.q21 = cell(1, total_time);
quad.q22 = cell(1, total_time);
ct = [];
decision_vars = {};
for i = 1:total_time
    if i == 1 || this_mult.quad_time_varying
        % Defining decision variables and constraints for box regions
        if strcmp(this_mult.region_type, 'box')
            q11j = cell(1, n_dels);
            q12j = cell(1, n_dels);
            for j = 1 : n_dels
                q11j{j} = sdpvar(this_mult.dim_outins(j, i));
                q12j{j} = sdpvar(size(q11j{j}, 1), size(q11j{j}, 2), 'skew');
                ct = ct + ((q11j{j} >= 0):['SLTV Repeated Multiplier, ',...
                                             this_mult.name,...
                                             ', q11 >= 0,',...
                                             ' Time: ', num2str(i),...
                                             ' Subblock: ', num2str(j)]);       %#ok<BDSCA>
            end
            decision_vars = [decision_vars, q11j];
            decision_vars = [decision_vars, q12j];
            q11 = blkdiag(q11j{:});
            q12 = blkdiag(q12j{:});
            q22 = -q11;
        % Defining decision variables for polytope regions
        elseif strcmp(this_mult.region_type, 'polytope')
            q11 = sdpvar(delta.dim_in(i));
            q12 = sdpvar(delta.dim_in(i), delta.dim_out(i));
            q22 = sdpvar(delta.dim_out(i));
            decision_vars = [decision_vars, {q11, q12, q22}];
        elseif strcmp(this_mult.region_type, 'ellipse')
            error('MultiplierSltvRepeated:MultiplierSltvRepeated',...
                  'region_type does not currently support "ellipse"')
        else
            error('MultiplierSltvRepeated:MultiplierSltvRepeated',...
                  'region_type must be "polytope", "box", or "ellipse"')
        end
        
    end
    % Defining constraints for polytope regions
    if strcmp(this_mult.region_type, 'polytope')
        ct = ct + ((q22 <= 0):['SLTV Repeated Multiplier Polytope, ',...
                               this_mult.name,...
                               ' q22 <= 0, ',...
                               'Time: ', num2str(i)]);                          %#ok<BDSCA>
%         I   = eye(delta.dim_in(i));
%         for j = 1 : n_dels
%             start_column = 1 + sum(this_mult.dim_outins(1 : j - 1, i));
%             end_column  = sum(this_mult.dim_outins(1 : j, i));
%             selector = I(:, start_column : end_column);
%             ct = ct + ((selector' * q22 * selector <= 0):...
%                        ['SLTV Repeated Multiplier Polytope, ',...
%                         this_mult.name,...
%                         ', col^T q22 col <= 0, ',...
%                         'Time: ', num2str(i), ' Delta: ', num2str(j)]);
%         end
        n_vertices = size(this_mult.vertices{i}, 2);
        for j = 1 : n_vertices
            vertex = this_mult.vertices{i}(:, j);
            del_vertex = [];
            for k = 1 : n_dels
                del_vertex = ...
                    blkdiag(del_vertex,...
                            vertex(k) * eye(this_mult.dim_outins(k, i)));
            end
            eye_delta = [eye(delta.dim_in(i)); del_vertex];
            q = [q11, q12; q12', q22];
            ct = ct + ((eye_delta' * q * eye_delta >= 0):...
                       ['SLTV Repeated Multiplier Polytope, ',...
                        this_mult.name,...
                        ', [I, del(vertex)] * quad * [I; del(vertex)] >= 0',...
                        ' Time: ', num2str(i), ' Vertex: ', num2str(j)]);       %#ok<BDSCA>
        end
    end
    quad.q11{i} = q11;
    quad.q12{i} = q12;
    quad.q21{i} = q12';
    quad.q22{i} = q22;
end
this_mult.quad = quad;
this_mult.constraints = ct;
this_mult.decision_vars = decision_vars;
end

%% Getter methods for dependent properties
function lower_bounds = get.lower_bounds(this_mult)
    lower_bounds = NaN;
    if strcmp(this_mult.region_type, 'box')
        lower_bounds = cellfun(@(array) array(:, 1),...
                               this_mult.region_data,...
                               'UniformOutput', false);
    end
end

function upper_bounds = get.upper_bounds(this_mult)
    upper_bounds = NaN;
    if strcmp(this_mult.region_type, 'box')
        upper_bounds = cellfun(@(array) array(:, 2),...
                               this_mult.region_data,...
                               'UniformOutput', false);
    end
end

function ellipses = get.ellipses(this_mult)
    ellipses = NaN;
    if strcmp(this_mult.region_type, 'ellipse')
        ellipses = this_mult.region_data;
    end
end

function vertices = get.vertices(this_mult)
    vertices = NaN;
    if strcmp(this_mult.region_type, 'polytope')
        vertices = this_mult.region_data;
    end
end
end

end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)