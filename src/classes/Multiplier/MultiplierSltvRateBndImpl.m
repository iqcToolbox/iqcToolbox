classdef MultiplierSltvRateBndImpl < MultiplierDelta
%% MULTIPLIERSLTVRATEBNDIMPLEMENTATION class for multiple static, linear, time-varying
%  uncertainties (DeltaSltvRateBndImpl). Extends the base class MultiplierDelta.
%
%  extended methods:
%     MultiplierSltvRateBndImpl(delta, varargin) :: Constructor method
%
%  extended properties:
%     region_type : string :: type of uncertainty region.  May be 'box', 
%                            'ellipse', or 'polytope'
%     region_data : cell of doubles :: specification of region, which will
%                                      differ by the region_type
%     dim_out : double :: output dimensions of uncertainty
%     dim_in : double :: input dimensions of uncertainty
%     upper_bound : double :: upper bound of uncertainty
%     lower_bound : double :: lower bound of uncertainty
%     upper_rate : double :: upper bound on rate-of-change of uncertainty
%     lower_rate : double :: lower bound on rate-of-change of uncertainty
%     ellipses : double :: half-length of axes for centered ellipse
%     vertices : double :: set of vertices describing polytope
%     quad_time_varying : logical :: true if desired to define a
%                         time-varying quadratic decision variable
%     discrete : logical :: true if discrete-time, false if continuous-time
%     basis_length : natural :: length of filter's basis_function
%     basis_poles : array of complex doubles :: list of basis_function poles
%     basis_function : tf :: a (basis_length+1 x 1) stable transfer
%                            function for defining filter
%     basis_realization : ss :: the ss realization of basis_function
%     block_realization : struct :: contains fields blk11 and blk22, which
%                             are ss of the upper-left and lower-right
%                             blocks of filter

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    dim_out double
    dim_in double
    region_type char
    region_data cell
    basis_length double
    basis_poles double
    basis_function tf
    basis_realization ss    
end

properties (Dependent)
    upper_bound double
    lower_bound double
    upper_rate double
    lower_rate double
    vertices double
    ellipses double
end

properties (SetAccess = immutable)
    quad_time_varying logical
end

properties (SetAccess = private)
    block_realization struct
end

methods
function this_mult = MultiplierSltvRateBndImpl(delta, varargin)
%% MULTIPLIERSLTVRATEBNDIMPLEMENTATION constructor
%
%  this_mult = MultiplierSltvRateBndImpl(delta, 'quad_time_varying', true, 'discrete', true, 'basis_length', 2, 'basis_pole', -0.5)
%  this_mult = MutliplierSltvRateBndImpl(delta) assumes the values provided above
%
%  Variables:
%  ---------
%    Input:
%       delta : The delta object characterized by this multiplier
%       varargin : Key-value pair of optional arguments for multiplier (see below for more details)
%    Output:
%       this_mult : MultiplierSltvRateBndImpl object
%
%  This constructor has four specification "groups":
%    delta : The delta object characterized by this multiplier
%    discrete : Whether or not this multiplier is in discrete-time
%    constraint : This group contains the following variables:
%       quad_time_varying : logical :: Whether or not the constraint
%                            on quad is time-varying (the former will possibly
%                            be less conservative, but more computationally
%                            challenging)
%    filter_construction : This group contains the following variables:
%       basis_poles : array of complex doubles :: list of basis_function poles
%       basis_function : tf :: a (basis_length+1 x 1) stable transfer
%                              function for defining filter
%       basis_realization : ss :: the ss realization of basis_function
%
%  A MultiplierSltvRateBndImpl can be constructed by providing inputs with varying levels of granularity:
%     least granular --------------------------------------------------- most granular
%        basis_pole    -------->    basis_function   --------->  basis_realization

%  If users provide less granular inputs, those inputs will construct the more granular ones.
%  If users provide more granular inputs, the less granular properties will be left unspecified.
%
%  See also MultiplierSltvRateBndImpl

% Inputs
input_parser = inputParser;
addRequired(input_parser,...
            'delta',...
            @(del) validateattributes(del, 'DeltaSltvRateBndImpl', {'nonempty'}))
addParameter(input_parser,...
             'quad_time_varying',...
             true,...
             @(quad) validateattributes(quad, 'logical', {'nonempty'}))
addParameter(input_parser,...
             'discrete',...
             true,...
             @(disc) validateattributes(disc, 'logical', {'nonempty'}))
addParameter(input_parser,...
             'basis_poles',...
             -0.5,...
             @(bp) isnumeric(bp))
addParameter(input_parser,...
             'basis_function',...
             [],...
             @(bf) isa(bf, 'tf'))
addParameter(input_parser,...
             'basis_realization',...
             [],...
             @(br) isa(br, 'ss'))

% Parsing inputs
parse(input_parser, delta, varargin{:})
delta               = input_parser.Results.delta;
quad_time_varying   = input_parser.Results.quad_time_varying;
discrete            = input_parser.Results.discrete;
basis_poles         = input_parser.Results.basis_poles;
basis_function      = input_parser.Results.basis_function;
basis_realization   = input_parser.Results.basis_realization;
for i = 1 : sum(delta.horizon_period)
    if strcmp(delta.region_type, 'box')
assert(all(abs(delta.upper_bound + delta.lower_bound) < 1e-8,'all'),...
       'MultiplierSltvRateBndImpl:MultiplierSltvRateBndImpl',...
       ['MultiplierSltvRateBndImpl currently does not support',...
       ' asymetric bounds for DeltaSltvRateBndImpl',...
       ' (i.e., lower_bound must equal -upper_bound). Normalize'])
assert(all(abs(delta.upper_rate + delta.lower_rate) < 1e-8,'all'),...
       'MultiplierSltvRateBndImpl:MultiplierSltvRateBndImpl',...
       ['MultiplierSltvRateBndImpl currently does not support',...
        ' asymetric rate-bounds for DeltaSltvRateBndImpl',...
        ' (i.e., lower_rate must equal upper_rate). Normalize'])
    end
end
basis_length = delta.basis_length;

this_mult.name               = delta.name;
this_mult.horizon_period     = delta.horizon_period;
this_mult.dim_out            = delta.dim_out;
this_mult.dim_in             = delta.dim_in;
this_mult.discrete           = discrete;
this_mult.region_type        = delta.region_type;
this_mult.region_data        = delta.region_data;
this_mult.quad_time_varying  = quad_time_varying;
this_mult.basis_length       = basis_length;

if ~isempty(basis_realization)
    this_mult.basis_realization = basis_realization;
    [this_mult.basis_function,...
     this_mult.basis_poles,...
     this_mult.basis_length] = deal([]);
elseif ~isempty(basis_function)
    this_mult.basis_function = basis_function;
    [this_mult.basis_poles,...
     this_mult.basis_length] = deal([]);
else
    if basis_length == 1
        basis_poles = [];
    end
    if isempty(basis_poles)
        basis_length = 1;
    end
    if discrete
        assert(all(abs(basis_poles) < 1, 'all'),...
               'MultiplierSltvRateBndImpl:MultiplierSltvRateBndImpl',...
               'basic_poles must be within the unit circle')
    else
        assert(all(real(basis_poles) < 0, 'all'),...
               'MultiplierSltvRateBndImpl:MultiplierSltvRateBndImpl',...
               'basic_poles must be in the left plane')
    end
    if size(basis_poles, 1) > 1
        assert(size(basis_poles, 1) == basis_length - 1,...
               'MultiplierSltvRateBndImpl:MultiplierSltvRateBndImpl',...
               ['if number of rows in basis_poles is more than one, ',...
                '(i.e., not following repeated pole scheme) ',...
                'its number of rows must be one less than basis_length'])
    end
    if any(~isreal(basis_poles), 'all')
        try 
            cplxpair(basis_poles);
        catch ME
            if strcmp(ME.identifier, 'MATLAB:cplxpair:ComplexValuesPaired')
                assert(false,...
                       'MultiplierSltvRateBndImpl:MultiplierSltvRateBndImpl',...
                       'complex poles must be in conjugate pairs')
            end
        end
    end
    this_mult.basis_poles  = basis_poles;
    this_mult.basis_length = basis_length;
    z_zpk = cell(basis_length, 1);
    p_zpk = cell(basis_length, 1);
    k_zpk = ones(basis_length, 1);
    for i = 2:basis_length
        if size(basis_poles, 1) > 1
            p_zpk{i, 1} = basis_poles(i - 1);
        else
            p_zpk{i, 1} = repmat(basis_poles(1,:), 1, i - 1);
        end
    end
    if discrete
        this_mult.basis_function = zpk(z_zpk, p_zpk, k_zpk, []);
    else
        this_mult.basis_function = zpk(z_zpk, p_zpk, k_zpk);
    end
end
end

%% Setter methods from less granular to more granular

function this_mult = set.basis_function(this_mult, basis_function)
    this_mult.basis_function = basis_function;
    if ~isempty(basis_function)
        assert(ifAndOnlyIf(this_mult.discrete, (basis_function.Ts~=0)),...
               'MultiplierSltvRateBndImpl:set:basis_function',...
               'time domain of basis_function must match multiplier')
        assert(isstable(basis_function),...
               'MultiplierSltvRateBndImpl:set:basis_function',...
               'basis_function must be stable')
        assert(size(basis_function, 2) == 1,...
               'MultiplierSltvRateBndImpl:set:basis_function',...
               'basis_function must be of size ? x 1')                
        this_mult.basis_realization = ss(basis_function, 'minimal');
    end
end

function this_mult = set.basis_realization(this_mult, basis_realization)
    this_mult.basis_realization = basis_realization;
    if ~isempty(basis_realization)
        assert(ifAndOnlyIf(this_mult.discrete, (basis_realization.Ts~=0)),...
               'MultiplierSltvRateBndImpl:set:basis_realization',...
               'time domain of basis_realization must match multiplier')
        assert(isstable(basis_realization),...
               'MultiplierSltvRateBndImpl:set:basis_realization',...
               'basis_realization must be stable')
        assert(size(basis_realization, 2) == 1,...
               'MultiplierSltvRateBndImpl:set:basis_realization',...
               'basis_realization must be of size ? x 1')
        assert(all(this_mult.dim_in == this_mult.dim_in(1)),...
               'MultiplierSltvRateBndImpl:set:basis_realization',...
               ['Cannot make time-invariant basis function if multiplier',...
                ' has time-varying dimensions']);
        basis    = basis_realization;
        dim_in   = this_mult.dim_in(1);
        if this_mult.discrete
            a_tilde = basis.a;
            b_tilde = basis.b;
        else
            a_tilde = eye(size(basis.a, 1));
            b_tilde = zeros(size(basis.b));
        end
        block11 =  ss(kron(basis.a, eye(dim_in)),...
                      kron(basis.b, eye(dim_in)),...
                      kron([basis.c; a_tilde], eye(dim_in)),...
                      kron([basis.d; b_tilde], eye(dim_in)),...
                      basis.Ts);
        block22 =  ss(kron(basis.a, eye(dim_in)),...
                      kron([basis.b, eye(size(basis.b, 1))], eye(dim_in)),...
                      kron([basis.c; zeros(size(basis.a, 2))], eye(dim_in)),... 
                      kron(blkdiag(basis.d, eye(size(basis.a, 1))), eye(dim_in)),...
                      basis.Ts);
        block_realization.block11 = block11;
        block_realization.block22 = block22;
        this_mult.block_realization = block_realization;
    end
end

function this_mult = set.block_realization(this_mult, block_realization)
    br = block_realization;
    blk11 = br.block11;
    blk22 = br.block22;
    this_mult.block_realization = br;

    % Define filter
    total_time = sum(this_mult.horizon_period);
    filter.a   = cell(1, total_time);
    filter.b1  = cell(1, total_time);
    filter.b2  = cell(1, total_time);
    filter.c1  = cell(1, total_time);
    filter.c2  = cell(1, total_time);
    filter.d11 = cell(1, total_time);
    filter.d12 = cell(1, total_time);
    filter.d21 = cell(1, total_time);
    filter.d22 = cell(1, total_time);
    [dim_out1, dim_in1] = size(blk11.d);
    [dim_out2, dim_in2] = size(blk22.d);
    dim_state1 = size(blk11.a, 1);
    dim_state2 = size(blk22.a, 1);
    for i = 1:total_time
        filter.a{i}   = blkdiag(blk11.a, blk22.a);
        filter.b1{i}  = [blk11.b; 
                         zeros(dim_state2, dim_in1)];
        filter.b2{i}  = [zeros(dim_state1, dim_in2); 
                         blk22.b];
        filter.c1{i}  = [blk11.c,         zeros(dim_out1, dim_state2)];
        filter.c2{i}  = [zeros(dim_out2, dim_state1),         blk22.c];
        filter.d11{i} = blk11.d;
        filter.d12{i} = zeros(dim_out1, dim_in2);
        filter.d21{i} = zeros(dim_out2, dim_in1);
        filter.d22{i} = blk22.d;
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
                q11_bnd = sdpvar(this_mult.dim_in(i) * this_mult.basis_length);
                q12_bnd = sdpvar(size(q11_bnd, 1), size(q11_bnd, 2), 'skew');
                ct = ct + ((q11_bnd >= 0):['SLTV-RB Multiplier, ',...
                                          this_mult.name,...
                                          ', q11 >= 0,',...
                                          ' Time: ', num2str(i),...
                                          ', magnitude bound']);                %#ok<BDSCA>
                q11_rte = sdpvar(this_mult.dim_in(i) * ...
                                 (this_mult.basis_length - 1));
                q12_rte = sdpvar(size(q11_rte, 1), size(q11_rte, 2), 'skew');
                ct = ct + ((q11_rte >= 0):['SLTV-RB Multiplier, ',...
                                          this_mult.name,...
                                          ', q11 >= 0,',...
                                          ' Time: ', num2str(i),...
                                          ', rate bound']);                     %#ok<BDSCA>
                decision_vars = [decision_vars,...
                                 {q11_bnd, q12_bnd, q11_rte, q12_rte}];
                q11 = blkdiag(this_mult.upper_bound(i) ^ 2 * q11_bnd,...
                              this_mult.upper_rate(i) ^ 2 * q11_rte);
                q12 = blkdiag(q12_bnd, q12_rte);
                q22 = -blkdiag(q11_bnd, q11_rte);
            % Defining decision variables for polytope regions
            elseif strcmp(this_mult.region_type, 'polytope')
                q11 = sdpvar(size(this_mult.block_realization.block11, 1));
                q22 = sdpvar(size(this_mult.block_realization.block22, 1));
                q12 = sdpvar(size(q11, 1), size(q22, 1), 'full');
                decision_vars = [decision_vars, {q11, q12, q22}];
            elseif strcmp(this_mult.region_type, 'ellipse')
                error('MultiplierSltvRateBndImpl:set:block_realization',...
                      'region_type does not currently support "ellipse"')
            else
                error('MultiplierSltvRateBndImpl:set:block_realization',...
                      'region_type must be "polytope", "box", or "ellipse"')
            end

        end
        % Defining polytope constraints that must hold at each timestep
        if strcmp(this_mult.region_type, 'polytope')
            dim_delta_mag = this_mult.basis_length * this_mult.dim_in;
            ct = ct + ((q22(1:dim_delta_mag, 1:dim_delta_mag) <= 0):...
                       ['SLTV-RB Multiplier Polytope, ',...
                            this_mult.name,...
                            ', q22_bound <= 0, ',...
                            'Time: ', num2str(i)]);                             %#ok<BDSCA>
            ct = ct + ((q22(dim_delta_mag+1:end,dim_delta_mag+1:end) <= 0):...
                       ['SLTV-RB Multiplier Polytope, ',...
                            this_mult.name,...
                            ', q22_rate <= 0, ',...
                            'Time: ', num2str(i)]);                             %#ok<BDSCA>
            n_vertices = size(this_mult.vertices{i}, 2);
            for j = 1 : n_vertices
                vertex = this_mult.vertices{i}(:, j);
                del_mag = vertex(1) * ...
                    eye(this_mult.basis_length * this_mult.dim_in);
                del_rate= vertex(2) * ...
                    eye((this_mult.basis_length - 1) * this_mult.dim_in);
                del_vertex = blkdiag(del_mag, del_rate);
                eye_delta = [eye(size(del_vertex, 1)); del_vertex];
                q = [q11, q12; q12', q22];
                ct = ct + ((eye_delta' * q * eye_delta >= 0):...
                           ['SLTV-RB Multiplier Polytope, ',...
                            this_mult.name,...
                            ', [I, del(vertex)] * quad * [I; del(vertex)] >= 0',...
                            ' Time: ', num2str(i), ' Vertex: ', num2str(j)]);   %#ok<BDSCA>
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
function lower_bound = get.lower_bound(this_mult)
    lower_bound = NaN;
    if strcmp(this_mult.region_type, 'box')
        lower_bound = cellfun(@(array) array(1, 1),...
                               this_mult.region_data);
    end
end

function upper_bound = get.upper_bound(this_mult)
    upper_bound = NaN;
    if strcmp(this_mult.region_type, 'box')
        upper_bound = cellfun(@(array) array(1, 2),...
                               this_mult.region_data);
    end
end

function lower_rate = get.lower_rate(this_mult)
    lower_rate = NaN;
    if strcmp(this_mult.region_type, 'box')
        lower_rate = cellfun(@(array) array(2, 1),...
                               this_mult.region_data);
    end
end

function upper_rate = get.upper_rate(this_mult)
    upper_rate = NaN;
    if strcmp(this_mult.region_type, 'box')
        upper_rate = cellfun(@(array) array(2, 2),...
                               this_mult.region_data);
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