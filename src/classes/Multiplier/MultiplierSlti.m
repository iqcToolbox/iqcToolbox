classdef MultiplierSlti < MultiplierDelta
%% MULTIPLIERSLTI class for static, linear, time-invariant uncertainties 
%  (DeltaSlti).
%  Extends the base class MultiplierDelta.
%
%  extended methods:
%     MultiplierSlti(delta, varargin) :: Constructor method
%
%  extended properties:
%     discrete : logical :: true if discrete-time, false if continuous-time
%     basis_length : natural :: length of filter's basis_function
%     basis_poles : array of complex doubles :: list of basis_function poles
%     basis_function : tf :: a (basis_length+1 x 1) stable transfer
%                            function for defining filter
%     basis_realization : ss :: the ss realization of basis_function
%     block_realization : ss :: the ss of the upper-left block defining
%                               filter = [bound * blk_realzation, 0;
%                                                  0             , blk_realization];
%     upper_bound : double :: upper bound of uncertainty
%     dim_outin : double :: dimensions of uncertainty
%     constraint_q11_kyp : logical :: true if a kyp-based constraint for
%                                     quad.q11 is desired
%     constraint_q12_kyp : logical :: true if a kyp-based constraint for
%                                     quad.q12 is desired
%
%  See also MultiplierSlti.MultiplierSlti

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    basis_length double
    basis_poles double
    basis_function tf
    basis_realization ss
    block_realization ss
    upper_bound double
    dim_outin double
end

properties (SetAccess = immutable)
    constraint_q11_kyp logical
    constraint_q12_kyp logical
end

methods
function this_mult = MultiplierSlti(delta, varargin)
%% MULTIPLIERSLTI constructor
%
%  this_mult = MultiplierSlti(delta, 'discrete', true, 'basis_length', 2, 'basis_pole', -0.5, 'constraint_q11_kyp', true, 'constraint_q12_kyp', false)
%  this_mult = MutliplierSlti(delta) assumes the values provided above
%
%  Variables:
%  ---------
%    Input:
%       delta : The delta object characterized by this multiplier
%       varargin : Key-value pair of optional arguments for multiplier (see below for more details)
%    Output:
%       this_mult : MultiplierSlti object
%
%  This constructor has four specification "groups":
%    delta : The delta object characterized by this multiplier
%    discrete : Whether or not this multiplier is in discrete-time
%    constraint : This group contains the following variables:
%       constraint_q11_kyp : logical :: Whether or not the constraint
%                            on quad.q11 will involve the kyp lemma
%                            (the former will possibly be less
%                            conservative, but more computationally
%                            challenging)
%       constraint_q12_kyp : logical :: Whether or not the constraint
%                            on quad.q12 will involve the kyp lemma
%    filter_construction : This group contains the following variables:
%       basis_length : natural :: length of filter's basis_function
%       basis_poles : array of complex doubles :: list of basis_function poles
%       basis_function : tf :: a (basis_length+1 x 1) stable transfer
%                              function for defining filter
%       basis_realization : ss :: the ss realization of basis_function
%       block_realization : ss :: the ss of the upper-left/lower-right block defining
%                                 filter = [bound * blk_realzation, 0;
%                                                  0             , blk_realization];    
%
%  A MultiplierSlti can be constructed by providing inputs with varying levels of granularity:
%     least granular --------------------------------------------------- most granular
%        [basis_length ->  basis_function -> basis_realization -> block_realization
%         basis_pole]
%  If users provide less granular inputs, those inputs will construct the more granular ones.
%  If users provide more granular inputs, the less granular properties will be left unspecified.
%
%  See also MultiplierSlti

% Inputs
input_parser = inputParser;
addRequired(input_parser,...
            'delta',...
            @(del) validateattributes(del, {'DeltaSlti'}, {'nonempty'}))
addParameter(input_parser,...
             'discrete',...
             true,...
             @(disc) validateattributes(disc, {'logical'}, {'nonempty'}))
addParameter(input_parser,...
             'constraint_q11_kyp',...
             true,...
             @(constraint) islogical(constraint))
addParameter(input_parser,...
             'constraint_q12_kyp',...
             false,...
             @(constraint) islogical(constraint))
addParameter(input_parser,...
             'basis_length',...
             2,...
             @(bl) validateattributes(bl,...
                                      {'numeric'},...
                                      {'integer', 'positive', 'nonempty'}))
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
addParameter(input_parser,...
             'block_realization',...
             [],...
             @(block) isa(block, 'ss'))
addParameter(input_parser,... % This parameter is not used. Only defined for compatibility with other Multiplier constructor calls
             'exponential',...
             [],...
             @(disc) validateattributes(disc, {'double'}, {'finite'})) 

% Parsing inputs
parse(input_parser, delta, varargin{:})
delta               = input_parser.Results.delta;
discrete            = input_parser.Results.discrete;
constraint_q11_kyp  = input_parser.Results.constraint_q11_kyp;
constraint_q12_kyp  = input_parser.Results.constraint_q12_kyp;
basis_length        = input_parser.Results.basis_length;
basis_poles         = input_parser.Results.basis_poles;
basis_function      = input_parser.Results.basis_function;
basis_realization   = input_parser.Results.basis_realization;
block_realization   = input_parser.Results.block_realization;

assert(all(abs(delta.upper_bound + delta.lower_bound) < 1e-8),...
       'MultiplierSlti:MultiplierSlti',...
       'MultiplierSlti currently does not support asymmetric bounds',...
       ' on DeltaSlti (lower_bound must equal -upper_bound). Normalize')
this_mult.name               = delta.name;
this_mult.horizon_period     = delta.horizon_period;
this_mult.upper_bound        = delta.upper_bound;  
this_mult.dim_outin          = delta.dim_out;
this_mult.discrete           = discrete;
this_mult.constraint_q11_kyp = constraint_q11_kyp;
this_mult.constraint_q12_kyp = constraint_q12_kyp;

if ~isempty(block_realization)
    this_mult.block_realization = block_realization;
    [this_mult.basis_realization,...
     this_mult.basis_function,...
     this_mult.basis_poles,...
     this_mult.basis_length] = deal([]);
elseif ~isempty(basis_realization)
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
        assert(all(all(abs(basis_poles) < 1)),...
               'MultiplierSlti:MultiplierSlti',...
               'basis_poles must be within the unit circle')
    else
        assert(all(all(real(basis_poles) < 0)),...
               'MultiplierSlti:MultiplierSlti',...
               'basis_poles must be in the left plane')
    end
    if size(basis_poles, 1) > 1
        assert(size(basis_poles, 1) == basis_length - 1,...
               'MultiplierSlti:MultiplierSlti',...
               ['if number of rows in basis_poles is more than one, ',...
               '(i.e., not following repeated pole scheme) ',...
               'its number of rows must be one less than basis_length'])
    end
    if any(any(~isreal(basis_poles)))
        try 
            cplxpair(basis_poles);
        catch ME
            if strcmp(ME.identifier, 'MATLAB:cplxpair:ComplexValuesPaired')
                assert(false,...
                       'MultiplierSlti:MultiplierSlti',...
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
        % Size is guaranteed to align with basis_length because of previous assert
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
               'MultiplierSlti:set:basis_function',...
               'time domain of basis_function must match multiplier')
        assert(isstable(basis_function),...
               'MultiplierSlti:set:basis_function',...
               'basis_function must be stable')
        assert(size(basis_function, 2) == 1,...
               'MultiplierSlti:set:basis_function',...
               'basis_function must be of size ? x 1')                
        this_mult.basis_realization = ss(basis_function, 'minimal');
    end
end

function this_mult = set.basis_realization(this_mult, basis_realization)
    this_mult.basis_realization = basis_realization;
    if ~isempty(basis_realization)
        assert(ifAndOnlyIf(this_mult.discrete, (basis_realization.Ts~=0)),...
               'MultiplierSlti:set:basis_realization',...
               'time domain of basis_realization must match multiplier')
        assert(isstable(basis_realization),...
               'MultiplierSlti:set:basis_realization',...
               'basis_realization must be stable')
        assert(size(basis_realization, 2) == 1,...
               'MultiplierSlti:set:basis_realization',...
               'basis_realization must be of size ? x 1')
        basis = basis_realization;
        dim_out = this_mult.dim_outin(1);
        this_mult.block_realization = ss(kron(basis.a, eye(dim_out)),...
                                         kron(basis.b, eye(dim_out)),...
                                         kron(basis.c, eye(dim_out)),...
                                         kron(basis.d, eye(dim_out)),...
                                         basis.Ts);
    end
end

function this_mult = set.block_realization(this_mult, block_realization)
    br = block_realization;
    assert(ifAndOnlyIf(this_mult.discrete, (block_realization.Ts~=0)),...
               'MultiplierSlti:set:block_realization',...
               'time domain of block_realization must match multiplier')
    assert(isstable(block_realization),...
           'MultiplierSlti:set:block_realization',...
           'block_realization must be stable')
    assert(size(block_realization, 2) == this_mult.dim_outin(1),...
           'MultiplierSlti:set:block_realization',...
           'block_realization must be of size ? x delta.dim_out')
    this_mult.block_realization = br;

    % Define filter and quadratic
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
    quad.q11 = cell(1, total_time);
    quad.q12 = cell(1, total_time);
    quad.q21 = cell(1, total_time);
    quad.q22 = cell(1, total_time);

    [dim_out1, dim_in1] = size(br.d);
    dim_state1 = size(br.a, 1);
    dim_state2 = dim_state1;
    dim_out2   = dim_out1;
    dim_in2    = dim_in1;
    q11 = sdpvar(dim_out1);
    if this_mult.constraint_q12_kyp
        q12 = sdpvar(dim_out1, dim_out2, 'full');
    else
        q12 = sdpvar(dim_out1, dim_out2, 'skew');
    end
    for i = 1:total_time
        filter.a{i}   = blkdiag(br.a, br.a);
        filter.b1{i}  = [br.b; zeros(dim_state2, dim_in1)];
        filter.b2{i}  = [zeros(dim_state1, dim_in2); br.b];
        filter.c1{i}  = this_mult.upper_bound(1) * ...
                        [br.c, zeros(dim_out1, dim_state2)];
        filter.c2{i}  = [zeros(dim_out2, dim_state1), br.c];
        filter.d11{i} = this_mult.upper_bound(1) * br.d;
        filter.d12{i} = zeros(dim_out1, dim_in2);
        filter.d21{i} = zeros(dim_out2, dim_in1);
        filter.d22{i} = br.d;

        quad.q11{i} = q11;
        quad.q12{i} = q12;
        quad.q21{i} = q12';
        quad.q22{i} = -q11;
    end
    this_mult.filter        = filter;
    this_mult.quad          = quad;
    this_mult.decision_vars = {q11, q12};

    % Define constraints
    constraints = [];
    if ~this_mult.constraint_q11_kyp
        c_q11 = (q11 >= 0):['SLTI Multiplier, ',...
                            this_mult.name,...
                            ', q11 >= 0'];                                      %#ok<BDSCA>
    else
        kyp_var_q11 = sdpvar(dim_state1);
        lmi_mat_q11 = kypLmiLti(br, q11, kyp_var_q11, this_mult.exponential);
        c_q11 = (lmi_mat_q11 >= 0):['SLTI Multiplier, ',...
                                    this_mult.name,...
                                    ', kyp(filter11, q11) >= 0'];               %#ok<BDSCA>
        this_mult.decision_vars{end + 1} = kyp_var_q11;
    end
    constraints = constraints + c_q11;
    if this_mult.constraint_q12_kyp
        kyp_var_q12 = sdpvar(dim_state1);
        lmi_mat_q12 = kypLmiLti(br, q12 + q12', kyp_var_q12);
        c_q12 = (lmi_mat_q12 == 0):['SLTI Multiplier, ',...
                                    this_mult.name,...
                                    ', kyp(filter11, q12 + q12^T) == 0'];       %#ok<BDSCA>
        this_mult.decision_vars{end + 1} = kyp_var_q12;
        constraints = constraints + c_q12;
    end
    this_mult.constraints = constraints;
end
end

% SET CONSTRAINTS TO BE MUTABLE (BUT DO SET.BLOCK_REALIZATION)
% SET Filter and constraints and quad to be private accessible
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)