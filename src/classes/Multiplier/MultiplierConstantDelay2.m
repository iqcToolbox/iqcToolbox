classdef MultiplierConstantDelay2 < MultiplierDelta
%% MULTIPLIERConstantDelay2 class for uncertainties of constant delay with a 
%  given maximum allowable delay (DeltaConstantDelay2).
%  Extends the base class MultiplierDelta.
%
%  extended methods:
%     MultiplierConstantDelay2(delta, varargin) :: Constructor method
%     
%  extended properties:
%     discrete : logical :: true if discrete-time, false if continuous-time
%     basis_length : natural :: length of filter's basis_function
%     basis_poles : array of complex doubles :: list of basis_function poles
%     basis_function : tf :: a (basis_length+1 x 1) stable transfer
%                            function for defining filter
%     basis_realization : ss :: the ss realization of basis_function
%     basis_delay : ss :: the ss realization of a basis function whose magnitude
%                         is greater than any permissible delay
%     dim_outin : double :: dimensions of uncertainty
%     constraint_q11_kyp : logical :: true if a kyp-based constraint for
%                                     quad.q11 is desired
%     constraint_q22_kyp : logical :: true if a kyp-based constraint for
%                                     quad.q22 is desired
%     delay_max : double :: maximum allowed delay, informs constraints for
%                           generalized, exponential IQCs
%
%  See also MultiplierConstantDelay2.MultiplierConstantDelay2

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    basis_length double
    basis_poles double
    basis_function tf
    basis_realization ss
    basis_delay ss
    delay_max double
    dim_outin double
    delay_filter_eps double
end

properties (SetAccess = immutable)
    constraint_q11_kyp logical
    constraint_q2_minus_q1_kyp logical
    constraint_q2_plus_q1_kyp logical
end

methods
function this_mult = MultiplierConstantDelay2(delta, varargin)
%% MULTIPLIERConstantDelay2 constructor
%
%  this_mult = MultiplierConstantDelay2(delta, 'discrete', true, 'basis_length', 2, 'basis_pole', -0.5, 'constraint_q11_kyp', true, 'constraint_q22_kyp', true, 'exponential', [], delay_filter_eps, 0.2)
%  this_mult = MutliplierConstantDelay2(delta) assumes the values provided above
%
%  Variables:
%  ---------
%    Input:
%       delta : The delta object characterized by this multiplier
%       varargin : Key-value pair of optional arguments for multiplier (see below for more details)
%    Output:
%       this_mult : MultiplierConstantDelay2 object
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
%       constraint_q2_minus_q1_kyp : logical :: Whether or not a constraint
%                            on quad.q22 will involve the kyp lemma
%                            (the former will possibly be less
%                            conservative, but more computationally
%                            challenging)
%       constraint_q2_plus_q1_kyp : logical :: Whether or not a constraint
%                            on quad.q22 will involve the kyp lemma
%                            (the former will possibly be less
%                            conservative, but more computationally
%                            challenging)
%    filter_construction : This group contains the following variables:
%       basis_length : natural :: length of filter's basis_function
%       basis_poles : array of complex doubles :: list of basis_function poles
%       basis_function : tf :: a (basis_length+1 x 1) stable transfer
%                              function for defining filter
%       basis_realization : ss :: the ss realization of basis_function  
%       exponential : double :: exponential rate bound used to define
%                               exponential IQCs
%       delay_filter_eps : double :: small, positive scalar used to define basis_delay,
%                                    the portion of the filter related to delay_max
%       the constructed filter appears as filter = [basis_realization * basis_delay, 0
%                                                              basis_realization,    0
%                                                                      0        , basis_realization]
%                                    
%
%  A MultiplierConstantDelay2 can be constructed by providing inputs with varying levels of granularity:
%     least granular --------------------------------------------------- most granular
%        [basis_length;     -->     basis_function      -->     basis_realization    
%         basis_pole]
%  If users provide less granular inputs, those inputs will construct the more granular ones.
%  If users provide more granular inputs, the less granular properties will be left unspecified.
%
%  See also MultiplierConstantDelay2

% Inputs
input_parser = inputParser;
addRequired(input_parser,...
            'delta',...
            @(del)validateattributes(del, {'DeltaConstantDelay2'}, {'nonempty'}))
addParameter(input_parser,...
             'discrete',...
             true,...
             @(disc) validateattributes(disc, {'logical'}, {'nonempty'}))
addParameter(input_parser,...
             'constraint_q11_kyp',...
             true,...
             @(constraint) islogical(constraint))
addParameter(input_parser,...
             'constraint_q2_minus_q1_kyp',...
             true,...
             @(constraint) islogical(constraint))
addParameter(input_parser,...
             'constraint_q2_plus_q1_kyp',...
             true,...
             @(constraint) islogical(constraint))
addParameter(input_parser,...
             'basis_length',...
             2,...
             @(bl) validateattributes(bl, {'numeric'}, {'integer', 'positive'}))
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
             'exponential',...
             [],...
             @(expo) validateattributes(expo, {'numeric'}, {'nonnegative'}))
addParameter(input_parser,...
             'delay_filter_eps',...
             0.2,...
             @(expo) validateattributes(expo, {'numeric'}, {'positive'}))

% Parsing inputs
parse(input_parser, delta, varargin{:})
delta                       = input_parser.Results.delta;
discrete                    = input_parser.Results.discrete;
constraint_q11_kyp          = input_parser.Results.constraint_q11_kyp;
constraint_q2_minus_q1_kyp  = input_parser.Results.constraint_q2_minus_q1_kyp;
constraint_q2_plus_q1_kyp   = input_parser.Results.constraint_q2_plus_q1_kyp;
basis_length                = input_parser.Results.basis_length;
basis_poles                 = input_parser.Results.basis_poles;
basis_function              = input_parser.Results.basis_function;
basis_realization           = input_parser.Results.basis_realization;
exponential                 = input_parser.Results.exponential;
delay_filter_eps            = input_parser.Results.delay_filter_eps;
if isempty(exponential)
% Set default value of exponential depending on continuous- or discrete-time
    if discrete
        exponential = 1;
    else
        exponential = 0;
    end
end
if discrete
    assert(0 < exponential && exponential <= 1,...
           'MultiplierConstantDelay2:MultiplierConstantDelay2',...
           'Exponential IQCs in discrete-time must have exponential in (0, 1]')
    assert(floor(delta.delay_max) == delta.delay_max,...
           'MultiplierConstantDelay2:MultiplierConstantDelay2',...
           ['Delta ', delta.name, ' must have delay_max be a natural number',...
            ' if the uncertain system is discrete-time'])
end

% Setting multiplier properties
this_mult@MultiplierDelta(exponential);
this_mult.name                          = delta.name;
this_mult.horizon_period                = delta.horizon_period;
this_mult.dim_outin                     = delta.dim_out;
this_mult.delay_max                     = delta.delay_max;
this_mult.discrete                      = discrete;
this_mult.constraint_q11_kyp            = constraint_q11_kyp;
this_mult.constraint_q2_minus_q1_kyp    = constraint_q2_minus_q1_kyp;
this_mult.constraint_q2_plus_q1_kyp     = constraint_q2_plus_q1_kyp;
this_mult.delay_filter_eps              = delay_filter_eps;

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
        assert(all(all(abs(basis_poles) < 1)),...
               'MultiplierConstantDelay2:MultiplierConstantDelay2',...
               'basis_poles must be within the unit circle')
    else
        assert(all(all(real(basis_poles) < 0)),...
               'MultiplierConstantDelay2:MultiplierConstantDelay2',...
               'basis_poles must be in the left plane')
    end
    if size(basis_poles, 1) > 1
        assert(size(basis_poles, 1) == basis_length - 1,...
               'MultiplierConstantDelay2:MultiplierConstantDelay2',...
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
                       'MultiplierConstantDelay2:MultiplierConstantDelay2',...
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
               'MultiplierConstantDelay2:set:basis_function',...
               'time domain of basis_function must match multiplier')
        assert(isstable(basis_function),...
               'MultiplierConstantDelay2:set:basis_function',...
               'basis_function must be stable')
        assert(size(basis_function, 2) == 1,...
               'MultiplierConstantDelay2:set:basis_function',...
               'basis_function must be of size ? x 1')                
        this_mult.basis_realization = ss(basis_function, 'minimal');
    end
end

function this_mult = set.basis_realization(this_mult, basis_realization)
    % Set basis_realization, check consistency
    this_mult.basis_realization = basis_realization;
    assert(ifAndOnlyIf(this_mult.discrete, (basis_realization.Ts~=0)),...
           'MultiplierConstantDelay2:set:basis_realization',...
           'time domain of basis_realization must match multiplier')
    assert(isstable(basis_realization),...
           'MultiplierConstantDelay2:set:basis_realization',...
           'basis_realization must be stable')
    assert(size(basis_realization, 2) == 1,...
           'MultiplierConstantDelay2:set:basis_realization',...
           'basis_realization must be of size ? x 1')
    basis    = basis_realization;
    dim_outin = this_mult.dim_outin(1);
    
    % Set basis_delay
    dmax = this_mult.delay_max;
    expo = this_mult.exponential;
    if this_mult.discrete
        omega = logspace(-4, log10(pi), 100);
        upper_bound = zeros(1, length(omega));
        z = tf('z');
        for delay = 1:dmax
            sys = 1/ expo ^ delay / z^delay - 1;
            sys_fr = frd(sys, omega);
            upper_bound = max([upper_bound;abs(squeeze(sys_fr.ResponseData))']);
        end
        upper_bound = upper_bound + sqrt(expo ^(-2 * dmax) - 1);
        ubnd_fr = frd(upper_bound, omega, sys.Ts);
        n_states = 4;
%         basis_delay = fitfrd(ubnd_fr, n_states);
        mag_constraint.UpperBound = [];
        mag_constraint.LowerBound = ubnd_fr;
        basis_delay = fitmagfrd(ubnd_fr, n_states, [], [], mag_constraint);
        shifted_delay = basis_delay;
        shifted_delay.a = expo * shifted_delay.a;
        shifted_delay.b = expo * shifted_delay.b;
        this_mult.basis_delay = shifted_delay;
    else
        s = tf('s');
        num = (s + 4/dmax/pi) * (s + this_mult.delay_filter_eps/dmax + expo);
        den = s^2 - pi * cos(pi^2 / 4) / dmax * s + pi^2 / dmax^2 / 4;
        gain = (1 + exp(expo * dmax));
        basis_delay = ss(gain * num / den);
        shifted_delay = ss(basis_delay + sqrt(exp(2 * expo * dmax) - 1));
        shifted_delay.a = shifted_delay.a - expo * eye(size(shifted_delay.a));
        this_mult.basis_delay = shifted_delay;
    end
    
    % Set filter blocks
    blk1 = [this_mult.basis_delay * this_mult.basis_realization;
             this_mult.basis_realization];
    blk11 = ss(kron(blk1.a, eye(dim_outin)), kron(blk1.b, eye(dim_outin)),...
               kron(blk1.c, eye(dim_outin)), kron(blk1.d, eye(dim_outin)),...
               basis.Ts);
    blk2 = this_mult.basis_realization;
    blk22 = ss(kron(blk2.a, eye(dim_outin)), kron(blk2.b, eye(dim_outin)),...
               kron(blk2.c, eye(dim_outin)), kron(blk2.d, eye(dim_outin)),...
               basis.Ts);
    
    
    % Define exponential stability variables
    if this_mult.discrete
        expo = 1 / this_mult.exponential ^ (2 * this_mult.delay_max);
    else
        expo = exp(2 * this_mult.exponential * this_mult.delay_max);
    end

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

    [dim_out1, dim_in1] = size(blk11.d);
    [dim_out2, dim_in2] = size(blk22.d);
    dim_state1 = size(blk11.a, 1);
    dim_state2 = size(blk22.a, 1);
    q1 = sdpvar(dim_out2);
    q2 = sdpvar(dim_out2);
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
%         filter.c1{i}  = expo_bound * [br.c, zeros(dim_out1, dim_state2)];
%         filter.d11{i} = expo_bound * br.d;

        quad.q11{i} = blkdiag(q1, zeros(dim_out2));
        quad.q12{i} = [zeros(dim_out2); q2];
        quad.q21{i} = quad.q12{i}';
        quad.q22{i} = q2 - q1;
    end
    this_mult.filter        = filter;
    this_mult.quad          = quad;
    this_mult.decision_vars = {q1, q2};

    % Define constraints
    constraints = [];
    if ~this_mult.constraint_q11_kyp
        c_q11 = (q1 >= 0):['Constant Delay Multiplier, ',...
                                   this_mult.name,...
                                   ', q1 >= 0'];                               %#ok<BDSCA>
    else
        kyp_var_q11 = sdpvar(dim_state2);
        lmi_mat_q11 = kypLmiLti(blk22, q1, kyp_var_q11, this_mult.exponential);
        c_q11 = (lmi_mat_q11 >= 0):['Constant Delay Multiplier, ',...
                                    this_mult.name,...
                                    ', kyp(filter22, q1) >= 0'];               %#ok<BDSCA>
        this_mult.decision_vars{end + 1} = kyp_var_q11;
    end
    if ~this_mult.constraint_q2_minus_q1_kyp
        c_q2_m_q1 = (q2 - q1 <= 0):['Constant Delay Multiplier, ',...
                                   this_mult.name,...
                                   ', q2 - q1 <= 0'];                               %#ok<BDSCA>
    else
        kyp_var_q2_m_q1 = sdpvar(dim_state2);
        lmi_mat_q2_m_q1 = kypLmiLti(blk22,...
                                    q2 - q1,...
                                    kyp_var_q2_m_q1,...
                                    this_mult.exponential);
        c_q2_m_q1 = (lmi_mat_q2_m_q1 <= 0):['Constant Delay Multiplier, ',...
                                      this_mult.name,...
                                      ', kyp(filter22, q2 - q1) <= 0'];               %#ok<BDSCA>
        this_mult.decision_vars{end + 1} = kyp_var_q2_m_q1;
    end
    if expo ~= 1
        if ~this_mult.constraint_q2_plus_q1_kyp
            c_q2_p_q1 = (q2 + q1 >= 0):['Constant Delay Multiplier, ',...
                                       this_mult.name,...
                                       ', q2 + q1 >= 0'];                               %#ok<BDSCA>
        else
            kyp_var_q2_p_q1 = sdpvar(dim_state2);
            lmi_mat_q2_p_q1 = kypLmiLti(blk22,...
                                        q2 + q1,...
                                        kyp_var_q2_p_q1,...
                                        this_mult.exponential);
            c_q2_p_q1 = (lmi_mat_q2_p_q1 >= 0):['Constant Delay Multiplier, ',...
                                          this_mult.name,...
                                          ', kyp(filter22, q2 + q1) >= 0'];               %#ok<BDSCA>
            this_mult.decision_vars{end + 1} = kyp_var_q2_p_q1;
        end
    else
        c_q2_p_q1 = [];
    end
    
    constraints = constraints + c_q11 + c_q2_m_q1 + c_q2_p_q1;
    this_mult.constraints = constraints;
end
end
end

%%  CHANGELOG
% Apr. 8, 2021: Added after v0.9.0 - Micah Fry (micah.fry@ll.mit.edu)