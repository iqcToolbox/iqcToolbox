function [result, valid, yalmip_report, lft_analyzed] = iqcAnalysis(lft_in, varargin) 
%% IQCANALYSIS function for conducting IQC analysis on Ulft objects.
%  result.performance expresses the optimal objective value pertaining to the
%  performance multipliers combined with amplification due to uncertain initial conditions
%  If lft_in.performance is simply PerformanceL2Induced,
%     lft_in.disturbance is simply DisturbanceL2,
%     and initial conditions are 0 (options.init_cond_ellipse == nan)
%  then result.performance is an upper bound on the worst-case l2-induced norm 
%     (or hinf norm for LTI systems) of lft_in
%
%  [result, valid, yalmip_report, lft_analyzed] = iqcAnalysis(lft_in,...
%                                                             'multipliers_delta', mults_del,...
%                                                             'multipliers_performance', mults_perf,...
%                                                             'multipliers_disturbance', mults_dis,...
%                                                             'analysis_options', options)
%  result = iqcAnalysis(lft_in) assumes default expressions for 'multipliers_delta', 'multipliers_performance', 'multipliers_disturbance', 'analysis_options'
%
%  Variables:
%  ---------
%     Input:
%       lft_in                  : Ulft object 
%                               :: lft for IQC analysis
%       multipliers_delta       : (1 x num_del) array of MultiplierDelta objects 
%                               :: default or specified multipliers for each Delta in lft_in.delta 
%                                  (respecting the same ordering)
%       multipliers_performance : (1 x num_del) array of MultiplierPerformance objects 
%                               :: default or specified multipliers for each Performance in lft_in.performance
%                                  (respecting the same ordering)
%       multipliers_disturbance : (1 x num_del) array of MultiplierDisturbance objects 
%                               :: default or specified multipliers for each disturbance in lft_in.disturbance
%                                  (respecting the same ordering)
%       analysis_options        : AnalysisOptions object 
%                               :: options for IQC analysis and solving the associated SDP
%     Output:
%       result        : struct      :: structure holding information on the feasibility and solution of the associated SDP
%       valid         : bool        :: true if SDP is feasible and the solution satisfies all constraints
%       yalmip_report : struct      :: structure holding output information from YALMIP regarding SDP
%       lft_analyzed  : Ulft object :: the actual LFT analyzed.  This may differ from lft_in if it has Deltas 
%                                      that require transforming lft_in into an analyzable lft first.
%
%  See also Ulft, generateReachabilityLft

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% PARSE INPUTS
input_parser = inputParser;

addRequired(input_parser, 'lft_in', @(lft) isa(lft, 'Ulft'));

addParameter(input_parser,...
             'multipliers_delta',...
             initMultiplierDelta(length(lft_in.delta.deltas)),...
             @(mult) isa(mult, 'MultiplierDelta') && ~isempty(mult));

addParameter(input_parser,...
             'multipliers_disturbance',...
             initMultiplierDisturbance(length(lft_in.disturbance.disturbances)),...
             @(mult) isa(mult, 'MultiplierDisturbance') && ~isempty(mult));

addParameter(input_parser,...
             'multipliers_performance',...
             initMultiplierPerformance(length(lft_in.performance.performances)),...
             @(mult) isa(mult, 'MultiplierPerformance') && ~isempty(mult));

addParameter(input_parser,...
             'analysis_options',...
             AnalysisOptions(),...
             @(opts) isa(opts, 'AnalysisOptions'));
         
parse(input_parser, lft_in, varargin{:})

lft_in = input_parser.Results.lft_in;
mults_del  = input_parser.Results.multipliers_delta;
mults_dis  = input_parser.Results.multipliers_disturbance;
mults_perf = input_parser.Results.multipliers_performance;
options = input_parser.Results.analysis_options;

is_discrete = ~isempty(lft_in.timestep) && any(lft_in.timestep);
if isempty(options.exponential)
    if is_discrete
        options.exponential = 1;
    elseif ~isempty(lft_in.timestep) && all(~lft_in.timestep) 
    % is continuous-time
        options.exponential = 0;
    end
end



%% Check nominal stability
if lft_in.uncertain && ~isempty(lft_in.timestep)
    % We normalize here to ensure nominal analysis is in the appropriate range
    % Normalization is not forced on the intended analysis of the system because the user-provided
    % Multipliers pertain to the Deltas of the given lft, and not the normalized ones
    lft_normalized = lft_in.removeDisturbance(1:length(lft_in.disturbance.names))...
                           .removePerformance(1:length(lft_in.performance.names))...
                           .addPerformance({PerformanceStable()})...
                           .normalizeLft;
    nominal_system = removeUncertainty(lft_normalized,...
                                       lft_normalized.delta.names(2:end));
    op = AnalysisOptions('verbose', options.yalmip_settings.verbose,...
                         'solver', options.yalmip_settings.solver,...
                         'lmi_shift', options.lmi_shift,...
                         'exponential', options.exponential);
    [res, valid_nominal, yr, sys] = iqcAnalysis(nominal_system,...
                                                'analysis_options', op);
    if ~valid_nominal
        result = res;
        valid = valid_nominal;
        yalmip_report = yr;
        lft_analyzed = sys;
        disp(['Nominal system is unstable,',...
              ' therefore, uncertain system is not robustly stable'])
        return
    end
end

%% Modify LFT for pertinent Delta
lft_analyzed = modifyLft(lft_in);

%% Define joint Delta/Disturbance multiplier
exponential = options.exponential;
total_time = sum(lft_analyzed.horizon_period);
delta = lft_analyzed.delta;
num_dels = length(delta.names);

% Define multipliers for each Delta
num_mults_del = length(mults_del);
assert(num_dels == num_mults_del,...
       'iqcAnalysis:iqcAnalysis',...
       'Number of Deltas and Multipliers are not the same')

for i=1:num_dels
    type_del = delta.types{i}(6:end);
    type_mult = class(mults_del(i));
    if ~isa(mults_del(i), 'MultiplierDeltaDefault')
        assert(strcmp(type_mult(11:end), type_del),...
               'iqcAnalysis:iqcAnalysis',...
               'A given Multiplier''s class does not match with its Delta')
        assert(strcmp(mults_del(i).name, delta.names(i)),...
               'iqcAnalysis:iqcAnalysis',...
               'A given Multiplier"s name does not match with its Delta')
        assert(all(all(mults_del(i).horizon_period == ...
                       delta.horizon_periods(i, :))),...
               'iqcAnalysis:iqcAnalysis',...
               ['A given Multiplier"s horizon_period does not',...
                'match w/ its Delta'])
    else
        mults_del(i) = deltaToMultiplier(delta.deltas{i},...
                                         'discrete', is_discrete,...
                                         'exponential', exponential);
    end
end

% Define multipliers for each Disturbance
dis = lft_analyzed.disturbance;
num_disturbances = length(dis.names);
if num_disturbances == 0
    mults_dis = initMultiplierDisturbance(0);
end
assert(num_disturbances == length(mults_dis),...
       'iqcAnalysis:iqcAnalysis',...
       'Number of Disturbances and Multipliers are not the same')

for i=1:num_disturbances
    type_dis = dis.types{i}(12:end);
    type_mult = class(mults_dis(i));
    if ~isa(mults_dis(i), 'MultiplierDisturbanceDefault')
        assert(strcmp(type_mult(11:end), type_dis),...
               'iqcAnalysis:iqcAnalysis',...
               'A given Multiplier"s class doesnt match with its Disturbance')
        assert(strcmp(mults_dis(i).name, dis.names(i)),...
               'iqcAnalysis:iqcAnalysis',...
               'A given Multiplier"s name does not match with its Disturbance')
        assert(all(all(mults_dis(i).horizon_period == ...
                       dis.horizon_periods(i, :))),...
               'iqcAnalysis:iqcAnalysis',...
               ['A given Multiplier"s horizon_period does not',...
                'match w/ its Disturbance'])
    else
        mults_dis(i) = ...
            disturbanceToMultiplier(dis.disturbances{i},...
                                    'discrete', is_discrete,...
                                    'dim_in_lft', size(lft_analyzed, 2));
    end
end

% Define multipliers for each Performance
perf = lft_analyzed.performance;
num_performances = length(perf.names);
if num_performances == 0
    mults_perf = initMultiplierPerformance(0);
end
assert(num_performances == length(mults_perf),...
       'iqcAnalysis:iqcAnalysis',...
       'Number of Performances and Multipliers are not the same')

for i=1:num_performances
    type_perf = perf.types{i}(12:end);
    type_mult = class(mults_perf(i));
    if ~isa(mults_perf(i), 'MultiplierPerformanceDefault')
        assert(strcmp(type_mult(11:end), type_perf),...
               'iqcAnalysis:iqcAnalysis',...
               'A given Multiplier"s class doesnt match with its Performance')
        assert(strcmp(mults_perf(i).name, perf.names(i)),...
               'iqcAnalysis:iqcAnalysis',...
               'A given Multiplier"s name does not match with its Performance')
        assert(all(all(mults_perf(i).horizon_period == ...
                       perf.horizon_periods(i, :))),...
               'iqcAnalysis:iqcAnalysis',...
               ['A given Multiplier"s horizon_period does not',...
                'match w/ its Performance'])
    else
        mults_perf(i) = ...
            performanceToMultiplier(perf.performances{i},...
                                    'discrete', is_discrete,...
                                    'dim_out_lft', size(lft_analyzed, 1),...
                                    'dim_in_lft', size(lft_analyzed, 2));
    end
end

% Combine multipliers and augment with performance metric
mult_del  = MultiplierDeltaCombined(mults_del);
mult_dis  = MultiplierDisturbanceCombined(mults_dis);
mult_perf = MultiplierPerformanceCombined(mults_perf); 

dim_out = size(lft_analyzed, 1);
mult = combineAllMultipliers(mult_del, mult_dis, mult_perf, dim_out);

%% Check that exponential rate is amenable with lft_analyzed and mult_del (might be done in stability analysis)
if ~isempty(mult_del.exponential) && ~isempty(options.exponential)
    assert(mult_del.exponential == options.exponential,...
           'iqcAnalysis:iqcAnalysis',...
           ['The exponential factor of every multiplier must match the ',...
            'exponential factor specified in analysis_options'])
    
end
mult_filt = mult_del.filter_lft;
if lft_analyzed.uncertain && ~isempty(mult_filt.timestep)
% Uncertain system, and multiplier has dynamics
    if ~isempty(options.exponential)
    % Analysis specifies a non-default exponential rate, check if filter has 
    % exponential convergence rate
        mult_op = AnalysisOptions('verbose', options.yalmip_settings.verbose,...
                                  'solver', options.yalmip_settings.solver,...
                                  'lmi_shift', options.lmi_shift,...
                                  'exponential', options.exponential);
        mult_filt = mult_filt.addPerformance({PerformanceStable()});
        m_results = iqcAnalysis(mult_filt, 'analysis_options', mult_op);
        assert(m_results.valid,...
               'iqcAnalysis:iqcAnalysis',...
               ['The multiplier filters are not exponentially stable at the',...
                ' specified rate: options.exponential = ',...
                num2str(options.exponential),'\nEither specify multipliers',...
                ' whose filters have the given rate of convergence, \nor',...
                ' specify options.exponential to be a slower rate of',...
                ' convergence \n(for discrete time, closer to 1; for',...
                ' continuous time, closer to 0)',])
    end    
%     mult_filt = lftToSs(mult_filt);
%         p_expo = sdpvar(size(mult_filt.a, 1));
%         p_constraint = [];%p_expo >= options.lmi_shift * eye(size(p_expo));
%         mult_quad = eye(size(mult_filt, 1));
%         lmi_expo = kypLmiLti(mult_filt, mult_quad, p_expo, options.exponential);
%         lmi_expo = kypLmiLti(mult_filt, mult_quad, p_expo, 0);
%         lmi_constraint = lmi_expo <= -options.lmi_shift * eye(size(lmi_expo));
%         mult_constraint = p_constraint + lmi_constraint;
%         mult_objective = 0;
%         mult_settings = sdpsettings(options.yalmip_settings);
%         mult_report = optimize(lmi_constraint, mult_objective, mult_settings);
%         mult_residual = check(mult_constraint);
%         valid_solution = all(mult_residual >= -options.lmi_shift);
%         assert(valid_solution,...
%                'iqcAnalysis:iqcAnalysis',...
%                ['The multiplier filters are not exponentially stable to the',...
%                 ' specified rate: options.exponential = ',...
%                 num2str(options.exponential),'\nEither specify multipliers',...
%                 ' whose filters have the given rate of convergence, \nor',...
%                 ' specify options.exponential to be a slower rate of',...
%                 ' convergence \n(for discrete time, closer to 1; for',...
%                 ' continuous time, closer to 0)',])
%     end
end

%% Formulate KYP lmis
[objective, kyp_constraints, state_amplification, ellipse, kyp_variables, exponential] = ...
                                    kypLmi(lft_analyzed, mult, options);

%% Gather constraints, set yalmip options, solve
constraints = mult.constraints + kyp_constraints;
settings = sdpsettings(options.yalmip_settings);
yalmip_report = optimize(constraints, objective, settings);

%% Check correctness, gather data
% Ignore ghost variable constraint if it exists and is "close enough" to equality
try 
    ghost_close = check(constraints('Ghost variable == 1')) > -2e-7;
    if ghost_close
        % Eliminate inaccuracy of near equality satisfaction from constraints
        constraints('Ghost variable == 1') = [];
    end
catch                                                                           %#ok<CTCH>
end

primal_residual = check(constraints);

% Check if the system_kyp constraints are truly negative definite
% (the primal_residuals may be negative because of the shift)
system_kyp_residuals = nan(total_time, 1);
for i = 1:total_time
    system_kyp_residuals(i) = check(constraints(['KYP LMI, ', num2str(i)]));
end
% If KYP matrix PD constraints exist, check if they're valid too
for i = 1:total_time
    try
        ct_tag = ['KYP matrix PD, ' num2str(i)];
        system_kyp_residuals(total_time + i) = check(constraints(ct_tag));
    catch                                                                       %#ok<CTCH>
    end
end
% If Initial Condition LMI constraint exists, check validity too
try
    system_kyp_residuals(end+1) = check(constraints('Initial Condition LMI'));
catch                                                                           %#ok<CTCH>
end
% How many system_kyp constraints fail YALMIP residuals?
system_kyp_failed_yalmip = sum(system_kyp_residuals < 0);
% How many YALMIP residuals are less than LMI shift?
system_kyp_neg_def =  sum((system_kyp_residuals(system_kyp_residuals < 0) ...
                           + options.lmi_shift) >= 0);
% How many YALMIP constraints failed?
failed_yalmip = sum(primal_residual < 0);
% Are all failing YALMIP-constraints the system_kyp-constraints that are
% truly negative-definite?
primal_satisfied = (failed_yalmip == system_kyp_failed_yalmip)...
                   && (system_kyp_failed_yalmip == system_kyp_neg_def);

valid = all(all(~isnan(primal_residual))) && ...
        (all(primal_residual >= 0) || primal_satisfied);        
if valid
    performance = sqrt(cellfun(@double, {mults_perf.objective}));
else
    performance = inf;
end

result.performance             = performance;
result.state_amplification     = double(state_amplification);
result.ellipse                 = double(ellipse);
result.multiplier_combined     = mult;
result.multipliers_delta       = mults_del;
result.multipliers_disturbance = mults_dis;
result.multipliers_performance = mults_perf;
result.kyp_variables           = kyp_variables;
result.exponential             = exponential;
result.debug.constraints       = constraints;
result.debug.yalmip_report     = yalmip_report;
result.valid                   = valid;
end

function [objective,...
          kyp_constraints,...
          state_amplification,...
          ellipse,...
          kyp_variables,...
          expo] = kypLmi(lft_in, mult, options)
%% KYPLMI function for generating variables and constraints associated with
%  the Kalman-yakubovich-popov lemma. 
%  Variables:
%     lft_in : Ulft object
%     mult : MultiplierPerformanceCombined object

total_time = sum(lft_in.horizon_period);

% Formulate state-space matrices for LMI
if isempty(lft_in.timestep)
    state_in = zeros(1, total_time);
    state_out = state_in;
else
    state_in = lft_in.delta.dim_outs(1,:);
    state_out = lft_in.delta.dim_ins(1,:);
end
a_lft = cell(1, total_time);
b_lft = cell(1, total_time);
c_lft = cell(1, total_time);
d_lft = cell(1, total_time);
for i = 1:total_time
    % Create system matrices for lft with only the state operator in the Delta block
    abcd = [lft_in.a{i}, lft_in.b{i};
            lft_in.c{i}, lft_in.d{i}];
    a_lft{i} = abcd(1:state_out(i), 1:state_in(i));
    b_lft{i} = abcd(1:state_out(i), state_in(i)+1:end);
    c_lft{i} = abcd(state_out(i)+1:end, 1:state_in(i));
    d_lft{i} = abcd(state_out(i)+1:end, state_in(i)+1:end);
end
if ~isempty(lft_in.timestep)
    delta_state = SequenceDelta(lft_in.delta.deltas{1});
    if lft_in.timestep(1) % Discrete-time system
        delta_state.deltas{1}.timestep = -1 * ones(1, total_time);
    end 
else
    delta_state = SequenceDelta();
end
lft_in_state = Ulft(a_lft, b_lft, c_lft, d_lft, delta_state,...
                    'horizon_period', lft_in.horizon_period);
dim_in = num2cell(size(lft_in_state, 2));

% Generate system: filter_lft_eye = mult.filter * [lft ; eye]
identity_cell = cellfun(@(dim) eye(dim), dim_in, 'UniformOutput', false);
if all(cellfun(@isempty, identity_cell))
    % This will corner case will occur if the combined multiplier is empty
    % (i.e., stability analysis on nominal systems)
    filt_lft_eye = lft_in_state;
else
    identity_lft = toLft(identity_cell, lft_in.horizon_period);
    filt_lft_eye = mult.filter_lft * [lft_in_state; identity_lft];
end

% Set up exponential rate and check for consistency
expo = options.exponential;
if ~isempty(filt_lft_eye.timestep)
% LFT has dynamics
    if isempty(options.exponential)
    % Exponential rate is unspecified
        if filt_lft_eye.timestep
        % LFT is discrete-time
            expo = 1;
        else
        % LFT is continuous-time
            expo = 0;
        end
    end
    default_exponential =    (all(logical(filt_lft_eye.timestep)) && expo == 1)...
                          || (all(~logical(filt_lft_eye.timestep)) && expo == 0);
    if ~default_exponential    
    % Non-default exponential rate is specified
    % for either discrete- or continuous-time system
        assert(isequal(filt_lft_eye.horizon_period, [0, 1]),...
               'iqcAnalysis:kypLmi',...
               'Cannot have non-default exponential rate for time-varying LFTs')
    end
end

% Reset state_in dimensions to pertain to filt_lft_eye, rather than lft_in
if ~isempty(filt_lft_eye.timestep)
    state_in = filt_lft_eye.delta.dim_outs(1,:);
else
    state_in = zeros(1, total_time);
end
% Reset abcd matrices to pertain to filt_lft_eye
abcd = cellfun(@(a, b, c, d) [a, b; c, d],...
              filt_lft_eye.a, filt_lft_eye.b, filt_lft_eye.c, filt_lft_eye.d,...
              'UniformOutput', false);

% Initialize lyapunov matrices
p = cell(1, total_time);
for i = 1:total_time
    p{i} = sdpvar(state_in(i));
end
if ~isempty(options.p0)
    p{1} = options.p0;
end
lmi_shift         = options.lmi_shift;
kyp_constraints   = [];

% Force P to be positive definite for nominal systems or non-default exponential rates
if ~lft_in.uncertain || ~default_exponential %&& ~isempty(lft_in.timestep)
    for i = 1:total_time
    mat_shift = lmi_shift * eye(size(p{i}, 2));

    if ~isempty(p{i})
    kyp_constraints = kyp_constraints ...
                      + ((p{i} >= mat_shift):['KYP matrix PD, ' num2str(i)]);   %#ok<BDSCA>
    end
    end
end

lmi_mat = cell(1, total_time);
if isempty(filt_lft_eye.timestep)
% LMIs for memoryless uncertain system
    for i = 1:total_time
        quad_now = [mult.quad.q11{i}, mult.quad.q12{i}; 
                    mult.quad.q21{i}, mult.quad.q22{i}];
        % Note that matrix dimensions conform because a, b, and c are empty
        lmi_mat{i} = abcd{i}' * quad_now * abcd{i};
    end
elseif any(filt_lft_eye.timestep)
% Discrete-time KYP LMIs
    % Initialize eventually periodic time indices
    time_indices = 1 : total_time + 1;
    time_indices(end) = time_indices(filt_lft_eye.horizon_period(1) + 1);
        
    for i = 1:total_time
        p_now = p{time_indices(i)};
        quad_now = [mult.quad.q11{i}, mult.quad.q12{i}; 
                    mult.quad.q21{i}, mult.quad.q22{i}];
        p_next = p{time_indices(i + 1)};
        lmi_mat{i} = [eye(state_in(i), size(abcd{i}, 2)); abcd{i}]' * ...
                     blkdiag(-expo^2 * p_now, p_next, quad_now) * ...
                     [eye(state_in(i), size(abcd{i}, 2)); abcd{i}];
    end
elseif ~(all(filt_lft_eye.timestep))
% Continuous-time KYP LMIs
    assert(isequal(filt_lft_eye.horizon_period, [0, 1]),...
           'iqcAnalysis:kypLmi',...
           'Can only conduct IQC analysis on time-invariant continuous-time systems')
    p_mat    = [2 * expo * p{1},  p{1};
                p{1},             zeros(size(p{1}))];
    quad_mat = [mult.quad.q11{1}, mult.quad.q12{1};
                mult.quad.q21{1}, mult.quad.q22{1}];
    lmi_mat{i} = [eye(state_in(1), size(abcd{1}, 2)); abcd{1}]' * ...
                 blkdiag(p_mat, quad_mat) * ...
                 [eye(state_in(1), size(abcd{1}, 2)); abcd{1}];
end
for i = 1:total_time
% Set negative-definiteness constraints on lmi_mat
    % Correct if floating operations make asymmetric
    if ~issymmetric(lmi_mat{i})
        lmi_mat{i} = (lmi_mat{i} + lmi_mat{i}')/2;
    end
    if isa(lmi_mat{i}, 'sdpvar')
    mat_shift = lmi_shift * eye(size(lmi_mat{i}, 2));
    kyp_constraints = kyp_constraints ...
                      + ((lmi_mat{i} <= -mat_shift):['KYP LMI, ' num2str(i)]);  %#ok<BDSCA>
    else
    % Make a ghost variable and a constraint that lmi eigs are negative.  This is done
    % in order to parse problems w/o decision vars as yalmip optimization problems
        if i == 1
            ghost = sdpvar(1);
            kyp_constraints = kyp_constraints + ...
                              (ghost == 1):'Ghost variable == 1';               %#ok<BDSCA>
        end
    lmi_eig = eig(lmi_mat{i});
    if all(all(lmi_eig == 0))
        lmi_eig = lmi_eig + 1e-8;
        lmi_shift = lmi_shift - 1e-8;
    end
    kyp_constraints = kyp_constraints ...
        + ((ghost * lmi_eig <= 0*lmi_eig - lmi_shift):['KYP LMI, ' num2str(i)]);%#ok<BDSCA>
    end
end
% Set constraints for initial conditions
if any(any(isnan(options.init_cond_ellipse)))
    % Initial condition at zero
    ellipse = nan;
    state_amplification = nan;
    objective_state = 0;
    p0_state = p{1};
elseif ~isempty(filt_lft_eye.timestep) && any(filt_lft_eye.timestep)
    % Handle uncertain initial conditions for discrete-time systems
    ellipse = options.init_cond_ellipse;
    states = options.init_cond_states;
    states = [false(1, size(mult.filter.a{1}, 1)), states];
    p0_state = p{1}(states,states);
    if any(any(isinf(ellipse)))
        ellipse = sdpvar(size(ellipse, 1)); 
        state_amplification = 1;
        ellipse_eigenvalues = sdpvar(size(ellipse,1),1);
        objective_state = sum(ellipse_eigenvalues'*ellipse_eigenvalues);
        ellipse = diag(ellipse_eigenvalues);
        kyp_constraints = kyp_constraints + (ellipse_eigenvalues >= 0);
        kyp_constraints = kyp_constraints + (p0_state <= ellipse);
    else
        objective_state = sdpvar(1);
        state_amplification = sqrt(objective_state);
        mat_shift = lmi_shift * eye(size(ellipse, 2));
        ic = (p0_state <= objective_state * ellipse - mat_shift):...
                                            'Initial Condition LMI';            %#ok<BDSCA>
        kyp_constraints = kyp_constraints + ic; 
    end
else
    error('iqcAnalysis:kypLmi',...
          ['Uncertain initial conditions can only be incorporated for',...
          ' discrete-time system'])
end
objective = mult.objective + options.scale_state_obj * objective_state;  
kyp_variables = horzcat(p, {p0_state});
end

function lft_out = modifyLft(lft_in)
%% MODIFYLFT function for generating modified LFTs which are analyzed in 
%  place of the original LFT.  This is only necessitated for specific
%  Deltas (as encoded in their class), such as DeltaRateBndSltv, and relies
%  on an extension of the Delta/Disturbance/Performance.modifyLft methods
%    Variables:
%    ---------
%      Input:
%         lft_in : Ulft object
%      Output:
%         lft_out : Ulft object

% For deltas
lft_out = lft_in;
deltas = lft_out.delta.deltas;
for i = 1:length(deltas)
    mod_lft_handle = deltas{i}.modifyLft();
    if ~isempty(mod_lft_handle)
        lft_out = mod_lft_handle(lft_out);
    end
end

% For Performances
perfs = lft_out.performance.performances;
for i = 1:length(perfs)
    mod_lft_handle = perfs{i}.modifyLft();
    if ~isempty(mod_lft_handle)
        lft_out = mod_lft_handle(lft_out);
    end
end
end    

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)