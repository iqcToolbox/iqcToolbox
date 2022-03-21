%% Requirements:
%  1. IQC analysis shall measure the degree to which uncertain initial 
%     conditions excite the upper-bound on an uncertain system's output

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for IQC analysis with uncertain initial conditions
classdef testIqcAnalysisUncertainInitialCondition < matlab.unittest.TestCase
methods (Test)

function testReachabilityPureDelay(testCase)
    gain = 2;
    dim_state = 3;
    a_delay = {gain * eye(dim_state)};
    b_delay = {zeros(dim_state, 1)};
    c_delay = {eye(dim_state)};
    d_delay = {zeros(dim_state, 1)};
    timestep = -1;
    horizon_period = [0, 1];
    lft_delay = Ulft(a_delay, b_delay, c_delay, d_delay,...
                     DeltaDelayZ(dim_state, timestep, horizon_period),...
                     'horizon_period', horizon_period);
    options = AnalysisOptions('verbose', false,...
                              'lmi_shift', 1e-6,...
                              'init_cond_ellipse', eye(dim_state),...
                              'init_cond_states', true(1, dim_state));
                          
    final_time = 2;
    reach_lft2 = generateReachabilityLft(lft_delay, final_time);                          
    [result2, valid] = iqcAnalysis(reach_lft2, 'analysis_options', options);
    assertTrue(testCase, valid)
    verifyLessThan(testCase,...
                   abs(result2.state_amplification - gain^final_time),...
                   1e-2)    
    
    final_time = 10;
    reach_lft10 = generateReachabilityLft(lft_delay, final_time);
    [result10, valid] = iqcAnalysis(reach_lft10, 'analysis_options', options);
    assertTrue(testCase, valid)
    verifyLessThan(testCase,...
                   abs(result10.state_amplification - gain^final_time),...
                   1e-2)      
end

function testReachabilityWithUncertainInitialCondition(testCase)
    % Create a pure state-delay system (inputs make no impact)
    dim_state = 3;
    dim_out = dim_state;
    dim_in = 1;
    a = eye(dim_state);
    b = zeros(dim_state, dim_in);
    c = eye(dim_out);
    d = zeros(dim_out, dim_in);
    lft = Ulft(a, b, c, d, DeltaDelayZ(3));
    % Uncertain initial conditions within unit ball
    ellipse = eye(dim_state);
    uncertain_states = true(dim_state, 1);
    options = AnalysisOptions('verbose', false,...
                              'lmi_shift', 1e-6,...
                              'init_cond_ellipse', ellipse,...
                              'init_cond_states', uncertain_states);
    % Exagerating disturbance input to try and inflate state_amplification objective
    % (regardless of the inflation, it should still be minimally larger than 1)
    bound_l2_norm_in = 1e4;
    % Make reachability LFT and analyze
    N = 10;
    for i = 1:N
        lft_reach = generateReachabilityLft(lft, i);
        [dim_out, dim_in] = size(lft_reach);
        mult_perf = MultiplierL2Induced(lft_reach.performance.performances{1},...
                                        dim_out, dim_in,...
                                        'objective_scaling', bound_l2_norm_in^2);
        result = iqcAnalysis(lft_reach,...
                             'analysis_options', options,...
                             'multipliers_performance', mult_perf);
        final_state_bound = result.state_amplification;
        % final_state_bound should be 1. Checking against inflated bound in case
        % platforms have poor solver performance
        verifyLessThan(testCase, final_state_bound, 1.01)        
    end   
    
    % Understating disturbance input to drive state_amplification object as low a possible
    % (regardless, state_amplification should never be less than 1)
    bound_l2_norm_in = 1e-4;
    % Make reachability LFT and analyze
    N = 10;
    for i = 1:N
        lft_reach = generateReachabilityLft(lft, i);
        [dim_out, dim_in] = size(lft_reach);
        mult_perf = MultiplierL2Induced(lft_reach.performance.performances{1},...
                                        dim_out, dim_in,...
                                        'objective_scaling', bound_l2_norm_in^2);
        result = iqcAnalysis(lft_reach,...
                             'analysis_options', options,...
                             'multipliers_performance', mult_perf);
        final_state_bound = result.state_amplification;
        % final_state_bound should be 1. Checking against inflated bound in case
        % platforms have poor solver performance
        verifyGreaterThan(testCase, final_state_bound, 1)        
    end        
    
    % Now try to make the initial ellipse an decision variable, should be close to eye(3)
    options = AnalysisOptions('verbose', true,...
                              'lmi_shift', 1e-6,...
                              'init_cond_ellipse', inf(3),...
                              'init_cond_states', uncertain_states);
    result = iqcAnalysis(lft_reach,...
                             'analysis_options', options,...
                             'multipliers_performance', mult_perf);
    init_cond_ellipse = value(result.kyp_variables{1});
    diff_ellipse = norm(init_cond_ellipse - eye(3), 2);
    verifyLessThan(testCase, diff_ellipse, 1e-2)       
end

function testErrorUncertainInitialConditions(testCase)
    lft = toLft(rss);
    dim_state = lft.delta.deltas{1}.dim_out;
    ellipse = eye(dim_state);
    uncertain_states = true(dim_state, 1);
    options = AnalysisOptions('verbose', false,...
                              'init_cond_ellipse', ellipse,...
                              'init_cond_states', uncertain_states);
    verifyError(testCase,...
                @() iqcAnalysis(lft, 'analysis_options', options),...
                'iqcAnalysis:kypLmi')
end

end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)