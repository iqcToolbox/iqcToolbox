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
                   abs(sqrt(result2.state_amplification) - gain^final_time),...
                   1e-2)    
    
    final_time = 10;
    reach_lft10 = generateReachabilityLft(lft_delay, final_time);
    [result10, valid] = iqcAnalysis(reach_lft10, 'analysis_options', options);
    assertTrue(testCase, valid)
    verifyLessThan(testCase,...
                   abs(sqrt(result10.state_amplification) - gain^final_time),...
                   1e-2)      
end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)