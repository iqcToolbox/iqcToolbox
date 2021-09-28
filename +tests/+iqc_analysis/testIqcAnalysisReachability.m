%% Requirements:
%  1. IQC analysis shall produce a worst-case reachability upper-bound 
%     for uncertain-systems that are robustly stable. 

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for IQC-based reachability analysis
classdef testIqcAnalysisReachability < matlab.unittest.TestCase
methods (Test)

function testReachabilityResults(testCase)
    zero = [];
    pole = -.5;
    gain = 1;
    timestep = -1;
    g = ss(zpk(zero, pole, gain, timestep));
    lft = toLft(g);
    
    options = AnalysisOptions('verbose', false);
    [result0, valid] = iqcAnalysis(lft, 'analysis_options', options);
    assertTrue(testCase, valid)
    
    final_time = 2;
    reach_lft2 = generateReachabilityLft(lft, final_time);
    [result2, valid] = iqcAnalysis(reach_lft2, 'analysis_options', options);
    assertTrue(testCase, valid)
    verifyLessThan(testCase, result2.performance, result0.performance)
    
    final_time = 100;
    reach_lft100 = generateReachabilityLft(lft, final_time);
    [result100, valid] = iqcAnalysis(reach_lft100, 'analysis_options',options);
    assertTrue(testCase, valid)
    verifyLessThan(testCase, result100.performance, result0.performance)
    verifyGreaterThan(testCase, result100.performance, result2.performance)
end

function testReachabilityDelayDynamics(testCase)
    gain = 1;
    a_delay = {0, gain};
    b_delay = {gain, 0};
    c_delay = {0, gain};
    d_delay = {0, 0};
    dim_state = cellfun(@(mats) size(mats, 1), a_delay);
    timestep = -1;
    horizon_period = [1, 1];
    lft_delay = Ulft(a_delay, b_delay, c_delay, d_delay,...
                     DeltaDelayZ(dim_state, timestep, horizon_period),...
                     'horizon_period', horizon_period);
    options = AnalysisOptions('verbose', false, 'lmi_shift', 1e-6);    
    
    % For gain = 1, performance should always be 1, regardless of final_time
    final_time = 2;
    reach_lft = generateReachabilityLft(lft_delay, final_time);
    [result, valid] = iqcAnalysis(reach_lft,'analysis_options',options); 
    assertTrue(testCase, valid)
    verifyLessThan(testCase, abs(result.performance - gain^3), 1e-2)
    
    final_time = 15;
    reach_lft = generateReachabilityLft(lft_delay, final_time);
    [result, valid] = iqcAnalysis(reach_lft,'analysis_options',options); 
    assertTrue(testCase, valid)
    verifyLessThan(testCase, abs(result.performance - gain^3), 1e-2)
    
    final_time = 100;
    reach_lft = generateReachabilityLft(lft_delay, final_time);
    [result, valid] = iqcAnalysis(reach_lft,'analysis_options',options); 
    assertTrue(testCase, valid)
    verifyLessThan(testCase, abs(result.performance - gain^3), 1e-2)
end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)