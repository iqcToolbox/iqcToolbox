%% Requirements:
%  1. IQC analysis shall be capable of producing worst-case upper-bounds on
%     uncertain-systems which have disturbances constrained to be constant in
%     pre-specified windows. 

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for IQC analysis with Constant Disturbances
classdef testIqcAnalysisDisturbanceConstant < matlab.unittest.TestCase
methods (Test)

function testReachabilityWithConstantSignal(testCase)
    zero = [];
    pole = -.5;
    gain = 1;
    timestep = -1;
    g = ss(zpk(zero, pole, gain, timestep));
%     g = drss(3, 1, 1)
%     g.a = g.a * 0.9;
    lft = toLft(g);
    final_time = 5;
    lft_reach = generateReachabilityLft(lft, final_time);
    window = 0:final_time;
    d = DisturbanceConstantWindow('dis',{[]}, window, lft_reach.horizon_period);
    lft_reach = lft_reach.addDisturbance({d});
    options = AnalysisOptions('verbose', true, 'lmi_shift', 1e-7);
    [result, valid] = iqcAnalysis(lft_reach, 'analysis_options', options);
    assertTrue(testCase, valid)

    u_sim = ones(1, final_time + 1);
    y_sim = simulate(lft_reach, u_sim);
    perf_sim = norm(y_sim, 'fro') / norm(u_sim, 'fro');
    perf_diff = result.performance - perf_sim;
    testCase.verifyGreaterThan(perf_diff, 0)
    testCase.verifyLessThan(abs(perf_diff)/perf_sim, 1e-3)
end

function test
end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)