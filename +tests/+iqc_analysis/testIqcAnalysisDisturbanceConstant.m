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
    % Check analysis result against simulation result
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
    window = 1 : final_time;
    d = DisturbanceConstantWindow('dis',{[]}, window, lft_reach.horizon_period);
    lft_reach = lft_reach.addDisturbance({d});
    options = AnalysisOptions('verbose', false, 'lmi_shift', 1e-5);
    [result, valid] = iqcAnalysis(lft_reach, 'analysis_options', options);
    assertTrue(testCase, valid)

    u_sim = ones(1, final_time + 1);
    y_sim = simulate(lft_reach, u_sim);
    perf_sim = norm(y_sim, 'fro') / norm(u_sim, 'fro');
    perf_diff = result.performance - perf_sim;
    testCase.verifyGreaterThan(perf_diff, 0)
    testCase.verifyLessThan(abs(perf_diff)/perf_sim, 1e-3)
    
    % No impact for memory-less systems (if constant portion irrelevant to reach window)
    rct_object = zeros(3);
    for i = 1:3
        var = randatom('ureal');
        base = rand(3);
        base(base < .5) = 0;
        base(base >= .5) = 1;
        rct_object = rct_object + var * base;
    end  
    rct_object = uss(rct_object);
    rct_result = wcgain(rct_object);
    assumeTrue(testCase, isfinite(rct_result.LowerBound))
    assumeTrue(testCase, isfinite(rct_result.UpperBound))
    lft = rctToLft(rct_object);
    lft = lft - zeros(3) * DeltaDelayZ(3);
    % Make a reachability LFT to allow specification of Disturbance to be constant over a subset of time
    final_time = 3;
    lft_reach = generateReachabilityLft(lft, final_time);
%     lft_reach = lft_reach.matchHorizonPeriod(lft_reach.horizon_period + [1, 0]);
    % This constant window should be before the reachability timestep, which means that it should not reduce the upper bound
    d = DisturbanceConstantWindow('d',...
                                  {[]},...
                                  1:(final_time-1),...
                                  lft_reach.horizon_period);
    lft_reach = lft_reach.addDisturbance({d});
    result = iqcAnalysis(lft_reach, 'analysis_options', options);
    verifyGreaterThan(testCase, result.performance, rct_result.LowerBound * .99)
    verifyLessThan(testCase, result.performance, rct_result.UpperBound * 1.01)
    
    % Now, if constant portion exists during the reachability timestep, upper bound will decrease
    lft_reach = lft_reach.removeDisturbance(1);
%     lft_reach = lft_reach.matchHorizonPeriod(lft_reach.horizon_period + [1, 0]);
    d = DisturbanceConstantWindow('d',...
                                  {[]},...
                                  final_time,...
                                  lft_reach.horizon_period);
    lft_reach = lft_reach.addDisturbance({d});
    result2 = iqcAnalysis(lft_reach, 'analysis_options', options);
    testCase.verifyLessThan(result2.performance, result.performance)
end

function testFullWindowIsReducedResult(testCase)
     % Check analysis result against simulation result
    zero = [];
    pole = -.5;
    gain = 1;
    timestep = -1;
    g = ss(zpk(zero, pole, gain, timestep));
%     g = drss(3, 1, 1)
%     g.a = g.a * 0.9;
    lft = toLft(g);
    lft = lft.addDisturbance({DisturbanceConstantWindow('d')});
    options = AnalysisOptions('verbose', false, 'lmi_shift', 1e-7);
    [result, valid] = iqcAnalysis(lft, 'analysis_options', options);
    testCase.assertTrue(valid)
    true_performance = norm(g, 'inf');
    testCase.verifyLessThan(result.performance, true_performance)  
    
    lft = lft.matchHorizonPeriod([7, 10]);
    [result2, valid] = iqcAnalysis(lft, 'analysis_options', options);
    testCase.assertTrue(valid)
    diff_perf = abs(result.performance - result2.performance);
    testCase.verifyLessThan(diff_perf/result.performance, 1e-3)  
end

function testPartiallyConstantThroughPeriodIsImprovement(testCase)
    filter_order = 4;
    cutoff_freq = 0.8;
    [z, p, k] = butter(filter_order, cutoff_freq, 'high');
    high_pass_ss = ss(zpk(z, p, k, -1));
    high_pass = toLft(high_pass_ss);
    options = AnalysisOptions('verbose', false, 'lmi_shift', 1e-7);
    [result_hpf, valid] = iqcAnalysis(high_pass, 'analysis_options', options);
    testCase.assertTrue(valid)
    
    % Add constant disturbance with short window, which should have little impact on upper-bound
    period = 20;
    high_pass = matchHorizonPeriod(high_pass, [0, period]);
    window = 1;
    d = DisturbanceConstantWindow('d', {[]}, window, high_pass.horizon_period);
    high_pass = high_pass.addDisturbance({d});
    result_hpf_fast = iqcAnalysis(high_pass, 'analysis_options', options);
    testCase.assertTrue(result_hpf_fast.valid)
    diff_perf = abs(result_hpf.performance - result_hpf_fast.performance);
    testCase.verifyLessThan(diff_perf / result_hpf.performance, 1e0);
    % Now give constant property over large window, which should significantly reduce upper-bound
    window = 1:(period - 1);
    d = DisturbanceConstantWindow('d', {[]}, window, high_pass.horizon_period);
    high_pass = high_pass.removeDisturbance(1).addDisturbance({d});
    result_hpf_slow = iqcAnalysis(high_pass, 'analysis_options', options);
    testCase.assertTrue(result_hpf_slow.valid)
    diff_perf = abs(result_hpf.performance - result_hpf_slow.performance);
    testCase.verifyGreaterThan(diff_perf / result_hpf.performance, 0.5); % at least 50% reduction
end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)