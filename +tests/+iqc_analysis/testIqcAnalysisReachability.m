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

function testReachabilityNonTrivialPeriod(testCase)
    % This test would have failed previous to hotfix-009
    timestep = rand + 1e-8; % Creata specific timestep (be sure it is not -1)
    period = randi([2, 10]); % Have greater-than-one period
    horizon_period = [randi([0, 10]), period];
    dim_state = randi([1, 3]);
    del = DeltaDelayZ(dim_state, timestep, horizon_period);
    lft = Ulft.random('num_deltas', 1, 'req_deltas', {del});
    final_time = horizon_period(1) + randi([1, 10]); % A final_time greater than the horizon
    lft_reach = generateReachabilityLft(lft, final_time);
    options = AnalysisOptions('verbose', false, 'lmi_shift', 1e-5);
    [result, valid] = iqcAnalysis(lft_reach, 'analysis_options', options);
    assertTrue(testCase, valid)
    
    % Create an equivalent reachability LFT
    hp_equiv = [final_time + 1, 1];
    hp_common = commonHorizonPeriod([hp_equiv; horizon_period]);
    ind = makeNewIndices(horizon_period, hp_common);
    total_time = sum(hp_equiv);
    ind = ind(1:total_time); % Truncate new indices to hp_equiv
    a = lft.a(1, ind);
    b = lft.b(1, ind);
    c = lft.c(1, ind);
    d = lft.d(1, ind);
    for i = 1:total_time
        if i < hp_equiv(1)
            c{i} = zeros(size(c{i}));
            d{i} = zeros(size(d{i}));
        elseif i > hp_equiv(1)
            a{i} = zeros(size(a{i}));
            b{i} = zeros(size(b{i}));
            c{i} = zeros(size(c{i}));
            d{i} = zeros(size(d{i}));
        end
    end
    del_equiv = DeltaDelayZ(dim_state, timestep, hp_equiv);
    lft_equiv = Ulft(a, b, c, d, del_equiv, 'horizon_period', hp_equiv);
    result_equiv = iqcAnalysis(lft_equiv, 'analysis_options', options);
    assertTrue(testCase, result_equiv.valid)
    perf_error = abs((result.performance - result_equiv.performance)/...
                     result_equiv.performance);
    verifyLessThan(testCase, perf_error, 1e-4)  
end
    
function testReachabilityNonTrivialPeriodShortHorizon(testCase)
    % This test would have failed previous to hotfix-009
    % While the previous tested final times somewhere beyond the horizon of hp,
    % this test is when the final time is less than the horizon of hp
    % (effectively truncating the system dynamics)
    timestep = -1;
    period = randi([2, 10]); % Have greater-than-one (non-trivial) period
    horizon_period = [randi([4, 10]), period]; % Have a non-zero (non-trivial) horizon
    dim_state = randi([1, 3]);
    del = DeltaDelayZ(dim_state, timestep, horizon_period);
    lft = Ulft.random('num_deltas', 1, 'req_deltas', {del});
    final_time = randi([0, horizon_period(1) - 2]); % Have a final timestep before the end of the horizon
    lft_reach = generateReachabilityLft(lft, final_time);
    options = AnalysisOptions('verbose', false, 'lmi_shift', 1e-5);
    [result, valid] = iqcAnalysis(lft_reach, 'analysis_options', options);
    assertTrue(testCase, valid)
    
    % Create an equivalent reachability LFT
    hp_equiv = [final_time + 1, 1];
    hp_common = commonHorizonPeriod([hp_equiv; horizon_period]);
    ind = makeNewIndices(horizon_period, hp_common);
    total_time = sum(hp_equiv);
    ind = ind(1:total_time); % Truncate new indices to hp_equiv
    a = lft.a(1, ind);
    b = lft.b(1, ind);
    c = lft.c(1, ind);
    d = lft.d(1, ind);
    for i = 1:total_time
        if i < hp_equiv(1)
            c{i} = zeros(size(c{i}));
            d{i} = zeros(size(d{i}));
        elseif i > hp_equiv(1)
            a{i} = zeros(size(a{i}));
            b{i} = zeros(size(b{i}));
            c{i} = zeros(size(c{i}));
            d{i} = zeros(size(d{i}));
        end
    end
    del_equiv = DeltaDelayZ(dim_state, timestep, hp_equiv);
    lft_equiv = Ulft(a, b, c, d, del_equiv, 'horizon_period', hp_equiv);
    result_equiv = iqcAnalysis(lft_equiv, 'analysis_options', options);
    assertTrue(testCase, result_equiv.valid)
    perf_error = abs((result.performance - result_equiv.performance)/...
                     result_equiv.performance);
    verifyLessThan(testCase, perf_error, 1e-4)
end
end

methods (Test, TestTags = {'RCT'})
function testReachabilityMemoryless(testCase)
    % Reachability analysis should be the same for memoryless systems
    rng(1, 'twister')
    dim_outin = 2;
    rct_object = zeros(dim_outin);
    for i = 1:3
        var = randatom('ureal');
        base = rand(dim_outin);
        base(base < .5) = 0;
        base(base >= .5) = 1;
        rct_object = rct_object + var * base;
    end  
    rct_object = uss(rct_object);
    rct_result = wcgain(rct_object);
    assumeTrue(testCase, isfinite(rct_result.LowerBound))
    assumeTrue(testCase, isfinite(rct_result.UpperBound))
    lft = rctToLft(rct_object);
    lft = lft + zeros(dim_outin) * DeltaDelayZ(dim_outin);
    options = AnalysisOptions('verbose', false, 'lmi_shift', 1e-7);
    result = iqcAnalysis(lft, 'analysis_options', options);
    % Ensure IQC analysis is returning correct results for LTI system
    verifyGreaterThan(testCase, result.performance, rct_result.LowerBound * .99)
    verifyLessThan(testCase, result.performance, rct_result.UpperBound * 1.01)
    
    lft_reach = generateReachabilityLft(lft, 0);
    result2 = iqcAnalysis(lft_reach, 'analysis_options', options);
    % Ensure IQC analysis is returning correct results for LTI system
    diff_perf = abs(result.performance - result2.performance);
    testCase.verifyLessThan(diff_perf / result.performance, 1e-4)
    
    lft_reach = generateReachabilityLft(lft, 15);
    result2 = iqcAnalysis(lft_reach, 'analysis_options', options);
    % Ensure IQC analysis is returning correct results for LTI system
    diff_perf = abs(result.performance - result2.performance);
    testCase.verifyLessThan(diff_perf / result.performance, 1e-4)
end
end

end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)