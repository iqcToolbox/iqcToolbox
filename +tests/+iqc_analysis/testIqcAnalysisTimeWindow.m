%% Requirements:
%  1. IQC analysis shall produce an upper-bound on worst-case performance
%     for many systems with disturbance signals that are windowed in time.
%     Producing a tight upper-bound for ALL uncertain systems with 
%     time-windowed disturbances is not expected.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for IQC analysis with time-windowed disturbances
classdef testIqcAnalysisTimeWindow < matlab.unittest.TestCase
methods (Test)

function testTimeWindowAllSteps(testCase)
    % Define nominal system 
    g = ss(zpk([], 0.9 * rand, 1, -1));
    norm_g = norm(g, 'inf');
    l = toLft(g);
    
    % Define disturbance and lft
    dis = DisturbanceTimeWindow('tw', {1}, 0, [0, 1]);
    lft_g = Ulft(l.a, l.b, l.c, l.d, l.delta,...
                 'disturbance', dis);
        
    opts = AnalysisOptions('verbose', false,...
                           'solver', 'sdpt3',...
                           'lmi_shift', 1e-7);
    [result, valid] = iqcAnalysis(lft_g, 'analysis_options', opts);
    assertTrue(testCase, valid)
    verifyLessThan(testCase, abs(result.performance - norm_g), 1e-4)
end

function testTimeWindowExtended(testCase)
    % Define nominal system 
    g = ss(zpk([], 0.9 * rand, 1, -1));
    norm_g = norm(g, 'inf');
    l = toLft(g);
    
    % Define disturbance and lft
    dis = DisturbanceTimeWindow('tw', {1}, 0, [0, 1]);
    lft_g = Ulft(l.a, l.b, l.c, l.d, l.delta,...
                 'disturbance', dis);
    % Extend horizon_period of lft
    lft_g = matchHorizonPeriod(lft_g, [randi(20), randi(20)]);
    
    opts = AnalysisOptions('verbose', false,...
                           'solver', 'sdpt3',...
                           'lmi_shift', 1e-7);
    [result, valid] = iqcAnalysis(lft_g, 'analysis_options', opts);
    perf_complete_window = result.performance;
    assertTrue(testCase, valid)
    verifyLessThan(testCase, abs(perf_complete_window - norm_g), 1e-4)
    
    % Remove signals during even time steps
    dis_odd = lft_g.disturbance.disturbances{1};
    dis_odd.window = dis_odd.window(1:2:end);
    lft_g_odd = Ulft(lft_g.a, lft_g.b, lft_g.c, lft_g.d, lft_g.delta,...
                     'horizon_period', lft_g.horizon_period,...
                     'disturbance', dis_odd);
    [result, valid] = iqcAnalysis(lft_g_odd, 'analysis_options', opts);
    assertTrue(testCase, valid)
    verifyLessThan(testCase, result.performance, perf_complete_window)    
end

function testTimeWindowReachability(testCase)
    % Define nominal system 
    g = ss(zpk([], 0.9 * rand, 1, -1));
    g = blkdiag(g, g / 2, g);
    % Define disturbance and lft
    dis = DisturbanceTimeWindow('tw', {[1; 3]}, 0, [0, 1]);
    l = generateReachabilityLft(toLft(g), randi(20));
    lft_dis = Ulft(l.a, l.b, l.c, l.d, l.delta,...
                   'horizon_period', l.horizon_period,...
                   'disturbance', dis);
    lft_no_dis = l;
        
    opts = AnalysisOptions('verbose', false,...
                           'solver', 'sdpt3',...
                           'lmi_shift', 1e-7);
    % Compare analysis results w/ and w/o disturbance (should be equal)
    [result, valid] = iqcAnalysis(lft_no_dis, 'analysis_options', opts);
    perf_no_dis = result.performance;
    assertTrue(testCase, valid)
    
    [result, valid] = iqcAnalysis(lft_dis, 'analysis_options', opts);
    perf_dis = result.performance;
    assertTrue(testCase, valid)
    verifyLessThanOrEqual(testCase, abs(perf_dis - perf_no_dis), 1e-4)
    
    % Remove latter part of signal
    dis_half = lft_dis.disturbance.disturbances{1};
    dis_half.window = dis_half.window(1 : floor(end / 2));
    lft_half = Ulft(lft_dis.a, lft_dis.b, lft_dis.c, lft_dis.d, lft_dis.delta,...
                    'horizon_period', lft_dis.horizon_period,...
                    'disturbance', dis_half);
    [result, valid] = iqcAnalysis(lft_half, 'analysis_options', opts);
    perf_dis_half = result.performance;
    assertTrue(testCase, valid)
    verifyLessThan(testCase, perf_dis_half, perf_dis * 1.001)                
    
    % Set disturbance multiplier to have time-invariant decision vars
    dim_in = size(lft_half, 2);
    dis_mult = disturbanceToMultiplier(dis_half,...
                                       'quad_time_varying', false,...
                                       'dim_in_lft', dim_in);
    [result, valid] = iqcAnalysis(lft_half,...
                                  'analysis_options', opts,...
                                  'multipliers_disturbance', dis_mult);
    perf_dis_half = result.performance;
    assertTrue(testCase, valid)
    verifyLessThan(testCase, perf_dis_half, perf_dis * 1.001)  
end

function testTimeWindowSomeChannels(testCase)
    % Generate random ss
    g = blkdiag(drss(3, 1, 1), 0);
    while(max(abs(eig(g.a))) >= 0.99)
        g = blkdiag(drss(3, 1, 1), 0);
    end
    g_norm = norm(g, 'inf');
    g_lft = toLft(g);
    
    % Make a LFT which has no disturbance impact for first timestep, but
    % on all subsequent ones
    g_lft = matchHorizonPeriod(g_lft, [1, 1]);
    g_lft.b{1} = zeros(3, 2);
    g_lft.d{1} = zeros(2);
    
    % Check that iqc analysis matches typical hinfnorm if dis at all times
    channel_first = {1};
    time_window_all = [0, 1];
    horizon_period = [1, 1];
    dis1 = DisturbanceTimeWindow('first at all times',...
                                channel_first,...
                                time_window_all,...
                                horizon_period);
    g_lft_d1 = addDisturbance(g_lft, {dis1});
    opts = AnalysisOptions('verbose', false);
    result = iqcAnalysis(g_lft_d1, 'analysis_options', opts);
    norm_percent_error = 100 * (g_norm - result.performance) / g_norm;
    verifyTrue(testCase, result.valid)
    verifyLessThan(testCase, norm_percent_error, 1e-1)

    % Check that iqc analysis matches typical hinfnorm if dis at all times
    channel_all = {[]};
    dis2 = DisturbanceTimeWindow('all at all times',...
                                channel_all,...
                                time_window_all,...
                                horizon_period);
    g_lft_d2 = addDisturbance(g_lft, {dis2});
    result = iqcAnalysis(g_lft_d2, 'analysis_options', opts);
    norm_percent_error = 100 * (g_norm - result.performance) / g_norm;
    verifyTrue(testCase, result.valid)
    verifyLessThan(testCase, norm_percent_error, 1e-1)
    
    % Check that performance is zero if window only allowed on first
    % timestep
    time_window_first = [0];
    dis3 = DisturbanceTimeWindow('all at first time',...
                                channel_all,...
                                time_window_first,...
                                horizon_period);
    g_lft_d3 = addDisturbance(g_lft, {dis3});
    result = iqcAnalysis(g_lft_d3, 'analysis_options', opts);
    verifyTrue(testCase, result.valid)
    verifyLessThan(testCase, result.performance, 1e-3)
    
    % Check that performance is zero if window only allowed on first
    % timestep
    time_window_first = [0];
    dis4 = DisturbanceTimeWindow('first at first time',...
                                channel_first,...
                                time_window_first,...
                                horizon_period);
    g_lft_d4 = addDisturbance(g_lft, {dis4});
    result = iqcAnalysis(g_lft_d4, 'analysis_options', opts);
    verifyTrue(testCase, result.valid)
    verifyLessThan(testCase, result.performance, 1e-3)

end


function testTimeWindowSomeChannelsWithL2Disturbance(testCase)
    % Generate random ss
    g = blkdiag(drss(3, 1, 1), 0);
    while(max(abs(eig(g.a))) >= 0.99)
        g = blkdiag(drss(3, 1, 1), 0);
    end
    g_norm = norm(g, 'inf');
    g_lft = toLft(g);
    
    % Make a LFT which has no disturbance impact for first timestep, but
    % on all subsequent ones
    g_lft = matchHorizonPeriod(g_lft, [1, 1]);
    g_lft.b{1} = zeros(3, 2);
    g_lft.d{1} = zeros(2);
    
    % Check that iqc analysis matches typical hinfnorm if dis at all times
    channel_first = {1};
    time_window_all = [0, 1];
    horizon_period = [1, 1];
    dis1 = DisturbanceTimeWindow('first at all times',...
                                channel_first,...
                                time_window_all,...
                                horizon_period);
    dis_l2_one = DisturbanceL2('def', {1});
    dis_l2_all = DisturbanceL2('def', {[]});
    g_lft_d1 = addDisturbance(g_lft, {dis1, dis_l2_one});
    opts = AnalysisOptions('verbose', false);
    result = iqcAnalysis(g_lft_d1, 'analysis_options', opts);
    norm_percent_error = 100 * (g_norm - result.performance) / g_norm;
    verifyTrue(testCase, result.valid)
    verifyLessThan(testCase, norm_percent_error, 1e-1)

    % Check that iqc analysis matches typical hinfnorm if dis at all times
    channel_all = {[]};
    dis2 = DisturbanceTimeWindow('all at all times',...
                                channel_all,...
                                time_window_all,...
                                horizon_period);
    g_lft_d2 = addDisturbance(g_lft, {dis2, dis_l2_all});
    result = iqcAnalysis(g_lft_d2, 'analysis_options', opts);
    norm_percent_error = 100 * (g_norm - result.performance) / g_norm;
    verifyTrue(testCase, result.valid)
    verifyLessThan(testCase, norm_percent_error, 1e-1)
    
    % Check that performance is zero if window only allowed on first
    % timestep
    time_window_first = [0];
    dis3 = DisturbanceTimeWindow('all at first time',...
                                channel_all,...
                                time_window_first,...
                                horizon_period);
    g_lft_d3 = addDisturbance(g_lft, {dis3, dis_l2_one});
    result = iqcAnalysis(g_lft_d3, 'analysis_options', opts);
    verifyTrue(testCase, result.valid)
    verifyLessThan(testCase, result.performance, 1e-3)
    
    % Check that performance is zero if window only allowed on first
    % timestep
    time_window_first = [0];
    dis4 = DisturbanceTimeWindow('first at first time',...
                                channel_first,...
                                time_window_first,...
                                horizon_period);
    g_lft_d4 = addDisturbance(g_lft, {dis4, dis_l2_all});
    result = iqcAnalysis(g_lft_d4, 'analysis_options', opts);
    verifyTrue(testCase, result.valid)
    verifyLessThan(testCase, result.performance, 1e-3)

end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)