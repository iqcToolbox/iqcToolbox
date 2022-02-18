%% Requirements:
%  1. PerformanceStable shall construct a valid stability performance metric,
%      given the horizon_period

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for PerformanceL2Induced
classdef testPerformanceStable < matlab.unittest.TestCase
methods (Test)

function testFullConstructorAndDisplay(testCase)
    horizon_period = [3, 2];
    total_time = sum(horizon_period);
    perf = PerformanceStable(horizon_period)
    verifyEqual(testCase, perf.name, 'stability');
    verifyEqual(testCase, perf.chan_out, repmat({0}, 1, total_time));
    verifyEqual(testCase, perf.chan_in, repmat({0}, 1, total_time));
    verifyEqual(testCase, perf.horizon_period, horizon_period)
end

function testOneArgumentConstructor(testCase)
    perf = PerformanceStable();
    verifyEqual(testCase, perf.name, 'stability');
    verifyEqual(testCase, perf.chan_out, {0});
    verifyEqual(testCase, perf.chan_in, {0});
    verifyEqual(testCase, perf.horizon_period, [0, 1])
end

function testPerformanceToMultiplier(testCase)
    p = PerformanceStable();
    verifyError(testCase,...
                @() p.performanceToMultiplier(),...
                'PerformanceStable:performanceToMultiplier')
end

function testMatchHorizonPeriod(testCase)
    perf = PerformanceStable();
    hp_new = [2, 2];
    total_time = sum(hp_new);
    perf = perf.matchHorizonPeriod(hp_new)
    testCase.verifyEqual(perf.horizon_period, hp_new)
    testCase.verifyEqual(perf.chan_in, repmat({0}, 1, total_time))
    testCase.verifyEqual(perf.chan_out, repmat({0}, 1, total_time))
    perf = PerformanceStable(hp_new);
    perf.chan_out = {0, 0};
    verifyError(testCase,...
                @() perf.matchHorizonPeriod(),...
                'PerformanceStable:matchHorizonPeriod')
    perf.chan_out = {0};
    perf.chan_in = {0, 0, 0};
    verifyError(testCase,...
                @() perf.matchHorizonPeriod(),...
                'PerformanceStable:matchHorizonPeriod')
end
end
end

%%  CHANGELOG
% Oct. 18, 2021: Added after v0.6.0 - Micah Fry (micah.fry@ll.mit.edu)