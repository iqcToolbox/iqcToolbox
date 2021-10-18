%% Requirements:
%  1. PerformancePassive shall construct a valid passive performance metric,
%      given the name, channels in, channels out, and horizon_period

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for PerformanceL2Induced
classdef testPerformancePassive < matlab.unittest.TestCase
methods (Test)
    
function testBadPerformanceConstruction(testCase)
    verifyError(testCase,...
                @() PerformancePassive(),...
                'PerformancePassive:PerformancePassive')

    verifyError(testCase,...
                @() PerformancePassive('test', {[1;2]}),...
                'PerformancePassive:PerformancePassive')
            
    verifyError(testCase,...
                @() PerformancePassive('test', {[1;2]}, {[1]}),...
                'PerformancePassive:PerformancePassive')
end

function testBadHorizonPeriod(testCase)
    horizon_period = [0, 3];
    good_chan = {1};
    bad_chan = {1, 1};
    verifyError(testCase,...
                @() PerformancePassive('test', good_chan, bad_chan, [0, 1]),...
                'PerformancePassive:matchHorizonPeriod')
    verifyError(testCase,...
                @() PerformancePassive('test', bad_chan, good_chan, [0, 1]),...
                'PerformancePassive:matchHorizonPeriod')
end

function testBadPerformanceMultiplierConstruction(testCase)
    dim_outin = 2;
    channel_out_of_dims = 3;
    perf = PerformancePassive('test', {channel_out_of_dims}, {1});
    verifyError(testCase,...
                @() MultiplierPerformancePassive(perf, dim_outin, dim_outin),...
                'MultiplierPerformancePassive:MultiplierPerformancePassive')
    perf = PerformancePassive('test', {1}, {channel_out_of_dims});
    verifyError(testCase,...
                @() MultiplierPerformancePassive(perf, dim_outin, dim_outin),...
                'MultiplierPerformancePassive:MultiplierPerformancePassive')
end

function testFullConstructorAndDisplay(testCase)
    name = 'test';
    chan_out = {1};
    chan_in = {2};
    horizon_period = [3, 2];
    total_time = sum(horizon_period);
    perf = PerformancePassive(name, chan_out, chan_in, horizon_period)
    verifyEqual(testCase, perf.name, name);
    verifyEqual(testCase, perf.chan_out, repmat(chan_out, 1, total_time));
    verifyEqual(testCase, perf.chan_in, repmat(chan_in, 1, total_time));
    verifyEqual(testCase, perf.horizon_period, horizon_period)
end

function testPerformanceToMultiplier(testCase)
    perf = PerformancePassive('test');
    verifyError(testCase,...
                @() performanceToMultiplier(perf),...
                'PerformancePassive:performanceToMultiplier')
    verifyError(testCase,...
                @() performanceToMultiplier(perf, 'dim_out_lft', 1),...
                'PerformancePassive:performanceToMultiplier')
end

end
end

%%  CHANGELOG
% Oct. 18, 2021: Added after v0.6.0 - Micah Fry (micah.fry@ll.mit.edu)