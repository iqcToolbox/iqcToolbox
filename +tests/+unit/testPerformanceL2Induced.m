% Progressive tests on PerformanceL2Induced

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for PerformanceL2Induced
classdef testPerformanceL2Induced < matlab.unittest.TestCase
methods (Test)

function testFullConstructor(testCase)
    name = 'test';
    chan_out = {[1; 2]};
    chan_in = {[1; 3; 2]};
    gain = [];
    horizon_period = [0, 1];
    perf = PerformanceL2Induced(name, chan_out, chan_in, gain, horizon_period);
    verifyEqual(testCase, perf.name, name)
    verifyEqual(testCase, perf.chan_in, chan_in)
    verifyEqual(testCase, perf.chan_out, chan_out)
    verifyEqual(testCase, perf.gain, gain)
    verifyEqual(testCase, perf.horizon_period, horizon_period)
    
    chan_out = {[]};
    chan_in = {[]};
    gain = 3;
    perf = PerformanceL2Induced(name, chan_out, chan_in, gain, horizon_period);
    verifyEqual(testCase, perf.name, name)
    verifyEqual(testCase, perf.chan_in, chan_in)
    verifyEqual(testCase, perf.chan_out, chan_out)
    verifyEqual(testCase, perf.gain, gain)
    verifyEqual(testCase, perf.horizon_period, horizon_period)

    chan_out = {};
    chan_in = {};
    perf = PerformanceL2Induced(name, chan_out, chan_in, gain, horizon_period);
    verifyEqual(testCase, perf.name, name)
    verifyEqual(testCase, perf.chan_in, {[]})
    verifyEqual(testCase, perf.chan_out, {[]})
    verifyEqual(testCase, perf.gain, gain)
    verifyEqual(testCase, perf.horizon_period, horizon_period)
end

function testDisplayDisturbance(testCase)
    perf = PerformanceL2Induced('test')
    perf = PerformanceL2Induced('test', {1}, {1}, 3)    
    perf = PerformanceL2Induced('test', {[1;4]}, {}, 2)    
end

function testConstructMultiplier(testCase)
    name = 'test';
    chan_out = {};
    chan_in = {[1; 3; 2]};
    gain = [];
    horizon_period = [0, 1];
    dim_out = 4;
    dim_in = 4;
    perf = PerformanceL2Induced(name, chan_out, chan_in, gain, horizon_period);
    mult = performanceToMultiplier(perf,...
                                   'dim_out_lft', dim_out,...
                                   'dim_in_lft', dim_in);
    verifyEqual(testCase, mult.name, name)
    verifyEqual(testCase, mult.chan_in, chan_in)
    verifyEqual(testCase, mult.chan_out, {[]})
    verifyEqual(testCase, mult.dim_in, dim_in)
    verifyEqual(testCase, mult.dim_out, dim_out)
    verifyClass(testCase, mult.gain, 'sdpvar')
    verifyEqual(testCase, mult.horizon_period, horizon_period) 

    mult = MultiplierL2Induced(perf, dim_out, dim_in);
    verifyEqual(testCase, mult.name, name)
    verifyEqual(testCase, mult.chan_in, chan_in)
    verifyEqual(testCase, mult.chan_out, {[]})
    verifyEqual(testCase, mult.dim_in, dim_in)
    verifyEqual(testCase, mult.dim_out, dim_out)
    verifyClass(testCase, mult.gain, 'sdpvar')
    verifyEqual(testCase, mult.horizon_period, horizon_period)
end

function testSequencePerformance(testCase)
    name = 'test';
    chan_in = {[2;3], [1:3]'};
    chan_out = {};
    gain = 5;
    horizon_period = [1, 1];
    perf = PerformanceL2Induced(name, chan_out, chan_in, gain, horizon_period);
    perfs = SequencePerformance(perf);
    verifyEqual(testCase, perfs.names{1}, name)
    verifyEqual(testCase, perfs.chan_ins(1, :), chan_in)
    verifyEqual(testCase, perfs.chan_outs(1, :), {[], []})
    verifyEqual(testCase, perfs.horizon_periods(1, :), horizon_period)
end

function testBadConstructor(testCase)
    verifyError(testCase,...
                @() PerformanceL2Induced(),...
                'PerformanceL2Induced:PerformanceL2Induced')
    verifyError(testCase,...
                @() PerformanceL2Induced('test', {1}),...
                'PerformanceL2Induced:PerformanceL2Induced')
    perf = PerformanceL2Induced('test', {1}, {1});
    verifyError(testCase,...
                @() performanceToMultiplier(perf),...
                'PerformanceL2Induced:performanceToMultiplier')
    verifyError(testCase,...
                @() performanceToMultiplier(perf, 'dim_out_lft', 1),...
                'PerformanceL2Induced:performanceToMultiplier')
end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)