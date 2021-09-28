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

function testFullConstructorDifferentHorizonPeriod(testCase)
    name = 'test';
    channel = {[1:5]'};
    horizon_period = [2, 4];
    total_time = sum(horizon_period);
    window = (1:total_time) - 1;
    dis = DisturbanceTimeWindow(name, channel, window, horizon_period);
    verifyEqual(testCase, dis.name, name)
    verifyEqual(testCase,...
                dis.chan_in,...
                repmat(channel, 1, total_time))
    verifyEqual(testCase,...
                dis.window,...
                window)
    verifyEqual(testCase, dis.horizon_period, horizon_period)
end

function testOneArgConstructor(testCase)
    name = 'test';
    dis = DisturbanceTimeWindow(name);
    verifyEqual(testCase, dis.name, name)
    verifyEqual(testCase, dis.chan_in, {[]})
    verifyEqual(testCase, dis.window, 0)
    verifyEqual(testCase, dis.horizon_period, [0, 1])
end

function testMatchHorizonPeriod(testCase)
   name = 'test';
   chan = {[1; 4; 10]};
   win = [0, 2, 4];
   horizon_period = [1, 4];
   total_time = sum(horizon_period);
   dis = DisturbanceTimeWindow(name, chan, win, horizon_period);
   
   assertEqual(testCase, dis.horizon_period, horizon_period);
   assertEqual(testCase, dis.chan_in, repmat(chan, 1, total_time));
   assertEqual(testCase, dis.window, win);
   
   % Checking horizon_period and making sure it fits for all properties
   dis2 = dis;
   horizon_period2 = [5, 8];
   total_time2 = sum(horizon_period2);
   win2 = [0, 2, 4, 6:2:13];
   
   dis = matchHorizonPeriod(dis, horizon_period2);
   verifyEqual(testCase, dis.horizon_period, horizon_period2)
   verifyEqual(testCase, dis.chan_in, repmat(chan, 1, total_time2))
   verifyEqual(testCase, dis.window, win2)
end

function testFailedName(testCase)
    verifyError(testCase, @() DisturbanceTimeWindow(), ?MException)
    verifyError(testCase, @() DisturbanceTimeWindow(1), ?MException)
end

function testFailedNumberOfArguments(testCase)
    verifyError(testCase, @() DisturbanceTimeWindow('test', {1}), ?MException)
    verifyError(testCase, @() DisturbanceTimeWindow('test', {1}, 0), ?MException)
end

function testFailedChannels(testCase)
    name = 'test';
    hp = [0, 1];
    win = 0;
    verifyError(testCase, @() DisturbanceTimeWindow(name, 1, win, hp), ?MException)
    verifyError(testCase, @() DisturbanceTimeWindow(name, {1.1}, win, hp), ?MException)
    verifyError(testCase, @() DisturbanceTimeWindow(name, {-1}, win, hp), ?MException)
    verifyError(testCase, @() DisturbanceTimeWindow(name, {0}, win, hp), ?MException)
    verifyError(testCase, @() DisturbanceTimeWindow(name, {[1, 2]}, win, hp), ?MException)
end

function testFailedWindow(testCase)
    name = 'test';
    hp = [0, 1];
    chan = {1};
    verifyError(testCase, @() DisturbanceTimeWindow(name, chan, -1, hp), ?MException)
    verifyError(testCase, @() DisturbanceTimeWindow(name, chan, 1.1, hp), ?MException)
    verifyError(testCase, @() DisturbanceTimeWindow(name, chan, [0, 1, 2], hp), ?MException)
    verifyError(testCase, @() DisturbanceTimeWindow(name, chan, [0, 0], hp), ?MException)
end

function testFailedRate(testCase)
    hp = [2, 2];
    total_time = sum(hp);
    lr = -1.2 * ones(1,total_time);
    ur = linspace(4, 5, total_time);
    lr(total_time) = ur(total_time) + 1;

    verifyError(testCase,...
                @() DisturbanceTimeWindow('test', 1, -1, 1, lr, ur, hp),...
                ?MException)

    lr(total_time) = -inf;
    verifyError(testCase,...
                @() DisturbanceTimeWindow('test', 1, -1, 1, lr, ur, hp),...
                ?MException)

    lr(total_time) = ur(total_time) - 1;
    ur(1) = nan;
    verifyError(testCase,...
                @() DisturbanceTimeWindow('test', 1, -1, 1, lr, ur, hp),...
                ?MException)
end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)