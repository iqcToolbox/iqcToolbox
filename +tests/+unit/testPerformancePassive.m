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

