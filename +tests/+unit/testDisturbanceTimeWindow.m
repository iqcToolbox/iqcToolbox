%% Requirements:
%  1. DisturbanceTimeWindow shall be defined by it's name, the input
%      channels of interest, the window of time-instances when the
%      signal of interest is non-zero, and the horizon_period of the
%      channels and time-window
%  2. Upon construction, and when queried by user, it shall display the
%      information described in (1).
%
%  3. If dimension and/or bound information is not provided by the user, by
%      default the channel shall be {1}, the window shall be [0], and the
%      horizon_period shall be [0, 1].
%
%  4. If the user provides no name, DisturbanceTimeWindow shall throw an 
%      exception
%  5. If the user provides a channel that is not a cell of vectors of
%      natural numbers DisturbanceTimeWindow shall throw an exception
%  6. If the user provides a window and horizon_period that are
%      inconsistent with each other, DisturbanceTimeWindow shall throw an
%      exception
%  7. If the user provides a window with duplicate time indices,
%      DisturbanceTimeWindow shall throw an exception
%
%  7. DisturbanceTimeWindow shall ensure that it's properties are consistent 
%      with its current horizon_period property
%  8. DisturbanceTimeWindow shall be capable of changing it's properties to 
%      match a newly input horizon_period, as long as the new 
%      horizon_period is consistent with the prior horizon_period
%
%  9. DisturbanceTimeWindow shall be capable of generating a
%      MultiplierTimeWindow from a DisturbanceTimeWindow object

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for DisturbanceTimeWindow.
classdef testDisturbanceTimeWindow < matlab.unittest.TestCase
methods (Test)

function testFullConstructor(testCase)
    name = 'test';
    channel = {[1:5]'};
    window = [0];
    horizon_period = [0, 1];
    dis = DisturbanceTimeWindow(name, channel, window, horizon_period);
    verifyEqual(testCase, dis.name, name)
    verifyEqual(testCase, dis.chan_in, channel)
    verifyEqual(testCase, dis.window, window)
    verifyEqual(testCase, dis.horizon_period, dis.horizon_period)
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

function testDisturbanceToMultiplier(testCase)
    dis = DisturbanceTimeWindow('test');
    dim_in = 1;
    mult = disturbanceToMultiplier(dis, 'dim_in_lft', dim_in);
    verifyClass(testCase, mult, 'MultiplierTimeWindow')
    verifyError(testCase,...
                @() disturbanceToMultiplier(dis),...
                ?MException)
end

function testDisplayDisturbance(testCase)
    dis = DisturbanceTimeWindow('test')
    dis = DisturbanceTimeWindow('test', {1,2,3,4,5,6,7}, [0:6], [0, 7])
end

function testFilterLftDisturbance(testCase)
    chan = {[1]};
    window = 1;
    horizon_period = [0, 2];
    dis = DisturbanceTimeWindow('test', chan, window, horizon_period);
    mult = MultiplierTimeWindow(dis, [2, 2]);
    filter = mult.filter_lft;
    verifyEqual(testCase, filter.d, {eye(2), zeros(2)})
end

end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)