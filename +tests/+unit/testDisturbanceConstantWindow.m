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

%% Test class for DisturbanceConstantWindow and MultiplierConstantWindow.
classdef testDisturbanceConstantWindow < matlab.unittest.TestCase
methods (Test)

function testDisturbanceFullConstructor(testCase)
    name = 'd';
    chan_in = {2};
    window = [0:2];
    horizon_period = [3, 1];
    override = true;
    d = DisturbanceConstantWindow(name,...
                                  chan_in,...
                                  window,...
                                  horizon_period,...
                                  override);
    testCase.verifyEqual(d.name, 'd')
    testCase.verifyEqual(d.chan_in, repmat(chan_in, 1, sum(horizon_period)))
    testCase.verifyEqual(d.window, window)
    testCase.verifyEqual(d.horizon_period, horizon_period)
    testCase.verifyEqual(d.override, override)
end

function testFourArgConstructor(testCase)
    name = 'd';
    chan_in = {[1;3]};
    window = 2;
    horizon_period = [2, 2];
    d = DisturbanceConstantWindow(name, chan_in, window, horizon_period);
    testCase.verifyEqual(d.name, 'd')
    testCase.verifyEqual(d.chan_in, repmat(chan_in, 1, sum(horizon_period)))
    testCase.verifyEqual(d.window, window)
    testCase.verifyEqual(d.horizon_period, horizon_period)
    testCase.verifyEqual(d.override, false)
end

function testOneArgConstructor(testCase)
    name = 'd';
    d = DisturbanceConstantWindow(name);
    testCase.verifyEqual(d.name, 'd')
    testCase.verifyEqual(d.chan_in, {[]})
    testCase.verifyEqual(d.window, 0)
    testCase.verifyEqual(d.horizon_period, [0, 1])
    testCase.verifyEqual(d.override, false)
end

function testMixingNonperiodicAndPeriodicTimesteps(testCase)
    % Should throw an error without override (because window includes periodic and non-periodic portions)
    window = [0:1];
    horizon_period = [2, 2];
    bad_d = @()DisturbanceConstantWindow('d', {[]}, window, horizon_period);
    testCase.verifyError(bad_d,...
                         'DisturbanceConstantWindow:DisturbanceConstantWindow')
    
    % Now with override, it will make a window that has a periodic portion
    override = true;
    d = DisturbanceConstantWindow('d', {[]}, window, horizon_period, override);
    dim_in = ones(1, sum(d.horizon_period));
    m = MultiplierConstantWindow(d, dim_in);
    q_class = cellfun(@class, m.quad.q, 'UniformOutput', false);
    q_true = {'double', 'sdpvar', 'sdpvar', 'double'};
    testCase.verifyEqual(q_class, q_true)
    
    % Periodicity can be seen by extending the horizon
    extend_horizon = 2;
    new_hp = d.horizon_period + [extend_horizon, 0];
    d = d.matchHorizonPeriod(new_hp);
    dim_in = ones(1, sum(d.horizon_period));
    m = MultiplierConstantWindow(d, dim_in);
    q_class = cellfun(@class, m.quad.q, 'UniformOutput', false);
    q_true = {'double', 'sdpvar', 'sdpvar', 'double', 'sdpvar', 'double'};
    testCase.verifyEqual(q_class, q_true)    
end

function testFullWindow(testCase)
    horizon_period = [0, 1];
    window = [0, 1];
    override = true;
    d = DisturbanceConstantWindow('d', {[]}, window, horizon_period, override);
    dim_lft_in = ones(1, sum(horizon_period));
    m = MultiplierConstantWindow(d, dim_lft_in);
    testCase.verifyClass(m.quad.q{1}, 'sdpvar')
    
    horizon_period = [0, 2];
    window = [0, 1, 2];
    override = true;
    d = DisturbanceConstantWindow('d', {[]}, window, horizon_period, override);
    dim_lft_in = ones(1, sum(horizon_period));
    m = MultiplierConstantWindow(d, dim_lft_in);
    testCase.verifyClass(m.quad.q{1}, 'sdpvar')
    testCase.verifyClass(m.quad.q{2}, 'sdpvar')
end

end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)