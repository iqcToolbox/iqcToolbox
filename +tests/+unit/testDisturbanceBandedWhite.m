%% Test class for DisturbanceBandedWhite and MultiplierBandedWhite
classdef testDisturbanceBandedWhite < matlab.unittest.TestCase
    
methods (TestMethodSetup)
function seedAndReportRng(testCase)
    seed = floor(posixtime(datetime('now')));
    rng(seed);
    diagnose_str = ...
        sprintf(['Random inputs may be regenerated by calling: \n',...
                 '>> rng(%10d) \n',...
                 'before running the remainder of the test''s body'],...
                seed);
    testCase.onFailure(@() fprintf(diagnose_str));
end    
end
    
methods (Test)
function testDisturbanceFullConstructor(testCase)
    name = 'test';
    chan_in = {2};
    omega = pi/2;
    horizon_period = [3, 1];
    d = DisturbanceBandedWhite(name, chan_in, omega, horizon_period)
    testCase.verifyEqual(d.name, name)
    testCase.verifyEqual(d.chan_in, repmat(chan_in, 1, sum(horizon_period)))
    testCase.verifyEqual(d.omega, omega)
    testCase.verifyEqual(d.horizon_period, horizon_period)
end

function testThreeArgConstructor(testCase)
    name = 'test';
    chan_in = {[1]};
    omega = pi;
    d = DisturbanceBandedWhite(name, chan_in, omega);
    testCase.verifyEqual(d.name, name)
    testCase.verifyEqual(d.chan_in, chan_in)
    testCase.verifyEqual(d.omega, omega)
    testCase.verifyEqual(d.horizon_period, [0, 1])
end

function testTwoArgConstructor(testCase)
    name = 'test';
    chan_in = {4};
    d = DisturbanceBandedWhite(name, chan_in);
    testCase.verifyEqual(d.name, name)
    testCase.verifyEqual(d.chan_in, chan_in)
    testCase.verifyEqual(d.omega, pi)
    testCase.verifyEqual(d.horizon_period, [0, 1])
end

function testOneArgConstructor(testCase)
    name = 'test';
    d = DisturbanceBandedWhite(name);
    testCase.verifyEqual(d.name, name)
    testCase.verifyEqual(d.chan_in, {1})
    testCase.verifyEqual(d.omega, pi)
    testCase.verifyEqual(d.horizon_period, [0, 1])
end

function testBadConstructorCalls(testCase)
    testCase.verifyError(@() DisturbanceBandedWhite(),...
                         'DisturbanceBandedWhite:DisturbanceBandedWhite')
    bad_chan = {1, 1};
    testCase.verifyError(@() DisturbanceBandedWhite('test', bad_chan, 1, [0, 2]),...
                         'DisturbanceBandedWhite:DisturbanceBandedWhite')
    bad_chan = {[]};
    testCase.verifyError(@() DisturbanceBandedWhite('test', bad_chan),...
                         'DisturbanceBandedWhite:DisturbanceBandedWhite')
    bad_chan = {[1; 2]};
    testCase.verifyError(@() DisturbanceBandedWhite('test', bad_chan),...
                         'DisturbanceBandedWhite:DisturbanceBandedWhite')
    bad_omega = 0;
    testCase.verifyError(@() DisturbanceBandedWhite('test', {1}, bad_omega),...
                         ?MException)
end

function testMatchHorizonPeriod(testCase)
    name = 'test';
    chan_in = {2};
    omega = pi/2;
    d = DisturbanceBandedWhite(name, chan_in, omega);
    new_hp = [2, 5];
    d = d.matchHorizonPeriod(new_hp);
    testCase.verifyEqual(d.chan_in, repmat(chan_in, 1, sum(new_hp)))
    testCase.verifyEqual(d.omega, omega)
    testCase.verifyEqual(d.horizon_period, new_hp)
end

function testMultiplierConstruction(testCase)
    name = 'test';
    d = DisturbanceBandedWhite(name);
    dim_in = 1;
    discrete = true;
    m = MultiplierBandedWhite(d, dim_in, discrete)
end

% function testMixingNonperiodicAndPeriodicTimesteps(testCase)
%     % Should throw an error without override (because window includes periodic and non-periodic portions)
%     window = [1:2];
%     horizon_period = [2, 2];
%     bad_d = @()DisturbanceConstantWindow('d', {[]}, window, horizon_period);
%     testCase.verifyError(bad_d,...
%                          'DisturbanceConstantWindow:DisturbanceConstantWindow')
%     
%     % Now with override, it will make a window that has a periodic portion
%     override = true;
%     d = DisturbanceConstantWindow('d', {[]}, window, horizon_period, override);
%     dim_in = ones(1, sum(d.horizon_period));
%     m = MultiplierConstantWindow(d, dim_in);
%     q_class = cellfun(@class, m.quad.q, 'UniformOutput', false);
%     q_true = {'double', 'sdpvar', 'sdpvar', 'double'};
%     testCase.verifyEqual(q_class, q_true)
%     
%     % Periodicity can be seen by extending the horizon
%     extend_horizon = 2;
%     new_hp = d.horizon_period + [extend_horizon, 0];
%     d = d.matchHorizonPeriod(new_hp);
%     dim_in = ones(1, sum(d.horizon_period));
%     m = MultiplierConstantWindow(d, dim_in);
%     q_class = cellfun(@class, m.quad.q, 'UniformOutput', false);
%     q_true = {'double', 'sdpvar', 'sdpvar', 'double', 'sdpvar', 'double'};
%     testCase.verifyEqual(q_class, q_true)    
% end
% 
% function testMatchHorizonPeriodCorrectness(testCase)
%     horizon_period = [2, 5];
%     window = [1, 2, 4];
%     override = true;
%     d = DisturbanceConstantWindow('dis', {[]}, window, horizon_period,override);
%     dim_in = ones(1, sum(horizon_period));
%     mult = MultiplierConstantWindow(d, dim_in);
%     q_class = cellfun(@class, mult.quad.q, 'UniformOutput', false);
%     q_true = {'double', 'sdpvar',...  % Non-periodic portion
%               'sdpvar', 'double', 'sdpvar', 'double', 'double'}; % Periodic portion
%     % Establish correctness before changing horizon_period
%     testCase.assertEqual(q_class, q_true)
%     
%     new_hp = [6, 10];
%     d_new_hp = d.matchHorizonPeriod(new_hp);
%     dim_in = ones(1, sum(new_hp));
%     mult = MultiplierConstantWindow(d_new_hp, dim_in);
%     q_class = cellfun(@class, mult.quad.q, 'UniformOutput', false);
%     q_true = {'double', 'sdpvar', 'sdpvar', 'double', 'sdpvar', 'double',...  % Non-periodic portion
%               'double', 'sdpvar', 'double', 'sdpvar', 'double', 'double', 'sdpvar', 'double', 'sdpvar', 'double'}; % Periodic portion
%     testCase.verifyEqual(q_class, q_true)
% end
% 
% function testFullWindow(testCase)
%     horizon_period = [0, 1];
%     window = 1;
%     override = true;
%     d = DisturbanceConstantWindow('d', {[]}, window, horizon_period, override);
%     dim_lft_in = ones(1, sum(horizon_period));
%     m = MultiplierConstantWindow(d, dim_lft_in);
%     % Check correct structure of multiplier
%     testCase.verifyClass(m.quad.q{1}, 'sdpvar')
%         
%     horizon_period = [0, 2];
%     window = 1:2;
%     override = true;
%     d = DisturbanceConstantWindow('d', {[]}, window, horizon_period, override);
%     dim_lft_in = ones(1, sum(horizon_period));
%     m = MultiplierConstantWindow(d, dim_lft_in);
%     testCase.verifyClass(m.quad.q{1}, 'sdpvar')
%     testCase.verifyClass(m.quad.q{2}, 'sdpvar')
% end
% 
% function testBadConstructorCalls(testCase)
%     testCase.verifyError(@() DisturbanceConstantWindow(),...
%                          'DisturbanceConstantWindow:DisturbanceConstantWindow')
%     testCase.verifyError(@() DisturbanceConstantWindow('d', {[]}),...
%                          'DisturbanceConstantWindow:DisturbanceConstantWindow')
% end
% 
% function testDisplayDisturbance(testCase)
%     lft = toLft(eye(2));
%     horizon_period = [3, 5];    
%     lft = lft.matchHorizonPeriod(horizon_period);
%     d1 = DisturbanceConstantWindow('d1', {1}, 4:6, horizon_period);
%     d2 = DisturbanceConstantWindow('d2', {2}, 1:2, horizon_period);
%     % Add disturbances and display
%     lft = lft.addDisturbance({d1, d2})
% end
% 
% function testBadDeltaToMultiplierCall(testCase)
%     d = DisturbanceConstantWindow('d');
%     testCase.verifyError(@()d.disturbanceToMultiplier,...
%                          'DisturbanceConstantWindow:disturbanceToMultiplier')
% end
% 
% function testMultiplierConstructor(testCase)
%     name = 'd';
%     chan_in = {2};
%     window = [1:3];
%     horizon_period = [3, 1];
%     override = true;
%     d = DisturbanceConstantWindow(name,...
%                                   chan_in,...
%                                   window,...
%                                   horizon_period,...
%                                   override);              
%     dim_in_lft = 2 * ones(1, sum(horizon_period));
%     m = MultiplierConstantWindow(d, dim_in_lft);
%     testCase.verifyEqual(m.name, 'd')
%     testCase.verifyEqual(m.chan_in, repmat(chan_in, 1, sum(horizon_period)))
%     testCase.verifyEqual(m.window, window)
%     testCase.verifyEqual(m.horizon_period, horizon_period)    
%     testCase.verifyEqual(m.dim_in, dim_in_lft)
% end
% 
% function testMultiplierBadConstructionCall(testCase)
%     testCase.verifyError(@() MultiplierConstantWindow([], []),...
%                          ?MException)
% end
% 
% function testTimeVaryingQuad(testCase)
%     g = drss(3);
%     g.a = g.a * 0.9;
%     lft = toLft(g);
%     % Extend lft to be time-varying
%     horizon_period = [0, 20];
%     lft = lft.matchHorizonPeriod(horizon_period);
%     % Create a disturbance which is constant for all time (therefore, 0)
%     dis = DisturbanceConstantWindow('d',...
%                                     {[]},...
%                                     1:(horizon_period(2) - 1),...
%                                     horizon_period);
%     lft = lft.addDisturbance({dis});
%     options = AnalysisOptions('verbose', false, 'lmi_shift', 1e-6);
%     % The upper bound using time-varying quads
%     result_lower = iqcAnalysis(lft, 'analysis_options', options);
%     testCase.assertTrue(result_lower.valid)
%     quads = cellfun(@double, result_lower.multipliers_disturbance.quad.q,...
%                     'UniformOutput', false);
%     testCase.verifyEqual(quads{1}, 0)
%     testCase.verifyFalse(all(cellfun(@(q) isequal(q, quads{2}), quads(2:end))))
%     
%     dis = lft.disturbance.disturbances{1};
%     dim_in = size(lft, 2);
%     mult = MultiplierConstantWindow(dis, dim_in, 'quad_time_varying', false);
%     % Generate results forcing quad to be time-invariant (should be greater than result_lower)
%     result_qti = iqcAnalysis(lft,...
%                              'analysis_options', options,...
%                              'multipliers_disturbance', mult);
%     testCase.assertTrue(result_qti.valid)
%     quads_ti = cellfun(@double, result_qti.multipliers_disturbance.quad.q,...
%                     'UniformOutput', false);
%     testCase.verifyEqual(quads_ti{1}, 0)
%     testCase.verifyTrue(all(cellfun(@(q) isequal(q, quads_ti{2}), quads_ti(2:end))))
% end
end
end

%%  CHANGELOG
% Nov. 18, 2021: Added after v0.6.0 - Micah Fry (micah.fry@ll.mit.edu)