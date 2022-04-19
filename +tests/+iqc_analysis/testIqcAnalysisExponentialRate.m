%% Requirements:
%  1. IQC analysis shall produce an "infeasible problem" result when analyzing
%     systems if their exponential decay rate is greater than that specified
%     by the analysis options.
%  2. IQC analysis shall produce a valid solution proving that an uncertain
%     system has a pre-specified exponential decay rate for many uncertain
%     systems which have such a decay rate. Producing a certificate for ALL
%     uncertain systems with such a decay rate is not expected.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for IQC analysis identifying Exponential Rate of Convergence
classdef testIqcAnalysisExponentialRate < matlab.unittest.TestCase

methods (TestMethodSetup)
function seedAndReportRng(testCase)
    seed = floor(posixtime(datetime('now')));
    rng(seed, 'twister');
    diagnose_str = ...
        sprintf(['Random inputs may be regenerated by calling: \n',...
                 '>> rng(%10d) \n',...
                 'before running the remainder of the test''s body'],...
                seed);
    testCase.onFailure(@() fprintf(diagnose_str));
end    
end

methods (Test)
function testNominalSystems(testCase)
    % Discrete-time
    g = drss;
    g.a = g.a * 0.95;
    exponential = max(abs(eig(g.a)));
    g_lft = toLft(g);
    g_lft = g_lft.addPerformance({PerformanceStable()});
    options = AnalysisOptions('verbose', false);
    % Check if decay rate may be slightly slower than known rate
    options.exponential = exponential * 1.01; 
    result = iqcAnalysis(g_lft, 'analysis_options', options);
    testCase.verifyTrue(result.valid)
    % Check if decay rate may be slightly faster than known rate
    options.exponential = exponential * 0.99; 
    result = iqcAnalysis(g_lft, 'analysis_options', options);
    testCase.verifyFalse(result.valid)

    % Continuous-time
    g = rss;
    g.a = g.a - 0.05 * eye(size(g.a));
    exponential = -max(real(eig(g.a)));
    g_lft = toLft(g);
    g_lft = g_lft.addPerformance({PerformanceStable()});
    options = AnalysisOptions('verbose', false);
    % Check if decay rate may be slightly slower than known rate
    options.exponential = exponential * 0.99;
    result = iqcAnalysis(g_lft, 'analysis_options', options);
    testCase.verifyTrue(result.valid)
    % Check if decay rate may be slightly faster than known rate
    options.exponential = exponential * 1.1;
    result = iqcAnalysis(g_lft, 'analysis_options', options);
    testCase.verifyFalse(result.valid)
end
   

function testUncertainSystemsExpoIndependent(testCase)
%     rng(1649681707) % Makes discrete system fail when P indefinite (even though G and Psi have small enough eigs)
                    % Also makes continuous system fail when P indefinite (even though G and Psi have "left enough" eigs)
  % These tests fail with some Slti, unless I make P > 0 (rather than indefinite)
  rng(1650397751)
    deltas = {'DeltaSlti','DeltaSltv','DeltaSltvRateBnd','DeltaSectorBounded'};
    del_ind = randi([1, length(deltas)]);
    del_type = deltas{del_ind};
%     del_type = 'DeltaSlti'; % Uncomment this is topline rng is uncommented

    % Discrete-time
    valid_exponential = false;
    while ~valid_exponential
        g = drss(randi([1,5]));
        g.a = g.a * 0.95;
        exponential = max(abs(eig(g.a)));
        valid_exponential = exponential > 0.3;
    end
    margin = (1 - exponential) / exponential;
    del_bnd = 1 + margin / 2; % Guaranteed to be greater than 1, will not destabilize
    switch del_type
        case 'DeltaSlti'
            del = DeltaSlti('del', size(g.a, 1), -del_bnd, del_bnd);
        case 'DeltaSltv'
            del = DeltaSltv('del', size(g.a, 1), -del_bnd, del_bnd);
        case 'DeltaSltvRateBnd'
            del = DeltaSltvRateBnd('del', size(g.a, 1), -del_bnd, del_bnd);
        case 'DeltaSectorBounded'
            del = DeltaSectorBounded('del', size(g.a, 1), -del_bnd, del_bnd);
    end
    a_del = g.a * del;
    g_lft_del = toLft(a_del, g.b, g.c, g.d, -1);
    g_lft_del = g_lft_del.addPerformance({PerformanceStable()});
    options = AnalysisOptions('verbose', false, 'lmi_shift', 1e-6);
    options.exponential = exponential * (1 + 2 * margin / 3);
    testCase.assertGreaterThan(options.exponential, del_bnd * exponential)
    result = iqcAnalysis(g_lft_del, 'analysis_options', options);
    testCase.verifyTrue(result.valid)
    options.exponential = exponential * (1 + margin / 3);
    testCase.assertLessThan(options.exponential, del_bnd * exponential)
    result = iqcAnalysis(g_lft_del, 'analysis_options', options);
    testCase.verifyFalse(result.valid)
    
    % Continuous-time
    valid_exponential = false;
    while ~valid_exponential
        g = rss;
        g.a = g.a - 0.05 * eye(size(g.a));
        exponential = -max(real(eig(g.a)));
        valid_exponential = exponential < 2.5;
    end
    del_bnd = exponential * 0.9;
    switch del_type
        case 'DeltaSlti'
            del = DeltaSlti('del', size(g.a, 1), -del_bnd, del_bnd);
        case 'DeltaSltv'
            del = DeltaSltv('del', size(g.a, 1), -del_bnd, del_bnd);
        case 'DeltaSltvRateBnd'
            del = DeltaSltvRateBnd('del', size(g.a, 1), -del_bnd, del_bnd);
        case 'DeltaSectorBounded'
            del = DeltaSectorBounded('del', size(g.a, 1), -del_bnd, del_bnd);        
    end
    a_del = g.a + del;
    g_lft_del = toLft(a_del, g.b, g.c, g.d);
    g_lft_del = g_lft_del.addPerformance({PerformanceStable()});
    options = AnalysisOptions('verbose', false, 'lmi_shift', 1e-6);
    options.exponential = exponential - 1.1 * del_bnd;
    testCase.assertLessThan(options.exponential, exponential - del_bnd)
    result = iqcAnalysis(g_lft_del, 'analysis_options', options);
    testCase.verifyTrue(result.valid)
    options.exponential = exponential - 0.9 * del_bnd;
    testCase.assertLessThan(exponential - del_bnd, options.exponential)
    result = iqcAnalysis(g_lft_del, 'analysis_options', options);
    testCase.verifyFalse(result.valid)
end

function testDiscreteTimeUncertainSystems(testCase) 
% Drawn from Bozcar, Ross. "Performance Guarantees in Learning and Robust Control," (2019)
    z = tf('z');
    num = -(z + 1) * (10 * z + 9);
    den = (2 * z - 1) * (5 * z - 1) * (10 * z - 1);
    g = toLft(ss(num / den));
    dim_outin = 1;
    bnd = 0.75;
    del_bnd = DeltaBounded('bnd', dim_outin, dim_outin, bnd);
    g_del = interconnect(toLft(del_bnd), g);
    g_del = g_del.addPerformance({PerformanceStable()});
    options = AnalysisOptions('verbose', false);
    options.exponential = 0.95;
    result = iqcAnalysis(g_del, 'analysis_options', options);
    testCase.verifyTrue(result.valid)
    options.exponential = 0.93;
    result = iqcAnalysis(g_del, 'analysis_options', options);
    testCase.verifyFalse(result.valid)
    
    bnd = 2;
    del_sb = DeltaSectorBounded('sb', dim_outin, 0, bnd);
    g_del = interconnect(toLft(del_sb), g);
    g_del = g_del.addPerformance({PerformanceStable()});
    options.exponential = 0.95;
    result = iqcAnalysis(g_del, 'analysis_options', options);
    testCase.verifyTrue(result.valid)
    options.exponential = 0.9;
    result = iqcAnalysis(g_del, 'analysis_options', options);
    testCase.verifyFalse(result.valid)
end

function testConstantDelayContinuousTime(testCase)
    % Check that robust stability is correctly concluded (no exponential rate bound)
    wn = 10; % natural freq rad/s
    zeta = 0.5; % damping coefficient
    s = tf('s');
    g = wn^2 / (s^2 + 2 * zeta * wn * s + wn^2);
    g = ss(g);
    delay_max = 0.1;
    del = DeltaConstantDelay2('del', 1, delay_max); % 0.1 second delay
    g_del = (1 + del) * g;
    g_del_cl = interconnect(-1, [1; 1] * g_del * [1, 1]);
    options = AnalysisOptions('verbose', false, 'lmi_shift', 1e-7);
    result = iqcAnalysis(g_del_cl, 'analysis_options', options);
    testCase.verifyTrue(result.valid); 
    % Include exponential rate bound
    expo = 0.1; % Should extend out to 0.7
    options.exponential = expo;
    m(2) = MultiplierConstantDelay2(del,...
                                    'discrete', false,...
                                    'basis_poles', -expo * 1.1,...
                                    'exponential', expo);
    result = iqcAnalysis(g_del_cl,...
                         'analysis_options', options,...
                         'multipliers_delta', m);
    testCase.verifyTrue(result.valid); 
    % Set exponential rate bound too fast, check that it correctly fails
    expo = 0.8; 
    options.exponential = expo;
    m(2) = MultiplierConstantDelay2(del,...
                                    'discrete', false,...
                                    'basis_poles', -expo * 1.1,...
                                    'exponential', expo);
    result = iqcAnalysis(g_del_cl,...
                         'analysis_options', options,...
                         'multipliers_delta', m);
    testCase.verifyFalse(result.valid); 
    
    % Check that robust stability correctly fails (with big enough delay, no exponential specification)
    delay_max = 0.2;
    del = DeltaConstantDelay2('del', 1, delay_max); 
    g_del = (1 + del) * g;
    g_del_cl = interconnect(-1, [1; 1] * g_del * [1, 1]);
    options.exponential = 0;
    result = iqcAnalysis(g_del_cl, 'analysis_options', options);
    testCase.verifyFalse(result.valid);
% These tests can be verified with checking the gain/phase margins of delayed plant:
% %     delay_max = 0.1;
% %     s = tf('s');
% %     delay = exp(-s * delay_max);
% %     margin(delay * g)
% %     expo = 0.8; % 0.7 shows a positive gain/phase margin
% %     g_expo = g;
% %     g_expo.a = g_expo.a + expo * eye(size(g_expo.a));
% %     delay_expo = exp(delay_max * expo) * exp(-s * delay_max);
% %     margin(delay_expo * g_expo)
end

function testConstantDelayDiscreteTime(testCase)
    % Check that robust stability is correctly concluded (no exponential rate bound)
    wn = 10; % natural freq rad/s
    zeta = 0.5; % damping coefficient
    s = tf('s');
    g = ss(wn^2 / (s^2 + 2 * zeta * wn * s + wn^2));
    delay_max_c = 0.1; % 0.1 sec delay
    dt = 0.01; % 100 Hz frequency
    gd = c2d(g, dt);
    delay_max = delay_max_c / dt;
    del = DeltaConstantDelay2('del', 1, delay_max); % 10 step delay
    g_del = (1 + del) * gd;
    g_del_cl = interconnect(-1, [1; 1] * g_del * [1, 1]);
    options = AnalysisOptions('verbose', false, 'lmi_shift', 1e-6);
    result = iqcAnalysis(g_del_cl, 'analysis_options', options);
    testCase.verifyTrue(result.valid); 
    
    % Check that robust stability correctly fails (with big enough delay, no exponential specification)
    delay_max = 20;
    del = DeltaConstantDelay2('del', 1, delay_max); 
    g_del = (1 + del) * gd;
    g_del_cl = interconnect(-1, [1; 1] * g_del * [1, 1]);
    result = iqcAnalysis(g_del_cl, 'analysis_options', options);
    testCase.verifyFalse(result.valid); 
    
    % Check that robustly stable system passes conforming exponential rate bound
    delay_max = 10;
    del = DeltaConstantDelay2('del', 1, delay_max); 
    g_del = (1 + del) * gd;
%     g_del = toLft(gd);
    g_del_cl = interconnect(-1, [1; 1] * g_del * [1, 1]);
    expo = 0.999; % Should extend down to 0.994
    options.exponential = expo;
    result = iqcAnalysis(g_del_cl, 'analysis_options', options);
    testCase.verifyTrue(result.valid); 
    % Set exponential rate bound too fast, check that it correctly fails
    expo = 0.993; 
    options.exponential = expo;
    result = iqcAnalysis(g_del_cl, 'analysis_options', options);
    testCase.verifyFalse(result.valid); 
% These tests can be verified with checking the gain/phase margins of delayed plant:
% %     delay_max = 10;
% %     z = tf('z');
% %     delay = 1/z^delay_max;
% %     margin(delay * gd)
% %     expo = .993; % 0.994 shows positive gain/phase margin
% %     g_expo = gd;
% %     g_expo.a = g_expo.a / expo;
% %     g_expo.b = g_expo.b / expo;
% %     delay_expo = 1/ expo^delay_max / z^delay_max;
% %     margin(delay_expo * g_expo)
% %     cl_expo = delay_expo * g_expo / (1 + delay_expo * g_expo);
% %     abs(eig(cl_expo.a))
end

function testExponentiallyStableMultipliers(testCase)
    del_slti = DeltaSlti('slti');
    del_sltvrb = DeltaSltvRateBnd('sltvrb');
%     del_delay = DeltaConstantDelay2('delay');
    
    % Check that multipliers with conforming convergence rate pass (continuous time)
    gc = toLft(ss(-100, 1, 1, 0)) * del_slti * del_sltvrb;% * (1 + del_delay);
    expo = 0.49;
    options = AnalysisOptions('verbose', false,...
                              'exponential', expo,...
                              'lmi_shift', 1e-5);
    result = iqcAnalysis(gc, 'analysis_options', options);
    m_del = MultiplierDeltaCombined(result.multipliers_delta);
    m_shifted = lftToSs(m_del.shiftMultiplier(expo).filter_lft);
    testCase.assertTrue(isstable(m_shifted))
    testCase.verifyTrue(result.valid)
    
    % Check that multipliers with slower convergence rate fail (continuous time)
    expo = 0.51;
    options.exponential = expo;
    testCase.verifyError(@() iqcAnalysis(gc, 'analysis_options', options),...
                         'iqcAnalysis:iqcAnalysis')
    
    % Check that multipliers with conforming convergence rate pass (discrete time)
    del_z = DeltaDelayZ();
    gd = del_z * del_slti * del_sltvrb;% * (1 + del_delay);
    expo = 0.51;
    options.exponential = expo;
    result = iqcAnalysis(gd, 'analysis_options', options);
    m_del = MultiplierDeltaCombined(result.multipliers_delta);
    m_shifted = lftToSs(m_del.shiftMultiplier(expo).filter_lft);
    testCase.assertTrue(isstable(m_shifted))
    testCase.verifyTrue(result.valid)
    
    % Check that multipliers with slower convergence rate fail (continuous time)
    expo = 0.49;
    options.exponential = expo;
    testCase.verifyError(@() iqcAnalysis(gd, 'analysis_options', options),...
                         'iqcAnalysis:iqcAnalysis')
end

function testDefaultExponentialRateForMemoryless(testCase)
    dim_outin = randi([1, 10]);
    delay_max = randi([1, 10]);
    d_delay = DeltaConstantDelay('delay', dim_outin, delay_max);
    g = randn(dim_outin, dim_outin);
    g = g * 0.99 / norm(g, 2);
    eye_mat = eye(dim_outin);
    g = [eye_mat; eye_mat] * g * [eye_mat, eye_mat];
%         if isempty(g.a)
%             g = g.d;
%         end
    lft_delay = interconnect(toLft(d_delay), g);
    m = MultiplierConstantDelay(d_delay, 'discrete', true);
    options = AnalysisOptions('lmi_shift', 1e-6, 'verbose', false);
    result = iqcAnalysis(lft_delay,...
                         'analysis_options', options,...
                         'multipliers_delta', m);
    testCase.verifyTrue(result.valid);
end
% function testConstantDelayDiscreteTime(testCase)
% Derived from Stabilization of Discrete-Time Systems with Input Delays - Bin Zhou
%     n = 4;
%     dt = 0.1;
%     a = zeros(n);
%     b = [];
%     for i = 1:n
%         term = dt ^ (i - 1) / factorial(i - 1);
%         a = a + diag(term * ones(n - (i - 1), 1), (i - 1));
%         b = [dt^i / factorial(i); b];
%     end
%     g_ol = toLft(a, b, eye(n), zeros(n, 1), -1);
%     g_aug = [eye(n); eye(n)] * g_ol * [1, 1];
%     gam = 0.025;
%     f = -[gam^4 / dt^4;
%           (8 - 3 * gam) * gam^3 / 2 / dt^3;
%           (11 * gam^2 - 48 * gam + 72) * gam^2 / 12 / dt^2;
%           (48 - 36 * gam - 3 * gam^3 + 16 * gam^2) * gam / 12 / dt]';
%     delay_max = 4;
%     del = DeltaConstantDelay('del', size(b, 2), delay_max);
%     g_cl = interconnect(del * f, g_aug);
    
    
%     wn = 10; % natural freq rad/s
%     zeta = 0.5; % damping coefficient
%     s = tf('s');
%     g = wn^2 / (s^2 + 2 * zeta * wn * s + wn^2);
%     delay_max = 0.1;
%     g = ss(g);
%     z = tf('z');
%     dt = 0.01;
%     gd = c2d(g, dt);
%     delay_max_d = delay_max / dt; % delay in discrete timesteps
%     del = DeltaConstantDelay('del', 1, delay_max); 
%     expo = .9;
%     options.exponential = expo;
%     delay = z^(-delay_max_d);
%     g_expo = gd;
%     g_expo.a = g_expo.a / expo; 
%     g_expo.b = g_expo.b / expo;
%     gd_tf = tf(gd);
%     expo_vec = [expo^2, expo, 1];
%     num = gd_tf.numerator; num{1} = num{1} .* expo_vec;
%     den = gd_tf.denominator; den{1} = den{1} .* expo_vec;
%     gd_tf_expo = tf(num, den, dt);
%     margin(delay * g_expo)
%     g_del = del * gd;
%     result = iqcAnalysis(g_del, 'analysis_options', options);
%     testCase.verifyTrue(result.valid);    
%     expo = 1.4;
%     options.exponential = expo;
%     delay = exp(-s * delay_max_d);
%     result = iqcAnalysis(g_del, 'analysis_options', options);
%     testCase.verifyFalse(result.valid);
    
% end
end
end

%%  CHANGELOG
% Apr. 18, 2022: Added after v0.9.0 - Micah Fry (micah.fry@ll.mit.edu)