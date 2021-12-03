%% Requirements:
%  1. Ulft.divide shall divide two LFTs by multiplying the numerator with the 
%     inverse of the denominator.
%  2. If the user provides a denominator which is non-invertible, Ulft shall
%     throw an error.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for Ulft.
classdef testUlftDivide < matlab.unittest.TestCase
    
methods (TestMethodSetup)
    function seedAndReportRng(testCase)
        seed = floor(posixtime(datetime('now')));
        rng('default');
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

function testDivideCorrectness(testCase)
    % Conformable dimensions
    ss_left = drss;
    ss_left.a = ss_left.a * 0.9;
    ss_right = drss(randi([1, 10]), size(ss_left, 2), size(ss_left, 2));
    ss_right.a = ss_right.a * 0.9;
    ss_right.d = ss_right.d + eye(size(ss_left, 2));
    ss_out = ss_left / ss_right;
    lft_left = toLft(ss_left);
    lft_right = toLft(ss_right);
    
    lft_out = lft_left / lft_right;
    % Check correctness of /
    sys_diff = norm(ss_out - lftToSs(lft_out), 'inf');
    testCase.verifyLessThan(sys_diff, 1e-6);
    
    % Make sure warning is thrown when doing typical mtimes
    warning_state = warning;
    warning('on', 'Ulft:rdivide')
    testCase.verifyWarning(@()  lft_left ./ lft_right, 'Ulft:rdivide');
    warning(warning_state)
    lft_out = lft_left ./ lft_right;
    % Check correctness of ./
    sys_diff = norm(ss_out - lftToSs(lft_out), 'inf');
    testCase.verifyLessThan(sys_diff, 1e-6);
        
    
    % When lft_right is scalar (need to compare via Linf norm)
    ss_left = drss(randi([1, 10]), randi([2, 10]), randi([2, 10]));
    ss_left.a = ss_left.a * 0.9;
    scalar_right = drss(1, 1, 1);
    scalar_right.a = scalar_right.a * 0.9;
    scalar_right.d = scalar_right.d + 1;
    ss_out = ss_left / scalar_right;
    lft_left = toLft(ss_left);
    lft_right = toLft(scalar_right);
    lft_out = lft_left ./ lft_right;
    sys_diff = norm(ss_out - lftToSs(lft_out), 'inf');
    testCase.verifyLessThan(sys_diff, 1e-6);
end

function testTimesNonLftRandom(testCase)
    % Define lfts
    
    % - LTI
    dim_outin = 2;
    lfti_l = Ulft.random('req_deltas', {DeltaIntegrator()},...
                         'dim_out', dim_outin,...
                         'dim_in', dim_outin);
    lfti_r = toLft(2 * eye(dim_outin));
    
    lfti_l_size = size(lfti_l, 2);
    % Convertible non-lft objects
    lft_double = eye(lfti_l_size);
    lft_delta = DeltaSltv('a', dim_outin);
    lft_delta_sqr = (DeltaSlti('a', dim_outin).^2 + 3 * eye(dim_outin))^2;
    lft_ss_cont = ss(eye(dim_outin), eye(dim_outin), eye(dim_outin), eye(dim_outin));
    lft_ss_disc = ss(eye(dim_outin), eye(dim_outin), eye(dim_outin), eye(dim_outin), -1);
    
    obj = 'non-convertible';
    
    % non-lft object is right
    lft_a1 = lfti_l ./ lft_double;
    lft_a3 = lfti_l / lft_delta_sqr;
    lft_a4 = lfti_l ./ lft_ss_cont;
    lft_a5 = lfti_l / lft_ss_cont;
    
    verifyClass(testCase, lft_a1, 'Ulft') 
    verifyClass(testCase, lft_a3, 'Ulft') 
    verifyClass(testCase, lft_a4, 'Ulft') 
    verifyClass(testCase, lft_a5, 'Ulft') 
    
    % non-lft object is left
    lft_b1 = lft_double / lfti_r;
    lft_b2 = lft_delta / lfti_r;
    lft_b3 = lft_delta_sqr ./ lfti_r;
    lft_b4 = lft_ss_disc / lfti_r;
    lft_b5 = lft_ss_cont ./ lfti_r;
    lft_b6 = lft_delta ./ lfti_r;
    
    verifyClass(testCase, lft_b1, 'Ulft') 
    verifyClass(testCase, lft_b2, 'Ulft') 
    verifyClass(testCase, lft_b3, 'Ulft') 
    verifyClass(testCase, lft_b4, 'Ulft') 
    verifyClass(testCase, lft_b5, 'Ulft') 
    verifyClass(testCase, lft_b6, 'Ulft') 
    
    % non-convertible non-lft object - right
    verifyError(testCase, @() rdivide(lfti_l, obj), ?MException)
    % non-convertible non-lft object - left
    verifyError(testCase, @() mrdivide(obj, lfti_r), ?MException)    
end

function testTimesCommonHp(testCase)
    z = DeltaDelayZ();
    a = {0};
    b = {1};
    c = {1};
    d = {2};
    a1 = {1, 1, 1, 1};
    b1 = {1, 1, 1, 1};
    c1 = {1, 1, 1, 1};
    d1 = {3, 1, 2, 5};
    lft1 = Ulft(a, b, c, d, z, 'horizon_period', [0 1]);
    lft2 = Ulft(a1, b1, c1, d1, z, 'horizon_period', [2 2]);
    correct_hp = [2 2];
    lftm = mrdivide(lft1, lft2);
    verifyEqual(testCase, lftm.horizon_period, correct_hp)
    
    lftr = rdivide(lft1, lft2);
    verifyEqual(testCase, lftr.horizon_period, correct_hp)
end
end
end

%%  CHANGELOG
% Dec. 02, 2021: Added after v0.6.0 - Micah Fry (micah.fry@ll.mit.edu)