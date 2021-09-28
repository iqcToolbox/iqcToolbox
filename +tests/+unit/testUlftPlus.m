%% Requirements:
% 1. Ulft.plus shall perform the addition operation between two LFTs, whose output LFT is as expressed in Section 2.4 of "A Review of LFTs, LMIs and Mu" (1991).
%
% 2. Ulft.plus shall throw an error if the two LFTs are not conformable. A pair of LFTs conformable to addition must:
%     - have the same size
%     - have the same specified disturbances
%     - have the same specified performances
%
% 3. Ulft.plus shall output an LFT that does not have duplicates of the same disturbance or performance.
%
% 4. Ulft.plus shall be capable of taking one input argument that is not a Ulft object, as long as the input argument can be converted to a Ulft. If the non-Ulft input is not convertible, an error shall be thrown. Objects convertible to Ulfts are:
%     - doubles
%     - Delta objects
%     - ss objects
%
% 5. If the Ulft.plus operands have different horizon_periods, Ulft.vertcat shall ensure that the output Ulft shall have a resulting horizon_period that is consistent with the horizon_periods of both operands.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for Ulft.plus
classdef testUlftPlus < matlab.unittest.TestCase
methods (Test)

function testPlus(testCase)
    % Iterate over randomized situations
    n = 10;
    for i = 1 : n
        % Some properties will need to be unified across the random LFTs
        dim_in = randi([1,10]);
        dim_out = randi([1,10]);
        if i < n/2
            req_deltas = {'DeltaDelayZ'};
        else
            req_deltas = {'DeltaIntegrator'};
        end
        % Create the random LFTs ensuring that they never have same-named deltas
        % (otherwise the gathering operation will yield discrepancy with Doyle)
        lfts = cell(1,2);
        for j = 1 : length(lfts)
            lfts{j} = Ulft.random('dim_in', dim_in,...
                                  'dim_out', dim_out,...
                                  'req_deltas', req_deltas);
        end
        if mod(i, 2)
            lfts{1} = removeUncertainty(lfts{1}, 1);
        end
        lfts{2} = removeUncertainty(lfts{2}, 1);
        % Verify that the sum's SS matrices are correct (implicitly verifies combined horizon_period)
        lft_sum = lfts{1} + lfts{2};
        lfts{1} = matchHorizonPeriod(lfts{1}, lft_sum.horizon_period);
        lfts{2} = matchHorizonPeriod(lfts{2}, lft_sum.horizon_period);
        for t = 1 : sum(lft_sum.horizon_period)
            verifyEqual(testCase, lft_sum.a{t}, blkdiag(lfts{1}.a{t}, lfts{2}.a{t}));
            verifyEqual(testCase, lft_sum.b{t}, vertcat(lfts{1}.b{t}, lfts{2}.b{t}));
            verifyEqual(testCase, lft_sum.c{t}, horzcat(lfts{1}.c{t}, lfts{2}.c{t}));
            verifyEqual(testCase, lft_sum.d{t}, lfts{1}.d{t} + lfts{2}.d{t});
        end
        % Verify that the sum's deltas are correct
        verifyEqual(testCase, lft_sum.delta.names, [lfts{1}.delta.names lfts{2}.delta.names]);
    end
end

function testPlusWithGather(testCase)
    % Iterate over randomized situations
    n = 10;
    for i = 1 : n
        % Create conformable random LFTs that share identical deltas, disturbances, and performances
        dim_in = randi([1,10]);
        dim_out = randi([1,10]);
        req_deltas = Ulft.random('dim_in', dim_in, 'dim_out', dim_out).delta.deltas;
        dists = {DisturbanceL2('test_dist')};
        perfs = {PerformanceL2Induced('test_perf')};
        lfts = cell(1,2);
        for j = 1 : length(lfts)
            lfts{j} = Ulft.random('dim_in', dim_in,...
                                  'dim_out', dim_out,...
                                  'num_deltas', length(req_deltas),...
                                  'req_deltas', req_deltas);
            lfts{j} = lfts{j}.addDisturbance(dists).addPerformance(perfs);
        end
        % Make sure the sum properly combines the deltas, disturbances, and performances
        lft_sum = lfts{1} + lfts{2};
        verifyEqual(testCase, lft_sum.delta.names, lfts{1}.delta.names);
        verifyEqual(testCase, lft_sum.disturbance.names, lfts{1}.disturbance.names);
        verifyEqual(testCase, lft_sum.performance.names, lfts{1}.performance.names);
    end
end

function testPlusWithGatherABCD(testCase)
    % Iterate over randomized situations
    n = 10;
    for i = 1 : n
        % Some properties will need to be unified across the random LFTs
        dim_in = randi([1,10]);
        dim_out = randi([1,10]);
        if i < n/2
            req_deltas = {'DeltaDelayZ'};
        else
            req_deltas = {'DeltaIntegrator'};
        end
        % Create the random LFTs ensuring there is only a duplicate DeltaDelayZ/Integrator
        lft1 = Ulft.random('dim_in', dim_in,...
                           'dim_out', dim_out,...
                           'req_deltas', req_deltas,...
                           'num_deltas', 1);
        lft2 = Ulft.random('dim_in', dim_in,...
                           'dim_out', dim_out,...
                           'req_deltas', {lft1.delta.deltas{1}});
        % Verify that the sum's SS matrices are correct (implicitly verifies combined horizon_period)
        lft_sum = lft1 + lft2;
        lft1 = matchHorizonPeriod(lft1, lft_sum.horizon_period);
        lft2 = matchHorizonPeriod(lft2, lft_sum.horizon_period);
        for t = 1 : sum(lft_sum.horizon_period)
            verifyEqual(testCase, lft_sum.a{t}, blkdiag(lft1.a{t}, lft2.a{t}));
            verifyEqual(testCase, lft_sum.b{t}, vertcat(lft1.b{t}, lft2.b{t}));
            verifyEqual(testCase, lft_sum.c{t}, horzcat(lft1.c{t}, lft2.c{t}));
            verifyEqual(testCase, lft_sum.d{t}, lft1.d{t} + lft2.d{t});
        end
        % Verify that the sum's deltas are correct
        verifyEqual(testCase, lft_sum.delta.names, lft2.delta.names);
    end
end

function testPlusNonLft(testCase)
    % Iterate over randomized situations
    n = 10;
    for i = 1 : n
        % Nonlft on the left
        verifyClass(testCase,...
                    1 + Ulft.random('dim_in', 1, 'dim_out', 1),...
                    'Ulft');
        verifyClass(testCase,...
                    DeltaSlti('f') + Ulft.random('dim_in', 1, 'dim_out', 1),...
                    'Ulft');
        verifyClass(testCase,...
                    rss(3) + Ulft.random('dim_in', 1, 'dim_out', 1, 'req_deltas', {'DeltaIntegrator'}),...
                    'Ulft');
        % Nonlft on the right
        verifyClass(testCase,...
                    Ulft.random('dim_in', 1, 'dim_out', 1) + 1,...
                    'Ulft');
        verifyClass(testCase,...
                    Ulft.random('dim_in', 1, 'dim_out', 1) + DeltaSlti('f'),...
                    'Ulft');
        verifyClass(testCase,...
                    Ulft.random('dim_in', 1, 'dim_out', 1, 'req_deltas', {'DeltaIntegrator'}) + rss(3),...
                    'Ulft');
    end
end

function testPlusWithErrors(testCase)
    % Try to add LFTs with different size
    lft1 = Ulft.random('dim_in', 1, 'dim_out', 1);
    lft2 = Ulft.random('dim_in', 2, 'dim_out', 2);
    verifyError(testCase,...
                @() lft1 + lft2,...
                ?MException);
    % Try to add LFTs with different disturbances or performances
    lft1 = Ulft.random('dim_in', 1, 'dim_out', 1);
    lft1 = Ulft(lft1.a,...
                lft1.b,...
                lft1.c,...
                lft1.d,...
                lft1.delta,...
                'horizon_period', lft1.horizon_period,...
                'disturbance', DisturbanceL2('test_dist'));
    lft2 = Ulft.random('dim_in', 1, 'dim_out', 1);
    lft2 = Ulft(lft2.a,...
                lft2.b,...
                lft2.c,...
                lft2.d,...
                lft2.delta,...
                'horizon_period', lft2.horizon_period,...
                'performance', PerformanceL2Induced('test_perf'));
    verifyError(testCase,...
                @() lft1 + lft2,...
                ?MException);
    lft3 = Ulft.random('dim_in', 1, 'dim_out', 1);
    verifyError(testCase,...
                @() lft1 + lft3,...
                ?MException);
    verifyError(testCase,...
                @() lft2 + lft3,...
                ?MException);
end

end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0): Removed tests for implicitly converting cells of doubles - Micah Fry (micah.fry@ll.mit.edu)
% Aug. 26, 2021 (v.0.5.0): Initial release - Jason Nezvadovitz (jason.nezvadovitz@ll.mit.edu)