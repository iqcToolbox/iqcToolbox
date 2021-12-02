%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for DeltaSectorBounded and MultiplierSectorBounded
classdef testDeltaSectorBounded < matlab.unittest.TestCase
    
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
function testFullConstructor(testCase)
    name = 'test';
    dim_outin = 2;
    lower_bound = 0;
    upper_bound = 2;
    horizon_period = [2, 3];
    d = DeltaSectorBounded(name,...
                           dim_outin,...
                           lower_bound,...
                           upper_bound,...
                           horizon_period)
    total_time = sum(horizon_period);
    testCase.verifyEqual(d.name, name)
    testCase.verifyEqual(d.dim_out, repmat(dim_outin, 1, total_time))
    testCase.verifyEqual(d.dim_in,  repmat(dim_outin, 1, total_time))
    testCase.verifyEqual(d.lower_bound, repmat(lower_bound, 1, total_time))
    testCase.verifyEqual(d.upper_bound, repmat(upper_bound, 1, total_time))
    testCase.verifyEqual(d.horizon_period, horizon_period)
end

function testFourArgConstructor(testCase)
    name = 'test3';
    dim_outin = 3;
    lower_bound = -3.1;
    upper_bound = 1.5;
    d = DeltaSectorBounded(name, dim_outin, lower_bound, upper_bound);
    testCase.verifyEqual(d.name, name)
    testCase.verifyEqual(d.dim_out, dim_outin)
    testCase.verifyEqual(d.dim_in,  dim_outin)
    testCase.verifyEqual(d.lower_bound, lower_bound)
    testCase.verifyEqual(d.upper_bound, upper_bound)
    testCase.verifyEqual(d.horizon_period, [0, 1])
end

function testTwoArgConstructor(testCase)
    name = 'test2';
    dim_outin = 3;
    d = DeltaSectorBounded(name, dim_outin);
    testCase.verifyEqual(d.name, name)
    testCase.verifyEqual(d.dim_out, dim_outin)
    testCase.verifyEqual(d.dim_in,  dim_outin)
    testCase.verifyEqual(d.lower_bound, -1)
    testCase.verifyEqual(d.upper_bound, 1)
    testCase.verifyEqual(d.horizon_period, [0, 1])
end

function testOneArgConstructor(testCase)
    name = 'test1';
    d = DeltaSectorBounded(name);
    testCase.verifyEqual(d.name, name)
    testCase.verifyEqual(d.dim_out, 1)
    testCase.verifyEqual(d.dim_in,  1)
    testCase.verifyEqual(d.lower_bound, -1)
    testCase.verifyEqual(d.upper_bound, 1)
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
    del = DeltaSectorBounded(name);
    m = MultiplierSectorBounded(del);
    testCase.verifyTrue(m.quad_time_varying)

    quad_time_varying = false;
    m = MultiplierSectorBounded(del, 'quad_time_varying', quad_time_varying);
    testCase.verifyEqual(m.quad_time_varying, quad_time_varying)
end
end
end

%%  CHANGELOG
% Nov. 23, 2021: Added after v0.6.0 - Micah Fry (micah.fry@ll.mit.edu)