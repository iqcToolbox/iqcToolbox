%% Requirements:
% 1. toLft shall be capable of taking Ulft-convertible objects as single inputs and outputting a Ulft object to represent the equivalent lft form of the input argument.
%    Ulft-convertible objects are:
%   - 2D double arrays (nonempty)
%   - 1 x 1 cell array of doubles
%   - ss objects (both continuous- and discrete-time)
%   - Delta objects
%
% 2. toLft shall be capable of taking a, b, c, and d matrices and outputting a continuous-time Ulft object equivalent to the state-space system represented by those matrices.
%    toLft shall also take a fifth argument representing the timestep of a discrete-time system and output a discrete-time Ulft object.
%
% 3. toLft shall be capable of taking uncertain a, b, c, and d matrices (represented with Ulft objects), and outputting a continuous-time Ulft object equivalent to the uncertain state-space system represented by those uncertain matrices.
%    toLft shall also take a fifth argument representing the timestep of a discrete-time system and output a discrete-time Ulft object.
%
% 4. Each uncertain matrix provided to toLft may only have memoryless uncertainties, otherwise toLft shall throw an error.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for toLft
classdef testToLft < matlab.unittest.TestCase
methods (Test)

function testK(testCase)
    dim_in = randi([1,10]);
    dim_out = randi([1,10]);
    K = rand(dim_out, dim_in);
    lft = toLft(K);
    verifyEmpty(testCase, lft.a{1});
    verifyEmpty(testCase, lft.b{1});
    verifyEmpty(testCase, lft.c{1});
    verifyEqual(testCase, lft.d{1}, K);
    verifyEmpty(testCase, lft.delta.deltas);
end

function testKCell(testCase)
    horizon_period = [randi([0, 10]), randi([1, 10])];
    total_time = sum(horizon_period);
    dim_in = num2cell(randi([1,10], 1, total_time));
    dim_out = num2cell(randi([1,10], 1, total_time));
    K = cellfun(@(out, in) randn(out, in), dim_out, dim_in,...
                'UniformOutput', false);
    lft = toLft(K, horizon_period);
    verifyEmpty(testCase, lft.a{1});
    verifyEmpty(testCase, lft.b{1});
    verifyEmpty(testCase, lft.c{1});
    verifyEqual(testCase, lft.d, K);
    verifyEmpty(testCase, lft.delta.deltas);
    verifyEqual(testCase, lft.horizon_period, horizon_period);    
end

function testKCellBad(testCase)
    horizon_period = [0, 1];
    K = {1, 1};
    verifyError(testCase, @() toLft(K, horizon_period), ?MException);
end

function testErrorBadK(testCase)
    verifyError(testCase, @() toLft([]), ?MException);
    verifyError(testCase, @() toLft(rand(3, 3, 3)), ?MException);
    verifyError(testCase, @() toLft(rand(4, 4, 4, 4)), ?MException);
    verifyError(testCase, @() toLft(rand(5, 5, 5, 5, 5)), ?MException);
    verifyError(testCase, @() toLft({}), ?MException);
    verifyError(testCase, @() toLft({1, 1}), ?MException);
    verifyError(testCase, @() toLft({rand(2,2), rand(2,2)}), ?MException);
end

function testContinuousStateSpace(testCase)
    ss = rss(randi([1,10]));
    lft = toLft(ss);
    verifyEqual(testCase, lft.delta.deltas{1}, DeltaIntegrator(length(ss.a)));
    verifyEqual(testCase, lft.a{1}, ss.a);
    verifyEqual(testCase, lft.b{1}, ss.b);
    verifyEqual(testCase, lft.c{1}, ss.c);
    verifyEqual(testCase, lft.d{1}, ss.d);
end

function testDiscreteStateSpace(testCase)
    ss = drss(randi([1,10])); ss.Ts = rand();
    lft = toLft(ss);
    verifyEqual(testCase, lft.delta.deltas{1}, DeltaDelayZ(length(ss.A), ss.Ts));
    verifyEqual(testCase, lft.a{1}, ss.a);
    verifyEqual(testCase, lft.b{1}, ss.b);
    verifyEqual(testCase, lft.c{1}, ss.c);
    verifyEqual(testCase, lft.d{1}, ss.d);
end

function testDelta(testCase)
    del = Ulft.random('num_deltas', 2).delta.deltas{end};
    lft = toLft(del);
    verifyEqual(testCase, lft.a{1}, zeros(del.dim_in(1), del.dim_out(1)));
    verifyEqual(testCase, lft.b{1}, eye(del.dim_in(1)));
    verifyEqual(testCase, lft.c{1}, eye(del.dim_out(1)));
    verifyEqual(testCase, lft.d{1}, zeros(del.dim_out(1), del.dim_in(1)));
    verifyEqual(testCase, lft.delta.deltas{1}, del);
end

function testContinuousScalarABCD(testCase)
    a = rand(); b = rand(); c = rand(); d = rand();
    lft = toLft(a, b, c, d);
    verifyEqual(testCase, lft.delta.deltas{1}, DeltaIntegrator(length(a)));
    verifyEqual(testCase, lft.a{1}, a);
    verifyEqual(testCase, lft.b{1}, b);
    verifyEqual(testCase, lft.c{1}, c);
    verifyEqual(testCase, lft.d{1}, d);
end

function testContinuousABCD(testCase)
    arb = Ulft.random('num_deltas', 1, 'req_deltas', {'DeltaIntegrator'});
    lft = toLft(arb.a{1}, arb.b{1}, arb.c{1}, arb.d{1});
    verifyEqual(testCase, lft.delta.deltas{1}, DeltaIntegrator(length(arb.a{1})));
    verifyEqual(testCase, lft.a{1}, arb.a{1});
    verifyEqual(testCase, lft.b{1}, arb.b{1});
    verifyEqual(testCase, lft.c{1}, arb.c{1});
    verifyEqual(testCase, lft.d{1}, arb.d{1});
end

function testDiscreteTimeInvariantABCD(testCase)
    arb = Ulft.random('num_deltas', 1,...
                      'req_deltas', {'DeltaDelayZ'},...
                      'horizon_period', [0, 1]);
    lft = toLft(arb.a{1}, arb.b{1}, arb.c{1}, arb.d{1}, arb.delta.deltas{1}.timestep(1));
    verifyEqual(testCase, lft.delta, arb.delta);
    verifyEqual(testCase, lft.a{1}, arb.a{1});
    verifyEqual(testCase, lft.b{1}, arb.b{1});
    verifyEqual(testCase, lft.c{1}, arb.c{1});
    verifyEqual(testCase, lft.d{1}, arb.d{1});
end

function testUncertainContinuousABCD(testCase)
    a = 1 + DeltaSltv('del_a');
    b = 1 + DeltaSlti('del_b');
    c = 1 + DeltaSltvRateBnd('del_c');
    d = 1 + DeltaSltvRepeated({'del_d'});
    lft = toLft(a, b, c, d);
    verifyEqual(testCase, lft.delta.deltas{1}, DeltaIntegrator(1));
    lft = lft.removeUncertainty({'del_a', 'del_b', 'del_c', 'del_d'});
    verifyEqual(testCase, lft.a{1}, 1);
    verifyEqual(testCase, lft.b{1}, 1);
    verifyEqual(testCase, lft.c{1}, 1);
    verifyEqual(testCase, lft.d{1}, 1);
end

function testUncertainDiscreteABCD(testCase)
    a = 1 + DeltaSltv('del_a', 1, -1, 1, [randi([1, 10]), randi([2, 10])]);%??? work around to make this pass is use hp [0,1]
    b = 1 + DeltaSlti('del_b');
    c = 1 + DeltaSltvRateBnd('del_c');
    d = 1 + DeltaSltvRepeated({'del_d'});
    T = rand();
    lft = toLft(a, b, c, d, T);
    verifyEqual(testCase,...
                lft.delta.deltas{1},...
                DeltaDelayZ(1, T, a.horizon_period));
    lft = lft.removeUncertainty({'del_a', 'del_b', 'del_c', 'del_d'});
    verifyEqual(testCase, lft.a{1}, 1);
    verifyEqual(testCase, lft.b{1}, 1);
    verifyEqual(testCase, lft.c{1}, 1);
    verifyEqual(testCase, lft.d{1}, 1);
end

function testErrorABCD(testCase)
    dim_in = randi([1,10]);
    dim_out = randi([1,10]);
    dim_state = randi([1,10]);
    a = Ulft.random('dim_in', dim_state, 'dim_out', dim_state, 'req_deltas', {'DeltaIntegrator'});
    b = toLft(rand(dim_state, dim_in));
    c = toLft(rand(dim_out, dim_state));
    d = toLft(rand(dim_out, dim_in));
    verifyError(testCase, @() toLft(a, b, c, d), ?MException);
    verifyError(testCase, @() toLft(zeros(dim_state, dim_state-1), b, c, d), ?MException);
end

function testNotCorrectNumberOfArguments(testCase)
    verifyError(testCase, @() toLft(1, 2, 3), ?MException)
end

end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0): Changed testKCell* and added incorrect nargin test - Micah Fry (micah.fry@ll.mit.edu)
% Aug. 26, 2021 (v.0.5.0): Initial release - Jason Nezvadovitz (jason.nezvadovitz@ll.mit.edu)