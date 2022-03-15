%% Requirements:
%  1. iqcAnalysis shall produce a valid upper-bound on worst-case performance
%     for any memory-less nominal system.
%  2. 
%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for IQC analysis for memoryless system
classdef testIqcAnalysisMemoryless < matlab.unittest.TestCase

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
function testErrorFreeAnalysis(testCase)
    options = AnalysisOptions('verbose', false, 'solver', 'sdpt3');
    del_types = {'DeltaSltv', 'DeltaBounded'};
    for i = 1:3
        req_deltas = cell(1, randi([1, length(del_types)]));
        for j = 1:length(req_deltas)
            req_deltas{j} = del_types{randi([1, length(del_types)])};
        end
        lft = Ulft.random('num_deltas', length(req_deltas),...
                          'req_deltas', req_deltas);
        [dim_out, dim_in] = size(lft);
        perf = lft.performance.performances{1};
        mult_perf = MultiplierL2Induced(perf, dim_out, dim_in);
        iqcAnalysis(lft,...
                    'analysis_options', options,...
                    'multipliers_performance', mult_perf);
    end
    no_error = true;
    verifyTrue(testCase, no_error)
end
end

methods (Test, TestTags = {'RCT'})
function testCorrectAnalysis(testCase)
for j = 1:5
    rct_object = zeros(3);
    for i = 1:3
        var = randatom('ureal');
        base = rand(3);
        base(base < .5) = 0;
        base(base >= .5) = 1;
        rct_object = rct_object + var * base;
    end    
    rct_result = wcgain(uss(rct_object));
    lft = rctToLft(rct_object);
    options = AnalysisOptions('verbose', false, 'solver', 'sdpt3');
    result = iqcAnalysis(lft, 'analysis_options', options);
    assumeTrue(testCase, isfinite(rct_result.LowerBound))
    assumeTrue(testCase, isfinite(rct_result.UpperBound))
    verifyGreaterThan(testCase, result.performance, rct_result.LowerBound * .99)
    verifyLessThan(testCase, result.performance, rct_result.UpperBound * 1.01)
end
end
end

end

%%  CHANGELOG
% Oct. 12, 2021: Added after v0.6.0 - Micah Fry (micah.fry@ll.mit.edu)