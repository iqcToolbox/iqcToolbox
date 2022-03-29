%% RUN_TESTS_SCRIPT for running tests
%  tests.run_tests_script() will run all tests

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Setup plugins and test suite
import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.selectors.HasTag
import matlab.unittest.constraints.IsEqualTo
import matlab.unittest.plugins.CodeCoveragePlugin
import matlab.unittest.plugins.codecoverage.CoberturaFormat
runner = TestRunner.withTextOutput;
top_path = mfilename('fullpath');
top_path(end - length(mfilename):end) =  [];
top_path = fullfile(top_path,'..','src');
quick_tests = {'tests.iqc_analysis.testIqcAnalysisDlti',...
               'tests.iqc_analysis.testIqcAnalysisPassive',...
               'tests.iqc_analysis.testIqcAnalysisSectorBounded',...
               'tests.iqc_analysis.testIqcAnalysisSltv',...
               'tests.iqc_analysis.testIqcAnalysisSltvRateBnd',...
               'tests.iqc_analysis.testIqcAnalysisUncertainInitialCondition',...
               'tests.iqc_analysis.testIqcAnalysisDisturbanceConstant',...
               'tests.iqc_analysis.testIqcAnalysisPerformanceStable/testPassiveTheorems',...
               'tests.iqc_analysis.testIqcAnalysisReachability'};
suite = testsuite(quick_tests);
% Check if user has Robust Control Toolbox, skip RCT-dependent tests if not
try
    ureal;
    rct_tbx = true;
catch
    rct_tbx = false;
end
if ~rct_tbx
    suite = suite.selectIf(~HasTag(IsEqualTo('RCT'))); 
end
% Check if user has Signal Processing Toolbox, skip SGT-dependent tests if not
try
    butter(2, 0.5);
    sgt_tbx = true;
catch
    sgt_tbx = false;
end
if ~sgt_tbx
    suite = suite.selectIf(~HasTag(IsEqualTo('SGT'))); 
end

%% Run all tests in suite, save results                  
result = runner.run(suite);
dtime = datetime;
save('+tests/quick_test_results', 'result', 'dtime')

%%  CHANGELOG
% Mar. 21, 2021: Added after v0.6.0 - Micah Fry (micah.fry@ll.mit.edu)