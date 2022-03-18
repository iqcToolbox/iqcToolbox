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

cov = CodeCoveragePlugin.forFolder(top_path,...
                                   'IncludingSubfolders', true,...
                                   'Producing', CoberturaFormat('+tests/coverage.xml'));
runner.addPlugin(cov)
suite = TestSuite.fromPackage('tests', 'IncludingSubpackages', true);
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
save('+tests/test_results', 'result', 'dtime')

%% To run a test from a specific class
% result = runner.run(TestSuite.fromClass(?tests.unit.testDeltaSltvRateBndImpl));
% result = runner.run(TestSuite.fromClass(?tests.iqc_analysis.testIqcAnalysisSltv));

%% To begin interactive testing
% testCase = matlab.unittest.TestCase.forInteractiveUse

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)
