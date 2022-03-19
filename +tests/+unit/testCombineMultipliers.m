%% Requirements:
%  1. combineAllMultipliers shall take MultiplierDeltaCombined,
%     MultiplierDisturbanceCombined, and MultiplierPerformanceCombined, and 
%     output a MultiplierPerformanceCombined which conglomerates all inputs

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for PerformanceL2Induced
classdef testCombineMultipliers < matlab.unittest.TestCase
methods (Test)
    function testCombiningMultipliers(testCase)
        dim_in  = 1;
        dim_out = 1;
        
        del_slti  = DeltaSlti('del');
        m_slti_d  = deltaToMultiplier(del_slti, 'discrete', true);
        m_slti_c  = deltaToMultiplier(del_slti, 'discrete', false);
        
        perf_l2   = PerformanceL2Induced('perf');
        m_l2      = performanceToMultiplier(perf_l2,...
                                            'dim_out_lft', dim_out,...
                                            'dim_in_lft', dim_in);
        
        dis_wh = DisturbanceBandedWhite('dis');
        m_wh   = disturbanceToMultiplier(dis_wh, 'dim_in_lft', dim_in);

        m_combined = combineAllMultipliers(MultiplierDeltaCombined(m_slti_d),...
                                           MultiplierDisturbanceCombined(m_wh),...
                                           MultiplierPerformanceCombined(m_l2),...
                                           dim_out);
        verifyEqual(testCase, m_combined.name, {m_slti_d.name, m_l2.name, m_wh.name})
        verifyError(testCase,...
                    @() combineAllMultipliers(MultiplierDeltaCombined(m_slti_c),...
                                              MultiplierDisturbanceCombined(m_wh),...
                                              MultiplierPerformanceCombined(m_l2),...
                                              dim_out),...
                    'combineAllMultipliers:combineAllMultipliers')
        
    end
end
end