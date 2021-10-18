%% Test script to check displays of all objects

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for DeltaSlti.
classdef testDisplays < matlab.unittest.TestCase
methods (Test)

function testDefault(testCase)
    z = toLft(DeltaDelayZ());
    dlti = toLft(DeltaDlti('dlti'));
    slti = toLft(DeltaSlti('slti'));
    sltv = toLft(DeltaSltv('sltv'));
    sltv_box = toLft(DeltaSltvRepeated('sltv_box'));
    sltv_ell = toLft(DeltaSltvRepeated('sltv_ell', 1, 'ellipse', {1}));
    sltv_polytope = ...
        toLft(DeltaSltvRepeated('sltv_poly', 1, 'polytope', {[-1, 1]}));
    sltv_rb = toLft(DeltaSltvRateBnd('sltv_rate_bound'));
    sltv_rbi_ell = toLft(DeltaSltvRateBndImpl('sltv_rbi_ell',...
                                                        2,...
                                                        'ellipse',...
                                                        {[1; 2]}));
    sltv_rbi_pol = toLft(DeltaSltvRateBndImpl('sltv_rbi_poly',...
                                                        2,...
                                                        'polytope',...
                                                        {[-1, -1, 1, 1;
                                                          -1,  1,-1, 1]}));
    sltv_rbi_box = toLft(DeltaSltvRateBndImpl('sltv_rbi_box',...
                                                        2,...
                                                        'box',...
                                                        {[-1, 1; -1, 1]}));
    passive = toLft(DeltaPassive('passive'));
    lft_discrete = z + dlti + slti + sltv + sltv_box + sltv_ell + ...
                       sltv_polytope + sltv_rb
    lft_discrete_2 = sltv_rbi_ell + sltv_rbi_pol + sltv_rbi_box
    
    integrator = toLft(DeltaIntegrator());
    lft_continuous = integrator + sltv + sltv_box + sltv_rb + passive
end

end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)