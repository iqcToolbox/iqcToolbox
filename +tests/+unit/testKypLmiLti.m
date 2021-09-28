%% Requirements:
%  1. kypLmiLti shall output the KYP-based LMI matrix associated with the 
%       terms: filter(jw)' * quad * filter(jw) 
%       (see "On the KYP lemma (Rantzer, 1997)")
%  1.1 kypLmiLti shall generate an LMI based off sdpvars from the input 
%       "kyp_var" and possibly from the input "quad"
%  1.2 kypLmiLti shall output the discrete-time or continuous-time LMI as
%       indicated by the input "filter"

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for kypLmiLti
classdef testKypLmiLti < matlab.unittest.TestCase
    methods (Test)
        function testDiscreteTime(testCase)
            dim_state = randi(10);
            dim_in = randi(10);
            dim_out = randi(10);
            lyap_matrix = sdpvar(dim_state);
            kap = sdpvar(1);
            
            g = drss(dim_state, dim_out, dim_in);
            g.a = 0.9 * g.a;
            g = g / norm(g, 'inf');
            gain = norm(g, 'inf');
                        
            lmi_mat = kypLmiLti([g; eye(dim_in)],...
                                blkdiag(eye(dim_out), -kap * eye(dim_in)),...
                                lyap_matrix);
            options = sdpsettings('showprogress', 0, 'verbose', 0);
            report = optimize([lmi_mat <= 1e-8], kap, options);
            kyp_gain = sqrt(double(kap));
            verifyLessThan(testCase, abs(kyp_gain - gain), 1e-3)
        end

        function testContinuousTime(testCase)
            dim_state = randi(10);
            dim_in = randi(10);
            dim_out = randi(10);
            lyap_matrix = sdpvar(dim_state);
            kap = sdpvar(1);
            
            g = rss(dim_state, dim_out, dim_in);
            g.a = g.a - 0.1 * eye(dim_state);
            g = g / norm(g, 'inf');
            gain = norm(g, 'inf');
            
            lmi_mat = kypLmiLti([g; eye(dim_in)],...
                                blkdiag(eye(dim_out), -kap * eye(dim_in)),...
                                lyap_matrix);
            options = sdpsettings('showprogress', 0, 'verbose', 0);
            optimize([lmi_mat <= 1e-8], kap, options);
            kyp_gain = sqrt(double(kap));
            verifyLessThan(testCase, abs(kyp_gain - gain), 1e-3)
        end
    end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0): Fixed creating of unstable test cases - Micah Fry
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)