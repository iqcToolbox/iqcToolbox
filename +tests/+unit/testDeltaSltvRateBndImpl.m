%% Requirements:
%  1. DeltaSltvRateBndImpl shall be defined by its name, input
%      dimensions, the type of region in which the uncertainty and its rate
%      reside, a specification of that region, the horizon_period, and the
%      basis length desired for its multiplier.
%  2. Upon construction, and when queried by user, it shall display the
%      information described in (1).
%
%  3. The admissible region types of DeltaSltvRateBndImpl shall be a box, a
%      polytope described by vertices of the polytope, or a origin-centered
%      axes-symmetric ellipse.
%  4. If region type, or region specification is not provided 
%      by the user, by default the uncertainty and its rate-of-change shall
%      reside in a unit box with a horizon_period [0, 1].
%
%  5. If the user provides no name, DeltaSltvRateBndImpl shall 
%      throw an exception
%  6. If the user provides input dimensions that is not a natural
%      number DeltaSltvRateBndImpl shall throw an exception
%  7. If the user provides inadmissible specifications of the uncertainty 
%      region, DeltaSltvRateBndImpl shall throw an exception
%  8. For a box region, if the user provides a lower bound which is greater
%      than the upper bound, DeltaSltvRateBndImpl shall throw an exception
%  9. For a polytope region, if the convex hull of the vertices does not 
%      contain the origin, DeltaSltvRateBndImpl shall throw an exception
%  10.For an ellipse region, if the radii are not well-defined,
%      DeltaSltvRateBndImpl shall throw an exception
%
%  11.The in/out dimensions of DeltaSltvRateBndImpl shall be equal.
%
%  12.DeltaSltvRateBndImpl shall ensure that it's properties are consistent 
%      with its current horizon_period property
%  13.DeltaSltv shall be capable of changing it's properties to match a
%      newly input horizon_period, as long as the new horizon_period is
%      consistent with the prior horizon_period

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for DeltaSlti.
classdef testDeltaSltvRateBndImpl < matlab.unittest.TestCase
methods (Test)

function testFullBoxConstructor(testCase)
	name = 'test';
    dim_in = 5;
    region_type = 'box';
    lower_bound = -5;
    upper_bound = -3;
    lower_rate = -1;
    upper_rate = 0.2;
    region_data = {[lower_bound, upper_bound; lower_rate, upper_rate]};
    horizon_period = [0, 1];
    basis_length = 3;
    del = DeltaSltvRateBndImpl(name,...
                                         dim_in,...
                                         region_type,...
                                         region_data,...
                                         horizon_period,...
                                         basis_length);
    verifyEqual(testCase, del.name, name)
    verifyEqual(testCase, del.dim_in, dim_in)
    verifyEqual(testCase, del.dim_out, dim_in * basis_length)
    verifyEqual(testCase, del.horizon_period, horizon_period)
    verifyEqual(testCase, del.region_type, region_type)
    verifyEqual(testCase, del.region_data, region_data)
    verifyEqual(testCase, del.basis_length, basis_length)
    verifyEqual(testCase, del.lower_bound, lower_bound)
    verifyEqual(testCase, del.upper_bound, upper_bound)
    verifyEqual(testCase, del.lower_rate, lower_rate)
    verifyEqual(testCase, del.upper_rate, upper_rate)
    verifyEqual(testCase, del.vertices, nan)
    verifyEqual(testCase, del.ellipses, nan)
end

function testFullEllipseConstructor(testCase)
    name = 'test';
    dim_in = 5;
    region_type = 'ellipse';
    axes = [1; 2;];
    region_data = {axes};
    horizon_period = [3, 2];
    basis_length = 2;
    del = DeltaSltvRateBndImpl(name,...
                                         dim_in,...
                                         region_type,...
                                         region_data,...
                                         horizon_period,...
                                         basis_length);
    total_time = sum(horizon_period);
    verifyEqual(testCase, del.name, name)
    verifyEqual(testCase, del.dim_in, repmat(dim_in, 1, total_time))
    verifyEqual(testCase, del.dim_out, repmat(dim_in * basis_length, 1, total_time))
    verifyEqual(testCase, del.horizon_period, horizon_period)
    verifyEqual(testCase, del.region_type, region_type)
    verifyEqual(testCase, del.region_data, repmat(region_data, 1, total_time))
    verifyEqual(testCase, del.basis_length, basis_length)
    verifyEqual(testCase, del.ellipses, repmat({axes}, 1, total_time))
    verifyEqual(testCase, del.lower_bound, nan)
    verifyEqual(testCase, del.upper_bound, nan)
    verifyEqual(testCase, del.lower_rate, nan)
    verifyEqual(testCase, del.upper_rate, nan)
    verifyEqual(testCase, del.vertices, nan)
end

function testFullPolytopeConstructor(testCase)
    name = 'x';
    dim_in = [2];
    region_type = 'polytope';
    vertices = [ 1,  1, -1, -1;
                 1, -1,  1, -1];
    region_data = {vertices};
    horizon_period = [5, 3];
    total_time = sum(horizon_period);
    basis_length = 4;
    del = DeltaSltvRateBndImpl(name,...
                                         dim_in,...
                                         region_type,...
                                         {vertices},...
                                         horizon_period,...
                                         basis_length);
    verifyEqual(testCase, del.name, name)
    verifyEqual(testCase, del.dim_in, repmat(dim_in, 1, total_time))
    verifyEqual(testCase, del.dim_out, repmat(dim_in * basis_length, 1, total_time))
    verifyEqual(testCase, del.horizon_period, horizon_period)
    verifyEqual(testCase, del.region_type, region_type)
    verifyEqual(testCase, del.region_data, repmat(region_data, 1, total_time))
    verifyEqual(testCase, del.basis_length, basis_length)
    verifyEqual(testCase, del.vertices, repmat({vertices}, 1, total_time))
    verifyEqual(testCase, del.lower_bound, nan)
    verifyEqual(testCase, del.upper_bound, nan)
    verifyEqual(testCase, del.lower_rate, nan)
    verifyEqual(testCase, del.upper_rate, nan)
    verifyEqual(testCase, del.ellipses, nan)
end

function testOneArgConstructor(testCase)
    name = 'test';
    del = DeltaSltvRateBndImpl(name);
    verifyEqual(testCase, del.name, name)
    verifyEqual(testCase, del.dim_in, 1)
    verifyEqual(testCase, del.dim_out, 2)
    verifyEqual(testCase, del.horizon_period, [0, 1])
    verifyEqual(testCase, del.region_type, 'box')
    verifyEqual(testCase, del.region_data, {[-1, 1; -1, 1]})
    verifyEqual(testCase, del.basis_length, 2)
    verifyEqual(testCase, del.lower_bound, -1)
    verifyEqual(testCase, del.upper_bound, 1)
    verifyEqual(testCase, del.lower_rate, -1)
    verifyEqual(testCase, del.upper_rate, 1)
    verifyEqual(testCase, del.vertices, nan)
    verifyEqual(testCase, del.ellipses, nan)
end

function testTwoArgConstructor(testCase)
    name = 'test';
    dim_in = 4;
    del = DeltaSltvRateBndImpl(name, dim_in);
    verifyEqual(testCase, del.name, name)
    verifyEqual(testCase, del.dim_in, dim_in)
    verifyEqual(testCase, del.dim_out, dim_in * 2)
    verifyEqual(testCase, del.horizon_period, [0, 1])
    verifyEqual(testCase, del.region_type, 'box')
    verifyEqual(testCase, del.region_data, {[-1, 1; -1, 1]})
    verifyEqual(testCase, del.basis_length, 2)
    verifyEqual(testCase, del.lower_bound, -1)
    verifyEqual(testCase, del.upper_bound, 1)
    verifyEqual(testCase, del.lower_rate, -1)
    verifyEqual(testCase, del.upper_rate, 1)
    verifyEqual(testCase, del.vertices, nan)
    verifyEqual(testCase, del.ellipses, nan)
end

function testFourArgConstructor(testCase)
    name = 'test';
    dim_in = 7;
    region_type = 'box';
    lower_bound = -5;
    upper_bound = -3;
    lower_rate = -1;
    upper_rate = 0.2;
    region_data = {[lower_bound, upper_bound; lower_rate, upper_rate]};
    del = DeltaSltvRateBndImpl(name,...
                                         dim_in,...
                                         region_type,...
                                         region_data);
    verifyEqual(testCase, del.name, name)
    verifyEqual(testCase, del.dim_in, dim_in)
    verifyEqual(testCase, del.dim_out, dim_in * 2)
    verifyEqual(testCase, del.horizon_period, [0, 1])
    verifyEqual(testCase, del.region_type, region_type)
    verifyEqual(testCase, del.region_data, region_data)
    verifyEqual(testCase, del.basis_length, 2)
    verifyEqual(testCase, del.lower_bound, lower_bound)
    verifyEqual(testCase, del.upper_bound, upper_bound)
    verifyEqual(testCase, del.lower_rate, lower_rate)
    verifyEqual(testCase, del.upper_rate, upper_rate)
    verifyEqual(testCase, del.vertices, nan)
    verifyEqual(testCase, del.ellipses, nan)
end

function testFiveArgConstructor(testCase)
    name = 'test';
    dim_in = 4;
    region_type = 'ellipse';
    region_data = {[3; 5]};
    horizon_period = [3, 1];
    del = DeltaSltvRateBndImpl(name,...
                                         dim_in,...
                                         region_type,...
                                         region_data,...
                                         horizon_period);
    total_time = sum(horizon_period);
    verifyEqual(testCase, del.name, name)
    verifyEqual(testCase, del.dim_in, repmat(dim_in, 1, total_time))
    verifyEqual(testCase, del.dim_out, repmat(dim_in * 2, 1, total_time))
    verifyEqual(testCase, del.horizon_period, horizon_period)
    verifyEqual(testCase, del.region_type, region_type)
    verifyEqual(testCase, del.region_data, repmat(region_data, 1, total_time))
    verifyEqual(testCase, del.basis_length, 2)
    verifyEqual(testCase, del.lower_bound, nan)
    verifyEqual(testCase, del.upper_bound, nan)
    verifyEqual(testCase, del.lower_rate, nan)
    verifyEqual(testCase, del.upper_rate, nan)
    verifyEqual(testCase, del.vertices, nan)
end

% function testDeltaToMultiplier(testCase)
%     name = 'test';
%     delta_sltv = DeltaSltv(name);
%     default_mult = deltaToMultiplier(delta_sltv);
%     verifyEqual(testCase, default_mult.quad_time_varying, true)
%     
%     quad_time_varying = false;
%     mult = deltaToMultiplier(delta_sltv,...
%                              'quad_time_varying',...
%                              quad_time_varying);
%     verifyEqual(testCase, mult.quad_time_varying, quad_time_varying)
% end

function testFailedName(testCase)
    verifyError(testCase, @() DeltaSltvRateBndImpl(), ?MException)
    verifyError(testCase, @() DeltaSltvRateBndImpl(1), ?MException)
end

function testFailedDimension(testCase)
    verifyError(testCase, @() DeltaSltvRateBndImpl('test', -4), ?MException)
    verifyError(testCase, @() DeltaSltvRateBndImpl('test', 5.1), ?MException)
end

function testFailedBox(testCase)
    hp = [2, 2];
    total_time = sum(hp);
    lb = -1.2 * ones(1, 2 * total_time);
    ub = linspace(4, 5, 2 * total_time);
    region_type = 'box';
    region_data = mat2cell([lb', ub'], 2 * ones(1, total_time), 2)';
%     assert
    try 
        DeltaSltvRateBndImpl('x', 1, region_type, region_data, hp);
    catch
        assertFail(testCase);
    end
    region_data{end}(1, 1) = region_data{end}(1, 2) + 1;
    verifyError(testCase,...
                @() DeltaSltvRateBndImpl('x', 1, region_type, region_data, hp),...
                ?MException)

    region_data{end}(1) = -inf;
    verifyError(testCase,...
                @() DeltaSltvRateBndImpl('x', 1, region_type, region_data, hp),...
                ?MException)

    region_data{total_time}(1) = ub(total_time) - 1;
    region_data{total_time}(2) = nan;
    verifyError(testCase,...
                @()DeltaSltvRateBndImpl('x', 1, region_type, region_data, hp),...
                ?MException)
            
    region_data = repmat({[-ones(3, 1), ones(3, 1)]}, 1, total_time);
    verifyError(testCase,...
            @()DeltaSltvRateBndImpl('x', 1, region_type, region_data, hp),...
            ?MException)

    region_data = repmat({[-ones(1, 1), ones(1, 1)]}, 1, total_time);
    verifyError(testCase,...
            @()DeltaSltvRateBndImpl('x', 1, region_type, region_data, hp),...
            ?MException)
end
% 
function testFailedEllipse(testCase)
    name = 'test';
    hp = [1, 3];
    total_time = sum(hp);
    dim_in = 3;
    region_type = 'ellipse';
    region_data = mat2cell([ones(1, total_time); 3 * ones(1, total_time)],...
                           2,...
                           ones(1, total_time));
    try
        DeltaSltvRateBndImpl(name,...
                                     dim_in,...
                                     region_type,...
                                     region_data,...
                                     hp);
    catch
        assertFail(testCase);
    end
        
    region_data{end}(1) = 0;
    verifyError(testCase,...
                @()DeltaSltvRateBndImpl(name,...
                                     dim_in,...
                                     region_type,...
                                     region_data,...
                                     hp),...
                ?MException)
    
    region_data{end} = [1; NaN];
    verifyError(testCase,...
                @()DeltaSltvRateBndImpl(name,...
                                     dim_in,...
                                     region_type,...
                                     region_data,...
                                     hp),...
                ?MException)
            
    region_data{end} = [inf; 2];
    verifyError(testCase,...
                @()DeltaSltvRateBndImpl(name,...
                                     dim_in,...
                                     region_type,...
                                     region_data,...
                                     hp),...
                ?MException)  
            
    region_data = repmat({ones(3, 1)}, 1, total_time);
    verifyError(testCase,...
            @()DeltaSltvRateBndImpl('x', 1, region_type, region_data, hp),...
            ?MException)

    region_data = repmat({ones(1, 1)}, 1, total_time);
    verifyError(testCase,...
            @()DeltaSltvRateBndImpl('x', 1, region_type, region_data, hp),...
            ?MException)
end

function testFailedPolytope(testCase)
    names = 'test';
    dim_in = 2;
    region_type = 'polytope';
    region_data = repmat({[1, 1, -1, -1; 1, -1,  1, -1]}, 1, 2);
    hp = [1, 1];
    try
        DeltaSltvRateBndImpl(names,...
                                     dim_in,...
                                     region_type,...
                                     region_data,...
                                     hp);
    catch
        assertFail()
    end
    
    region_data = repmat({[1, 1, -1, -1; 1, 3,  1, 3]}, 1, 2);
    verifyError(testCase,...
                @()DeltaSltvRateBndImpl(names,...
                                     dim_in,...
                                     region_type,...
                                     region_data,...
                                     hp),...
                ?MException)
end    

function testHorizonPeriod(testCase)
   name = 'test';
   del = DeltaSltvRateBndImpl(name);
   assertEqual(testCase, del.horizon_period, [0, 1])
   
   % Checking horizon_period and making sure it fits for all properties
   hp2 = [4, 3];
   total_time2 = sum(hp2);
   region_type = 'polytope';
   bnd = linspace(1, 2, total_time2);
   for i = 1:total_time2
       region_data{i} = [[-bnd(i); -bnd(i)], [bnd(i); bnd(i)]];
   end
   dim_in = 3;
   del = DeltaSltvRateBndImpl(name, dim_in, region_type, region_data, hp2);
   del = matchHorizonPeriod(del);
   assertEqual(testCase, del.horizon_period, hp2)
   assertEqual(testCase, del.vertices, region_data)
   assertEqual(testCase, del.dim_in, repmat(dim_in, 1, total_time2))

   % Resetting horizon_period and making sure it fits for all properties
   horizon_period3 = [total_time2, hp2(2) * 2];
   total_time3 = sum(horizon_period3);
   del = matchHorizonPeriod(del, horizon_period3);
   verifyEqual(testCase, del.horizon_period, horizon_period3)
   verifyEqual(testCase,...
               del.vertices,...
               [region_data,...
                region_data((1:(horizon_period3(2) / 2)) + hp2(1)),...
                region_data((1:(horizon_period3(2) / 2)) + hp2(1))])
   verifyEqual(testCase, del.dim_in,  repmat(sum(dim_in), 1, total_time3))
end

function testFailedHorizonPeriod(testCase)
   name = 'test';
   del = DeltaSltvRateBndImpl(name);
   assertEqual(testCase, del.horizon_period, [0, 1])
   
   % Resetting horizon_period and making sure it fits for all properties
   horizon_period2 = [4, 7];
   del.horizon_period = horizon_period2;
   del = matchHorizonPeriod(del);
   assertEqual(testCase, del.horizon_period, horizon_period2)
   
   % Resetting horizon_period and incorrectly trying to force a fit with
   % other properties
   horizon_period3 = [5, 14];
   horizon_period3 = commonHorizonPeriod([horizon_period2; horizon_period3]);
   del.horizon_period = horizon_period3;
   verifyError(testCase, @() matchHorizonPeriod(del), ?MException)
end

function testDeltaToMultiplier(testCase)
    name = 'test';
    del = DeltaSltvRateBndImpl(name);
    default_mult = deltaToMultiplier(del);
    verifyEqual(testCase, default_mult.quad_time_varying, true)
    verifyEqual(testCase, default_mult.discrete, true)
    
    quad_time_varying = false;
    is_discrete = false;
    mult = deltaToMultiplier(del,...
                            'discrete', is_discrete,...
                            'quad_time_varying', quad_time_varying);
    verifyEqual(testCase, mult.discrete, is_discrete)
end

function testRecastMatricesAndDelta(testCase)
    name = 'test';
    del = DeltaSltvRateBndImpl(name);
    [recastA, recastB, recastC, newDelta] = recastMatricesAndDelta(del);
    verifyEmpty(testCase, recastA)
    verifyEmpty(testCase, recastB)
    verifyEmpty(testCase, recastC)
    verifyEmpty(testCase, newDelta)
end

end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)