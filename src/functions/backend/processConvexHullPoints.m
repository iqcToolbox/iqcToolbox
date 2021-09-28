function [valid, points] = processConvexHullPoints(points)
%% PROCESSCONVEXHULLPOINTS function which determines if the convex hull of 
%  a set of points is valid, and removes any interior or extraneous points
%
%  [valid, points] = processConvexHullPoints(points)
%
%  Variables:
%  ---------
%     Input:
%       points : numeric array :: an n x m array of m points in R^n
%     Output:
%       valid : boolean :: indicates if convex hull contains the origin
%       points : numeric array :: an n x (m - p) array of points in R^n, where
%                                 any points interior to the convex hull are
%                                 removed from the set
%
%  See also DeltaSltvRepeated

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

validateattributes(points, 'numeric', {'nonempty',...
                                      'finite',...
                                      '2d'})
dim_points = size(points, 1);
valid = false;
containsOrigin = @(points) pointInConvexHull(zeros(dim_points, 1), points);
if containsOrigin(points)
    valid = true;
    points_reduced = points;
    i = 1;
    while i <= size(points_reduced, 2)
        point = points_reduced(:, i);
        points_trial = points_reduced(:, [1:i-1,i+1:end]);
        if pointInConvexHull(point, points_trial)
            points_reduced = points_trial;
        else
            i = i + 1;
        end
        points = points_reduced;
    end
end

end

function inside = pointInConvexHull(point, points)
[dim_points, n_points] = size(points);
validateattributes(point, 'numeric', {'nonempty',...
                                      'finite',...
                                      'size', [dim_points, 1]})
validateattributes(points, 'numeric', {'nonempty', 'finite'})

Aeq = [points; ones(1, n_points)];
beq = [point; 1];
lb = zeros(n_points, 1);
x = sdpvar(n_points, 1);
ct = [Aeq * x == beq,  lb <= x];
options = sdpsettings('solver', 'lpsolve', 'verbose', false);
yalmip_report = optimize(ct, [], options);

inside = yalmip_report.problem == 0;
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0): Replaced dependence on Optimization Toolbox with LPSOLVE - Micah Fry
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)