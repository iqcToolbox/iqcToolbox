function vertices_time = boxToVertices(box)
%% BOXTOVERTICES function which creates vertices whose convex hull is
%  equivalent to the box input.
%
%  vertices_time = boxToVertices(box)
%
%  Variables:
%  ---------
%     Input:
%       box : cell array of matrices :: the time-varying upper/lower bounds
%                                         of each dimension describing a box
%     Output:
%       vertices_time : (1 x total_time) cell array of matrices :: 
%                          each column of each matrix represents a vertex describing the box
%
%  See also DeltaSltvRepeated

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

validateattributes(box, 'cell', {'nonempty'})
total_time = length(box);
dims = size(box{1}, 1);

for i = 1 : total_time
    validateattributes(box{i},...
                       'numeric',...
                       {'size', [dims, 2], 'nonnan', 'finite'})
    assert(all(box{i}(:, 1) <= box{i}(:, 2), 'all'),...
           'boxToVertices:boxToVertices',...
           ['region_data for "box" must have lower bounds (:, 1)',...
            ' less than upper bounds (:, 2)'])
    assert(all(abs(box{i}(:, 1) + box{i}(:, 2)) < 1e-8,'all'),...
           'boxToVertices:boxToVertices',...
           'boxToVertices currently does not support asymetric bounds',...
           ' (i.e., lower_bound must equal -upper_bound). Normalize')
end

signs = ones(dims, 1);
flips = 2 .^ (1:dims) / 2;
vertices_time = cell(1, total_time);
for i = 1 : total_time
    vertices = nan(dims, 2 ^ dims);
    for j = 1 : 2 ^ dims
        for k = 1 : dims
            if ~mod(j, flips(k))
                signs(k) = -signs(k);
            end
        end
        vertex = box{i}(:, 2) .* signs;
        vertices(:, j) = vertex;
    end
    vertices_time{i} = circshift(vertices, [0, 1]);
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)