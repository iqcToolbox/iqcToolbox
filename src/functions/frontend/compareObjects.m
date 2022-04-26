function compareObjects(obj1, obj2)
%% COMPAREOBJECTS function to help quickly determine where two objects are not equivalent
%  by printing which properties are equivalent, and which are not
%
%  compareObjects(obj1, obj2)
%
%  Variables:
%  ---------
%     Input:
%       obj1 : any MATLAB object
%       obj2 : any MATLAB object (as long as it is of the same class as obj1)

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

assert(isa(obj1, class(obj2)), 'Comparing objects of different classes!')

fnames = fieldnames(obj1);
sobj1 = struct(obj1);
sobj2 = struct(obj2);
for i = 1:length(fnames)
    nt = [];
    fname = fnames{i};
    if ~isequaln(getfield(sobj1, fname), getfield(sobj2, fname))
        nt = 'n''t';
    end
    disp(['Property ', fname, ' is',nt, ' equal for both objects'])
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)
