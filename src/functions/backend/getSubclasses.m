function subclasses = getSubclasses(superclass)
%% GETSUBCLASSES returns an array of each subclass of the requested superclass as strings
%
%  subclasses = getSubclasses(superclass)
%
%  Variables:
%  ---------
%     Input:
%        superclass : string :: name of the class to get the subclasses of
%     Output:
%        subclasses : cell array of strings :: names of the subclasses

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology
%  SPDX-License-Identifier: GPL-2.0
%%

validateattributes(superclass, {'char'}, {'nonempty'});

directory = fileparts(which(superclass));
files = struct2cell(dir(directory));
files = files(1,:);

subclasses = cell(0);
for i = 1 : length(files)
    [~, name, ~] = fileparts(files{i});
    if exist(name, 'class') && ismember(superclass, superclasses(name))
        subclasses{end+1} = name;
    end
end

end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Jason Nezvadovitz (jason.nezvadovitz@ll.mit.edu)