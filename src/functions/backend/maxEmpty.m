function out = maxEmpty(in)
%% MAXEMPTY a quick-fix utility function to return 0 if calling the built-in
%  max() function on an empty array
%
%  out = maxEmpty(in)
%
%  Variables:
%  ---------
%     Input:
%       in : numeric array 
%     Output:
%       out : numeric scalar :: maximum of input
%
%  See also Ulft.Ulft

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

assert(isnumeric(in), 'maxEmpty:maxEmpty', 'Input to maxEmpty must be numeric')
if isempty(in)
    out = 0;
else
    out = max(in);
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)