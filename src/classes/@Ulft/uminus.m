function minus_lft = uminus(this_lft)
%% - UMINUS overloaded function for Ulft objects.
%
%     this_lft = -this_lft
%     this_lft = uminus(this_lft)
%
%     Variables:
%     ---------
%       Input:
%         this_lft : Ulft object :: the lft on the right side of the "-" symbol
%       Output:
%         minus_lft : Ulft object :: the negated lft
%
%     See also Ulft, minus, plus, uplus, mtimes.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

minus_lft = this_lft;
minus_lft.c = cellfun(@uminus, minus_lft.c, 'UniformOutput', false);
minus_lft.d = cellfun(@uminus, minus_lft.d, 'UniformOutput', false);
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)