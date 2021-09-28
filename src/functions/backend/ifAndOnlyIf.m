function true_false = ifAndOnlyIf(bool1, bool2)
%% IFANDONLYIF basic boolean-valued function for iff
%
%  true_false = ifAndOnlyIf(bool1, bool2)
%
%  Variables:
%  ---------
%     Input:
%        bool1 : boolean 
%        bool2 : boolean
%     Output:
%        true_false : boolean :: true if bool1 <==> bool2, false otherwise

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

validateattributes(bool1, 'logical', {'nonempty'})
validateattributes(bool2, 'logical', {'nonempty'})

true_false = (~bool1 || bool2) && (~bool2 || bool1);
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)