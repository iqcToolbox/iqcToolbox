function ss_out = lftToSs(in)
%% lftToSs for converting Ulft objects to ss objects.
%
%     ss_out = lftToSs(ulft_input)
%
%     Variables:
%     ---------
%       Input:
%         in : ulft object
%       Output:
%         ss_out : ss object :: the resultant state space object
%
%     See also Ulft, Ulft.Ulft, toLft, lftToSs

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

validateattributes(in, {'Ulft'}, {'nonempty'})
good_horizon = @(ulft) isequal(in.horizon_period, [0 1]);
assert(good_horizon(in),...
       'lftToSs:lftToSs',...
       strcat('ss objects do not support eventually-periodic systems. ',...
       ' Ulft objects must have a [0 1] horizon period.'))
assert(~in.uncertain,...
       'lftToSs:lftToSs',...
       ['Uncertain Ulft objects cannot be converted to ss objects. ',...
        'Consider using lftToRct.'])

if isempty(in.timestep)
    ss_out = ss(in.d{1});
elseif isfinite(in.timestep) && (in.timestep >= 0 || in.timestep == -1)
    ss_out = ss(in.a{1}, in.b{1}, in.c{1}, in.d{1}, in.timestep);
else
    error('lftToSs:lftToSs',...
          ['The timestep of the given LFT is inadmissible. ',...
           'It must either be empty, -1, or nonnegative'])
end        
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0): Added after v0.5.0 - Micah Fry (micah.fry@ll.mit.edu)