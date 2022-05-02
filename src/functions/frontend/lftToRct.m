function [uss_out, delta_map] = lftToRct(ulft_input)
%% lftToRct for converting Ulft objects to uss objects.
%
%     [uss_out, delta_map] = lftToRct(ulft_input)
%
%     Variables:
%     ---------
%       Input:
%         ulft_input : ulft object
%       Output:
%         uss_out : uss_object :: the resultant uncertain state space object
%         delta_map: containers.Map object :: a mapping between names of
%                       of uncertainties, and the original Delta object.
%                       This is only needed when wanting to preserve information
%                       on uncertainties that are not part of the Robust Control
%                       Toolbox, because the uss object will merely create
%                       nondescript udyn objects as
%
%     See also Ulft, Ulft.Ulft, rctToLft, toLft

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

    % generates a m matrix to allow for manipulation to match uss expected
    % input
    validateattributes(ulft_input, {'Ulft'}, {'nonempty'})
    good_horizon = @(ulft) isequal(ulft_input.horizon_period, [0 1]);
    assert(good_horizon(ulft_input),...
           'lftToRct:lftToRct',...
           strcat('RCT does not support eventually-periodic systems. ',...
           ' Ulft systems must have a [0 1] horizon period.'))
    m = [ulft_input.a, ulft_input.b; ulft_input.c, ulft_input.d];
    m = cell2mat(m);
    
    delta = ulft_input.delta;
    delta = delta.deltas;
    uncertainties = cell(size(delta, 2),1);
    for i = 1:size(delta, 2)
        uncertainties{i} = toRct(delta{i});
    end

    uncertainties = blkdiag(uncertainties{:});
    
    % uss A matrix only holds the information the state deltas - if there
    % are none, the dimension will be 0
    dim_in = 0;
    dim_out = 0;
    if ~isempty(delta) && isa(delta{1}, 'DeltaDelayZ') ||...
       ~isempty(delta) && isa(delta{1}, 'DeltaIntegrator')
            delta_integrator = delta{1};
            dim_in = delta_integrator.dim_in;
            dim_out = delta_integrator.dim_out;
    end
    
    ts = 0;
    
    if ~isempty(delta) && isa(delta{1}, 'DeltaDelayZ')
        ts = delta{1}.timestep;
    end
    
    a = m(1:dim_in, 1:dim_out);
    b = m(1:dim_in, dim_out+1:end);
    c = m(dim_in+1:end, 1:dim_out);
    d = m(dim_in+1:end, dim_out+1:end);
    
    ss_lft = ss(a, b, c, d, ts);
    delta_map = delta_mapper(ulft_input);
    uss_out = lft(uncertainties, ss_lft);
end

function mapped_delta = delta_mapper(in_lft)
    mapped_delta = containers.Map;
    for i = 1:length(in_lft.delta.deltas)
        mapped_delta(in_lft.delta.deltas{i}.name) = in_lft.delta.deltas{i};
    end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Philip Mulford (philip.mulford@ll.mit.edu)