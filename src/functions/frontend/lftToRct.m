function uss_out = lftToRct(ulft_input)
%% lftToRct for converting Ulft objects to uss objects.
%
%     uss_out = lftToRct(ulft_input)
%
%     Variables:
%     ---------
%       Input:
%         ulft_input : ulft object
%       Output:
%         uss_out : uss_object :: the resultant uncertain state space
%         object
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
    
    uss_out = lft(uncertainties, ss_lft);
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Philip Mulford (philip.mulford@ll.mit.edu)