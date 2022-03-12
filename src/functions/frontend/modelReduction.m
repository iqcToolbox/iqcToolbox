function [lft_red, error_out] = modelReduction(lft_in, desired_error)
%% modelReduction for reducing the dimension of deltas in ulft objects.
%
%     lft_red = modelReduction(lft_in, desired_error)
%     
%
%     Variables:
%     ---------
%       Input:
%         lft_in : input Ulft object with deltas to be reduced
%
%         desired_error : acceptable max error post-reduction
%             
%       Output:
%         lft_red : Ulft object :: the reduced lft
%
%       Note: we can only reduce LFTs with DeltaSlti, DeltaSltv,
%       DeltaSltvRateBounded, and DeltaDelayZ deltas
%       For stable systems, we use a balanced truncation approach. 
%       For unstable systems, we use Coprime Factors Reduction
%
%     See also Ulft, Ulft.Ulft.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

    %% Process Inputs
    input_parser = inputParser;
    addRequired(input_parser,...
               'lft_in',...
               @(lft_in) validateattributes(lft_in, {'Ulft'}, {'nonempty'}));
    
    % Error Bound Checking
    addOptional(input_parser,...
               'desired_error',...
               0.01,...
               @(err) validateattributes(err, {'numeric'}, {'nonnegative'}));
           
    parse(input_parser, lft_in, desired_error)
    lft_in        = input_parser.Results.lft_in;
    desired_error = input_parser.Results.desired_error;
    
    
    %Establishes iqcToolbox-specific constraints
    assert(isempty(lft_in.timestep) || ~isequal(lft_in.timestep, 0),...
           'modelReduction:modelReduction',...
           'The given LFT must be memoryless or discrete-time');
    isSquare = @(delta) isequal(delta.dim_in, delta.dim_out);
    goodDelta = @(delta) isa(delta, 'DeltaSlti') ...
                         || isa(delta, 'DeltaSltv') ...
                         || isa(delta, 'DeltaSltvRateBounded') ...
                         || isa(delta, 'DeltaDelayZ');
                     
    for i = 1:size(lft_in.delta.deltas, 1)
        if ~isa(lft_in.delta.deltas{i}, 'DeltaDelayZ')
            assert(isSquare(lft_in.delta.deltas{i}),...
               'modelReduction:modelReduction',...
               'All deltas must be square for balanced truncation');
        end
        assert(goodDelta(lft_in.delta.deltas{i}),...
               'modelReduction:modelReduction',...
               strcat('Delta must be one of a DeltaDelayZ, DeltaSlti,',...
                      ' DeltaSltv, DeltaSltvRateBounded'));
    end
    
    lft_norm = normalizeLft(lft_in);
    if ~isequaln(lft_norm, lft_in)
        lft_in = lft_norm;
        warning('modelReduction:modelReduction',...
                strcat('Input LFT is not normalized. Normalzing LFT to',...
                       ' apply model reduction algorithm'))
    end
    
    %% Check if strongly stable, get Right Coprime Factorization if not
    
    is_stable = checkStronglyStable(lft_in);
    
    if is_stable
        lft_trunc_ind = lft_in;
    else
        % If unstable, truncation indices will be determined by RCF
        lft_trunc_ind = rightCoprimeFactorization(lft_in);
        assert(checkStronglyStable(lft_trunc_ind),...
               'modelReduction:modelReduction',...
               ['The produced coprime factorization is not stable, cannot ',...
                'reduce model'])
    end
    
    %% Balance LFT
    [x_big, y_big] = computeGramians(lft_trunc_ind);
    [lft_trunc_ind, x_bal] = balance(lft_trunc_ind, x_big, y_big);
    
    %% Get truncation indices
    [ind_trunc, error_out] = truncationIndices(lft_trunc_ind, x_bal, desired_error);
    
    %% Truncate LFT
    if is_stable
        lft_bal = lft_trunc_ind;
    else
        % If unstable, must also balance the original LFT by RCF grammians
        lft_bal = balance(lft_in, x_big, y_big);
        error_out = -error_out;
        warning('modelReduction:modelReduction',...
                ['Error bound pertains to right coprime factor of LFT, ',...
                 'and cannot be used to bound the truncation error for ',...
                 'the original LFT. Error bound is returned as a negative ',...
                 'value to reflect its inapplicability.'])
    end
    lft_red = truncate(lft_bal, ind_trunc);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% First tier subfunctions %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [is_stable] = checkStronglyStable(lft_in)
%% CHECKSTRONGLYSTABLE determines if an lft is strongly stable or not
    [a, ~, ~, ~] = getABCD(lft_in);
    total_time = sum(lft_in.horizon_period);
    num_deltas = size(lft_in.delta.deltas, 2);
    dim_ins = lft_in.delta.dim_ins;
    
    % Formulate sub-block and full p's
    ct_stable = [];
    p = cell(num_deltas, total_time);
    for j = 1:total_time
        for i = 1:num_deltas
            p{i, j} = sdpvar(dim_ins(i, j));
            ct_stable = ct_stable + ((p{i, j} >= 1e-8 * eye(size(p{i, j}))):...
                                     ['PD stability, time: ', num2str(j),...    
                                      'delta: ', num2str(i)]);                  %#ok<BDSCA>
        end
    end
    
    p_big = cell(1, total_time);
    for j = 1:total_time
        p_big{j} = blkdiag(p{:, j});
    end
    p_big_plus = shiftForward(p, lft_in.timestep, lft_in.horizon_period);

    % Formulate LMIs and solve
    for j = 1:total_time
        lmi_mat = a{j} * p_big{j} * a{j}' - p_big_plus{j};
        ct_stable = ct_stable + ((lmi_mat <= -1e-8 * eye(size(lmi_mat))):...
                                 ['Stability LMI, time: ', num2str(j)]);        %#ok<BDSCA>
    end
    yalmip_settings = sdpsettings('verbose', false);
    diagnostic = optimize(ct_stable, [], yalmip_settings);
    
    % System is strongly stable if optimizer returns valid solution
    primal_residual = check(ct_stable);
    is_stable = diagnostic.problem ~= 1 ...
                && all(~isnan(primal_residual), 'all') ...
                && all(primal_residual >= -1e-8);
end

function lft_rcf = rightCoprimeFactorization(lft_in)
%% RIGHTCOPRIMEFACTORIZATION formulates the right coprime factor of an unstable 
%  lft by checking if it is stabilizable, and then constructing a stabilizing
%  feedback matrix
    [a, b, c, d] = getABCD(lft_in);
    total_time   = sum(lft_in.horizon_period);
    num_deltas   = length(lft_in.delta.deltas);
    dim_ins      = lft_in.delta.dim_ins;
    
    % Generate stabilizability matrix P
    p = cell(num_deltas, total_time);
    ct_p = [];
    
    for j = 1:total_time
        for i = 1:num_deltas
            p{i, j} = sdpvar(dim_ins(i, j));
            ct_p = ct_p + ((p{i, j} >= 1e-7 * eye(size(p{i, j}))):...
                           ['PD p time: ' num2str(j), ' delta: ' num2str(i)]);  %#ok<BDSCA>
        end
    end
    p_big = cell(1, total_time);
    for j = 1:total_time
        p_big{j} = blkdiag(p{:, j});
    end
    p_big_plus = shiftForward(p, lft_in.timestep, lft_in.horizon_period);
    
    for j = 1:total_time
        lmi_mat = a{j} * p_big{j} * a{j}' - p_big_plus{j} - b{j} * b{j}';
        ct_p = ct_p + ((lmi_mat <= -1e-7 * eye(size(lmi_mat))):...
                       ['P LMI constraints, timestep: ' num2str(j)]);           %#ok<BDSCA>
    end
    
    yalmip_settings = sdpsettings('verbose', false);
    diagnostic = optimize(ct_p, [], yalmip_settings);
    primal_residual = check(ct_p);
    is_stabilizable = diagnostic.problem ~= 1 ...
                      && all(~isnan(primal_residual), 'all') ...
                      && all(primal_residual >= -1e-7);
    assert(is_stabilizable,...
           'modelReduction:rightCoprimeFactorization',...
           'Unstable system is not stabilizable, unable to reduce model')
            
    f = cell(1, total_time);    
    for j = 1:total_time
        p_plus = value(p_big_plus{j});
        f{j} = -pinv(b{j}' / p_plus * b{j}) * b{j}' / p_plus * a{j};
    end
    
    %Generate A, C, D matrices for the RCf realization 
    a_rcf = cell(1, total_time);
    c_rcf = cell(1, total_time);
    d_rcf = cell(1, total_time);
    for j = 1:total_time
        a_rcf{j} = a{j} + b{j} * f{j};
        c_rcf{j} = [c{j} + d{j} * f{j}; f{j}];
        d_rcf{j} = [d{j}; eye(size(f{j}, 1))];
    end
    
    lft_rcf = Ulft(a_rcf, b, c_rcf, d_rcf, lft_in.delta,...
                   'horizon_period', lft_in.horizon_period,...
                   'disturbance', lft_in.disturbance,...
                   'performance', lft_in.performance);
end

function [x, y] = computeGramians(lft_in)
%% COMPUTEGRAMIANS takes a stable LFT and returns unbalanced controllability
%  and observability gramians for the LFT
    lmi_shift = 5e-7;
    [a, b, c] = getABCD(lft_in);
    total_time = sum(lft_in.horizon_period);
    dim_ins = lft_in.delta.dim_ins;
    num_deltas = length(lft_in.delta.deltas);
    
    % Formulate gramians
    x = cell(num_deltas, total_time);
    y = cell(num_deltas, total_time);
    ct_x = [];
    ct_y = [];
    for j = 1:total_time
        for i = 1:num_deltas
            x{i, j} = sdpvar(dim_ins(i, j));
            y{i, j} = sdpvar(dim_ins(i, j));
            ct_x = ct_x + ((x{i, j} >= lmi_shift * eye(size(x{i, j}))):...
                           ['PD x time: ' num2str(j), ' delta: ' num2str(i)]);  %#ok<BDSCA>
            ct_y = ct_y + ((y{i, j} >= lmi_shift * eye(size(y{i, j}))):...
                           ['PD y time: ' num2str(j), ' delta: ' num2str(i)]);  %#ok<BDSCA>
        end
    end
    x_big = cell(1, total_time);
    y_big = cell(1, total_time);
    obj_x = 0;
    obj_y = 0;
    for j = 1:total_time
        x_big{j} = blkdiag(x{:, j});
        y_big{j} = blkdiag(y{:, j});
    end
    x_big_plus = shiftForward(x, lft_in.timestep, lft_in.horizon_period);
    y_big_plus = shiftForward(y, lft_in.timestep, lft_in.horizon_period);
    
    % Formulate LMIs
    for j = 1:total_time
        lmi_x_mat = a{j} * x_big{j} * a{j}' - x_big_plus{j} + b{j} * b{j}';
        ct_x = ct_x + ((lmi_x_mat <= -lmi_shift * eye(size(x_big{j}))):...
                       ['X constraints, timestep: ' num2str(j)]);               %#ok<BDSCA>
        obj_x = obj_x + trace(x_big{j});
        
        lmi_y_mat = a{j}' * y_big_plus{j} * a{j} - y_big{j} + c{j}' * c{j};
        ct_y = ct_y + (( lmi_y_mat <= -lmi_shift * eye(size(y_big{j}))):...
                       ['Y constraints, timestep: ' num2str(j)]);               %#ok<BDSCA>
        obj_y = obj_y + trace(y_big{j});
    end
    
    % Solve for gramians
    yalmip_settings = sdpsettings('verbose', false);
    optimize(ct_x, obj_x, yalmip_settings);
    optimize(ct_y, obj_y, yalmip_settings);
    primal_x = check(ct_x);
    primal_y = check(ct_y);
    valid_x = all(~isnan(primal_x), 'all') &&  all(primal_x >= -lmi_shift);
    valid_y = all(~isnan(primal_y), 'all') &&  all(primal_y >= -lmi_shift);
    assert(valid_x && valid_y,...
           'modelReduction:computeGramian',...
           'Could not find admissible Gramians');
end

function [lft_bal, x_hat] = balance(lft_in, x, y)
%% BALANCE takes an LFT and gramians associated with that LFT (either of the LFT
%  directly, or of the LFT's RCF) and balances the LFT with those gramians.  It 
%  also outputs the balanced gramian associated with the given gramians.
    [a, b, c, d] = getABCD(lft_in);
    total_time = sum(lft_in.horizon_period);
    num_deltas = length(lft_in.delta.deltas);    
    
    % Calculate transformation matrices
    t = cell(num_deltas, total_time);
    t_inv = cell(num_deltas, total_time);
    for j = 1:total_time
        for i = 1:num_deltas
            x_chol = chol(value(x{i, j}));
            y_chol = chol(value(y{i, j}));
            [~, s, v] = svd(y_chol * x_chol');
            t_inv{i, j} = x_chol' * v * s^(-1/2);
            t{i, j} = inv(t_inv{i, j}); %s^(-1/2)*u*y_chol;
        end
    end    
    t_inv_big = cell(1, total_time);
    for j = 1:total_time
        t_inv_big{j} = blkdiag(t_inv{:, j});
    end
    t_big_plus = shiftForward(t, lft_in.timestep, lft_in.horizon_period);
    
    a_hat = cell(1, total_time);
    b_hat = cell(1, total_time);
    c_hat = cell(1, total_time); 
    x_hat = cell(num_deltas, total_time);
    for j = 1:total_time
        a_hat{j} = t_big_plus{j} * a{j} * t_inv_big{j};
        b_hat{j} = t_big_plus{j} * b{j};
        c_hat{j} = c{j} * t_inv_big{j};
        for i = 1:num_deltas
            x_hat{i, j} = t{i, j} * value(x{i, j}) * t{i, j}';
        end
    end
    lft_bal = Ulft(a_hat, b_hat, c_hat, d, lft_in.delta,...
                   'horizon_period', lft_in.horizon_period,...
                   'performance', lft_in.performance,...
                   'disturbance', lft_in.disturbance);
end

function [trunc_ind, error_out] = truncationIndices(lft_in, x_bal, desired_err)
%% TRUNCATIONINDICES takes an LFT, a balanced gramian, and a desired bound on
%  the error, and returns the indices of each Delta that should be truncated,
%  along with the derived upper bound on the truncation error.  If the LFT is stable
%  this bound applies to the error induced by truncating that LFT, if the LFT is
%  unstable, this bound applies to the error induced by truncating the RCF of 
%  that LFT.
    total_time = sum(lft_in.horizon_period);
    num_del = length(lft_in.delta.deltas);
    total_error = zeros(1, total_time);
    trunc_ind = cell(num_del, total_time);
    dim_outs = lft_in.delta.dim_outs;
    
    for j = 1:total_time
        for i = 1:num_del
            error_bound = 2 * diag(x_bal{i, j});
            if i > 1
                index_shift = 0;
                for k = 1:i - 1
                    index_shift = index_shift + dim_outs(k, j);
                end
            else
                index_shift = 0;
            end
            temp_index = find(error_bound < desired_err);
            trunc_ind{i, j} = index_shift + temp_index;
            total_error(j) = total_error(j) + sum(error_bound(temp_index));
        end
    end
    error_out = total_error;
end

function lft_red = truncate(lft_in, trunc_ind)
%% TRUNCATE takes an LFT and the desired truncation indices, and returns the
%  truncated LFT.
    [a, b, c, d] = getABCD(lft_in);
    total_time = sum(lft_in.horizon_period);
    num_del = length(lft_in.delta.deltas);
    
    for j = 1:total_time
        col_drop_ind = [];
        row_drop_ind = [];
        for k = 1:num_del
            col_drop_ind = [col_drop_ind; trunc_ind{k, j}];
            if ~isequal(j, total_time)
                row_drop_ind = [row_drop_ind; trunc_ind{k, j + 1}];
            else
                row_drop_ind = [row_drop_ind;
                                trunc_ind{k, lft_in.horizon_period(1) + 1}];
            end
        end

        a{j}(:, col_drop_ind) = [];
        a{j}(row_drop_ind, :) = [];
        b{j}(row_drop_ind, :) = [];
        c{j}(:, col_drop_ind) = [];

    end

    new_deltas = lft_in.delta.deltas;
    new_dim_outins = cellfun(@(dim_oi, inds) dim_oi - length(inds),...
                             num2cell(lft_in.delta.dim_outs),...
                             trunc_ind);
       
    dim_state = zeros(1, total_time);
    for j = 1:total_time
        dim_state(j) = new_deltas{1}.dim_out(j) - size(trunc_ind{1, j}, 1);
    end
    
    removeable_deltas = [];
    for i = 1:num_del
        if isa(new_deltas{i}, 'DeltaDelayZ') && ~all(dim_state == 0)
            new_deltas{i} = DeltaDelayZ(dim_state,...
                                    new_deltas{i}.timestep,...
                                    new_deltas{i}.horizon_period);
        elseif isa(new_deltas{i}, 'DeltaDelayZ') && all(dim_state == 0)
            removeable_deltas = [removeable_deltas; i];
        else
            if all(new_dim_outins(i, :) == 0 )
                removeable_deltas = [removeable_deltas; i];
            else
                new_deltas{i}.dim_in = new_dim_outins(i, :);
                new_deltas{i}.dim_out = new_dim_outins(i, :);
            end
        end
    end
    
    new_deltas(removeable_deltas) = [];
    
    new_deltas = SequenceDelta(new_deltas);
    lft_red = Ulft(a, b, c, d, new_deltas,...
                   'horizon_period', lft_in.horizon_period,...
                   'performance', lft_in.performance,...
                   'disturbance', lft_in.disturbance);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Second tier subfunctions %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [a, b, c, d] = getABCD(lft_in)
%% GETABCD helpful getter function for LFTs.
a = lft_in.a;
b = lft_in.b;
c = lft_in.c;
d = lft_in.d;
end

function p_plus = shiftForward(p, timestep, horizon_period)
%% SHIFTFORWARD helpful function to shift structured sequences of matrices
%  forward in time.  If the lft is discrete-time, the p matrix pertaining to
%  the DeltaDelayZ is shifted forward, while all the rest are not. If the lft
%  is not discrete, no matrices are shifted forward.
    total_time = sum(horizon_period);
    delay = ~isempty(timestep);
    p_plus = cell(1, size(p, 2));
    for j = 1:total_time
        if delay
            if j < total_time
                next_time = j + 1;
            else
                next_time = horizon_period(1) + 1;
            end
        else
            next_time = j;
        end
        p_plus{j} = blkdiag(p{1, next_time}, p{2:end, j});
    end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0): Added after v0.5.0 - Philip Mulford and Micah Fry ({philip.mulford, micah.fry}@ll.mit.edu)