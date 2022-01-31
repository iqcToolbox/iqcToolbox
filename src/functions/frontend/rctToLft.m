function lft_out = rctToLft(rct_obj, varargin) 
%% RCTTOLFT for converting ureal, ultidyn, umat, uss objects to Ulft objects.
%
%     lft_out = rctToLft(rct_obj)
%     lft_out = rctToLft(ureal)
%     lft_out = rctToLft(ultidyn)
%     lft_out = rctToLft(umat)
%     lft_out = rctToLft(uss(a,b,c,d))
%
%     Variables:
%     ---------
%       Input:
%         rct_obj : ureal, ultidyn, umat, or uss object
%         varagin : optional containers.Map object of delta names to rct
%                   deltas. 
%       Output:
%         lft_out : Ulft object :: the resultant lft
%
%     See also Ulft, Ulft.Ulft, toLft, lftToRct

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Calling appropriate functions for different input arguments
    if isa(rct_obj, 'uss')
        lft_out = ussToLft(rct_obj, varargin);
    elseif isa(rct_obj, 'umat')
        lft_out = umatToLft(rct_obj);
    elseif isa(rct_obj, 'ultidyn')
        lft_out = ultidynToLft(rct_obj);
    elseif isa(rct_obj, 'ureal')
        lft_out = urealToLft(rct_obj);
    else
        error('rctToLft:rctToLft', ['rctToLft only accepts uss, umat, ureal, or '...
              'ultidyn objects']);
    end
end
function lft_out = umatToLft(umat_in)

    % umats are only functions of ureal/ucomplex/ucomplexm objects. As we
    % do not currently support complex Deltas in iqcToolbox, this function
    % only deals with converting umats made of ureal uncertainties
    [m, delta_std, ~, delta_norm] = lftdata(umat_in);    
    delta_std = struct2cell(delta_std.Uncertainty);
    %Adding warning that we only use normalized uncertainty
    for i = 1:size(delta_norm,1)                         
        if ~isequal(delta_std{i}.Range(1), -1) ||...
           ~isequal(delta_std{i}.Range(2), 1) ||...
           ~isequal(delta_std{i}.NominalValue, 0)
            warning('rctToLft:rctToLft',...
                    strcat('toLft will use a normalized delta to convert to',...
                           ' a Ulft object. Be aware that the Delta bounds of ',... 
                           ' the Ulft object may not match the bounds of the',...
                           ' uss object'));
            break
        end
    end
    
    % Building the representative deltas for ureal, ultidyn
    delta_cell = cell(size(delta_norm, 1), 1);
    for i = 1:size(delta_norm)
        delta_cell{i} = toDelta(delta_norm{i});
    end
    
    %Finding and removing empty cells that represent repeated deltas
    delta_object = SequenceDelta(delta_cell);
    dim_in = sum(delta_object.dim_ins);
    dim_out = sum(delta_object.dim_outs);
    a = m(1:dim_in, 1:dim_out);
    b = m(1:dim_in, dim_out+1:end);
    c = m(dim_in+1:end, 1:dim_out);
    d = m(dim_in+1:end, dim_out+1:end);
    
    lft_out = Ulft(a, b, c, d, delta_object);
end

function lft_out = ussToLft(uss_in, varargin)

    [m, delta_std, ~, delta_norm] = lftdata(uss_in);
    m = [m.A, m.B; m.C, m.D];
    dim_state = size(uss_in.A, 1);
    % Checking if the system has a state - if not, time_delta is a flag set
    % to -1
    if dim_state <= 0   
        time_delta = -1;
    elseif isequal(uss_in.Ts, 0)
        time_delta = DeltaIntegrator(dim_state);
    elseif isequal(uss_in.Ts, -1)
        time_delta = DeltaDelayZ(dim_state);
    else
        time_delta = DeltaDelayZ(dim_state, uss_in.Ts);
    end
    
    uncertainties = struct2cell(uss_in.Uncertainty);
    
    for i = 1:size(uncertainties,1)
        std_uncertainty = struct2cell(delta_std.Uncertainty);
        std_uncertainty = std_uncertainty{i};
        type = class(std_uncertainty);
        
        % generate warnings that normalized values will be used instead of
        % the given nominal values
        switch type
            case 'ultidyn'
                if ~isequal(std_uncertainty.Bound, 1)
                    warning('rctToLft:rctToLft',...
                            strcat('rctToLft will used a normalized delta to',...
                                   ' convert to a Ulft object. Be aware that',...
                                   ' the upper bound of the DeltaDlti may',...
                                   ' not match the bounds of the ultidyn',...
                                   ' uncertainty'))
                end
            otherwise
            if ~isequal(std_uncertainty.Range(1), -1) ||...
               ~isequal(std_uncertainty.Range(2), 1) ||...
               ~isequal(std_uncertainty.NominalValue, 0)
                warning('rctToLft:rctToLft',...
                        strcat('rctToLft will use a normalized delta to',...
                               ' convert to a Ulft object. Be aware that',...
                               ' the Delta bounds of the DeltaSlti may',...
                               ' not match the bounds of the ureal uncertainty'));
            end    
        end
    end
    
    delta_cell = cell(size(delta_norm,1), 1);
    map = varargin{1}{1};
    if nargin == 1
        for i = 1:size(delta_norm)
            delta_cell{i} = toDelta(delta_norm{i});
        end
    elseif nargin == 2
        for i = 1:size(delta_norm)
            delta_cell{i} = toDelta(delta_norm{i}, map);
        end
    end
    
    % check to see if there are any states/whether any time deltas were
    % generated
    if ~isequal(time_delta, -1)
        delta_object = SequenceDelta(time_delta, delta_cell);
    else
        delta_object = SequenceDelta(delta_cell);
    end
    
    dim_out =  sum(delta_object.dim_outs);
    dim_in = sum(delta_object.dim_ins);
    a = m(1:dim_in, 1:dim_out);
    b = m(1:dim_in, dim_out + 1:end);
    c = m(dim_in + 1:end, 1:dim_out);
    d = m(dim_in + 1:end, dim_out + 1:end);
    
    lft_out = Ulft(a, b, c, d, delta_object);
    
end

function lft_out = ultidynToLft(ultidyn_in)  
    lft_out = toLft(toDelta(ultidyn_in));
end

function lft_out = urealToLft(ureal_in)
    lft_out = toLft(toDelta(ureal_in));
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Philip Mulford (philip.mulford@ll.mit.edu)