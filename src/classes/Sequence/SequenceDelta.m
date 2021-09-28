classdef SequenceDelta
%% SEQUENCEDELTA class for collecting different Delta objects. This class
% assists in conducting the same operation to all the Delta objects in its
% sequence
%
%   methods:
%     SequenceDelta(varargin) :: Constructor method
%     disp(this_deltas) :: Display method
%     matchHorizonPeriod(this_deltas, horizon_period) :: Matches deltas properties to new horizon_period
%
%   properties:
%     deltas : a cell array of Delta objects
%
%   dependent properties:
%     dim_outs : (num_del x total_time) array of delta output dimensions
%     dim_ins  : (num_del x total_time) array of delta input dimensions
%     horizon_periods : (num_del x 2) array of delta horizon_periods
%     names : (1 x num_del) cell array of delta names
%     types : (1 x num_del) cell array of each delta's class
%
%   See also SequenceDelta.SequenceDelta, SequenceDisturbance, SequencePerformance, Ulft.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    deltas
end

properties (Dependent)
    dim_outs
    dim_ins
    horizon_periods
    names
    types
end

methods
    function this_deltas = SequenceDelta(varargin)
    %% SEQUENCEDELTA constructor
    %
    %  this_deltas = SequenceDelta(delta1, delta2, delta3)
    %  this_deltas = SequenceDelta(cell_of_deltas)
    %  this_deltas = SequenceDelta(cell_of_deltas, delta_4, cell_of_deltas_2)
    %
    %  Variables:
    %  ---------
    %     Input:
    %       any combination of Delta objects and cells of Delta objects
    %     Output:
    %       this_deltas : SequenceDelta object :: collection of given Deltas
    %
    %  See also SequenceDelta.
    
        % Allocate empty cell of appropriate length
        len_all = sum(cellfun(@(arg) length(arg), varargin));
        deltas = cell(1, len_all);
        
        % Check input arguments and fill
        isDelta = @(del) isa(del, 'Delta');
        isCellOfDeltas = @(del) iscell(del) && all(cellfun(isDelta, del));
        ind_del = 1;
        for i = 1:nargin
            del = varargin{i};
            assert(isDelta(del) || isCellOfDeltas(del),...
                   'SequenceDelta:SequenceDelta',...
                   'Input objects must be Deltas and/or Cells of Deltas')
            len_del = length(del);
            if isDelta(del)
                deltas{ind_del} = del;
            else
                deltas(ind_del : ind_del + len_del - 1) = horzcat(del);
            end
            ind_del = ind_del + len_del;
        end
        
        notDiscreteAndContinuous = ...
            ~(any(cellfun(@isa, deltas, repmat({'DeltaDelayZ'}, size(deltas))))...
              && ...
              any(cellfun(@isa, deltas, repmat({'DeltaIntegrator'}, size(deltas)))));
        assert(notDiscreteAndContinuous,...
               'SequenceDelta:SequenceDelta',...
               'Cannot combine DeltaDelayZ and DeltaIntegrator in a Sequence')
        
        % Check if like-named deltas are indeed equivalent
        names = cellfun(@(del) del.name, deltas, 'UniformOutput', false);
        [~, unique_indices] = unique(names);
        for i = unique_indices'
            delta = deltas{i};
            delta_group = deltas(strcmp(names{i}, names));
            for j = 1:length(delta_group)
                % Remove checks on dim_in, dim_out
                warn_state = warning('query', 'MATLAB:structOnObject');
                warning('off', 'MATLAB:structOnObject');
                del_no_dim = rmfield(struct(delta), {'dim_in', 'dim_out'});
                del_grp_no_dim = rmfield(struct(delta_group{j}),...
                                         {'dim_in', 'dim_out'});
                warning(warn_state);
                % Check all other properties
                assert(isequaln(del_no_dim, del_grp_no_dim),...                        
                       'SequenceDelta:SequenceDelta',...
                       ['There cannot be multiple Deltas that share the',...
                        ' same name but have different properties'])
            end
        end
            
        horizon_periods = cell2mat(cellfun(@(del) del.horizon_period,...
                                           deltas',...
                                           'UniformOutput', false));
        horizon_period = commonHorizonPeriod(horizon_periods);
        deltas = cellfun(@(del) matchHorizonPeriod(del, horizon_period),...
                         deltas,...
                         'UniformOutput', false);
        
        % Assign to class property
        this_deltas.deltas = deltas;
    end
    
    function disp(this_deltas)
    %% DISP function for SequenceDelta object
    %
    %  disp(this_deltas) (e.g., disp(SequenceDelta(DeltaSlti('dis')) )
    %
    %  Variables:
    %  ---------
    %     Input:
    %       this_dis : SequenceDelta object  
    %
    %  See also Ulft.disp, Delta.disp
        
        fprintf([repmat('*', [1,80]), '\n',...
                 repmat('*', [1,34]),'   DELTAS   ', repmat('*', [1,34]), '\n']);
        cellfun(@(del) disp(del), this_deltas.deltas);
    end
    
    function this_deltas = matchHorizonPeriod(this_deltas, horizon_period)
    %% MATCHHORIZONPERIOD function to apply matchHorizonPeriod to all Deltas
    %  in SequenceDelta object
    %
    %  this_deltas = matchHorizonPeriod(this_deltas, horizon_period)
    %
    %  Variables:
    %  ---------
    %     Input:
    %       this_deltas : SequenceDelta object
    %       new_horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of new sequence      
    %     Output:
    %       this_deltas : SequenceDelta object
    %
    %  See also Ulft.matchHorizonPeriod
        this_deltas.deltas = ...
            cellfun(@(del, hp) matchHorizonPeriod(del, hp),...
                    this_deltas.deltas,...
                    repmat({horizon_period}, size(this_deltas.deltas)),...
                    'UniformOutput', false);
    end    
       
    %% Getter functions
    function dim_outs = get.dim_outs(this_deltas)
        dim_outs = cell2mat(cellfun(@(del) del.dim_out,...
                                    this_deltas.deltas',...
                                    'UniformOutput', false));
    end
    
    function dim_ins = get.dim_ins(this_deltas)
        dim_ins = cell2mat(cellfun(@(del) del.dim_in,...
                                    this_deltas.deltas',...
                                    'UniformOutput', false));
    end
    
    function names = get.names(this_deltas)
        names = cellfun(@(del) del.name,...
                        this_deltas.deltas,...
                        'UniformOutput', false);
    end
    
    function types = get.types(this_deltas)
        types = cellfun(@(del) class(del),...
                        this_deltas.deltas,...
                        'UniformOutput', false);
    end
    
    function horizon_period = get.horizon_periods(this_deltas)
        horizon_period = cell2mat(cellfun(@(del) del.horizon_period,...
                                          this_deltas.deltas',...
                                          'UniformOutput', false));
    end

end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)