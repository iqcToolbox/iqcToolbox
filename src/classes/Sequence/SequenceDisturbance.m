classdef SequenceDisturbance
%% SEQUENCEDISTURBANCE class for collecting different Disturbance objects.
% This class assists in conducting the same operation to all the 
% Disturbance objects in its sequence
%
%   methods:
%     SequenceDisturbance(varargin) :: Constructor method
%     disp(this_disturbances) :: Display method
%     matchHorizonPeriod(this_disturbances, horizon_period) :: Matches disturbances properties to new horizon_period
%
%   properties:
%     disturbances : a cell array of Disturbance objects
%
%   dependent properties:
%     channels : (num_del x 1) cell array of each disturbances channel property
%     horizon_periods : (num_del x 2) array of disturbance horizon_periods
%     names : (1 x num_del) cell array of disturbance names
%     types : (1 x num_del) cell array of each disturbance's class
%
%   See also SequenceDisturbance.SequenceDisturbance, SequenceDelta, SequencePerformance, Ulft.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    disturbances (1, :) cell
end

properties (Dependent)
    chan_ins
    horizon_periods
    names
    types
end

methods
    function this_disturbances = SequenceDisturbance(varargin)
    %% SEQUENCEDISTURBANCE constructor
    %
    %  this_disturbances = SequenceDisturbance(dis1, dis2, dis3)
    %  this_disturbances = SequenceDisturbance(cell_of_diss)
    %  this_disturbances = SequenceDisturbance(cell_of_diss, dis_4, cell_of_diss_2)
    %
    %  Variables:
    %  ---------
    %     Input:
    %       any combination of Disturbance objects and cells of Disturbance objects
    %     Output:
    %       this_disturbances : SequenceDisturbance object :: collection of given Disturbances
    %
    %  See also SequenceDisturbance.
    
        % Allocate empty cell of appropriate length
        len_all = sum(cellfun(@(arg) length(arg), varargin));
        disturbances = cell(1, len_all);
        
        % Check input arguments and fill
        isDisturbance = @(dis) isa(dis, 'Disturbance');
        isCellOfDisturbances = @(dis) iscell(dis) ...
                                       && ...
                                       all(cellfun(isDisturbance, dis));
        ind_dis = 1;
        for i = 1:nargin
            dis = varargin{i};
            assert(isDisturbance(dis) || isCellOfDisturbances(dis),...
                   'SequenceDisturbance:SequenceDisturbance',...
                   ['Input objects must be Disturbances',...
                    'and/or Cells of Disturbances'])
            len_dis = length(dis);
            if isDisturbance(dis)
                disturbances{ind_dis} = dis;
            else
                disturbances(ind_dis : ind_dis + len_dis - 1) = horzcat(dis);
            end
            ind_dis = ind_dis + len_dis;
        end
        
        % Check if like-named disturbances are indeed equivalent
        names = cellfun(@(dis) dis.name, disturbances, 'UniformOutput', false);
        [~, unique_indices] = unique(names);
        for i = unique_indices'
            dis = disturbances{i};
            dis_group = disturbances(strcmp(names{i}, names));
            for j = 1:length(dis_group)
                % Remove checks on channel
                warn_state = warning('query', 'MATLAB:structOnObject');
                warning('off', 'MATLAB:structOnObject');
                dis_no_chan = rmfield(struct(dis), {'chan_in'});
                dis_grp_no_dim = rmfield(struct(dis_group{j}),...
                                         {'chan_in'});
                warning(warn_state);
                % Check all other properties
                assert(isequaln(dis_no_chan, dis_grp_no_dim),...                        
                       'SequenceDisturbance:SequenceDisturbance',...
                       ['There cannot be multiple Disturbances that share',...
                        'the same name but have different properties'])
            end
        end
        
        horizon_periods = cell2mat(cellfun(@(dis) dis.horizon_period,...
                                           disturbances',...
                                           'UniformOutput', false));
        horizon_period = commonHorizonPeriod(horizon_periods);
        disturbances = cellfun(@(dis) matchHorizonPeriod(dis, horizon_period),...
                               disturbances,...
                               'UniformOutput', false);
        
        % Assign to class property
        this_disturbances.disturbances = disturbances;
    end     
    
    function disp(this_disturbances)
    %% DISP function for SequenceDisturbance object
    %
    %  disp(this_disturbances) (e.g., disp(SequenceDisturbance(DisturbanceL2('l2')) )
    %
    %  Variables:
    %  ---------
    %     Input:
    %       this_dis : SequenceDisturbance object   
    %
    %  See also Ulft.disp, Disturbance.disp
    
        fprintf([repmat('*', [1,80]), '\n',...
                 repmat('*', [1,31]),'   DISTURBANCES   ', repmat('*', [1,31]), '\n']);
        cellfun(@(dis) disp(dis), this_disturbances.disturbances);
    end
    
    function this_dis = matchHorizonPeriod(this_dis, horizon_period)
    %% MATCHHORIZONPERIOD function to apply matchHorizonPeriod to all Disturbances
    %  in SequenceDisturbance object
    %
    %  this_disturbances = matchHorizonPeriod(this_disturbances, horizon_period)
    %
    %  Variables:
    %  ---------
    %     Input:
    %       this_disturbances : SequenceDisturbance object
    %       new_horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of new sequence   
    %     Output:
    %       this_disturbances : SequenceDisturbance object
    %
    %  See also Ulft.matchHorizonPeriod
    
        this_dis.disturbances = ...
            cellfun(@(dis, hp) matchHorizonPeriod(dis, hp),...
                    this_dis.disturbances,...
                    repmat({horizon_period}, size(this_dis.disturbances)),...
                    'UniformOutput', false);
    end    
    
    %% Getter functions
    function chan_ins = get.chan_ins(this_disturbances)
        channel_cell = cellfun(@(dis) dis.chan_in,...
                               this_disturbances.disturbances',...
                               'UniformOutput', false);
        chan_ins = vertcat(channel_cell{:});
%         channels = cell2mat(cellfun(@(dis) dis.chan_in,...
%                                     this_disturbances.disturbances',...
%                                     'UniformOutput', false));
    end
       
    function names = get.names(this_disturbances)
        names = cellfun(@(dis) dis.name,...
                        this_disturbances.disturbances,...
                        'UniformOutput', false);
    end
    
    function types = get.types(this_disturbances)
        types = cellfun(@(dis) class(dis),...
                        this_disturbances.disturbances,...
                        'UniformOutput', false);
    end
    
    function horizon_period = get.horizon_periods(this_disturbances)
        horizon_period = cell2mat(cellfun(@(dis) dis.horizon_period,...
                                          this_disturbances.disturbances',...
                                          'UniformOutput', false));
    end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)