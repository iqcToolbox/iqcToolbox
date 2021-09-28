classdef (Abstract) Delta 
%% DELTA abstract base class for uncertainties
% This class should not be directly accessed.  Subclasses of Delta (such as
% DeltaSlti, DeltaBounded, etc.) extend and concretize this class.
%
%   concrete methods:
%     Delta(name, dim_out, dim_in) :: Constructor method
%     disp(this_delta) :: Display method
%     recastMatricesAndDelta(this_delta) :: Method for changing initial LFT to the analyzed one
%     normalizeDelta(this_delta) :: Method for changing unnormalized LFT to normalized one
%     toRct(this_delta) :: Method for converting Delta object into robust control toolbox objects
%     sample(this_delta, discrete) :: Randomly generates a valid realization of this delta given whether this LFT is discrete-time or not
%     validateSample(this_delta, value, discrete) :: Throws an error if the given value is not a valid realization of this delta for the given time-type
%     
%   abstract methods:
%     matchHorizonPeriod(this_delta, new_horizon_period) 
%                        :: Matches delta properties to new horizon_period
%     deltaToMultiplier(this_delta, varargin)
%                        :: Method for constructing a multiplier from a delta
%
%   properties:
%     name : char array :: unique ID of the uncertainty (ex. 'dmass' or 'dX')
%     dim_out : natural :: output dimensions of uncertainty
%     dim_in : natural :: input dimensions of uncertainty
%     horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) 
%                                                 of Delta properties
%
%   See also Delta.Delta.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    name
    dim_out
    dim_in
    horizon_period
end

methods (Abstract)
    this_delta = matchHorizonPeriod(this_delta, total_time)
    multiplier = deltaToMultiplier(this_delta, varargin)
end

methods 
    function this_delta = Delta(name, dim_out, dim_in, horizon_period)
    %% DELTA constructor, must be called from a subclass constructor
    %
    %     this_delta = Delta(name, dim_out, dim_in, horizon_period)
    %
    %     Variables:
    %     ---------
    %       Input:
    %          name : string :: unique ID of the uncertainty (ex. 'dmass')
    %          dim_in : natural :: input dimensions of uncertainty
    %          dim_out : natural :: output dimensions of uncertainty
    %          horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) 
    %                                                      of Delta properties
    %       Output:
    %          this_delta : Delta object :: the base class Delta
    %
    %    See also Delta.
    validateattributes(name, 'char', {'nonempty'});
    notAllZeroDim = @(dim) any(dim > 0);
    validateattributes(dim_out, 'numeric', {'integer', 'nonnegative'});
    assert(notAllZeroDim(dim_out),...
           'Delta:Delta',...
           'Delta must have at least one non-zero output dimension')
    validateattributes(dim_in, 'numeric', {'integer', 'nonnegative'});
    assert(notAllZeroDim(dim_in),...
           'Delta:Delta',...
           'Delta must have at least one non-zero input dimension')
    validateattributes(horizon_period, 'numeric', {'size', [1,2],...
                                                   'integer',...
                                                   'nonnegative'});
    validateattributes(horizon_period(2), 'numeric', {'positive'})                                               

    % Setting properties of Delta
    this_delta.name           = name;
    this_delta.dim_out        = dim_out;
    this_delta.dim_in         = dim_in;
    this_delta.horizon_period = horizon_period;
    end

    function disp(this_delta, type)
    %% DISP function for Delta object, must be called from a subclass disp method.
    %
    %     disp(this_delta, type)
    %
    %     Variables:
    %     ---------
    %       Input:
    %          this_delta : Delta object 
    %          type : char array :: type of uncertainty ('SLTI', 'Delay', etc.')
    %
    %     See also Ulft.disp, SequenceDelta.disp.
       fprintf('%4s %-7s is a %3d x %-3d %18s \n',...
               '',...
               this_delta.name,...
               this_delta.dim_out(1),...
               this_delta.dim_in(1),...
               type)                      
    end
    
    function [recastA, recastB, recastC, recastD, recastDelta] = recastMatricesAndDelta(this_delta) %#ok<MANU>
    %% RECASTMATRICESANDDELTA method for creating a modified LFT for IQC analysis.
    %  this method should be extended for any subclass of Delta whereby IQC
    %  analysis is conducted on an analyzable, but different LFT (see, for
    %  example, DeltaSltvRateBnd and DeltaSltvRateBndImpl). The output
    %  arguments are function handles for modifying the initial LFT a, b,
    %  c, d matrices and Delta objects. These handles are used in the
    %  sub-function iqcAnalysis/modifyLft.
    %
    %    [newA, newB, newC, newD, newDelta] = recastMatricesAndDelta(this_delta)
    %
    %    Variables:
    %    ---------
    %      Input:
    %         this_delta : Delta object
    %      Output:
    %         recastA : function_handle :: function to transform a matrices of LFT
    %         recastB : function_handle :: function to transform b matrices of LFT
    %         recastC : function_handle :: function to transform c matrices of LFT
    %         recastD : function_handle :: function to transform d matrices of LFT
    %         recastDelta : Delta object :: new Delta object for modified LFT
    %
    %    See also iqcAnalysis.modifyLft, DeltaSltvRateBnd.recastMatricesAndDelta.
        recastA     = function_handle.empty();
        recastB     = function_handle.empty();
        recastC     = function_handle.empty();
        recastD     = function_handle.empty();
        recastDelta = [];
    end
    
    function [del_diff, del_ave, del_scale, del_norm] = normalizeDelta(this_delta)
    %% NORMALIZEDELTA function for Delta object. This method should be extended
    %  for any subclass of Delta that intends to conduct the normalization procedure.
    %  The method for the Delta superclass returns operations that, with the 
    %  Ulft.normalizeLft method, would not result in a normalized LFT. I.e., without
    %  extension of this method in Delta subclasses, the method in the Delta
    %  superclass defines the default behavior as "no normalization".
    %  If the superclass Delta method is called, it returns del_norm as empty,
    %  to indicate that the default method (i.e., no normalization) method was
    %  called.
    %
    %    [del_diff, del_ave, del_scale, del_norm] = normalizeDelta(this_delta)
    %
    %    Variables:
    %    ---------
    %      Input:
    %         this_delta : Delta object
    %      Output:
    %         del_diff : double matrix :: Matrix for modifying portions of LFT a, b, c, d
    %                                       matrices pertaining to this Delta
    %         del_ave : double matrix :: Matrix for modifying portions of LFT a, b, c, d
    %                                       matrices pertaining to this Delta
    %         del_scale : double matrix :: Matrix for modifying portions of LFT a, b, c, d
    %                                       matrices pertaining to this Delta
    %         del_norm : Delta object (or empty) :: Normalized Delta object. If empty, indicates
    %                                      that no normalization took place and a warning is thrown.
    %
    %    See also Ulft.normalizeLft, DeltaDlti.normalizeDelta, DeltaSltvRateBnd.normalizeDelta.
        total_time = sum(this_delta.horizon_period);
        del_diff   = cell(1, total_time);
        del_ave    = cell(1, total_time);
        del_scale  = cell(1, total_time);
        for i = 1:total_time
            del_diff{i}  = eye(this_delta.dim_out(i));
            del_ave{i}   = zeros(this_delta.dim_out(i), this_delta.dim_in(i));
            del_scale{i} = eye(this_delta.dim_in(i));
        end
        del_norm = [];
    end
    function rct_obj = toRct(this_delta)
    %% TORCT function for Delta object. This method is extended for specific 
    %  subclasses of Delta (such as DeltaDlti and DeltaSlti), which can be
    %  converted to Robust Control Toolbox objects.  Most Delta subclasses
    %  do not have a direct conversion, so the default conversion is to a
    %  udyn object with matching dimensions. Delta subclasses that do have
    %  direct conversion to RCT objects may implement those by extending
    %  this method.
    %    
    %    rct_obj = toRct(this_delta)
    % 
    %    Variables:
    %    ---------
    %      Input:
    %        this_delta : Delta object to be converted
    %      Output:
    %        rct_obj : udyn object (from the Robust Control Toolbox)
        rct_obj = udyn(this_delta.name, [this_delta.dim_out,this_delta.dim_in]);
    end

    function value = sample(this_delta, ~)                                      %#ok<STOUT>
        %% Randomly generates a valid realization of this delta.
        error('Delta:sample',...
              ['The class for "',this_delta.name,'" (',class(this_delta),...
               ') has no sampling method implemented yet.']);
    end

    function validateSample(this_delta, value, ~, validate_base)
        %% Throws an error if the given value is not a valid realization of this delta.
        if ~exist('validate_base', 'var') || ~validate_base
            error('Delta:validateSample',...
                  ['The class for"',this_delta.name,'" (',class(this_delta),...
                   ') has no sample validation method implemented yet.']);
        else
            % Check uncertainties
            assert(~value.uncertain,...
                   [class(this_delta),':validateSample'],...
                   ['Specified value for delta "',this_delta.name,...
                    '" must have no uncertainties.']);
            % Check horizon period
            assert(all(value.horizon_period == this_delta.horizon_period),...
                   [class(this_delta),':validateSample'],...
                   ['Specified value has invalid horizon_period for delta "',...
                    this_delta.name,'".']);
            % Check dimensions
            assert(isequal(size(value, 1), this_delta.dim_out) &&...
                   isequal(size(value, 2), this_delta.dim_in),...
                   [class(this_delta),':validateSample'],...
                   ['Specified value has invalid dimensions for delta "',...
                   this_delta.name,'".']);
        end
    end
end

methods (Sealed)
%% These methods are to allow the addition/multiplication/etc. of Delta objects 
%  with Ulft objects.
    function out = plus(varargin)
        del_ind =  find(cellfun(@(a) isa(a, 'Delta'), varargin), 1, 'first');
        varargin{del_ind} = toLft(varargin{del_ind});
        out = plus(varargin{:});
    end

    function out = mtimes(varargin)
        del_ind =  find(cellfun(@(a) isa(a, 'Delta'), varargin), 1, 'first');
        varargin{del_ind} = toLft(varargin{del_ind});
        out = mtimes(varargin{:});
    end
    
    function out = minus(varargin)
        del_ind =  find(cellfun(@(a) isa(a, 'Delta'), varargin), 1, 'first');
        varargin{del_ind} = toLft(varargin{del_ind});
        out = minus(varargin{:});
    end
    
    function out = vertcat(varargin)
        del_ind =  find(cellfun(@(a) isa(a, 'Delta'), varargin), 1, 'first');
        varargin{del_ind} = toLft(varargin{del_ind});
        out = vertcat(varargin{:});
    end
    
    function out = horzcat(varargin)
        del_ind =  find(cellfun(@(a) isa(a, 'Delta'), varargin), 1, 'first');
        varargin{del_ind} = toLft(varargin{del_ind});
        out = horzcat(varargin{:});
    end
    
    function out = blkdiag(varargin)
        del_ind =  find(cellfun(@(a) isa(a, 'Delta'), varargin), 1, 'first');
        varargin{del_ind} = toLft(varargin{del_ind});
        out = blkdiag(varargin{:});
    end
end


end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0): Added sample and validateSample methods - Jason Nezvadovitz (jason.nezvadovitz@ll.mit.edu)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry and Philip Mulford ({micah.fry, philip.mulford}@ll.mit.edu)