function delta_out = toDelta(rct_obj, del_map)
%% TODELTA - Function for converting uncertain RCT objects to iqcToolbox Deltas.
%
%     delta_out = toDelta(rct_obj, del_map)
%
%     If input is ultidyn, output is a DeltaDlti.
%     If input is ureal, output is a DeltaSlti. 
%     If input is udyn, output is a Delta object as defined by 
%                                 the del_map(rct_obj.Name) 
%     Note: rct_obj.NominalValue is lost when converting to Delta objects.
% 
%     Variables:
%     ---------
%       Input:
%         rct_obj : ureal, ultidyn, or udyn (from Robust Control Toolbox)
%         del_map : map of delta names and objects from lftToRct, only necessary
%                    when rct_obj is a udyn
%       output:
%         delta_out : Delta object :: the resultant delta
%
%     See also rctToLft, lftToRct

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Calling appropriate functions for different input arguments

if isa(rct_obj, 'ureal')
    delta_out = urealToDelta(rct_obj);
elseif isa(rct_obj, 'ultidyn')
    delta_out = ultidynToDelta(rct_obj);
elseif isa(rct_obj, 'udyn')
    delta_out = udynToDelta(rct_obj, del_map);
    
else
    error('toDelta:toDelta', ['Can only convert udyn',...
                              ' (accompanied with a map to the pertinent Delta',...
                              ', ureal, or ultidyn objects.']);
end
end

function delta_out = ultidynToDelta(ulti_in)
    goodDelta = @(delta) isa(delta, 'ultidyn') &&...
                         isequal(delta.Type, 'GainBounded');
    assert(goodDelta(ulti_in),...
           'toDelta:toDelta',...
           strcat('Cannot make a DeltaDlti obect where uncertainties are not',...
                  ' GainBounded ultidyns'))
    delta_out = DeltaDlti(ulti_in.Name,...
                          size(ulti_in.OutputUnit,1),...
                          size(ulti_in.InputUnit,1),...
                          ulti_in.Bound);
end

function delta_out = urealToDelta(ureal_in)
    nom_check = @(delta) isequal(delta.NominalValue, 0);
    
    if ~nom_check(ureal_in)
       warning('toDelta:toDelta',...
                strcat('Nominal value information will not be retained',...
                       ' in the DeltaSlti'));
    end
    
    delta_out = DeltaSlti(ureal_in.Name,...
                          1,...
                          ureal_in.Range(1),...
                          ureal_in.Range(2));
end

function delta_out = udynToDelta(udyn_in, del_map)
    name = erase(udyn_in.Name, 'Normalized');
    assert(isKey(del_map, name), 'toDelta:toDelta',...
           ['The container "del_map" does not have a key-value pair for the',...
            ' given udyn uncertainty: ', name, '. Note that the name should',...
            ' not be preceded with the string "Normalized"'])
    delta_out = toLft(del_map(name)).delta.deltas{1};
end
%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Philip Mulford (philip.mulford@ll.mit.edu)