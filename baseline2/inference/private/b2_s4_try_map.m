function [map, exception] = b2_s4_try_map(channel, fisher, prior, initial, active, cfg)
% Only the original solver's explicit geometry-domain failure is recoverable.
% Invalid input, indefinite Fisher and unknown program/toolbox errors propagate.
exception=struct('identifier',"",'message',"");
try
    map=refine_position_map(channel,fisher,prior,initial,active,cfg);
catch caught
    if ~strcmp(caught.identifier,'refine_position_map:PositionAtRisCenter')
        rethrow(caught);
    end
    exception.identifier=string(caught.identifier);
    exception.message=string(caught.message);
    map=struct('positionEstimate',nan(6,1),'staticPositionEstimate',nan(3,1), ...
        'mobilePositionEstimate',nan(3,1),'objective',nan,'exitflag',nan,'iterations',nan);
end
end
