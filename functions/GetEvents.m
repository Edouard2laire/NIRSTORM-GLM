function [events, missing_events] = GetEvents(sEvents, eventNames)
% Return the events based on their name
    event_idx = zeros(1, length(eventNames));
    missing_events = {};
    for iEvt = 1:length(eventNames)
        idx = find(strcmp({sEvents.label}, eventNames{iEvt}));
        if isempty(idx)
            missing_events{end+1} =  eventNames{iEvt};
            continue;
        end
        event_idx(iEvt) = idx;
    end
    event_idx = event_idx(event_idx > 0);
    events = sEvents(event_idx);
end