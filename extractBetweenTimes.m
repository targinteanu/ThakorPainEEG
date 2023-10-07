function outEEG = extractBetweenTimes(inEEG, bound)
% Extract EEG between the time bounds specified as a new EEG object. 
% 
% Inputs: 
%   inEEG: EEG object input (must contain times in bounds +- 0.5 second buffer
%   bound: 2-element vector as [start time, end time] - SECONDS
% Outputs: 
%   outEEG: EEG object output only between bound times 

    buf = .5; % s
    bound(1) = max(inEEG.xmin, bound(1)-buf);
    bound(2) = min(inEEG.xmax, bound(2)+buf);
    outEEG = pop_select(inEEG, 'time', bound);
end