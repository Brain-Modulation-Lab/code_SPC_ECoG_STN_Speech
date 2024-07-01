function  [trialsMatrix, trialsVect] = getTrialsMatrix(events, data, dtTime, tWidth, tOffset)
% GETTRIALSMATRIX splits continuous data into epochs around events with a width 
% tWidth and an offset tOffset (positive numbers shift the epoch to start
% earlier than the event, negative numbers shift it to start after the
% event by the set offset), works both for time-frequency and time-data only
% OUTPUTS 
% trialsMatrix - matrix with epochs (time*freq*trials) if data is a
%               (time*freq) matrix, or (time*trials) if data is a time-vector      

events     = events(~isnan(events)); % discard any nans
nEvents    = length(events);

nWidth  = round(tWidth/dtTime);  % convert tWidth from seconds in samples    
nOffset = round(tOffset/dtTime); 

trialsMatrix = nan(size(data,1), nWidth + 1, nEvents);
trialsVect   = nan(size(data));

for i = 1:nEvents
    if isnan(events(i)), continue; end
    currtime = round(events(i)/dtTime);
    n1 = currtime - nOffset;
    n2 = n1 + nWidth;         
    trialsMatrix(:, :, i) = data(:, n1:n2);
    if size(data,1) == 1  % only compute this if you receive a vector instead
        % of a time-frequency matrix
        trialsVect(n1:n2) = data(:, n1:n2);
    end
end

trialsMatrix = squeeze(trialsMatrix);