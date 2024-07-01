function [curr_spks, store] = calc_binLen_selectSpikes(idcs, w, spikes_dataLen, nSpikes, store, halves, trialHalf)
% this function calculates the bin width around each bin centre to make
% sure that an equal number of spikes (nSpikes) is contained in each window


spk = [];
% take equally as many spikes to the right...
idcs1 = bsxfun(@plus, idcs(:,w), repmat(0:2000, size(idcs,1), 1));
spk(:,:,1) = spikes_dataLen(idcs1);
% ...and to the left of the current window center (but need to 
% flip it, so that the counting starts from the time of the event)
idcs2 = bsxfun(@plus, idcs(:,w), repmat(-2000:0, size(idcs,1), 1));

spk(:,:,2)   = fliplr(spikes_dataLen(idcs2));
sumSpk       = cumsum(sum(sum(spk, 3)));
winHalfWidth = find(sumSpk >= nSpikes, 1);  % find the very first index where the sum spikes is larger than our target number of spikes (nSpikes)

if iscell(halves)
    store.(halves{trialHalf}).winHalfWidths(w) = winHalfWidth;
else
    store.all.winHalfWidths(w) = winHalfWidth;
end

win = -winHalfWidth:winHalfWidth;
spkSearch = bsxfun(@plus, idcs(:,w), repmat(win, size(idcs,1), 1)); % get the indices of the windows that have the required length (to include the target nr of spikes)
curr_spks = nan(size(spikes_dataLen));                              % prepare a nan vector
curr_spks(spkSearch(:)) = spikes_dataLen(spkSearch(:));   % copy only the spike information from the selected windows to the nan vector. this format allows an easy extraction of the phase later on.

if iscell(halves)
    store.(halves{trialHalf}).nSpikes          = nSpikes;
    store.(halves{trialHalf}).nSpikes_win(w)   = sum(spikes_dataLen(spkSearch(:)));
else
    store.nSpikes        = nSpikes;
    store.nSpikes_win(w) = sum(spikes_dataLen(spkSearch(:)));
end
