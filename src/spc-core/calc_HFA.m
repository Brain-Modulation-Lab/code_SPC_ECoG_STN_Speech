function HFA = calc_HFA(data, FS, cfg)
% calculate HFA (Fischer et al 2020 Elife)

Fhp = cfg.Fhp;
Flp = cfg.Flp;
ordFhp = cfg.ordFhp;
ordFlp = cfg.ordFlp;

FN = FS/2;

% 2. High pass filter at F_hp (e.g., 300 Hz)
[bHp, aHp] = butter(ordFhp, Fhp/FN ,'high');
data = filtfilt(bHp, aHp, data);

% 3. Full-wave rectification 
data = abs(data);

% 4. Low-pass filter at F_lp (e.g., 100 Hz)
[bLp, aLp] = butter(ordFlp, Flp/FN ,'low');
HFA = filtfilt(bLp, aLp, data);