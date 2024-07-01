function x_ci = mybootci(y, nbtsp,varargin); 

if nargin > 2
    func = varargin{1}; % function handle
else
    func = @(x) mean(x,'omitnan');
end
x_btsp = bootstrp(nbtsp,func,y);
x_ci = [prctile(x_btsp,5); func(y); prctile(x_btsp,95)];

