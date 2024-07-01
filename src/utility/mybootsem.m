function x_sem = mybootsem(y, nbtsp,varargin); 

if nargin > 2
    func = varargin{1}; % function handle
else
    func = @(x) mean(x,'omitnan');
end
x_btsp = bootstrp(nbtsp,func,y);
x_sem =  func(y) + [-std(x_btsp,'omitnan'); zeros(size(func(y))); std(x_btsp,'omitnan')];

