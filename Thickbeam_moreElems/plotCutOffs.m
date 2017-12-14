function [hnd] = plotCutOffs(x,y,COvector,colors,lines)
% Plots one figure with the y vector sectioned according to the COvector
% This code will work if the last 3 arguments are not given, but you cant
% skip arguments. So if you declare 'lines', COvector can be [], but colors
% cannot.
if (isempty(COvector))
    COvector = length(y);
end
n = length(COvector);
ini = 1;

if (nargin<4)
    for ii=1:n
        endi = COvector(ii);
        plot(x(ini:endi),y(ini:endi))
        hold on
        ini = COvector(ii);
    end
    plot(x(ini:end),y(ini:end))
else
    if nargin<5
        for ii=1:n
            endi = COvector(ii);
            plot(x(ini:endi),y(ini:endi),colors{ii})
            hold on
            ini = COvector(ii);
        end
        plot(x(ini:end),y(ini:end),colors{ii+1})
    else
        for ii=1:n
            endi = COvector(ii);
            plot(x(ini:endi),y(ini:endi),colors{ii},'LineWidth',lines(ii))
            hold on
            ini = COvector(ii);
        end
        plot(x(ini:end),y(ini:end),colors{ii+1},'LineWidth',lines(ii+1))
    end
end
% hold off
hnd = gcf;