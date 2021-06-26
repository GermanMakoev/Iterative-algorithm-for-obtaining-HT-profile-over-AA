function [hFig, hAxes] = ssuCreateMultiAxesFig(n1, n2, name, maximize, squared)

if (~exist('name', 'var'))
    name = 'Default Figure';
end
if (~exist('maximize', 'var'))
    maximize = true;
end
if (~exist('squared', 'var'))
    squared = false;
end

hFig = figure('Name', name);

Wx = 1;
Wy = 0.98;
dx = 0.01;
dy = 0.01;
wx = (Wx-dx)/n2 - dx;
wy = (Wy-dy)/n1 - dx;

for k1 = 1:n1
    for k2 = 1:n2
        hAxes(k1, k2) = axes('NextPlot', 'add');
        if squared
            set(hAxes(k1, k2), 'DataAspectRatio', [1 1 1]);
        end
        drawnow
    end
end

if (maximize)
    set(hFig, 'Units', 'normalized', 'Position', [0 0 1 1]);
end

for k1 = 1:n1
    for k2 = 1:n2
        set(hAxes(k1, k2), 'OuterPosition', [dx+(k2-1)*wx dy+(n1-k1)*wy wx wy]);
        drawnow
    end
end

end
