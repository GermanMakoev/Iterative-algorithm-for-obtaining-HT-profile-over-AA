function ssuAxImSurf(hAx, imdata, scale, offset, cmap, colbar)

if ~exist('scale', 'var') || isempty(scale)
    scale = [1 1];
end
if ~exist('offset', 'var') || isempty(offset)
    offset = [0 0];
end
if ~exist('colbar', 'var')
    colbar = true;
end

mainMapSize = 231;
centered = false;
nandata = mainMapSize;
if ~exist('cmap', 'var') || isempty(cmap)
    cmap = ssuSymmColormap(mainMapSize);
    centered = true;
    nandata = 130;
end

imdata(isnan(imdata)) = nandata;

absmax = max(max(abs(imdata)));
if (centered)
    imdata(1, 1) = -absmax;
    imdata(end, end) = absmax;
end

xrange = (0:size(imdata, 2)-1)*scale(2) + offset(2);
yrange = (0:size(imdata, 1)-1)*scale(1) + offset(1);
surf(hAx, xrange, yrange, zeros(size(imdata)), imdata, 'LineStyle','none');
set(hAx,'YDir','normal', 'Layer', 'top');
set(hAx,'XLimMode','manual', 'XLim', [min(xrange) max(xrange)]);
set(hAx,'YLimMode','manual', 'YLim', [min(yrange) max(yrange)]);
view(0, 90);

set(hAx,'DataAspectRatio', [1,1,1]);
colormap(hAx, cmap)
if colbar
    colorbar(hAx)
end

end
