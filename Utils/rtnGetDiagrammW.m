function dwidth = rtnGetDiagrammW(freq, mode, c, b)

if mode == 3
    dwidth = 0.2  + 9.4  *3e10./freq;
elseif mode == 2
    dwidth = 4.38 + 6.87 *3e10./freq;
elseif mode == 1
    dwidth = 8.5         *3e10./freq;
else
    dwidth = c    + b    *3e10./freq;
end

end
