function B = fpuField2XYZ(Bany)

if (~isstruct(Bany)) % suppose B3/B4
    if ndims(Bany) == 4
        B.x = Bany(:, :, :, 1);
        B.y = Bany(:, :, :, 2);
        B.z = Bany(:, :, :, 3);
    else
        B.x = Bany(:, 1);
        B.y = Bany(:, 2);
        B.z = Bany(:, 3);
    end
else
    if isfield(Bany, 'x')
        B = Bany;
    else
        if isfield(Bany, 'BX')
            B.x = Bany.BX;
            B.y = Bany.BY;
            B.z = Bany.BZ;
        else
            disp('Wrong Structure!')
            return
        end
    end
end

end
