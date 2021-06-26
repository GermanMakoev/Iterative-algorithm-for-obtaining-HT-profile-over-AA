function [Babs, Bincl, Bazim, transv, Bz] = fpuFieldVal(Bc)
% Bc - field in B.CC (Cartesian coordinates) form (B.x, B.y, B.z) or B4 form (index should exist or default 1)
% B - observed field in B.ME (Milne-Eddington, MERLIN) form (B.absB, B.incl, B.azim)
% NOTE: B.azim is true (-180...180)

Babs = [];
Bincl = [];
Bazim = [];
transv = [];
Bz = [];

Bc = fpuField2XYZ(Bc);

transv = sqrt(Bc.x.^2 + Bc.y.^2);
Bz = Bc.z;
Babs = sqrt(transv.^2 + Bc.z.^2);
Bincl = acos(Bc.z./Babs)*180/pi;
Bazim = atan2(Bc.x, Bc.y)*180/pi;

end
