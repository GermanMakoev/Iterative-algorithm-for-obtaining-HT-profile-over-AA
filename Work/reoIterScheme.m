% �������� ���������� ���� (������� mfoLoadField �������� ����� ������������ iouLoadField)
mfoData = mfoLoadField('C:\Users\Alan Makoev\Desktop\Matlab � 2\SDO\12470_hmi.M_720s.20151218_082209.W85N13CR.CEA.NAS_1000_sst.sav'); % for az=+12, +10
if isempty(mfoData)
    return
end

% �������� ������ ����� (��. �������� ��������� ������). ��� ������� �����
% �� ���� ������ � �������� ������ � ��������� ratan
ratan = iouLoadRATANData('C:\Users\Alan Makoev\Desktop\Matlab � 2\RATAN\RATAN_AR12470_20151218_082624_az+10_SPECTRA__flocculae-included_stille_appr.dat');

% � ������ ����� ������� � ������ ������ �����, ������� ����� (5) �������������
% �������� ��������� ����, ����� ��
select = 5;

% �������� ������� � �������
freqs = ratan.freqs;
Robs = ratan.right(select, :);
Lobs = ratan.left(select, :);

% ����� �������� � � ��������� �������, ��, ����� ������ ���, �����
% ����������� ���������������� ������� ���������� ��������.
Robsappr = asmAsym2SigOpt(ratan.freqs*1e-9, Robs);
Lobsappr = asmAsym2SigOpt(ratan.freqs*1e-9, Lobs);
% ��� ����������� � ���, ��� ������������� ���������, ����� ���������,
% ������, Robs � Robsappr �� ����� �������

% ����������� ���� ����� ��� ���������� � ������, ������� ������
% ������������� ��� � mfoData
mfoData.posangle = ratan.header.RATAN_P;

% ����� (���������) ������� �������� � ������������� ����������, � ���������, ��������
% ����:
QT = 1;
freefree = 0;

hLib = reoInitLibrary(QT, freefree);
if hLib == 0
    return
end

[M, base, Bc] = reoSetField(hLib, mfoData, step);

% ������ ����� �������� ������� ����� ����� (� ���. ���.), ������� �������
% ����� ����� �� ����������� � �������� ������������ ���� ����� ��������
% ��������� �������:
pos = floor((ratan.pos(select) - base(2))/step(1)) + 1;

% ��� ��� ����� ���� ����� ������ �����, ������ ���� �� ��������

% � ������ SDO ������ ������� ������ ���������� � ����������, �� �������
% ����� ��������� ������� ����� ��� ����/��������/��������, ��� ���� ��
% ��������� �����.
