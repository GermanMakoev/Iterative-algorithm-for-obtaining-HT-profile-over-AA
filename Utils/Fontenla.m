[S,P,A]=getFontenla;
Ht1 = [0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.5 4.0 4.5 5.0 5.5  6 10 15 20 25]*1e8;
Ht2 = [1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.3 2.5 2.7 2.9 3.1 3.7 4.2 4.7 5.2 5.7 6.5 10 15 20 25 30]*1e8;
Hc = (Ht1+Ht2)/2;
T1=[5.891e4 1.13e5 1.915e5 2.578e5 3.198e5 3.978e5 4.527e5 5.237e5 5.966e5 6.679e5 7.489e5 8.264e5 9.166e5 1.032e6 1.152e6 1.31e6 1.423e6 1.569e6 1.721e6 1.868e6 2.001e6 2.111e6 2.193e6 2.233e6 2.232e6 2.232e6 2.232e6 2.232e6];
T2=[6e4 8.323e5 8.228e5 8.04e5 7.63e5 7.477e5 7.422e5 7.931e5 7.508e5 7.729e5 7.897e5 8.187e5 8.621e5 9.129e5 9.612e5 1.021e6 1.08e6 1.142e6 1.241e6 1.389e6 1.586e6 1.82e6 2.078e6 2.332e6 2.466e6 2.472e6 2.472e6 2.472e6];    
plot(S(:,1),10.^S(:,2));
title('��������� ���������� ������������-�������� �������� � ��������� Fontenla');
xlabel('������, ��');
ylabel('�����������, �');
hold on;
plot(P(:,1),10.^P(:,2));
hold on;
plot(A(:,1),10.^A(:,2));
hold on;
plot(Hc,T1);
hold on;
plot(Hc,T2)
legend('��������� ������','����','��������','AR11312','AR12740');

