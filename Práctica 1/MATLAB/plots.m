clc;
clear;

data = load('sizes10x10.txt');
dim = size(data);
a = data(:, 1);
coo = data(:, 2);
csr = data(:, 3);

figure(1)
subplot(1,2,1);
plot([1:99], a);
hold on
plot([1:99], coo);
plot([1:99], csr);
title("Dimensi�n de la matriz de 10x10");
ylabel("Tama�o (bytes)");
xlabel("N� de ceros");
legend('Matriz sin codificar','COO', 'CSR', 'Location', 'northeast');
hold off

data = load('sizes100x100.txt');
dim = size(data);
a = data(:, 1);
coo = data(:, 2);
csr = data(:, 3);

figure(1)
subplot(1,2,2);
plot([1:9999], a);
hold on
plot([1:9999], coo);
plot([1:9999], csr);
title("Dimensi�n de la matriz de 100x100");
ylabel("Tama�o (bytes)");
xlabel("N� de ceros");
legend('Matriz sin codificar','COO', 'CSR', 'Location', 'northeast');
hold off