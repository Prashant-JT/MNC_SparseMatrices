clc;
clear;

data = load('P3_times_Blas-CSR-COO.csv');
dim = size(data);
a = data(:, 1);
coo = data(:, 3);
csr = data(:, 2);

figure(1)
plot([1:100:30000], a);
hold on
plot([1:100:30000], coo);
plot([1:100:30000], csr);
title("Dimensión de la matriz de 1000x1000");
ylabel("Tiempo (segundos)");
xlabel("Nº de elementos distinto de 0");
legend('Matriz sin codificar','COO', 'CSR', 'Location', 'northeast');
hold off