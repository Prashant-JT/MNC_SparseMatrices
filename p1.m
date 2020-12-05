clc;
clear;

A = [1,7,0,0;0,2,8,0,;5,0,3,9;0,6,0,4];

[row, col, val] = COO(A);
disp("COO a la matriz A:");
disp("Vector de filas: ");
disp(row);
disp("Vector de columnas: ");
disp(col);
disp("Vector de valores: ");
disp(val);
disp("");

[row, col, val] = CSR(A);
disp("CSR a la matriz A:");
disp("Vector de filas: ");
disp(row);
disp("Vector de columnas: ");
disp(col);
disp("Vector de valores: ");
disp(val);


function [row, col, val] = COO(A)
    n = nnz(A); % number of nonzero
    row = zeros(1,n);
    col = zeros(1,n);
    val = zeros(1,n);
    
    k = 1;
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            if A(i,j) ~= 0
                row(k) = i;
                col(k) = j;
                val(k) = A(i,j);
                k = k + 1;
            end
        end
    end
end

function [rowOff, col, val] = CSR(A)
    n = nnz(A); % number of nonzero
    nz = 0;
    rowOff = zeros(1,size(A,1)+1);
    col = zeros(1,n);
    val = zeros(1,n);
    
    k = 1;
    l = 1;
    if A(1,1) ~= 0
        rowOff(l) = 0;
        l = 2;
    end
    
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            if A(i,j) ~= 0
                col(k) = j;
                val(k) = A(i,j);
                k = k + 1;
                nz = nz + 1;
            end
        end
        rowOff(l) = nz;
        l = l + 1;
    end
end

function X = sparseMatrix(m, n)
    X = sparse(m,n);
end