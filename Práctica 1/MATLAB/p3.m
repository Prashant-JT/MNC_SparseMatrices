clc;
clear;

for i=1:(10*10)-1
    disp("Matriz 10x10 con nº de ceros = " + i);
    A = generateSparse(10,10,i,1,50);
    [rowCOO, co1COO, valCOO] = COO(A);
    [rowCSR, colCSR, valCSR] = CSR(A);
    a = whos('A').bytes;
    coo = whos('rowCOO').bytes + whos('co1COO').bytes + whos('valCOO').bytes;
    csr = whos('rowCSR').bytes + whos('colCSR').bytes + whos('valCSR').bytes;
    disp("Bytes A = " + a);
    disp("Bytes COO = " + coo);
    disp("Bytes CSR = " + csr);
    disp("-----------------------------------");
end


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

function [A] = generateSparse(nRow, nCol, zDensity, vMin, vMax)
    if zDensity > nRow * nCol
        disp("La densidad de ceros debe ser menor que nFilas * nColumnas");
        A = zeros(1,1);
        return
    end
    
    if zDensity <= 0
        disp("La densidad de ceros debe ser mayor que 0");
        A = zeros(1,1);
        return
    end
    
    if nRow <= 0 || nCol <= 0
        disp("La dimensión debe ser al menos 1x1");
        A = zeros(1,1);
        return
    end
    
    if vMin <= 0 || vMax <= vMin
        disp("Los valores deben ser mayores que 0");
        A = zeros(1,1);
        return
    end
    
    A = randi([vMin vMax],nRow,nCol);
    
    for i = 1:zDensity
        posRow = randi([1 nRow],1,1); % position row in matrix
        posCol = randi([1 nCol],1,1); % position col in matrix
        
        if A(posRow, posCol) ~= 0
            A(posRow, posCol) = 0;
        else
            while A(posRow, posCol) == 0
                posRow = randi([1 nRow],1,1);
                posCol = randi([1 nCol],1,1);
            end
            A(posRow, posCol) = 0;
        end
    end
end