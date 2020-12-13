clc;
clear;

% method(dimension, number of nonzeros values, lower range, upper range)
A = generateSymmetricSparse(4,3,1,10)


function [A] = generateSymmetricSparse(n, zDensity, vMin, vMax)
    if zDensity > n * n
        disp("La densidad de ceros debe ser menor que nFilas * nColumnas");
        A = zeros(1,1);
        return
    end
    
    if zDensity <= 0
        disp("La densidad de ceros debe ser mayor que 0");
        A = zeros(1,1);
        return
    end
    
    if n <= 0
        disp("La dimensión debe ser al menos 1x1");
        A = zeros(1,1);
        return
    end
    
    A = zeros(n,n);
    
    for i = 1:zDensity
        value = randi([vMin vMax],1,1); % value 
        posRow = randi([1 n],1,1); % position row in matrix
        posCol = randi([1 n],1,1); % position col in matrix
        
        if A(posRow, posCol) == 0
            A(posRow, posCol) = value;
            A(posCol, posRow) = value;
        else
            while A(posRow, posCol) ~= 0
                posRow = randi([1 n],1,1);
                posCol = randi([1 n],1,1);
            end
            A(posRow, posCol) = value;
            A(posCol, posRow) = value;
        end
    end
end