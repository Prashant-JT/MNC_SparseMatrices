clc;
clear;

M = generateSparse(4,4,9,1,50)

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