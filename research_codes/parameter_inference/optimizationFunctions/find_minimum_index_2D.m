function [indi, indj] = find_minimum_index_2D(mse_vol)
indi = [];
indj = [];
min_v = min(mse_vol(:));
for i = 1 : size(mse_vol,1)
    for j = 1 : size(mse_vol,2)
        if mse_vol(i,j) == min_v
            indi = [indi; i];
            indj = [indj; j];
        end
    end
end