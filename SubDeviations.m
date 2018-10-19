function [ std_frame ] = SubDeviations(N,M,Matrix)
% Gabriel Moreira
% Given a matrix A with dimensions bigger than N and M,
% this function calculates a matrix with the standard deviations of the 
% submatrices of A
[nrows, ncols] = size(Matrix);

std_frame = zeros(floor(nrows/N), floor(ncols/M));

for i = 1:N:nrows
    for j = 1:M:ncols
        if i == 1
            ii = 1;
        else
            ii = floor((i-1)/N);
        end
        if j == 1
            jj = 1;
        else
            jj = floor((j-1)/M);
        end
        i_limit = i+N-1;
        j_limit = j+N-1;
        if i+N-1 > nrows
            i_limit = nrows;
        end
        if j+N-1 > ncols
            j_limit = ncols;
        end
        var_sample = reshape(Matrix(i:i_limit, j:j_limit), [1, ...
            (i_limit-i+1)*(j_limit-j+1)]);
        std_frame(ii, jj) = sqrt(var(var_sample));
    end  
end
end

