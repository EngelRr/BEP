clear
close all
clc
rng(1);

data = rand(1000,5);
C = parallel.pool.Constant(data);

x = eye(5);
for ii = 1:5
    parfor jj = 1:5
        x(ii,jj) = C.Value(ii,jj);
    end
end
x