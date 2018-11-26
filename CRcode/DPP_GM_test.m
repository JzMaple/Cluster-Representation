clear K_ker;
histo = zeros(10,10);
K_ker = decompose_kernel(KK_full);

for i = 1:1000
    K_idx = sample_dpp(K_ker,10);
    if K_idx(1) == -1
        continue;
    end
    [x, y] = ind2sub([10 10],K_idx);
    for j = 1:10
        histo(x(j),y(j)) = histo(x(j),y(j)) + 1;
    end
end