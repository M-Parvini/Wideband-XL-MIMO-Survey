function [H, c, lambda] = generate_ULA_parallel_new(fc, Nt, Nr,D_list)
%GENERATE_ULA_PARALLEL 此处显示有关此函数的摘要
%   此处显示详细说明
c = 3e8;
lambda = c/fc;

[~, d_len] = size(D_list);

tx_start = -(Nt/2-0.5)*lambda/2;
tx_loc = tx_start:lambda/2:-tx_start;

rx_start = -(Nr/2-0.5)*lambda/2;
rx_loc = rx_start:lambda/2:-rx_start;

H = zeros(d_len, Nt, Nr);
for i_d = 1:d_len
    D = D_list(i_d);
    dis_matrix = zeros(Nt, Nr);
    for i_tx = 1:Nt
        for i_rx = 1:Nr
            dis_matrix(i_tx, i_rx) = sqrt(D^2 + (tx_loc(i_tx) - rx_loc(i_rx))^2);
        end
    end
    H_tmp = exp(-1i*2*pi/lambda*dis_matrix).';
%     H_tmp = exp(-1i*2*pi/lambda*dis_matrix).'/D;
%     H(i_d, :, :) = H_tmp/sqrt(Nt*Nr);
    H(i_d, :, :) = H_tmp;
end


