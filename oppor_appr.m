Csel = zeros(K, 1);
pathloss_store = zeros(K, 1);
for k = 1 : K
    pathloss_store(k, 1) = max(beta(k, :));
end
[~, Csel(:, 1)] = sort(pathloss_store, 'ascend');
C_len = K - sum(ill_cond_s(:, 1));
Csel_temp = zeros(C_len, 1);
k_xx = 0;
for k_x = 1 : K
    if ill_cond_s(Csel(k_x, 1), 1) > 0
    else
        k_xx = k_xx + 1;
        Csel_temp(k_xx, 1) = Csel(k_x, 1);
    end
end
Csel = Csel_temp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H1B_brev = H_Bb_t(:, Csel);
HM_brev = zeros(Mm, P*length(Csel_temp));
for k_xx = 1 : length(Csel_temp)
    HM_brev(:, (k_xx-1)*P+1:k_xx*P) = H_Mb_t(:, (Csel_temp(k_xx,1)-1)*P+1:Csel_temp(k_xx,1)*P);
end
Csel_temp = Csel;
if cond(H1B_brev'*H1B_brev) > 1e3
    for k_xx = 2 : C_len
        H1B_brev_temp = H_Bb_t(:, Csel(k_xx:C_len, 1));
        Csel_temp = Csel(k_xx:C_len, 1);
        if cond(H1B_brev_temp'*H1B_brev_temp) > 1e3
        else
            break;
        end
    end
    H1B_brev = H1B_brev_temp;
    HM_brev = zeros(Mm, P*length(Csel_temp));
    for k_xx = 1 : length(Csel_temp)
        HM_brev(:, (k_xx-1)*P+1:k_xx*P) = H_Mb_t(:, (Csel_temp(k_xx,1)-1)*P+1:Csel_temp(k_xx,1)*P);
    end
    Csel = Csel_temp;
else
end
gk_matr = inv(H1B_brev' * H1B_brev);
%%%%%%%%%%%%%%%%%%%%%%%%%%
sum_nomi = 0;
for nc = 1 : length(Csel)
    sum_nomi = sum_nomi + gk_matr(nc, nc) * abs(beta(Csel(nc, 1), p_store(Csel(nc, 1), 1)))^(-2);
end
for nc = 1 : length(Csel)
    if block_store(Csel(nc, 1), 1) > 0
        rk_matr = inv(HM_brev(:, (nc-1)*P+2:nc*P)' * HM_brev(:, (nc-1)*P+2:nc*P));
        rk = rk_matr(p_store(Csel(nc, 1), 1)-1,p_store(Csel(nc, 1), 1)-1);
    else
        rk_matr = inv(HM_brev(:, (nc-1)*P+1:nc*P)' * HM_brev(:, (nc-1)*P+1:nc*P));
        rk = rk_matr(p_store(Csel(nc, 1), 1),p_store(Csel(nc, 1), 1));
    end
    capacity(variable_n, 1) = capacity(variable_n, 1) + log2(1 + PB/(sum_nomi * rk));
    if 1==ii && length(variable_s)==variable_n
        capacity_cdf(Csel(nc, 1), 1) = capacity_cdf(Csel(nc, 1), 1) + log2(1 + PB/(sum_nomi * rk));
    else
    end
end


