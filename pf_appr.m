capacity_cumu = zeros(K, 1);
capa_cu_s = zeros(K, 1);
for k = 1 : K
    H1B_brev = H_Bb_t(:, k);
    gk_matr = inv(H1B_brev' * H1B_brev);
    sum_nomi = gk_matr(1, 1) * abs(beta(k, p_store(k, 1)))^(-2);
    HM_brev = H_Mb_t(:, (k-1)*P+1:k*P);
    if ill_cond_s(k, 1) < 1
        if block_store(k, 1) > 0
            rk_matr = inv(HM_brev(:, 2:P)' * HM_brev(:, 2:P));
            rk = rk_matr(p_store(k, 1)-1,p_store(k, 1)-1);
        else
            rk_matr = inv(HM_brev(:, 1:P)' * HM_brev(:, 1:P));
            rk = rk_matr(p_store(k, 1),p_store(k, 1));
        end
        capa_cu_s(k, 1) = log2(1 + PB/(sum_nomi * rk));
    else
    end
end
for tfr = 1 : time_fre_resources
    pf_ratio_store = zeros(K, 1);
    for k = 1 : K
        if ill_cond_s(k, 1) < 1
            if capacity_cumu(k, 1) > 0
                pf_ratio_store(k, 1) = capa_cu_s(k, 1) / capacity_cumu(k, 1);
            else
                pf_ratio_store(k, 1) = capa_cu_s(k, 1) * 1e8;%%%?
            end
        else
        end
    end
    [~, Csel] = sort(pf_ratio_store, 'ascend');
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
        capacity(variable_n, 3) = capacity(variable_n, 3) + log2(1 + PB/(sum_nomi * rk));
        if 1==ii && length(variable_s)==variable_n
            capacity_cdf(Csel(nc, 1), 3) = capacity_cdf(Csel(nc, 1), 3) + log2(1 + PB/(sum_nomi * rk));
        else
        end
        capacity_cumu(Csel(nc, 1), 1) = capacity_cumu(Csel(nc, 1), 1) + log2(1 + PB/(sum_nomi * rk));
    end
end


