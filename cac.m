%%%%grouping%%%%%%%%%%%%%%%%%%%%%%%
group_num = zeros(K, 3);
group_count = zeros(n_arr, 1);
flag_group = zeros(n_arr, 1);
number_of_group = 0;
for k = 1 : K
    if ill_cond_s(k, 1) > 0
        continue;
    else
    end
    [~, p] = max(beta(k, :));
    um_p1 = cos(theta(k, p)) * sin(phi(k, p));
    vm_p1 = sin(theta(k, p));
    m1 = (um_p1+1) * 0.5 * naz;
    m2 = (vm_p1+1) * 0.5 * nel;
    if m1 - floor(m1) <= 0.5
        m1 = floor(m1);
    else
        m1 = ceil(m1);
    end
    if m2 - floor(m2) <= 0.5
        m2 = floor(m2);
    else
        m2 = ceil(m2);
    end
    m_all = m2 * naz + m1 + 1;
    group_num(k, :) = [m1, m2, m_all];
    group_count(m_all, 1) = group_count(m_all, 1) + 1;
    if flag_group(m_all, 1) < 1
        flag_group(m_all, 1) = 1;
        number_of_group = number_of_group + 1;
    else
    end
end
group_store = zeros(n_arr, max(group_count));
for k = 1 : K
    if ill_cond_s(k, 1) > 0
        continue;
    else
    end
    for nn_g = 1 : max(group_count)
        if group_store(group_num(k, 3), nn_g) < 1
            break;
        else
        end
    end
    group_store(group_num(k, 3), nn_g) = k;
end
[~, group_order] = sort(group_count, 'ascend');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
capacity_cumu = zeros(K, 1);
capacity_group = zeros(n_arr, 1);
for tfr = 1 : time_fre_resources
    Csel = zeros(number_of_group, 1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%MS scheduling%%%%%%%%%%%%%%
    number_of_nonscheduled_group = 0;
    for nn_gn = n_arr-number_of_group+1 : n_arr
        group_serial = group_order(nn_gn, 1);
        capacity_temp = zeros(group_count(group_serial, 1), 2);
        if n_arr-number_of_group+1 == nn_gn
            for nn_g = 1 : group_count(group_serial, 1)
                k = group_store(group_serial, nn_g);
                Csel_temp = k;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                H1B_brev = H_Bb_t(:, Csel_temp);
                HM_brev = zeros(Mm, P*length(Csel_temp));
                for k_xx = 1 : length(Csel_temp)
                    HM_brev(:, (k_xx-1)*P+1:k_xx*P) = H_Mb_t(:, (Csel_temp(k_xx,1)-1)*P+1:Csel_temp(k_xx,1)*P);
                end
                gk_matr = inv(H1B_brev' * H1B_brev);
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                sum_nomi = 0;
                for nc = 1 : length(Csel_temp)
                    sum_nomi = sum_nomi + gk_matr(nc, nc) * abs(beta(Csel_temp(nc, 1), p_store(Csel_temp(nc, 1), 1)))^(-2);
                end
                for nc = 1 : length(Csel_temp)
                    if block_store(Csel_temp(nc, 1), 1) > 0
                        rk_matr = inv(HM_brev(:, (nc-1)*P+2:nc*P)' * HM_brev(:, (nc-1)*P+2:nc*P));
                        rk = rk_matr(p_store(Csel_temp(nc, 1), 1)-1,p_store(Csel_temp(nc, 1), 1)-1);
                    else
                        rk_matr = inv(HM_brev(:, (nc-1)*P+1:nc*P)' * HM_brev(:, (nc-1)*P+1:nc*P));
                        rk = rk_matr(p_store(Csel_temp(nc, 1), 1),p_store(Csel_temp(nc, 1), 1));
                    end
                    if capacity_cumu(k, 1) > 0
                        capacity_temp(nn_g, 1) = capacity_temp(nn_g, 1) + log2(1 + PB/(sum_nomi * rk)) / capacity_cumu(k, 1);
                        capacity_temp(nn_g, 2) = capacity_temp(nn_g, 2) + log2(1 + PB/(sum_nomi * rk));
                    else
                        capacity_temp(nn_g, 1) = capacity_temp(nn_g, 1) + log2(1 + PB/(sum_nomi * rk)) * 1e8;%%%?
                        capacity_temp(nn_g, 2) = capacity_temp(nn_g, 2) + log2(1 + PB/(sum_nomi * rk));%%%?
                    end
                end
            end
            [~, nn_g] = max(capacity_temp(:, 1));
            capacity_check = capacity_temp(nn_g, 2);
            Csel(nn_gn-n_arr+number_of_group, 1) = group_store(group_serial, nn_g);
        else
            n_temp_Csel = 0;
            Csel_temp = zeros(nn_gn - n_arr+number_of_group - number_of_nonscheduled_group, 1);
            for mm_gm = 1 : nn_gn - n_arr+number_of_group-1
                if Csel(mm_gm, 1)<1
                else
                    n_temp_Csel = n_temp_Csel + 1;
                    Csel_temp(n_temp_Csel, 1) = Csel(mm_gm, 1);
                end
            end
            for nn_g = 1 : group_count(group_serial, 1)
                k = group_store(group_serial, nn_g);
                Csel_temp(nn_gn - n_arr+number_of_group - number_of_nonscheduled_group, 1) = k;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                H1B_brev = H_Bb_t(:, Csel_temp);
                if cond(H1B_brev'*H1B_brev) > 1e3
                    continue;
                else
                end
                HM_brev = zeros(Mm, P*length(Csel_temp));
                for k_xx = 1 : length(Csel_temp)
                    HM_brev(:, (k_xx-1)*P+1:k_xx*P) = H_Mb_t(:, (Csel_temp(k_xx,1)-1)*P+1:Csel_temp(k_xx,1)*P);
                end
                gk_matr = inv(H1B_brev' * H1B_brev);
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                sum_nomi = 0;
                for nc = 1 : length(Csel_temp)
                    sum_nomi = sum_nomi + gk_matr(nc, nc) * abs(beta(Csel_temp(nc, 1), p_store(Csel_temp(nc, 1), 1)))^(-2);
                end
                for nc = 1 : length(Csel_temp)
                    if block_store(Csel_temp(nc, 1), 1) > 0
                        rk_matr = inv(HM_brev(:, (nc-1)*P+2:nc*P)' * HM_brev(:, (nc-1)*P+2:nc*P));
                        rk = rk_matr(p_store(Csel_temp(nc, 1), 1)-1,p_store(Csel_temp(nc, 1), 1)-1);
                    else
                        rk_matr = inv(HM_brev(:, (nc-1)*P+1:nc*P)' * HM_brev(:, (nc-1)*P+1:nc*P));
                        rk = rk_matr(p_store(Csel_temp(nc, 1), 1),p_store(Csel_temp(nc, 1), 1));
                    end
                    if capacity_cumu(k, 1) > 0
                        capacity_temp(nn_g, 1) = capacity_temp(nn_g, 1) + log2(1 + PB/(sum_nomi * rk)) / capacity_cumu(k, 1);
                        capacity_temp(nn_g, 2) = capacity_temp(nn_g, 2) + log2(1 + PB/(sum_nomi * rk));
                    else
                        capacity_temp(nn_g, 1) = capacity_temp(nn_g, 1) + log2(1 + PB/(sum_nomi * rk)) * 1e8;%%%%%?
                        capacity_temp(nn_g, 2) = capacity_temp(nn_g, 2) + log2(1 + PB/(sum_nomi * rk));%%%%%?
                    end
                end
            end
            [~, nn_g] = max(capacity_temp(:, 1));
            capa_nn_g = capacity_temp(nn_g, 2);
            if capa_nn_g > capacity_check
                Csel(nn_gn-n_arr+number_of_group, 1) = group_store(group_serial, nn_g);
                capacity_check = capa_nn_g;
            else
                number_of_nonscheduled_group = number_of_nonscheduled_group + 1;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%Csel as the set of selected MSs%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    n_temp_Csel = 0;
    if number_of_nonscheduled_group < 1
    else
        Csel_temp = zeros(number_of_group - number_of_nonscheduled_group, 1);
        for mm_gm = 1 : number_of_group
            if Csel(mm_gm, 1)<1
            else
                n_temp_Csel = n_temp_Csel + 1;
                Csel_temp(n_temp_Csel, 1) = Csel(mm_gm, 1);
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    H1B_brev = H_Bb_t(:, Csel_temp);
    HM_brev = zeros(Mm, P*length(Csel_temp));
    for k_xx = 1 : length(Csel_temp)
        HM_brev(:, (k_xx-1)*P+1:k_xx*P) = H_Mb_t(:, (Csel_temp(k_xx,1)-1)*P+1:Csel_temp(k_xx,1)*P);
    end
    gk_matr = inv(H1B_brev' * H1B_brev);
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    sum_nomi = 0;
    for nc = 1 : length(Csel_temp)
        sum_nomi = sum_nomi + gk_matr(nc, nc) * abs(beta(Csel_temp(nc, 1), p_store(Csel_temp(nc, 1), 1)))^(-2);
    end
    for nc = 1 : length(Csel_temp)
        if block_store(Csel_temp(nc, 1), 1) > 0
            rk_matr = inv(HM_brev(:, (nc-1)*P+2:nc*P)' * HM_brev(:, (nc-1)*P+2:nc*P));
            rk = rk_matr(p_store(Csel_temp(nc, 1), 1)-1,p_store(Csel_temp(nc, 1), 1)-1);
        else
            rk_matr = inv(HM_brev(:, (nc-1)*P+1:nc*P)' * HM_brev(:, (nc-1)*P+1:nc*P));
            rk = rk_matr(p_store(Csel_temp(nc, 1), 1),p_store(Csel_temp(nc, 1), 1));
        end
        capacity(variable_n, 4) = capacity(variable_n, 4) + log2(1 + PB/(sum_nomi * rk));
        if 1==ii && length(variable_s)==variable_n
            capacity_cdf(Csel_temp(nc, 1), 4) = capacity_cdf(Csel_temp(nc, 1), 4) + log2(1 + PB/(sum_nomi * rk));
        else
        end
        capacity_cumu(Csel_temp(nc, 1), 1) = capacity_cumu(Csel_temp(nc, 1), 1) + log2(1 + PB/(sum_nomi * rk));
        capacity_group(group_num(Csel_temp(nc, 1), 3), 1) = capacity_group(group_num(Csel_temp(nc, 1), 3), 1) + log2(1 + PB/(sum_nomi * rk));
    end
    order_temp = 0;
    group_order_temp = zeros(number_of_group, 1);
    for nn_x = 1 : number_of_group
        if capacity_group(group_order(n_arr-number_of_group+nn_x, 1), 1) > 0
        else
            order_temp = order_temp + 1;
            group_order_temp(order_temp, 1) = group_order(n_arr-number_of_group+nn_x, 1);
        end
    end
    capa_group_store = zeros(number_of_group-order_temp, 2);
    order_temp_temp = 0;
    for nn_x = 1 : number_of_group
        if capacity_group(group_order(n_arr-number_of_group+nn_x, 1), 1) > 0
            order_temp_temp = order_temp_temp  + 1;
            capa_group_store(order_temp_temp, 1) = capacity_group(group_order(n_arr-number_of_group+nn_x, 1), 1);
            capa_group_store(order_temp_temp, 2) = group_order(n_arr-number_of_group+nn_x, 1);
        else
        end
    end
    [~, capa_group_store_ind] = sort(capa_group_store(:,1), 'ascend');
    group_order_temp(order_temp+1:number_of_group, 1) = capa_group_store(capa_group_store_ind, 2);
    group_order(n_arr-number_of_group+1:n_arr, 1) =  group_order_temp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


