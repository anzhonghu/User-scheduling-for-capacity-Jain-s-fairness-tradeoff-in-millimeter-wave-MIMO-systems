clear;
close all;
naz = 32;
nel = 151;
n_arr = naz * nel;
maz = 3;
mel = 4;
m_arr = maz * mel;
Rmin = 10;
Rmax = 100;
beam_numberforMS = 2;
beam_numberforBS = 8;
P = 4;
PB = 1e11;
variable_s = [20; 40; 60;80; 100; ];
height = 10;
time_fre_resources = 50;
f = 80*10^9;%1G bandwidth
lambda = 3 * 10^8 / f;
miu = 0.5;
beta_m = 10.3 * pi / 180;%-phi_m in the EL paper
Nite = 1e3;
capacity = zeros(length(variable_s),  4);
capacity_cdf = zeros(variable_s(length(variable_s), 1),  4);
angle = zeros(1, 2);
U = zeros(n_arr, n_arr);
for nx = 1 : naz
    for ny = 1 : nel
        angle(1, 2) = (-1+2*ny/nel);%el
        angle(1, 1) = (-1+2*nx/naz);%az
        n = (ny - 1) * naz + nx;
        for mx = 0 : naz-1
            for my = 0 : nel-1
                m = my * naz + 1 + mx;
                U(m, n) = exp(-1i * 2 * pi * miu * ((mx-0.5*(naz-1)) * angle(1, 1) + (my-0.5*(nel-1)) * angle(1, 2))) / sqrt(n_arr);
            end
        end
    end
end
Um = zeros(m_arr, m_arr);
for nx = 1 : maz
    for ny = 1 : mel
        angle(1, 2) = (-1+2*ny/mel);%el
        angle(1, 1) = (-1+2*nx/maz);%az
        n = (ny - 1) * maz + nx;
        for mx = 0 : maz-1
            for my = 0 : mel-1
                m = my * maz + 1 + mx;
                Um(m, n) = exp(-1i * 2 * pi * miu * ((mx-0.5*(maz-1)) * angle(1, 1) + (my-0.5*(mel-1)) * angle(1, 2))) / sqrt(m_arr);
            end
        end
    end
end
for variable_n = 1 : length(variable_s)
    K = variable_s(variable_n, 1);
    %H = zeros(K*m_arr, n_arr);
    H1B = zeros(n_arr, K);
    HM = zeros(m_arr, P*K);
    Mb = beam_numberforMS * K;
    Mm = beam_numberforBS;
    for ii = 1 : Nite
        pos = zeros(K, 3);
        theta = zeros(K, P);
        phi = zeros(K, P);
        beta = zeros(K, P);
        theta_m = zeros(K, P);
        phi_m = zeros(K, P);
        p_store = zeros(K, 1);
        block_store = zeros(K, 1);
        ill_cond_s = zeros(K, 1);
        for k = 1 : K
            for p = 1 : P
                if 1==p
                    if rand(1,1)>0.1
                        pos_temp = zeros(1, 2);
                        while norm(pos_temp) < Rmin || norm(pos_temp) > Rmax || abs(atan(pos_temp(1, 2) / pos_temp(1, 1))) > pi / 3
                            pos_temp(1, 1) = rand(1, 1) * Rmax;
                            pos_temp(1, 2) = (rand(1, 1) * 2 - 1) * Rmax;
                        end
                        pos(k, 1:2) = pos_temp;
                        pos(k, 3) = norm(pos_temp);
                        pos(k, 3) = norm([pos(k, 3), height]);%distance
                        phi(k, p) = asin(pos_temp(1, 2) / sqrt(pos_temp(1, 2)^2 + (pos_temp(1, 1) * cos(beta_m) + height * sin(beta_m))^2));%az
                        theta(k, p) = asin((pos_temp(1, 1) * sin(beta_m) - height * cos(beta_m)) / pos(k, 3));%el
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        phi_m(k, p) = pi - phi(k, p);
                        theta_m(k, p) = pi - theta(k, p);
                    else
                        block_store(k, 1) = 1;
                        beta(k, p) = 0;
                        continue;
                    end
                else
                    %take the reflector inside the cell
                    pos_temp = zeros(1, 2);
                    while norm(pos_temp) < Rmin || norm(pos_temp) > Rmax || abs(atan(pos_temp(1, 2) / pos_temp(1, 1))) > pi / 3
                        pos_temp(1, 1) = rand(1, 1) * Rmax;
                        pos_temp(1, 2) = (rand(1, 1) * 2 - 1) * Rmax;
                    end
                    pos(k, 1:2) = pos_temp;
                    pos(k, 3) = norm(pos_temp);
                    pos(k, 3) = norm([pos(k, 3), height]);%distance
                    phi(k, p) = asin(pos_temp(1, 2) / sqrt(pos_temp(1, 2)^2 + (pos_temp(1, 1) * cos(beta_m) + height * sin(beta_m))^2));%az
                    theta(k, p) = asin((pos_temp(1, 1) * sin(beta_m) - height * cos(beta_m)) / pos(k, 3));%el
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    phi_m(k, p) = pi - phi(k, p);
                    theta_m(k, p) = pi - theta(k, p);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Aaz = -min(12 * phi(k, p)^2 / (70/180*pi)^2, 25);
                Ael = -min(12 * theta(k, p)^2 / (7/180*pi)^2, 20);
                D0 = -min(-Aaz-Ael, 25);
                D0 = 10^(D0*0.1);
                if 1==p
                    beta(k, p) = sqrt(D0 * lambda^2 / (16 * pi^2 * pos(k, 3)^2)) * exp(1i * rand(1,1) * 2 * pi);
                else
                    beta(k, p) = sqrt(D0 * lambda^2 / (16 * pi^2 * pos(k, 3)^2)) * exp(1i * rand(1,1) * 2 * pi) * 10^((-rand(1,1) * 5 - 15)*0.05);%-15~-20dB loss
                end
                %                 hb = zeros(n_arr, 1);
                %                 for nx = 0 : naz-1
                %                     for ny = 0 : nel-1
                %                         n = ny * naz+ 1 + nx;
                %                         hb(n, 1) = exp(-1i * 2 * pi * miu * ((nx-0.5*(naz-1)) * cos(theta(k, p)) * sin(phi(k, p)) + (ny-0.5*(nel-1)) * sin(theta(k, p))));
                %                     end
                %                 end
                hm = zeros(m_arr, 1);
                for nx = 0 : maz-1
                    for ny = 0 : mel-1
                        n = ny * maz+ 1 + nx;
                        hm(n, 1) = exp(-1i * 2 * pi * miu * ((nx-0.5*(maz-1)) * cos(theta_m(k, p)) * sin(phi_m(k, p)) + (ny-0.5*(mel-1)) * sin(theta_m(k, p))));
                    end
                end
                %H((k-1)*m_arr+1: k*m_arr, :) = H((k-1)*m_arr+1: k*m_arr,:) + hm * hb' * beta(k, p);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                HM(:, (k-1)*P+p) = hm;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [~, p] = max(abs(beta(k, :)));
            p_store(k, 1) = p;
            for nx = 0 : naz-1
                for ny = 0 : nel-1
                    n = ny * naz+ 1 + nx;
                    H1B(n, k) = exp(-1i * 2 * pi * miu * ((nx-0.5*(naz-1)) * cos(theta(k, p)) * sin(phi(k, p)) + (ny-0.5*(nel-1)) * sin(theta(k, p))));
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Db = zeros(Mb, 1);
        diff = 1e4 * ones(n_arr, 1);
        H_Bb =  U' * H1B;
        [~, uhein] = sort(abs(H_Bb), 'descend');
        for k = 1 : K
            counttemp = 0;
            for n = 1 : n_arr
                if diff(uhein(n,k), 1) > 0
                    diff(uhein(n,k), 1) = 0;
                    counttemp = counttemp + 1;
                    Db((k-1)*beam_numberforMS+counttemp, 1) = uhein(n,k);
                else
                end
                if counttemp < beam_numberforMS
                else
                    break;
                end
            end
        end
        H_Bb_t = H_Bb(Db, :);%%Mb*K
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Dm = zeros(Mm, K);
        H_Mb =  Um' * HM;%m_arr, P*K
        H_Mb_t = zeros(Mm, P*K);%Mm, P*K
        for k = 1 : K
            [~, uhein] = sort(abs(H_Mb(:,(k-1)*P+p_store(k, 1))), 'descend');
            Dm(:, k) = uhein(1:Mm,1);
            H_Mb_t(:, (k-1)*P+1:k*P) = H_Mb(Dm(:,k), (k-1)*P+1:k*P);
            if block_store(k, 1) > 0
                rk_matr_cond = H_Mb_t(:, (k-1)*P+2:k*P)' * H_Mb_t(:, (k-1)*P+2:k*P);
            else
                rk_matr_cond = H_Mb_t(:, (k-1)*P+1:k*P)' * H_Mb_t(:, (k-1)*P+1:k*P);
            end
            if cond(rk_matr_cond) > 1e3
                ill_cond_s(k, 1) = 1;
            else
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cac;%4
        oppor_appr;%1
        maxmi_appr;%2
        pf_appr;%3
        disp([variable_n, ii])
    end
end
capacity(:, 2:4) = capacity(:, 2:4) / Nite / time_fre_resources;
capacity(:, 1) = capacity(:, 1) / Nite;
capacity_cdf(:, 2:4) = capacity_cdf(:, 2:4) / time_fre_resources;
capacity_cdf = sort(capacity_cdf, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1 = subplot(1,2,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(variable_s, capacity(1:length(variable_s), 1), 'k--s','LineWidth',1,'MarkerSize',10)
hold on
plot(variable_s, capacity(1:length(variable_s), 2), 'k--*','LineWidth',1,'MarkerSize',10)
plot(variable_s, capacity(1:length(variable_s), 3), 'k--o','LineWidth',1,'MarkerSize',12)
plot(variable_s, capacity(1:length(variable_s), 4), 'k-^','LineWidth',1,'MarkerSize',10)
xlim([min(variable_s), max(variable_s)])
le = legend('Opportunistic','max-min','PF','Proposed', 'Location', 'northwest');
set(le,'Fontname','Times')
set(gca,'XTick',variable_s)
xlabel('Number of MSs','Fontname','Times')
ylabel('Sum capacity (bps/Hz)','Fontname','Times')
grid on
%%%%%%%%%%%%%%%%%%%%%%%%
h2 = subplot(1,2,2);
capacity_plot = zeros(10,  4);
for plot_nn = 1 : 10
    capacity_plot(plot_nn, :) = capacity_cdf(variable_s(length(variable_s), 1)*plot_nn*0.1, :);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F1 = cdfplot(capacity_cdf(:, 1));
% set(F1,'LineWidth',2,'Color','k','LineStyle', '--')
% hold on
% plot(capacity_plot(:,1),(1:10)*0.1, 'ks','MarkerSize',10)
% F2 = cdfplot(capacity_cdf(:, 2));
% set(F2,'LineWidth',2,'Color','k','LineStyle', '--')
% plot(capacity_plot(:,2),(1:10)*0.1, 'k*','MarkerSize',10)
% F3 = cdfplot(capacity_cdf(:, 3));
% set(F3,'LineWidth',2,'Color','k','LineStyle', '--')
% plot(capacity_plot(:,3),(1:10)*0.1, 'ko','MarkerSize',12)
% F4 = cdfplot(capacity_cdf(:, 4));
% set(F4,'LineWidth',2,'Color','k','LineStyle', '-')
% plot(capacity_plot(:,4),(1:10)*0.1, 'k^','MarkerSize',10)
% xlabel('Capacity (bps/Hz)','Fontname','Times')
% ylabel('CDF','Fontname','Times')
% grid on
% title('')


plot(capacity_plot(:,1),(1:10)*0.1, 'k-s','MarkerSize',10)
hold on
plot(capacity_plot(:,2),(1:10)*0.1, 'k-*','MarkerSize',10)
plot(capacity_plot(:,3),(1:10)*0.1, 'k-o','MarkerSize',12)
plot(capacity_plot(:,4),(1:10)*0.1, 'k-^','MarkerSize',10)
xlabel('Capacity (bps/Hz)','Fontname','Times')
ylabel('CDF','Fontname','Times')
grid on
title('')

