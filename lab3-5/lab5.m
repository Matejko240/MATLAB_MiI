%ZADANIE 1 - Generowanie Danych
% Parametry
N = 500;         % liczba próbek
mu = 1;          % średnia
sigma = 1;       % odchylenie standardowe

% Generowanie danych z rozkładu normalnego N(1,1)
X = mu + sigma * randn(N, 1);

% Podgląd histogramu
figure;
histogram(X, 'Normalization', 'pdf','NumBins', 50);
title('Histogram próbek N(1,1)');
xlabel('x'); ylabel('gęstość');

%ZADANIE 2 – Estymator jądrowy, jądro prostokątne
% Parametry
hN_values = [0.1, 0.5, 1];   % różne wartości wygładzania
x = linspace(-3, 5, 1000);   % zakres oceny gęstości

% Funkcja jądra prostokątnego
K_rect = @(u) 0.5 * (abs(u) <= 1);   % prostokątne: K(u) = 0.5 dla |u| <= 1

figure;
hold on;
for hN = hN_values
    f_hat = zeros(size(x));
    for i = 1:length(x)
        u = (X - x(i)) / hN;
        f_hat(i) = sum(K_rect(u)) / (N * hN);
    end
    plot(x, f_hat, 'DisplayName', ['h = ' num2str(hN)]);
end
title('Estymator jądrowy - jądro prostokątne');
xlabel('x'); ylabel('f̂_N(x)');
legend show;
hold off;

%ZADANIE 3 – Porównanie różnych funkcji jądra

hN = 0.5; % stałe hN
x = linspace(-3, 5, 1000);

% Przykładowe jądra
kernels = {
    @(u) 0.5 * (abs(u) <= 1), 'Prostokątne';
    @(u) (3/4)*(1 - u.^2) .* (abs(u) <= 1), 'Epanechnikov';
    @(u) (1/sqrt(2*pi)) * exp(-0.5*u.^2), 'Gaussowskie'
};

figure;
hold on;
for k = 1:size(kernels,1)
    K = kernels{k,1};
    label = kernels{k,2};

    f_hat = zeros(size(x));
    for i = 1:length(x)
        u = (X - x(i)) / hN;
        f_hat(i) = sum(K(u)) / (N * hN);
    end
    plot(x, f_hat, 'DisplayName', label);
end
title(['Estymator dla różnych jąder (h = ' num2str(hN) ')']);
xlabel('x'); ylabel('f̂_N(x)');
legend show;
hold off;

% ZADANIE 4 – Błąd empiryczny (L prób)
% Parametry
N = 500; L = 10; M = 100;
x_vals = linspace(-3, 5, M);
hN_vec = linspace(0.1, 1.5, 20);
true_pdf = @(x) normpdf(x, 1, 1);  % gęstość rozkładu N(1,1)

% Jądro Gaussowskie
K = @(u) (1/sqrt(2*pi)) * exp(-0.5 * u.^2);

% Błąd empiryczny dla różnych hN
Err = zeros(size(hN_vec));

for j = 1:length(hN_vec)
    hN = hN_vec(j);
    total_error = 0;
    for l = 1:L
        X = 1 + randn(N,1); % próbka

        for m = 1:M
            u = (X - x_vals(m)) / hN;
            f_hat = sum(K(u)) / (N * hN);
            total_error = total_error + (f_hat - true_pdf(x_vals(m)))^2;
        end
    end
    Err(j) = total_error / (L * M);
end

figure;
plot(hN_vec, Err, '-o');
title('Błąd empiryczny estymatora jądrowego');
xlabel('h_N'); ylabel('Err[f̂_N]');

% CROSS-VALIDATION – AUTOMATYCZNY DOBÓR hN
N = 500;
X = 1 + randn(N, 1); % próbka z N(1,1)
hN_vec = linspace(0.1, 1.5, 30);
J_vals = zeros(size(hN_vec));

for j = 1:length(hN_vec)
    hN = hN_vec(j);
    % Estymator z całej próbki
    x_grid = linspace(-3, 5, 200);
    f_hat_sq = 0;
    for x = x_grid
        u = (X - x) / hN;
        f_hat = mean(normpdf(u)); % Gauss
        f_hat_sq = f_hat_sq + f_hat^2;
    end
    f_hat_sq = f_hat_sq * (x_grid(2) - x_grid(1)); % całkowanie prostokątne

    % Leave-one-out
    loo_sum = 0;
    for k = 1:N
        X_loo = X([1:k-1, k+1:end]);
        u = (X_loo - X(k)) / hN;
        f_hat_k = mean(normpdf(u));
        loo_sum = loo_sum + f_hat_k;
    end
    J_vals(j) = f_hat_sq - 2 * loo_sum / N;
end

figure;
plot(hN_vec, J_vals, '-o');
xlabel('hN');
ylabel('Ĵ(hN)');
title('Dobór parametru wygładzania hN metodą cross-validation');
grid on;
