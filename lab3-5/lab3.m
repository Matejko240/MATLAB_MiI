% daszek - estymata
% L - liczba powtórzeń eksperymentu
% N - ilość próbek

clear workspace;
close all;
% Parametry rozkładu normalnego
mu = 0.5; %średnia wart oczekiwana
sigma = 1.5; %odchylenie standardowe
N = 1000; % liczba próbek
%próbki z rozkładu normalnego
X = mu + sigma * randn(1,N);
%estymator średniej wartości oczekiwanej
mu_hat=mean(X);
%estymator wariancji populacyjny
sigma2_hat=var(X);
%estymator wariancji próbkowy
S2_hat=var(X,1);
fprintf('est śrL %.4f\n',mu_hat);
fprintf('est war pop %.4f\n',sigma2_hat);
fprintf('est war prób %.4f\n',S2_hat);
%==============================================
%zad2
L_values = 1:10:1000; % liczba realizacji estymatora
N = 1000; % różne liczby próbek N
errors_mu = zeros(length(L_values), 1);
errors_sigma2 = zeros(length(L_values), 1);
errors_S2 = zeros(length(L_values), 1);

for i = 1:length(L_values)
L = L_values(i);
mu_estimates = zeros(1, L);
sigma2_estimates = zeros(1, L);
S2_estimates = zeros(1, L);

for l = 1:L
X = mu + sigma * randn(1, N); % generowanie próby
mu_estimates(l) = mean(X); % estymator średniej
sigma2_estimates(l) = var(X); % estymator wariancji (populacyjny)
S2_estimates(l) = var(X, 1); % estymator wariancji (próbowy)
end

% błąd empiryczny
errors_mu(i) = mean((mu_estimates - mu).^2);
errors_sigma2(i) = mean((sigma2_estimates - sigma^2).^2);
errors_S2(i) = mean((S2_estimates - sigma^2).^2);
end

% Wykresy
figure;
subplot(1,3,1);
plot(L_values, errors_mu);
xlabel('L');
ylabel('Błąd empiryczny \mu');
title('Błąd \mu');

subplot(1,3,2);
plot(L_values, errors_sigma2);
xlabel('L');
ylabel('Błąd empiryczny \sigma^2');
title('Błąd \sigma^2');

subplot(1,3,3);
plot(L_values, errors_S2);
xlabel('L');
ylabel('Błąd empiryczny S^2');
title('Błąd S^2');

% ANALIZA WPŁYWU LICZBY PRÓBEK N NA BŁĄD EMPIRYCZNY
N_values = [10, 50, 100, 500, 1000, 5000];
L = 1000;
errors_mu_N = zeros(length(N_values), 1);
errors_sigma2_N = zeros(length(N_values), 1);
errors_S2_N = zeros(length(N_values), 1);

for i = 1:length(N_values)
    N = N_values(i);
    mu_estimates = zeros(1, L);
    sigma2_estimates = zeros(1, L);
    S2_estimates = zeros(1, L);

    for l = 1:L
        X = mu + sigma * randn(1, N);
        mu_estimates(l) = mean(X);
        sigma2_estimates(l) = var(X);
        S2_estimates(l) = var(X, 1);
    end

    errors_mu_N(i) = mean((mu_estimates - mu).^2);
    errors_sigma2_N(i) = mean((sigma2_estimates - sigma^2).^2);
    errors_S2_N(i) = mean((S2_estimates - sigma^2).^2);
end

figure;
semilogx(N_values, errors_mu_N, '-o', 'DisplayName', 'Błąd μ');
hold on;
semilogx(N_values, errors_sigma2_N, '-s', 'DisplayName', 'Błąd σ²');
semilogx(N_values, errors_S2_N, '-^', 'DisplayName', 'Błąd S²');
xlabel('Liczba próbek N');
ylabel('Błąd empiryczny');
legend;
title('Błąd estymatorów w funkcji liczby próbek N');
grid on;


% Parametry rozkładu Cauchy'ego
mu_cauchy = 0; % środek
sigma_cauchy = 1; % skala
N = 1000; % liczba próbek

% Generowanie próbek z rozkładu Cauchy'ego
X_cauchy = mu_cauchy + sigma_cauchy * tan(pi * (rand(1, N) - 0.5));

% Estymatory dla rozkładu Cauchy'ego
mu_hat_cauchy = mean(X_cauchy);
sigma2_hat_cauchy = var(X_cauchy);
S2_hat_cauchy = var(X_cauchy, 1);

% Wyświetlanie wyników dla rozkładu Cauchy'ego
fprintf('Estymator średniej (Cauchy): %.4f\n', mu_hat_cauchy);
fprintf('Estymator wariancji (populacyjny, Cauchy): %.4f\n', sigma2_hat_cauchy);
fprintf('Estymator wariancji (próbowy, Cauchy): %.4f\n', S2_hat_cauchy);
