%% ZADANIE 1 i 2 - Generowanie danych i chmura punktów
N = 500;
a = 1;
sigma_Z = 0.1;
X = -2 + 4 * rand(N, 1);
Z = sigma_Z * randn(N, 1);
m = @(x) atan(a * x);
Y = m(X) + Z;

figure;
scatter(X, Y, 10, 'filled');
hold on;
fplot(m, [-2.5, 2.5], 'r', 'LineWidth', 2);
title('Charakterystyka systemu i chmura pomiarów');
xlabel('x'); ylabel('y');
grid on;
legend('Punkty pomiarowe', 'Prawdziwa funkcja m(x)', 'Location', 'best');

%% ZADANIE 3 - Estymator jądrowy (jądro prostokątne)
hN = 0.2;
x_eval = linspace(-2.5, 2.5, 200)';
K_rect = @(u) 0.5 * (abs(u) <= 1);
m_hat = zeros(size(x_eval));

for i = 1:length(x_eval)
    u = (X - x_eval(i)) / hN;
    w = K_rect(u);
    if sum(w) ~= 0
        m_hat(i) = sum(Y .* w) / sum(w);
    else
        m_hat(i) = NaN;
    end
end

figure;
plot(x_eval, m_hat, 'b', 'LineWidth', 2);
hold on;
fplot(m, [-2.5, 2.5], 'r--', 'LineWidth', 2);
title('Estymator m̂_N(x) i funkcja m(x)');
xlabel('x'); ylabel('m̂_N(x)');
legend('Estymator', 'Prawdziwa funkcja');
grid on;

%% ZADANIE 4 - Różne funkcje jądra
different_kernels = {@(u) 0.5*(abs(u)<=1), @(u) exp(-0.5*u.^2)/sqrt(2*pi), @(u) 0.75*(1-u.^2).*(abs(u)<=1)};
kernel_names = {'prostokątne', 'Gaussa', 'Epanechnikova'};
hN = 0.2;

figure;
for j = 1:3
    for i = 1:length(x_eval)
        u = (X - x_eval(i)) / hN;
        w = different_kernels{j}(u);
        if sum(w) ~= 0
            m_hat(i) = sum(Y .* w) / sum(w);
        else
            m_hat(i) = NaN;
        end
    end
    subplot(3,1,j);
    plot(x_eval, m_hat, 'b', 'LineWidth', 1.5);
    hold on;
    fplot(m, [-2.5, 2.5], 'r--');
    title(['Jądro ' kernel_names{j}]);
    grid on;
end

%% ZADANIE 5 - Optymalizacja h
Q = 100;
x_q = linspace(-1, 1, 2*Q+1);
m_true = m(x_q);
h_space = linspace(0.05, 1, 100);
errors = zeros(size(h_space));

for k = 1:length(h_space)
    hN = h_space(k);
    m_hat_q = zeros(size(x_q));
    for i = 1:length(x_q)
        u = (X - x_q(i)) / hN;
        w = K_rect(u);
        if sum(w) ~= 0
            m_hat_q(i) = sum(Y .* w) / sum(w);
        else
            m_hat_q(i) = NaN;
        end
    end
    errors(k) = mean((m_hat_q - m_true).^2, 'omitnan');
end

[~, idx_opt] = min(errors);
h_opt = h_space(idx_opt);

figure;
plot(h_space, errors, 'k', 'LineWidth', 2);
title('Błąd walidacyjny w funkcji h');
xlabel('h'); ylabel('Błąd');
grid on;

%% ZADANIE 6 - Estymator z optymalnym h
hN = h_opt;
for i = 1:length(x_eval)
    u = (X - x_eval(i)) / hN;
    w = K_rect(u);
    if sum(w) ~= 0
        m_hat(i) = sum(Y .* w) / sum(w);
    else
        m_hat(i) = NaN;
    end
end

figure;
plot(x_eval, m_hat, 'b', 'LineWidth', 2);
hold on;
fplot(m, [-2.5, 2.5], 'r--', 'LineWidth', 2);
title(['Estymator z optymalnym h = ' num2str(h_opt)]);
xlabel('x'); ylabel('m̂_N(x)');
legend('Estymator', 'Prawdziwa funkcja');
grid on;

%% ZADANIE 7 - Zakłócenie z rozkładu Cauchyego
cauchy_rng = @(n, gamma) gamma * tan(pi * (rand(n,1) - 0.5));
gamma_vals = [0.01, 0.1, 0.5];

figure;
for g = 1:length(gamma_vals)
    gamma = gamma_vals(g);
    Z_cauchy = cauchy_rng(N, gamma);
    Yc = m(X) + Z_cauchy;
    for i = 1:length(x_eval)
        u = (X - x_eval(i)) / h_opt;
        w = K_rect(u);
        if sum(w) ~= 0
            m_hat(i) = sum(Yc .* w) / sum(w);
        else
            m_hat(i) = NaN;
        end
    end
    subplot(3,1,g);
    plot(x_eval, m_hat, 'b', 'LineWidth', 1.5);
    hold on;
    fplot(m, [-2.5, 2.5], 'r--');
    title(['Zakłócenie Cauchyego, \gamma = ' num2str(gamma)]);
    grid on;
end