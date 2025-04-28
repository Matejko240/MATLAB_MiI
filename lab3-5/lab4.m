clc;clear;close all;
% Zadanie 1: Generowanie liczb losowych
N = 1000; % liczba próbek
X = rand(1, N); % generowanie N liczb losowych z rozkładu U(0,1)
% to że histogram jest równy jeszcze nie oznacza że jest dobry
% Zastosowanie funkcji gęstości prawdopodobieństwa f(x) = 20x dla x in [0,1]
X = sqrt(X / 10) ; % transformacja do rozkładu o gęstości 20x
% Zadanie 2: Dystrybuanta empiryczna
x_vals = linspace(0, 1, 100); % punkty, w których będziemy obliczać dystrybuantę
usmiech=sqrt(1000/5);
F_empirical = zeros(1, length(x_vals)); % wektor na dystrybuantę empiryczną

for i = 1:length(x_vals)
    F_empirical(i) = sum(X <= x_vals(i)) / N; % obliczenie dystrybuanty empirycznej
end

% Obliczenie teoretycznej dystrybuanty F(x) = x^2
F_theoretical = 10*x_vals.^2;
for i = 1:length(x_vals)
    if F_theoretical(i) >=1
        F_theoretical(i)=1;
    end
end
% Wykres
figure;
plot(x_vals, F_theoretical, 'r-', 'LineWidth', 2); % dystrybuanta teoretyczna
hold on;
plot(x_vals, F_empirical, 'b--', 'LineWidth', 2); % dystrybuanta empiryczna
xlabel('x');
ylabel('F(x)');
legend('Dystrybuanta teoretyczna', 'Dystrybuanta empiryczna');
title('Dystrybuanta empiryczna i teoretyczna');
grid on;
% Zadanie 3: Błąd estymatora
D_N = zeros(1, N); % wektor na błąd estymatora
for i = 1:N
    F_empirical_i = sum(X <= X(i)) / N; % dystrybuanta empiryczna dla X(i)
    D_N(i) = abs(F_empirical_i - X(i)^2); % obliczenie błędu
end

% Wykres błędu estymatora
figure;
plot(1:N, D_N, 'g-', 'LineWidth', 2);
xlabel('Indeks próbki');
ylabel('Błąd estymatora');
title('Błąd estymatora dystrybuanty');
grid on;
data = sqrt(rand(1,1000)); % dla f(x)=2x, jak wcześniej
x = linspace(0,1,1000);
F_teor = x.^2;

[f_emp, x_emp] = ecdf(data);

% Usunięcie powtórzeń w x_emp (próbkach ECDF)
[x_emp_unique, idx_unique] = unique(x_emp);
f_emp_unique = f_emp(idx_unique);

% Interpolacja empirycznej dystrybuanty w punktach X
F_interp = interp1(x_emp_unique, f_emp_unique, data, 'linear', 'extrap');


blad = abs(F_interp - data.^2); % |F_N(X_i) - F(X_i)|

scatter(data, blad, 10, 'filled')
xlabel('X_i');
ylabel('Błąd estymatora |F_N(X_i) - F(X_i)|');
title('Błąd estymatora dystrybuanty empirycznej jako funkcja X_i');
grid on;


% Zadanie 4: Analiza danych z pliku
data = load('ModelowanieLab4Data.txt'); % wczytanie danych
N_values = [10, 100, 1000];
x_vals = linspace(min(data), max(data), 1000); % ustalona siatka punktów

% Dystrybuanta empiryczna dla różnych wartości N
figure;
hold on;
for idx = 1:length(N_values)
    N = N_values(idx);
    F_empirical_data = zeros(1, length(x_vals));
    for i = 1:length(x_vals)
        F_empirical_data(i) = sum(data(1:N) <= x_vals(i)) / N;
    end
    plot(x_vals, F_empirical_data, 'LineWidth', 2, 'DisplayName', ['N = ', num2str(N)]);
end
xlabel('Dane');
ylabel('F(x)');
legend;
title('Dystrybuanta empiryczna dla różnych N');
grid on;
% Dodanie dystrybuant teoretycznych do porównania
x_grid = linspace(min(data), max(data), 1000);
xlim([-50 50]);

F_norm_11 = normcdf(x_grid, 1, 1);
F_norm_05 = normcdf(x_grid, 0, sqrt(5));
F_cauchy = cdf(makedist('tLocationScale', 0, 1, 1), x_grid);


figure;
hold on;
[f_emp, x_emp] = ecdf(data);
plot(x_emp, f_emp, 'k-', 'LineWidth', 2, 'DisplayName', 'Empiryczna');

x = linspace(min(data), max(data), 1000);
plot(x, normcdf(x,1,1), '--r', 'LineWidth', 1.5, 'DisplayName', 'N(1,1)');
plot(x, normcdf(x,0,sqrt(5)), '--b', 'LineWidth', 1.5, 'DisplayName', 'N(0,5)');
plot(x, cdf(makedist('tLocationScale',0,1,1), x), '--g', 'LineWidth', 1.5, 'DisplayName', 'Cauchy(0,1)');

xlim([-50 50]); % lub inny sensowny zakres
xlabel('x');
ylabel('F(x)');
title('Porównanie dystrybuanty empirycznej z teoretycznymi');
legend('show');
grid on;


% Zadanie 5: Estymacja wariancji dystrybuanty empirycznej
x_vals = linspace(0, 1, 100); % przedział [0,1] bo rozkład teoretyczny 2x na [0,1]
L = 100; % liczba powtórzeń
N = 1000;
F_matrix = zeros(L, length(x_vals)); % macierz: wiersze - realizacje, kolumny - punkty x

for l = 1:L
    X = sqrt(rand(1, N) / 10); % losowanie z gęstości 2x
    for i = 1:length(x_vals)
        F_matrix(l, i) = sum(X <= x_vals(i)) / N;
    end
end

% Wariancja w każdym punkcie x
var_F = var(F_matrix, 0, 1); % wariancja wzdłuż wierszy (czyli między próbami)

% Wykres wariancji
figure;
plot(x_vals, var_F, 'm-', 'LineWidth', 2);
xlabel('x');
ylabel('Wariancja dystrybuanty empirycznej');
title('Wariancja dystrybuanty empirycznej (L=100, N=1000)');
grid on;

