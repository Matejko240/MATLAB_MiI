% Zadanie 1: Generowanie liczb losowych
N = 1000; % liczba próbek
X = rand(1, N); % generowanie N liczb losowych z rozkładu U(0,1)

% Zastosowanie funkcji gęstości prawdopodobieństwa f(x) = 2x dla x in [0,1]
X = sqrt(X); % transformacja do rozkładu o gęstości 2x
% Zadanie 2: Dystrybuanta empiryczna
x_vals = linspace(0, 1, 100); % punkty, w których będziemy obliczać dystrybuantę
F_empirical = zeros(1, length(x_vals)); % wektor na dystrybuantę empiryczną

for i = 1:length(x_vals)
    F_empirical(i) = sum(X <= x_vals(i)) / N; % obliczenie dystrybuanty empirycznej
end

% Obliczenie teoretycznej dystrybuanty F(x) = x^2
F_theoretical = x_vals.^2;

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
% Zadanie 4: Analiza danych z pliku
data = load('ModelowanieLab4Data.txt'); % wczytanie danych
N_values = [10, 100, 1000]; % różne wartości N do analizy

% Dystrybuanta empiryczna dla różnych wartości N
figure;
hold on;
for N = N_values
    F_empirical_data = zeros(1, N); % wektor na dystrybuantę empiryczną
    for i = 1:N
        F_empirical_data(i) = sum(data <= data(i)) / N;
    end
    plot(data, F_empirical_data, 'LineWidth', 2);
end
xlabel('Dane');
ylabel('F(x)');
legend(arrayfun(@(x) ['N = ', num2str(x)], N_values, 'UniformOutput', false));
title('Dystrybuanta empiryczna dla różnych N');
grid on;
% Zadanie 5: Estymacja wariancji dystrybuanty empirycznej
x_vals = linspace(0, 1, 100); % punkty, w których będziemy obliczać dystrybuantę
F_empirical = zeros(1, length(x_vals)); % wektor na dystrybuantę empiryczną
var_F = zeros(1, length(x_vals)); % wektor na wariancję

for i = 1:length(x_vals)
    F_empirical(i) = sum(X <= x_vals(i)) / N; % obliczenie dystrybuanty empirycznej
    var_F(i) = var(F_empirical(i)); % obliczenie wariancji
end

% Wykres wariancji
figure;
plot(x_vals, var_F, 'm-', 'LineWidth', 2);
xlabel('x');
ylabel('Wariancja dystrybuanty empirycznej');
title('Wariancja dystrybuanty empirycznej');
grid on;
