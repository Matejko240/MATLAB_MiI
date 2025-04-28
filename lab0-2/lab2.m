clc; clear; close all;

N = 100000; % Liczba generowanych próbek
edges = 100;

%% Zadanie 1: Histogram dla f(x)
samples = zeros(N,1);
for i = 1:N
    accept = false;
    while ~accept
        x = 2 * rand - 1;  % Generowanie x w (-1,1)
        u = rand;
        if x <= 0
            f_x = x + 1;
        else
            f_x = -x + 1;
        end
        if u <= f_x
            accept = true;
            samples(i) = x;
        end
    end
end

figure;
histogram(samples, edges, 'Normalization', 'pdf');
hold on;

% Teoretyczny rozkład
x = linspace(-1, 1, 1000);
y = arrayfun(@(x) (x + 1) .* (x <= 0) + (-x + 1) .* (x > 0), x);
plot(x, y, 'r', 'LineWidth', 2);

title('Rozkład f(x) = x+1 dla x \in (-1,0] oraz -x+1 dla x \in (0,1]');
xlabel('Wartości');
ylabel('Gęstość prawdopodobieństwa');
legend('Histogram', 'Funkcja teoretyczna');
grid on;
saveas(gcf, 'Rozklad_Zadanie1.png');

%% Zadanie 2: Wyznaczenie stałej c i implementacja generatora
syms c
c = solve(50 * (1/100) + c * (99/100) == 1, c);
c = double(c);

samples = zeros(N,1);
for i = 1:N
    accept = false;
    while ~accept
        x = rand;  % Generowanie x w (0,1)
        u = rand;
        if x <= 1/100
            f_x = 50;
        else
            f_x = c;
        end
        if u <= f_x
            accept = true;
            samples(i) = x;
        end
    end
end

figure;
histogram(samples, edges, 'Normalization', 'pdf');
title('Histogram dla zadania 2');
xlabel('Wartości');
ylabel('Gęstość prawdopodobieństwa');
grid on;
saveas(gcf, 'Rozklad_Zadanie2.png');

%% Zadanie 3: Generacja liczb z rozkładu półokręgu
samples = zeros(N,1);
for i = 1:N
    accept = false;
    while ~accept
        x = 2 * rand - 1;  % Przedział (-1,1)
        u = rand;
        f_x = sqrt(1 - x^2);  % Półokrąg
        if u <= f_x
            accept = true;
            samples(i) = x;
        end
    end
end

figure;
histogram(samples, edges, 'Normalization', 'pdf');
hold on;

% Teoretyczny rozkład półokręgu
x = linspace(-1, 1, 1000);
y = (2/pi) * sqrt(1 - x.^2);  % Normalizacja przez (2/pi) dla pełnej gęstości prawdopodobieństwa

plot(x, y, 'r', 'LineWidth', 2);

title('Rozkład półokręgu');
xlabel('Wartości');
ylabel('Gęstość prawdopodobieństwa');
legend('Histogram', 'Funkcja teoretyczna');
grid on;
saveas(gcf, 'Rozklad_Zadanie3.png');

%% Zadanie 4: Generacja liczb z rozkładu normalnego metodą odrzucania
lambda = 1; % Współczynnik rozkładu wykładniczego
samples = zeros(N,1);

for i = 1:N
    accept = false;
    while ~accept
        y = -log(rand) / lambda; % Rozkład wykładniczy
        u = rand;
        if u <= exp(-0.5 * (y - 1)^2)
            accept = true;
            if rand < 0.5
                samples(i) = y;
            else
                samples(i) = -y;
            end
        end
    end
end

figure;
histogram(samples, edges, 'Normalization', 'pdf');
hold on;

% Teoretyczny rozkład normalny
x = linspace(-4, 4, 1000);
y = (1 / sqrt(2 * pi)) * exp(-0.5 * x.^2);
plot(x, y, 'r', 'LineWidth', 2);

title('Rozkład normalny N(0,1)');
xlabel('Wartości');
ylabel('Gęstość prawdopodobieństwa');
legend('Histogram', 'Funkcja teoretyczna');
grid on;
saveas(gcf, 'Rozklad_Zadanie4.png');

