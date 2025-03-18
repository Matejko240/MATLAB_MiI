%% Modelowanie i identyfikacja - Laboratorium 2
% Generacja liczb losowych metodą odrzucania

clc; clear; close all;

N = 1000; % Liczba generowanych próbek

%% Zadanie 1: Generacja liczb z rozkładu o gęstości f(x)
% Rozkład składa się z dwóch funkcji liniowych na przedziałach (-1,0] i (0,1].
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
histogram(samples, 50);
title('Histogram dla f(x)'); xlabel('Wartości'); ylabel('Liczność');

%% Zadanie 2: Wyznaczenie stałej c i implementacja generatora
% Wyznaczanie c poprzez całkowanie funkcji gęstości do 1
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
histogram(samples, 50);
title('Histogram dla zadania 2'); xlabel('Wartości'); ylabel('Liczność');

%% Zadanie 3: Generacja liczb z rozkładu półokręgu
% Gęstość prawdopodobieństwa ma kształt półokręgu
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
histogram(samples, 50);
title('Histogram dla półokręgu'); xlabel('Wartości'); ylabel('Liczność');

%% Zadanie 4: Generacja liczb z rozkładu normalnego metodą odrzucania
% Generacja liczb z rozkładu normalnego N(0,1) wykorzystując metodę odrzucania
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
histogram(samples, 50);
title('Histogram dla rozkładu normalnego'); xlabel('Wartości'); ylabel('Liczność');

%% Zadanie 5: Podstawowe własności metody odrzucania
% Analiza efektywności metody poprzez wyznaczenie współczynnika akceptacji
accepted = 0;
attempts = 0;
for i = 1:N
    attempts = attempts + 1;
    if rand <= 0.5  % Przykładowy warunek akceptacji
        accepted = accepted + 1;
    end
end
acceptance_rate = accepted / attempts;
disp(['Współczynnik akceptacji: ', num2str(acceptance_rate)]);
