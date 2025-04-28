clc; clear; close all;

N = 100000; 
edges = 100;

%% Zadanie 1: f(x) = 2x, x ∈ [0,1]
U = rand(1, N);
X1 = sqrt(U);

figure;
histogram(X1, edges, 'Normalization', 'pdf');
hold on;

% Teoretyczna gęstość prawdopodobieństwa f(x) = 2x
x = linspace(0, 1, 100);
y = 2 * x;
plot(x, y, 'r', 'LineWidth', 2);

title('Rozkład f(x) = 2x na przedziale [0,1]');
xlabel('Wartości');
ylabel('Gęstość prawdopodobieństwa');
grid on;
legend('Histogram', 'Funkcja teoretyczna');
saveas(gcf, 'Rozklad_2x.png');

%% Zadanie 2: f(x) = x+1 dla x ∈ (-1,0), -x+1 dla x ∈ [0,1)
U = rand(1, N);
X2 = zeros(1, N);
X2(U < 0.5) = sqrt(2*U(U < 0.5)) - 1;
X2(U >= 0.5) = 1 - sqrt(2*(1 - U(U >= 0.5)));

figure;
histogram(X2, edges, 'Normalization', 'pdf');
hold on;

% Teoretyczna gęstość prawdopodobieństwa
x1 = linspace(-1, 0, 100);
x2 = linspace(0, 1, 100);
y1 = x1 + 1;
y2 = -x2 + 1;
plot(x1, y1, 'r', 'LineWidth', 2);
plot(x2, y2, 'r', 'LineWidth', 2);

title('Rozkład f(x) = x+1 dla x ∈ (-1,0) oraz -x+1 dla x ∈ [0,1)');
xlabel('Wartości');
ylabel('Gęstość prawdopodobieństwa');
grid on;
legend('Histogram', 'Funkcja teoretyczna');
saveas(gcf, 'Rozklad_x1.png');

%% Zadanie 3: Rozkład wykładniczy
X3 = -log(1 - U);

figure;
histogram(X3, edges, 'Normalization', 'pdf');
hold on;

% Teoretyczna funkcja gęstości dla rozkładu wykładniczego
x = linspace(0, 6, 100);
y = exp(-x);
plot(x, y, 'r', 'LineWidth', 2);

title('Rozkład wykładniczy f(x) = e^{-x}');
xlabel('Wartości');
ylabel('Gęstość prawdopodobieństwa');
grid on;
legend('Histogram', 'Funkcja teoretyczna');
saveas(gcf, 'Rozklad_wykladniczy.png');

%% Zadanie 4: Rozkład Laplace'a
X4 = zeros(1, N);
X4(U < 0.5) = log(2 * U(U < 0.5));
X4(U >= 0.5) = -log(2 * (1 - U(U >= 0.5)));

figure;
histogram(X4, edges, 'Normalization', 'pdf');
hold on;

% Teoretyczna funkcja gęstości dla rozkładu Laplace'a
x = linspace(-6, 6, 1000);
y = 0.5 * exp(-abs(x));
plot(x, y, 'r', 'LineWidth', 2);

title('Rozkład Laplace’a');
xlabel('Wartości');
ylabel('Gęstość prawdopodobieństwa');
grid on;
legend('Histogram', 'Funkcja teoretyczna');
saveas(gcf, 'Rozklad_Laplace.png');


