clc; clear; close all;

%% Generator piłokształtny (Zadanie 1)
X0 = 0.123456789; %wartosc poczatkowa
N = 10000; %Ilość iteracji
edges = 100;

z_values = [3, 15, 30];  % Poprawione - teraz jest to zmienna

% Układ wykresów z większymi odstępami
figure;
tiledlayout(3, 1, 'Padding', 'loose', 'TileSpacing', 'loose');

for i = 1:3
    z = z_values(i);  % Indeksujemy tablicę zmiennych z poprawnie
    Xn = X0;
    A = zeros(1,N);
    
    for j = 1:N
        A(j) = Xn;
        Xn = Xn * z - floor(Xn * z);
    end
    
    nexttile;
    histogram(A, edges);
    title(['Generator piłokształtny, z = ', num2str(z)]);
    xlabel('Wygenerowane liczby');
    ylabel('Częstość występowania');
    grid on;
end

saveas(gcf, 'Generator_piloksztaltny.png');



%% Generator rekurencyjny (Zadanie 2)
X0 = 0.31;
N = 1000; % Liczba generowanych próbek
edges = 100; % Liczba przedziałów histogramu

% Różne zestawy parametrów
params = [
    111, 777, 555, 2^32;
    111, 777, 0, 2^32;
    111, 777, 555, 2^28]

% Układ wykresów z większymi odstępami
figure;
tiledlayout(3, 1, 'Padding', 'loose', 'TileSpacing', 'loose');

for i = 1:size(params, 1)
    % Ustawianie parametrów dla danego przypadku
    a0 = params(i, 1);
    a1 = params(i, 2);
    c = params(i, 3);
    m = params(i, 4);

    % Inicjalizacja zmiennych
    Xn = floor(X0 * m); 
    Xn1 = floor(0.4353453 * m); 

    B = zeros(1, N);
    for j = 1:N
        B(j) = Xn / m; 
        Xn2 = mod(a0 * Xn + a1 * Xn1 + c, m); 
        Xn1 = Xn; 
        Xn = Xn2; 
    end

    % Wykres dla danego zestawu parametrów
    nexttile;
    histogram(B, edges, 'Normalization', 'pdf');
    title(['Generator rekurencyjny: a0=', num2str(a0), ', a1=', num2str(a1), ', c=', num2str(c), ', m=', num2str(m)]);
    xlabel('Wygenerowane liczby');
    ylabel('Gęstość prawdopodobieństwa');
    grid on;
end

saveas(gcf, 'Generator_rekurencyjny_porownanie.png');
