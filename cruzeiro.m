
clear; clc; close all;

%% 1. Dados do Boeing 787
area_asa = 325;        % m^2
massa_min = 126000;    % kg
massa_max = 253000;    % kg
gravidade = 9.8;       % m/s^2

%% 2. Condicoes de cruzeiro
rho = 0.3639;          % kg/m^3 (11 km)
V = 266.111;           % m/s

% Pressao dinamica
q = 0.5 * rho * V^2;

%% 3. Coeficientes aerodinamicos
% https://github.com/junzis/openap/blob/master/openap/data/dragpolar/b789.yml
% aircraftinvestigation.info/airplanes/Boeing_787-9.html
CL0 = 0.16;
CLalpha = 5.12;        % 1/rad

% Arrasto (OpenAP)
% github.com/junzis/openap
CD0 = 0.022;
k = 0.042;

% Assumindo alpha = theta
CL = @(theta) CL0 + CLalpha * theta;
CD = @(theta) CD0 + k * CL(theta).^2;

%% 4. Forcas aerodinamicas
L = @(theta) q * area_asa * CL(theta);
D = @(theta) q * area_asa * CD(theta);

%% 5. Equacoes de equilibrio

% Equilibrio horizontal
eq_horizontal = @(T, theta) T*cos(theta) - D(theta);

% Equilibrio vertical
eq_vertical = @(T, theta, m) T*sin(theta) + L(theta) - m*gravidade;

% Sistema das duas equacoes
% x(1) = T
% x(2) = theta

sistema = @(x, m) [
    eq_horizontal(x(1), x(2));
    eq_vertical(x(1), x(2), m)
];

%% 6. Calculo para diferentes massas

massa = linspace(massa_min, massa_max, 200);

theta = zeros(size(massa));
T = zeros(size(massa));

% Estimativa inicial [T; theta]
x0 = [100000; 0.05];

for i = 1:length(massa)

    % Resolver as duas equacoes simultaneamente
    solucao = fsolve(@(x) sistema(x, massa(i)), x0);

    T(i) = solucao(1);
    theta(i) = solucao(2);

    %x0 = solucao;

end

%% 7. Conversao de unidades

massa_ton = massa / 1000;
T_kN = T / 1000;
theta_graus = theta * 180/pi;

%% 8. Graficos

figure;

subplot(2,1,1)
plot(massa_ton, T_kN, 'b', 'LineWidth', 2)
grid on
xlabel('Massa [toneladas]')
ylabel('Empuxo T [kN]')
title('Empuxo necessario para cruzeiro')

subplot(2,1,2)
plot(massa_ton, theta_graus, 'r', 'LineWidth', 2)
grid on
xlabel('Massa [toneladas]')
ylabel('Theta [graus]')
title('Angulo necessario para cruzeiro')
