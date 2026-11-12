clear; clc; close all;

%% 1. Dados do Boeing 787
area_asa = 325;        % m^2
m = 180000;            % kg
gravidade = 9.80665;       % m/s^2
T0 = 288.15;           % K
Taxa_L = 0.0065;      % K/m 
P0 = 101325;           % Pa
M  = 0.0289644;        % kg/mol
R  = 287.05;           % J/(kg*K)
Rho0 = 1.225;           %
h = 11000;               % m
phi0 = 0; 

% cálculo da densidade do ar

massa_min = 126000;    % kg
massa_max = 253000;    % kg

Rho = @(h) Rho0 * ( 1 - Taxa_L * h / T0 ) ^ (gravidade / (R * Taxa_L) - 1);


%% 2. Condicoes de 

V_vento = [0; 
           20];           % m/s

V_aviao = [100;
           0];

% phi = arctan(Vaz / Vax) 

V_aparente = V_vento - V_aviao;

V_norm = norm(V_aparente);

phi  = pi - atan2(V_aparente(2),V_aparente(1))

% Pressao dinamica
q = 0.5 * Rho(h) * V_norm^2;

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
eq_horizontal = @(T, theta) T*cos(theta) + L * sin(phi) - D(theta); % Vx 

% Equilibrio vertical
eq_vertical = @(T, theta, m) T*sin(theta) + L(theta) + D * sin(phi) - m*gravidade; % Vz 

% Sistema das duas equacoes
% x(1) = T
% x(2) = theta

sistema = @(x, m) [
    eq_horizontal(x(1), x(2));
    eq_vertical(x(1), x(2), m)
    ];

x0 = [100000; 0.05];
solucao = fsolve(@(x) sistema(x, m), x0)

T = solucao(1)
tetha = solucao(2)

