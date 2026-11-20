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
           0];           % m/s

V_aviao = [250;
           0];

% phi = arctan(Vaz / Vax) 

V_aparente = V_vento - V_aviao;

V_norm = norm(V_aparente);

phi  = pi - atan2(V_aparente(2),V_aparente(1));

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

% alpha 
alpha = @(theta, phi) theta + phi; 

% Assumindo alpha = theta
CL = @(theta, phi) CL0 + CLalpha * alpha(theta,phi);
CD = @(theta, phi) CD0 + k * CL(theta,phi).^2;

%% 4. Forcas aerodinamicas
L = @(theta, phi) q * area_asa * CL(theta,phi);
D = @(theta, phi) q * area_asa * CD(theta,phi);

%% 5. Equacoes de equilibrio

% Equilibrio horizontal
eq_horizontal = @(T, theta) T*cos(theta) + L(theta, phi) * sin(phi) - D(theta, phi) * cos(phi); % Vx 

% Equilibrio vertical
eq_vertical = @(T, theta, m) T*sin(theta) + L(theta, phi) * cos(phi) + D(theta, phi) * sin(phi) - m*gravidade; % Vz 

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
theta = solucao(2)

%% 6. Equacoes da dinamica
% m * dVx/dt = T*cos(theta) + L*sin(phi) - D*cos(phi)
% m * dVz/dt = T*sin(theta) + L*cos(phi) + D*sin(phi) - m*g

V_aparente_din = @(Vx, Vz) V_vento - [Vx; Vz];
phi_din = @(Vx, Vz) atan2(V_vento(2) - Vz, Vx - V_vento(1));   % mesmo phi de cima
q_din   = @(Vx, Vz) 0.5 * Rho(h) * norm(V_aparente_din(Vx, Vz))^2;
L_din   = @(theta, Vx, Vz) q_din(Vx, Vz) * area_asa * CL(theta, phi_din(Vx, Vz));
D_din   = @(theta, Vx, Vz) q_din(Vx, Vz) * area_asa * CD(theta, phi_din(Vx, Vz));

% Forca horizontal (Vx)
Fx = @(T, theta, Vx, Vz) T*cos(theta) + L_din(theta, Vx, Vz) * sin(phi_din(Vx, Vz)) - D_din(theta, Vx, Vz) * cos(phi_din(Vx, Vz));
% Forca vertical (Vz)
Fz = @(T, theta, Vx, Vz) T*sin(theta) + L_din(theta, Vx, Vz) * cos(phi_din(Vx, Vz)) + D_din(theta, Vx, Vz) * sin(phi_din(Vx, Vz)) - m*gravidade;
% estados:  x = [Vx; Vz]
% entradas: u = [T; theta]
dinamica = @(x, u) [
    Fx(u(1), u(2), x(1), x(2)) / m;
    Fz(u(1), u(2), x(1), x(2)) / m
    ];

%% 7. Simulacao: equilibrio e depois degrau (modelo nao linear)

% Ponto de equilibrio
x_eq = V_aviao;            % [Vx; Vz] de equilibrio
u_eq = solucao;            % [T; theta] de equilibrio

% Tamanho dos degraus
dT = 0.05 * abs(T);        % degrau de 5% no empuxo
dtheta = deg2rad(1);       % degrau de 1 grau em theta

% Tempo
t_degrau = 50;             % s  (instante do degrau)
t_final  = 600;            % s
opcoes = odeset('MaxStep', 0.5);

% Entradas: equilibrio ate t_degrau, depois degrau
u_T  = @(t) u_eq + (t >= t_degrau) * [dT; 0];
u_th = @(t) u_eq + (t >= t_degrau) * [0; dtheta];

% Simulacao
[t1, x1] = ode45(@(t,x) dinamica(x, u_T(t)),  [0 t_final], x_eq, opcoes);
[t2, x2] = ode45(@(t,x) dinamica(x, u_th(t)), [0 t_final], x_eq, opcoes);

% Degrau em T
figure;
subplot(3,1,1); plot(t1, (T + dT*(t1 >= t_degrau))/1e3, 'LineWidth', 1.5); grid on; ylabel('T [kN]'); title('Degrau em T');
subplot(3,1,2); plot(t1, x1(:,1), 'LineWidth', 1.5); grid on; ylabel('V_x [m/s]');
subplot(3,1,3); plot(t1, x1(:,2), 'LineWidth', 1.5); grid on; ylabel('V_z [m/s]'); xlabel('t [s]');

% Degrau em theta
figure;
subplot(3,1,1); plot(t2, rad2deg(theta + dtheta*(t2 >= t_degrau)), 'LineWidth', 1.5); grid on; ylabel('\theta [graus]'); title('Degrau em \theta');
subplot(3,1,2); plot(t2, x2(:,1), 'LineWidth', 1.5); grid on; ylabel('V_x [m/s]');
subplot(3,1,3); plot(t2, x2(:,2), 'LineWidth', 1.5); grid on; ylabel('V_z [m/s]'); xlabel('t [s]');