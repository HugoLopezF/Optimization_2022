%                   MAIN TOCDA

clc; clear; close all
warning('off','all')

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultLineLineWidth',2)

%% Initial values

load("Estructura_datos_materiales.mat")
load("Estructura_datos_prop.mat")
    
M = 2; % Mach de vuelo para crucero
M_t = 1.25; % Mach de vuelo para giro
AOA = 2; % Ángulo de ataque para crucero [deg]
AOA_t = AOA; % Ángulo de ataque para giro [deg]
Sweep = 50; % Flecha [deg]
Wingspan = 12; % Envergadura [m]
hft = 30000; % Altitud para crucero [ft]
hft_t = 3300; % Altitud para giro [ft]
Fm = 8000; % Masa combustible [kg]
E = 75E9;   % Modulo de Young

% Restricciones de desigualdad

A = [1 0 0 0 0 0 0 0 0 0; ...
     -1 0 0 0 0 0 0 0 0 0; ...
     0 1 0 0 0 0 0 0 0 0; ...
     0 -1 0 0 0 0 0 0 0 0; ...
     0 0 0 0 1 0 0 0 0 0; ...
     0 0 0 0 -1 0 0 0 0 0; ...
     0 0 0 0 0 1 0 0 0 0; ...
     0 0 0 0 0 -1 0 0 0 0; ...
     0 0 0 0 0 0 0 0 0 1; ...
     0 0 0 0 0 0 0 0 0 -1];

b = [2.5 -1.5 1.5 -1 70 -50 25 -10 135 -70];

% Restricciones de igualdad

Aeq = [0 0 1 0 0 0 0 0 0 0; ...
       0 0 0 1 0 0 0 0 0 0 ; ...
       0 0 0 0 0 0 1 0 0 0; ...
       0 0 0 0 0 0 0 1 0 0; ...
       0 0 0 0 0 0 0 0 1 0];

beq = [AOA AOA_t hft hft_t Fm/1E3];

%$ Valores iniciales

vect_init_mono = [M, M_t, AOA, AOA_t, Sweep, Wingspan, hft, hft_t, Fm/1E3, E/1E9];

valores_T = Datos_prop.T; 
valores_T_af = Datos_prop.Taf;
valores_SFC = Datos_prop.SFC; 
valores_SFC_af = Datos_prop.SFCaf; 
valores_E = Datos_mat.moduloYoung;

% Precios
precio_motores = Datos_prop.precio; 
precio_materiales = Datos_mat.precio;
precio_combustible = 1;

%% Simulator test
%  
% % Aerodinámica
% [CL, ~, CD, Lift, Drag, F_beam, c, S, v] = SupersonicAerodynamics(M, AOA, Sweep, Wingspan, hft);
% t = 0.07*c;     % Espesor del ala
% 
% [CL_t, ~, CD_t, Lift_t, Drag_t, F_beam_t, c_t, S_t, v_t] = SupersonicAerodynamics(M_t, AOA_t, Sweep, Wingspan, hft_t);
% t_t = 0.07*c_t;     % Espesor del ala
% 
% % Estructura
% u = structure(E, Wingspan, c, t, F_beam);
% u_t = structure(E, Wingspan, c_t, t_t, F_beam_t);
% 
% % Actuaciones y propulsión
% [R, Aut, To, Rt, To_t, gforce, mu_deg, To_min, SFC] = Actuaciones(hft, hft_t, v, CL, CD, S, Fm, v_t, Lift_t, Drag_t, Lift, valores_T,...
%                                                                   valores_SFC, valores_T_af, valores_SFC_af);
% 
% % Costes
% [coste_motor, precio_mat, coste_mat, coste_comb, coste_total] = costes(valores_T, valores_T_af, precio_motores, To_min, valores_E,...
%                                                                 precio_materiales, E, t, S, precio_combustible);


%% OPTIMIZACIÓN MONO-OBJETIVO

options = optimoptions('fmincon','Display','iter','MaxIterations',6, 'Algorithm','sqp');
[X, F, exitflag, output, lambda, grad, hessian] = fmincon(@coste_monoobjetivo, vect_init_mono, A, b', Aeq, beq', [], [], [], options);

% FUNCIONES 

function coste_total = coste_monoobjetivo(x)

load("Estructura_datos_materiales.mat")
load("Estructura_datos_prop.mat")

M = x(1);
M_t = x(2);
AOA = x(3);
AOA_t = x(4);
Sweep = x(5);
Wingspan = x(6);
hft = x(7);
hft_t = x(8);
Fm = x(9);
E = x(10);

% Valores

valores_T = Datos_prop.T; 
valores_T_af = Datos_prop.Taf;
valores_SFC = Datos_prop.SFC; 
valores_SFC_af = Datos_prop.SFCaf; 
valores_E = Datos_mat.moduloYoung;

% Precios
precio_motores = Datos_prop.precio; 
precio_materiales = Datos_mat.precio;
precio_combustible = 1;

% Aerodinámica
[CL, ~, CD, Lift, ~, F_beam, c, S, v] = SupersonicAerodynamics(M, AOA, Sweep, Wingspan, hft);
t = 0.07*c;     % Espesor del ala

[~, ~, ~, Lift_t, Drag_t, F_beam_t, c_t, ~, v_t] = SupersonicAerodynamics(M_t, AOA_t, Sweep, Wingspan, hft_t);
t_t = 0.07*c_t;     % Espesor del ala

% Estructura :((
u = structure(E, Wingspan, c, t, F_beam);
u_t = structure(E, Wingspan, c_t, t_t, F_beam_t);          

% Actuaciones y propulsión
[R, Aut, ~, ~, ~, ~, ~, To_min, SFC] = Actuaciones(hft, hft_t, v, CL, CD, S, Fm*1E3, v_t, Lift_t, Drag_t, Lift, valores_T,...
                                                   valores_SFC, valores_T_af, valores_SFC_af);

% Costes
[coste_motor, ~, coste_mat, coste_comb, coste_total] = costes(valores_T, valores_T_af, precio_motores, To_min, Fm*1E3, valores_E,...
                                                             precio_materiales, E, S, precio_combustible);
                    
salidas = [Aut R SFC coste_motor, coste_mat, coste_comb];
disp('    Autonomía   Alcance   SFC')
disp(['    ' num2str(salidas(1)) ' h ' num2str(salidas(2)) ' km ' num2str(salidas(3)) ' kg/kN h '])
disp(' ')
disp('    Motor   Material   Combustible')
disp(['    ' num2str(salidas(4)) ' € ' num2str(salidas(5)) ' € ' num2str(salidas(6)) ' € '])
disp(' ')

end

