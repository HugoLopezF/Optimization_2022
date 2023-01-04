%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ACTUACIONES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R, Aut, To, Rt, To_t, gforce, mu_deg, T_min, SFC] = Actuaciones(hft, hft_t, v, CL, CD, S, Fm, v_t, Lift_t, Drag_t, Lift,...
                                                                          valores_T, valores_SFC, valores_T_af, valores_SFC_af)

    % PARÁMETROS 
    
    %   SFC = Specific Fuel Consumption [kg/kNh]
    %   R = Alcance [km]
    %   v = Velocidad [m/s]
    %   CL = Coef. sustentación
    %   CD = Coef. Resistencia
    %   Aut = Autonomía [h]
    %   T = Empuje [N]
    %   hft = Altura [ft]
    %   M = Número de Mach
    %   E = Eficiencia aerodinámica
    %   S = Superficie alar [m^2]
    %   Fm = Masa combustible [kg]
    %   v_t = Velocidad [m/s]
    %   Lift_t = Sustentación [N]
    %   Drag_t = Resistencia [N]
    %   Lift = Sustentación crucero (es el peso del avión) [N]
    %   mu_deg = Ángulo de balance de velocidad [deg] !!!!!!!!!!
    %   Rt = Radio del viraje [m]
    %   T_t = Empuje necesario [N]
    %   gforce = Fuerza g
    
    % AIRE
    
    [~, ~, ~, rho] = atmosisa(hft/3.28084); % Atmósfera estándar
    g = 9.80665; % Gravedad [kg/m*s^2]
    
    % AERONAVE %
    
    W = 0.5*rho*S*v^2*CL; % Peso aeronave [N]
    FW = Fm*g; % Peso combustible [N]
    
    
%% CÁLCULO VIRAJE PLANO HORIZONTAL %%

%%% Cálculo de fuerzas
            
    n = Lift_t/Lift; % Factor de carga
    
    mu = acos(1/n); % Ángulo de balance de velocidad      
    mu_deg = mu*(180/pi); % Ángulo de balance de velocidad [deg]
    Rt = v_t^2/(g*tan(mu)); % Radio de viraje [m]
    xidot = v_t/Rt; % Velocidad angular de viraje [rad/s]

    T_t = Drag_t; % Empuje viraje [N]
    
    gforce = v_t^2/(g*Rt); % Fuerza g del viraje
    
%% CÁLCULO ALCANCE AUTONOMÍA
    
    E = CL/CD; % Eficiencia aerodinámica
    T = 0.5*rho*S*v^2*CD; % Empuje crucero[N]
    To = T/(1-(2E-5)*hft); % Empuje crucero según altura [N]
    To_t = T_t/(1-(2E-5)*hft_t); % Empuje viraje según altura [N]
    T_min = max([To_t,To]); % Empuje mínimo necesario [N]
    
%%% Cálculo SFC

    if T_min < max(valores_T)

        x = linspace(0, max( valores_T ));
    
        coefs = polyfit( valores_T, valores_SFC, 2 ); % coefs polinomio para aproximar
        
        funcion_SFC = coefs(1)*x.^2 + coefs(2)*x + coefs(3);
        SFC_o = coefs(1)*T_min^2 + coefs(2)*T_min + coefs(3);

%         disp('No hace falta postcombustor')

%         figure
%         plot(valores_T, valores_SFC, '*')
%         hold on
%         plot(x, funcion_SFC, 'LineWidth',1)
%         title('Sin postcombustor')
%         xlabel('T [N]')
%         ylabel('SFC [kg/kNh]')
%         hold off

    else
        
        x = linspace(0, max( valores_T_af ));
    
        coefs = polyfit( valores_T_af, valores_SFC_af, 2 ); % coefs polinomio para aproximar
        
        funcion_SFC = coefs(1)*x.^2 + coefs(2)*x + coefs(3);
        SFC_o = coefs(1)*T_min^2 + coefs(2)*T_min + coefs(3);

%         disp('Es necesario postcombustor')

%         figure
%         plot(valores_T_af, valores_SFC_af, '*')
%         hold on
%         plot(x, funcion_SFC, 'LineWidth',1)
%         title('Con postcombustor')
%         xlabel('T [N]')
%         ylabel('SFC [kg/kNh]')
%         hold off

    end
    
%     if To_t > To
%     
%         SFC = SFC_o*(1-(5E-6)*hft_t); % Consumo según altura [kg/kNh]
%     
%     else
%     
%         SFC = SFC_o*(1-(5E-6)*hft); % Consumo según altura [kg/kNh]
%     
%     end
    
    SFC = SFC_o;
    
%%% Alcance

    R = (3600/(g*SFC))*(sqrt(CL)/CD)*sqrt(2/(rho*S))*(sqrt(W) - sqrt(W-FW)); % Alcance máximo [km]
    
%%% Autonomía

    Aut = (1000*E/(g*SFC))*log(W/(W-FW)); % Autonomía máxima [h]
    
end
