
%  Nome: kepler_nr.m
%  Número da versão actual: 1.0
%  Número do grupo: 44
%  Nomes dos autores: João Caires, nº113627
%                     Afonso Cardoso, nº113555
%                     Eduardo Monteiro, nº114191
%  Data da primeira versão: 05/Nov/2026
%  Breve descrição:
%     Resolve a equação de Kepler elíptica, f(E) = E - e*sin(E) - M = 0,
%     recorrendo ao método iterativo de Newton-Raphson, devolvendo a
%     aproximação final da anomalia excêntrica E, o número de iterações
%     realizadas, um indicador lógico de convergência, o histórico completo
%     das aproximações sucessivas (E_0, E_1, ..., E_k) e o resíduo final.
%  Lista de algoritmos utilizados e referências:
%     - Método de Newton-Raphson aplicado à equação de Kepler elíptica,
%       conforme formulado em:
%       [1] P. J. S. Gil, "Elementos de Mecânica Orbital", IST, 2015
%           (vs. 1.5.1, 2025-09-08).
%       [4] J. M. A. Danby, "Fundamentals of Celestial Mechanics", 2nd ed.,
%           Willmann-Bell, 1988.
%     - Alterações efectuadas relativamente ao algoritmo base:
%       1. Dupla condição de paragem explícita: |E_{k+1} - E_k| < tol e
%          |f(E_{k+1})| < tol.
%       2. Protecção contra divisão por derivada numericamente nula
%          (|1 - e*cos(E_k)| < tol_deriv) e contra valores NaN/Inf.
%       3. Preservação do tipo numérico de E0 (double ou single) para
%          permitir testes em precisão simples (alínea 3.c).
%  Intervalos de validade esperados e efectivamente testados:
%     - Intervalos de validade esperados: 0 <= e < 1 (órbitas elípticas e
%       circulares), M em R (rad), E0 em R (rad), tol > 0, Nmax >= 1.
%     - Intervalos efectivamente testados: e em [0.5, 0.995715],
%       M em [0.001635, 3.0] rad, E0 em [0, 2.24] rad,
%       tol em [1e-12, 1e-6], Nmax em [20, 100].
%  Lista de inputs:
%     - e    : Excentricidade orbital (adimensional, 0 <= e < 1)
%     - M    : Anomalia média (rad)
%     - E0   : Estimativa inicial da anomalia excêntrica (rad)
%     - tol  : Tolerância numérica de paragem, epsilon (rad)
%     - Nmax : Número máximo de iterações permitido (inteiro positivo)
%  Lista de outputs:
%     - E         : Aproximação final da anomalia excêntrica E_k (rad)
%     - n_iter    : Número de iterações efectivamente realizadas
%     - convergiu : Indicador booleano (true = convergiu; false = falha)
%     - E_seq     : Vector linha [E_0, E_1, ..., E_k] com a sequência (rad)
%     - r         : Resíduo final r = |E - e*sin(E) - M| (rad)
%  Lista de dependências:
%     Nenhuma dependência externa (usa apenas funções elementares nativas
%     do MATLAB: sin, cos, abs, eps, isnan, isinf).
%  Versão do software / ambiente de desenvolvimento e teste:
%     MATLAB R2024b / R2025a
%  Lista de alterações efectuadas e respectivas datas:
%     - 05/Nov/2026 (v1.0): Versão inicial desenvolvida para o TPC-1.
%     [Espaço reservado para registo de futuras alterações após a entrega]

function [E, n_iter, convergiu, E_seq, r] = kepler_nr(e, M, E0, tol, Nmax)

% Validação elementar dos parâmetros de entrada
if e < 0 || e >= 1
    error('kepler_nr:ExcentricidadeInvalida', ...
        'A excentricidade deve satisfazer 0 <= e < 1.');
end
if tol <= 0 || Nmax < 1
    error('kepler_nr:ParametrosInvalidos', ...
        'A tolerância deve ser positiva e Nmax >= 1.');
end

% Preservar o tipo de dados (double ou single) definido por E0
tipo = class(E0);
e_val = cast(e, tipo);
M_val = cast(M, tipo);
Ek    = cast(E0, tipo);

% Limiar de protecção contra divisão por valor numericamente nulo
tol_deriv = cast(100, tipo) * eps(cast(1, tipo));

% Pré-alocação do vector de histórico (inclui E_0 na posição 1)
E_seq = zeros(1, Nmax + 1, tipo);
E_seq(1) = Ek;

n_iter = 0;

% Ciclo iterativo de Newton-Raphson
for k = 1:Nmax
    % Cálculo da função e da derivada em E_k
    f_k  = Ek - e_val * sin(Ek) - M_val;
    df_k = cast(1, tipo) - e_val * cos(Ek);

    % Protecção contra divisão por derivada numericamente nula ou NaN/Inf
    if abs(df_k) < tol_deriv || isnan(df_k) || isinf(df_k)
        convergiu = false;
        E_seq = E_seq(1:k);
        E = Ek;
        r = abs(f_k);
        return;
    end

    % Passo de Newton-Raphson: E_{k+1} = E_k - f(E_k)/f'(E_k)
    passo = f_k / df_k;
    Ek1   = Ek - passo;

    % Protecção contra geração de NaN ou Inf
    if isnan(Ek1) || isinf(Ek1)
        convergiu = false;
        E_seq = E_seq(1:k);
        E = Ek;
        r = abs(f_k);
        return;
    end

    % Registo da nova aproximação e actualização do contador
    n_iter = k;
    E_seq(k + 1) = Ek1;

    % Avaliação explícita do resíduo em E_{k+1}
    f_k1 = Ek1 - e_val * sin(Ek1) - M_val;

    % Verificação explícita das duas condições de paragem (eq. 4)
    if (abs(Ek1 - Ek) < tol) && (abs(f_k1) < tol)
        convergiu = true;
        E_seq = E_seq(1:(k + 1));
        E = Ek1;
        r = abs(f_k1);
        return;
    end

    % Preparar próxima iteração
    Ek = Ek1;
end

% Se o ciclo terminou sem atingir o critério de convergência em Nmax
E_seq = E_seq(1:(n_iter + 1));
E = Ek;
r = abs(E - e_val * sin(E) - M_val);
convergiu = false;
