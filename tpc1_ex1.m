
%  Nome: tpc1_ex1.m
%  Número da versão actual: 1.0
%  Número do grupo: 42
%  Nomes dos autores: João Caires, nº113627
%                     Afonso Cardoso, nº113555
%                     Eduardo Monteiro, nº114191
%  Data da primeira versão: 05/Nov/2026
%  Breve descrição:
%     Programa de teste que executa todos os casos numéricos pedidos na
%     Pergunta 1 (Equação de Kepler) do TPC-1 para x = 0.635 (alíneas 2.a,
%     2.b, 3.b, 3.c, 3.f, 3.g e 3.h) e imprime as tabelas de resultados.
%  Lista de algoritmos utilizados e referências:
%     - Método de Newton-Raphson implementado em kepler_nr.m [1].
%     - Estimativa inicial segura E_0 = min(M + e, pi) (eq. 7 do enunciado).
%     - Estimativa inicial de Danby E_0 = M + 0.85*e [3, 4].
%  Intervalos de validade esperados e efectivamente testados:
%     - Testado para x = 0.635: Caso A (e_A = 0.5175, M_A = 0.5675*pi rad),
%       Caso B (e_B = 0.995715, M_B = 0.001635 rad) e Caso 3.c (e = 0.5,
%       M = 3.0 rad, E0 = 0 em precisão simples).
%  Lista de inputs:
%     Parâmetro de grupo x = 0.635 definido internamente no script.
%  Lista de outputs:
%     Tabelas formatadas na consola com E, n_iter, estado de convergência,
%     sequências de iteradas E_k e resíduos |f(E_k)|.
%  Lista de dependências:
%     Requer a função kepler_nr.m no mesmo directório.
%  Versão do software / ambiente de desenvolvimento e teste:
%     MATLAB R2024b / R2025a
%  Lista de alterações efectuadas e respectivas datas:
%     - 05/Nov/2026 (v1.0): Versão inicial desenvolvida para o TPC-1.
% =========================================================================


clear; clc;
format long;

% Parâmetro adimensional exacto do grupo
x = 0.635;

fprintf('=================================================================\n');
fprintf('  TPC-1 MECÂNICA ORBITAL (2026/27) - PERGUNTA 1 (x = %.3f)\n', x);
fprintf('=================================================================\n\n');

%% ------------------------------------------------------------------------
%  PERGUNTA 2: Caso regular (Caso A)
%  ------------------------------------------------------------------------
eA = 0.20 + 0.50 * x;
MA = pi * (0.25 + 0.50 * x);
E0_A = MA;
Nmax_A = 50;
tols = [1e-6, 1e-9, 1e-12];

fprintf('--- PERGUNTA 2.a: Caso A (e_A = %.6f, M_A = %.15f rad) ---\n', eA, MA);
fprintf('%-12s | %-20s | %-8s | %-10s | %-15s\n', ...
        'Tolerancia', 'E final (rad)', 'N_iter', 'Convergiu', 'Residuo r (rad)');
fprintf('-----------------------------------------------------------------------------\n');

for i = 1:length(tols)
    tol = tols(i);
    [E_f, n_it, conv, E_seq, res] = kepler_nr(eA, MA, E0_A, tol, Nmax_A);
    fprintf('%-12.0e | %-20.15f | %-8d | %-10d | %-15.6e\n', ...
            tol, E_f, n_it, conv, res);
    if tol == 1e-12
        E_seq_2b = E_seq;
    end
end

fprintf('\n--- PERGUNTA 2.b: Sequência completa E_k para tol = 1e-12 ---\n');
for k = 0:(length(E_seq_2b) - 1)
    fk = abs(E_seq_2b(k+1) - eA * sin(E_seq_2b(k+1)) - MA);
    fprintf('k = %d : E_%d = %.15f rad | |f(E_%d)| = %.6e rad\n', ...
            k, k, E_seq_2b(k+1), k, fk);
end
fprintf('\n');

%% ------------------------------------------------------------------------
%  PERGUNTA 3.b: Caso B com E0 = M_B
%  ------------------------------------------------------------------------
eB = 0.990 + 0.009 * x;
MB = 1e-3 * (1 + x);
tol_B = 1e-12;
Nmax_B = 100;

[EB_1, nit_B1, conv_B1, E_seq_B1, res_B1] = kepler_nr(eB, MB, MB, tol_B, Nmax_B);

fprintf('--- PERGUNTA 3.b: Caso B com E_0 = M_B (e_B = %.6f, M_B = %.6f rad) ---\n', eB, MB);
fprintf('Convergiu: %d | Número total de iterações: %d | E_final = %.15f rad | r = %.3e rad\n', ...
        conv_B1, nit_B1, EB_1, res_B1);
fprintf('Sequência completa (incluindo E_0 e as primeiras 5 iteradas produzidas):\n');
for k = 0:(length(E_seq_B1) - 1)
    fk = abs(E_seq_B1(k+1) - eB * sin(E_seq_B1(k+1)) - MB);
    fprintf('k = %d : E_%d = %.15f rad | |f(E_%d)| = %.6e rad\n', ...
            k, k, E_seq_B1(k+1), k, fk);
end
fprintf('\n');

%% ------------------------------------------------------------------------
%  PERGUNTA 3.c: Caso e = 0.5, M = 3.0 rad, E0 = 0 em precisão simples
%  ------------------------------------------------------------------------
e_3c    = single(0.5);
M_3c    = single(3.0);
E0_3c   = single(0.0);
tol_3c  = single(1e-12);
Nmax_3c = 20;

[E_3c, nit_3c, conv_3c, E_seq_3c, res_3c] = kepler_nr(e_3c, M_3c, E0_3c, tol_3c, Nmax_3c);

fprintf('--- PERGUNTA 3.c: e = 0.5, M = 3.0 rad, E_0 = 0 (Precisão Simples) ---\n');
fprintf('Convergiu: %d | Iterações realizadas: %d\n', conv_3c, nit_3c);
fprintf('%-4s | %-18s | %-18s\n', 'k', 'E_k (rad)', '|f(E_k)| (rad)');
fprintf('------------------------------------------------\n');
for k = 0:(length(E_seq_3c) - 1)
    Ek_val = E_seq_3c(k+1);
    fk_val = abs(Ek_val - e_3c * sin(Ek_val) - M_3c);
    fprintf('%-4d | %-18.7f | %-18.7f\n', k, Ek_val, fk_val);
end
fprintf('\n');

%% ------------------------------------------------------------------------
%  PERGUNTAS 3.f, 3.g e 3.h: Comparação de estimativas iniciais no Caso B
%  ------------------------------------------------------------------------
E0_B2 = min(MB + eB, pi);
[EB_2, nit_B2, conv_B2, E_seq_B2, res_B2] = kepler_nr(eB, MB, E0_B2, tol_B, Nmax_B);

fprintf('--- PERGUNTA 3.f: Caso B com E_0 = min(M_B + e_B, pi) = %.6f rad ---\n', E0_B2);
fprintf('Convergiu: %d | Número de iterações: %d | E_final = %.15f rad | r = %.3e rad\n', ...
        conv_B2, nit_B2, EB_2, res_B2);
for k = 0:(length(E_seq_B2) - 1)
    fprintf('k = %d : E_%d = %.15f rad\n', k, k, E_seq_B2(k+1));
end
fprintf('\n');

E0_B3 = MB + 0.85 * eB;
[EB_3, nit_B3, conv_B3, E_seq_B3, res_B3] = kepler_nr(eB, MB, E0_B3, tol_B, Nmax_B);

fprintf('--- PERGUNTA 3.g: Caso B com E_0 = M_B + 0.85*e_B = %.8f rad ---\n', E0_B3);
fprintf('Convergiu: %d | Número de iterações: %d | E_final = %.15f rad | r = %.3e rad\n', ...
        conv_B3, nit_B3, EB_3, res_B3);
for k = 0:(length(E_seq_B3) - 1)
    fprintf('k = %d : E_%d = %.15f rad\n', k, k, E_seq_B3(k+1));
end
fprintf('\n');

fprintf('--- PERGUNTA 3.h: Tabela Resumo (Caso B) ---\n');
fprintf('%-25s | %-16s | %-8s | %-18s\n', ...
        'Estimativa Inicial E_0', 'Valor E_0 (rad)', 'N_iter', 'E final (rad)');
fprintf('----------------------------------------------------------------------------\n');
fprintf('%-25s | %-16.8f | %-8d | %-18.15f\n', 'E_0 = M_B', MB, nit_B1, EB_1);
fprintf('%-25s | %-16.8f | %-8d | %-18.15f\n', 'E_0 = min(M_B + e_B, pi)', E0_B2, nit_B2, EB_2);
fprintf('%-25s | %-16.8f | %-8d | %-18.15f\n', 'E_0 = M_B + 0.85*e_B', E0_B3, nit_B3, EB_3);
