function [v1, v2, flag, info] = lambert_dir(r1, r2, dt, M, sentido, ramo, mu, tol)
%LAMBERT_DIR  Interface para o Lambert solver de R. Oldenhuis (lambert.m),
%   escolhendo o SENTIDO DO MOVIMENTO em vez de "short/long way".
%   [v1, v2, flag, info] = lambert_dir(r1, r2, dt, M, sentido, ramo, mu, tol)
%   r1, r2  [km] vectores 3x1; dt [s] tempo de voo; M n.o de revolucoes
%   sentido +1 directo (h segundo +z), -1 retrogrado (h segundo -z)
%   ramo    +1 ou -1: sinal passado a m (so relevante para M > 0)
%   tol     (opcional) tolerancia de Izzo -> usa lambert_tol.m
%   No solver: tf < 0 pede a solucao "long way"; o versor do momento
%   angular e sign(tf)*(r1 x r2)/|r1 x r2|. Logo, para obter o sentido
%   pedido: sign(tf) = sentido*sign((r1 x r2)_z). O solver usa tf em DIAS.
%   ---
%   Versao 1.0 | Autores: [nomes do grupo] | 1.a versao: 2026-11-08
%   Dependencias: lambert.m e lambert_tol.m (R. Oldenhuis [17])
%   Ambiente: MATLAB R2024b; testado em GNU Octave 8.4
%   Alteracoes: v1.0 2026-11-08 versao inicial (registar aqui as seguintes)
cz = r1(1)*r2(2) - r1(2)*r2(1);
s = sign(cz); if s == 0, s = 1; end     % caso degenerado (r1 x r2 = 0)
tf = sentido*s*dt/86400;                % [dias], com sinal
m  = ramo*M;
info = struct('alg','Izzo/Lancaster (lambert.m original)','iter',NaN,'err',NaN);
if nargin < 8 || isempty(tol)
    [V1, V2, ~, flag] = lambert(r1(:).', r2(:).', tf, m, mu);
else
    [V1, V2, ~, flag, info] = lambert_tol(r1(:).', r2(:).', tf, m, mu, tol);
end
v1 = V1(:); v2 = V2(:);
end
