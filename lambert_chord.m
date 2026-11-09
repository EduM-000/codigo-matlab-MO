function [v1, v2, alpha, beta, res_t, x] = lambert_chord(r1, r2, dt, a, M, mu, sentido)
%LAMBERT_CHORD  Velocidades terminais de uma orbita eliptica de Lambert
%   (sentido = +1 directo [omissao] ou -1 retrogrado; angulo de
%   transferencia medido no sentido do movimento, em [0, 2pi)) a partir do
%   semi-eixo maior a, usando a corda e os angulos de Lagrange:
%     v1 = (B + A) uc + (B - A) u1,   v2 = (B + A) uc - (B - A) u2,
%     A = sqrt(mu/(4a)) cot(alpha/2),  B = sqrt(mu/(4a)) cot(beta/2),
%   uc = (r2 - r1)/c, u1 = r1/|r1|, u2 = r2/|r2|,
%   sin^2(alpha/2) = s/(2a), sin^2(beta/2) = (s - c)/(2a),
%   s = (r1 + r2 + c)/2 (semi-perimetro), c = |r2 - r1| (corda).
%   O quadrante de alpha (alpha0 ou 2pi - alpha0) e o sinal de beta
%   (negativo se o angulo de transferencia > pi) sao escolhidos de modo a
%   satisfazer a equacao do tempo de Lagrange
%     dt = sqrt(a^3/mu) [2 pi M + (alpha - sin alpha) - (beta - sin beta)];
%   res_t e o residuo [s] dessa equacao para a escolha feita;
%   x = cos(alpha/2) e a variavel de Izzo (permite identificar o ramo).
%   ---
%   Versao 1.0 | Autores: [nomes do grupo] | 1.a versao: 2026-11-08
%   Dependencias: nenhuma
%   Ambiente: MATLAB R2024b; testado em GNU Octave 8.4
%   Alteracoes: v1.0 2026-11-08 versao inicial (registar aqui as seguintes)
if nargin < 7 || isempty(sentido), sentido = 1; end
r1 = r1(:); r2 = r2(:);
n1 = norm(r1); n2 = norm(r2);
cv = r2 - r1; c = norm(cv);
s  = (n1 + n2 + c)/2;
al0 = 2*asin(sqrt(s/(2*a)));
be0 = 2*asin(sqrt((s - c)/(2*a)));
dth = mod(sentido*atan2(r1(1)*r2(2) - r1(2)*r2(1), dot(r1, r2)), 2*pi);
if dth > pi, be = -be0; else, be = be0; end
cand = [al0, 2*pi - al0];
tt = sqrt(a^3/mu)*(2*pi*M + (cand - sin(cand)) - (be - sin(be)));
[res_t, k] = min(abs(tt - dt));
alpha = cand(k); beta = be;
A = sqrt(mu/(4*a))*cot(alpha/2);
B = sqrt(mu/(4*a))*cot(beta/2);
uc = cv/c; u1 = r1/n1; u2 = r2/n2;
v1 = (B + A)*uc + (B - A)*u1;
v2 = (B + A)*uc - (B - A)*u2;
x = cos(alpha/2);
end
