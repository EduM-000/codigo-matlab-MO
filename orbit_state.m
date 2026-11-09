function [r, v] = orbit_state(a, e, theta, mu)
%ORBIT_STATE  Posicao e velocidade numa orbita eliptica no plano xy,
%   movimento directo, periapside no semi-eixo +x (referencial perifocal).
%   [r, v] = orbit_state(a, e, theta, mu)
%   a [km], e [-], theta anomalia verdadeira [rad], mu [km^3/s^2]
%   r [km] e v [km/s] sao vectores coluna 3x1.
%   Formulas: equacao da orbita r = p/(1+e cos(theta)), p = a(1-e^2), e
%   v = sqrt(mu/p) [-sin(theta); e+cos(theta); 0] (EMO, componentes da
%   velocidade no referencial perifocal).
%   ---
%   Versao 1.0 | Autores: [nomes do grupo] | 1.a versao: 2026-11-08
%   Dependencias: nenhuma
%   Ambiente: MATLAB R2024b; testado em GNU Octave 8.4
%   Alteracoes: v1.0 2026-11-08 versao inicial (registar aqui as seguintes)
p = a*(1 - e^2);
rn = p/(1 + e*cos(theta));
r = rn*[cos(theta); sin(theta); 0];
v = sqrt(mu/p)*[-sin(theta); e + cos(theta); 0];
end
