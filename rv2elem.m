function [a, e, theta, evec, hvec] = rv2elem(r, v, mu)
%RV2ELEM  Semi-eixo maior, excentricidade e anomalia verdadeira a partir
%   de (r, v). [a, e, theta, evec, hvec] = rv2elem(r, v, mu)
%   a pela equacao vis-viva, evec = (v x h)/mu - r/|r| (vector
%   excentricidade, EMO), theta medido de evec para r no sentido do
%   movimento (sinal dado por (evec x r).h), em [0, 2*pi).
%   ---
%   Versao 1.0 | Autores: [nomes do grupo] | 1.a versao: 2026-11-08
%   Dependencias: nenhuma
%   Ambiente: MATLAB R2024b; testado em GNU Octave 8.4
%   Alteracoes: v1.0 2026-11-08 versao inicial (registar aqui as seguintes)
r = r(:); v = v(:);
rn = norm(r);
hvec = cross(r, v);
a = 1/(2/rn - dot(v, v)/mu);
evec = cross(v, hvec)/mu - r/rn;
e = norm(evec);
theta = atan2(dot(cross(evec, r), hvec)/norm(hvec), dot(evec, r));
theta = mod(theta, 2*pi);
end
