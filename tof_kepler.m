function dt = tof_kepler(a, e, th1, th2, mu)
%TOF_KEPLER  Tempo de voo numa orbita eliptica entre as anomalias
%   verdadeiras th1 e th2 (movimento directo, 0 <= th1 < th2 < 2*pi).
%   dt = tof_kepler(a, e, th1, th2, mu)   [s]
%   Usa tan(E/2) = sqrt((1-e)/(1+e)) tan(theta/2) e a equacao de Kepler
%   M = E - e sin(E), com n = sqrt(mu/a^3) (EMO).
%   ---
%   Versao 1.0 | Autores: [nomes do grupo] | 1.a versao: 2026-11-08
%   Dependencias: nenhuma
%   Ambiente: MATLAB R2024b; testado em GNU Octave 8.4
%   Alteracoes: v1.0 2026-11-08 versao inicial (registar aqui as seguintes)
n  = sqrt(mu/a^3);
E1 = 2*atan2(sqrt(1-e)*sin(th1/2), sqrt(1+e)*cos(th1/2));
E2 = 2*atan2(sqrt(1-e)*sin(th2/2), sqrt(1+e)*cos(th2/2));
M1 = E1 - e*sin(E1);
M2 = E2 - e*sin(E2);
dt = mod(M2 - M1, 2*pi)/n;
end
