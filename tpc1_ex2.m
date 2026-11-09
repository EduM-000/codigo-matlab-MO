%% TPC1_EX2  Mecanica Orbital 2026/27 - TPC-1, Problema 2 (Problema de Lambert)
%
% (a) Nome:        tpc1_ex2
% (b) Versao:      1.0
% (c) Autores:     [nomes dos elementos do grupo]
% (d) Data da 1.a versao: 2026-11-08
% (e) Descricao:   Script que executa todos os testes do Problema 2 com o
%                  Lambert solver lambert.m de R. Oldenhuis (MATLAB Central
%                  File Exchange #26348, v1.4; ref. [17] do enunciado).
% (f) Algoritmos:  Lambert: Izzo (secante em x, com fallback para Lancaster-
%                  Blanchard com melhorias de Gooding), via lambert.m [17];
%                  lambert_tol.m = copia de lambert.m com tolerancia e n.o
%                  maximo de iteracoes de Izzo como entradas (pergunta 3).
%                  Orbita conhecida: equacao da orbita + equacao de Kepler
%                  (EMO). Verificacao: formulas da corda com cot(alpha/2) e
%                  cot(beta/2) (EMO; Battin 1999).
% (g) Validade:    orbitas elipticas; testado com chi = 0.635, e em
%                  [0.227, 0.89525], M = 0 e 1, transferencias 42.7 a 270 graus
%                  e quase 180 graus.
% (h) Inputs:      chi, mu (parametros do enunciado).
% (i) Outputs:     tabelas na consola; figura das 4 orbitas (pergunta 5).
% (j) Dependencias: lambert.m, lambert_tol.m, lambert_dir.m, lambert_chord.m,
%                  orbit_state.m, tof_kepler.m, rv2elem.m (mesma pasta).
% (k) Ambiente:    MATLAB R2024b; testado tambem em GNU Octave 8.4.
% (l) Historico de alteracoes: v1.0 2026-11-08 versao entregue.

clear; clc; format long g;
chi = 0.635;            % chi do grupo (exato)
mu  = 3.986e5;          % mu Terra [km^3/s^2]
Re  = 6378;             % raio da Terra [km]

%% ------------------------------------------------------------ Pergunta 2
aA = 7600 + 3000*chi;   eA = 0.10 + 0.20*chi;
th1 = 30*pi/180;        th2 = (60 + 20*chi)*pi/180;
fprintf('Orbita A: a = %.10g km, e = %.10g, th1 = %.10g deg, th2 = %.10g deg\n', ...
        aA, eA, th1*180/pi, th2*180/pi);

% 2(a) solucao Kepleriana
[r1, v1] = orbit_state(aA, eA, th1, mu);
[r2, v2] = orbit_state(aA, eA, th2, mu);
dtA = tof_kepler(aA, eA, th1, th2, mu);
fprintf('\n2(a)\n r1 = [%.10f %.10f %.1f] km\n r2 = [%.10f %.10f %.1f] km\n', r1, r2);
fprintf(' v1 = [%.12f %.12f %.1f] km/s\n v2 = [%.12f %.12f %.1f] km/s\n', v1, v2);
fprintf(' dt = %.10f s\n', dtA);

% 2(b) Lambert, M = 0, directo
[v1L, v2L, flag] = lambert_dir(r1, r2, dtA, 0, +1, +1, mu);
fprintf('\n2(b) exitflag = %d\n v1L = [%.12f %.12f %.1f] km/s\n v2L = [%.12f %.12f %.1f] km/s\n', flag, v1L, v2L);
fprintf(' |v1L - v1| = %.3e km/s (relativo %.3e)\n', norm(v1L - v1), norm(v1L - v1)/norm(v1));
fprintf(' |v2L - v2| = %.3e km/s (relativo %.3e)\n', norm(v2L - v2), norm(v2L - v2)/norm(v2));

% 2(c) elementos a partir de (r1, v1L) e (r2, v2L)
[aL, eL, t1L] = rv2elem(r1, v1L, mu);
[aL2, eL2, t2L] = rv2elem(r2, v2L, mu);
fprintf('\n2(c) a = %.10f km (erro %.3e km)\n     e = %.12f (erro %.3e)\n', aL, aL - aA, eL, eL - eA);
fprintf('     th1 = %.10f deg (erro %.3e deg)\n     th2 = %.10f deg (erro %.3e deg)\n', ...
        t1L*180/pi, (t1L - th1)*180/pi, t2L*180/pi, (t2L - th2)*180/pi);
fprintf('     (a, e a partir de r2, v2L: %.10f km, %.12f)\n', aL2, eL2);

% 2(d) verificacao pelas formulas da corda (com a = aA, imposto)
[v1c, v2c, al, be, rt] = lambert_chord(r1, r2, dtA, aA, 0, mu);
c = norm(r2 - r1); s = (norm(r1) + norm(r2) + c)/2;
fprintf('\n2(d) c = %.10f km, s = %.10f km, alpha = %.12f rad, beta = %.12f rad\n', c, s, al, be);
fprintf('     residuo da eq. do tempo de Lagrange = %.3e s\n', rt);
fprintf('     v1c = [%.12f %.12f] km/s, |v1c - v1L| = %.3e km/s\n', v1c(1:2), norm(v1c - v1L));
fprintf('     v2c = [%.12f %.12f] km/s, |v2c - v2L| = %.3e km/s\n', v2c(1:2), norm(v2c - v2L));
[v1c2, v2c2] = lambert_chord(r1, r2, dtA, aL, 0, mu);      % com a do solver
fprintf('     (com a = aL do solver: |v1c-v1L| = %.3e, |v2c-v2L| = %.3e km/s)\n', ...
        norm(v1c2 - v1L), norm(v2c2 - v2L));

% 2(e) problema inverso no tempo: r2 -> r1, mesmo dt
[w1r, w2r, fr] = lambert_dir(r2, r1, dtA, 0, -1, +1, mu);   % retrogrado
[w1d, w2d, fd] = lambert_dir(r2, r1, dtA, 0, +1, +1, mu);   % directo (long way)
fprintf('\n2(e) retrogrado: flag %d, v1 = [%.12f %.12f %.1f], v2 = [%.12f %.12f %.1f]\n', fr, w1r, w2r);
fprintf('     |v1 + v2L| = %.3e, |v2 + v1L| = %.3e km/s\n', norm(w1r + v2L), norm(w2r + v1L));
[ad, ed] = rv2elem(r2, w1d, mu);
fprintf('     directo:    flag %d, v1 = [%.12f %.12f %.1f], v2 = [%.12f %.12f %.1f]\n', fd, w1d, w2d);
fprintf('                 a = %.6f km, e = %.8f\n', ad, ed);

%% ------------------------------------------------------------ Pergunta 3
fprintf('\n3(b) Orbita A, tolerancia de Izzo variavel (lambert_tol.m)\n');
fprintf('%8s %11s %11s %11s %11s %5s  %s\n', 'tol', 'err v1', 'err v2', 'err a', 'err e', 'iter', 'algoritmo');
tols = [1e-2 1e-4 1e-6 1e-8 1e-10 1e-12 1e-14 1e-15 1e-16];
for tol = tols
    [u1, u2, fl, inf3] = lambert_dir(r1, r2, dtA, 0, +1, +1, mu, tol);
    [aT, eT] = rv2elem(r1, u1, mu);
    fprintf('%8.0e %11.3e %11.3e %11.3e %11.3e %5d  %s\n', tol, norm(u1 - v1), ...
        norm(u2 - v2), abs(aT - aA), abs(eT - eA), inf3.iter, inf3.alg);
end

%% ------------------------------------------------------------ Pergunta 4
rpA = aA*(1 - eA);
eB = 0.80 + 0.15*chi;   aB = rpA/(1 - eB);
fprintf('\n4    rpA = %.10f km; Orbita B: a = %.10f km, e = %.10g\n', rpA, aB, eB);
[R1, V1] = orbit_state(aB, eB, th1, mu);
[R2, V2] = orbit_state(aB, eB, th2, mu);
dtB = tof_kepler(aB, eB, th1, th2, mu);
fprintf('4(a) r1 = [%.10f %.10f 0] km, r2 = [%.10f %.10f 0] km\n', R1(1:2), R2(1:2));
fprintf('     v1 = [%.12f %.12f 0] km/s, v2 = [%.12f %.12f 0] km/s\n', V1(1:2), V2(1:2));
fprintf('     dt = %.10f s\n', dtB);
[V1L, V2L, flB] = lambert_dir(R1, R2, dtB, 0, +1, +1, mu);
[aBL, eBL] = rv2elem(R1, V1L, mu);
fprintf('4(b) flag %d, v1L = [%.12f %.12f], v2L = [%.12f %.12f] km/s\n', flB, V1L(1:2), V2L(1:2));
fprintf('     err v1 = %.3e (rel %.3e), err v2 = %.3e (rel %.3e) km/s\n', ...
    norm(V1L - V1), norm(V1L - V1)/norm(V1), norm(V2L - V2), norm(V2L - V2)/norm(V2));
fprintf('     err a = %.3e km (rel %.3e), err e = %.3e\n', aBL - aB, (aBL - aB)/aB, eBL - eB);
fprintf('     (orbita A: err a rel = %.3e, err e = %.3e)\n', (aL - aA)/aA, eL - eA);
[~, ~, ~, infB] = lambert_dir(R1, R2, dtB, 0, +1, +1, mu, 1e-14);
[~, ~, ~, infA] = lambert_dir(r1, r2, dtA, 0, +1, +1, mu, 1e-14);
fprintf('     iteracoes de Izzo: orbita A %d, orbita B %d\n', infA.iter, infB.iter);

%% ------------------------------------------------------------ Pergunta 5
P1 = [7000; 0; 0]; P2 = [0; 8000; 0]; dt5 = 10000;
fprintf('\n5(a,b) M = 1\n');
nomeS = {'directo', 'retrogrado'}; sent = [+1 -1]; ramos = [+1 -1];
sol = struct('v1',{},'v2',{},'a',{},'e',{},'lab',{});
for i = 1:2
    for j = 1:2
        [q1, q2, fl] = lambert_dir(P1, P2, dt5, 1, sent(i), ramos(j), mu);
        [a5, e5] = rv2elem(P1, q1, mu);
        [~, ~, ~, ~, rt5, x5] = lambert_chord(P1, P2, dt5, a5, 1, mu, sent(i));
        lab = sprintf('%s, m = %+d', nomeS{i}, ramos(j));
        fprintf(' %-18s flag %d  v1 = [%9.6f %9.6f] v2 = [%9.6f %9.6f] km/s\n', lab, fl, q1(1:2), q2(1:2));
        fprintf('                    a = %.6f km, e = %.8f, rp = %.3f km, x = %+.6f (res. t = %.1e s)\n', ...
            a5, e5, a5*(1 - e5), x5, rt5);
        sol(end+1) = struct('v1', q1, 'v2', q2, 'a', a5, 'e', e5, 'lab', lab); %#ok<SAGROW>
    end
end
fprintf('\n5(c) M = 0\n');
for i = 1:2
    for j = 1:2
        [q1, q2, fl] = lambert_dir(P1, P2, dt5, 0, sent(i), ramos(j), mu);
        [a5, e5] = rv2elem(P1, q1, mu);
        fprintf(' %-10s (sinal %+d) flag %d v1 = [%9.6f %9.6f] v2 = [%9.6f %9.6f] a = %.4f km, e = %.6f\n', ...
            nomeS{i}, ramos(j), fl, q1(1:2), q2(1:2), a5, e5);
    end
end
% tempo minimo para M = 1 (minimo de t(x) na eq. de Lagrange) em cada sentido
for i = 1:2
    c5 = norm(P2 - P1); s5 = (norm(P1) + norm(P2) + c5)/2; am = s5/2;
    dth = mod(sent(i)*atan2(P1(1)*P2(2) - P1(2)*P2(1), dot(P1, P2)), 2*pi);
    xs = linspace(-0.999999, 0.999999, 200001); aa = am./(1 - xs.^2);
    alf = 2*acos(xs); bet = 2*asin(sqrt((s5 - c5)/2./aa)); if dth > pi, bet = -bet; end
    tx = sqrt(aa.^3/mu).*(2*pi + (alf - sin(alf)) - (bet - sin(bet)));
    [tmin, k] = min(tx);
    fprintf(' %s: tempo minimo para M = 1: T_min = %.2f s (em x = %.4f)\n', nomeS{i}, tmin, xs(k));
end

% figura das quatro orbitas M = 1
figure; hold on; grid on;
tt = linspace(0, 2*pi, 400);
fill(Re*cos(tt), Re*sin(tt), [0.85 0.85 0.85], 'EdgeColor', [0.4 0.4 0.4]);   % Terra
cores   = {[0 0.35 0.8], [0 0.35 0.8], [0.85 0.2 0.1], [0.85 0.2 0.1]};
tracos  = {'-', '--', '-', '--'};
lims = [0 0 0 0];
for k = 1:4
    [~, ~, ~, ev, hv] = rv2elem(P1, sol(k).v1, mu);
    Pu = ev/norm(ev); Qu = cross(hv/norm(hv), Pu);
    E = linspace(0, 2*pi, 600); b = sol(k).a*sqrt(1 - sol(k).e^2);
    X = sol(k).a*(cos(E) - sol(k).e)'*Pu.' + (b*sin(E))'*Qu.';
    plot(X(:,1), X(:,2), 'LineStyle', tracos{k}, 'Color', cores{k}, 'LineWidth', 1.5);
    lims = [min(lims(1), min(X(:,1))), max(lims(2), max(X(:,1))), ...
            min(lims(3), min(X(:,2))), max(lims(4), max(X(:,2)))];
end
axis equal; axis(lims + 500*[-1 1 -1 1]);
plot([P1(1) P2(1)], [P1(2) P2(2)], 'ko', 'MarkerFaceColor', 'k');
legend(['Terra', {sol.lab}], 'Location', 'southoutside');
xlabel('x [km]'); ylabel('y [km]');
title('Pergunta 5: as quatro solucoes com M = 1 (dt = 10000 s)');

%% ------------------------------------------------------------ Pergunta 6
raA = aA*(1 + eA); TA = 2*pi*sqrt(aA^3/mu);
Q1 = [rpA; 0; 0]; Q2 = [-raA; 0; 0];
fprintf('\n6(a) rp = %.10f km, ra = %.10f km, TA/2 = %.10f s\n', rpA, raA, TA/2);
fprintf('     |r1 x r2| = %.3e km^2\n', norm(cross(Q1, Q2)));
[z1, z2, fz] = lambert_dir(Q1, Q2, TA/2, 0, +1, +1, mu);
fprintf('     solver: flag %d, v1 = [%g %g %g], v2 = [%g %g %g]\n', fz, z1, z2);
[~, vp] = orbit_state(aA, eA, 0, mu); [~, va] = orbit_state(aA, eA, pi, mu);
fprintf('     (orbita A: |vp| = %.10f km/s, |va| = %.10f km/s)\n', norm(vp), norm(va));

th2b = pi - 1e-3;
[S1, U1] = orbit_state(aA, eA, 0, mu);
[S2, U2] = orbit_state(aA, eA, th2b, mu);
dt6 = tof_kepler(aA, eA, 0, th2b, mu);
fprintf('\n6(b) r2 = [%.10f %.10f 0] km, dt = %.10f s\n', S2(1:2), dt6);
fprintf('     |r1 x r2| = %.6e km^2, angulo de transferencia = %.10f deg\n', ...
    norm(cross(S1, S2)), atan2(norm(cross(S1, S2)), dot(S1, S2))*180/pi);
[y1, y2, fy] = lambert_dir(S1, S2, dt6, 0, +1, +1, mu);
[a6, e6] = rv2elem(S1, y1, mu);
fprintf('     flag %d, v1L = [%.12f %.12f %.3e], v2L = [%.12f %.12f %.3e]\n', fy, y1, y2);
fprintf('     err v1 = %.3e, err v2 = %.3e km/s; err a = %.3e km, err e = %.3e\n', ...
    norm(y1 - U1), norm(y2 - U2), a6 - aA, e6 - eA);
