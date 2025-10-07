clc;
clear;

// ===== 0) UŽIVATELSKÁ NASTAVENÍ =============================================
// (B) Volby identifikace (= řád modelu z B(s)/A(s))
// A0 vždy začíná 1, optimalizují se A0(2:$) a všechny koeficienty B0
A0 = [1, 0.5, 2];          // póly  ( A0(end)*s^(n) + A0(3)*s^(n-1) + A0(2)*s^(n-2) + A0(1)
B0 = [1, 1, 0.2];             // nuly  ( B0(end)*s^(n) + B0(3)*s^(n-1) + B0(2)*s^(n-2) + B0(1)
ftol = 1e-8;            // Tolerance na změnu reziduí – když se hodnota cílové funkce mezi iteracemi změní méně než ftol, algoritmus končí.
xtol = 1e-8;            // Tolerance na změnu parametrů – když se parametry x mezi iteracemi změní méně než xtol, končí.
gtol = 1e-8;            // Tolerance na gradient – pokud velikost gradientu (nebo Jacobián) klesne pod gtol, považuje se za konvergentní.
maxite = 200;           // Maximální počet iterací
file_name = 'data.sav'; // Název souboru s daty u, y a t

// ===== 1) DATA Z WORKSPACE ==================================================
load(file_name);

if ~exists("u","local") | ~exists("y","local") | ~exists("t","local") then
    error("Proměnné u, t a y musí být v workspace (Variable Browser).");
end

u = u(:); y = y(:); t = t(:);     // sloupcové vektory
N = size(u,1);                    // vektory musí být stejně dlouhé

if (size(y,1)<>N) | (size(t,1)<>N) then
    error("u, t a y musí mít stejný počet vzorků.");
end

if or(diff(t) <= 0) then
    error("Vektor t musí být striktně rostoucí.");
end

// ===== 2) PŘIPRAVA PARAMETRŮ K IDENTIFIKACI =================================

if A0(1)<>1 then
    error("A0 musí začínat 1 (monický jmenovatel).");   // Jmenovatel A0 zůstává konstantní pro jednoznačnost optimalizace
end

na = length(A0) - 1;              // počet a-ček (bez A(1))
nb = length(B0);                  // počet b-ček
theta0 = [A0(2:$), B0]';          // vektor počátečních parametrů

// ===== 3) REZIDUÁLNÍ FUNKCE PRO OPTIMALIZACI (lsqrsolve) ====================

function e = resid_fun(theta, m)
    a_tail = theta(1:na);                 // a1..a_na
    b      = theta(na+1:na+nb);           // b0..b_(nb-1)

    A_id = [1; a_tail];                      // [1, a1, a2, ...] jmenovatel
    B_id = b';                                // [b0, b1, ...] čitatel
    
    // polynomy z koeficientů (s^n ... s^0)
    Apoly = poly(A_id, 's', 'coeff');
    Bpoly = poly(B_id, 's', 'coeff');

    G = syslin('c', Bpoly, Apoly);        // spojitý systém
    yhat = csim(u', t', G);                 // simulace odezvy
    yhat = yhat(:);                       // sloupcový vektor
    
    e = y - yhat;                         // rezidualní vektor (N×1)
endfunction

// ===== 4) OPTIMALIZACE (NELINEÁRNÍ LS) =====================================

stop = [ftol, xtol, gtol, maxite, 0, 100];
[theta_opt, v, info] = lsqrsolve(theta0, resid_fun, N, stop);
if info<>1 then
    mprintf("Upozornění: lsqrsolve info=%d (konvergence nemusí být ideální)\n", info);
end

// ===== 5) REKONSTRUKCE OPTIMALIZOVAN0HO MODELU ==============================

a_tail_opt = theta_opt(1:na);
b_opt      = theta_opt(na+1:na+nb);

A_id = [1, a_tail_opt'];
B_id = b_opt';

s = poly(0,'s');
Apoly = poly(A_id, 's', 'coeff');
Bpoly = poly(B_id, 's', 'coeff');

G = syslin('c', Bpoly, Apoly);
yhat = csim(u', t', G);
yhat = yhat(:);

// ===== 6) METRIKA KVALITY ===================================================

e_full = y - yhat;
rmse   = sqrt(mean(e_full.^2));
fit    = 100*(1 - norm(e_full)/norm(y - mean(y)));
poles  = roots(Apoly);
zeros_ = roots(Bpoly);

// ===== 7) VÝSTUP A ZOBRAZENÍ ================================================

mprintf("\n============================ IDENTIFIKOVANÝ MODEL ============================\n")
mprintf("\nRMSE = %.6f,  FIT = %.2f %%,  sum(E^2) = %.4f\n", rmse, fit, sum(e_full.^2));
mprintf("Póly (s):\n"); disp(poles);
mprintf("Nuly (s):\n"); disp(zeros_);

close();
scf(); clf();
plot(t, y, t, yhat);
legend("y (měřené)", "yhat (model)");
xlabel("time [s]"); ylabel("Amplitude");
xgrid();

I = length(string(Apoly));
mprintf("\n        %s\n", string(Bpoly));
mprintf(" G = ---");
for i = 1:I
    mprintf("-")
end
printf("---\n")
mprintf("        %s\n", string(Apoly));




