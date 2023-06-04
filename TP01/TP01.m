clear all;
close all;
clc;

%% Page 1 : les differentes commandes
syms f(x) x ;
f(x) = (1/2) * (x^2) - sin(x);

vpa(f(1))
vpa(f(1), 2)
vpa(f(1), 3)
vpa(f(1), 5)

gradient(f(x))
diff(f(x))
taylor(f(x))
taylor(f(x),x,pi/2,'Order',3)
%subs(f, x, y)


%% Page 2 : Le script fait la methode de Steepest Sescent
clear all;
syms x1 x2 x3 f(x1, x2, x3) alpha;
f(x1, x2, x3)=(x1-4)^4+ (x2-3)^2+ 4*(x3+5)^4;
x=[4; 2; -1];

g=gradient(f(x1, x2, x3));
grad=vpa(subs(g, [x1, x2, x3], x.'));

a=0.004200;
phi=subs(f(x1, x2, x3), [x1;x2;x3],x-alpha*grad);
Q_1=vpa(subs(diff(phi),alpha,a));
Q_2=vpa(subs(diff(diff(phi)),alpha,a));
a=a-Q_1/Q_2;
x=x-a*grad;

%% Page 3/1 : Résolution par méthode de Newton
clear all;
syms x f(x) q(x)  ;
f(x) = 1/2 * x^2 - sin(x);
x0 = 0.5;
f_1 = diff(f(x));
f_2 = diff(diff(f(x)));
q(x) = f(x0)+f_1*(x-x0) + 1/2 * f_2*(x-x0);
xk = x0 - subs(f_1,x,x0)/subs(f_2,x,x0);
err = abs(xk - x0);
xknew = 0;
while err >= 10^(-5)
    xknew = xk - subs(f_1,x,xk)/subs(f_2,x,xk);
    err = abs(xknew - xk);
    xk = xknew;
end

%% Page 3/2 : Résolution par méthode de Sécante
clear all;
syms x g(x);

g(x) = (2*x-1)^2 + 4*(4-1024*x)^4;
a = 0.0042 ;
b = 0.004;
diffg = diff(g(x));
err = 1;
while err >=10^(-4)
    g_1 = subs(diffg,x,b);
    g_2 = (subs(diffg,x,b) - subs(diffg,x,a))/(b - a);
    xkn = b - g_1/g_2;
    err = abs(xkn - b);
    b = xkn;
    a = b;
end

%% Page 3/3 : Résolution par méthode de Steepest Descent
clear all;
syms x1 x2 x3 f(x1, x2, x3) alpha;
f(x1, x2, x3)=(x1-4)^4+ (x2-3)^2+ 4*(x3+5)^4;
x=[4; 2; -1];
g=gradient(f(x1, x2, x3));
err = 1;
aerr = 1;
while err >=10^(-5)
    grad=vpa(subs(g, [x1, x2, x3], x.'));
    a=0.004200;
    phi=subs(f(x1, x2, x3), [x1;x2;x3],x-alpha*grad);
    while aerr >=10^(-5)
        Q_1=vpa(subs(diff(phi),alpha,a));
        Q_2=vpa(subs(diff(diff(phi)),alpha,a));
        na=a-Q_1/Q_2;
        aerr = abs(na - a);
        a = na;
    end
    nx = x - a * grad;
    err = abs(nx - x);
    x = nx;
end
%% Page 3/4 : Résolution par méthode du gradient conjugué
clear all;
syms x1 x2 x3 f(x1,x2,x3) alpha ;
f(x1, x2, x3) = (3/2)*x1^2 + 2*x2^2 + (3/2)*x3^2 + x1*x3 + 2*x2*x3 - 3*x1 -x3;
x = [0;0;0];
g=gradient(f(x1, x2, x3));
err = 1;
aerr = 1;
nx = x;
bita = 0.9;
while err >=10^(-5)
    g = -subs(g,[x1, x2, x3],nx.') + bita*subs(g,[x1, x2, x3],x.');
    grad=vpa(subs(g, [x1, x2, x3], x.'));
    a=0.004200;
    phi=subs(f(x1, x2, x3),[x1;x2;x3],x-alpha*grad);
    while aerr >=10^(-5)
        Q_1=vpa(subs(diff(phi),alpha,a));
        Q_2=vpa(subs(diff(diff(phi)),alpha,a));
        na=a-Q_1/Q_2;
        aerr = abs(na - a);
        a = na;
    end
    nx = x - a * grad;
    err = abs(nx - x);
    x = nx;
end

