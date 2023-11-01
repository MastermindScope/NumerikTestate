%% Numerik Testat Aufgabe 1

f = 50; %Hz
T = 1/f; %s
U0 = 325; %V

% Integral x(t) dt von 0 bis T für x(t) = sin(x*f*2pi)
% (-cos(2piTf)+cos(0))/(2pif) = 0
% -> Analytisch: Arithmetischer MW = 0
%Numerisch wie folgt:

t = linspace(0,T,100000);
s = sin(2*pi*f*t)*U0;

arithMedian = 1/length(s)*sum(s);
disp('Der Arithmetische Mittelwert für eine Sinusschwingung von 325 V beläuft sich numerisch zu: ')
disp(num2str(arithMedian))

% Integral |x(t)| dt von 0 bis T für x(t) = sin(x*f*2pi)
% da gilt |sin(x)| = |sin(x+pi)| kann das Integral zu 2 + Integral
% sin(x*f*2pi) von 0 bis T/2 vereinfacht werden
% f*2pi*T/2 = pi = pi+f*2pi*T
% aus sin(x) >= 0 für 0<=x<=pi folgt |sin(x)| = sin(x) für diesen Bereich
% Gleichrichtwert für die Funktion also:
% (50 Hz * 2 * 325 V *(-cos(pi)+cos(0))/(2pi*50 Hz)) = 2 * 325 V / pi
% numerisch:

gleichricht = 1/length(s)*sum(abs(s));
disp('Der Gleichrichtwert für eine Sinusschwingung von 325 V beläuft sich numerisch zu: ')
disp(gleichricht)

% Integral x(t)^2 dt von 0 bis T für x(t) = sin(x*f*2pi)
% nach Formelsammlung Merziger I(sin(ax)^2,0,2pi) = x/2 - sin(2ax)/(4a)
% Effektivwert also analytisch:
% sqrt(U0^2*(T/2-sin(4pi)/(4*2pi*f)-0/2+sin(0)/(4*2pi*f)))/sqrt(T) =
% sqrt(U0^2*T/2)/sqrt(T) = 325 V / sqrt(2)
% numerisch:

eff = sqrt(1/length(s)*sum(s.^2));
disp('Der Effektivwert für eine Sinusschwingung von 325 V beläuft sich numerisch zu: ')
disp(eff)

%% Aufgabe 2










%% Aufgabe 3



% (a)
t0 = 0.2;
t = linspace(0,1,100000);
s = sin(t);
stan = sinTangent(t, t0);
figure;
plot(t,s,t,stan)
title("sin(x) und die Tangente in t0 über t")

% Momentangeschwindigkeit in t0 = 0.2 von s(t)=sin(t) 
% v(t0) = cos(t0) = s'(t0)

vt0 = cos(t0);

% (b)

dt = [10^-1, 10^-2, 10^-3, 10^-4, 10^-5, 10^-6, 10^-7, 10^-8, 10^-9];

e = abs(vd(t0,@sine,dt)-vt0)/abs(vt0);


% (c)
figure;
loglog(dt,e)
title("Fehler über dt")

disp("Ab einem dt von 10^-9 steigt der Fehler wieder")

% (d)


e = abs(vz(t0,@sine,dt)-vt0)/abs(vt0);

figure;
loglog(dt,e)
title("Fehler über dt mit 2dt approx")

disp("auch hier steigt der Fehler wieder, er sinkt aber zuerst schneller")


%% Aufgabe 4



% (a)
eps = logspace(-8,-2,7);


diffDoub = cos(double(eps))-1;

% (b)



diffSing = (cos(single(eps))-1);


% (c)


err = abs(diffDoub-diffSing)./diffDoub;

figure;
loglog(eps,err);
title("Relativer Fehler über eps")

% (d)

disp("Ein Fehler taucht erst bei ausreichend großen epsilon auf. " + ...
    "Sonst ist Epsilon so klein, dass es keinen Unterschied mehr macht.")




%% Funktionen

function s = sinTangent(x,t0)
    s = sin(t0)+cos(t0)*(x-t0);
end

function res = vd(t0, func, dt)
    res = (func(t0+dt)-func(t0))./(dt);
end

function res = vz(t0, func, dt)
    res = (func(t0+dt)-func(t0-dt))./(2*dt);
end

function res = sine(x)
    res = sin(x);
end

