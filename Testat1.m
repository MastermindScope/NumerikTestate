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


s2e = 1.41421356237309504880168872420969807856967187537694807317667;
s2m = 0.41421356237309504880168872420969807856967187537694807317667;

% sqrt2 in Binär Stelle vor dem Komma: 1
% sqrt2 in Binär nach dem Komma:
l = s2m;


m=l;
% implementierung meines Algorithmus für die Berechnung der Bits für der
% Mantisse
for i = 1:50
    mprime = 2*m;
    mprimeprime = 2*(m-0.5);
    if mprime < 1
        dig(i)=0;
        m = mprime;
    elseif mprimeprime < 1
        m = mprimeprime;
        dig(i)=1;
    end
end

% IEEE single point precision ist 32 bits 
% (-1)^b31 * 2^(b30...b23-127)*(1.b23....b0)
% bit 31 is sign, 0 in our case
% bit 30 to 23 are the exponent, (digits before comma) = e - 127 in binary ->
% 127 -> 01111111
% bit 22 to 0 are digits after the comma -> 1 from before the comma and
% then 22 bits from dig
s2b = 1;
exponent = 0;
expBits = [1,1,1,1,1,1,1,0];
for i=0:7
    exponent = exponent+2^i*expBits(i+1);
end



for i = 1:23
    s2b = s2b+dig(i)*2^(-i);
end

s2b = (-1)^0*2^(exponent-127)*s2b;

s2 = single(sqrt(2));

err2 = abs(s2e-s2b);
disp("Absoluter Fehler: ")
disp(err2)

err2 = abs(s2e-s2b)/s2e;
disp("Relativer Fehler: ")
disp(err2)






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
eps = logspace(-10,-2,9);


diffDoub = cos(double(eps))-1;

% (b)



diffSing = (cos(single(eps))-1);


% (c)


err = abs((diffDoub-double(diffSing))./diffDoub);
errabs = abs(diffDoub-double(diffSing));

figure;
loglog(eps,err);
title("Relativer Fehler über eps")

figure;
loglog(eps,errabs);
title("Absoluter Fehler über eps")

% (d)

disp("Der Fehler steigt mit sinkendem Epsilon, " + ...
    "da durch die schlechte Konditionierung der Subtraktion" + ...
    "die Konditionszahl so stark ansteigt, dass der Relative Fehler in die Höhe schießt")




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

