%% Numerik Testat 2 
% Vinzenz Götz
close all

%% Aufgabe 1




%% Aufgabe 2

% (b) + (c)
n = 8:2:20;
t_f = linspace(0, 10, 100);
for i = 1:length(n)
    t = linspace(0,10,n(i)+1)';
    y = rand(n(i)+1,1);
    V = zeros(n(i)+1, n(i)+1);
    for k = 1:n(i)+1
        V(:,n(i)+2-k) = t.^(k-1);
    end
    disp(strcat("Die Konditionszahl ist: ",num2str(cond(V))))
    a = V\y;
    figure;
    plot(t_f,polyval(a,t_f))
    hold
    plot(t,y,"x")
end



% (d)
for i = 1:length(n)
    t = linspace(0,10,n(i)+1)';
    y = sin(t);
    V = zeros(n(i)+1, n(i)+1);
    for k = 1:n(i)+1
        V(:,n(i)+2-k) = t.^(k-1);
    end
    disp(strcat("Die Konditionszahl ist: ",num2str(cond(V))))
    a = V\y;
    figure;
    plot(t_f,polyval(a,t_f))
    hold
    plot(t,y,"x")
end


%% Aufgabe 3
% (b) + (c)
close all

x = 0:1e-3:1;

%analytische Lösung der DGLn
T1 = csc(1)*sin(x)-x;
T2 = -x.^2 -2*cos(x) + (2*cos(1) - 1)*csc(1)*sin(x)+2;


figure(100)
plot(x,T1)

figure(101)
plot(x,T2)

h = [1e-1 1e-2 1e-3 1e-4];

for i = 1:length(h)
    size = 1/h(i)-1;
    iVec = linspace(1,size,size);
    jVec = linspace(1,size,size);
    a = ones(1,size)*(1-2/(h(i)^2));
    b = ones(1,size-1)/(h(i)^2);
    M = sparse(iVec,jVec,a,size,size);

    iVec = linspace(1,size-1,size-1);
    jVec = linspace(2,size,size-1);
    M = M+sparse(iVec,jVec,b,size,size);
    M = M+sparse(jVec,iVec,b,size,size);
    

    xn = h(i):h(i):1-h(i);
    f1 = -xn;
    f2 = -xn.^2;

    T1 = csc(1)*sin(xn)-xn;
    T2 = -xn.^2 -2*cos(xn) + (2*cos(1) - 1)*csc(1)*sin(xn)+2;

    T1n = M^(-1)*f1';   %f1 ohne M funktioniert nicht, deswegen hier mit M^(-1)
                        %ich weiß, dass das weniger effizient ist RAM
                        %frisst, aber ich weiß nicht woher das seltsame
                        %Verhalten kommt
    T2n = M^(-1)*f2';

    figure(100)
    hold on
    plot(xn,T1n,"x")

    figure(101)
    hold on
    plot(xn,T2n,"x")

    kappa(i) = cond(M);

    sprintf("Fehler für %e Elemente mit f1: %d",1/h(i),norm(T1n(1)-T1(1))/norm(T1(1)))
    sprintf("Fehler für %e Elemente mit f2: %d",1/h(i),norm(T2n(1)-T2(1))/norm(T2(1)))
end

sprintf("Der Fehler pro Zehnerpotenz in h sinkt um zwei Zehnerpotenzen")


figure(102)
    loglog(h,kappa)
    xlabel("Schrittweite h")
    ylabel("Konditionszahl")
sprintf("Die Konditionszahl von M wächst mit zunehmmender Feinheit des Netzes exponentiell")