%% NEWTON    sucht minimum in Funktion
% n = newton(g, x0, eps)
% takes symbolic expression g, column vector x0 for starting value 
% (same amount of starting values as g has variables),
% eps termination criterion when norm of two following iteration steps is
% smaller than this tolerance eps > 0
function n = newton(g, x0, eps)


    % variablen aus g
    X = symvar(g);
    % gradient aus g
    grad = gradient(g,X);
    % gradient zwischenspeichern
    grad_ = grad;
    % jacobi-matrix von
    jacobi = jacobian(grad,X);
    % jacobi-matrix zwischenspeichern
    jacobi_ = jacobi;
    % initialisieren des vorigen x-vektors
    xlast = [inf; inf];
    % erstellen der deltax variablen
    D = sym('del',[1 2]);
    % sketchy Lösung um den Variablen-Vektor zu transponieren weil ' nicht
    % geht
    D_ = diag(D);
    D = diag(D_);
    % ersten x-Vektor (Iterationsschritt 0) berechnen, sehr sketchy weil
    % double(solution(...)) nicht ging
    x = x0+vpa(table2array(struct2table(solve(subs(jacobi,X,x0')*D==-subs(grad,X,x0'),D))))';
    % jacobi-matrix und gradient wieder auf symbolische Werte setzen weil
    % subs die expressions in "zahlen" wandelt
    jacobi = jacobi_;
    grad = grad_;

    while norm(xlast - x) > eps
        % letzten iterationsschritt speichern
        xlast = x;
        % s.o.
        x = x+vpa(table2array(struct2table(solve(subs(jacobi,X,x')*D==-subs(grad,X,x'),D))))';
        jacobi = jacobi_;
        grad = grad_;
    end
    
    % Lösung ausgeben
    n = x;



end