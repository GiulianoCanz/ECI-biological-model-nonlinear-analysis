% ECI_model.m = script principale per l'analisi del modello;
% ECI_equations.m = equazioni del modello;
% hill_function.m = una funzione presente nelle equazioni del modello;
% p_hill_function.m = versione pointwise della precedente;
% parameters.mat = parametri da utilizzare;
% intersections.m = funzione usata per trovare gli stati di equilibrio.
% stability.m = funzione usata nel diagramma di biforcazione per stabilire
% la stabilità del punto

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Abbiamo 3 equazioni differenziali rappresentanti la dinamica
%delle densità di 3 tipologie di cellule: cancerose(C),dendritiche(D) e killer(K)

%i parametri delle equazioni sono contenuti in lambda,comp,n,g,k.
%rho è un ulteriore parametro, dal quale dipendono 2 valori di lambda.
%rho modula la risposta immunitaria dell'organismo.

clear all , close all
load('parameters.mat')

%settiamo la immune recognition (rho): 0.2(bassa), 1(nella norma), 2(migliorata)
rho = 2;
lambda.KCp = 1+2*rho;  %perchè lambdaKC0 è sempre 3
%settiamo i valori delle interazioni specifiche
%scegliere il caso da analizzare , e commentare gli altri


%CASO 1
lambdaDC0 = 2;
lambda.DCp = 1+rho*(lambdaDC0 - 1);
lambda.CD1m = 0.1;


%{
%CASO 2
lambdaDC0 = 1;
lambda.DCp = 1+rho*(lambdaDC0 - 1);
lambda.CD1m = 1;
%}

%{
%CASO 3
lambdaDC0 = 2;
lambda.DCp = 1+rho*(lambdaDC0 - 1);
lambda.CD1m = 1;
%}

tic

tspan = [0 1000];
ini_conds = [1000 1000 1000]; %partiamo dalla densità di un migliaio di cellule/microlitro
[t,x] = ode45(@(t,x) ECI_equations(t,x,g,k,comp,lambda,n),tspan,ini_conds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% SPAZIO DELLE FASI %%%%%%%%%%%%%%%%%%%%%%%%%%
figure, plot3(x(:,2),x(:,1),x(:,3)),title('phase space'),xlabel('dendritic cells (D)'),ylabel('cancer cells (C)'),zlabel('killer cells (K)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TRAIETTORIE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure,plot(t,x(:,1)),title('cancer cells (C)'),xlabel('time'),ylabel('C');
figure,plot(t,x(:,2)),title('dendritic cells (D)'),xlabel('time'),ylabel('D');
figure,plot(t,x(:,3)),title('killer cells (K)'),xlabel('time'),ylabel('K');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMPO VETTORIALE %%%%%%%%%%%%%%%%%%%%%%%%%%%
syms K 
s = 50;
C = [0:s:1800];
D = [0:s:1800];
[D,C] = meshgrid(D,C);

dxdt3 = g.k0.*D.*p_hill_function(C,comp.CK10,n.ck1,lambda.CK1p).*p_hill_function(C,comp.CK20,n.ck2,lambda.CK2m) - k.k.*K.*p_hill_function(C,comp.CK30,n.ck3,lambda.CK3p);

for i=1:size(dxdt3,2)
    for j=1:size(dxdt3)
       K1(i,j)=solve(dxdt3(i,j)==0,K); %mi trovo la K dalla derivata dK/dt=0 per ogni coppia (D,C)
    end
end

K=K1; %i valori di K ottenuti li inserisco nella variabile K per sostituirli nelle altre due derivate

dxdt1 = g.c0.*p_hill_function(C,comp.CC0,n.cc,lambda.CCp) - k.c.*C.*p_hill_function(D,comp.DC0,n.dc,lambda.DCp).*p_hill_function(K,comp.KC0,n.kc,lambda.KCp);
dxdt2 = g.d0.*p_hill_function(D,comp.DD0,n.dd,lambda.DDp).*p_hill_function(K,comp.KD0,n.kd,lambda.KDp).*p_hill_function(C,comp.CD10,n.cd1,lambda.CD1m).*p_hill_function(C,comp.CD20,n.cd2,lambda.CD2p)-k.d.*D;
figure
q = quiver(D,C,dxdt2,dxdt1,5,'k'),title('vector field');
annotation('textbox',[.65 .75 .05 .15],'BackgroundColor','w','String','dD = 0 & dK = 0','FitBoxToText','on','Color','blue');
annotation('textbox',[.65 .75 .005 .01],'BackgroundColor','w','String','dC = 0 & dK = 0','FitBoxToText','on','Color','red');
xlabel('dendritic cells'),ylabel('cancer cells');

hold on

%%%%%%%%%%%%%%%%%%%%%%%%% NULLCLINES %%%%%%%%%%%%%%%%%%%%%%%%%%%
syms C D K

dxdt1 = g.c0.*p_hill_function(C,comp.CC0,n.cc,lambda.CCp) - k.c.*C.*p_hill_function(D,comp.DC0,n.dc,lambda.DCp).*p_hill_function(K,comp.KC0,n.kc,lambda.KCp);
dxdt2 = g.d0.*p_hill_function(D,comp.DD0,n.dd,lambda.DDp).*p_hill_function(K,comp.KD0,n.kd,lambda.KDp).*p_hill_function(C,comp.CD10,n.cd1,lambda.CD1m).*p_hill_function(C,comp.CD20,n.cd2,lambda.CD2p)-k.d.*D;
dxdt3 = g.k0.*D.*p_hill_function(C,comp.CK10,n.ck1,lambda.CK1p).*p_hill_function(C,comp.CK20,n.ck2,lambda.CK2m) - k.k.*K.*p_hill_function(C,comp.CK30,n.ck3,lambda.CK3p);

%PRIMA NULLCLINE --->  dD = 0 & dK = 0
Sk = solve(dxdt3==0,K); % valori di K per cui dK = 0
eq = subs(dxdt2,'K',Sk); %sostituisco in dxdt2, che quindi conterrà i valori di K per cui dK = 0
fun1 = matlabFunction(eq); %converto il simbolico in funzione
C = [0:1800];
D = [0:1800];
[D,C] = meshgrid(D,C);
fun1 = fun1(C,D);
c1 = contour(fun1,[0 0],'b','LineWidth',3); %plotto la curva per cui fun = 0
hold on
contourTable1 = getContourLineCoordinates(c1);

%SECONDA NULLCLINE ---> dC = 0 & dK = 0
eq=subs(dxdt1,'K',Sk);  % sostituisco in dxdt1
fun2 = matlabFunction(eq); % converto il simbolico in equazione
fun2 = fun2(C,D);          
c2  = contour(fun2,[0 0],'r','LineWidth',3); %plotto la curva per la quale fun = 0
hold on
contourTable2 = getContourLineCoordinates(c2);

%%%%%%%%%%%%%%%%%%%%%%% FIXED POINTS %%%%%%%%%%%%%%%%%%%%%%%%% 
 
%le intersezioni tra le nullclines sono gli stati di equilibrio del sistema 
[xf,yf] = intersections(contourTable1.X,contourTable1.Y,contourTable2.X,contourTable2.Y,0);  
 
for i=1:length(xf)
plot(xf(i),yf(i),'x','LineWidth',3); 
hold on  
end
 
%%%%%%%%%%%%%%%%%%%%%%% STABILITY OF THE FIXED POINTS %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
% studiamo la stabilità di questi stati 
%per ognuno di essi studiamo gli autovalori della jacobiana 
% per rho = 2 : primo punto e terzo stabili, secondo punto instabile
syms C D K 
Sk = solve(dxdt3==0,K); %valori di K per cui dK/dt = 0
J = jacobian([dxdt1; dxdt2; dxdt3],[C D K]); % matrice jacobiana  
eigvi = zeros(3,length(xf));

for i = 1:length(xf)
K = subs(Sk,{'C','D'},{yf(i),xf(i)});
Ji = subs(J,{'C','D','K'},{yf(i), xf(i),K}); 
eigvi(:,i) = double(eig(Ji));
end

%%%%%%%%%%%%%%%%%%%%%%% BIFURCATION DIAGRAMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

syms rho C D K
%scriviamo dxdt1 in funzione di rho (solo questa equazione dipende da rho)
dxdt1 = g.c0.*p_hill_function(C,comp.CC0,n.cc,lambda.CCp) - k.c.*C.*p_hill_function(D,comp.DC0,n.dc,1+rho*(lambdaDC0-1)).*p_hill_function(K,comp.KC0,n.kc,1+rho*2);
%dxdt2 = g.d0.*p_hill_function(D,comp.DD0,n.dd,lambda.DDp).*p_hill_function(K,comp.KD0,n.kd,lambda.KDp).*p_hill_function(C,comp.CD10,n.cd1,lambda.CD1m).*p_hill_function(C,comp.CD20,n.cd2,lambda.CD2p)-k.d.*D;
%dxdt3 = g.k0.*D.*p_hill_function(C,comp.CK10,n.ck1,lambda.CK1p).*p_hill_function(C,comp.CK20,n.ck2,lambda.CK2m) - k.k.*K.*p_hill_function(C,comp.CK30,n.ck3,lambda.CK3p);
J = jacobian([dxdt1; dxdt2; dxdt3],[C D K]);
figure
for r=0:0.05:3 %facciamo variare rho ad ogni ciclo
%ad ogni ciclo calcolo i punti di intersezione delle nullclines e li plotto
%col giusto colore in base alla stabilità
c_eq = subs(dxdt1,'rho',r); %aggiorno dxdt1 col nuovo valore di rho
Sk = solve(dxdt3==0,K); % valori di K per cui dK = 0
eq = subs(dxdt2,'K',Sk); %sostituisco in dxdt2, che quindi conterrà solo C e D
fun1 = matlabFunction(eq); %converto il simbolico in funzione
C = [0:1800];
D = [0:1800];
[D,C] = meshgrid(D,C);
fun1 = fun1(C,D);
c1 = contourc(fun1,[0 0]); %prima nullcline --> dxdt2 = 0 & dxdt3 = 0
contourTable1 = getContourLineCoordinates(c1);
eq=subs(c_eq,'K',Sk);  % sostituisco anche in dxdt1
fun2 = matlabFunction(eq);
if r==0
    fun2=fun2(C);  
else
    fun2 = fun2(C,D);
end
c2  = contourc(fun2,[0 0]);%seconda nullcline --> dxdt1 = 0 & dxdt3 = 0
contourTable2 = getContourLineCoordinates(c2);
[xf,yf] = intersections(contourTable1.X,contourTable1.Y,contourTable2.X,contourTable2.Y,0); %calcoliamo le C e le D dei punti di equilibrio
for i = 1:size(yf,1)
    Sk1 = solve(dxdt3==0,K); %valori di K per cui dK/dt = 0
    Ki = subs(Sk1,{'C','D'},{yf(i),xf(i)}); %sostituiamo le coordinate (C e D) dei punti fissi, in questo modo troviamo il valore di K corrispondente
    Jr = subs(J,{'C','D','K','rho'},{yf(i), xf(i),Ki,r}); %per ogni punto calcoliamo la jacobiana
    eigv = double(eig(Jr)) %estraiamo gli autovalori
    %pause
    stab = stability(eigv);
    if stab == 1 %se stabile, plotta blu
        subplot(3,1,1),plot(r,yf(i),'.b'),xlabel('rho'),ylabel('cancer cells');  % C --> yf(i)
        hold on
        subplot(3,1,2),plot(r,xf(i),'.b'),xlabel('rho'),ylabel('dendritic cells');  % D --> xf(i)
        hold on
        subplot(3,1,3),plot(r,Ki,'.b'),xlabel('rho'),ylabel('killer cells');  % K --> Ki
        hold on
    end
    if stab == 0  %se instabile, plotta rosso
        subplot(3,1,1),plot(r,yf(i),'.r'),xlabel('rho'),ylabel('cancer cells');  % instabile --> rosso
        hold on
        subplot(3,1,2),plot(r,xf(i),'.r'),xlabel('rho'),ylabel('dendritic cells');  % D --> xf(i)
        hold on
        subplot(3,1,3),plot(r,Ki,'.r'),xlabel('rho'),ylabel('killer cells');  % K --> Ki
        hold on
    end
 end %end del for 
end %end del for principale
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
