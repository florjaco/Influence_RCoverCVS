% The events function in the ODE solvers in Matlab was used to detect when
% a valve should be opened or closed and break out of the sinll11ation. The state
% vector of the systmn is then changed to one of the alternatives in Figure 4.1 and
% the ODE solver is restarted. For the other variables in the state vector that are
% not affected, the final value before the solver was stopped is used as the initial
% condition for the resumed silnulation. This method of simulation nleans that
% instead of sinllllating one nlodel for a predefined tinle period, a series of models
% are sinlulated in succession each heartbeat.

function [P,V,Qentrada,Qsalida,T]=Simulacion_cardiac_chamber3
%le paso Vinicial en m^3, la condicion inicial para la ODE
%tambien le paso Tf, es decir cuanto tiempo quiero simular
%% en este caso estoy simulando un modelo de single chamber (ventrículo derecho) con inercia muy pequeña (uso la inercia que da el autor en el paper Minimal Haemodynamic system including ventricular interaction and valve dynamics en la sección 3.Results pag 136), uso las ctes que usa el autor en su tesis doctoral (tabla 5.1)
%% ver pag 74 y 76-77 ch 5 de la tesis que muestra parámetros y resultados con esos parámetros
global HR V0 lamda P0 R1 R2 Pin Pout Ees Vd L1 L2 tol

HR=80; %beats per minute
Fs=3*HR;
Ts=1/Fs;
Tspan=60/HR; %tiempo a simular: duracion de un latido[s]
Vinicial=25*10^(-6); %[m^3]=10mL %con cond inic por debajo de 23 mL no anda
% Qinicial=250*10^(-6); %[mL/s]
Ti=0;
tol=10^(-10); %tolerancia para las restas. si la resta es menor que esto que de cero.

% t=0:Ts:2; %2seg
V0=0*10^(-6); %[10^(-6)m^3]=[mlitro] %volume at zero pressure 

Ees=54e6; %elastancia de fin de sistole [10^6 N/m^5]=[kPa/litro]==>54 kPa/litro=54e6 N/m^5
Vd=0*10^(-6); %[10^(-6)m^3]=[mlitro] %unstressed chamber volume en l ;unstressed volume is the volume in a chamber that does not
% contribute to an increase in pressure, or the relaxed volume of a chalnber
lamda=23000; %[m^(-3)]%cte de la curva de fin de diastole 
P0=10; %[N/m^2]=[Pa]%condicion inicial de presion 

%en este caso que la camara que simulo es el ventrículo derecho:

R1=1*10^6;%resistencia de la válvula tricúspide [Ns/m^5]
R2=1*10^6;%resistencia de la válvula de la arteria pulmonar [Ns^2/m^5]

% Voy a usar inertancias muy pequeñas. El modelo este debería dar similar al que no tiene inertancia
% L=10^(-6) mmHg s^2/mL 
% [mmHg s^2/mL] = 133.332Pa s^2/(10^(-6) m^3)=133.332(N/(m^2)) s^2/(10^(-6) m^3)=1.3332x10^8 Ns^2/m^5

L1=10^(-6)*1.3332*10^8; %inertancia de la válvula tricúspide [Ns^2/m^5] 
L2=10^(-6)*1.3332*10^8; %inertancia de la válvula de la arteria pulmonar [Ns^2/m^5]

Pin=2*133.322; %2mmHg pasado a Pa ; pressure source: right atrium 
Pout=16.5*133.322; %16.5mmHg pasado a Pa ; pressure sink that represents the pressure in the pulmonary artery.


% Qin=(Pin-Pout)/R;
% Qout=0;
% V(1)=V0;

% [t,x]=ode45(@sistema,0:Ts:7,V0);
% Vinicial=10*10^(-6); %[m^3]=10ml

options1=odeset('Events',@event1);
options2=odeset('Events',@event2);
options3=odeset('Events',@event3);
options4=odeset('Events',@event4);
T=[];
V=[];
Qentrada=[];
Qsalida=[];

for k=1:5 %cuantos latidos quiero simular

%% llenado
Tf=Ti+Tspan;
[t1,x1,te,xe,~]=ode15s(@sistema1,[Ti Tf],[Vinicial 0],options1); %uso ode15s en vez de ode45 xq la señal tiene menos ruidito. Uso ode15s xq en la tesis del autor usa esa herramienta
T=[T;t1]; %vector temporal
V=[V;x1(:,1)]; %volumen en m^3
Qentrada=[Qentrada;x1(:,2)]; %flujo de entrada
Qsalida=[Qsalida;zeros(size(t1))]; %flujo de salida
elastancia=e(T);
P=elastancia.*Ees.*(V-Vd)+(1-elastancia).*P0.*(exp(lamda.*(V-V0))-1);

%voy ploteando las variables y veo como evolucionan
close all
figure();
grid on;
subplot(311),plot(T,P/133.322,T,repmat(Pin/133.322,size(T)),T,repmat(Pout/133.322,size(T))),xlabel('tiempo'), ylabel('P(t) [mmHg]'),legend('P2(t)','Pin','Pout');
subplot(312),plot(T,V*10^6),xlabel('tiempo'), ylabel('V(t) [mL]');
subplot(313),plot(T,Qentrada*10^6,T,Qsalida*10^6,'r'),xlabel('tiempo'), ylabel('flujo[mL/s]'),legend('Qin','Qout');

%diagrama P-V
figure();
plot(V*10^6,P/133.322),ylabel('P(t) [mmHg]'),xlabel('V(t) [mL]'),title('Diagrama P-V');
%% contraccion isovolumetrica
Tf=te+Tspan;
% Tf=t1(end-1)+Tspan;
% [t2,x2,te,xe,~]=ode15s(@sistema3,t1(end-1):Ts:Tf,x1(end-1,1),options3); %uso ode15s en vez de ode45 xq la señal tiene menos ruidito. Uso ode15s xq en la tesis del autor usa esa herramienta

[t2,x2,te,xe,~]=ode15s(@sistema3,[te Tf],xe(1),options3); %uso ode15s en vez de ode45 xq la señal tiene menos ruidito. Uso ode15s xq en la tesis del autor usa esa herramienta
T=[T;t2]; %voy armandome un vector con todos los tiempos
V=[V;x2(:,1)];
Qentrada=[Qentrada; zeros(size(t2))];
Qsalida=[Qsalida;zeros(size(t2))];
elastancia=e(T);
P=elastancia.*Ees.*(V-Vd)+(1-elastancia).*P0.*(exp(lamda.*(V-V0))-1);

%voy ploteando las variables y veo como evolucionan
close all
figure();
grid on
subplot(311),plot(T,P/133.322,T,repmat(Pin/133.322,size(T)),T,repmat(Pout/133.322,size(T))),xlabel('tiempo'), ylabel('P(t) [mmHg]'),legend('P2(t)','Pin','Pout');
subplot(312),plot(T,V*10^6),xlabel('tiempo'), ylabel('V(t) [mL]');
subplot(313),plot(T,Qentrada*10^6,T,Qsalida*10^6,'r'),xlabel('tiempo'), ylabel('flujo[mL/s]'),legend('Qin','Qout');
%diagrama P-V
figure();
plot(V*10^6,P/133.322),ylabel('P(t) [mmHg]'),xlabel('V(t) [mL]'),title('Diagrama P-V');
%% eyeccion
Tf=te+Tspan;
% Tf=t2(end-1)+Tspan;
% [t3,x3,te,xe,~]=ode15s(@sistema2,t2(end-1):Ts:Tf,[x2(end-1) 0],options2); %uso ode15s en vez de ode45 xq la señal tiene menos ruidito. Uso ode15s xq en la tesis del autor usa esa herramienta
[t3,x3,te,xe,~]=ode15s(@sistema2,[te Tf],[xe 0],options2); %uso ode15s en vez de ode45 xq la señal tiene menos ruidito. Uso ode15s xq en la tesis del autor usa esa herramienta
T=[T;t3]; %voy armandome un vector con todos los tiempos
V=[V;x3(:,1)];
Qentrada=[Qentrada; zeros(size(t3))];
Qsalida=[Qsalida;x3(:,2)];
elastancia=e(T);
P=elastancia.*Ees.*(V-Vd)+(1-elastancia).*P0.*(exp(lamda.*(V-V0))-1);

%voy ploteando las variables y veo como evolucionan
close all
figure();
grid on
subplot(311),plot(T,P/133.322,T,repmat(Pin/133.322,size(T)),T,repmat(Pout/133.322,size(T))),xlabel('tiempo'), ylabel('P(t) [mmHg]'),legend('P2(t)','Pin','Pout');
subplot(312),plot(T,V*10^6),xlabel('tiempo'), ylabel('V(t) [mL]');
subplot(313),plot(T,Qentrada*10^6,T,Qsalida*10^6,'r'),xlabel('tiempo'), ylabel('flujo[mL/s]'),legend('Qin','Qout');
%diagrama P-V
figure();
plot(V*10^6,P/133.322),ylabel('P(t) [mmHg]'),xlabel('V(t) [mL]'),title('Diagrama P-V');
%% expansion isovolumetrica
Tf=te+Tspan;
% Tf=t3(end-1)+Tspan;
% [t4,x4,te,xe,~]=ode15s(@sistema3,t3(end-1):Ts:Tf,x3(end-1,1),options4); %uso ode15s en vez de ode45 xq la señal tiene menos ruidito. Uso ode15s xq en la tesis del autor usa esa herramienta
[t4,x4,te,xe,~]=ode15s(@sistema3,[te Tf],xe(1),options4); %uso ode15s en vez de ode45 xq la señal tiene menos ruidito. Uso ode15s xq en la tesis del autor usa esa herramienta
%aca como que no corta cuando la P1>P2 para poder abrir la valvula
%tricuspide, por lo tanto, el volumen me queda negativo!!
T=[T;t4]; %voy armandome un vector con todos los tiempos
V=[V;x4(:,1)];
Qentrada=[Qentrada; zeros(size(t4))];
Qsalida=[Qsalida;zeros(size(t4))];
elastancia=e(T);
P=elastancia.*Ees.*(V-Vd)+(1-elastancia).*P0.*(exp(lamda.*(V-V0))-1);

% voy ploteando las variables y veo como evolucionan
close all
figure();
grid on
subplot(311),plot(T,P/133.322,T,repmat(Pin/133.322,size(T)),T,repmat(Pout/133.322,size(T))),xlabel('tiempo'), ylabel('P(t) [mmHg]'),legend('P2(t)','Pin','Pout');
hold on
plot(T,elastancia*30,'m');
hold off
grid on
subplot(312),plot(T,V*10^6),xlabel('tiempo'), ylabel('V(t) [mL]');
grid on
subplot(313),plot(T,Qentrada*10^6,T,Qsalida*10^6,'r'),xlabel('tiempo'), ylabel('flujo[mL/s]'),legend('Qin','Qout');
grid on
%diagrama P-V
figure();
plot(V*10^6,P/133.322),ylabel('P(t) [mmHg]'),xlabel('V(t) [mL]'),title('Diagrama P-V');

Vinicial=xe;
Ti=te;
% Vinicial=x4(end-1);
% Ti=t4(end-1);
% Qinicial=0;

end
% end

elastancia=e(T);
P=elastancia.*Ees.*(V-Vd)+(1-elastancia).*P0.*(exp(lamda.*(V-V0))-1);

%voy ploteando las variables y veo como evolucionan
close all
figure();
grid on
subplot(311),plotyy(T,P/133.322,T,e(T)),xlabel('tiempo'), ylabel('P(t) [mmHg]');
hold on;
plot(T,repmat(Pin/133.322,size(T)),'r',T,repmat(Pout/133.322,size(T)),'k'),legend('P2(t)','Pin','Pout','e(t)');
subplot(312),plot(T,V*10^6),xlabel('tiempo'), ylabel('V(t) [mL]');
subplot(313),plot(T,Qentrada*10^6,T,Qsalida*10^6,'r'),xlabel('tiempo'), ylabel('flujo[mL/s]'),legend('Qin','Qout');
%diagrama P-V
figure();
plot(V*10^6,P/133.322),ylabel('P(t) [mmHg]'),xlabel('V(t) [mL]'),title('Diagrama P-V');

% T=[t1;t2;t3;t4];
% V=[x1(1,:);x2(1,:);x3(1,:);x4(1,:)];

% V=x;
% 
% elastancia=e(t);
% P=elastancia.*Ees.*(x-Vd)+(1-elastancia).*P0.*exp(lamda.*(x-V0)-1);
% 

% figure();
% plot(t,elastancia),xlabel('tiempo'),ylabel('driver function e(t)');
% figure();
% plot(t,P/133.322),xlabel('tiempo'), ylabel('P(t) [mmHg]'); %esta en Pa y lo paso a mmHg
% figure();
% plot(t,V*10^6),xlabel('tiempo'), ylabel('V(t) [mL]'); %esta en m^3 lo paso a mL
% figure();
% plot(t,Qin*10^6),xlabel('tiempo'), ylabel('flujo de entrada [mL/s]'); %esta en m^3/s y lo paso a mL/s
% figure();
% plot(t,Qout*10^6),xlabel('tiempo'), ylabel('flujo de salida [mL/s]'); %esta en m^3/s y lo paso a mL/s

function et=e(t)%funcion elastancia variante en tiempo
    %tengo HR como variable global
          taux=mod(t,60/HR); %para que se repita con los latidos
%         et=exp(-80.*(taux-60/HR/2).^2); %centro la gaussiana en el latido
          et=exp(-80.*(taux-0.27).^2); 
end

%% llenado ventricular
function dxdt=sistema1(t,x) %la variable de estado x es un vector columna con el volumen,el flujo de entrada
  %este sería el caso de llenado de la camara (valvula tricuspide abierta y pulmonar cerrada)
%x=[V, Qin]
V_p=x(2); %derivada del volumen=Qin-Qout (como en este caso Qout=0==>V_p=Qin)
P2=e(t)*Ees*(x(1)-Vd)+(1-e(t))*P0*(exp(lamda*(x(1)-V0))-1); %presion en la camara durante el llenado
deltaP=Pin-P2;
% if deltaP<=tol
%     deltaP=0;
% end
Qin_p=(deltaP-x(2)*R1)/L1; %considero la inertancia, Qin_p es la derivada del flujo de entrada Qin
dxdt=[V_p;Qin_p]; %la variacion neta de volumen (Vpunto se debe al flujo de entrada solamente, xq no hay flujo de salida)
end

function [cierre stop direction]=event1(t,x) %condicion1 (en que momento se cierra la valvula tricuspide)
    cierre= x(2); %busca en que momento Qin (la segunda componente de la var de estado) se anula
%     P2=e(t)*Ees*(x(1)-Vd)+(1-e(t))*P0*exp(lamda*(x(1)-V0)-1); %presion en la camara durante el llenado
%     cierre= P2-Pin; %busca en que momento el gradiente de presion se vuelve positivo, es decir, P2>P1 (en que momento se cierra la valvula)
    stop= 1; %corta cuando sucede el evento, o sea cuando Qin se anula
    direction=-1; %agarro solo si se anula pero tiene derivada negativa, o sea pasa por cero para disminuir (no xq se aumenta)
%     direction=1; %agarro solo si se anula pero tiene derivada positiva, o sea pasa por cero para aumentar
end

%% eyeccion 
function dxdt=sistema2(t,x) %la variable de estado x es un vector columna con el volumen,el flujo de entrada y el flujo de salida
 %este seria el caso de vaciado de la camara (valvula pulmonar abierta y tricuspide cerrada)
    %x=[V, Qout]
V_p=-x(2); %la variacion del volumen se debe al flujo de salida!
P2=e(t)*Ees*(x(1)-Vd)+(1-e(t))*P0*(exp(lamda*(x(1)-V0))-1);
deltaP=P2-Pout;
% if deltaP<=tol
%     deltaP=0;
% end
Qout_p=(deltaP-R2*x(2))/L2; %considero la inertancia, Qout_p es la derivada del flujo de entrada Qin

dxdt=[V_p;Qout_p];%la variacion neta de volumen (Vpunto se debe al flujo de salida solamente, xq no hay flujo de entrada)
end

function [cierre stop direction]=event2(t,x) %condicion2 (en que momento se cierra la valvula pulmonar)
    cierre= x(2); %busca en que momento Qout (la segunda componente de la var de estado) se anula
%     cierre= Pout-e(t)*Ees*(x(1)-Vd)+(1-e(t))*P0*exp(lamda*(x(1)-V0)-1); %busca en que momento el gradiente de presion se vuelve positivo, es decir, P2>P1 (en que momento se cierra la valvula)
    stop= 1; %corta cuando sucede el evento, o sea cuando Qout se anula
    direction=-1; %agarro solo si se anula pero tiene derivada negativa, o sea pasa por cero para disminuir (no xq se aumenta)
end

%% contraccion o expansion isovolumetrica
function dxdt=sistema3(t,x) %la variable de estado x es el volumen 
    %corresponde a contraccion o expansion isovolumetrica en la que no hay
    %flujo de entrada o salida de la camara==> No hay variación del Volumen
    %en la camara
dxdt=0;

end

function [cierre stop direction]=event3(t,x) %condicion3 (en que momento se abre la valvula pulmonar)
    P2=e(t)*Ees*(x-Vd)+(1-e(t))*P0*(exp(lamda*(x-V0))-1);
    cierre= P2-Pout; %busca en que momento el gradiente de presion se vuelve negativo, es decir, P3<P2 (en que momento se abre la valvula)
    stop= 1; %corta cuando sucede el evento, o sea cuando P2-P3 se anula
    direction=1; %agarro solo si se anula pero tiene derivada positiva, o sea pasa por cero porque va a aumentar 
end

function [cierre stop direction]=event4(t,x) %condicion4 (en que momento se abre la valvula tricuspide)
    P2=e(t)*Ees*(x-Vd)+(1-e(t))*P0*(exp(lamda*(x-V0))-1);
    cierre= Pin-P2; %busca en que momento el gradiente de presion se vuelve negativo es decir, P2<P1 (en que momento se abre la valv)
    stop= 1; %corta cuando sucede el evento, o sea cuando P1-P2 se anula
    direction=1;%1; %agarro solo si se anula pero tiene derivada negativa, o sea pasa por cero porque va a disminuir 
end

end