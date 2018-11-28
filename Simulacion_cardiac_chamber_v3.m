function [P,V,Qin,Qout,t]=Simulacion_cardiac_chamber_v3(Vinicial,Tf)
%le paso Vinicial en m^3, la condicion inicial para la ODE
%tambien le paso Tf, es decir cuanto tiempo quiero simular
%% en este caso estoy simulando un modelo de single chamber (ventrículo derecho) sin inercia, uso las ctes que usa el autor en su tesis doctoral (tabla 5.1)
%% anda 10 points!
global HR V0 lamda P0 R1 R2 Pin Pout Ees Vd 

HR=80; %beats per minute
Fs=3*HR;
Ts=1/Fs;

% t=0:Ts:2; %2seg
V0=0*10^(-6); %[10^(-6)m^3]=[mlitro] %volume at zero pressure 

Ees=54e6; %elastancia de fin de sistole [10^6 N/m^5]=[kPa/litro]==>54 kPa/litro=54e6 N/m^5
Vd=0*10^(-6); %[10^(-6)m^3]=[mlitro] %unstressed chamber volume en l ;unstressed volume is the volume in a chamber that does not
% contribute to an increase in pressure, or the relaxed volume of a chalnber
lamda=23000; %[m^(-3)]%cte de la curva de fin de diastole 
P0=10; %[N/m^2]=[Pa]%condicion inicial de presion 

%en este caso que la camara que simulo es el ventrículo derecho:

R1=1*10^6;%resistencia de la válvula tricúspide [Ns/m^5]
R2=1*10^6;%resistencia de la válvula de la arteria pulmonar [Ns/m^5]

Pin=2*133.322; %2mmHg pasado a Pa ; pressure source: right atrium 
Pout=16.5*133.322; %16.5mmHg pasado a Pa ; pressure sink that represents the pressure in the pulmonary artery.


% Qin=(Pin-Pout)/R;
% Qout=0;
% V(1)=V0;

% [t,x]=ode45(@sistema,0:Ts:7,V0);
% Vinicial=10*10^(-6); %[m^3]=10ml
[t,x]=ode15s(@sistema,0:Ts:Tf,Vinicial); %uso ode15s en vez de ode45 xq la señal tiene menos ruidito. Uso ode15s xq en la tesis del autor usa esa herramienta
V=x;

elastancia=e(t);
P=elastancia.*Ees.*(x-Vd)+(1-elastancia).*P0.*exp(lamda.*(x-V0)-1);
% Pin=10*0.1333; %10mmHg pasados a kPa
% Pout=100*0.1333; %100mmHg pasados a kPa
% R1=19;%resistencia de la circ pulmonar[kPa s/l]
% R2=140;%resistencia periferica[kPa s/l]
Qin=(Pin-P)/R1;
Qin(Qin<=0)=0;
% ind= Qin<=0;
% Qin(ind)=0;
Qout=(P-Pout)/R2;
Qout(Qout<=0)=0;
% ind=find (Qout<=0);
% Qout(ind)=0;

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

function dxdt=sistema(t,x) %la variable de estado x es el volumen
P2=e(t)*Ees*(x-Vd)+(1-e(t))*P0*(exp(lamda*(x-V0))-1);
Qin=(Pin-P2)/R1;
if Qin<=0 
    Qin=0;
end
Qout=(P2-Pout)/R2;
if Qout<=0 
    Qout=0;
end
dxdt=Qin-Qout;

end
end