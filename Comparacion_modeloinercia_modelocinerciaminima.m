
%%modelo sin inercia mínima
% Voy a usar inertancias muy pequeñas. El modelo este debería dar similar al que no tiene inertancia
% L=10^(-6) mmHg s^2/mL 
% [mmHg s^2/mL] = 133.332Pa s^2/(10^(-6) m^3)=133.332(N/(m^2)) s^2/(10^(-6) m^3)=1.3332x10^8 Ns^2/m^5

[P1,V1,Qin1,Qout1,t1]=Simulacion_cardiac_chamber3; %simu con cond inicial 25mL

%% modelo sin inercia
Vinicial=25*10^(-6);
Tspan=60/HR; %duracion de un latido[s]
Tf=5*Tspan; %simulo 5 latidos
Pin=2*133.322; %2mmHg pasado a Pa ; pressure source: right atrium 
Pout=16.5*133.322; %16.5mmHg pasado a Pa ; pressure sink that represents the pressure in the pulmonary artery.

[P2,V2,Qin2,Qout2,t2]=Simulacion_cardiac_chamber_v3(Vinicial,Tf);

subplot(211),plot(t1,V1*10^6,'r',t2,V2*10^6,'b .'),xlabel('tiempo [s]'),ylabel('V(t) [mL]');
subplot(212),plot(t1,P1/133.322,'r',t2,P2/133.322,'b .',t1,repmat(Pin/133.322,size(t1)),t1,repmat(Pout/133.322,size(t1))),xlabel('tiempo [s]'),ylabel('P(t) [mmHg]');