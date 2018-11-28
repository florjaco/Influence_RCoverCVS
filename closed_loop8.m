function [PRV, Vrv, Qtc, Qpv, PLV, Vlv, Qmt, Qav, Vtot, T, Vvc, Vpa, Vpu, Vao, PVC, PPA, PPU,PAO]=closed_loop8(ploteo)
    global HR FR Rpul Ees_vc Ees_pa Ees_pu Ees_ao Rsys Vi_rv Vi_lv Vi_pu Vi_ao Vi_vc Vi_pa Qi_tc Qi_mt Qi_pul Qi_sys V0_spt V0_rvf V0_lvf V0_pcd Vd_spt Vd_lvf Vd_rvf lambda_spt lambda_lvf lambda_rvf lambda_pcd Ees_spt Ees_rvf Ees_lvf P0_spt P0_lvf P0_rvf P0_pcd Rtc Ltc Lpv Rpv Rmt Lmt Rav Lav Pra Ppa Ppu Pla Pao Pvc
%     INPUT: si le mando un uno a ploteo plotea paso a paso. si no plotea al
%     final
    if nargin<1
        ploteo=0;
    end
    FR=12; %frec respiratoria
    HR=80;                                      %beats per minute
    % Fs=3*HR;
    % Ts=1/Fs;
    Tspan=60/HR;                                %tiempo a simular: duracion de un latido[s]
%     Ti=0;

%     Vi_spt=2.4*10^(-6);                         %10^(-6)[m^3]=6.5mL
    Vi_rv= 120*10^(-6);%25*10^(-6);                          %10^(-6)[m^3]=25mL
    Vi_lv= 140*10^(-6);%25*10^(-6);                          %10^(-6)[m^3]=25mL
    
    Qi_tc=0;%50*10^(-6)/(300*10^(-3)); %50mL/300seg: vol de llenado ventr/dur del llenado;[m^3/s]%250*10^(-6);                         %[mL/s]
    Qi_mt=0;%50*10^(-6)/(300*10^(-3));
    Qi_pul=0;
    Qi_sys=0;

    V0_spt=2*10^(-6);                            %[10^(-6)m^3]=[mlitro], volume at zero pressure 
    V0_rvf=0*10^(-6);                            %[10^(-6)m^3]=[mlitro], volume at zero pressure 
    V0_lvf=0*10^(-6);                            %[10^(-6)m^3]=[mlitro], volume at zero pressure 
    V0_pcd=200*10^(-6);                          %[10^(-6)m^3]=[mlitro], volume at zero pressure 

    Vd_spt=2*10^(-6); 
    Vd_lvf=0*10^(-6); 
    Vd_rvf=0*10^(-6); %[10^(-6)m^3]=[mlitro]     %unstressed chamber volume en l, unstressed volume is the volume in a chamber that does not contribute to an increase in pressure, or the relaxed volume of a chalnber

    lambda_spt=435000; %[m^(-3)]%cte de la curva de fin de diastole 
    lambda_lvf=33000;  %[m^(-3)]%cte de la curva de fin de diastole            
    lambda_rvf=23000;  %[m^(-3)]%cte de la curva de fin de diastole                           
    lambda_pcd=30000;  %[m^(-3)]%cte de la curva de fin de diastole

    Ees_spt=6500e6;                              %elastancia de fin de sistole del SPT [10^6 N/m^5]=[kPa/litro]==>6500 kPa/litro=6500e6 N/m^5
    Ees_rvf=2*54e6;%123.12e6;%2*54e6;%54e6;                                %elastancia de fin de sistole del VD [10^6 N/m^5]=[kPa/litro]==>54 kPa/litro=54e6 N/m^5
    Ees_lvf=228e6;%2*100e6;%100e6;                              %elastancia de fin de sistole del VI [10^6 N/m^5]=[kPa/litro]==>100 kPa/litro=100e6 N/m^5
 % en este caso es loop cerrado asi que uso las elastancias de vena cava,
 % arteria pulmonar, vena pulmonar y aorta (se agregan esas 4 camaras al cerrar el loop)
    Ees_vc=1.3e6;
    Ees_pa=72e6;
    Ees_pu=1.9e6;
    Ees_ao=98e6;
 
    P0_spt=148;                                  %[N/m^2]=[Pa]%condicion inicial de presion
    P0_lvf=10;                                   %[N/m^2]=[Pa]%condicion inicial de presion
    P0_rvf=10;                                   %[N/m^2]=[Pa]%condicion inicial de presion
    P0_pcd=66.7;                                 %[N/m^2]=[Pa]%condicion inicial de presion 
%     Pth=-4*133.322;                              %-4mmHg pasado a Pa; esto deberia ser oscilatorio para simular inspiracion y espiracion; presion promedio del torax. en ventilacion respiratoria mecanica Pth esta aumentado: Pth=0 mmHg
%Pth tiene que ser variable
    % A_rvf 
    % A_lvf

    Rtc=1*10^6;                                 %resistencia de la válvula tricúspide [Ns/m^5]
    Ltc=1.3*10^4;                               %inertancia de la válvula tricúspide [Ns/m^5]
    Rpv=1*10^6;                                 %resistencia de la válvula pulmonar [Ns^2/m^5]
    Lpv=2*10^4;                                 %inertancia de la válvula pulmonar [Ns^2/m^5]
    Rmt=1.6*10^6; %6.1*10^6; 
    Lmt=1.3*10^4; 
    Rav=2.75*10^6;
    Lav=5*10^4;
    % en este caso es loop cerrado asi que uso Rpul y Rsys
    Rpul=0.5*9.4*10^6;%9.4*10^6;
    Rsys=0.7*170*10^6;%170*10^6;

    Pra=2*133.322;                              %2mmHg pasado a Pa ; pressure source: right atrium 
    Ppa=16.5*133.322;                           %16.5mmHg pasado a Pa ; pressure sink that represents the pressure in the pulmonary artery.
    Pla=12*133.322;  %Pmedia en Left auricle    %12mmHg pasado a Pa ; pressure source: left atrium 
    Pao=65*133.322;%80*133.322;  %Pmedia en aorta           %70mmHg pasado a Pa ; pressure sink that represents the pressure in the aorta.

    Vi_pu=Pla/(Ees_pu);%Pra/(Ees_pu*e(0)); %lo despejo de la ecuacion de calculo de P de la camara elastica en fn de V: Pes(v)=Ees*e(t)*V
    Vi_ao=Pao/(Ees_ao);%Pao/(Ees_ao*e(0));
    Vi_vc=Pra/(Ees_vc);%Pla/(Ees_vc*e(0));
    Vi_pa=Ppa/(Ees_pa);%Ppa/(Ees_pa*e(0));
    Vtot=Vi_rv+Vi_lv+Vi_pu+Vi_vc+Vi_ao+Vi_pa; %volumen total circulante en el sistema cerrado
    
    Ti=0;
    T=[];
    Vlv=[];
    Vrv=[];
    Qtc=[];
    Qpv=[];
    Qmt=[];
    Qav=[];
    Vpu=[];
    Vvc=[];
    Vao=[];
    Vpa=[];
    
    options1=odeset('Events',@event1);
    options2=odeset('Events',@event2);
    options3=odeset('Events',@event3);
    options4=odeset('Events',@event4);
    options5=odeset('Events',@event5);
    options6=odeset('Events',@event6);
    options7=odeset('Events',@event7);
    options8=odeset('Events',@event8);
    %% FUNCION ELASTANCIA VARIANTE EN EL TIEMPO
    function et=e(t)                            %funcion elastancia variante en tiempo
        taux=mod(t,60/HR);                      %para que se repita con los latidos
        et=exp(-80.*(taux-0.27).^2);   
    end
    %% FUNCION PTH (MODIFICA LA PRESION PLEURAL EN FUNCION DEL TIEMPO->INSPIRACION Y ESPIRACION)
    function Ppl=Pth(t)
        Tresp=60/FR; %duracion del ciclo resp
        taux=mod(t,Tresp);   %para que se repita en cada respiracion
        %duración de la fase inspiratoria: 1/3 del ciclo resp
        Tinsp=Tresp/3;
        %duracion de la fase espitaroria: 2/3 del ciclo resp
        Tesp=Tresp-Tinsp;
        if taux<=Tinsp
            Ppl=2.2*exp(-taux/(Tinsp/5.39))-5.88; 
        elseif taux>Tinsp && taux<=Tresp
        % Pth(i)=2.2*(1-exp(-(taux(i)-Tinsp)/(Tesp/13.8)));
        % Pth(i)=2.2*(1-exp(-(taux(i)-Tinsp)/(Tesp/6.9)));
        Ppl=2.2*(1-exp(-(taux-Tinsp)/(Tesp/4.6)))-5.88;
        end
        %lo calcule en mmHg
        Ppl=Ppl*133.322; %lo paso a Pa
    end

for latido=1:15 %cuantos latidos quiero simular
    %Variable de estado general:
    %x=[Vpu, Qmt, Vlv, Qav, Vao, Vvc, Qtc, Vrv, Qpv, Vpa]
%     Vpu: Volumen en la vena pulmonar (se trata como una cámara elastica)
%     Qmt: Flujo a través de la válvula mitral
%     Vlv: Volumen en el ventriculo izquierdo
%     Qav: Flujo a través de la válvula aórtica
%     Vao: Volumen en la aorta (se trata como una cámara elastica)
%     Vvc: Volumen en la vena cava (se trata como una cámara elastica)
%     Qtc: Flujo a través de la válvula tricúspide
%     Vrv: Volumen en el ventrículo derecho
%     Qpv: Flujo a través de la válvula pulmonar
%     Vpa: Volumen en la arteria pulmonar (se trata como una cámara elastica)
%% 1) ARRANCO CON VI y VD EN LLENADO
    %por lo tanto: Qpv=Qav=0 (no hay flujos de salida de los ventrículos)
%     Qsys=(Pao-Pvc)/Rsys
%     Qpul=(Ppa-Ppu)/Rpul  
    % por lo tanto la variable de estado va a ser:
%     Vpu_p=Qpul-Qmt
%     Qmt_p=(Ppu-Plv-Qmt*Rmt)/Lmt
%     Vlv_p=Qmt
%     Vao_p=-Qsys
%     Vvc_p=Qsys-Qtc
%     Qtc_p=(Pvc-Prv-Qtc*Rtc)/Ltc
%     Vrv_p=Qtc
%     Vpa_p=-Qpul
%     x=[Vpu, Qmt, Vlv, Vao, Vvc, Qtc, Vrv, Vpa]
    Tf=Ti+Tspan;
    [t1,x1,te,xe,~]=ode15s(@sistema1,[Ti Tf],[Vi_pu Qi_mt Vi_lv Vi_ao Vi_vc Qi_tc Vi_rv Vi_pa],options1); %uso ode15s en vez de ode45 xq la señal tiene menos ruidito. Uso ode15s xq en la tesis del autor usa esa herramienta
    T=[T;t1]; %vector temporal
    Vpu=[Vpu; x1(:,1)]; %volumen en m^3 
    Qmt=[Qmt; x1(:,2)]; %flujo de entrada al VI
    Vlv=[Vlv; x1(:,3)]; %volumen en m^3   
    Vao=[Vao; x1(:,4)]; %volumen en m^3
    Vvc=[Vvc; x1(:,5)]; %volumen en m^3
    Qtc=[Qtc; x1(:,6)]; %flujo de entrada al VD
    Vrv=[Vrv; x1(:,7)]; %volumen en m^3
    Vpa=[Vpa; x1(:,8)]; %volumen en m^3
    Qpv=[Qpv; zeros(size(t1))]; %flujo de salida del VD
    Qav=[Qav; zeros(size(t1))]; %flujo de salida del VI
%     elastancia=e(T);
%     Plv=elastancia.*Ees.*(Vlv-Vd_lv)+(1-elastancia).*P0.*(exp(lamda.*(V-V0))-1);

if(ploteo)
for k=1:length(T)
[PRV(k,1), PLV(k,1), VSPT(k,1)]=PresionesVentr(T(k),Vrv(k),Vlv(k));
% PVC(k,1)=Ees_vc*e(T(k))*Vvc(k);
% PPA(k,1)=Ees_pa*e(T(k))*Vpa(k)+Pth;
% PPU(k,1)=Ees_pu*e(T(k))*Vpu(k)+Pth;
% PAO(k,1)=Ees_ao*e(T(k))*Vao(k);
%la elastancia no es variable==> es un comportamiento pasivo de la camara
%elastica
PVC(k,1)=Ees_vc*Vvc(k);
PPA(k,1)=Ees_pa*Vpa(k)+Pth(T(k));
PPU(k,1)=Ees_pu*Vpu(k)+Pth(T(k));
PAO(k,1)=Ees_ao*Vao(k);
end
figure();
subplot(321),plot(T,Qtc*10^6,T,Qpv*10^6,'r'),legend('Flujo de entrada Qtc','Flujo de salida Qpv'),xlabel('tiempo [s]'),ylabel('flujos del VD [mL/s]');
subplot(322),plot(T,Qmt*10^6,T,Qav*10^6,'r'),legend('Flujo de entrada Qmt','Flujo de salida Qav'),xlabel('tiempo [s]'),ylabel('flujos del VI [mL/s]');
subplot(323),plot(T,Vrv*10^6),xlabel('tiempo [s]'),ylabel('volumen VD [mL]');
subplot(324),plot(T,Vlv*10^6),xlabel('tiempo [s]'),ylabel('volumen VI [mL]');
subplot(325),plot(T,PRV/133.322,T,PVC/133.322,T,PPA/133.322),legend('Presion ventricular','Presion vena cava','Presion arteria pulmonar'),xlabel('tiempo [s]'),ylabel('presion VD [mL]');
subplot(326),plot(T,PLV/133.322,T,PPU/133.322,T,PAO/133.322),legend('Presion ventricular','Presion vena pulmonar','Presion aorta'),xlabel('tiempo [s]'),ylabel('presion VI [mmHg]');
% subplot(421),plot(T,Qtc*10^6,T,Qpv*10^6,'r'),legend('Flujo de entrada Qtc','Flujo de salida Qpv'),xlabel('tiempo [s]'),ylabel('flujos del VD [mL/s]');
% subplot(422),plot(T,Qmt*10^6,T,Qav*10^6,'r'),legend('Flujo de entrada Qmt','Flujo de salida Qav'),xlabel('tiempo [s]'),ylabel('flujos del VI [mL/s]');
% subplot(423),plot(T,Vrv*10^6),xlabel('tiempo [s]'),ylabel('volumen VD [mL]');
% subplot(424),plot(T,Vlv*10^6),xlabel('tiempo [s]'),ylabel('volumen VI [mL]');
% subplot(425),plot(T,PRV/133.322,T,PVC/133.322,T,PPA/133.322),legend('Presion ventricular','Presion aurícula der','Presion arteria pulmonar'),xlabel('tiempo [s]'),ylabel('presion VD [mL]');
% subplot(426),plot(T,PLV/133.322,T,PPU/133.322,T,PAO/133.322),legend('Presion ventricular','Presion aurícula izq','Presion aorta'),xlabel('tiempo [s]'),ylabel('presion VI [mmHg]');
% subplot(427),plot(T,Vvc*10^6,T,Vao*10^6),xlabel('tiempo [s]'),ylabel('volumenes [mL]'),legend('Volumen vena cava','Volumen aorta');
% subplot(428),plot(T,Vpa*10^6,T,Vpu*10^6),xlabel('tiempo [s]'),ylabel('volumenes [mL]'),legend('Volumen art pulmonar','Volumen vena pulmonar');
figure();
plot(T,VSPT*10^6), xlabel('tiempo [s]'),ylabel('Vol del septum [mL]');
suptitle('fase 1');    
end
    %% 2) ARRANCA ISOVOL SISTOLICA DEL VI
    % MIENTRAS EL VD SIGUE EN LLENADO
    %por lo tanto: Qpv=Qav=Qmt=0 (no hay flujos de salida de los ventrículos ni de entrada al ventriculo izquierdo)
%     Qsys=(Pao-Pvc)/Rsys
%     Qpul=(Ppa-Ppu)/Rpul  
    % por lo tanto la variable de estado va a ser:
%     Vpu_p=Qpul
%     Vlv_p=0
%     Vao_p=-Qsys
%     Vvc_p=Qsys-Qtc
%     Qtc_p=(Pvc-Prv-Qtc*Rtc)/Ltc
%     Vrv_p=Qtc
%     Vpa_p=-Qpul
%     x=[Vpu, Vlv, Vao, Vvc, Qtc, Vrv, Vpa]
    Tf=Ti+Tspan;
    [t1,x1,te,xe,~]=ode15s(@sistema2,[te Tf],[xe(1) xe(3) xe(4) xe(5) xe(6) xe(7) xe(8)],options2); %uso ode15s en vez de ode45 xq la señal tiene menos ruidito. Uso ode15s xq en la tesis del autor usa esa herramienta
    T=[T;t1]; %vector temporal
    Vpu=[Vpu; x1(:,1)]; %volumen en m^3 
    Qmt=[Qmt; zeros(size(t1))]; %flujo de entrada al VI
    Vlv=[Vlv; x1(:,2)]; %volumen en m^3   
    Vao=[Vao; x1(:,3)]; %volumen en m^3
    Vvc=[Vvc; x1(:,4)]; %volumen en m^3
    Qtc=[Qtc; x1(:,5)]; %flujo de entrada al VD
    Vrv=[Vrv; x1(:,6)]; %volumen en m^3
    Vpa=[Vpa; x1(:,7)]; %volumen en m^3
    Qpv=[Qpv; zeros(size(t1))]; %flujo de salida del VD
    Qav=[Qav; zeros(size(t1))]; %flujo de salida del VI
%     elastancia=e(T);
%     Plv=elastancia.*Ees.*(Vlv-Vd_lv)+(1-elastancia).*P0.*(exp(lamda.*(V-V0))-1);

if(ploteo)
for k=1:length(T)
[PRV(k,1), PLV(k,1), VSPT(k,1)]=PresionesVentr(T(k),Vrv(k),Vlv(k));
% PVC(k,1)=Ees_vc*e(T(k))*Vvc(k);
% PPA(k,1)=Ees_pa*e(T(k))*Vpa(k)+Pth;
% PPU(k,1)=Ees_pu*e(T(k))*Vpu(k)+Pth;
% PAO(k,1)=Ees_ao*e(T(k))*Vao(k);
%la elastancia no es variable==> es un comportamiento pasivo de la camara
%elastica
PVC(k,1)=Ees_vc*Vvc(k);
PPA(k,1)=Ees_pa*Vpa(k)+Pth(T(k));
PPU(k,1)=Ees_pu*Vpu(k)+Pth(T(k));
PAO(k,1)=Ees_ao*Vao(k);
end
figure();
subplot(321),plot(T,Qtc*10^6,T,Qpv*10^6,'r'),legend('Flujo de entrada Qtc','Flujo de salida Qpv'),xlabel('tiempo [s]'),ylabel('flujos del VD [mL/s]');
subplot(322),plot(T,Qmt*10^6,T,Qav*10^6,'r'),legend('Flujo de entrada Qmt','Flujo de salida Qav'),xlabel('tiempo [s]'),ylabel('flujos del VI [mL/s]');
subplot(323),plot(T,Vrv*10^6),xlabel('tiempo [s]'),ylabel('volumen VD [mL]');
subplot(324),plot(T,Vlv*10^6),xlabel('tiempo [s]'),ylabel('volumen VI [mL]');
subplot(325),plot(T,PRV/133.322,T,PVC/133.322,T,PPA/133.322),legend('Presion ventricular','Presion vena cava','Presion arteria pulmonar'),xlabel('tiempo [s]'),ylabel('presion VD [mL]');
subplot(326),plot(T,PLV/133.322,T,PPU/133.322,T,PAO/133.322),legend('Presion ventricular','Presion vena pulmonar','Presion aorta'),xlabel('tiempo [s]'),ylabel('presion VI [mmHg]');
figure();
plot(T,VSPT*10^6), xlabel('tiempo [s]'),ylabel('Vol del septum [mL]');
suptitle('fase 2');    
end
    %% 3) ARRANCA CON ISOVOL SISTOLICA DEL VD
    % MIENTRAS EL VD TAMBIEN ESTA EN ISOVOL SISTOLICA
    %por lo tanto: Qpv=Qav=Qmt=Qtc=0 (no hay flujos de salida ni de entrada a los ventrículos)
%     Qsys=(Pao-Pvc)/Rsys
%     Qpul=(Ppa-Ppu)/Rpul  
    % por lo tanto la variable de estado va a ser:
%     Vpu_p=Qpul
%     Vlv_p=0
%     Vao_p=-Qsys
%     Vvc_p=Qsys
%     Vrv_p=0
%     Vpa_p=-Qpul
%     x=[Vpu, Vlv, Vao, Vvc, Vrv, Vpa]
    Tf=te+Tspan;
    [t1,x1,te,xe,~]=ode15s(@sistema3,[te Tf],[xe(1) xe(2) xe(3) xe(4) xe(6) xe(7)],options3); %uso ode15s en vez de ode45 xq la señal tiene menos ruidito. Uso ode15s xq en la tesis del autor usa esa herramienta
    T=[T;t1]; %vector temporal
    Vpu=[Vpu; x1(:,1)]; %volumen en m^3 
    Qmt=[Qmt; zeros(size(t1))]; %flujo de entrada al VI
    Vlv=[Vlv; x1(:,2)]; %volumen en m^3   
    Vao=[Vao; x1(:,3)]; %volumen en m^3
    Vvc=[Vvc; x1(:,4)]; %volumen en m^3
    Qtc=[Qtc; zeros(size(t1))]; %flujo de entrada al VD
    Vrv=[Vrv; x1(:,5)]; %volumen en m^3
    Vpa=[Vpa; x1(:,6)]; %volumen en m^3
    Qpv=[Qpv; zeros(size(t1))]; %flujo de salida del VD
    Qav=[Qav; zeros(size(t1))]; %flujo de salida del VI
%     elastancia=e(T);
%     Plv=elastancia.*Ees.*(Vlv-Vd_lv)+(1-elastancia).*P0.*(exp(lamda.*(V-V0))-1);

if(ploteo)
for k=1:length(T)
[PRV(k,1), PLV(k,1), VSPT(k,1)]=PresionesVentr(T(k),Vrv(k),Vlv(k));
% PVC(k,1)=Ees_vc*e(T(k))*Vvc(k);
% PPA(k,1)=Ees_pa*e(T(k))*Vpa(k)+Pth;
% PPU(k,1)=Ees_pu*e(T(k))*Vpu(k)+Pth;
% PAO(k,1)=Ees_ao*e(T(k))*Vao(k);
%la elastancia no es variable==> es un comportamiento pasivo de la camara
%elastica
PVC(k,1)=Ees_vc*Vvc(k);
PPA(k,1)=Ees_pa*Vpa(k)+Pth(T(k));
PPU(k,1)=Ees_pu*Vpu(k)+Pth(T(k));
PAO(k,1)=Ees_ao*Vao(k);
end
figure();
subplot(321),plot(T,Qtc*10^6,T,Qpv*10^6,'r'),legend('Flujo de entrada Qtc','Flujo de salida Qpv'),xlabel('tiempo [s]'),ylabel('flujos del VD [mL/s]');
subplot(322),plot(T,Qmt*10^6,T,Qav*10^6,'r'),legend('Flujo de entrada Qmt','Flujo de salida Qav'),xlabel('tiempo [s]'),ylabel('flujos del VI [mL/s]');
subplot(323),plot(T,Vrv*10^6),xlabel('tiempo [s]'),ylabel('volumen VD [mL]');
subplot(324),plot(T,Vlv*10^6),xlabel('tiempo [s]'),ylabel('volumen VI [mL]');
subplot(325),plot(T,PRV/133.322,T,PVC/133.322,T,PPA/133.322),legend('Presion ventricular','Presion vena cava','Presion arteria pulmonar'),xlabel('tiempo [s]'),ylabel('presion VD [mL]');
subplot(326),plot(T,PLV/133.322,T,PPU/133.322,T,PAO/133.322),legend('Presion ventricular','Presion vena pulmonar','Presion aorta'),xlabel('tiempo [s]'),ylabel('presion VI [mmHg]');
figure();
plot(T,VSPT*10^6), xlabel('tiempo [s]'),ylabel('Vol del septum [mL]');
suptitle('fase 3');    
end
    %% 4) ARRANCA EYECCIÓN DEL VD
    % MIENTRAS EL VI SIGUE EN ISOVOL SISTOLICA
    %por lo tanto: Qav=Qmt=Qtc=0 (no hay flujos de entrada a los ventrículos ni de salida del VI)
%     Qsys=(Pao-Pvc)/Rsys
%     Qpul=(Ppa-Ppu)/Rpul  
    % por lo tanto la variable de estado va a ser:
%     Vpu_p=Qpul
%     Vlv_p=0
%     Vao_p=-Qsys
%     Vvc_p=Qsys
%     Vrv_p=0-Qpv
%     Qpv_p=(Prv-Ppa-Qpv.Rpv)/Lpv
%     Vpa_p=Qpv-Qpul
%     x=[Vpu, Vlv, Vao, Vvc, Vrv, Qpv, Vpa]
    Tf=te+Tspan;
    [t1,x1,te,xe,~]=ode15s(@sistema4,[te Tf],[xe(1) xe(2) xe(3) xe(4) xe(5) 0 xe(6)],options4); %uso ode15s en vez de ode45 xq la señal tiene menos ruidito. Uso ode15s xq en la tesis del autor usa esa herramienta
    T=[T;t1]; %vector temporal
    Vpu=[Vpu; x1(:,1)]; %volumen en m^3 
    Vlv=[Vlv; x1(:,2)]; %volumen en m^3   
    Vao=[Vao; x1(:,3)]; %volumen en m^3
    Vvc=[Vvc; x1(:,4)]; %volumen en m^3
    Vrv=[Vrv; x1(:,5)]; %volumen en m^3
    Qpv=[Qpv; x1(:,6)]; %flujo de salida del VD
    Vpa=[Vpa; x1(:,7)]; %volumen en m^3
    Qav=[Qav; zeros(size(t1))]; %flujo de salida del VI
    Qmt=[Qmt; zeros(size(t1))]; %flujo de entrada al VI
    Qtc=[Qtc; zeros(size(t1))]; %flujo de entrada al VD
%     elastancia=e(T);
%     Plv=elastancia.*Ees.*(Vlv-Vd_lv)+(1-elastancia).*P0.*(exp(lamda.*(V-V0))-1);

if(ploteo)
for k=1:length(T)
[PRV(k,1), PLV(k,1), VSPT(k,1)]=PresionesVentr(T(k),Vrv(k),Vlv(k));
% PVC(k,1)=Ees_vc*e(T(k))*Vvc(k);
% PPA(k,1)=Ees_pa*e(T(k))*Vpa(k)+Pth;
% PPU(k,1)=Ees_pu*e(T(k))*Vpu(k)+Pth;
% PAO(k,1)=Ees_ao*e(T(k))*Vao(k);
%la elastancia no es variable==> es un comportamiento pasivo de la camara
%elastica
PVC(k,1)=Ees_vc*Vvc(k);
PPA(k,1)=Ees_pa*Vpa(k)+Pth(T(k));
PPU(k,1)=Ees_pu*Vpu(k)+Pth(T(k));
PAO(k,1)=Ees_ao*Vao(k);
end
figure();
subplot(321),plot(T,Qtc*10^6,T,Qpv*10^6,'r'),legend('Flujo de entrada Qtc','Flujo de salida Qpv'),xlabel('tiempo [s]'),ylabel('flujos del VD [mL/s]');
subplot(322),plot(T,Qmt*10^6,T,Qav*10^6,'r'),legend('Flujo de entrada Qmt','Flujo de salida Qav'),xlabel('tiempo [s]'),ylabel('flujos del VI [mL/s]');
subplot(323),plot(T,Vrv*10^6),xlabel('tiempo [s]'),ylabel('volumen VD [mL]');
subplot(324),plot(T,Vlv*10^6),xlabel('tiempo [s]'),ylabel('volumen VI [mL]');
subplot(325),plot(T,PRV/133.322,T,PVC/133.322,T,PPA/133.322),legend('Presion ventricular','Presion vena cava','Presion arteria pulmonar'),xlabel('tiempo [s]'),ylabel('presion VD [mL]');
subplot(326),plot(T,PLV/133.322,T,PPU/133.322,T,PAO/133.322),legend('Presion ventricular','Presion vena pulmonar','Presion aorta'),xlabel('tiempo [s]'),ylabel('presion VI [mmHg]');
figure();
plot(T,VSPT*10^6), xlabel('tiempo [s]'),ylabel('Vol del septum [mL]');
suptitle('fase 4');    
end
    %% 5) ARRANCA EYECCION DEL VI
    % MIENTRAS EL VD SIGUE EN EYECCION
    %por lo tanto: Qmt=Qtc=0 (no hay flujos de entrada a los ventrículos)
%     Qsys=(Pao-Pvc)/Rsys
%     Qpul=(Ppa-Ppu)/Rpul  
    % por lo tanto la variable de estado va a ser:
%     Vpu_p=Qpul
%     Vlv_p=0-Qav
%     Qav_p=(Plv-Pao-Qav*Rav)/Lav;
%     Vao_p=Qav-Qsys
%     Vvc_p=Qsys
%     Vrv_p=0-Qpv
%     Qpv_p=(Prv-Ppa-Qpv.Rpv)/Lpv
%     Vpa_p=Qpv-Qpul
%     x=[Vpu, Vlv, Qav, Vao, Vvc, Vrv, Qpv, Vpa]
    Tf=te+Tspan;
    [t1,x1,te,xe,~]=ode15s(@sistema5,[te Tf],[xe(1) xe(2) 0 xe(3) xe(4) xe(5) xe(6) xe(7)],options5); %uso ode15s en vez de ode45 xq la señal tiene menos ruidito. Uso ode15s xq en la tesis del autor usa esa herramienta
    T=[T;t1]; %vector temporal
    Vpu=[Vpu; x1(:,1)]; %volumen en m^3 
    Vlv=[Vlv; x1(:,2)]; %volumen en m^3  
    Qav=[Qav; x1(:,3)]; %flujo de salida del VI
    Vao=[Vao; x1(:,4)]; %volumen en m^3
    Vvc=[Vvc; x1(:,5)]; %volumen en m^3
    Vrv=[Vrv; x1(:,6)]; %volumen en m^3
    Qpv=[Qpv; x1(:,7)]; %flujo de salida del VD
    Vpa=[Vpa; x1(:,8)]; %volumen en m^3
    Qmt=[Qmt; zeros(size(t1))]; %flujo de entrada al VI
    Qtc=[Qtc; zeros(size(t1))]; %flujo de entrada al VD
%     elastancia=e(T);
%     Plv=elastancia.*Ees.*(Vlv-Vd_lv)+(1-elastancia).*P0.*(exp(lamda.*(V-V0))-1);
if(ploteo)
for k=1:length(T)
[PRV(k,1), PLV(k,1), VSPT(k,1)]=PresionesVentr(T(k),Vrv(k),Vlv(k));
% PVC(k,1)=Ees_vc*e(T(k))*Vvc(k);
% PPA(k,1)=Ees_pa*e(T(k))*Vpa(k)+Pth;
% PPU(k,1)=Ees_pu*e(T(k))*Vpu(k)+Pth;
% PAO(k,1)=Ees_ao*e(T(k))*Vao(k);
%la elastancia no es variable==> es un comportamiento pasivo de la camara
%elastica
PVC(k,1)=Ees_vc*Vvc(k);
PPA(k,1)=Ees_pa*Vpa(k)+Pth(T(k));
PPU(k,1)=Ees_pu*Vpu(k)+Pth(T(k));
PAO(k,1)=Ees_ao*Vao(k);
end
figure();
subplot(321),plot(T,Qtc*10^6,T,Qpv*10^6,'r'),legend('Flujo de entrada Qtc','Flujo de salida Qpv'),xlabel('tiempo [s]'),ylabel('flujos del VD [mL/s]');
subplot(322),plot(T,Qmt*10^6,T,Qav*10^6,'r'),legend('Flujo de entrada Qmt','Flujo de salida Qav'),xlabel('tiempo [s]'),ylabel('flujos del VI [mL/s]');
subplot(323),plot(T,Vrv*10^6),xlabel('tiempo [s]'),ylabel('volumen VD [mL]');
subplot(324),plot(T,Vlv*10^6),xlabel('tiempo [s]'),ylabel('volumen VI [mL]');
subplot(325),plot(T,PRV/133.322,T,PVC/133.322,T,PPA/133.322),legend('Presion ventricular','Presion vena cava','Presion arteria pulmonar'),xlabel('tiempo [s]'),ylabel('presion VD [mL]');
subplot(326),plot(T,PLV/133.322,T,PPU/133.322,T,PAO/133.322),legend('Presion ventricular','Presion vena pulmonar','Presion aorta'),xlabel('tiempo [s]'),ylabel('presion VI [mmHg]');
figure();
plot(T,VSPT*10^6), xlabel('tiempo [s]'),ylabel('Vol del septum [mL]');
suptitle('fase 5');  
end
 %% 6) ARRANCA ISOVOL DIAST DEL VI 
    % MIENTRAS EL VD SIGUE EN EYECCION
    %por lo tanto: Qmt=Qtc=Qav=0 (no hay flujos de entrada a los ventrículos ni de salida del VI)
%     Qsys=(Pao-Pvc)/Rsys
%     Qpul=(Ppa-Ppu)/Rpul  
    % por lo tanto la variable de estado va a ser:
%     Vpu_p=Qpul
%     Vlv_p=0
%     Vao_p=-Qsys
%     Vvc_p=Qsys
%     Vrv_p=-Qpv
%     Qpv_p=(Prv-Ppa-Qav*Rpv)/Lpv;
%     Vpa_p=Qpv-Qpul
%     x=[Vpu, Vlv, Vao, Vvc, Vrv, Qpv, Vpa]
    Tf=te+Tspan;
    [t1,x1,te,xe,~]=ode15s(@sistema6,[te Tf],[xe(1) xe(2) xe(4) xe(5) xe(6) xe(7) xe(8)],options6); %uso ode15s en vez de ode45 xq la señal tiene menos ruidito. Uso ode15s xq en la tesis del autor usa esa herramienta
    T=[T;t1]; %vector temporal
    Vpu=[Vpu; x1(:,1)]; %volumen en m^3 
    Vlv=[Vlv; x1(:,2)]; %volumen en m^3  
    Qav=[Qav; zeros(size(t1))]; %flujo de salida del VI
    Vao=[Vao; x1(:,3)]; %volumen en m^3
    Vvc=[Vvc; x1(:,4)]; %volumen en m^3
    Vrv=[Vrv; x1(:,5)]; %volumen en m^3
    Qpv=[Qpv; x1(:,6)]; %flujo de salida del VD
    Vpa=[Vpa; x1(:,7)]; %volumen en m^3
    Qmt=[Qmt; zeros(size(t1))]; %flujo de entrada al VI
    Qtc=[Qtc; zeros(size(t1))]; %flujo de entrada al VD
%     elastancia=e(T);
%     Plv=elastancia.*Ees.*(Vlv-Vd_lv)+(1-elastancia).*P0.*(exp(lamda.*(V-V0))-1);
if(ploteo)
for k=1:length(T)
[PRV(k,1), PLV(k,1), VSPT(k,1)]=PresionesVentr(T(k),Vrv(k),Vlv(k));
% PVC(k,1)=Ees_vc*e(T(k))*Vvc(k);
% PPA(k,1)=Ees_pa*e(T(k))*Vpa(k)+Pth;
% PPU(k,1)=Ees_pu*e(T(k))*Vpu(k)+Pth;
% PAO(k,1)=Ees_ao*e(T(k))*Vao(k);
%la elastancia no es variable==> es un comportamiento pasivo de la camara
%elastica
PVC(k,1)=Ees_vc*Vvc(k);
PPA(k,1)=Ees_pa*Vpa(k)+Pth(T(k));
PPU(k,1)=Ees_pu*Vpu(k)+Pth(T(k));
PAO(k,1)=Ees_ao*Vao(k);
end
figure();
subplot(321),plot(T,Qtc*10^6,T,Qpv*10^6,'r'),legend('Flujo de entrada Qtc','Flujo de salida Qpv'),xlabel('tiempo [s]'),ylabel('flujos del VD [mL/s]');
subplot(322),plot(T,Qmt*10^6,T,Qav*10^6,'r'),legend('Flujo de entrada Qmt','Flujo de salida Qav'),xlabel('tiempo [s]'),ylabel('flujos del VI [mL/s]');
subplot(323),plot(T,Vrv*10^6),xlabel('tiempo [s]'),ylabel('volumen VD [mL]');
subplot(324),plot(T,Vlv*10^6),xlabel('tiempo [s]'),ylabel('volumen VI [mL]');
subplot(325),plot(T,PRV/133.322,T,PVC/133.322,T,PPA/133.322),legend('Presion ventricular','Presion vena cava','Presion arteria pulmonar'),xlabel('tiempo [s]'),ylabel('presion VD [mL]');
subplot(326),plot(T,PLV/133.322,T,PPU/133.322,T,PAO/133.322),legend('Presion ventricular','Presion vena pulmonar','Presion aorta'),xlabel('tiempo [s]'),ylabel('presion VI [mmHg]');
figure();
plot(T,VSPT*10^6), xlabel('tiempo [s]'),ylabel('Vol del septum [mL]');
suptitle('fase 6'); 
end

%% 7) ARRANCA ISOVOL DIAST DEL VD
    % MIENTRAS EL VI SIGUE EN isovol sist
   %por lo tanto: Qmt=Qtc=Qpv=Qav=0 (no hay flujos de entrada ni de salida de los ventrículos)
%     Qsys=(Pao-Pvc)/Rsys
%     Qpul=(Ppa-Ppu)/Rpul  
    % por lo tanto la variable de estado va a ser:
%     Vpu_p=Qpul
%     Vlv_p=0
%     Vao_p=-Qsys
%     Vvc_p=Qsys
%     Vrv_p=0
%     Vpa_p=0-Qpul
%     x=[Vpu, Vlv, Vao, Vvc, Vrv, Vpa]
    Tf=te+Tspan;
    [t1,x1,te,xe,~]=ode15s(@sistema7,[te Tf],[xe(1) xe(2) xe(3) xe(4) xe(5) xe(7)],options7); %uso ode15s en vez de ode45 xq la señal tiene menos ruidito. Uso ode15s xq en la tesis del autor usa esa herramienta
    T=[T;t1]; %vector temporal
    Vpu=[Vpu; x1(:,1)]; %volumen en m^3 
    Vlv=[Vlv; x1(:,2)]; %volumen en m^3  
    Qav=[Qav; zeros(size(t1))]; %flujo de salida del VI
    Vao=[Vao; x1(:,3)]; %volumen en m^3
    Vvc=[Vvc; x1(:,4)]; %volumen en m^3
    Vrv=[Vrv; x1(:,5)]; %volumen en m^3
    Qpv=[Qpv; zeros(size(t1))]; %flujo de salida del VD
    Vpa=[Vpa; x1(:,6)]; %volumen en m^3
    Qmt=[Qmt; zeros(size(t1))]; %flujo de entrada al VI
    Qtc=[Qtc; zeros(size(t1))]; %flujo de entrada al VD
%     elastancia=e(T);
%     Plv=elastancia.*Ees.*(Vlv-Vd_lv)+(1-elastancia).*P0.*(exp(lamda.*(V-V0))-1);
if(ploteo)
for k=1:length(T)
[PRV(k,1), PLV(k,1), VSPT(k,1)]=PresionesVentr(T(k),Vrv(k),Vlv(k));
% PVC(k,1)=Ees_vc*e(T(k))*Vvc(k);
% PPA(k,1)=Ees_pa*e(T(k))*Vpa(k)+Pth;
% PPU(k,1)=Ees_pu*e(T(k))*Vpu(k)+Pth;
% PAO(k,1)=Ees_ao*e(T(k))*Vao(k);
%la elastancia no es variable==> es un comportamiento pasivo de la camara
%elastica
PVC(k,1)=Ees_vc*Vvc(k);
PPA(k,1)=Ees_pa*Vpa(k)+Pth(T(k));
PPU(k,1)=Ees_pu*Vpu(k)+Pth(T(k));
PAO(k,1)=Ees_ao*Vao(k);
end
figure();
subplot(321),plot(T,Qtc*10^6,T,Qpv*10^6,'r'),legend('Flujo de entrada Qtc','Flujo de salida Qpv'),xlabel('tiempo [s]'),ylabel('flujos del VD [mL/s]');
subplot(322),plot(T,Qmt*10^6,T,Qav*10^6,'r'),legend('Flujo de entrada Qmt','Flujo de salida Qav'),xlabel('tiempo [s]'),ylabel('flujos del VI [mL/s]');
subplot(323),plot(T,Vrv*10^6),xlabel('tiempo [s]'),ylabel('volumen VD [mL]');
subplot(324),plot(T,Vlv*10^6),xlabel('tiempo [s]'),ylabel('volumen VI [mL]');
subplot(325),plot(T,PRV/133.322,T,PVC/133.322,T,PPA/133.322),legend('Presion ventricular','Presion vena cava','Presion arteria pulmonar'),xlabel('tiempo [s]'),ylabel('presion VD [mL]');
subplot(326),plot(T,PLV/133.322,T,PPU/133.322,T,PAO/133.322),legend('Presion ventricular','Presion vena pulmonar','Presion aorta'),xlabel('tiempo [s]'),ylabel('presion VI [mmHg]');
figure();
plot(T,VSPT*10^6), xlabel('tiempo [s]'),ylabel('Vol del septum [mL]');
suptitle('fase 7');  
end

%% 8) ARRANCA LLENADO DEL VD 
    % MIENTRAS EL VI SIGUE EN isovol diast
   %por lo tanto: Qmt=Qpv=Qav=0 (no hay flujo de salida de los ventrículos ni de entrada al VI)
%     Qsys=(Pao-Pvc)/Rsys
%     Qpul=(Ppa-Ppu)/Rpul  
    % por lo tanto la variable de estado va a ser:
%     Vpu_p=Qpul
%     Vlv_p=0
%     Vao_p=-Qsys
%     Vvc_p=Qsys-Qtc
%     Qtc_p=(Pvc-Prv-Qtc*Rtc)/Ltc
%     Vrv_p=Qtc-0
%     Vpa_p=0-Qpul
%     x=[Vpu, Vlv, Vao, Vvc, Qtc, Vrv, Vpa]
    Tf=te+Tspan;
    [t1,x1,te,xe,~]=ode15s(@sistema8,[te Tf],[xe(1) xe(2) xe(3) xe(4) 0 xe(5) xe(6)],options8); %uso ode15s en vez de ode45 xq la señal tiene menos ruidito. Uso ode15s xq en la tesis del autor usa esa herramienta
    T=[T;t1]; %vector temporal
    Vpu=[Vpu; x1(:,1)]; %volumen en m^3 
    Vlv=[Vlv; x1(:,2)]; %volumen en m^3  
    Qav=[Qav; zeros(size(t1))]; %flujo de salida del VI
    Vao=[Vao; x1(:,3)]; %volumen en m^3
    Vvc=[Vvc; x1(:,4)]; %volumen en m^3
    Qtc=[Qtc; x1(:,5)]; %flujo de entrada al VD
    Vrv=[Vrv; x1(:,6)]; %volumen en m^3
    Qpv=[Qpv; zeros(size(t1))]; %flujo de salida del VD
    Vpa=[Vpa; x1(:,7)]; %volumen en m^3
    Qmt=[Qmt; zeros(size(t1))]; %flujo de entrada al VI
%     elastancia=e(T);
%     Plv=elastancia.*Ees.*(Vlv-Vd_lv)+(1-elastancia).*P0.*(exp(lamda.*(V-V0))-1);
if(ploteo)
for k=1:length(T)
[PRV(k,1), PLV(k,1), VSPT(k,1)]=PresionesVentr(T(k),Vrv(k),Vlv(k));
% PVC(k,1)=Ees_vc*e(T(k))*Vvc(k);
% PPA(k,1)=Ees_pa*e(T(k))*Vpa(k)+Pth;
% PPU(k,1)=Ees_pu*e(T(k))*Vpu(k)+Pth;
% PAO(k,1)=Ees_ao*e(T(k))*Vao(k);
%la elastancia no es variable==> es un comportamiento pasivo de la camara
%elastica
PVC(k,1)=Ees_vc*Vvc(k);
PPA(k,1)=Ees_pa*Vpa(k)+Pth(T(k));
PPU(k,1)=Ees_pu*Vpu(k)+Pth(T(k));
PAO(k,1)=Ees_ao*Vao(k);
end
figure();
subplot(321),plot(T,Qtc*10^6,T,Qpv*10^6,'r'),legend('Flujo de entrada Qtc','Flujo de salida Qpv'),xlabel('tiempo [s]'),ylabel('flujos del VD [mL/s]');
subplot(322),plot(T,Qmt*10^6,T,Qav*10^6,'r'),legend('Flujo de entrada Qmt','Flujo de salida Qav'),xlabel('tiempo [s]'),ylabel('flujos del VI [mL/s]');
subplot(323),plot(T,Vrv*10^6),xlabel('tiempo [s]'),ylabel('volumen VD [mL]');
subplot(324),plot(T,Vlv*10^6),xlabel('tiempo [s]'),ylabel('volumen VI [mL]');
subplot(325),plot(T,PRV/133.322,T,PVC/133.322,T,PPA/133.322),legend('Presion ventricular','Presion vena cava','Presion arteria pulmonar'),xlabel('tiempo [s]'),ylabel('presion VD [mL]');
subplot(326),plot(T,PLV/133.322,T,PPU/133.322,T,PAO/133.322),legend('Presion ventricular','Presion vena pulmonar','Presion aorta'),xlabel('tiempo [s]'),ylabel('presion VI [mmHg]');
figure();
plot(T,VSPT*10^6), xlabel('tiempo [s]'),ylabel('Vol del septum [mL]');
suptitle('fase 8');  
end
x=[Vpu, Vlv, Vao, Vvc, Qtc, Vrv, Vpa]
Vi_pu=xe(1);
Qi_mt=0;
Vi_lv=xe(2);
Vi_ao=xe(3);
Vi_vc=xe(4);
Qi_tc=xe(5);
Vi_rv=xe(6);
Vi_pa=xe(7);
Ti=te;

end
for k=1:length(T)
[PRV(k,1), PLV(k,1),Vspt(k,1)]=PresionesVentr(T(k),Vrv(k),Vlv(k));
% PVC(k,1)=Ees_vc*e(T(k))*Vvc(k);
% PPA(k,1)=Ees_pa*e(T(k))*Vpa(k)+Pth;
% PPU(k,1)=Ees_pu*e(T(k))*Vpu(k)+Pth;
% PAO(k,1)=Ees_ao*e(T(k))*Vao(k);
%la elastancia no es variable==> es un comportamiento pasivo de la camara
%elastica
PVC(k,1)=Ees_vc*Vvc(k);
PPA(k,1)=Ees_pa*Vpa(k)+Pth(T(k));
PPU(k,1)=Ees_pu*Vpu(k)+Pth(T(k));
PAO(k,1)=Ees_ao*Vao(k);
PTH(k,1)=Pth(T(k));
end
VPCD=Vlv+Vrv;
figure();
subplot(321),plot(T,Qtc*10^6,T,Qpv*10^6,'r'),legend('Flujo de entrada Qtc','Flujo de salida Qpv'),xlabel('tiempo [s]'),ylabel('flujos del VD [mL/s]');
subplot(322),plot(T,Qmt*10^6,T,Qav*10^6,'r'),legend('Flujo de entrada Qmt','Flujo de salida Qav'),xlabel('tiempo [s]'),ylabel('flujos del VI [mL/s]');
subplot(323),plot(T,Vrv*10^6),xlabel('tiempo [s]'),ylabel('volumen VD [mL]');
subplot(324),plot(T,Vlv*10^6),xlabel('tiempo [s]'),ylabel('volumen VI [mL]');
subplot(325),plot(T,PRV/133.322,T,PVC/133.322,T,PPA/133.322),legend('Presion ventricular','Presion vena cava','Presion arteria pulmonar'),xlabel('tiempo [s]'),ylabel('presion VD [mmHg]');
subplot(326),plot(T,PLV/133.322,T,PPU/133.322,T,PAO/133.322),legend('Presion ventricular','Presion pulmonar','Presion aorta'),xlabel('tiempo [s]'),ylabel('presion VI [mmHg]');
% suptitle('fase 8');  
figure();
subplot(121),plot(Vlv*10^6,PLV/133.322),xlabel('volumen VI [mL]'),ylabel('Presion VI [mmHg]');
subplot(122),plot(Vrv*10^6,PRV/133.322),xlabel('volumen VD [mL]'),ylabel('Presion VD [mmHg]');
figure();
[AX,~,~]=plotyy(T,[Vrv*10^6 Vlv*10^6],T,PTH/133.322);xlabel('tiempo [s]');%,ylabel('volumenes ventriculares [mL]'); %legend('VD','VI'),grid on;
set(get(AX(1),'Ylabel'),'String','Volumenes ventriculares [mL]'),legend('VD','VI','Pth'); 
set(get(AX(2),'Ylabel'),'String','Presion pleural [mmHg]');
figure();
[AX,~,~]=plotyy(T,Vspt*10^6,T,VPCD*10^6);xlabel('tiempo [s]');
set(get(AX(1),'Ylabel'),'String','Volumen del septum [mL]'); 
set(get(AX(2),'Ylabel'),'String','Volumen del pericardio [mL]');
figure();
Vtot=Vrv+Vlv+Vpu+Vvc+Vao+Vpa; %volumen total en cada iteracion
plot(T,Vtot*10^6),xlabel('tiempo [s]'),ylabel('volumen total del circuito [mL]');
%% 1) ARRANCO CON VI y VD EN LLENADO
        function dxdt=sistema1(t,x) %x=[Vpu, Qmt, Vlv, Vao, Vvc, Qtc, Vrv, Vpa]
    %por lo tanto: Qpv=Qav=0 (no hay flujos de salida de los ventrículos)
%     Qsys=(Pao-Pvc)/Rsys-->if Qsys<=0==>Qsys=0
%     Qpul=(Ppa-Ppu)/Rpul--> if Qpul<=0==>Qpul=0  
    % por lo tanto la variable de estado va a ser:
%     Vpu_p=Qpul-Qmt
%     Qmt_p=(Ppu-Plv-Qmt*Rmt)/Lmt
%     Vlv_p=Qmt
%     Vao_p=-Qsys
%     Vvc_p=Qsys-Qtc
%     Qtc_p=(Pvc-Prv-Qtc*Rtc)/Ltc
%     Vrv_p=Qtc
%     Vpa_p=-Qpul
%     x=[Vpu, Qmt, Vlv, Vao, Vvc, Qtc, Vrv, Vpa]
    
%calculo las presiones
            [Prv, Plv, ~]=PresionesVentr(t,x(7),x(3));
%             Pvc=Ees_vc*e(t)*x(5);
%             Ppa=Ees_pa*e(t)*x(8)+Pth;
%             Ppu=Ees_pu*e(t)*x(1)+Pth;
%             Pao=Ees_ao*e(t)*x(4);
%VC, PA, PU y AO son camaras elasticas con comportamiento pasivo
            Pvc=Ees_vc*x(5);
            Ppa=Ees_pa*x(8)+Pth(t);
            Ppu=Ees_pu*x(1)+Pth(t);
            Pao=Ees_ao*x(4);
%-------------------------------------------------------------------------
            %calculo las derivadas
            Qsys=(Pao-Pvc)/Rsys; %flujo a través de capilares sistémicos
            if Qsys<=0
                Qsys=0;
            end
            Qpul=(Ppa-Ppu)/Rpul; %flujo a través de capilares pulmonares
            if Qpul<=0
                Qpul=0;
            end
            
            Vpu_p=Qpul-x(2);
            Qmt_p=(Ppu-Plv-x(2)*Rmt)/Lmt;
            Vlv_p=x(2);
            Vao_p=-Qsys;
            Vvc_p=Qsys-x(6);
            Qtc_p=(Pvc-Prv-x(6)*Rtc)/Ltc;
            Vrv_p=x(6);
            Vpa_p=-Qpul;

            dxdt=[Vpu_p; Qmt_p; Vlv_p; Vao_p; Vvc_p; Qtc_p; Vrv_p; Vpa_p]; 

        end

        function [cierre stop direction]=event1(t,x) %condicion1 (en que momento se cierra la valvula mitral y el VI entra en contraccion isovol)
            cierre= x(2); %busca en que momento Qmt se anula
            stop= 1; %corta cuando sucede el evento, o sea cuando Qmt se anula
            direction=-1; %agarro solo si se anula pero tiene derivada negativa, o sea pasa por cero para disminuir (no xq se aumenta)
        end

    %% 2) ISOVOL SISTOLICA DEL VI + LLENADO DE VD
        function dxdt=sistema2(t,x) %x=[Vpu, Vlv, Vao, Vvc, Qtc, Vrv, Vpa]
 %por lo tanto: Qpv=Qav=Qmt=0 (no hay flujos de salida de los ventrículos ni de entrada al ventriculo izquierdo)
%     Qsys=(Pao-Pvc)/Rsys-->if Qsys<=0==>Qsys=0
%     Qpul=(Ppa-Ppu)/Rpul-->if Qpul<=0==>Qpul=0  
    % por lo tanto la variable de estado va a ser:
%     Vpu_p=Qpul
%     Vlv_p=0
%     Vao_p=-Qsys
%     Vvc_p=Qsys-Qtc
%     Qtc_p=(Pvc-Prv-Qtc*Rtc)/Ltc
%     Vrv_p=Qtc
%     Vpa_p=-Qpul
%     x=[Vpu, Vlv, Vao, Vvc, Qtc, Vrv, Vpa]

            %calculo las presiones
            [Prv, Plv, ~]=PresionesVentr(t,x(6),x(2));
%             Pvc=Ees_vc*e(t)*x(4);
%             Ppa=Ees_pa*e(t)*x(7)+Pth;
%             Ppu=Ees_pu*e(t)*x(1)+Pth;
%             Pao=Ees_ao*e(t)*x(3);
%VC, PA, PU y AO son camaras elasticas con comportamiento pasivo
            Pvc=Ees_vc*x(4);
            Ppa=Ees_pa*x(7)+Pth(t);
            Ppu=Ees_pu*x(1)+Pth(t);
            Pao=Ees_ao*x(3);
%-------------------------------------------------------------------------
            %calculo las derivadas
            Qsys=(Pao-Pvc)/Rsys; %flujo a través de capilares sistémicos
            if Qsys<=0
                Qsys=0;
            end
            Qpul=(Ppa-Ppu)/Rpul; %flujo a través de capilares pulmonares
            if Qpul<=0
                Qpul=0;
            end
            
            Vpu_p=Qpul;
            Vlv_p=0;
            Vao_p=-Qsys;
            Vvc_p=Qsys-x(5);
            Qtc_p=(Pvc-Prv-x(5)*Rtc)/Ltc;
            Vrv_p=x(5);
            Vpa_p=-Qpul;

            dxdt=[Vpu_p; Vlv_p; Vao_p; Vvc_p; Qtc_p; Vrv_p; Vpa_p]; 
        end

        function [cierre stop direction]=event2(t,x) %condicion1 (en que momento se cierra la valvula tricuspide y el VD entra en contraccion isovol)
            cierre= x(5); %busca en que momento Qtc se anula
            stop= 1; %corta cuando sucede el evento, o sea cuando Qtc se anula
            direction=-1; %agarro solo si se anula pero tiene derivada negativa, o sea pasa por cero para disminuir (no xq se aumenta)
        end
    %% 3) ISOVOL SISTOLICA DEL VI y VD
        function dxdt=sistema3(t,x) %x=[Vpu, Vlv, Vao, Vvc, Vrv, Vpa]
    %por lo tanto: Qpv=Qav=Qmt=Qtc=0 (no hay flujos de salida ni de entrada a los ventrículos)
%     Qsys=(Pao-Pvc)/Rsys
%     Qpul=(Ppa-Ppu)/Rpul  
    % por lo tanto la variable de estado va a ser:
%     Vpu_p=Qpul
%     Vlv_p=0
%     Vao_p=-Qsys
%     Vvc_p=Qsys
%     Vrv_p=0
%     Vpa_p=-Qpul
%     x=[Vpu, Vlv, Vao, Vvc, Vrv, Vpa]

            %calculo las presiones
            [Prv, Plv, ~]=PresionesVentr(t,x(5),x(2));
%             Pvc=Ees_vc*e(t)*x(4);
%             Ppa=Ees_pa*e(t)*x(6)+Pth;
%             Ppu=Ees_pu*e(t)*x(1)+Pth;
%             Pao=Ees_ao*e(t)*x(3);
%VC, PA, PU y AO son camaras elasticas con comportamiento pasivo            
            Pvc=Ees_vc*x(4);
            Ppa=Ees_pa*x(6)+Pth(t);
            Ppu=Ees_pu*x(1)+Pth(t);
            Pao=Ees_ao*x(3);
%-------------------------------------------------------------------------
            %calculo las derivadas
            Qsys=(Pao-Pvc)/Rsys; %flujo a través de capilares sistémicos
            if Qsys<=0
                Qsys=0;
            end
            Qpul=(Ppa-Ppu)/Rpul; %flujo a través de capilares pulmonares
            if Qpul<=0
                Qpul=0;
            end
            
            Vpu_p=Qpul;
            Vlv_p=0;
            Vao_p=-Qsys;
            Vvc_p=Qsys;
            Vrv_p=0;
            Vpa_p=-Qpul;

            dxdt=[Vpu_p; Vlv_p; Vao_p; Vvc_p; Vrv_p; Vpa_p]; 
        end

        function [cierre stop direction]=event3(t,x) %condicion2 (en que momento se abre la valvula pulmonar y el VD empieza a eyectar)
            [Prv, ~, ~]=PresionesVentr(t,x(5),x(2));
            %VC, PA, PU y AO son camaras elasticas con comportamiento pasivo  
%             Ppa=Ees_pa*e(t)*x(6)+Pth;
            Ppa=Ees_pa*x(6)+Pth(t);
            cierre= Prv-Ppa; %busca en que momento Prv>Partpulmonar
            stop= 1; %corta cuando sucede el evento
            direction=1; %agarro solo si se anula pero tiene derivada positiva
        end
    %% 4) EYECCIÓN DEL VD + ISOVOL SIST DEL VI
        function dxdt=sistema4(t,x) %x=[Vpu, Vlv, Vao, Vvc, Vrv, Qpv, Vpa]
    %por lo tanto: Qav=Qmt=Qtc=0 (no hay flujos de entrada a los ventrículos ni de salida del VI)
%     Qsys=(Pao-Pvc)/Rsys
%     Qpul=(Ppa-Ppu)/Rpul  
    % por lo tanto la variable de estado va a ser:
%     Vpu_p=Qpul
%     Vlv_p=0
%     Vao_p=-Qsys
%     Vvc_p=Qsys
%     Vrv_p=0-Qpv
%     Qpv_p=(Prv-Ppa-Qpv.Rpv)/Lpv
%     Vpa_p=Qpv-Qpul
%     x=[Vpu, Vlv, Vao, Vvc, Vrv, Qpv, Vpa]

            %calculo las presiones
            [Prv, Plv, ~]=PresionesVentr(t,x(5),x(2));
%             Pvc=Ees_vc*e(t)*x(4);
%             Ppa=Ees_pa*e(t)*x(7)+Pth;
%             Ppu=Ees_pu*e(t)*x(1)+Pth;
%             Pao=Ees_ao*e(t)*x(3);
            %VC, PA, PU y AO son camaras elasticas con comportamiento pasivo  
            Pvc=Ees_vc*x(4);
            Ppa=Ees_pa*x(7)+Pth(t);
            Ppu=Ees_pu*x(1)+Pth(t);
            Pao=Ees_ao*x(3);
%-------------------------------------------------------------------------
            %calculo las derivadas
            Qsys=(Pao-Pvc)/Rsys; %flujo a través de capilares sistémicos
            if Qsys<=0
                Qsys=0;
            end
            Qpul=(Ppa-Ppu)/Rpul; %flujo a través de capilares pulmonares
            if Qpul<=0
                Qpul=0;
            end
            
            Vpu_p=Qpul;
            Vlv_p=0;
            Vao_p=-Qsys;
            Vvc_p=Qsys;
            Vrv_p=-x(6);
            Qpv_p=(Prv-Ppa-x(6)*Rpv)/Lpv;
            Vpa_p=x(6)-Qpul;

            dxdt=[Vpu_p; Vlv_p; Vao_p; Vvc_p; Vrv_p; Qpv_p; Vpa_p]; 
        end

        function [cierre stop direction]=event4(t,x) %condicion3 (en que momento se abre la valvula aortica y el VI empieza a eyectar)
            [~, Plv, ~]=PresionesVentr(t,x(5),x(2));
            %VC, PA, PU y AO son camaras elasticas con comportamiento pasivo  
%             Pao=Ees_ao*e(t)*x(3);
            Pao=Ees_ao*x(3);
            cierre= Plv-Pao; %busca en que momento Plv>Paorta
            stop= 1; %corta cuando sucede el evento
            direction=1; %agarro solo si se anula pero tiene derivada positiva
        end
    %% 5) ARRANCA EYECCION DEL VI
    % MIENTRAS EL VD SIGUE EN EYECCION
        function dxdt=sistema5(t,x) %x=[Vpu, Vlv, Qav, Vao, Vvc, Vrv, Qpv, Vpa]
%por lo tanto: Qmt=Qtc=0 (no hay flujos de entrada a los ventrículos)
%     Qsys=(Pao-Pvc)/Rsys
%     Qpul=(Ppa-Ppu)/Rpul  
    % por lo tanto la variable de estado va a ser:
%     Vpu_p=Qpul
%     Vlv_p=0-Qav
%     Qav_p=(Plv-Pao-Qav*Rav)/Lav;
%     Vao_p=Qav-Qsys
%     Vvc_p=Qsys
%     Vrv_p=0-Qpv
%     Qpv_p=(Prv-Ppa-Qpv.Rpv)/Lpv
%     Vpa_p=Qpv-Qpul
%     x=[Vpu, Vlv, Qav, Vao, Vvc, Vrv, Qpv, Vpa]

             %calculo las presiones
            [Prv, Plv, ~]=PresionesVentr(t,x(6),x(2));
%             Pvc=Ees_vc*e(t)*x(5);
%             Ppa=Ees_pa*e(t)*x(8)+Pth;
%             Ppu=Ees_pu*e(t)*x(1)+Pth;
%             Pao=Ees_ao*e(t)*x(4);
            %VC, PA, PU y AO son camaras elasticas con comportamiento pasivo 
            Pvc=Ees_vc*x(5);
            Ppa=Ees_pa*x(8)+Pth(t);
            Ppu=Ees_pu*x(1)+Pth(t);
            Pao=Ees_ao*x(4);
%-------------------------------------------------------------------------
            %calculo las derivadas
            Qsys=(Pao-Pvc)/Rsys; %flujo a través de capilares sistémicos
            if Qsys<=0
                Qsys=0;
            end
            Qpul=(Ppa-Ppu)/Rpul; %flujo a través de capilares pulmonares
            if Qpul<=0
                Qpul=0;
            end
            
            Vpu_p=Qpul;
            Vlv_p=0-x(3);
            Qav_p=(Plv-Pao-x(3)*Rav)/Lav;
            Vao_p=x(3)-Qsys;
            Vvc_p=Qsys;
            Vrv_p=-x(7);
            Qpv_p=(Prv-Ppa-x(7)*Rpv)/Lpv;
            Vpa_p=x(7)-Qpul;

            dxdt=[Vpu_p; Vlv_p; Qav_p; Vao_p; Vvc_p; Vrv_p; Qpv_p; Vpa_p]; 
        end

        function [cierre stop direction]=event5(t,x) %condicion4 (en que momento termina la eyeccion del VI)
%             [Prv, ~, ~]=PresionesVentr(t,x(6),x(2));
            cierre= x(3); %busca en que momento Qav (flujo en la valv aortica) se hace nulo
            stop= 1; %corta cuando sucede el evento, o sea cuando Qpv se anula
            direction=-1; %agarro solo si se anula pero tiene derivada negativa
        end
    %% 6) ARRANCA ISOVOL DIAST DEL VI
    % MIENTRAS EL VD SIGUE EN EYECCION
function dxdt=sistema6(t,x) %x=[Vpu, Vlv, Vao, Vvc, Vrv, Qpv, Vpa]
    %por lo tanto: Qmt=Qtc=Qav=0 (no hay flujos de entrada a los ventrículos ni de salida del VI)
%     Qsys=(Pao-Pvc)/Rsys
%     Qpul=(Ppa-Ppu)/Rpul  
    % por lo tanto la variable de estado va a ser:
%     Vpu_p=Qpul
%     Vlv_p=0
%     Vao_p=-Qsys
%     Vvc_p=Qsys
%     Vrv_p=-Qpv
%     Qpv_p=(Prv-Ppa-Qav*Rpv)/Lpv;
%     Vpa_p=Qpv-Qpul
%     x=[Vpu, Vlv, Vao, Vvc, Vrv, Qpv, Vpa]
            
            %calculo las presiones
            [Prv, Plv, ~]=PresionesVentr(t,x(5),x(2));
%             Pvc=Ees_vc*e(t)*x(4);
%             Ppa=Ees_pa*e(t)*x(7)+Pth;
%             Ppu=Ees_pu*e(t)*x(1)+Pth;
%             Pao=Ees_ao*e(t)*x(3);
             %VC, PA, PU y AO son camaras elasticas con comportamiento pasivo 
            Pvc=Ees_vc*x(4);
            Ppa=Ees_pa*x(7)+Pth(t);
            Ppu=Ees_pu*x(1)+Pth(t);
            Pao=Ees_ao*x(3);
%-------------------------------------------------------------------------
            %calculo las derivadas
            Qsys=(Pao-Pvc)/Rsys; %flujo a través de capilares sistémicos
            if Qsys<=0
                Qsys=0;
            end
            Qpul=(Ppa-Ppu)/Rpul; %flujo a través de capilares pulmonares
            if Qpul<=0
                Qpul=0;
            end
            
            Vpu_p=Qpul;
            Vlv_p=0;
            Vao_p=-Qsys;
            Vvc_p=Qsys;
            Vrv_p=-x(6);
            Qpv_p=(Prv-Ppa-x(6)*Rpv)/Lpv;
            Vpa_p=x(6)-Qpul;

            dxdt=[Vpu_p; Vlv_p; Vao_p; Vvc_p; Vrv_p; Qpv_p; Vpa_p]; 
        end

        function [cierre stop direction]=event6(t,x) %condicion4 (en que momento termina la eyeccion del VD)
%             [Prv, ~, ~]=PresionesVentr(t,x(1),x(3));
            cierre= x(6); %busca en que momento Qpv (flujo en la valv pulmonar) se hace nulo
            stop= 1; %corta cuando sucede el evento
            direction=-1; %agarro solo si se anula pero tiene derivada negativa
        end
    %% 7) ARRANCA ISOVOL DIAST DEL VI
    % MIENTRAS EL VD SIGUE EN isovol sist
function dxdt=sistema7(t,x) %x=[Vpu, Vlv, Vao, Vvc, Vrv, Vpa]
%por lo tanto: Qmt=Qtc=Qpv=Qav=0 (no hay flujos de entrada ni de salida de los ventrículos)
%     Qsys=(Pao-Pvc)/Rsys
%     Qpul=(Ppa-Ppu)/Rpul  
    % por lo tanto la variable de estado va a ser:
%     Vpu_p=Qpul
%     Vlv_p=0
%     Vao_p=-Qsys
%     Vvc_p=Qsys
%     Vrv_p=0
%     Vpa_p=0-Qpul
%     x=[Vpu, Vlv, Vao, Vvc, Vrv, Vpa]
            
           %calculo las presiones
            [Prv, Plv, ~]=PresionesVentr(t,x(5),x(2));
%             Pvc=Ees_vc*e(t)*x(4);
%             Ppa=Ees_pa*e(t)*x(6)+Pth;
%             Ppu=Ees_pu*e(t)*x(1)+Pth;
%             Pao=Ees_ao*e(t)*x(3);
            %VC, PA, PU y AO son camaras elasticas con comportamiento pasivo 
            Pvc=Ees_vc*x(4);
            Ppa=Ees_pa*x(6)+Pth(t);
            Ppu=Ees_pu*x(1)+Pth(t);
            Pao=Ees_ao*x(3);
%-------------------------------------------------------------------------
            %calculo las derivadas
            Qsys=(Pao-Pvc)/Rsys; %flujo a través de capilares sistémicos
            if Qsys<=0
                Qsys=0;
            end
            Qpul=(Ppa-Ppu)/Rpul; %flujo a través de capilares pulmonares
            if Qpul<=0
                Qpul=0;
            end
            
            Vpu_p=Qpul;
            Vlv_p=0;
            Vao_p=-Qsys;
            Vvc_p=Qsys;
            Vrv_p=0;
            Vpa_p=0-Qpul;

            dxdt=[Vpu_p; Vlv_p; Vao_p; Vvc_p; Vrv_p; Vpa_p]; 
        end

        function [cierre stop direction]=event7(t,x) %condicion4 (en que momento termina la isovol diast del VD)
            [Prv, ~, ~]=PresionesVentr(t,x(5),x(2));
            %VC, PA, PU y AO son camaras elasticas con comportamiento pasivo 
%             Pvc=Ees_vc*e(t)*x(4);
            Pvc=Ees_vc*x(4);
            cierre= Pvc-Prv; %busca en que momento Pvc>Prv
            stop= 1; %corta cuando sucede el evento
            direction=1; %agarro solo si se anula pero tiene derivada positiva
        end
%% 8) ARRANCA LLENADO DEL VD 
    % MIENTRAS EL VI SIGUE EN isovol diast
function dxdt=sistema8(t,x) %x=[Vpu, Vlv, Vao, Vvc, Qtc, Vrv, Vpa]
%por lo tanto: Qmt=Qpv=Qav=0 (no hay flujo de salida de los ventrículos ni de entrada al VI)
%     Qsys=(Pao-Pvc)/Rsys
%     Qpul=(Ppa-Ppu)/Rpul  
    % por lo tanto la variable de estado va a ser:
%     Vpu_p=Qpul
%     Vlv_p=0
%     Vao_p=-Qsys
%     Vvc_p=Qsys-Qtc
%     Qtc_p=(Pvc-Prv-Qtc*Rtc)/Ltc
%     Vrv_p=Qtc-0
%     Vpa_p=0-Qpul
%     x=[Vpu, Vlv, Vao, Vvc, Qtc, Vrv, Vpa]
            
            %calculo las presiones
            [Prv, Plv, ~]=PresionesVentr(t,x(6),x(2));
%             Pvc=Ees_vc*e(t)*x(4);
%             Ppa=Ees_pa*e(t)*x(7)+Pth;
%             Ppu=Ees_pu*e(t)*x(1)+Pth;
%             Pao=Ees_ao*e(t)*x(3);
%VC, PA, PU y AO son camaras elasticas con comportamiento pasivo 
            Pvc=Ees_vc*x(4);
            Ppa=Ees_pa*x(7)+Pth(t);
            Ppu=Ees_pu*x(1)+Pth(t);
            Pao=Ees_ao*x(3);
%-------------------------------------------------------------------------
            %calculo las derivadas
            Qsys=(Pao-Pvc)/Rsys; %flujo a través de capilares sistémicos
            if Qsys<=0
                Qsys=0;
            end
            Qpul=(Ppa-Ppu)/Rpul; %flujo a través de capilares pulmonares
            if Qpul<=0
                Qpul=0;
            end
            
            Vpu_p=Qpul;
            Vlv_p=0;
            Vao_p=-Qsys;
            Vvc_p=Qsys-x(5);
            Qtc_p=(Pvc-Prv-x(5)*Rtc)/Ltc;
            Vrv_p=x(5)-0;
            Vpa_p=0-Qpul;

            dxdt=[Vpu_p; Vlv_p; Vao_p; Vvc_p; Qtc_p; Vrv_p; Vpa_p]; 
        end

        function [cierre stop direction]=event8(t,x) %condicion4 (en que momento termina isovol diast del VI e inicia el llenado del VI)
            [~, Plv, ~]=PresionesVentr(t,x(6),x(2));
            %VC, PA, PU y AO son camaras elasticas con comportamiento pasivo 
%             Ppu=Ees_pu*e(t)*x(1)+Pth;
            Ppu=Ees_pu*x(1)+Pth(t);
            cierre= Ppu-Plv; %busca en que momento Ppu>Plv
            stop= 1; %corta cuando sucede el evento
            direction=1; %agarro solo si se anula pero tiene derivada positiva
        end
    %% calculo las presiones ventriculares
    function [Prv, Plv, Vspt]=PresionesVentr(t,VRV,VLV)
                %1) calculo el volumen del septum
                %voy a calcular el Vseptum en mL por una cuestion de
                %resolucion del solver; si le paso los volumenes en m^3 (me van a quedar valores *10^(-6)==> va a cortar mas rapido el solver, xq su resolucion es ~10^(-16))
                f = @(Vs) VolumenSeptum(Vs,t,VRV*10^6,VLV*10^6); % function of dummy variable y
                options = optimset('TolX',sqrt(eps),'TolFun',sqrt(eps)); %le seteo las tolerancias para que corte con un cero aproximado de ~10^(-16)-->eso es eps
                Vspt = fsolve(f,0,options); %el resultado me lo da en mL
%                 %quiero ver que me da el resultado de evaluar el supuesto
%                 %cero en la funcion
%                 pepe=Funcionpepe(Vspt,t,VRV*10^6,VLV*10^6);
                Vspt=Vspt*10^(-6); %paso el resultado a m^3
                %2) Calculo Vrvf y Vlvf (volumenes free wall)

                Vrvf= VRV + Vspt;
                Vlvf= VLV - Vspt;

                %3) Calculo las presiones ES y ED para los dos ventriculos free wall

                Pes_rvf=Ees_rvf*(Vrvf-Vd_rvf);
                Ped_rvf=P0_rvf*(exp(lambda_rvf*(Vrvf-V0_rvf))-1);

                Pes_lvf=Ees_lvf*(Vlvf-Vd_lvf);
                Ped_lvf=P0_lvf*(exp(lambda_lvf*(Vlvf-V0_lvf))-1);

                %4) Calculo las presiones de los ventriculos free wall

                Prvf=e(t)*Pes_rvf+(1-e(t))*Ped_rvf;
                Plvf=e(t)*Pes_lvf+(1-e(t))*Ped_lvf;

                %5) Calculo el volumen del pericardio

                Vpcd= VRV + VLV;

                %6) Calculo la presion en la pared del pericardio
                Ppcd= P0_pcd*(exp(lambda_pcd*(Vpcd-V0_pcd))-1);

                %7) Calculo la presion en el pericardio

                Pperi= Ppcd + Pth(t); %Pth= presion en el torax

                %8) Calculo la presion en los ventriculos

                Prv= Prvf + Pperi;
                Plv= Plvf + Pperi;

        end
    %% para calcular numericamente el Vspt
        function F= VolumenSeptum(Vs,t,VRV,VLV) %calculo el volumen del septum
            %voy a usar todos los volumenes en mL para que calcule mejor,
            %si lo hago en m^3, los valores son del orden de 10^(-6)==> va
            %a cortar antes la cuenta por un tema de resolucion==> ajusto
            %lambdas y elastancias para que todas las unidades me den bien
%                     x(1)=Vrv y x(3)= Vlv
                    aux1=e(t)*(Ees_spt*10^(-6))*(Vs-(Vd_spt*10^6));
                    aux2=(1-e(t))*P0_spt*(exp((lambda_spt*10^(-6))*(Vs-(V0_spt*10^6)))-1);
                    aux3=e(t)*(Ees_lvf*10^(-6))*(VLV-Vs-Vd_lvf);
                    aux4=(1-e(t))*P0_lvf*(exp((lambda_lvf*10^(-6))*(VLV-Vs-(V0_lvf*10^6)))-1);
                    aux5=e(t)*(Ees_rvf*10^(-6))*(VRV+Vs-(Vd_rvf*10^6));
                    aux6=(1-e(t))*P0_rvf*(exp((lambda_rvf*10^(-6))*(VRV+Vs-(V0_rvf*10^6)))-1);
%                     aux5=aux51+aux52;
                    F= aux1+aux2-aux3-aux4+aux5+aux6; 
% %                     F=e(t)*Ees_spt*(Vs-Vd_spt)+(1-e(t))*P0_spt*(exp(lambda_spt*(Vs-V0_spt))-1)-e(t)*Ees_lvf(x(3)-Vs-Vd_lvf)+(1-e(t))*P0_lvf*(exp(lambda_lvf*(x(3)-Vs-V0_lvf))-1)-(e(t)*Ees_rvf*(x(1)-Vs-Vd_rvf)+(1-e(t))*P0_rvf*(exp(lambda_rvf*(x(1)-Vs-V0_rvf))-1));
        end
    
%         function Resulpepe= Funcionpepe(Vs,t,VRV,VLV) %calculo el volumen del septum
% %                     x(1)=Vrv y x(3)= Vlv
%                     aux1=e(t)*Ees_spt*10^(-6)*(Vs-Vd_spt*10^6);
%                     aux2=(1-e(t))*P0_spt*(exp(lambda_spt*10^(-6)*(Vs-V0_spt*10^6))-1);
%                     aux3=e(t)*Ees_lvf*10^(-6)*(VLV-Vs-Vd_lvf);
%                     aux4=(1-e(t))*P0_lvf*(exp(lambda_lvf*10^(-6)*(VLV-Vs-V0_lvf*10^6))-1);
%                     aux51=e(t)*Ees_rvf*10^(-6)*(VRV-Vs-Vd_rvf*10^6);
%                     aux52=(1-e(t))*P0_rvf*(exp(lambda_rvf*10^(-6)*(VRV-Vs-V0_rvf*10^6))-1);
%                     aux5=aux51+aux52;
%                     Resulpepe= aux1+aux2-aux3+aux4-aux5; 
% % %                     F=e(t)*Ees_spt*(Vs-Vd_spt)+(1-e(t))*P0_spt*(exp(lambda_spt*(Vs-V0_spt))-1)-e(t)*Ees_lvf(x(3)-Vs-Vd_lvf)+(1-e(t))*P0_lvf*(exp(lambda_lvf*(x(3)-Vs-V0_lvf))-1)-(e(t)*Ees_rvf*(x(1)-Vs-Vd_rvf)+(1-e(t))*P0_rvf*(exp(lambda_rvf*(x(1)-Vs-V0_rvf))-1));
%         end

end