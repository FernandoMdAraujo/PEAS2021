%Fernando Marques de Araujo - NUP 6801112

%Creditos pelo Codigo:
%Este script foi feito com base na Aula #10 da Professora Maria D. Miranda
%Aula ministrada em 21/set/2021
%Curso PTC3451 – Processamento Estatístico e Adaptativo de Sinais

%Teste #3

var_s=0.25; %s(n) eh um ruido gaussiano de media nula
omega=2*pi/40;
A=5; %Amplitude
theta=pi/6;
M=2;% numeo de coeficientes

N=4500;
n=0:N-1;

phi_v=2*pi*rand;
phi_u=2*pi*rand;

s=sqrt(var_s)*randn(N,1);
x=sin(omega*n+theta+phi_v)';
u=A*sin(omega*n+phi_v)';

d=s+x;

var_d=cov(d)

r=xcorr(u,M-1,'unbiased');
ru=r(M:end);
R=toeplitz(ru)
pdu=xcorr(d,u,M-1,'unbiased');
p=pdu(M:end)

wotimo=R\p
Jmin=var_d-p'*wotimo

yotimo=filter(wotimo,1,u);
eotimo=d-yotimo;
figure(1)
N2=200;
subplot(311)
plot([0:N2-1],u(1:N2),'g'); hold on;
plot([0:N2-1],d(1:N2),'k'); hold off;
grid; legend('u(n)','d(n)');
axis([0 N2-1 -5.1 5.1]);
title('Visualização dos sinais de entrada correlacionados u(n) e d(n)');
subplot(312)
plot([0:N2-1],x(1:N2),'g'); hold on;
plot([0:N2-1],d(1:N2),'k');
plot([0:N2-1],yotimo(1:N2),'m.'); hold off;
grid; legend('x(n)','d(n)','yotimo(n)');
axis([0 N2-1 -2.2 2.2]);
title('Visualização da correlação dos sinais x(n) e yotimo');
subplot(313)
plot([0:N2-1],d(1:N2),'k'); hold on;
plot([0:N2-1],s(1:N2),'g'); plot([0:N2-1],eotimo(1:N2),'m.'); hold off;
grid; legend('d(n)','s(n)','eotimo(n)');
xlabel('n'); axis([0 N2-1 -2.2 2.2]);
title('Visualização de como o erro ótimo "persegue" o sinal desejado s(n)');




%B
N1=1500; mu=0.01;
wSD=zeros(M,N1); ySD=zeros(N1,1); eSD=zeros(N1,1);
JSD=zeros(N1,1); uv=zeros(M,1);
for i=1:N1
wSD(:,i+1)=wSD(:,i)+mu*(p-R*wSD(:,i));
uv=[u(i); uv(1:M-1)];
ySD(i)=uv'*wSD(:,i);
eSD(i)=d(i)-ySD(i);
JSD(i)=var_d-2*wSD(:,i)'*p+wSD(:,i)'*R*wSD(:,i);
end
Jmin=var_d-2*wotimo'*p+wotimo'*R*wotimo;
wotimoN1=kron(wotimo,ones(1,N1))';
figure(2)
subplot(311)
plot(wotimoN1(:,1),'m'); hold on;
plot(wSD(1,1:N1)','b');
plot(wotimoN1(1:N1,2),'m');
plot(wSD(2,1:N1)','b'); hold off;
grid; legend('wotimo','wSD(n)');
title('Coeficientes de Steepest Descent vs Coeficientes de Wiener');
axis([0 N1-1 -1.2 1.2]);
subplot(312)
plot([0:N1-1],eSD(1:N1),'b'); hold on;
plot([0:N1-1],s(1:N1),'g'); hold off;
legend('eSD(n)','s(n)'); grid;
axis([0 N1-1 -2.5 2.5]);
title('sinais erro, eSD(n), e sinal desejado, s(n)');
subplot(313)
plot(JSD,'b'); hold on;
plot(Jmin*ones(N1,1),'m'); hold off;
grid; legend('JSD(n)','Jmin');
axis([0 N1-1 0 2]); xlabel('n');
title('Erro Quadrático Médio (SD) vs Erro Quadrático Médio Mínimo (Wiener)');


%C
N1=1500; n=0:N1-1; mu=0.01;
EQM=zeros(N1,1); EQME=zeros(N1,1);
EQMdB=zeros(N1,3);
EQMEdB=zeros(N1,3);
EQMEteodB=zeros(3,N1);
K=500; % Número de realizações
for j=3:-1:1
    for k=1:K
    phiV=2*pi*rand(1);
    s=sqrt(var_s)*randn(N1,1);
    x=sin(omega*n+theta+phiV)';
    u=A*sin(omega*n+phiV)';
    d=s+x;
    %[yLMS,eLMS,wLMS,eaLMS]=LMS(M,mu,u,d,wo);
    N=length(u);
    wLMS=zeros(M,N); yLMS=zeros(N,1);
    eLMS=zeros(N,1); uv=zeros(M,1);
    eaLMS=zeros(N,1);
        for i=1:N
        uv=[u(i); uv(1:M-1)];
        yLMS(i)=uv'*wLMS(:,i);
        eLMS(i)=d(i)-yLMS(i);
        wLMS(:,i+1)=wLMS(:,i)+j*mu*eLMS(i)*uv;
        eaLMS(i) = uv'*(wotimo - wLMS(:,i));
        end
    EQM=EQM+eLMS.*eLMS;
    EQME=EQME+eaLMS(1:N).*eaLMS(1:N);
    end
    EQM=EQM/K; EQME=EQME/K;
    EQMdB(:,j)=10*log10(EQM); EQMEdB(:,j)=10*log10(EQME);


EQMEteo(j,:)=Jmin*((trace(R)*j*mu)/(2-(trace(R)*j*mu)));
EQMEteodB(j,:)=10*log10(EQMEteo(j,:));
end
figure(4)
subplot(211); plot(n,EQM,'b'); hold on;
plot(Jmin*ones(N1,1),'m'); hold off;
grid; legend('EQM(n)','Jmin');
axis([0 N1-1 0 2]); xlabel('n'); title('EQM');
subplot(212); plot(n,EQME,'b'); hold on;
plot(n,EQMEteo(1,:)*ones(N1,1),'m'); hold off;
grid; legend('EQME(n)','EQMEmin');
axis([0 N1-1 0 2]); xlabel('n'); title('EQME');

figure(5)
subplot(211)
plot(n,EQMdB(:,1),'r'); hold on
plot(n,EQMdB(:,2),'b'); plot(n,EQMdB(:,3),'g');
plot(n,10*log10(Jmin)*ones(N1,1),'k--'); hold off
grid; legend('mu=0.01','mu=0.02','mu=0.03','Jmin');
ylim([-6.5 -1.5]); ylabel('EQM dB')
title('EQM [dB] para \mu=0,01,  \mu=0,02, e \mu=0,03 ');
subplot(212)
plot(n,EQMEdB(:,1),'r'); hold on
plot(n,EQMEteodB(1,:),'r--'); plot(n,EQMEdB(:,2),'b');
plot(n,EQMEteodB(2,:),'b--'); plot(n,EQMEdB(:,3),'g');
plot(n,EQMEteodB(3,:),'g--'); hold off
grid;
legend('mu=0.01','mu=0.01 teo','mu=0.02','mu=0.02 teo','mu=0.03','mu=0.03 teo');
ylim([-15 -4]);ylabel('EQME dB')
title('EQME [dB], simulado vs teórico, para \mu=0,01,  \mu=0,02, e \mu=0,03 ');


%C1
figure(7)
subplot(2,2,1:2)
%plot([0:N1-1],eSD(1:N1),'b'); hold on;
plot([0:N1-1],s(1:N1),'g');  hold on;
plot([0:N1-1],eLMS(1:N1),'r');
hold off;
%legend('eSD(n)','s(n)', 'eLMS(n)'); grid;
legend('s(n)', 'eLMS(n)'); grid;
axis([0 N1-1 -2.5 2.5]);
title('sinal de erro, eLMS(n), e sinal desejado, s(n)');
subplot(2,2,3)
%plot([0:200],eSD(1:201),'b'); hold on;
plot([0:200],s(1:201),'g');  hold on;
plot([0:200],eLMS(1:201),'r');
hold off;
%legend('eSD(n)','s(n)', 'eLMS(n)'); grid;
legend('s(n)', 'eLMS(n)'); grid;
axis([0 100 -2.5 2.5]);
title('Zoom do início');
subplot(2,2,4)
%plot([1400:1499],eSD(1401:1500),'b'); hold on;
plot([1400:1499],s(1401:1500),'g');  hold on;
plot([1400:1499],eLMS(1401:1500),'r');
hold off;
legend('eSD(n)','s(n)', 'eLMS(n)'); grid;
legend('s(n)', 'eLMS(n)'); grid;
axis([1400 1499 -2.5 2.5]);
title('Zoom do final');
figure(6)
plot(wotimoN1(:,1),'m'); hold on;
plot(wSD(1,1:N1)','b');
plot(wLMS(1,1:N1)','r');
plot(wotimoN1(1:N1,2),'m');
plot(wSD(2,1:N1)','b'); 
plot(wLMS(2,1:N1)','r');hold off;
grid; legend('wotimo','wLMS(n)' );
title('Coeficientes de Wiener vs Coeficientes do LMS');
axis([0 N1-1 -1.2 1.2]);

