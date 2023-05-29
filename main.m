% Student: Willian Leandro dos Santos Pinto

clear all
close all
clc

mu = 20;
e = -1/10;
a = -1/5;
r = sqrt(6/25);
sigma = 2/5;
V_Vec = linspace(1e-9,3,500);

[VF, p] = pkmethod(mu,e,a,r,sigma, V_Vec, 1e-4,1000);

figure(1)

hold on
plot(V_Vec, max(real(p)),'k', 'LineWidth',1.5)
plot(V_Vec, min(real(p)),'k','LineWidth',1.5)
hold off
grid on
grid minor
box on
set(gcf,'Color',[1 1 1])
xlabel('V')
ylabel('\Gamma/\omega_\theta')

figure(2)
hold on
plot(V_Vec, max(imag(p)),'k','LineWidth',1.5)
plot(V_Vec, min(imag(p)),'k','LineWidth',1.5)
hold off
grid on
grid minor
box on
set(gcf,'Color',[1 1 1])
xlabel('V')
ylabel('\Omega/\omega_\theta')
