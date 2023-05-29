function [VF,p] = pkmethod(mu,e,a,r,sigma,V_Vec,tol,k_max)
%PKMETHOD Encontra a velocidade de flutter e autovalores usando o método pk

j=0;
i=sqrt(-1);
for V=V_Vec
    j = j+1;
    k= 1;
    iter = 1;
    erro = 1;
    while erro >=tol && iter < k_max
        k_old = k;
        Ck = (0.01365 +0.2808*i*k - k^2/2)/(0.01365 + 0.3455*i*k - k^2); 
        p11 = [1 0 sigma^2/V^2- k^2/mu+(2*i*k*Ck)/mu];
        p12 = [(e-a) 0 (k*(i+a*k)+(2+i*k*(1-2*a))*Ck)/mu];
        p21 = [(e-a) 0 (a*k^2-i*k*(1+2*a)*Ck)/mu];
        p22 = [r^2 0 (8*mu*r^2/V^2+4*i*((1+2*a)*(2*i-k*(1-2*a))*Ck)-k*(k-4*i+8*a*(i+a*k)))/(8*mu)];
        DFlutter = conv(p11,p22) - conv(p12,p21);
        p(:,j) = roots(DFlutter);
        k = min(abs(imag(p(:,j))));
        iter = iter + 1;
        erro = abs(k-k_old)/k;
    end
    k_vec(j) = k;
    p(:,j) = V*p(:,j);
end
k = imag(p);
re = p(k>=0);
clear p
for i = 1:size(k,2)
    if round(imag(re(2*i-1)),4) == round(max(k(:,i)),4)
    p(1,i)=re(2*i-1);
    p(2,i)=re(2*i);
    else
    p(1,i) = re(2*i);
    p(2,i) = re(2*i-1);
    end
end

[~,i] = max(max(sign(real(p))));
gamma_max = max(real(p));
if i == 1
    VF = NaN;
else
VF = V_Vec(i-1) + (gamma_max(i)-gamma_max(i-1))/(V_Vec(i)-V_Vec(i-1)) *(-gamma_max(i-1));
end

end

