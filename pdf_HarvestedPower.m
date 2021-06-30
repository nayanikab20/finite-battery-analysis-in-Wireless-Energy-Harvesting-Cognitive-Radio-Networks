clear all;
s = 100000;
M=3;
PU_tx = db2pow(0);
eta = 0.8;
alpha= 0.5;%3;
T = 1;%(1/1)*10^-9;
threshold1 = 10;   % for battery for SS

c = threshold1; %capacity
x = 0:0.05:c;
d1=1; d2=sqrt(5);  d4=sqrt(5); d6=2;
n= -3;
v1 = d1^n;

m=3;
b=v1;

l =  eta*alpha*PU_tx*T;

%simulation
f11 = nkg_sq2(b,m,x,s);nkg_sq2(b,m,x,s); f12 = nkg_sq2(b,m,x,s); f13 = nkg_sq2(b,m,x,s); 

for i = 1:size(f11)
            
            Ehs(i) = eta*alpha*T*PU_tx*( f11(i) + f12(i) + f13(i) );                     
            %var(i) = f11(i) + f12(i) + f13(i);            
            phs(i) = min(Ehs(i), threshold1);
end


% analytical

for i= 1:length(0:0.05:(threshold1-0.05))
    
    p(i) = (x(i))^(3*m-1);
    q(i) = exp(-(x(i))*(m/(l*b)));
    r(i) = 1/gamma(3*m);
    t(i) =(m/(l*b))^3*m;
    
    phs_a(i) = p(i)*q(i)*r(i)*t(i);
end
phs_a(length(0:0.05:threshold1)) = gammainc((threshold1*m/(l*b)), 3*m, 'upper');

figure, 
histogram(phs,'Normalization', 'pdf')
hold on;
plot( 0:0.05:threshold1, phs_a, 'LineWidth', 2)
legend('simulated', 'analytical')
xlabel('Ehs')
ylabel('PDF')

for i= 1:length(0:0.05:(threshold1))
        
    cdf_a(i) = gammainc((x(i))*m/(l*b), 3*m);
end

for i= 1:length(threshold1:0.05:c)
        
    cdf_a( length(0:0.05:(threshold1))+i ) = 1;
end


figure, histogram(phs, 'Normalization', 'cdf')
hold on;
plot(0:0.05:c,cdf_a, 'LineWidth', 2)
xlim([1 c])
ylim([0 1.5])
legend('simulated', 'analytical')
xlabel('Ehs')
ylabel('CDF')