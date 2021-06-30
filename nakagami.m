% simulation of outage probability for a cognitive relay network with energy harvesting

% number of transmitters = 3
% range of PI =-20:1:20
% range of gamma_th =-10:5:0
% RV with rayleign distribution are f(1,j), f(2,j), f(3,j), h1, h2
% number of samples 10,000
% assumptions - 1. f(1,j) same for all j; similarly for f(2,j) and f(3,j)
%               2. alpha =0.5
%               3. eta =0.8

%clear all;
%close all;

s = 100000; % sample size
M=3;
PU_tx = db2pow(0);
eta = 0.8;
alpha= 0.5;
T = (1/1)*10^-9;

%range of values for channel coeff rv
x = 0:0.05:10 ;

d1=1; d2=sqrt(2); d3=sqrt(5); d4=sqrt(5); d5=sqrt(2); d6=1; d7=1;
n= -3;
v1 = d1^n; v2 = d2^n; v3 = d3^n;
w1 = d4^n; w2 = d5^n;
y1 = d6^n; y2 = d7^n;

m = 2;

% link gain realisation rv and their pdf
h1 = nkg_sq2(y1,m,x,s);
h2 = nkg_sq2(y2,m,x,s);
f11 = nkg_sq2(v1,m,x,s); f12 = nkg_sq2(v1,m,x,s); f13 = nkg_sq2(v1,m,x,s); 
f21 = nkg_sq2(v2,m,x,s); f22 = nkg_sq2(v2,m,x,s); f23 = nkg_sq2(v2,m,x,s);
f31 = nkg_sq2(v3,m,x,s); f32 = nkg_sq2(v3,m,x,s); f33 = nkg_sq2(v3,m,x,s);
g11 = nkg_sq2(w1,m,x,s); g12 = nkg_sq2(w1,m,x,s); g13 = nkg_sq2(w1,m,x,s);
g21 = nkg_sq2(w2,m,x,s); g22 = nkg_sq2(w2,m,x,s); g23 = nkg_sq2(w2,m,x,s);


size = length(f11);

count = 0;
P_out =zeros(3, 41);

for th = 1:3
     gamma_th = db2pow( -10+(th-1)*5) 
    for interference = 1:41
        PI = db2pow(-20 + (interference -1)*1);
        for i = 1:size
            
            Ehs(i) = eta*alpha*T*PU_tx*( f11(i) + f12(i) + f13(i) );
            Ehr(i) = eta*alpha*T*PU_tx*( f21(i) + f22(i) + f23(i) );
            
            Phs = 2*Ehs(i)/(( 1 - alpha )*T);
            g1 =[g11(i); g12(i); g13(i)];
            PIs = PI/max(g1);
            
            p = 2*eta*alpha/(1-alpha);
            
            gammaR1(i,interference,th) = p*h1(i)*( f11(i) + f12(i) + f13(i) )/( f21(i) +f22(i) +f23(i) );
            gammaR2(i,interference,th) = PI*h1(i)/(max(g1)* PU_tx*( f21(i) +f22(i) +f23(i) ) );
            
            if(Phs<PIs)
                gammaR(i,interference,th) =gammaR1(i,interference,th);
                
            else
                
                 gammaR(i,interference,th) =gammaR2(i,interference,th);
            end
            
            Phr = 2*Ehr(i)/(( 1 - alpha )*T);
            g2 =[g21(i); g22(i); g23(i)];
            PIr = PI/max(g2);
            
            
            gammaD1(i,interference,th)= p*h2(i)*( f21(i) + f22(i) + f23(i) )/( f31(i) +f32(i) +f33(i));
            gammaD2(i,interference,th) =  PI*h2(i)/(max(g2)*PU_tx*( f31(i) +f32(i) +f33(i) ) );
          
            
            if(Phr<PIr)
                 gammaD(i,interference,th) = gammaD1(i,interference,th);
                
            else
                  gammaD(i,interference,th) = gammaD2(i,interference,th);
                  
            end
            
            if(gammaR(i,interference,th)>= gamma_th && gammaD(i,interference,th) >= gamma_th)
                count= count+1;
            end
            
        end
        
        P_out(th, interference) = (s -count)/s;
        count = 0;
    end
    
end

figure,       
semilogy(-20:1:20, P_out(1,:), '.-r')
hold on; 
semilogy(-20:1:20, P_out(2,:), '.-g')
semilogy(-20:1:20, P_out(3,:), '.-b')
xlabel('P_I (dBW)')
ylabel('P out')
legend('gamma th= -10', 'gamma th = -5', 'gamma th = 0')
title('Nakagami-m distributed channels')


figure, histogram(Ehs)
title('Nakagami-m distributed channels')

figure, histogram(Ehr)
title('Nakagami-m distributed channels')



