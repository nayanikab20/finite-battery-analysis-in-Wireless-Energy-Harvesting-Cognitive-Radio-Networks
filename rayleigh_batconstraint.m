% simulation of outage probability for a cognitive relay network with
% energy harvesting with battery constraint

% number of transmitters = 3
% range of PI =-20:1:20
% range of gamma_th =-10:5:0
% RV with rayleign distribution are f(1,j), f(2,j), f(3,j), h1, h2
% number of samples 10,000
% assumptions - 1. f(1,j) same for all j; similarly for f(2,j) and f(3,j)
%               2. alpha =0.5
%               3. eta =0.8

clear all;
%close all;

s = 100000;
M=3;
PU_tx = db2pow(0);
eta = 0.8;
alpha= 0.5;
T = (1/1)*10^-9;
threshold1 = db2pow(25);    % for battery for SS in watts
threshold2 = db2pow(25);   % for battery of SR in watts

d1=1; d2=sqrt(2); d3=sqrt(5); d4=sqrt(5); d5=sqrt(2); d6=1; d7=1;
m= -4;
v1 = d1^m; v2 = d2^m; v3 = d3^m;
w1 = d4^m; w2 = d5^m;
y1 = d6^m; y2 = d7^m;

% rv and their pdf
h1 = exprnd(y1,s,1);
h2 = exprnd(y2,s,1);
f11 = exprnd(v1,s,1); f12 = exprnd(v1,s,1); f13 = exprnd(v1,s,1);
f21 = exprnd(v2,s,1); f22 = exprnd(v2,s,1); f23 = exprnd(v2,s,1);
f31 = exprnd(v3,s,1); f32 = exprnd(v3,s,1); f33 = exprnd(v3,s,1);
g11 = exprnd(w1,s,1); g12 = exprnd(w1,s,1); g13 = exprnd(w1,s,1);
g21 = exprnd(w2,s,1); g22 = exprnd(w2,s,1); g23 = exprnd(w2,s,1);

count = 0;
P_out =zeros(3, 41);

for th = 1:3
     gamma_th = db2pow( -10+(th-1)*5) 
    for interference = 1:41
        PI = db2pow(-20 + (interference -1)*1);
        for i = 1:s
            
            Ehs(i) = eta*alpha*T*PU_tx*( f11(i) + f12(i) + f13(i) );
            Ehr(i) = eta*alpha*T*PU_tx*( f21(i) + f22(i) + f23(i) );
            
            Phs(i) = 2*Ehs(i)/(( 1 - alpha )*T);
            g1 =[g11(i); g12(i); g13(i)];
            PIs = PI/max(g1);
            p = 2*eta*alpha/(1-alpha);
            
            if(Phs(i)<min(PIs, threshold1))
                Phs(i)= Phs(i);
                gammaR(i,interference,th) = p*h1(i)*( f11(i) + f12(i) + f13(i) )/( f21(i) +f22(i) +f23(i) );
                
            elseif(threshold1<min(Phs(i), PIs))
                Phs(i) = threshold1;
                gammaR(i,interference,th) = threshold1*h1(i)/(PU_tx*( f21(i) +f22(i) +f23(i) ));
                
            elseif(PIs<min(threshold1, Phs(i)))
                
                 gammaR(i,interference,th) = PI*h1(i)/(max(g1)* PU_tx*( f21(i) +f22(i) +f23(i) ) );
            end
            
            
            Phr(i) = 2*Ehr(i)/(( 1 - alpha )*T);
            g2 =[g21(i); g22(i); g23(i)];
            PIr = PI/max(g2);
 
            if(Phr(i)<min(PIr, threshold2))
                Phr(i)= Phr(i);
                 gammaD(i,interference,th) = p*h2(i)*( f21(i) + f22(i) + f23(i) )/( f31(i) +f32(i) +f33(i));
                 
            elseif(threshold2<min(Phr(i), PIr))
                Phr(i) = threshold2;
                gammaD(i,interference,th) = threshold2*h2(i)/(PU_tx*( f31(i) +f32(i) +f33(i) ) );
                
            elseif(PIr<min(threshold2, Phr(i)))
                  gammaD(i,interference,th) = PI*h2(i)/(max(g2)*PU_tx*( f31(i) +f32(i) +f33(i) ) );
                  
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
title('Rayleigh Fading Relay network- Battery constraint condition')

% figure, histogram(Ehs)
% title('Rayleigh distributed channels')
% 
% figure, histogram(Ehr)
% title('Rayleigh distributed channels')
% 
% figure, histogram(Phs)
% title('Rayleigh distributed channels')
% 
% figure, histogram(Phr)
% title('Rayleigh distributed channels')
