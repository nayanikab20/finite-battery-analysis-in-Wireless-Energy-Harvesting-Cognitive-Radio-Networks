 % simulation of outage probability for a cognitive relay network with
% energy harvesting with battery constraint

% number of transmitters = 3
% range of PI =-20:1:20
% range of gamma_th =-10:5:0
% RV with nakagami distribution are f(1,j), f(2,j), f(3,j), h1, h2
% number of samples 10,000
% assumptions - 1. f(1,j) same for all j; similarly for f(2,j) and f(3,j)
%               2. alpha =0.5
%               3. eta =0.8

clear all;
%close all;

s = 100000;
M=3;
PU_tx = db2pow(0);
PI = db2pow(10);
eta = 0.8;
alpha= 0.5;
T =1;% (1/1)*10^-9;
%Bmax = 2.323;   % for battery for SS


x = 0:0.05:10;
d1=1; d2=sqrt(5);  d4=sqrt(5); d6=2;
n= -3;
v1 = d1^n; v2 = d2^n; 
w1 = d4^n; 
y1 = d6^n; 

m=3;

h1 = nkg_sq2(y1,m,x,s);
f11 = nkg_sq2(v1,m,x,s); f12 = nkg_sq2(v1,m,x,s); f13 = nkg_sq2(v1,m,x,s); 
f21 = nkg_sq2(v2,m,x,s); f22 = nkg_sq2(v2,m,x,s); f23 = nkg_sq2(v2,m,x,s);
g11 = nkg_sq2(w1,m,x,s); g12 = nkg_sq2(w1,m,x,s); g13 = nkg_sq2(w1,m,x,s);

count = 0;
P_out =zeros(3, 41);
step_size = 0.5;

size=size(f11);

for th = 1:3
     gamma_th = db2pow( -10+(th-1)*5) 
    for Bmax = 1:41
        Bmax_E = db2mag(0 + (Bmax -1)*step_size);
        for i = 1:size
            
            Ehs(i) = eta*alpha*T*PU_tx*( f11(i) + f12(i) + f13(i) );
                                  
%             Phs(i) = 2*Ehs(i)/(( 1 - alpha )*T);
            Phs(i) = Ehs(i)/(( 1 - alpha )*T);
            threshold = Bmax_E/(( 1 - alpha )*T);
            
            g1 =[g11(i); g12(i); g13(i)];
            PIs = PI/max(g1);
            
            phs(i) = min(Phs(i), threshold); %put battery value in terms of power
            
%             p = 2*eta*alpha/(1-alpha);
            
            p = eta*alpha/(1-alpha);
            
            if(Phs(i)<PIs & Phs(i)< threshold)
                gammaR(i,Bmax,th) =p*h1(i)*( f11(i) + f12(i) + f13(i) )/( f21(i) +f22(i) +f23(i) );
                
            elseif(threshold<Phs(i) & threshold< PIs)
                
                 gammaR(i,Bmax,th) =threshold*h1(i)/(PU_tx*( f21(i) +f22(i) +f23(i) ));
                 
            elseif(PIs<threshold & PIs<Phs(i))
                
                gammaR(i,Bmax,th) = PI*h1(i)/(max(g1)* PU_tx*( f21(i) +f22(i) +f23(i) ) );
            end
            
           
            if(gammaR(i,Bmax,th)>= gamma_th)
                count= count+1;
                
            end
            
        end
        
        P_out(th, Bmax) = (s -count)/s;
        count = 0;
    end
    
end
%%
%figure,        
semilogy(0:step_size:20, P_out(1,:), '-r')
hold on; 
%semilogy(0:step_size:20, P_out(2,:), '-g')
semilogy(0:step_size:20, P_out(3,:), '-b')
ylim([ 0.01 1])
xlabel('Bmax (dBW)')
ylabel('P out')
%legend('gamma th= -10', 'gamma th = -5', 'gamma th = 0')
%title('Nakagami-m non-relay network- Battery constraint condition PI=10dB')

%%
%figure, histogram(phs)
%hold on;
%l=2*eta*alpha*PU_tx/(1-alpha);
%ideal = l*nkg_sq2(3*v1,3*m,x,s);
%histogram(ideal)
%title('Nakagami-m distributed channels')

 