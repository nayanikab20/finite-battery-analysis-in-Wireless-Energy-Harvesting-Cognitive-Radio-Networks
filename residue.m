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
eta = 0.8;
alpha= 0.5;
T =1;% (1/1)*10^-9;
Bmax = 2.323;   % for battery for SS
%threshold2 =     % for battery of SR

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

size=size(f11);

for th = 1:3
     gamma_th = db2pow( -10+(th-1)*5) 
    for interference = 1:41
        PI = db2pow(-20 + (interference -1)*1);
        for i = 1:size
            
            Ehs(i) = eta*alpha*T*PU_tx*( f11(i) + f12(i) + f13(i) );
                                  
%             Phs(i) = 2*Ehs(i)/(( 1 - alpha )*T);
            Phs(i) = Ehs(i)/(( 1 - alpha )*T);
            threshold1 = Bmax/(( 1 - alpha )*T);
            
            g1 =[g11(i); g12(i); g13(i)];
            PIs = PI/max(g1);
            
            phs(i) = min(Phs(i), threshold1); %put battery value in terms of power
            
%             p = 2*eta*alpha/(1-alpha);
            
            p = eta*alpha/(1-alpha);
            
            if(Phs(i)<PIs & Phs(i)< threshold1)
                gammaR(i,interference,th) =p*h1(i)*( f11(i) + f12(i) + f13(i) )/( f21(i) +f22(i) +f23(i) );
                
            elseif(threshold1<Phs(i) & threshold1< PIs)
                
                 gammaR(i,interference,th) =threshold1*h1(i)/(PU_tx*( f21(i) +f22(i) +f23(i) ));
                 
            elseif(PIs<threshold1 & PIs<Phs(i))
                
                gammaR(i,interference,th) = PI*h1(i)/(max(g1)* PU_tx*( f21(i) +f22(i) +f23(i) ) );
            end
            
           
            if(gammaR(i,interference,th)>= gamma_th)
                count= count+1;
                %residual_e(th, interference, i) = Phs(i) -min( min(threshold1, Phs(i)), PIs);%(threshold1 - min(Ehs(i), Bmax))*(( 1 - alpha )*T);
            end
            
            residual_e(th, interference, i) = phs(i) -min( phs(i), PIs);
            
        end
        
        P_out(th, interference) = (s -count)/s;
        count = 0;
    end
    
end
%%
figure,        
semilogy(-20:1:20, P_out(1,:), '*r')
hold on; 
semilogy(-20:1:20, P_out(2,:), '*g')
semilogy(-20:1:20, P_out(3,:), '*b')
ylim([ 0.01 1])
xlabel('P_I (dBW)')
ylabel('P out')
legend('gamma th= -10', 'gamma th = -5', 'gamma th = 0')
title('Nakagami-m non-relay network- Battery constraint condition')

%%
%figure, histogram(phs)
%hold on;
%l=2*eta*alpha*PU_tx/(1-alpha);
%ideal = l*nkg_sq2(3*v1,3*m,x,s);
%histogram(ideal)
%title('Nakagami-m distributed channels')

%%
for i=1:41
    
    figure, 
    histogram(residual_e(1,i, :))
    title(-20+i)
end
%%

%Bmax = threshold1;
I3=0;
I1=0;
I2=0;

for th = 1:3
     gamma_th = db2pow( -10+(th-1)*5) 
    for interference = 1:41
        PI = db2pow(-20 + (interference -1)*1);
        for i = 1:size
            
            Ehs(i) = eta*alpha*T*PU_tx*( f11(i) + f12(i) + f13(i) );
            V(i) = min(Bmax, Ehs(i))/((1-alpha)*T);
            P(i) = PI/max([g11(i); g12(i); g13(i)]);
            C(i) = gamma_th*PU_tx*(f21(i) + f22(i) + f23(i))/h1(i);
            
            if( V(i)<C(i) && P(i)<C(i) )
                I3 = I3+1;
            else
               I3=I3;
            end
             
            if(V(i)<C(i))
                I1 =I1+1;
            else
                I1=I1;
            end
            
            if(P(i)<C(i))
                I2 =I2+1;
            else
                I2=I2;
            end
                
        end
        I3_norm(th, interference) = I3/s;
        I1_norm(th, interference) = I1/s;
        I2_norm(th, interference) = I2/s;
        I3= 0;
        I1=0;
        I2=0;
    end
end

%%
figure, semilogy(-20:1:20, I3_norm(1,:), '*r')
hold on;
semilogy(-20:1:20, I3_norm(2,:), '*b')
semilogy(-20:1:20, I3_norm(3,:), '*g')
xlabel('P_I (dBW)')
ylabel('I3')
legend('gamma th= -10', 'gamma th = -5', 'gamma th = 0')
title('I3 simulated')
    
figure, semilogy(-20:1:20, I1_norm(1,:), '-r')
hold on;
semilogy(-20:1:20, I1_norm(2,:), '-b')
semilogy(-20:1:20, I1_norm(3,:), '-g')
xlabel('P_I (dBW)')
ylabel('I1')
legend('gamma th= -10', 'gamma th = -5', 'gamma th = 0')
title('I1 simulated')

figure, semilogy(-20:1:20, I2_norm(1,:), '.r')
hold on;
semilogy(-20:1:20, I2_norm(2,:), '.b')
semilogy(-20:1:20, I2_norm(3,:), '.g')
xlabel('P_I (dBW)')
ylabel('I2')
legend('gamma th= -10', 'gamma th = -5', 'gamma th = 0')
title('I2 simulated')
            
P_out_sim(1,:) = I1_norm(1,:)+I2_norm(1,:)-I3_norm(1,:);
P_out_sim(2,:) = I1_norm(2,:)+I2_norm(2,:)-I3_norm(2,:);
P_out_sim(3,:) = I1_norm(3,:)+I2_norm(3,:)-I3_norm(3,:);

figure, semilogy(-20:1:20, P_out_sim(1,:), 'r')
hold on;
semilogy(-20:1:20, P_out_sim(2,:), 'b')
semilogy(-20:1:20, P_out_sim(3,:), 'g')
xlabel('P_I (dBW)')
ylabel('P_out')
legend('gamma th= -10', 'gamma th = -5', 'gamma th = 0')
title('P_out simulated')
semilogy(-20:1:20, P_out(1,:), '*r')
hold on; 
semilogy(-20:1:20, P_out(2,:), '*g')
semilogy(-20:1:20, P_out(3,:), '*b')
ylim([ 0.01 1])
xlabel('P_I (dBW)')
ylabel('P out')
legend('gamma th= -10', 'gamma th = -5', 'gamma th = 0')
title('Nakagami-m non-relay network- Battery constraint condition')