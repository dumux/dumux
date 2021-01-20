%% Konstanten
phi=0.4;            %Porosität
rho_w=1000;         %Wasserdichte
mu_w=10^(-3);		%dynamische Viskostät des Wassers bei T=285.15k
K_int=5*10^(-12);	%intrinsische Permeabilität
g=[0;0;-9.81];		%Gravitationsbeschleunigung
p_n=10^5;           %Luftdruck

pcsweLow=10^4;
pcsweHigh=1000;
m=-5*10^5;
sThres=0.01;

%% Linearmodell-Parameter
EntryPc=0;       %[Pa]
MaxPc=10^10;    %[Pa]


%% Brooks-corey-Modell-Parameter
p_d=1000;           %Eindringdruck [Pa]
lambda=2;           %Verteilung der Porengröße 

%% Variablen
x=sym('x');                                 %Raum
globalPos1=sym('globalPos1');               %Raum y
z=sym('z');                                 %Raum
time=sym('time');                           %Zeit
pwBottom=sym('pwBottom');
pwTop=sym('pwTop');


%syms x y z t

%% Lösung Wasserdruck
%p_w=pwBottom+(pwTop-pwBottom)*(globalPos1*globalPos1*globalPos1*time/100)/64;  %einfach
p_w=pwBottom+0.5*(pwTop-pwBottom)*(1+tanh(5*globalPos1-15+time/10));            %komplex

%% Sättigung
%S_w=0.7;
%S_w=(p_n-p_w-MaxPc)/(EntryPc-MaxPc);        %Lineares Modell
%S_w=((p_n-p_w)/p_d)^(-lambda));              %Brooks-Corey-Modell nicht regularisiert

% S_w=@(p_w) ( ...
%     (p_w<=10000) .* (p_w) + ...                                       % regularization towards zero
%     and((p_w>10000),(p_w < 30000)) .* (((p_n-p_w)/p_d)^(-lambda)) + ... % Brooks-Corey part
%     (p_w >= 30000) .* (p_w) ... % regularization towards infinity
%     );             % regularisiertes Brooks-Corey-Modell


% p_c=p_n-p_w;
% S_w=@(p_c) ( ...
%             (sThres+(p_c-pcsweLow)/m) .* (p_c>=pcsweLow) + ...                       % regularization towards zero
%             (((p_c)/p_d)^(-lambda)) .* ((p_c < pcsweLow) && (p_c > pcsweHigh)) +...         % Brooks-Corey part
%             (1+(p_c-pcsweHigh)/m)   .* (p_c <=pcsweHigh)   ...                     % regularization towards infinity
%            ); 
  
p_c=p_n-p_w;
S_w=@(p_c) ( ...
            (sThres+(p_c-pcsweLow)/m)*heaviside(p_c-pcsweLow) + ...                       % regularization towards zero
            (((p_c)/p_d)^(-lambda))*(heaviside(p_c-pcsweHigh)-heaviside(p_c-pcsweLow))+...         % Brooks-Corey part
            (1+(p_c-pcsweHigh)/m)   .* heaviside((-p_c)+pcsweHigh)   ...                     % regularization towards infinity
           );            




ezplot(S_w,[0,5000]);
%% Permeabilität
%Krw=max(min(S_w,1.0),0.0);
%K=K_int*S_w;                          %Lineares Modell
K=K_int*S_w(p_c)^(2/lambda+3);              %Brooks-Corey-Modell

%% Richards-Quellterm
q_w=phi*diff((rho_w*S_w(p_c)),time)-divergence((rho_w/mu_w*K*([diff(p_w,x);diff(p_w,globalPos1);diff(p_w,z)]-rho_w*g)),[x,globalPos1,z]);
