clc
clear all
close all

syms s


k=0.01;

omg=k*(1/s);


c3232_1=43*10^9;
c3232_2=45.3*10^9;
c3131_1=c3232_1;
c3131_2=c3232_2;
c3333_1=162*10^9;
c3333_2=269.5*10^9;

e113_1=11.6;
e113_2=0;
e223_1=e113_1;
e223_2=e113_2;
e333_1=18.6;
e333_2=0;

q113_1=0;
q113_2=550;   
q223_1=q113_1;
q223_2=q113_2;
q333_1=0;
q333_2=-699.7;    


kappa11_1=11.1*10^(-9);  
kappa11_2=0.08*10^(-9);  
kappa22_1=kappa11_1;
kappa22_2=kappa11_2;
kappa33_1=12.6*10^(-9); 
kappa33_2=0.093*10^(-9);  

alpha11_1=0;
alpha11_2=0;
alpha22_1=0;
alpha22_2=0;
alpha33_1=0;
alpha33_2=0;

mu11_1=5*10^(-6);  
mu11_2=-590*10^(-6);  
mu22_1=mu11_1;
mu22_2=mu11_2;
mu33_1=10*10^(-6);  
mu33_2=157*10^(-6);  

rho_1=5800;
rho_2=5300;




global V;
global H;
global L;

L=0.2;
V=0.5;
H=0.05;
eps=0.05;


rho=prom(rho_1,rho_2);

M11_1=[c3131_1 e113_1 q113_1; e113_1 -kappa11_1 -alpha11_1; q113_1 -alpha11_1 -mu11_1];
M11_2=[c3131_2 e113_2 q113_2; e113_2 -kappa11_2 -alpha11_2; q113_2 -alpha11_2 -mu11_2];

M22_1=[c3232_1 e223_1 q223_1; e223_1 -kappa22_1 -alpha22_1; q223_1 -alpha22_1 -mu22_1];
M22_2=[c3232_2 e223_2 q223_2; e223_2 -kappa22_2 -alpha22_2; q223_2 -alpha22_2 -mu22_2];


dN1_1= @(x) ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\( ( prom( (rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x))\eye(3),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\eye(3) ) )...
    \ prom( ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\(rho1(x)*M11_1),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\(rho1(x)*M11_2) ) -rho1(x)*M11_1  );
dN1_2= @(x) ( rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x) )\( ( prom( (rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x))\eye(3),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\eye(3) ) )...
    \prom( ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\(rho1(x)*M11_1),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\(rho1(x)*M11_2) ) -rho1(x)*M11_2  );

dN2_1= @(x) ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\( ( prom( (rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x))\eye(3),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\eye(3) ) )...
    \prom( ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\(rho2(x)*M22_1),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\(rho2(x)*M22_2) ) -rho2(x)*M22_1  );
dN2_2= @(x) ( rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x) )\( ( prom( (rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x))\eye(3),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\eye(3) ) )...
    \prom( ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\(rho2(x)*M22_1),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\(rho2(x)*M11_2) ) -rho2(x)*M22_2  );



calM1= @(x) prom(M22_1*drho2(x)*dN1_1(x),M22_2*drho2(x)*dN1_2(x));
calM2= @(x) prom(M22_1*drho2(x)*dN2_1(x),M22_2*drho2(x)*dN2_2(x));



Q1_1=@(x) (-M22_1*drho2(x)*dN1_1(x)+rho_1/rho*calM1(x))*V^2/2;
Q1_2=@(x) (-M22_1*drho2(x)*dN1_1(x)+rho_1/rho*calM1(x))*V*(1-V)+(-M22_2*drho2(x)*dN1_2(x)+rho_2/rho*calM1(x))*(1/2-V^2/2-V*(1-V));

Q2_1=@(x) (-M22_1*drho2(x)*dN2_1(x)+rho_1/rho*calM2(x))*V^2/2;
Q2_2=@(x) (-M22_1*drho2(x)*dN2_1(x)+rho_1/rho*calM2(x))*V*(1-V)+(-M22_2*drho2(x)*dN2_2(x)+rho_2/rho*calM2(x))*(1/2-V^2/2-V*(1-V));



A1=@(x)  -( prom((M11_1+rho2(x)*M22_1*rho2(x))\eye(3),(M11_2+rho2(x)*M22_2*rho2(x))\eye(3) ) )\ prom( (M11_1+rho2(x)*M22_1*rho2(x))\Q1_1(x), ...
         (M11_2+rho2(x)*M22_2*rho2(x))\Q1_2(x) );

A2=@(x)  -( prom((M11_1+rho2(x)*M22_1*rho2(x))\eye(3),(M11_2+rho2(x)*M22_2*rho2(x))\eye(3) ) )\ prom( (M11_1+rho2(x)*M22_1*rho2(x))\Q2_1(x), ...
         (M11_2+rho2(x)*M22_2*rho2(x))\Q2_2(x) );



dcalN1_1= @(x) (M11_1+rho2(x)*M22_1*rho2(x))\Q1_1(x)+(M11_1+rho2(x)*M22_1*rho2(x))\A1(x);
dcalN1_2= @(x) (M11_2+rho2(x)*M22_2*rho2(x))\Q1_2(x)+(M11_2+rho2(x)*M22_2*rho2(x))\A1(x);


dcalN2_1= @(x) (M11_1+rho2(x)*M22_1*rho2(x))\Q2_1(x)+(M11_1+rho2(x)*M22_1*rho2(x))\A2(x);
dcalN2_2= @(x) (M11_2+rho2(x)*M22_2*rho2(x))\Q2_2(x)+(M11_2+rho2(x)*M22_2*rho2(x))\A2(x);



% 
M11= @(x) prom(M11_1*rho1(x)*dN1_1(x)+M11_1,M11_2*rho1(x)*dN1_2(x)+M11_2);
M22= @(x) prom(M22_1*rho2(x)*dN2_1(x)+M22_1,M22_2*rho2(x)*dN2_2(x)+M22_2);

M12= @(x) prom(M11_1*rho1(x)*dN2_1(x),M11_2*rho1(x)*dN2_2(x));
M21= @(x) prom(M22_1*rho2(x)*dN1_1(x),M22_2*rho2(x)*dN1_2(x));

sfM1= @(x) prom(M22_1*drho2(x)*dcalN1_1(x),M22_2*drho2(x)*dcalN1_2(x));
sfM2= @(x) prom(M22_1*drho2(x)*dcalN2_1(x),M22_2*drho2(x)*dcalN2_2(x));

calM1= @(x) prom(M22_1*drho2(x)*dN1_1(x),M22_2*drho2(x)*dN1_2(x));
calM2= @(x) prom(M22_1*drho2(x)*dN2_1(x),M22_2*drho2(x)*dN2_2(x));


B_mat=@(x)  sfM1(x)*(calM1(x)\eye(3))+sfM2(x)*(calM2(x)\eye(3));  




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  first order BC

% BC at x=0

N1_1_BC0= @(x) ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\( ( prom( (rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x))\eye(3),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\eye(3) ) )...
    \prom( ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\(rho1(x)*M11_1),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\(rho1(x)*M11_2) ) - rho1(x)*M11_1  )*V^2/2 ...
    + ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\( ( prom( (rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x))\eye(3),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\eye(3) ) )...
    \prom( ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\(rho1(x)*M11_1),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\(rho1(x)*M11_2) ) - rho1(x)*M11_1  )*V*(1-V) ...
    + ( rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x) )\( ( prom( (rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x))\eye(3),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\eye(3) ) )...
    \prom( ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\(rho1(x)*M11_1),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\(rho1(x)*M11_2) ) - rho1(x)*M11_2  )*((1/2-V^2/2)-V*(1-V)) ;

% BC at x=1

N1_2_BC1= @(x) ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\( ( prom( (rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x))\eye(3),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\eye(3) ) )...
    \prom( ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\(rho1(x)*M11_1),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\(rho1(x)*M11_2) ) - rho1(x)*M11_1  )*V...
     + ( rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x) )\( ( prom( (rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x))\eye(3),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\eye(3) ) )...
    \prom( ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\(rho1(x)*M11_1),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\(rho1(x)*M11_2) ) - rho1(x)*M11_2  )*(1-V)- N1_1_BC0(x); 



B1b=@(x) N1_1_BC0(x);

B2b=@(x) ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\( ( prom( (rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x))\eye(3),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\eye(3) ) )...
    \prom( ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\(rho2(x)*M22_1),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\(rho2(x)*M22_2) ) - rho2(x)*M22_1  )*V^2/2 ...
    + ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\( ( prom( (rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x))\eye(3),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\eye(3) ) )...
    \prom( ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\(rho2(x)*M22_1),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\(rho2(x)*M22_2) ) - rho2(x)*M22_1  )*V*(1-V) ...
    + ( rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x) )\( ( prom( (rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x))\eye(3),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\eye(3) ) )...
    \prom( ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\(rho2(x)*M22_1),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\(rho2(x)*M22_2) ) - rho2(x)*M22_2  )*((1/2-V^2/2)-V*(1-V)) ;




tau1_1= @(x) dN1_1(x);   
tau1_2= @(x) dN1_2(x); 

Q11_1_coef=@(x)  ( -M11_1*rho1(x)*tau1_1(x)-M11_1+rho_1/rho*M11(x) ); 
Q11_2_coef=@(x)  ( -M11_2*rho1(x)*tau1_2(x)-M11_2+rho_2/rho*M11(x) ); 

Q11_1=@(x) Q11_1_coef(x)*V^2/2;

Q11_2=@(x) Q11_1_coef(x)*V*(1-V)+ Q11_2_coef(x)*((1/2-V^2/2)-V*(1-V));

tau2_1= @(x) dN2_1(x);    
tau2_2= @(x) dN2_2(x); 

Q22_1_coef=@(x)  ( -M22_1*rho2(x)*tau2_1(x)-M22_1+rho_1/rho*M22(x) ); 
Q22_2_coef=@(x)  ( -M22_2*rho2(x)*tau2_2(x)-M22_2+rho_2/rho*M22(x) ); 

Q22_1=@(x) Q22_1_coef(x)*V^2/2;

Q22_2=@(x) Q22_1_coef(x)*V*(1-V)+ Q22_2_coef(x)*((1/2-V^2/2)-V*(1-V));


Q12_1_coef=@(x)  ( -M11_1*rho1(x)*tau2_1(x)+rho_1/rho*M12(x) ); 
Q12_2_coef=@(x)  ( -M11_2*rho1(x)*tau2_2(x)+rho_2/rho*M12(x) ); 

Q12_1=@(x) Q12_1_coef(x)*V^2/2;

Q12_2=@(x) Q12_1_coef(x)*V*(1-V)+ Q12_2_coef(x)*((1/2-V^2/2)-V*(1-V));


Q21_1_coef=@(x)  ( -M22_1*rho2(x)*tau1_1(x)+rho_1/rho*M21(x) ); 
Q21_2_coef=@(x)  ( -M22_2*rho2(x)*tau1_2(x)+rho_2/rho*M21(x) ); 

Q21_1=@(x) Q21_1_coef(x)*V^2/2;

Q21_2=@(x) Q21_1_coef(x)*V*(1-V)+ Q21_2_coef(x)*((1/2-V^2/2)-V*(1-V));



N1_1_coef=@(x) ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\( ( prom( (rho1(x)*M11_1*rho1(x) ...
        +rho2(x)*M22_1*rho2(x))\eye(3),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\eye(3) ) )\prom( ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\(rho1(x)*M11_1), ...
        (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\(rho1(x)*M11_2) ) - rho1(x)*M11_1);

N1_2_coef=@(x) ( rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x) )\( ( prom( (rho1(x)*M11_1*rho1(x) ...
        +rho2(x)*M22_1*rho2(x))\eye(3),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\eye(3) ) )\prom( ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\(rho1(x)*M11_1), ...
        (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\(rho1(x)*M11_2) ) - rho1(x)*M11_2);

N1_1=@(x)  N1_1_coef(x)*V^2/2 - B1b(x)*V ;

N1_2=@(x) N1_1_coef(x)*V*(1-V)+N1_2_coef(x)*((1/2-V^2/2)-V*(1-V))- B1b(x)*(1-V);


N2_1_coef=@(x) ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\( ( prom( (rho1(x)*M11_1*rho1(x) ...
        +rho2(x)*M22_1*rho2(x))\eye(3),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\eye(3) ) )\prom( ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\(rho2(x)*M22_1), ...
        (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\(rho2(x)*M22_2) ) - rho2(x)*M22_1);

N2_2_coef=@(x) ( rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x) )\( ( prom( (rho1(x)*M11_1*rho1(x) ...
        +rho2(x)*M22_1*rho2(x))\eye(3),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\eye(3) ) )\prom( ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\(rho2(x)*M22_1), ...
        (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\(rho2(x)*M22_2) ) - rho2(x)*M22_2);

N2_1=@(x)  N2_1_coef(x)*V^2/2 - B2b(x)*V ;

N2_2=@(x) N2_1_coef(x)*V*(1-V)+N2_2_coef(x)*((1/2-V^2/2)-V*(1-V))- B2b(x)*(1-V);




calA_11=@(x) -( prom( (rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x))\eye(3),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\eye(3) ) ) ...
        \prom( (rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x))\( Q11_1(x)-rho1(x)*M11_1*N1_1(x) ), (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\( Q11_2(x)-rho1(x)*M11_2*N1_2(x) ) );

calA_22=@(x) -( prom( (rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x))\eye(3),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\eye(3) ) ) ...
        \prom( (rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x))\( Q22_1(x)-rho2(x)*M22_1*N2_1(x) ), (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\( Q22_2(x)-rho2(x)*M22_2*N2_2(x) ) );

calA_12=@(x) -( prom( (rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x))\eye(3),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\eye(3) ) ) ...
        \prom( (rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x))\( Q12_1(x) -rho1(x)*M11_1*N2_1(x) ), (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\( Q12_2(x) -rho1(x)*M11_2*N2_2(x) ) );

calA_21=@(x) -( prom( (rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x))\eye(3),  (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\eye(3) ) ) ...
        \prom( (rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x))\( Q21_1(x) -rho2(x)*M22_1*N1_1(x) ), (rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x))\( Q21_2(x) -rho2(x)*M22_2*N1_2(x) ) );




dN11_1=@(x) ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\(Q11_1(x)-rho1(x)*M11_1*N1_1(x) )+ ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\calA_11(x); 

dN11_2=@(x) ( rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x) )\(Q11_2(x)-rho1(x)*M11_2*N1_2(x) )+ ( rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x) )\calA_11(x); 

dN22_1=@(x) ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\(Q22_1(x)-rho2(x)*M22_1*N2_1(x) )+ ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\calA_22(x); 

dN22_2=@(x) ( rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x) )\(Q22_2(x)-rho2(x)*M22_2*N2_2(x) )+ ( rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x) )\calA_22(x); 

dN12_1=@(x) ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\(Q12_1(x)-rho1(x)*M11_1*N2_1(x) )+ ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\calA_12(x); 

dN12_2=@(x) ( rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x) )\(Q12_2(x)-rho1(x)*M11_2*N2_2(x) )+ ( rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x) )\calA_12(x); 

dN21_1=@(x) ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\(Q21_1(x)-rho2(x)*M22_1*N1_1(x) )+ ( rho1(x)*M11_1*rho1(x)+rho2(x)*M22_1*rho2(x) )\calA_21(x); 

dN21_2=@(x) ( rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x) )\(Q21_2(x)-rho2(x)*M22_2*N1_2(x) )+ ( rho1(x)*M11_2*rho1(x)+rho2(x)*M22_2*rho2(x) )\calA_21(x); 



calM11=@(x) prom(M22_1*drho2(x)*dN11_1(x)+M11_1*rho1(x)*dcalN1_1(x), ...
                   M22_2*drho2(x)*dN11_2(x)+M11_2*rho1(x)*dcalN1_2(x) );

calM22=@(x) prom(M22_1*drho2(x)*dN22_1(x)+M22_1*rho2(x)*dcalN2_1(x), ...
                   M22_2*drho2(x)*dN22_2(x)+M22_2*rho2(x)*dcalN2_2(x) );

calM12=@(x) prom(M22_1*drho2(x)*dN12_1(x)+M11_1*rho1(x)*dcalN2_1(x), ...
                   M22_2*drho2(x)*dN12_2(x)+M11_2*rho1(x)*dcalN2_2(x) );

calM21=@(x) prom(M22_1*drho2(x)*dN21_1(x)+M22_1*rho2(x)*dcalN1_1(x), ...
                   M22_2*drho2(x)*dN21_2(x)+M22_2*rho2(x)*dcalN1_2(x) );


P=[rho 0 0; 0 0 0; 0 0 0];



m1=0.1;


kk=1;

for theta=0:1:360

   eqn = ( M11(m1)*k^2*cosd(theta)^2+M12(m1)*k^2*cosd(theta)*sind(theta)+M21(m1)*k^2*cosd(theta)*sind(theta)+M22(m1)*k^2*sind(theta)^2) ...
           -1i*(calM1(m1)*cosd(theta)+calM2(m1)*sind(theta))*k ...
                   - omg^2*(P+eps*(calM11(m1)*(M11(m1)\eye(3))+calM12(m1)*(M12(m1)\eye(3))+calM21(m1)*(M21(m1)\eye(3))+calM22(m1)*(M22(m1)\eye(3)))*P);
 

exp=det(eqn);

rr=real(simplify(solve(exp,s)));

if rr(1)>rr(2)
h(kk)=rr(1);
else
    h(kk)=rr(2);
end

kk=kk+1;

end

theta=linspace(0,2*pi,361);

h1=polar( theta, h);


set(h1, 'color', 'g','linewidth',2,'linestyle', '--');
title('')


hold on

%%


syms s


k=0.01;

omg=k*(1/s);


m1=0.125;


kk=1;

for theta=0:1:360

   eqn = ( M11(m1)*k^2*cosd(theta)^2+M12(m1)*k^2*cosd(theta)*sind(theta)+M21(m1)*k^2*cosd(theta)*sind(theta)+M22(m1)*k^2*sind(theta)^2) ...
           -1i*(calM1(m1)*cosd(theta)+calM2(m1)*sind(theta))*k ...
                   - omg^2*(P+eps*(calM11(m1)*(M11(m1)\eye(3))+calM12(m1)*(M12(m1)\eye(3))+calM21(m1)*(M21(m1)\eye(3))+calM22(m1)*(M22(m1)\eye(3)))*P);
 

 
exp=det(eqn);

rr=real(simplify(solve(exp,s)));

if rr(1)>rr(2)
h(kk)=rr(1);
else
    h(kk)=rr(2);
end

kk=kk+1;

end


theta=linspace(0,2*pi,361);

h2=polar( theta, h);


set(h2, 'color', 'b','linewidth',2,'linestyle','--');
title('')


hold on



%%


syms s


k=0.01;

omg=k*(1/s);


m1=0.15;


kk=1;

for theta=0:1:360

   eqn = ( M11(m1)*k^2*cosd(theta)^2+M12(m1)*k^2*cosd(theta)*sind(theta)+M21(m1)*k^2*cosd(theta)*sind(theta)+M22(m1)*k^2*sind(theta)^2) ...
           -1i*(calM1(m1)*cosd(theta)+calM2(m1)*sind(theta))*k ...
                   - omg^2*(P+eps*(calM11(m1)*(M11(m1)\eye(3))+calM12(m1)*(M12(m1)\eye(3))+calM21(m1)*(M21(m1)\eye(3))+calM22(m1)*(M22(m1)\eye(3)))*P);
 
    
 
exp=det(eqn);

rr=real(simplify(solve(exp,s)));

if rr(1)>rr(2)
h(kk)=rr(1);
else
    h(kk)=rr(2);
end

kk=kk+1;

end


theta=linspace(0,2*pi,361);

h3=polar( theta, h);


set(h3, 'color', 'r','linewidth',2,'linestyle','--');
title('')


hold on


%%


syms s


k=0.01;

omg=k*(1/s);


m1=0.175;


kk=1;

for theta=0:1:360

   eqn = ( M11(m1)*k^2*cosd(theta)^2+M12(m1)*k^2*cosd(theta)*sind(theta)+M21(m1)*k^2*cosd(theta)*sind(theta)+M22(m1)*k^2*sind(theta)^2) ...
           -1i*(calM1(m1)*cosd(theta)+calM2(m1)*sind(theta))*k ...
                   - omg^2*(P+eps*(calM11(m1)*(M11(m1)\eye(3))+calM12(m1)*(M12(m1)\eye(3))+calM21(m1)*(M21(m1)\eye(3))+calM22(m1)*(M22(m1)\eye(3)))*P);
 

 
exp=det(eqn);

rr=real(simplify(solve(exp,s)));

if rr(1)>rr(2)
h(kk)=rr(1);
else
    h(kk)=rr(2);
end

kk=kk+1;

end


theta=linspace(0,2*pi,361);

h4=polar( theta, h);


set(h4, 'color', 'm','linewidth',2,'linestyle','--');
title('')


hold on


%%


syms s


k=0.01;

omg=k*(1/s);


m1=0.2;


kk=1;

for theta=0:1:360

    eqn = ( M11(m1)*k^2*cosd(theta)^2+M12(m1)*k^2*cosd(theta)*sind(theta)+M21(m1)*k^2*cosd(theta)*sind(theta)+M22(m1)*k^2*sind(theta)^2) ...
           -1i*(calM1(m1)*cosd(theta)+calM2(m1)*sind(theta))*k ...
                   - omg^2*(P+eps*(calM11(m1)*(M11(m1)\eye(3))+calM12(m1)*(M12(m1)\eye(3))+calM21(m1)*(M21(m1)\eye(3))+calM22(m1)*(M22(m1)\eye(3)))*P);
 

 
exp=det(eqn);

rr=real(simplify(solve(exp,s)));

if rr(1)>rr(2)
h(kk)=rr(1);
else
    h(kk)=rr(2);
end

kk=kk+1;

end


theta=linspace(0,2*pi,361);

h5=polar( theta, h);


set(h5, 'color', 'k','linewidth',2,'linestyle','--');
title('')


hold on
