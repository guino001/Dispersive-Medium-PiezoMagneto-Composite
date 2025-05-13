
clc
close all
clear all

warning('off')

% Parameters
    Nx = 80; % Number of grid points in x-direction
    Ny = 80; % Number of grid points in y-direction
    Lx = 1;  % Length of domain in x-direction
    Ly = 1;  % Length of domain in y-direction
    dx = Lx / (Nx + 1); % Grid spacing in x-direction
    dy = Ly / (Ny + 1); % Grid spacing in y-direction
    
% Grid points
    x = linspace(0, Lx, Nx);
    y = linspace(0, Ly, Ny);
    [X, Y] = meshgrid(x, y);


% Material properties

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


% Control parameters

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


 
M11= @(x) prom(M11_1*rho1(x)*dN1_1(x)+M11_1,M11_2*rho1(x)*dN1_2(x)+M11_2);
M22= @(x) prom(M22_1*rho2(x)*dN2_1(x)+M22_1,M22_2*rho2(x)*dN2_2(x)+M22_2);

M12= @(x) prom(M11_1*rho1(x)*dN2_1(x),M11_2*rho1(x)*dN2_2(x));
M21= @(x) prom(M22_1*rho2(x)*dN1_1(x),M22_2*rho2(x)*dN1_2(x));

sfM1= @(x) prom(M22_1*drho2(x)*dcalN1_1(x),M22_2*drho2(x)*dcalN1_2(x));
sfM2= @(x) prom(M22_1*drho2(x)*dcalN2_1(x),M22_2*drho2(x)*dcalN2_2(x));

calM1= @(x) prom(M22_1*drho2(x)*dN1_1(x),M22_2*drho2(x)*dN1_2(x));
calM2= @(x) prom(M22_1*drho2(x)*dN2_1(x),M22_2*drho2(x)*dN2_2(x));


B_mat=@(x)  sfM1(x)*(calM1(x)\eye(3))*P+sfM2(x)*(calM2(x)\eye(3))*P;  %(sfM1(x)+sfM2(x))*((calM1(x)+calM2(x))\eye(3));




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



calM11=@(x) prom(M22_1*drho2(x)*dN11_1(x)+M11_1*rho1(x)*dcalN1_1(x), M22_2*drho2(x)*dN11_2(x)+M11_2*rho1(x)*dcalN1_2(x) );

calM22=@(x) prom(M22_1*drho2(x)*dN22_1(x)+M22_1*rho2(x)*dcalN2_1(x),  M22_2*drho2(x)*dN22_2(x)+M22_2*rho2(x)*dcalN2_2(x) );

calM12=@(x) prom(M22_1*drho2(x)*dN12_1(x)+M11_1*rho1(x)*dcalN2_1(x),  M22_2*drho2(x)*dN12_2(x)+M11_2*rho1(x)*dcalN2_2(x) );

calM21=@(x) prom(M22_1*drho2(x)*dN21_1(x)+M22_1*rho2(x)*dcalN1_1(x), M22_2*drho2(x)*dN21_2(x)+M22_2*rho2(x)*dcalN1_2(x) );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % construction of A matrix

A = sparse(Nx * Ny * 3, Nx * Ny * 3);


for j = 1:Ny
    
    for i = 1:Nx

        idx = (i-1)*Ny + j;

            % Diagonal term

            A(3*idx-2:3*idx, 3*idx-2:3*idx)=-(M11(y(j)))*2/dx^2 - ( M22(y(j)) )*2 /dy^2;
                        
            % Neighboring points in x1-direction
            if i > 1
                 left_idx = (i-2)*Ny + j;  % Index of left neighbor
                A(3*idx-2:3*idx, 3*left_idx-2:3*left_idx) = M11(y(j))/dx^2 - calM1(y(j))/(2*dx);

            end
            if i < Nx

                 right_idx = i*Ny + j;  % Index of right neighbor
            A(3*idx-2:3*idx, 3*right_idx-2:3*right_idx) = M11(y(j))/dx^2 + calM1(y(j))/(2*dx);

            end

            % Neighboring points in x2-direction
            if j > 1
               
                  bottom_idx = (i-1)*Ny + (j-1);  % Index of bottom neighbor
            A(3*idx-2:3*idx, 3*bottom_idx-2:3*bottom_idx) =  M22(y(j))/dy^2 - calM2(y(j))/(2*dy);
                
                          
            end
            if j < Ny
                
                top_idx = (i-1)*Ny + (j+1);  % Index of top neighbor
            A(3*idx-2:3*idx, 3*top_idx-2:3*top_idx) = M22(y(j))/dy^2 + calM2(y(j))/(2*dy);
                
            end

            % Mixed derivative terms

            if i > 1 && j > 1

                left_idx = (i-2)*Ny + j-1;  % Index of left neighbor

            A(3*idx-2:3*idx, 3*left_idx-2:3*left_idx) =( M12(y(j)) + M21(y(j)) )/(4*dx*dy);

                
            end
            if i < Nx && j > 1

                 left_idx = i*Ny + j-1;  % Index of left neighbor

            A(3*idx-2:3*idx, 3*left_idx-2:3*left_idx) = -(M12(y(j))+M21(y(j)))/(4*dx*dy);


            end
            if i > 1 && j < Ny

                 left_idx = (i-2)*Ny + j+1;  % Index of left neighbor

            A(3*idx-2:3*idx, 3*left_idx-2:3*left_idx) =-(M12(y(j))+M21(y(j)))/(4*dx*dy);


            end
            if i < Nx && j < Ny

                   left_idx = i*Ny + (j+1);  % Index of left neighbor

            A(3*idx-2:3*idx, 3*left_idx-2:3*left_idx) = (M12(y(j)) + M21(y(j)) )/(4*dx*dy);

            end

            % %  Robin BCs


            if i==1 

                left_idx =  j;  % Index of left neighbor
                A(3*idx-2:3*idx, 3*left_idx-2:3*left_idx) = M11(y(j))/dx^2 - calM1(y(j))/(2*dx) ...
                                                            +(M11(y(j))/dx^2 - calM1(y(j))/(2*dx))*( (((eye(3)+eps*N1_1_BC0(y(j))/dx)\(eps*N1_1_BC0(y(j))/dx ) ) ) );

                if j>1

                 left_idx =  j-1;  % Index of left neighbor

                A(3*idx-2:3*idx, 3*left_idx-2:3*left_idx) =( M12(y(j)) + M21(y(j)) )/(4*dx*dy)...
                                                            +(( M12(y(j)) + M21(y(j)) )/(4*dx*dy))*( (((eye(3)+eps*N1_1_BC0(y(j))/dx)\(eps*N1_1_BC0(y(j))/dx ) ) ) );
                end


                if j<Ny

                  left_idx =  j+1;  % Index of left neighbor

                A(3*idx-2:3*idx, 3*left_idx-2:3*left_idx) =-( M12(y(j)) + M21(y(j)) )/(4*dx*dy)...
                                                            -(( M12(y(j)) + M21(y(j)) )/(4*dx*dy))*( (((eye(3)+eps*N1_1_BC0(y(j))/dx)\(eps*N1_1_BC0(y(j))/dx ) ) ) );

                end

            end


             if i==Nx 

                left_idx =  (Nx-1)*Ny + j;  % Index of left neighbor
                A(3*idx-2:3*idx, 3*left_idx-2:3*left_idx) = M11(y(j))/dx^2 + calM1(y(j))/(2*dx) ...
                                                        + (M11(y(j))/dx^2 + calM1(y(j))/(2*dx))*(- ( (eps/dx)* N1_2_BC1(y(j)) ) \( eye(3)-eps*N1_2_BC1(y(j))/dx ) );


                if j>1

                 left_idx = (Nx-1)*Ny + j-1;  % Index of left neighbor

                A(3*idx-2:3*idx, 3*left_idx-2:3*left_idx) =-( M12(y(j)) + M21(y(j)) )/(4*dx*dy)...
                                                            -(( M12(y(j)) + M21(y(j)) )/(4*dx*dy))*(- ( (eps/dx)* N1_2_BC1(y(j)) ) \( eye(3)-eps*N1_2_BC1(y(j))/dx ) );

                end


                if j<Ny

                 left_idx = (Nx-1)*Ny + j+1;  % Index of left neighbor

                A(3*idx-2:3*idx, 3*left_idx-2:3*left_idx) =( M12(y(j)) + M21(y(j)) )/(4*dx*dy)...
                                                            +(( M12(y(j)) + M21(y(j)) )/(4*dx*dy))*(- ( (eps/dx)* N1_2_BC1(y(j)) ) \( eye(3)-eps*N1_2_BC1(y(j))/dx ) );

                end

             end
                   
    end
end


M1=sparse(Nx * Ny * 3, Nx * Ny * 3);



P=[rho 0 0; 0 0 0; 0 0 0];




for i = 1:Nx
    for j = 1:Ny
             k = (i-1)*Ny + j;
             idx_start = 3*(k-1) + 1;
            idx_end = 3*k;
           
            M1(idx_start:idx_end, idx_start:idx_end) = P;
           
    end
end

M2=sparse(Nx * Ny * 3, Nx * Ny * 3);



for i = 1:Nx
    for j = 1:Ny
            idx = (i-1)*Ny + j;
                M2(3*idx-2:3*idx, 3*idx-2:3*idx) = (calM11(y(j))*(M11(y(j))\eye(3))+calM22(y(j))*(M22(y(j))\eye(3)) ...
                                                    +calM12(y(j))*(M12(y(j))\eye(3))+calM21(y(j))*(M21(y(j))\eye(3)))*P;

     end
end



M=-(M1-eps*M2);


    [eigenfunctions, eigenvalues] = eigs(A, M, 40, 10);


  % %   Extract the eigenvalues and eigenvectors
eigenvalues = diag(eigenvalues);

% % Sort the eigenvalues in ascending order
eigenvalues1=real(eigenvalues);
eigenvalues1=eigenvalues1(eigenvalues1>0);
[sorted_eigenvalues, sorted_indices] = sort(real(eigenvalues1));

frq=(sqrt(sorted_eigenvalues))/(2*pi)

sorted_eigenfunctions =real( eigenfunctions(:, sorted_indices));


sorted_eigenfunctions_u = sorted_eigenfunctions(1:3:end,:);
sorted_eigenfunctions_phi = sorted_eigenfunctions(2:3:end,:);
sorted_eigenfunctions_psi = sorted_eigenfunctions(3:3:end,:);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;


n=1;

 subplot(2,3,1)
 u = reshape(sorted_eigenfunctions_u(:, n)/max(sorted_eigenfunctions_u(:, n)), Nx, Ny);

     [c,h]=contourf(X, Y, u,15);
    set(h, 'edgecolor','none');
    xlabel('$x_1$','interpreter','latex', 'FontSize',15);
    ylabel('$x_2$','interpreter','latex', 'FontSize',15);
     title(['$\hat{u}_3$, Frequency: ', num2str(frq(kf)), ' Hz'],'interpreter','latex', 'FontSize',15);
     colorbar;

    subplot(2,3,2)
     u = reshape(sorted_eigenfunctions_phi(:, n)/max(sorted_eigenfunctions_phi(:, n)), Nx, Ny);

     [c,h]=contourf(X, Y, u,15);
    set(h, 'edgecolor','none');
    xlabel('$x_1$','interpreter','latex', 'FontSize',15);
    ylabel('$x_2$','interpreter','latex', 'FontSize',15);
     title(['$\hat{\phi}$, Frequency: ', num2str(frq(kf)), ' Hz'], 'interpreter','latex', 'FontSize',15);
     colorbar;

        subplot(2,3,3)
     u = reshape(sorted_eigenfunctions_psi(:, n)/min(sorted_eigenfunctions_psi(:, n)), Nx, Ny);

      [c,h]=contourf(X, Y, u,15);
    set(h, 'edgecolor','none');
    xlabel('$x_1$','interpreter','latex', 'FontSize',15);
    ylabel('$x_2$','interpreter','latex', 'FontSize',15);
     title(['$\hat{\psi}$, Frequency: ', num2str(frq(n)), ' Hz'],'interpreter','latex', 'FontSize',15);
     colorbar;



n=3;

 subplot(2,3,4)
 u = reshape(sorted_eigenfunctions_u(:, n)/max(sorted_eigenfunctions_u(:, n)), Nx, Ny);

      [c,h]=contourf(X, Y, u,15);
    set(h, 'edgecolor','none');
    xlabel('$x_1$','interpreter','latex', 'FontSize',15);
    ylabel('$x_2$','interpreter','latex', 'FontSize',15);
     title(['$\hat{u}_3$, Frequency: ', num2str(frq(2)), ' Hz'],'interpreter','latex', 'FontSize',15);
    colorbar;

    subplot(2,3,5)
   u = reshape(sorted_eigenfunctions_phi(:, n)/max(sorted_eigenfunctions_phi(:, n)), Nx, Ny);

      [c,h]=contourf(X, Y, u,15);
    set(h, 'edgecolor','none');
    xlabel('$x_1$','interpreter','latex', 'FontSize',15);
    ylabel('$x_2$','interpreter','latex', 'FontSize',15);
     title(['$\hat{\phi}$, Frequency: ', num2str(frq(2)), ' Hz'], 'interpreter','latex', 'FontSize',15);
    colorbar;

        subplot(2,3,6)
   u = reshape(sorted_eigenfunctions_psi(:, n)/max(sorted_eigenfunctions_psi(:, n)), Nx, Ny);

     [c,h]=contourf(X, Y, u,15);
    set(h, 'edgecolor','none');
    xlabel('$x_1$','interpreter','latex', 'FontSize',15);
    ylabel('$x_2$','interpreter','latex', 'FontSize',15);
    title(['$\hat{\psi}$, Frequency: ', num2str(frq(2)), ' Hz'],'interpreter','latex', 'FontSize',15);
     colorbar;
