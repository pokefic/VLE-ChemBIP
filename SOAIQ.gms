$title Optimization of Binary Parameters NRTL

$onText
This is a GAMS-code to estimate Binary interaction parameters for VLE.
$offText

$onExternalInput
Set
    i 'Compounds'
    / ET    'Ethanol'
      WA    'Water' /
$offExternalInput
    PTexp 'Presión and Temperature'
    / Texp  'Temperature'
      Pexp  'Presure' /
    
    PsPar 'Parameters of Antoine Equation'
    / 1st
      2nd
      3rd
      4th
      5th /;
      
Alias (i,j,h);

Set       k 'Number of experimental points';


$onExternalInput
Parameter Table xpx(k<,i) 'Composition of experimental Data'
*       Liq-1        Liq-2  
          ET          WA         
    1   1e-16       1.0000      
    2   0.0170      0.9830      
    3   0.0287      0.9713      
    4   0.0516      0.9484      
    5   0.0655      0.9345      
    6   0.1042      0.8958      
    7   0.1355      0.8645      
    8   0.2243      0.7757      
    9   0.3695      0.6305      
    10  0.4924      0.5076      
    11  0.6529      0.3471      
    12  1.0000      1e-16 ;      

Table xpy(k,i) 'Composition of experimental Data'
*        Vap-1      Vap-2
          ET          WA    
    1   1e-16      1.0000
    2   0.1543      0.8457
    3   0.2627      0.7373
    4   0.3008      0.6992
    5   0.3770      0.6230
    6   0.4424      0.5576
    7   0.4846      0.5154
    8   0.5504      0.4496
    9   0.6052      0.3948
    10  0.6473      0.3527
    11  0.7212      0.2788
    12  1.0000      1e-16;
      


Table PT(k,PTexp) 'Experimental Condition'
*    Presure(Pa)   Temperature(K)
          Pexp        Texp     
    1    101.03e3    373.05
    2    101.03e3    369.80 
    3    101.03e3    366.62 
    4    101.03e3    364.36 
    5    101.03e3    362.07 
    6    101.03e3    359.64 
    7    101.03e3    358.07 
    8    101.03e3    355.92 
    9    101.03e3    353.95 
    10   101.03e3    353.24 
    11   101.03e3    352.45 
    12   101.03e3    351.38;
    
Table Psat(i,PsPar) 'Parameters of Antoine Equation for compounds'
           1st         2nd       3rd          4th         5th
    ET   73.304     -7122.3    -7.1424     2.8853e-6      2
    WA   73.649     -7258.2    -7.3037     4.1653e-6      2;
$offExternalInput

Positive Variable
    P(k)              'Presure Model'
    T(k)              'Temperature of Model'
    Ps(k,i)           'Saturation Presure of i'
    x(k,i)            'Liquid fraction of component i'
    y(k,i)            'Vapor fraction of component i'
$onExternalOutput    
    alpha(i,j)        'Binary interaction alpha'
$offExternalOutput
    Gex(k,i,j)        'G variable of NRTL'
    gama(k,i)        'Gamma variable of NRTL (activity)';

Variable
    tau(k,i,j)        'Tau variable of NRTL'
$onExternalOutput
    a(i,j)            'Binary interaction Parameter a'
    b(i,j)            'Binary interaction Parameter b'
    Obj               'Objective of optimization';
$offExternalOutput
Equation
    Ant(k,i)          'Antoine equation'
    taue(k,i,j)       'Tau equation of NRTL'
    Gexe(k,i,j)       'G equation of NRTL'
    Gam(k,i)          'Gamma Equation of NRTL'
    Eq(k,i)           'Equilibrium Equation'
    Res1(k)           'Balance y'
    Res2(k)           'Balance x'
    Res3(i,j)         'alpha verification'
    Res4(i,j)         'alpha 0'
    Res5(i)           'Diagonal'
    Res6(i)           'Diagonal'
    Obje              'Objective function';


Ant(k,i)..      log(Ps(k,i)) =e= Psat(i,'1st') + Psat(i,'2nd')/T(k) +
                Psat(i,'3rd')*log(T(k)) + Psat(i,'4th')*(T(k) ** Psat(i,'5th'));

taue(k,i,j)..   tau(k,i,j) =e= a(i,j) + b(i,j)/T(k);
Gexe(k,i,j)..   Gex(k,i,j) =e= exp(-alpha(i,j) * tau(k,i,j));

Gam(k,i)..      gama(k,i) =e= exp((sum(j, tau(k,j,i) * Gex(k,j,i) * x(k,j))/
                sum(h, Gex(k,h,i) * x(k,h))) + sum(j, (x(k,j) * Gex(k,i,j)/
                sum(h, Gex(k,h,j) * x(k,h))) * (tau(k,i,j) -
                sum(h, x(k,h) * tau(k,h,j) * Gex(k,h,j)) /
                sum(h, x(k,h) * Gex(k,h,j)))));

Eq(k,i)..       P(k) * y(k,i) =e= gama(k,i) * x(k,i) * Ps(k,i);

Res1(k)..       sum(i, y(k,i)) =e= 1;
Res2(k)..       sum(i, x(k,i)) =e= 1;
Res3(i,j)..     alpha(i,j) =e= alpha(j,i);
Res4(i,j)..     alpha(i,i) =e= alpha(i,j);
Res5(i)..       a(i,i)     =e= 0;
Res6(i)..       b(i,i)     =e= 0;   

Obje..          obj =e= sum((k,i), power(x(k,i) - xpx(k,i), 2)) +
                sum((k,i), power(y(k,i) - xpy(k,i), 2)) +
                sum(k, power(P(k) - PT(k,'Pexp'), 2)) +
                sum(k, power(T(k) - PT(k,'Texp'), 2));

* Initilization parameters
    alpha.l(i,j)  =  0.2;

    a.l(i,j)      =  5;
    b.l(i,j)      =  1000;
    a.l(i,i)      =  0;
    b.l(i,i)      =  0;
    

    P.l(k)        =  PT(k,'Pexp');
    T.l(k)        =  PT(k,'Texp');
    x.l(k,i)      =  xpx(k,i);    
    y.l(k,i)      =  xpy(k,i) ;      
    Ps.l(k,i)     =  exp(Psat(i,'1st') + Psat(i,'2nd')/T.l(k) +
                     Psat(i,'3rd')*log(T.l(k)) + Psat(i,'4th')*(T.l(k) ** Psat(i,'5th')));
    tau.l(k,i,j)  =  a.l(i,j) + b.l(i,j)/T.l(k);
    Gex.l(k,i,j)  =  exp(-alpha.l(i,j) * tau.l(k,i,j));
    
    gama.l(k,i)   =  exp((sum(j, tau.l(k,j,i) * Gex.l(k,j,i) * x.l(k,j))/
                     sum(h, Gex.l(k,h,i) * x.l(k,h))) + sum(j, (x.l(k,j) * Gex.l(k,i,j)/
                     sum(h, Gex.l(k,h,j) * x.l(k,h))) * (tau.l(k,i,j) -
                     sum(h, x.l(k,h) * tau.l(k,h,j) * Gex.l(k,h,j)) /
                     sum(h, x.l(k,h) * Gex.l(k,h,j)))));

    a.lo(i,j)     =  -10;
    a.up(i,j)     =   10;
    
    b.lo(i,j)     =  -2500;
    b.up(i,j)     =   2500;
    
    alpha.lo(i,j) =   0.2;
    alpha.up(i,j) =   0.426;
    
Model Opti / all /;

*Options NLP = IPOPT;

solve Opti using NLP minimizing obj;

display a.l,b.l,alpha.l,obj.l;

            
    