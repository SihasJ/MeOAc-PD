$title CENG0013 Project - Process & Cost Analysis for Level 4 Douglas

* ---- Defining all the components and streams ---- *
Sets
$ontext
* ============================================== *
*    Components are  CO   - Carbon Monoxide      *
*                    H2   - Hydrogen             *
*                   DME   - Dimethyl Ether       *
*                   MeOH  - Methanol             *
*                   MeOAc - Methyl Acetate       *
* ============================================== *
$offText
c      set of components   /H2,CO,DME,MeOAc,MeOH/

* For the distillation column
cl(c)  subset of components lighter than the light key  /CO,H2/

$ontext
* ====================================================== *
*  Streams are   F1     -  CO Feed (+H2)                 *
*                F2     -  DME Feed (+MeOH)              *
*                Rin  - Stream to reactor                *
*                Rout - Reactor effluent                 *
*                V  - Vapour stream from flash drum      *
*                L  - Liquid stream from flash drum      *
*                Recycle - Recycle stream from splitter  *
*                Purge  - Purge Stream                   *
*                Product - Product Stream                *
*                DCTop  - Top distillation product       *
* ====================================================== *
$offText
s       set of streams
        /F1,F2,Rin,Rout,V,L,Recycle,Purge,Product,DCTop/;

Positive Variables
n(s,c)  molar flow rate of component c in stream s (kmol per h)
xi      extent of reaction
r       distillation column recovery amount
FR(s)   Total Stream Flow Rate of stream s
x(c)    Molar composition of liquid stream from the flash drum
y(c)    Molar composition of vapour stream from the flash drum
z(c)    Molar composition of stream to the flash drum;

Parameters
Reactor_P     Pressure of the reactor [bar]         /26/
Flash_T       System temperature [K]                /460/
Flash_P       System pressure [bar]                 /1/
sf            Split fraction to the recycle stream  /0.9/
spc           Single pass conversion (SPC) (%)
nu(c)         Stoichiometry for reaction ( DME + CO -> MeOAc )
              /DME   -1
               CO    -1
               MeOAc  1
               H2     0
               MeOH   0/
pstar(c)    Vapour pressure for each component [bar] in the flash drum
gfunc       Gamma function for compressor costing;

* Antoine coefficients for log_10 (P [bar]) = A - B / (T[K] + C)
* to determine flash drum distributions

Sets
coeff  /A, B, C/;

Table Ant(c,coeff) Antoine coefficients
              A          B           C
H2         3.54314     99.395       7.726
CO         3.36515    230.272     -13.15
DME        4.11475    894.669     -30.604
MeOAc      4.20364   1164.426     -52.69
MeOH       5.15853   1569.613     -34.846;

pstar(c) = 10**(Ant(c,'A') - Ant(c,'B') / ((Flash_T) + Ant(c,'C')));

Equations

* ======== Constraint & Specification Equations ======== *

f1comp           F1 composition specifications
f2comp           F2 composition specifications
TSFR(s)          Total Stream Flow Rate [kmol per h]
purity           Purity Constraint
minco            Minimum carbon monoxide entering reactor constraint;

* Stream F1
n.fx('F1','DME') = 0;      
n.fx('F1','MeOH') = 0;     
n.fx('F1','MeOAc') = 0;    
* Stream F2
FR.fx('F2') = 200;
n.fx('F2','CO') = 0;
n.fx('F2','H2') = 0;
n.fx('F2','MeOAc') = 0;

minco..     n('Rin','CO') =g= n('Rin','DME');
f1comp..    n('F1','H2') =e= 0.02 * FR('F1');
f2comp..    n('F2','MeOH') =e= 0.001 * FR('F2');

* Total stream flow rate equations

TSFR(s)..   sum(c,n(s,c)) =e= FR(s);

* =========== Single Pass Conversion equations =========== *
Equations
spcdef           Single pass conversion definition equation;

spcdef..     n('Rin','DME') * spc =e= xi;

spc = (1.3373 * Reactor_P + 4.0682)/100;


* ============ Mass balance equations ============ *
Equations
mbmix(c)         Mass balance of the mixer equation
mbreac(c)        Mass balance of the reactor equation
mbsplit(c)       Mass balance of the splitter equation
mbfd(c)          Mass balance of the flash drum
mbdc(c)          Mass balance of the distillation column;

mbmix(c)..   n('F1',c) + n('F2',c) + n('Recycle',c) + n('DCTop',c) =e= n('Rin',c);
mbreac(c)..  n('Rout',c) =e= n('Rin',c) + nu(c) * xi;
mbsplit(c).. n('V',c) =e= n('purge',c) + n('recycle',c);
mbfd(c)..    FR('Rout') * z(c) =e= FR('L') * x(c) + FR('V') * y(c);
mbdc(c)..    n('L',c) =e= n('Product',c) + n('DCTop',c);

* ============ Flash Drum Equations ============= *
Equations
xcomp(c)         Composition equation of stream L
ycomp(c)         Composition equation of stream V
zcomp(c)         Composition equation of stream Rout
xt               Total composition balance of x
yt               Total composition balance of y
zt               Total composition balance of z
raoults_law(c)   Raoult's law for each component;

xcomp(c)..   x(c) * FR('L') =e= n('L',c);
ycomp(c)..   y(c) * FR('V') =e= n('V',c);
zcomp(c)..   z(c) * FR('Rout') =e= n('Rout',c);

xt..         sum(c,x(c)) =e= 1;
yt..         sum(c,y(c)) =e= 1;
zt..         sum(c,z(c)) =e= 1;

raoults_law(c).. y(c) * Flash_P =e= x(c) * pstar(c);

n.l('L',c) = 0.5 * n.l('Rout',c);
n.l('V',c) = 0.5 * n.l('Rout',c);

* ============ Distillation column balances ============= *
Equations
dclkt            Distillation column recovery of the light key to tops
dclkb            Distillation column recovery of the light key to bottoms
dchkb            Distillation column recovery of the heavy key to bottoms
dchkt            Distillation column recovery of the heavy key to tops

dcllk            Distillation column for components lighter than the light key
dchhk            Distillation column for components heavier than the heavy key;

* Recovery fraction constraint
r.lo = 0.9;
r.up = 0.998;

dclkt..     n('DCTop','DME') =e= r * n('L','DME');
dclkb..     n('Product','DME') =e= (1-r) * n('L','DME');
dchkb..     n('Product','MeOAc') =e= r * n('L','MeOAc');
dchkt..     n('DCTop','MeOAc') =e= (1-r) * n('L','MeOAc');
dcllk(cl).. n('DCTop',cl) =e= n('L',cl);  
dchhk..     n('Product','MeOH') =e= n('L','MeOH');

purity..    n('Product','MeOAc') =g= 0.99 * FR('Product');

n.fx('Product','H2') = 0;
n.fx('Product','CO') = 0;
n.fx('DCTop','MeOH') = 0;

* ============ Splitter Equations ============ *
Equations
splitfrac(c)     Split fraction equations of recycled proportion for splitter;

splitfrac(c).. n('V',c) * sf =e= n('recycle',c);

$onText
* =============== Cost Function =============== *
*  This section is for calculating the cost for *
* the purpose of maximising the profit function *
*  to determine the F1 flowrate, distillation   *
*  recovery and the temperatures and pressures  *
*       of the flash drum and reactor           *
* ============================================= *
$offText

Equation
profiteq      overall profit equation
matnetprof    material net profit equation
DCCosteq      Cost equation of the distillation column
FDCosteq      Cost equation of the flash drum
RCosteq       Cost equation of the reactor
CCosteq       Cost from flow to the compressor
CWCosteq      Cost from work of the compressor
C1Workeq      Work equation of the compressor of the F1 stream
C2Workeq      Work equation of the compressor of the F2 stream
C3Workeq      Work equation of the compressor of the Recycle stream
C4Workeq      Work equation of the compressor of the DCTop stream
unitnetcost   total unit costs;

Parameters
* ---- Value of the streams in USD per kmol ---- *
CostFG        Feed gas cost                         /13.2/
CostFL        Feed liquid cost                      /78.2/
PriceP        Product price                         /135/
OPTime        Operational Time in hours per year    /8150/
RV            relative volatilities;

Variables
DCCost        cost of distillation column in MUSD per annum
FDCost        cost of flash drum in in MUSD per annum
RCost         cost of reactor in in MUSD per annum
CCost         cost of compressors in MUSD per annum
profit        total profit
MAT           material net profit
UNIT          total unit cost
CWCost        work cost for the compressor in MUSD per annum
C1Work        Work equation of the compressor of the F1 stream
C2Work        Work equation of the compressor of the F2 stream
C3Work        Work equation of the compressor of the Recycle stream
C4Work        Work equation of the compressor of the DCTop stream;

* The temperature of the feed to the distillation column will be assumed to be
* the same as the flash drum to calculate a reasonable relative volatility

RV = pstar('DME')/pstar('MeOAc');

* Costs for the major units (except heat exchanger)
DCCosteq.. DCCost =e= sqrt(FR('L'))/((RV-1)*(100*(1-r)));
FDCosteq.. FDCost =e= 0.001 * FR('Rout');
RCosteq .. RCost  =e= FR.l('Rin') * 13248*10**(-6);
CCosteq .. CCost  =e= 0.02 * (FR('F1') + FR('F2') + FR('Recycle') + FR('DCTop'));

* Balances for the components costs/revenue and unit costs 
matnetprof.. MAT =e= (PriceP * FR('Product') - CostFL * FR('F2') - CostFG * FR('F1')) * OPtime * 1e-6;
unitnetcost.. UNIT =e= DCCost + FDCost + RCost + CCost + CWCost;
profiteq.. profit =e= MAT - UNIT;

* Gamma Function for compressor work estimation
gfunc = 1.4/(1.4 - 1);

* Compressor Costs from Work
* including unit conversions, actual power requirements for the compressors
* and 10% more than destination pressure
CWCosteq.. CWCost =e= (C1Work + C2Work + C3Work + C4Work) * 3600 * OPTime * 1e-9 * 1e-6 * 17;
C1Workeq.. C1Work =e= 1.562 * gfunc * 8.314 * (25 + 273.15) * FR('F1')  * 1000/3600 * ((Reactor_P*1.1/1)**(1/gfunc) - 1);
C2Workeq.. C2Work =e= 1.562 * gfunc * 8.314 * (25 + 273.15) * FR('F2')  * 1000/3600 * ((Reactor_P*1.1/1)**(1/gfunc) - 1);
C3Workeq.. C3Work =e= 1.562 * gfunc * 8.314 * (Flash_T) * FR('Recycle') * 1000/3600 * (((Reactor_P*1.1)/Flash_P)**(1/gfunc) - 1);
C4Workeq.. C4Work =e= 1.562 * gfunc * 8.314 * (Flash_T) * FR('DCTop')   * 1000/3600 * (((Reactor_P*1.1)/Flash_P)**(1/gfunc) - 1);

Model process /all/;

* Creates an output file for easy data maniplulation
file results /outputfile/;
put results;

put "    profit      F1 Flow Rate      r         sf       Flash_P    Flash_T     Reactor_P     Distillate Flowrate"/;

* Iterative loop that goes over split fraction, reactor pressure,
* flash drum temperature and flash drum pressure

for (sf = 0.93 to 0.99 by 0.01,
    for (Reactor_P = 20 to 30 by 2,
        for (Flash_T = 310 to 410 by 20,
            for (Flash_P = 20 to 29 by 1,
            
*               After each iteration, the values for spc, pstar and RV are recalculated

                spc = (1.3373 * Reactor_P + 4.0682)/100;
                pstar(c) = 10**(Ant(c,'A') - Ant(c,'B') / ((Flash_T) + Ant(c,'C')));
                RV = pstar('DME')/pstar('MeOAc');
                
                solve process using NLP maximizing profit;
                
*               The if-statement below filters the results to only put
*               optimal solutions into the data file

                if (process.modelstat eq 2,
*                   The output contains the required design variables
                    put profit.l FR.l('F1') r.l sf Flash_P Flash_T Reactor_P FR.l('Product')/;
                );
            );
        );
    );
);

$ontext
* =================================================================================== *
*  After running the above code and finding the maximum profit, the design variables  *
*      are reapplied to find the molecular flow rates of each individual stream       *
* =================================================================================== *

sf = 0.95;
Reactor_P = 30;
Flash_P = 29;
Flash_T = 310;

spc = (1.3373 * Reactor_P + 4.0682)/100;
pstar(c) = 10**(Ant(c,'A') - Ant(c,'B') / ((Flash_T) + Ant(c,'C')));
RV = pstar('DME')/pstar('MeOAc');

Solve process using NLP maximizing profit;
Display FR.l, n.l;
$offText
