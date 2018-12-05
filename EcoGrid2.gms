$eolcom //
SETS
P        DSO service type /
$include Sets\P.inc
/
I        Conventional Units /
$include Sets\I.inc
/
T        Time Periods /
$include Sets\T.inc
/
D        DR Blocks offered by each DR aggregator /
$include Sets\D.inc
/
C        DR aggregator /
$include Sets\C.inc
/
;


Parameter PP(p) /
Sched 1
Cond1 0.3
Cond2 0.45
/;

parameter QQ(p);
QQ('Cond2') = 0.15;


Parameters PReg(p,t) Required Response for DR service/
$include Params\PReg.inc
/
PReb(p,t) Allowed rebound for DR service/
$include Params\PReb.inc
/
PUpCon(p,i) Limit of load reduction (KW) for conventional unit/
$include Params\PUpCon.inc
/
B_DR(p,c,d) Asymmetric block limit/
$include Params\B_DR.inc
/
Q_DR(p,c,d,t) Asymmetric block definition (KW per t) for each block /
$include Params\Q_DR.inc
/
CResCon(p,i,t) Reserve cost for conventional unit (€ per KW per t) /
$include Params\CResCon.inc
/
CDispCon(p,i,t)  Dispatch cost for conventional unit (€ per KW per t)/
$include Params\CDispCon.inc
/
CResDR(p,c,d) Reserve cost for asymmetric block aggregator (€ per Block)/
$include Params\CResDR.inc
/
CDispDR(p,c,d) Dispatch cost for asymmetric block aggregator (€ per Block)/
$include Params\CDispDR.inc
/
CResReb(p,t) Reserve rebound cost for DR service (€ per KW per t)/
$include Params\CResReb.inc
/
CDispReb(p,t) Dispatch rebound cost for DR service (€ per KW per t)/
$include Params\CDispReb.inc
/
CResDSO(p) Bid-value to reserve DR service (€) to satisfy requirement p/
$include Params\CResDSO.inc
/
CDispDSO(p) Bid-value to dispatch DR service (€) to satisfy requirement p/
$include Params\CDispDSO.inc
/
;

PARAMETERS
R_DR(p,c,d)      Stores Rup
P_Con(p,i,t)     Stores Pup
S_DSO(p,t)       Stores Sup
PregZ_DSO(p,t)   Preg for chosen DR service
Ract(p,c,t)      Dispatched asymmetric block in terms of power reduction at t
Pbal_pi(p,t)     Clearing price

T_Con(p,i)    Payment from the DSO to the conventional unit
T_DR(p,c)     Payment from the DSO to the demand response unit

W_Con(p,i)    Disbenefit for each conventional unit
W_DR(p,c)     Disbenefit for each demand response unit
W_DSO(p)      Disbenefit for the DSO

// Payment adjusted parameters

h_T_Con(p,i)     Payment from the DSO to the conventional unit
h_T_DR(p,c)      Payment from the DSO to the demand response unit

h_W_Con(p,i)    Disbenefit for each conventional unit
h_W_DR(p,c)     Disbenefit for each demand response unit
h_W_DSO(p)      Disbenefit for the DSO
;

alias(t,tau);

INTEGER VARIABLE
Rup(p,c,d) Block DR block number
;

POSITIVE VARIABLES
Pup(p,i,t) Conventional DR
Sup(p,t)   DR Rebound
;

BINARY VARIABLE
z(p)       Chosen DR service
m(p,c,d)   Chosen DR block
;

VARIABLES
BC     Peak shave cost (objective function value)
RC(p)  Cost of meeting each reserve requirement (excluding operation)
OC(p)  Cost of activating reserve when they are required
;


EQUATIONS
PeakShaveCost    Balancing cost (objective function)
ReserveCost(p)   Product reserve cost
OperationCost(p) Cost to activate reserve
MaxPup(p,i,t)    Maximum up-reg - bid i time t
MaxRup(p,c,d)    Limit on DR blocks
MaxM(p,c)        Exclusivity in DR blocks
MaxSup(p,t)      Rebound limit
oneDSO           Only one DSO service dispatched
Pbal(p,t)        Service satisfied
;

// (1) Objective function
PeakShaveCost..     BC =e= sum(p, RC(p) + PP(p) * OC(p) ) ;

// (2) Defines the reserve cost
ReserveCost(p)..     RC(p) =e=  sum((i,t), CResCon(p,i,t) * Pup(p,i,t))
                    + sum((c,d), CResDR(p,c,d) * Rup(p,c,d) )
                    + sum(t, CResReb(p,t) * Sup(p,t))
                    - CResDSO(p) *   z(p) ;

// (3) Defines the dispatch cost
OperationCost(p)..   OC(p) =e=  sum((i,t), CDispCon(p,i,t) * Pup(p,i,t))
                    + sum((c,d), CDispDR(p,c,d) * Rup(p,c,d) )
                    + sum(t, CDispReb(p,t) * Sup(p,t))
                    - CDispDSO(p) * z(p) ;

// (4) Limit on conventional unit regulation
MaxPup(p,i,t).. Pup(p,i,t) =l= PUpCon(p,i) * z(p) ;

// (5) Maximum number of DR blocks dispatched
MaxRup(p,c,d).. Rup(p,c,d) =l= B_DR(p,c,d) * m(p,c,d);

// (6) Only one type of DR block is dispatched
MaxM(p,c)..     sum(d,m(p,c,d)) =l= z(p);

// (7) Limits rebound
MaxSup(p,t).. Sup(p,t) =l= Preb(p,t) * z(p);

// (8) Ensures only 1 DSO service clears
oneDSO.. sum(p, z(p)) =l= 1;

// (9) Ensures the cleared DR service is satisfied by DR units
Pbal(p,t)..         sum(i,Pup(p,i,t))
                    + sum((c,d), Q_DR(p,c,d,t) * Rup(p,c,d) )
                    =g= Preg(p,t) * z(p) - Sup(p,t);

model Ecogrid2 / PeakShaveCost
                 ReserveCost
                 OperationCost
                 MaxPup
                 MaxRup
                 MaxM
                 MaxSup
                 oneDSO
                 Pbal /;

*Solver options
FILE OPT cplex OPTION FILE /cplex.opt/
put opt;
put 'threads 0'/;
putclose;


Rup.lo(p,c,d) = 0;
solve Ecogrid2 minimizing BC using mip;
z.fx(p)       = z.l(p);
m.fx(p,c,d)   = m.l(p,c,d);
Rup.fx(p,c,d) = Rup.l(p,c,d);
solve Ecogrid2 minimizing BC using rmip;


// Storing dispatch  decisions into parameters
R_DR(p,c,d)  = Rup.l(p,c,d)  ;    //
P_Con(p,i,t) = Pup.l(p,i,t)  ;
S_DSO(p,t)       = Sup.l(p,t)  ;
PregZ_DSO(p,t)   = Preg(p,t)*z.l(p);
Ract(p,c,t)      = sum(d, Q_DR(p,c,d,t) * Rup.l(p,c,d));
Pbal_pi(p,t)     = Pbal.m(p,t);

// The payment to each aggregator
T_Con(p,i) = sum(t, Pbal.m(p,t) * Pup.l(p,i,t));
T_DR(p,c)  = max(sum(d, sum(t, Pbal.m(p,t) * Rup.l(p,c,d) * Q_DR(p,c,d,t))) ,
                 sum(d, (CResDR(p,c,d) + PP(p) * CDispDR(p,c,d)) * Rup.l(p,c,d))
                 );


// The profit of each agent for each service
W_Con(p,i) = T_Con(p,i)
           - sum(t, (CResCon(p,i,t) + QQ(p) * CDispCon(p,i,t)) * Pup.l(p,i,t));
W_DR(p,c)  = T_DR(p,c)
           - sum(d, (CResDR(p,c,d) + QQ(p) * CDispDR(p,c,d)) * Rup.l(p,c,d));
W_DSO(p)   = (CResDSO(p) + QQ(p) * CDispDSO(p)) * z.l(p)
           - sum(t, (CResReb(p,t) + QQ(p) * CDispReb(p,t)) * Sup.l(p,t) )
           - sum(i,T_Con(p,i))
           - sum(c,T_DR(p,c));

// The payment to each aggregator including dispatch payments
h_T_Con(p,i) = T_Con(p,i)
             + (QQ(p)-PP(p)) * sum(t, CDispCon(p,i,t) * Pup.l(p,i,t));
h_T_DR(p,c)  = T_DR(p,c)
             + (QQ(p)-PP(p)) * sum(d, CDispDR(p,c,d) * Rup.l(p,c,d));

// The profit of each participant for each service including dispatch payments
h_W_Con(p,i) = h_T_Con(p,i)
             - sum(t, (CResCon(p,i,t) + QQ(p) * CDispCon(p,i,t)) * Pup.l(p,i,t));
h_W_DR(p,c)  = h_T_DR(p,c)
             - sum(d, (CResDR(p,c,d) + QQ(p) * CDispDR(p,c,d)) * Rup.l(p,c,d));
h_W_DSO(p)   = sum(t,(CResReb(p,t) + QQ(p) * CDispReb(p,t)) * Sup.l(p,t))
             + (CResDSO(p) + QQ(p) * CDispDSO(p)) * z.l(p)
             - sum(i,h_T_Con(p,i)) - sum(c,h_T_DR(p,c));



Execute_Unload 'DSO_Market_Output.gdx',
                         P,I,T,D,C,
                         PReg,PUpCon,Q_DR,
                         CResCon,CDispCon,CResDR,CDispDR,CResDSO,CDispDSO,
                         PP,R_DR,P_Con,S_DSO,PregZ_DSO,Ract,
                         Pbal_pi,
                         T_Con,T_DR,
                         W_Con,W_DR,W_DSO,
                         h_T_Con,h_T_DR,
                         h_W_Con,h_W_DR,h_W_DSO;




