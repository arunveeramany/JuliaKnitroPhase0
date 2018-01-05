using Complementarity;
using JuMP;

include("def.jl");
include("dProc.jl");
include("buildMod.jl");


#-------------------------------------------------------------------------------
#INPUTS (FOR OFFLINE TESTING)
#-------------------------------------------------------------------------------
base = "/home/svcarpacomp/data/105/Phase_0_RTS96/scenario_1/"
#base = "/data/105/Phase_0_RTS96/scenario_45/"
base = "/data/141982/Phase_0_IEEE14_1Scenario/scenario_1/"

rawFile =base * "powersystem.raw";
genFile =base * "generator.csv";
contFile=base * "contingency.csv";
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#LOAD DATA FROM INPUT FILES
#comment 'function...end' for offline testing
#-------------------------------------------------------------------------------

#function MyJulia1(rawFile, genFile, contFile) 

# written by Haoxiang Yang, Northwestern University, 404-421-0638, haoxiangyang2019@u.northwestern.edu
# modified by Stephen Elbert, Jesse Holzer and Arun Veeramany, PNNL; Miles Lubin, Google
   
    println("RAW file: ",rawFile)
    println("gen file: ",genFile)
    println("con file: ",contFile)

    println("Data generation")

    baseMVA,busSeg,loadSeg,shuntSeg,genSeg,branchSeg,transformerSeg = preProc(rawFile);
    busList,busDList,genList,genDList = busgenProc(baseMVA,busSeg,shuntSeg,loadSeg,genSeg,genFile);
    brList,brListSingle,brDList = branchProc(baseMVA,branchSeg,transformerSeg);
    contList,contDList = contProc(contFile);

    fData = fixedData(baseMVA,busList,busDList,genList,genDList,brList,brListSingle,brDList);
    uData = uncertainData(contList,contDList);

    mp = buildMod(fData,uData, contDList);

  # Input:
  # fData - grid data with all parameters
  # uData - contingency data

  # readin the data from the structure
  baseMVA = fData.baseMVA;
  bList = fData.busList;
  bData = fData.busDList;

  gList = fData.genList;
  gData = fData.genDList;

  brList = fData.brList;
  brData = fData.brDList;

  S = uData.contList;
  contData = uData.contDList;

#-------------------------------------------------------------------------------
#LOOP THROUGH CONTINGENCY CASES
#-------------------------------------------------------------------------------

contingency_cases = Bool[false, true];
v0_bak = [];

for contingency_case in contingency_cases

  # set up the model
   mp = Model(solver = KnitroSolver(bar_initmu=0.12)); 
   			
  # create the variables for the base case
  @variable(mp,bData[i].Vmin <= v0[i in bList] <= bData[i].Vmax);
  @variable(mp,gData[l].Pmin <= sp0[l in gList] <= gData[l].Pmax);
  @variable(mp,gData[l].Qmin <= sq0[l in gList] <= gData[l].Qmax);
  @variable(mp,p0[k in brList]);
  @variable(mp,q0[k in brList]);
  @variable(mp,psh0[i in bList]);
  @variable(mp,qsh0[i in bList]);
  @variable(mp,theta0[i in bList]);
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#REMOVE CONTINGENCY ELEMENTS (BRANCHES FOR NOW) - FEASIBLE BUT INSIGNIFICANT
#-------------------------------------------------------------------------------
	if contingency_case == true
		for c in values(contDList)
				println("Deleting branch...",c)
				for i in 1:2 delete!( brData, [c.Loc[i]] ); end
				for i in 1:2 deleteat!( brList,findin(brList,[c.Loc[i]] )); end
		end
	end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# EQNS 22, 24	OHM'S LAW
# EQNS 23, 25 	REVERSE OF ABOVE CAUSES INFEASIBILITY, UNIMPLEMENTED
# EQ 29			THETA BOUNDS
#-------------------------------------------------------------------------------
@NLconstraint(mp,pFlow0[k in brList], 
		p0[k] == 
			(1/brData[k].tau^2) * brData[k].g* v0[brData[k].From]^2 -
            (1/brData[k].tau)* v0[brData[k].From]*v0[brData[k].To]*
             (
             	brData[k].g* cos(theta0[brData[k].From] - theta0[brData[k].To] - brData[k].thetatr) +
             	brData[k].b* sin(theta0[brData[k].From] - theta0[brData[k].To] - brData[k].thetatr)
             )
       );


 @NLconstraint(mp,qFlow0[k in brList], 
 			q0[k] == 
 			-(1.0/(brData[k].tau^2))*(brData[k].b+brData[k].bc/2)*v0[brData[k].From]^2 -
 			 (1/brData[k].tau)* v0[brData[k].From]*v0[brData[k].To]*
             (
             	brData[k].g* cos(theta0[brData[k].From] - theta0[brData[k].To] - brData[k].thetatr) -
             	brData[k].b* sin(theta0[brData[k].From] - theta0[brData[k].To] - brData[k].thetatr)
             )
 			
        );

@constraint(mp,theta_c[k in brList], 
		-pi/3.0 <= theta0[brData[k].From]-theta0[brData[k].To] <= pi/3.0);

#------------------------------------------------------------------------------------------




#------------------------------------------------------------------------------------------
#EQNS 26-27	KIRCHOFF'S LAW
#------------------------------------------------------------------------------------------
@NLconstraint(mp,pBalance0[i in bList], sum(sp0[l] for l in bData[i].gen) ==  psh0[i]*v0[i]^2 + bData[i].Pd + sum(p0[k] for k in brList if (brData[k].From == i)));
@NLconstraint(mp,qBalance0[i in bList], sum(sq0[l] for l in bData[i].gen) == -qsh0[i]*v0[i]^2 + bData[i].Qd + sum(q0[k] for k in brList if (brData[k].From == i)));
#------------------------------------------------------------------------------------------



#------------------------------------------------------------------------------------------
#APPARENT FLOW BOUND --- UNDO baseMVA^2 - FEASIBLE 
#------------------------------------------------------------------------------------------
@NLconstraint(mp,flowBound0[k in brList],p0[k]^2 + q0[k]^2 <= fData.baseMVA^2 * brData[k].t^2);
#------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------
#SHUNTS --- CAUSES INFEASIBILITY
#------------------------------------------------------------------------------------------
#@NLconstraint(mp,pShunt0[i in bList],psh0[i] == bData[i].gsh*v0[i]^2);
#@NLconstraint(mp,qShunt0[i in bList],qsh0[i] == -bData[i].bsh*v0[i]^2);
#------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------
#COMPLEMENTARITY ---- CONTINGENCY CASE ONLY --- FEASIBLE BUT INSIGNIFICANT 
#------------------------------------------------------------------------------------------
if contingency_case == true
	@variable(mp, bData[i[1]].Vmin <= vx[i in gList] <= bData[i[1]].Vmax);
	@variable(mp, v1[i in gList] >= 0);
	@variable(mp, v2[i in gList] >= 0);
	@variable(mp, gData[l].Qmin <= sq[l in gList] <= gData[l].Qmax);
	
	for i in gList
	  @complements(mp, v1[i] - v2[i], gData[i].Qmin <= sq[i] <= gData[i].Qmax)
	end
	
    @constraint(mp, vConstr[i in gList], v0_bak[i[1]] + vx[i] ==  v1[i] - v2[i]);
end
#------------------------------------------------------------------------------------------



#------------------------------------------------------------------------------------------
# OBJECTIVE FUNCTION
#------------------------------------------------------------------------------------------
@variable(mp, CF); #STE
@NLconstraint(mp, CF == sum(sum(gData[l].cParams[n]*(sp0[l]*fData.baseMVA)^n for n in gData[l].cn )        for l in gList if gData[l].Pmin<gData[l].Pmax)     #STE
                      + sum(sum(gData[l].cParams[n]*(gData[l].Pmin*fData.baseMVA)^n for n in gData[l].cn ) for l in gList if gData[l].Pmin==gData[l].Pmax) 
                        ); 
@NLobjective(mp, Min, CF);  #STE needed to avoid error for RTS96
#------------------------------------------------------------------------------------------



#------------------------------------------------------------------------------------------
#SOLVE THE OPTIMIZATION PROBLEM
#------------------------------------------------------------------------------------------
 solve(mp);
#------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------
#EXTRACT
#------------------------------------------------------------------------------------------
    sphat = getvalue(mp[:sp0]);
    sqhat = getvalue(mp[:sq0]);
    costhat = getobjectivevalue(mp);
    
    v0_bak = getvalue(mp[:v0]);
#------------------------------------------------------------------------------------------

end	#contigency_case for loop   
    
#------------------------------------------------------------------------------------------
#WRITE SOLUTION FILES
#------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------


#end		#end build function, commented for offline testing
