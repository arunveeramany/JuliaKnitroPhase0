using Complementarity;
using JuMP;
#using Ipopt;

include("def.jl");
include("dProc.jl");
include("buildMod.jl");


#-------------------------------------------------------------------------------
#INPUTS (FOR OFFLINE TESTING)
#-------------------------------------------------------------------------------
basedata = "/home/svcarpacomp/data"
basedata = "/data"
base = basedata * "/105/Phase_0_RTS96/scenario_45/"
#base = basedata * "/141982/Phase_0_IEEE14_1Scenario/scenario_1/"

#rawFile =base * "powersystem.raw";
#genFile =base * "generator.csv";
#contFile=base * "contingency.csv";
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#LOAD DATA FROM INPUT FILES
#comment 'function...end' for offline testing
#-------------------------------------------------------------------------------

function MyJulia1(rawFile, genFile, contFile) 

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

  # reading data from the structure
  baseMVA = fData.baseMVA;
  bList = fData.busList;
  bData = fData.busDList;

  gList = fData.genList;
  gData = fData.genDList;

  brList = fData.brList;
  brData = fData.brDList;

  S = uData.contList;
  contData = uData.contDList;
  

sol2file = open("solution2.txt","w") 
sol2gen = String[]
sol2bus = String[]
sol2branch = String[]

#-------------------------------------------------------------------------------
#LOOP THROUGH CONTINGENCY CASES
#-------------------------------------------------------------------------------

contingency_cases = collect(0:length(contDList));
v0_base = [];
p0_base = [];
q0_base = [];
sp0_base = [];
sq0_base = [];
theta0_base = [];
costhat_base = 0;

brData_backup = copy(brData);
brList_backup = copy(brList);

for contingency_case in contingency_cases

  # set up the model
   mp = Model(solver = KnitroSolver(bar_initmu=0.12)); 
   #mp = Model(solver = IpoptSolver());
   			
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

	if contingency_case > 0
	
		brData = brData_backup;
		brList = brList_backup;
	
		c = contDList[ contingency_case ];
		println("Deleting branch...",c)
		for i in 1:2 delete!( brData, [c.Loc[i]] ); end
		#for i in 1:2 deleteat!( brList,findin(brList,[c.Loc[i]] )); end
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
#EQNS 43,44 MERGED - APPARENT FLOW BOUND --- UNDO baseMVA^2 - FEASIBLE 
#------------------------------------------------------------------------------------------
@NLconstraint(mp,flowBound0[k in brList],p0[k]^2 + q0[k]^2 <= fData.baseMVA^2 * brData[k].t^2);
#------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------
#SHUNTS --- CAUSES INFEASIBILITY --- gsh ZERO
#------------------------------------------------------------------------------------------
#@NLconstraint(mp,pShunt0[i in bList],psh0[i] == bData[i].gsh*v0[i]^2);
#@NLconstraint(mp,qShunt0[i in bList],qsh0[i] == -bData[i].bsh*v0[i]^2);
#------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------
#COMPLEMENTARITY ---- CONTINGENCY CASE ONLY --- FEASIBLE BUT INSIGNIFICANT 
#MODIFY TO CONSIDER ONE CONTINGENCY AT A TIME
#------------------------------------------------------------------------------------------
if contingency_case == -1 # > 0
	@variable(mp, bData[i[1]].Vmin <= vx[i in gList] <= bData[i[1]].Vmax);
	@variable(mp, v1[i in gList] >= 0);
	@variable(mp, v2[i in gList] >= 0);
	@variable(mp, gData[l].Qmin <= sq[l in gList] <= gData[l].Qmax);
	
	for i in gList
	  @complements(mp, v1[i] - v2[i], gData[i].Qmin <= sq[i] <= gData[i].Qmax)
	end
	
    #@constraint(mp, vConstr[i in gList], v0_base[i[1]] + vx[i] ==  v1[i] - v2[i]);
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
    v0hat = 	getvalue(mp[:v0]);
#    theta0hat = getvalue(mp[:theta0]/pi*180.0);

if contingency_case == 0
	  open("solution1.txt","w") do f
	      write(f, "--generation dispatch \nbus id,unit id,pg(MW),qg(MVar) \n");
	      for i in fData.genList
	        loc = fData.genDList[i].Loc;
	        name = fData.genDList[i].Name;
	        spTemp = sphat[i]*fData.baseMVA;
	        sqTemp = sqhat[i]*fData.baseMVA;
	        write(f, "$loc,$name,$spTemp,$sqTemp \n");
	      end
	      write(f,"--end of generation dispatch \n");
	    end
	        
    	v0_base = getvalue(mp[:v0]);		#------backup base case solution
    	sp0_base = getvalue(mp[:sp0]);
		sq0_base = getvalue(mp[:sq0]);
		theta0_base = getvalue(mp[:theta0]);
		costhat_base = getobjectivevalue(mp);
		p0_base = getvalue(mp[:p0]);
		q0_base = getvalue(mp[:q0]);
		
else
	counter=0
	 for i in fData.genList
          counter += 1;
          loc = fData.genDList[i].Loc;
          idTemp = fData.genDList[i].ID;
          name = fData.genDList[i].Name;
          sqTemp = sqhat[i]*fData.baseMVA;
          push!(sol2gen, "$contingency_case,l_$counter,$loc,$name,$sqTemp \n")
     end
     
     for i in fData.busList
          id = fData.busDList[i].ID;
          name = fData.busDList[i].Name;
          vTemp = v0hat[i];
          thetaTemp = getvalue(mp[:theta0][i]/pi*180)
          push!(sol2bus,"$contingency_case,$id,$vTemp,$thetaTemp \n");
     end
     
     for i in fData.brListSingle
        idTemp = fData.brDList[i].ID;
        revidTemp = fData.brDList[i].revID;
        fromTemp = fData.brDList[i].From;
        toTemp = fData.brDList[i].To;
        name = fData.brDList[i].CKT;
        pTemp = getvalue(mp[:p0][i])*fData.baseMVA;
        qTemp = getvalue(mp[:q0][i])*fData.baseMVA;
        pRevTemp = 0 #getvalue(mp[:p][revidTemp,s])*fData.baseMVA;
        qRevTemp = 0 #getvalue(mp[:q][revidTemp,s])*fData.baseMVA;
        push!(sol2branch,"$contingency_case,$name,$fromTemp,$toTemp,i_$(fromTemp)_$(toTemp)_$(name),$pTemp,$qTemp,$pRevTemp,$qRevTemp \n");
     end
end	#if contingency_case
   
#------------------------------------------------------------------------------------------

end	#contigency_case for loop

write(sol2file, "--contingency generator \nconID,genID,busID,unitID,q(MW) \n");

for line in sol2gen
	write(sol2file, line);
end

write(sol2file,"--end of contingency generator \n--bus \ncontingency id,bus id,v(pu),theta(deg) \n");

for i in fData.busList
        id = fData.busDList[i].ID;
        name = fData.busDList[i].Name;
        vTemp = (v0_base[i]);
        thetaTemp = (theta0_base[i]/pi*180);
        write(sol2file,"0,$id,$vTemp,$thetaTemp \n");
end
      
for line in sol2bus
	write(sol2file, line);
end

write(sol2file,"--end of bus \n--Delta \ncontingency id,Delta(MW) \n");

for s in 1:length(contingency_cases)-1
        pdeltaTemp = 0 #getvalue(mp[:pdelta][s])*fData.baseMVA;
        write(sol2file,"$s,$pdeltaTemp \n");
end

write(sol2file,"--end of Delta \n--line flow \ncontingency id,line id,origin bus id,destination bus id,circuit id,p_origin(MW),q_origin(MVar),p_destination(MW),q_destination(MVar) \n");

for i in fData.brListSingle
        idTemp = fData.brDList[i].ID;
        revidTemp = fData.brDList[i].revID;
        fromTemp = fData.brDList[i].From;
        toTemp = fData.brDList[i].To;
        name = fData.brDList[i].CKT;
        pTemp = p0_base[i]*fData.baseMVA;
        qTemp = q0_base[i]*fData.baseMVA;
        pRevTemp = 0 #p0_base[revidTemp]*fData.baseMVA;
        qRevTemp = 0 #getvalue(mp[:q0][revidTemp])*fData.baseMVA;
        write(sol2file,"0,$name,$fromTemp,$toTemp,i_$(fromTemp)_$(toTemp)_$(name),$pTemp,$qTemp,$pRevTemp,$qRevTemp \n");
end

for line in sol2branch
	write(sol2file, line);
end


write(sol2file,"--end of line flow \n")

close(sol2file)   

end		#end build function, commented for offline testing
