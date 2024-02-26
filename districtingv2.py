import homodistrictingv2 as d
#import laparcgisx6btested as d
import time,sys,random
#parameters: filename1 filename2 number_of_facilities pop_size timelimit spp seed max_loops_solution_not_improved
feas=0
first=0
seed=-1
psize=10 #population size
maxloops=1000 #max_loops_solution_not_improved
t=250 #timelimit for searching, SPP modeling
spp=0 #set partitioning
n=0 # number of facilities, number of service areas
if len(sys.argv) >= 2:
    fn=sys.argv[1]
if len(sys.argv) >= 3:
    fn2=sys.argv[2]
if len(sys.argv) >= 4:
    feas=int(sys.argv[3])
if len(sys.argv) >= 5:
    first=int(sys.argv[4])
if len(sys.argv) >= 6:
    n=int(sys.argv[5])
if len(sys.argv) >= 7:
    psize=int(sys.argv[6])
if len(sys.argv) >= 8:
    maxloops=int(sys.argv[7])
if seed<0:
    seed=random.randint(0,100)
random.seed(seed)
d.seed=seed

def save_sol_location_problem():
    fn2=fn+"_sol"+str(n)+"_"+str(int(d.biobjective))+".txt"
    f = open(fn2,"w")
    idx=0
    f.write("Idx\tRID\tCenter\n")
    for i in range(d.num_units):
        k=d.node_groups[i]
        s=str(i)+"\t"+str(k)+'\t'+str(d.centersID[k])+"\n"
        f.write(s)
    f.close()


def save_homoDP_sol():
    fn2=fn+"_sol"+str(n)+"_"+str(int(d.biobjective))+".txt"
    f = open(fn2,"w")
    idx=0
    f.write("FID\tX\tY\txidx\tyidx\tRID\n")
    for i in range(d.num_units):
        for x in d.nodes[i][:5]: 
            f.write(str(x)+"\t")
        f.write(str(d.node_groups[i])+"\n")
    f.close()

#problem definition				GAP		SAP		SSCFLP	SSCKFLP	CFLSAP	CKFLSAP	pdp		PMP		CPMP
#d.location_problem 			0		0		1		1		1		1		2		3		1
#d.fixed_cost_obj				0		0		1		1		1		1		0		0		0
#d.spatial_contiguity			0		1		0		0		1		1		1		0/1/2	0/1
#d.adaptive_number_of_facilities0		0		1		0		1		0		0		0		0
#d.pop_dis_coeff=				10000	10000	10000	10000	10000	10000	10000	10000	10000 penalty on facility overload/district balance
#d.pop_deviation=				0.0		0.0		0.0		0.0		0.0		0.0		0.05	0.0		0.0

d.location_problem=3
d.all_units_as_candadate_locations=0
d.fixed_cost_obj=0
d.spatial_contiguity=02
#0 non-connected, 1 strictly-connected, 2 connected with several parts 
d.spatial_contiguity_minimum_percentage=10.0
d.adaptive_number_of_facilities=0
d.pop_dis_coeff=100000#10000#1000000 #100000 for FLP 1-5 for PDP
d.pop_deviation=0.00 #for PDP
d.max_loops_solution_not_improved=maxloops #for search termination

t0=time.time()
#read instance file(s)
#d.read_pmp_climate_instance(fn,fn2,0,0) #PMP benchmark
d.read_pmp_climate_instance(fn,fn2,feas,first) #PMP benchmark

for i in range(d.num_units):
    if d.node_neighbors[i]==[]:
        print "node", i, "has no any neighbor"
#solution method			GAP		SAP		SSCFLP	SSCKFLP	FLSAP	CKFLSAP	pdp
#d.initial_solution_method=	0-construction,1-LP (solver needed), 8 sampling, 9 LR (all except PDP)
#d.mip_solver="cplex"		"gurobi","cplex", "cbc", ""
#d.multi_start_count	solver parameter: multistart or population size
#d.is_spp_modeling		solver parameter: set partitioning 
#d.heuristic_time_limit	solver parameter: time in second
#d.operators_selected=[0,1]		assignment operators: 0 one-unit move, 1 two-unit move, 2 three-unit move
#d.solution_similarity_limit=5 	5~15

d.initial_solution_method=0
    #0 construction method, region growth, k-medoids, weighted k-medoids
    #1 MIP model with relaxation, TP for GAP/SAP, CFLP for SSCFLP/PDP, greedy for pmp  
    #2 w_k_medoids for edp, PDP;
    #3 random for PMP
    #8 sampling for PMP; x,y relaxation
    #9 Lagrangian relaxation, for GAP/SAP/SSCFLP/CKFLP/PDP
d.location_operators_selected=[]#[0,4,4,4,4] #[9]#
    #0 swap
    #1 drop, 2 add  (DO NOT used for CKFLP/FLSAP/CKFLSAP!!!)
    #4 multi-exchange
    #5 k-location-exchange (PMP only), 7 TB (PMP only), 
    #8 PMP submpdel (PMP only,not for CPMP) 
    #9 CFLP submodel (for CKFLP/FLSAP/CKFLSAP)
d.assignment_operators_selected=[0] 
    #0 for 1-unit move, 1 for 2-unit move, 2 for 3-unit move(slow!!!) 
d.ruin_oprators=[3,4]#[0,1,2,3,4,5,6]#[0,1,2,3]#[0,1,2,3,4,5,5,5,5] #[0,1,3,4,6]
    #0 ruin a region around a location 
    #1 ruin a few districts randomly
    #2 ruin a few connected districts 
    #3 ruin some boundary units
    #4 move a few boundary units
    #6: path_relink  (GAP,SAP only)
    #7: assign with tp model 
    #8: assign with gap model
    #9: mode a few locations(for PMP, PDP)

d.solution_similarity_limit=5.0 #max(10.0,100.0/n)
d.is_spp_modeling=01
d.acceptanceRule="ils"
d.heuristic_time_limit=t 
d.multi_start_count=psize  #population_size
d.ruin_percentage=5.0
d.max_loops_solution_not_improved=maxloops #used for search termination
#d.location_cover_check()
d.seed=seed
d.solver_message=01

d.ils_large_homoDP(n,psize)

#
algorithm=3000
if d.location_problem==3: #pmp, spatial clustering
    if algorithm==0: d.ils_pmp(n,psize,t,seed) #good
    if algorithm==1: d.mip(n,d.location_problem,d.spatial_contiguity,t)
    if algorithm==2: d.multi_exchange(n,psize,t,seed)	#heuristic
    if algorithm==3: d.ils_large_homoDP(n,psize,t,seed)  #sampling for large instances,#best?
    if algorithm==4: d.location_sub_mip(n,psize,t,seed) 
#d.ils(n,psize,t,spp,seed)
#d.ils_lr_pdp(n,psize,t,spp,seed)

print "=========================Final results========================="
print "objective:",d.biobjective
print "facility cost",d.objective_fcost
print "transportation cost:",d.objective
print "srrvice overload", d.objective_overload
print "pool size",len(d.region_pool)
print "total time",time.time()-t0
print "facilities selected",[x for x in d.centersID if x>=0]
#print "demand assignment:", d.node_groups
#print "service area stat:"
#for i in range(d.num_districts):
#    if d.district_info[i][0]==0: continue
#    print d.facilityCandidate[i],d.district_info[i], (d.district_info[i][2]+d.facilityCost[i])/d.district_info[i][1],
#    print "continuality?",d.check_continuality_feasibility(d.node_groups,i)
print "solution",
for x in d.node_groups: print x,
print 
save_homoDP_sol()

