# -*- coding: utf-8 -*-
## An unified algorithm framework for the Facility Location, Allocation and Service Area Problems. 
## yfkong@henu.edu.cn, Apr.,2021
## to solve problems such as:
## 1 GAP,SAP
## 2 SSCFLP, SSCKFLP, CFLSAP, CKFLSAP
## 3 PMP/CPMP
## 4 PDP
## 5 Eaqul-capacity PMP (ECPMP)
## 6 SSCFLP_R:CFLP with covering radius and covering percentage

import sys,os,random,time,copy,math,tempfile
#ArcGIS
has_arcpy=0
try:
    import arcpy
    has_arcpy=1
except:
    has_arcpy=0
#mip solver
mip_solvers=[] #MIP solvers supported 
mip_solver=''  #MIP solver, "cplex", "cbc" or ""
mip_file_path=tempfile.gettempdir()
#os.chdir(mip_file_path)  #used in arcgis
try:
    import cplex
    #mip_solvers.append('cplex')
except: 
    pass
try:
    import pulp
    s=pulp.apis.GUROBI_CMD()
    if s.available()==False:
        pass
    else:
        mip_solvers.append('gurobi')
    s=pulp.apis.cplex_api.CPLEX_CMD()
    if s.available()==False:
        pass
    else:
        mip_solvers.append('cplex')
    s=pulp.apis.COIN_CMD()
    if s.available()==False:
        pass
    else:
        mip_solvers.append('cbc')
except: 
    pass
if len(mip_solvers)>0: mip_solver=mip_solvers[0]

#constant
MAXNUMBER=1.0e+20
MINNUMBER=1.0e-10
#instance info
nodes=[]
nodes_std=[] #for homoDP only
weight_features=[] #for homoDP only
num_units=-1
nodedij=[]
nodedik=[]  #weighted cost from i to k, =nodedij*nodes[][3] 
nodendij=[] #network distance
node_neighbors=[]
facility_neighbors=[]
total_pop=0
avg_pop=0
total_supply=0
all_units_as_candadate_locations=0
facilityCandidate=[]
facilityCapacity=[]
facilityCost=[]
num_facilityCandidate=-1
num_districts=-1 # number of service areas/facilities
avg_dis_min=0.0
potential_facilities=[]
NearFacilityList=[]
nearCustomer=[]
nearCustomers=[]
geo_instance=1
pmp_I_eaquls_J=1

#parameters for districting
location_problem=1
max_num_facility=999
adaptive_number_of_facilities=1
fixed_cost_obj=1
spatial_contiguity=1 # 0 no, 1 yes, 2 yes with multi-parts
spatial_contiguity_minimum_percentage=5
pop_dis_coeff=10000.0 #used in the objective function
pop_deviation=0.00 #for pdp, 5%

#current solution
centersID=[]
node_groups=[]
district_info=[] #[[0,0,0.0] for x in range(num_districts)] # solution
objective_overload=0
obj_balance=MAXNUMBER
objective=MAXNUMBER
objective_fcost=MAXNUMBER
biobjective=MAXNUMBER
objective_supply=0.0

given_solution=0 #reserved
all_solutions=[]

#best solution in each start
best_solution =[] # node_groups[:]
best_centersID=[]
best_biobjective=MAXNUMBER
best_objective=MAXNUMBER
best_objective_overload = MAXNUMBER
best_objective_fcost = MAXNUMBER
#global best solution 
#best_centers_global=[]
best_solution_global=[]
best_centersID_global=[]
best_biobjective_global = MAXNUMBER
best_objective_global = MAXNUMBER
best_objective_fcost_global = MAXNUMBER
best_overload_global = MAXNUMBER

#search statistics
time_check=0
time_check_edge_unit=0
time_spp=0.0
time_update_centers=0.0
time_op=[0.0 for x in range(10)]
time_ruin_recreate=[0.0 for x in range(10)]
time_location=[0.0 for x in range(10)]
time_pmp_re_location=0.0
time_Whitaker=0.0
time_repair=0
count_op=[0.0 for x in range(10)]
check_count=0
improved=0
move_count=0

#search histry
region_pool = []
pool_index=[]


#local search
acceptanceRule="hc" #solver name
assignment_operators_selected=[0,1] #0 one-unit move, 1 two-unit move, 2 three-unit move
location_operators_selected=[0,1,2,3,4] #0 swap, 1 drop, 2 add, 3 add+drop, 4 me
ruin_oprators=[0,1,2,3,4] #ruin0, ruin1, 9 mip assign
multi_start_count=6 #population size for GA, ILS, VNS, LIS+VND
initial_solution_method=0 #0 construction, 1 LP
assign_method=0 #not used
assign_or_Location_search_method=0
large_facility_cost=0
maxloops=1000
max_loops_solution_not_improved=-1
SA_maxloops = 100 # maximum number of search loops for GA
SA_temperature=1.0
op_random = 1 # operators used sequentially (0) or randomly(1)
last_loop_obj=0.0
ruin_percentage=20
mainloop=0
mutation_rate=0.01 
cross_methods=-1
adj_pop_loops=5000
solution_similarity_limit=10.0

#mip modelling for inititail solution, spp and full model
is_spp_modeling=1 #0 no spp modelling, 1 modelling at the end, 2 modelling in case of local optimum
linear_relaxation=0
spp_loops=400
solver_time_limit=7200 #used for mip modeling
solver_mipgap=0.000000000001
solver_message=0
heuristic_time_limit=300
seed =random.randint(0,10000)
random.seed(seed)
locTabuLength=100  #nf*psize
locTabuList=[]
locTabuList2=[]

def update_locTabuList(locs):
    global locTabuList
    if locs not in locTabuList:
        locTabuList.append(locs)
    if len(locTabuList) >locTabuLength:
        del locTabuList[0]

def arcpy_print(s):
    if has_arcpy==1: 
        arcpy.AddMessage(s)
    else:
        print s

#record a region in current solution
def update_region_pool(rid):
    global pool_index
    global time_spp
    global region_pool
    t=time.time()
    if is_spp_modeling<=0: return
    if centersID[rid]<0: return
    #if spatial_contiguity==1 and check_continuality_feasibility(node_groups,rid)<1:
    #    return
    ulist=[x for x in  range(num_units) if node_groups[x]==rid]
    if ulist==[]:
        #print "empty area:",rid,node_groups
        return
    cost1=district_info[rid][2]
    cost2=sum(nodedik[x][rid] for x in ulist)
    if abs(cost1-cost2)>0.001: print rid,cost1,cost2
    obj=district_info[rid][2]+district_info[rid][4]*pop_dis_coeff
    idx=int(obj*100000)
    idx+=sum(x*x for x in ulist)
    if idx not in pool_index[rid]:
        pool_index[rid].append(idx)
        region_pool.append([ulist,district_info[rid][2],district_info[rid][1],district_info[rid][4],rid])
    time_spp+=time.time()-t
    return

#record all regions in current solution
def update_region_pool_all():
    if is_spp_modeling<=0:
        return
    for rid in range (num_districts):
        if centersID[rid]<0: continue
        update_region_pool(rid)
def delete_region_pool_all():
    global region_pool
    global pool_index
    region_pool=[]
    for i in range(num_districts):
        pool_index[i]=[]

#check continuality of a solution (sol)
def check_solution_continuality_feasibility(sol):
    if spatial_contiguity==0: return 9
    feasible = spatial_contiguity
    for i in range (num_districts):
        if centersID[i]<0: continue
        if check_continuality_feasibility(sol,i) <= 0:
            feasible=0  #infeas.
            break
    return feasible

def check_current_solution_continuality_feasibility():
    feaslist=[]
    for i in range (num_districts):
        if centersID[i]<0: continue
        sta=check_continuality_feasibility(node_groups,i)
        feaslist.append(sta)
    return feaslist

#check continuality of a region (rid) in solution (sol)
def check_continuality_feasibility(sol,rid):
    global time_check
    global check_count
    if spatial_contiguity==0: return 9
    #if geo_instance==0: return -1
    u=facilityCandidate[rid]
    check_count+=1
    t=time.time()
    ulist1=[x for x in range(num_units) if sol[x]==rid and x!=u]
    ulist2=[u]
    #ulist2.append(ulist1.pop())
    for x in ulist2:
        for i in range(len(ulist1)):
            j=ulist1[i]
            if j in node_neighbors[x]:
                ulist2.append(j)
                ulist1[i]=-1
        ulist1=[x for x in ulist1 if x>=0]
    if len(ulist1)==0:          
        time_check+=time.time()-t
        return 1  #feasible
    if spatial_contiguity==2:
        ulist=[x for x in range(num_units) if sol[x]==rid]
        flist=frag_unit_minority(ulist)
        if flist==[]: 
            time_check+=time.time()-t
            return 2
    time_check+=time.time()-t
    return 0    #infeasible

#check continuality of a list of units
def check_ulist_continuality(ulist):
    if spatial_contiguity==0: return 1
    global time_check
    global check_count
    t=time.time()
    ulist1=ulist[:]
    ulist2=[]
    ulist2.append(ulist1.pop())
    check_count+=1
    for x in ulist2:
        for i in range(len(ulist1)):
            #if ulist1[i]==-1: continue
            if ulist1[i] in node_neighbors[x]:
                ulist2.append(ulist1[i])
                ulist1[i]=-1
        ulist1=[x for x in ulist1 if x>=0]         
    #ulist3=[x for x in ulist1 if x!=-1]
    if len(ulist1)==0:          
        time_check+=time.time()-t
        return 1  #feasible
    if spatial_contiguity==2:
        flist=frag_unit_minority(ulist)
        time_check+=time.time()-t
        if flist==[]: return 2
    return 0    #infeasible

#return a list of boundary units
def find_edge_units():
    if spatial_contiguity==0 and random.random()>0.5:
        return find_tail_units()
    ulist=[]
    for x in range(num_units):
        if geo_instance==1 and spatial_contiguity==1 and x in centersID:
            continue  #bug for benckmark instances
        k=node_groups[x]
        for y in node_neighbors[x]:
            if node_groups[y] != k:
                ulist.append(x)
                break
    random.shuffle(ulist)
    if objective_overload==0: 
        return ulist
    ulist=[[x, district_info[node_groups[x]][4]] for x in ulist]
    ulist.sort(key=lambda x:-x[1])
    ulist=[ x[0] for x in ulist]
    return ulist

def find_tail_units():
    #p=1.0
    #nf=sum(1 for x in centersID if x>=0)
    #nu=num_units*1.0/nf
    #if nu>7:
    #    p=(4.0*math.sqrt(nu)-4)/nu
    tlist=[]
    for k in range(num_districts):
        if centersID[k]<0: continue
        ulist=[[x,nodedij[x][k]] for x in range(num_units) if node_groups[x]==k]
        ulist.sort(key=lambda x:-x[1])
        #n=int(len(ulist)*p)
        n=len(ulist)
        if n>7: n=int(4*math.sqrt(n))-4
        if n<1: n=1
        ulist=ulist[:n]
        tlist+=[x[0] for x in ulist]
    random.shuffle(tlist)
    #return tlist 
    if objective_overload==0: #!!!
        return tlist
    #ulist=[ [x, district_info[node_groups[x]][4]] for x in tlist]
    ulist=[[x,nodedij[x][node_groups[x]]*(1+district_info[node_groups[x]][4]*10)] for x in tlist]
    ulist.sort(key=lambda x:-x[1])
    ulist=[ x[0] for x in ulist]
    #print len(tlist)
    return ulist

#return a list of edge units that having three or more neighor regions
def find_edge_units_2():
    ulist=[]
    for x in range(num_units):
        if spatial_contiguity==1 and x in centersID: continue #bug for benckmark instances
        rlist=[node_groups[y] for y in node_neighbors[x]]
        rlist.append(node_groups[x])
        rset=set(rlist)
        if len(rset)>2:
            ulist.append(x)
    random.shuffle(ulist)
    return ulist


#check an edge unit (uid), reserved
def is_edge_unit(uid):
    #global time_check_edge_unit
    if spatial_contiguity==0: return True
    t=time.time()
    rlist = [node_groups[x] for x in node_neighbors[uid]]
    rlist+=node_groups[uid]
    #time_check_edge_unit+=time.time()-t
    if len(set(rlist))==1:
        return False
    return True

#update region information of the current solution
def update_district_info():
    global objective_overload
    global objective
    global biobjective
    global objective_fcost
    global district_info
    global move_count
    global obj_balance
    global centersID
    global objective_supply
    global avg_dis_min
    for k in range(num_districts):
        district_info[k][0] = 0
        district_info[k][1] = 0.0
        district_info[k][2] = 0.0
        district_info[k][3] = 0.0
        district_info[k][4] = 0.0
    for k in range(num_districts):
        if centersID[k]<0 and k in node_groups:
            arcpy_print("debug: a facility not selected but used: " + str(k))
            centersID[k]=facilityCandidate[k]
    for k in range(num_districts):
        if centersID[k]<0:
            continue
        ulist=[x for x in range(num_units) if node_groups[x]==k]
        if len(ulist)==0:
            if location_problem==3: continue
            if adaptive_number_of_facilities==1:
                supply=sum(facilityCapacity[x] for x in range(num_districts) if centersID[x]>=0)
                centersID[k]=-1
                continue
        district_info[k][0] = len(ulist)
        district_info[k][1] = sum(nodes[x][3] for x in ulist)
        district_info[k][2] = sum(nodedik[x][k] for x in ulist)
        district_info[k][3] = facilityCapacity[k] 
        district_info[k][4] = max(0.0,district_info[k][1]-facilityCapacity[k]) # -district_info[k][3]
        if location_problem==3: district_info[k][4]=0 #pmp
        if location_problem==2:
            bal=0.0
            dev=pop_deviation*total_pop/max_num_facility
            if district_info[k][1]>district_info[k][3]+dev: bal=district_info[k][1]-district_info[k][3]-dev
            if district_info[k][1]<district_info[k][3]-dev: bal=district_info[k][3]-district_info[k][1]-dev
            district_info[k][4]=bal
        #print centersID,node_groups
    bal=sum(x[4] for x in district_info)
    objective=sum([x[2] for x in district_info])
    objective_overload=bal
    #if objective/total_pop<avg_dis_min:
    #    avg_dis_min=objective/total_pop
    avg_dis_min=objective/total_pop
    biobjective=objective+objective_overload*avg_dis_min*pop_dis_coeff


    objective_supply=sum(facilityCapacity[x] for x in range(num_districts) if centersID[x] >=0)
    #biobjective=objective+objective_overload*avg_dis_min*1000000
    #biobjective=bal2*avg_dis_min*1000000
    if fixed_cost_obj==1:
        fcost=sum(facilityCost[x] for x in range(num_districts) if centersID[x] >=0)
        objective_fcost=fcost
        biobjective+=fcost
    move_count+=1

def find_frag_unit():
    if spatial_contiguity==0: return []
    global time_check_edge_unit
    t=time.time()    
    frag_units=[]
    if spatial_contiguity==2:
        for k in range(num_districts):
            if centersID[k]==-1: continue
            ulist=[x for x in range(num_units) if node_groups[x]==k]
            nflist=frag_unit_minority(ulist)
            frag_units+=nflist
    if spatial_contiguity!=2:
        for k in range(num_districts):
            if centersID[k]==-1: continue
            ulist2=[centersID[k]]
            ulist1=[x for x in range(num_units) if node_groups[x]==k and x!=centersID[k]]
            for x in ulist2:
                for i in range(len(ulist1)):
                    if ulist1[i]==-1: continue
                    if ulist1[i] in node_neighbors[x]:
                        ulist2.append(ulist1[i])
                        ulist1[i]=-1
                ulist1=[x for x in ulist1 if x>=0]
            frag_units+=ulist1

    random.shuffle(frag_units)
    time_check_edge_unit+=time.time()-t
    #print frag_units
    return frag_units    

def frag_unit_minority(ulist):
    final_list=[]
    ulist2=ulist[:1]
    ulist1=ulist[1:]
    total_area=sum(x[3] for x in nodes)
    while 1:
        for x in ulist2:
            for i in range(len(ulist1)):
                if ulist1[i]==-1: continue
                if ulist1[i] in node_neighbors[x]:
                    ulist2.append(ulist1[i])
                    ulist1[i]=-1
            ulist1=[x for x in ulist1 if x>=0]
        final_list.append([len(ulist2),ulist2[:]])
        if len(ulist1)<=1:
           if len(ulist1)==1:
                final_list.append([1,ulist1[:]])
           break
        u=ulist1[0]
        ulist2=[u]
        del ulist1[0]
    if len(final_list)==1: return []
    final_list.sort(key=lambda x:x[0])
    #del final_list[-1]
    flist=[]
    n=total_area*spatial_contiguity_minimum_percentage/max_num_facility/100
    for x in final_list: 
        area=sum(nodes[i][3] for i in x[1])
        if area>n:continue
        flist+=x[1]
    #print [len(flist),len(ulist)],
    return flist

def district_with_multiparts(ulist):
    final_list=[]
    ulist2=ulist[:1]
    ulist1=ulist[1:]
    while 1:
        for x in ulist2:
            for i in range(len(ulist1)):
                if ulist1[i]==-1: continue
                if ulist1[i] in node_neighbors[x]:
                    ulist2.append(ulist1[i])
                    ulist1[i]=-1
            ulist1=[x for x in ulist1 if x>=0]
        final_list.append([len(ulist2),ulist2[:]])
        if len(ulist1)<=1:
           if len(ulist1)==1:
                final_list.append([1,ulist1[:]])
           break
        u=ulist1[0]
        ulist2=[u]
        del ulist1[0]
    final_list.sort(key=lambda x:-x[0])
    #del final_list[-1]
    flist=[x[0] for x in final_list]
    return flist

def repair_fragmented_solution_large_instance():
    global node_groups
    global centersID
    global time_repair
    if spatial_contiguity==0: return
    t=time.time()
    frag_units=find_frag_unit()
    #print "frag_units",len(frag_units),
    if len(frag_units)==0: return
    #sol=node_groups[:]
    for x in frag_units:
        node_groups[x]=-1
    update_district_info()
    #print len(frag_units),
    while len(frag_units)>0:
        cands=[]
        for x in frag_units:
            for y in node_neighbors[x]:
                k=node_groups[y]
                if k<0: continue
                cands.append([x,k,nodedik[x][k]])
        if cands==[]: break
        cands.sort(key=lambda x:x[2])
        n=len(cands)/10
        if n<1: n=1
        for x in cands[:n]:
            nid=x[0]
            if node_groups[nid]>=0: continue
            node_groups[nid]=x[1]
            if nid in frag_units: frag_units.remove(x[0])
    for x in frag_units:
        for k in NearFacilityList[x]:
            if centersID[k]>=0:
                 node_groups[x]=k
                 break
    #print len(frag_units)
    #print int(objective),
    update_district_info()
    #print int(objective),
    time_repair+=time.time()-t
    repair_fragmented_solution()
    #print objective,
	
#repair the fragmented solution
def repair_fragmented_solution():
    if spatial_contiguity==0: return
    #if num_units/max_num_facility>100:
    #    repair_fragmented_solution_large_instance()
    #    return
    global node_groups
    global centersID
    global time_repair
    t=time.time()
    for k in range(num_districts):
        if centersID[k]<0: continue
        u=nearCustomer[k]
        node_groups[u]=k
    update_district_info()
    frag_units=find_frag_unit()
    #print "frag_units",frag_units,
    if len(frag_units)==0: return
    sol=node_groups[:]
    for x in frag_units:
        node_groups[x]=-1
    # if location_problem>=3:
        # for k in range(num_districts):
            # if centersID[k]<0: continue
            # if node_groups[k]>=0: continue
            # c,ulist=update_center(k)
            # centersID[k]=-1
            # centersID[c]=c
            # for x in ulist: node_groups[x]=c
    update_district_info()
    #print len(frag_units),
    while len(frag_units)>0:
        newk=-1
        nid=-1
        cost=MAXNUMBER
        for x in frag_units:
            for y in node_neighbors[x]:
                k=node_groups[y]
                if k>=0:
                    gap=max(district_info[k][1]+ nodes[x][3] - facilityCapacity[k],0)
                    if location_problem==3: gap=0
                    cost2=gap*avg_dis_min*pop_dis_coeff + nodedik[x][k]
                    if cost2<cost:
                        nid=x
                        newk=k
                        cost=cost2
        if newk>=0:
            node_groups[nid]=newk
            update_district_info()
            frag_units.remove(nid)
        else:
            break
    #print "frag_units", frag_units
    for x in frag_units:
        node_groups[x]=sol[x]
    #print len(frag_units)
    #print int(objective),
    update_district_info()
    #print check_current_solution_continuality_feasibility()
    time_repair+=time.time()-t
    #print int(objective),


def assign_ruin_recreate(idx):
    global time_ruin_recreate
    #if location_problem==4: return -1,1
    t=time.time()
    if idx<0 and len(ruin_oprators)<1: return -1
    ruin_idx=idx
    if idx<0:
        r=random.randint(0,len(ruin_oprators)-1)
        ruin_idx=ruin_oprators[r]
    #ruin_idx=0~4, diversification
    if ruin_idx==0: r_r_perb_location()
    if ruin_idx==1: r_r_perb_district()
    if ruin_idx==2: r_r_large_region()
    if ruin_idx==3: r_r_perb_edge_units()
    if ruin_idx==4: r_r_perb_mutation()
    if ruin_idx==9: #move a few locations for PMP, PDP
        r_r_perb_center_locations()
    if spatial_contiguity>=1: 
        if location_problem==3 and num_units>2000:
            repair_fragmented_solution_large_instance()
        else:
            repair_fragmented_solution()
    time_ruin_recreate[ruin_idx]+=time.time()-t
    return ruin_idx,1

#r&r method
#remove the edge units and reassign them to nearest facilities
def r_r_perb_edge_units():
    global node_groups
    ulist=find_edge_units() 
    num=int(len(ulist)*ruin_percentage/100)
    if num<1: num=1
    ulist=ulist[:num]
    for x in ulist:
        for k in NearFacilityList[x]:
            if centersID[k]>=0:
                node_groups[x]=k
                break
    update_district_info()
    if spatial_contiguity>=1: repair_fragmented_solution()

def r_r_large_region():
    global node_groups
    dlist=select_region(-1)
    nf=sum(1 for x in range(num_districts) if centersID[x]>=0)
    num=int(nf*ruin_percentage/100)
    num=max(num,3)
    num=min(num,5)
    dlist=dlist[:num]
    #if len(dlist)>3: dlist=dlist[:3]
    ulist=[x for x in range(num_units) if node_groups[x] in dlist]
    for x in ulist:
        for k in NearFacilityList[x]:
            if centersID[k]>=0:
                node_groups[x]=k
                break
    update_district_info()
    if spatial_contiguity>=1: repair_fragmented_solution()


#r&r method
#remove 1/40 units arround an edge unit
def r_r_perb_location():
    global node_groups
    nf =sum(1 for x in centersID if x>=0)
    ulist=find_edge_units_2()
    if len(ulist)==0:
        ulist=find_edge_units()
    if len(ulist)==0:
        return 
    ulist=[ulist[0]]
    pop=0
    num=int(num_units*ruin_percentage/100)
    #if num<2: num=2
    for x in ulist:
        for y in node_neighbors[x]:
            if y not in ulist:
                ulist.append(y)
        if len(ulist)>=num:
            break
    for x in ulist:
        for k in NearFacilityList[x]:
            if centersID[k]>=0:
                node_groups[x]=k
                break
    update_district_info()
    if spatial_contiguity>=1: repair_fragmented_solution()

def r_r_perb_mutation():
    global node_groups
    rate=ruin_percentage
    ulist=find_edge_units()
    nf=sum(1 for x in range(num_districts) if centersID[x]>=0)
    #num=int(num_units*rate)
    num=nf/2
    #if location_problem==3: num=max(5,nf/4)

    if num<1: num=1
    for x in ulist:
        if node_groups[x]==x and location_problem==3: continue
        for y in node_neighbors[x]:
            if node_groups[x]!=node_groups[y]:
                node_groups[x]=node_groups[y]
                num-=1
                break
        if num<=0:
            break
    update_district_info()
    #if spatial_contiguity>=1: repair_fragmented_solution()


def select_region(seed):
    nf=sum(1 for x in range(num_districts) if centersID[x]>=0)
    n=100*nf/num_units  #areas with 100 units
    if nf<=5: n=nf
    #if nf>=7: n=random.randint(nf/2+1,nf)
    if nf>=6 and nf<=11: n=random.randint(nf/2+1,9)
    if nf>=12 and nf<=16: 
        n=random.randint(7,10)
    if nf>=17: 
        n=random.randint(7,10)
    if n*num_units/nf<80: 
        n=min(10,80*nf/num_units)
    if location_problem==3: n=min(128,num_districts/10)
    #clist=[]
    #u=random.randint(0,num_units-1)
    #if r>=0: 
    #    u=nearCustomer[r]
    #    clist=[r]
    #for k in NearFacilityList[u]:
    #    if centersID[k]<0: continue
    #    if k not in clist: clist.append(k)
    #    if len(clist)==n: break
    #return clist
    #if objective_overload>0: ???
    u=seed
    if u<0: u=random.randint(0,num_units-1)
    r=node_groups[u]
    if location_problem==0 and objective_overload>0: #SAP
        rlist=[k for k in range(num_districts) if district_info[k][4]>0]
        random.shuffle(rlist)
        r=rlist[0]
        u=nearCustomer[r]
        return NearFacilityList[u][:n]

    clist=[r]
    if random.random()>-0.5:
        for k in NearFacilityList[u]:
            if centersID[k]<0: continue
            if k==r: continue
            clist.append(k)
            if len(clist)==n: break
        #clist.sort()
        return clist

    for i in facility_neighbors[r]:
        if centersID[i]<0: continue
        clist.append(i)
        if len(clist)>=n: break
    #clist.sort()
    return clist

#update the best and the global best solution
#if the current solution is better than the best
def update_best_solution():
    global best_solution
    global best_centersID
    global best_biobjective
    global best_objective
    global best_objective_fcost
    global best_overload
    global best_objective_overload
    global best_centersID_global
    global best_solution_global
    global best_biobjective_global
    global best_objective_global
    global best_objective_fcost_global
    global best_overload_global    
    global improved_loop
    global improved
    global avg_dis_min
    improve =0
    if location_problem==1 and adaptive_number_of_facilities==0:
        nf=sum(1 for x in centersID if x>=0)
        if nf!=max_num_facility: return 0
    #if spatial_contiguity==1 and check_solution_continuality_feasibility(node_groups)==0:
    #    ##noprint "check_solution_continuality_feasibility!!!"
    #    return improve
    biobj=biobjective
    biobj_best=best_biobjective
    if biobj<=biobj_best:
        best_biobjective=biobj
        best_objective = objective
        best_objective_fcost=objective_fcost
        best_objective_overload = objective_overload
        best_solution = node_groups[:]
        best_centersID=centersID[:]
        improved_loop=mainloop
        improve =1
        improved+=1
    if biobj<best_biobjective_global:
        best_biobjective_global=biobj
        best_centersID_global=centersID[:]
        best_overload_global = objective_overload
        best_solution_global =node_groups[:]
        best_objective_global = objective
        best_objective_fcost_global=objective_fcost
        avg_dis_min=biobj/total_pop
    return improve

def pmp_VND_search():
    global node_groups
    global time_op
    global count_op
    improved = 0
    t=time.time()
    if spatial_contiguity==0: return 0
    while 1:
        improve = 0
        #if spatial_contiguity==0:
        ulist=find_edge_units()
        for u in ulist:
            k=node_groups[u]
            for u2 in node_neighbors[u]:
                if node_groups[u2]==k: continue
                k2=node_groups[u2]
                if nodedik[u][k2]>=nodedik[u][k]: continue
                #if spatial_contiguity>=1: 
                #    ulist2=[x for x in range(num_units) if node_groups[x]==k and x!=u]
                #    if ulist2==[]: continue
                #    if check_ulist_continuality(ulist2)==0: continue
                node_groups[u] = k2
                improve+=1
                improved+=1
        if improve == 0: break
    time_op[0]+=time.time()-t
    #print "improved",improved,
    update_district_info()
    return improved
def RRT_local_search():
    #if location_problem==3:
    #    pmp_VND_search()
    #    return
    global improved
    #global node_neighbors
    improved=0
    #operators=[0,1,2,3,4]
    operators=assignment_operators_selected[:]
    #if op_random == 1 and random.random()>0.5:
    #    random.shuffle(operators)
    for op in operators:
        if op == 0:
            one_unit_move()
            #update_region_pool_all()
        if op == 1:
            two_unit_move()
            #update_region_pool_all()
        if op == 2:
            #if random.random()<0.1: three_unit_move()
            three_unit_move()
            #update_region_pool_all()
        if op == 3:
            three_unit_move_simple()
    return


def read_pmp_climate_instance(f1,f2,num_feat,num_first): #climate data
	global num_units
	global nodes
	global nodes_std
	global node_neighbors
	global nodedij
	global nodedik
	global facilityCandidate
	global facilityCapacity
	global facilityCost
	global num_facilityCandidate
	global num_districts
	global district_info
	global potential_facilities
	global geo_instance
	global total_pop
	global weight_features
	geo_instance=1
	f = open(f1)
	line = f.readline()
	line=line[:-1]
	items=line.split("\t")
	num_units=int(items[0])
	num_features=int(items[1])
	num_first_feature=5
	if len(items)>=3:
		num_first_feature=int(items[2])
	data_standardlize=1
	if len(items)>=4:
		data_standardlize=int(items[3])
	if num_feat>=1:
		num_features=num_feat
		num_first_feature=num_first
	print "instance size, features, first columme and std.",num_units,num_features,num_first_feature,data_standardlize

	weight_features=[1 for x in range(num_features)]
	line = f.readline()
	line=line[:-1]
	items=line.split("\t")
	for i in range(num_features):
		weight_features[i]=float(items[num_first_feature+i])
	#print weight_features[:10],"..."
	line = f.readline()
	line=line[:-1]
	idx=0
	for i in range(num_units):
		if i%100==0: print ".",
		line = f.readline()
		line=line[:-1]
		items = line.split('\t')
		node=[int(items[0]),float(items[1]),float(items[2]),float(items[3])]
		#items[:num_first_feature]
		for j in range(len(items)):
			if items[j]=="": items[j]="0"
		if data_standardlize<=1:
			attr=[float( items[x+num_first_feature] ) for x in range(num_features)]
		else:
			attr=[items[x+num_first_feature] for x in range(num_features)]
		nodes.append(node+attr)
	f.close()
	num_districts=num_units
	facilityCandidate=range(num_districts)
	facilityCapacity=[num_units for x in range(num_districts)]
	print "num_units, features",num_units, num_features
	print "weights",weight_features
	#print "nodes", nodes[0]
	#print "nodes", nodes[num_units-1]
	nodes_std=[[0.0 for x in range(num_features)] for y in range(num_units)]
	if data_standardlize==0:
		for k in range(num_features):
			if weight_features[k]<0.0000001: continue
			for i in range(num_units):
				nodes_std[i][k]=nodes[i][k+4]*math.sqrt(weight_features[k])
	if data_standardlize==1:
		print "data standardlizing..."
		area_total=sum(x[3] for x in nodes)
		for k in range(num_features):
			if weight_features[k]<0.0000001: continue
			attr=[x[3]*x[k+4] for x in nodes]
			avg=sum(attr)/area_total
			var=sum((x-avg)*(x-avg) for x in attr)
			if var<0.0000001: continue
			var/=num_units
			dev=pow(var,0.5)
			for i in range(num_units):
				nodes_std[i][k]=(nodes[i][k+4]-avg)/dev*math.sqrt(weight_features[k])

							   
											
							 
															   
	nodedik=[[0.0 for x in range(num_units)] for y in range(num_units)]
	#if spatial_contiguity==1:
	#	nodedij=[[0 for x in range(num_units)] for y in range(num_units)]
	#else: nodedij=nodedik
	nodedij=nodedik
	#print "nodes_std",nodes_std
	if data_standardlize==1:
		n=0
		for x in nodes_std:
			for y in x:
				if y>2 or y<-2: n+=1
		print "large values (>2, <-2):",n
		print "calculating the matrix ..."
	for i in range(num_units):
		nodedik[i][i]=0
				 
		for j in range(i,num_units):
			if j<=i: continue
			d=0.0
			if data_standardlize==2:
				for k in range(num_features):
					if weight_features[k]<0.00000001: continue
					ai=nodes[i][k+4]
					aj=nodes[j][k+4]
					if ai!=aj: 
						d+=math.sqrt(weight_features[k])
			if data_standardlize<=1:
				d=sum((nodes_std[i][x]-nodes_std[j][x])*(nodes_std[i][x]-nodes_std[j][x]) for x in range(num_features)) #/num_features
			nodedik[i][j]=nodes[i][3]*d
			nodedik[j][i]=nodes[j][3]*d
			#if spatial_contiguity==1:
			#	d=abs(nodes[i][1]-nodes[j][1])+abs(nodes[i][2]-nodes[j][2])
			#	nodedij[i][j]=d
			#	nodedij[j][i]=d
		if i%10==0: print ".",
		#if i>0 and i%(num_units/100) == 0: print str(i)+"/"+str(num_units)
	print
	#print nodedij
	print "creating node neighbors ..."
	#if spatial_contiguity!=1:
	nodedij=nodedik
	node_neighbors = [[]for x in range(num_units)]
	# for i in range(num_units):
		# id1= nodes[i][3]
		# id2= nodes[i][4]
		# neighb=[]
		# for j in range(num_units):
			# if j==i: continue
			# if abs(id1-nodes[j][3])<=1 and abs(id2-nodes[j][4])<=1:
				# neighb.append(j)
		# node_neighbors[i]=neighb[:]
		# if i%10==0: print ".",
	f = open(f2)
	line = f.readline()
	line = f.readline()
	links=0
	ids=[x[0] for x in nodes]
	while line:
		items = line.split(',')
		if len(items)<=2:
			items = line.split('\t')
		if int (items[1]) != int (items[2]):
			id1=int (items[1])
			id2=int (items[2])
			idx1=-1
			idx2=-1
			for i in range(num_units):
				if nodes[i][0]==id1:
					idx1=i
				if nodes[i][0]==id2:
					idx2=i
				if idx1>=0 and idx2>0:
					break
			if idx1>=0 and idx2>0:
			    node_neighbors[idx1].append(idx2)
		line = f.readline()
	f.close()
	print 
	#print node_neighbors
	#create_facility_neighbors()
	#find_NearFacilityList(num_districts)
	#find_near_customer()
	potential_facilities=[x for x in range(num_districts)]
	print "data prepaird!!"
	total_pop=num_units
	return
def create_facility_neighbors():
    return 
    global facility_neighbors
    mindij=[[MAXNUMBER for x in range(num_districts)] for y in range(num_districts)]
    for i in range(num_districts):
        for j in range(num_districts):
            if j<=i: continue
            dlist=[nodedij[x][i]-nodedij[x][j] for x in range(num_units)]
            d=sum(x*x for x in dlist)
            mindij[i][j]=d
            mindij[j][i]=d
    facility_neighbors = [[]for x in range(num_districts)]
    for i in range(num_districts):
        dlist=[[x, mindij[i][x]] for x in range(num_districts) ]
        dlist.sort(key=lambda x:x[1])
        nghrs=[x[0] for x in dlist]
        facility_neighbors[i]=nghrs[:]
def find_near_customer():
    global nearCustomer
    global nearCustomers
    if location_problem>=2 and pmp_I_eaquls_J==1: 
        nearCustomers=NearFacilityList
        nearCustomer=[x for x in range(num_districts)]
        return
    nearCustomer=[-1 for x in range(num_districts)]
    nearCustomers=[[] for x in range(num_districts)]
    dis=[]
    for k in range(num_districts):
        dis=[ [x,nodedij[x][k]] for x in range(num_units)]
        dis.sort(key=lambda x: x[1])
        nearCustomer[k]=dis[0][0]
        nearCustomers[k]=[x[0] for x in dis]


def initialize_instance():
    global num_units
    global num_districts
    global num_facilityCandidate
    global centersID
    global node_groups
    global facilityCost
    global facilityCandidate
    global facilityCapacity
    global nodedik
    global avg_pop
    global total_pop
    global avg_dis_min
    global total_supply
    global fixed_cost_obj
    global max_num_facility
    global adaptive_number_of_facilities
    #solution obj 
    global district_info
    global objective_overload
    global objective
    global biobjective
    global all_solutions
    #best solution 
    global best_solution
    global best_centersID
    global best_biobjective
    global best_objective
    global best_objective_overload
    #global best solution 
    global best_solution_global
    global best_centersID_global
    global best_biobjective_global
    global best_objective_global
    global best_overload_global
    global potential_facilities
    global max_exclusion_list
    num_districts=len(facilityCandidate)
    #num_units=len(nodes)
    total_pop=sum(x[3] for x in nodes)
    #print total_pop,nodes[:10]
	#sum(nodes[x][3] for x in range(num_units))
    node_groups=[-1 for x in range(num_units)]
    if location_problem==0:
        fixed_cost_obj=0
        max_num_facility=num_districts
    if fixed_cost_obj==0:
        facilityCost=[0 for x in range(num_districts)]
    if location_problem==1 and max_num_facility<1:
        max_num_facility=num_districts
        adaptive_number_of_facilities=1
    if location_problem==2:
        if all_units_as_candadate_locations==1:
            facilityCandidate=[x for x in range(num_districts)]
            facilityCost=[0.0 for x in range(num_districts)]
            popa=total_pop*1.0/max_num_facility
            facilityCapacity=[popa for x in range(num_districts)]
        if all_units_as_candadate_locations==0:
            facilityCost=[0.0 for x in range(num_districts)]
            popa=total_pop*1.0/max_num_facility
            facilityCapacity=[popa for x in range(num_districts)]
    if location_problem==3: #pmp
        #num_districts=num_units
        #facilityCandidate=[x for x in range(num_districts)]
        facilityCost=[0.0 for x in range(num_districts)]
        facilityCapacity=[total_pop for x in range(num_districts)]
    centersID=facilityCandidate[:]
    num_facilityCandidate=len(facilityCandidate)
    district_info = [[0,0.0,0.0,0.0,0.0] for x in range(num_districts)]
    total_supply=sum(facilityCapacity)
    #arcpy_print("total demand: "+str(total_pop))
    #arcpy_print("total supply: "+str(total_supply))
    #arcpy_print("avg. distance to nearest facility: "+str(avg_dis_min))

    objective_overload=MAXNUMBER
    obj_balance=MAXNUMBER
    objective=MAXNUMBER
    biobjective=MAXNUMBER
    all_solutions=[]

    #best solution in each start
    best_solution =[] # node_groups[:]
    best_centersID=[]
    best_biobjective=MAXNUMBER
    best_objective=MAXNUMBER
    best_objective_overload = MAXNUMBER

    #global best solution 
    best_solution_global=[]
    best_centersID_global=[]
    best_biobjective_global = MAXNUMBER
    best_objective_global = MAXNUMBER
    best_overload_global = MAXNUMBER
    #if geo_instance==1:
    #    nodedik=[[nodedij[y][facilityCandidate[x]]*nodes[y][3] for x in range(num_districts)] for y in range(num_units)]
    avg_dis_min =sum(nodedik[x][0] for x in range(num_units))/total_pop
    if spatial_contiguity>=1:
        find_near_customer()
    find_NearFacilityList(num_districts)
    if linear_relaxation==1:
        lplocs,sol=location_model_linear_relexation(max_num_facility,0,heuristic_time_limit,0.0001)	
        potential_facilities=[x for x in range(num_districts) if lplocs[x]>0.0001]
        print "Potential facilities by Linear Relax",potential_facilities    
    max_exclusion_list=[0.0 for x in range(num_districts)]

def find_NearFacilityList(nnn):
    global NearFacilityList
    if len(NearFacilityList)>0: return
    NearFacilityList=[]
    n=nnn#num_districts
    if n>num_districts: n=num_districts
    dis=0.0
    print "find_NearFacilityList()",
    for i in range(num_units):
        if i%100==0: print ".",
        fdlist=[ [x,nodedik[i][x]] for x in range(num_districts)]
        fdlist.sort(key=lambda x:x[1])
        flist=[x[0] for x in fdlist[0:n]]
        NearFacilityList.append(flist[:])
    print
    if geo_instance==0:
        return
TB_tabu_list=[]
def pmp_TB():
    global node_groups
    global centersID
    global TB_tabu_list
    #if random.random()>0.5: pmp_TB2()
    best_exch=[] #klist,save
    best_obj=biobjective
    #if centersID in TB_tabu_list: 
    #    return -1
    clist=[x for x in range(num_districts) if centersID[x]>=0]
    random.shuffle(clist)
    klist=[x for x in range(num_districts) if centersID[x]<0]
    random.shuffle(klist)
    for k in clist:
        #u=nearCustomer[k]
        #klist=NearFacilityList[u][num_districts/20:num_districts/5]
        random.shuffle(klist)
        for nk in klist:
            #if nk in clist: continue
            obj=0.0
            new_centersID=centersID[:]
            new_centersID[k]=-1
            new_centersID[nk]=nk
            if new_centersID in TB_tabu_list: continue
            for i in range(num_units):
                j=node_groups[i]
                if j!=k:
                    obj+=min(nodedik[i][j],nodedik[i][nk])
                else:
                    for j in NearFacilityList[i]:
                        if new_centersID[j]>=0:
                            obj+=nodedik[i][j]
                            break
            if obj<best_obj: 
                best_obj=obj
                best_exch.append([k,nk,obj])
                #break
            #print obj,
            if len(best_exch)>=1: break
        if len(best_exch)>=1: break
    if best_exch==[]: 
        #TB_tabu_list.append(centersID[:])
        #print "tabu len=",len(TB_tabu_list)
        return -9
    #print best_exch,best_obj,biobjective
    best_exch.sort(key=lambda x:x[2])
    k=best_exch[0][0]
    centersID[k]=-1
    nk=best_exch[0][1]
    centersID[nk]=nk
    obj=biobjective
    for i in range(num_units):
        j=node_groups[i]
        if j!=k:
            if nodedik[i][j]>nodedik[i][nk]: node_groups[i]=nk
        else:
            for j in NearFacilityList[i]:
                if centersID[j]>=0:
                    node_groups[i]=j
                    break
    if centersID not in TB_tabu_list: TB_tabu_list.append(centersID[:])
    update_district_info()
    #if spatial_contiguity>=1: repair_fragmented_solution()
    #update_centers()
    #print "tb",obj,biobjective,best_exch
    if centersID not in TB_tabu_list: TB_tabu_list.append(centersID[:])
    return 1
def pmp_find(fi,facility2): #find fr for fi
    savings=0.0
    loss=[[x,0.0] for x in range(num_districts)]
    if fixed_cost_obj==1: 
        savings-=facilityCost[fi]
        loss=[[x,-facilityCost[x]] for x in range(num_districts)]
    for i in range(num_units):
        k=node_groups[i]
        if nodedik[i][fi]<nodedik[i][k]: 
            savings+=nodedik[i][k]-nodedik[i][fi]
        else:
            loss[k][1]+=min(nodedik[i][fi],nodedik[i][facility2[i]])-nodedik[i][k]
    loss=[x for x in loss if centersID[x[0]]>=0]
    loss.sort(key=lambda x:x[1])
    fr=loss[0][0]
    profit=savings-loss[0][1]
    return fr,profit

#On the implementation of a swap-based local search procedure for the p-median problem
def pmp_Whitaker():
    global node_groups
    global centersID
    global time_Whitaker
    #global TB_tabu_list
    #if centersID in TB_tabu_list: return -1
    t=time.time()
    facility2=[0 for x in range(num_districts)]
    for i in range(num_units):
        k1=node_groups[i]
        for k2 in NearFacilityList[i]:
            if centersID[k2]>=0 and k2!=k1: 
                facility2[i]=k2
                break
    klist=[x for x in range(num_districts) if centersID[x]<0]
    random.shuffle(klist)
    improved=-1
    best_swap=[-1,-1,0.0]
    for nk in klist:
        k,profit=pmp_find(nk,facility2)
        if profit<=0: continue
        if profit>best_swap[2]:
            best_swap[0]=k
            best_swap[1]=nk
            best_swap[2]=profit
            break

    if best_swap[0]>=0:
        k=best_swap[0]
        nk=best_swap[1]
        centersID[k]=-1
        centersID[nk]=nk
        for i in range(num_units):
            j=node_groups[i]
            if j!=k:
                if nodedik[i][j]>nodedik[i][nk]: node_groups[i]=nk
            else:
                for j in NearFacilityList[i]:
                    if centersID[j]>=0:
                        node_groups[i]=j
                        break
        update_district_info()
        #print [x for x in centersID if x>=0]
        improved=1
    time_Whitaker+=time.time()-t
    #if improved==-1: 
    #    if centersID not in TB_tabu_list: TB_tabu_list.append(centersID[:])
    return improved

def pmp_Whitaker_first():
    global node_groups
    global centersID
    #global TB_tabu_list
    #if centersID in TB_tabu_list: return -1
    facility2=[0 for x in range(num_units)]
    for i in range(num_units):
        k1=node_groups[i]
        for k2 in NearFacilityList[i]:
            if centersID[k2]>=0 and k2!=k1: 
                facility2[i]=k2
                break
    klist=[x for x in range(num_districts) if centersID[x]<0]
    random.shuffle(klist)
    improved=-1
    for nk in klist:
        k,profit=pmp_find(nk,facility2)
        if profit<=0: continue
        #print [x for x in centersID if x>=0]
        centersID[k]=-1
        centersID[nk]=nk
        for i in range(num_units):
            j=node_groups[i]
            if j!=k:
                if nodedik[i][j]>nodedik[i][nk]: node_groups[i]=nk
            else:
                for j in NearFacilityList[i]:
                    if centersID[j]>=0:
                        node_groups[i]=j
                        break
        update_district_info()
        #print [x for x in centersID if x>=0]
        print k,nk,profit,
        improved=1
        break
    #if improved==-1: 
    #    if centersID not in TB_tabu_list: TB_tabu_list.append(centersID[:])
    return improved


def r_r_reselect_location_homoDP(centers,candlist,ulist):
    global node_groups
    global centersID
    global tabu_list_homoDP
    global time_pmp_re_location
    t=time.time()
    best_exch=[] #klist,save
    sol=[-1 for x in ulist]
    best_obj=0.0
    if len(ulist)<=num_units or ulist==[]:
        for i in range(len(ulist)):
            u=ulist[i]
            bestk=-1
            bestd=MAXNUMBER
            for k in centers:
                if nodedik[u][k]<bestd:
                    bestk=k
                    bestd=nodedik[u][k]
            sol[i]=bestk
            best_obj+=bestd
    else:
        sol=node_groups[:]
        best_obj=biobjective
    best_obj-=0.0000001
    clist=centers[:]
    random.shuffle(clist)
    for k in clist:
        klist=candlist
        if len(clist)==len(candlist):
            idx=clist.index(k)
            klist=candlist[idx]
        for nk in klist:
            if nk in tabu_list_homoDP[k]: continue
            if nk in clist: continue
            obj=0.0
            nklist=clist[:]
            nklist.remove(k)
            nklist.append(nk)
            ulist1=[ulist[x] for x in range(len(ulist)) if sol[x]==k]
            #print ulist1,clist
            for i in ulist1:
                obj+=min(nodedik[i][x] for x in nklist)
            ulist2=[x for x in range(len(ulist)) if sol[x]!=k]
            for i in ulist2:
                u=ulist[i]
                c=sol[i]
                obj+=min(nodedik[u][c], nodedik[u][nk])
            #print obj1,obj2,obj,best_obj
            if obj<best_obj: 
                best_obj=obj
                best_exch=[k,nk]
            #else:
            tabu_list_homoDP[k].append(nk)
            if best_exch!=[]: break
    if best_exch==[]: 
        return -1,best_obj,centers
    klist=centers[:]
    k1=best_exch[0]
    k2=best_exch[1]
    klist.remove(k1)
    klist.append(k2)
    #print k1,k2,NearFacilityList[k1].index(k2)
    if len(ulist)<num_units: return 1,best_obj,klist
    node_groups=sol
    ulist=[x for x in range(num_units) if sol[x]==k1]
    for i in ulist:
        bestd=MAXNUMBER
        bestk=-1
        for k in klist:
            if nodedik[i][k]<bestd:
                 bestk=k
                 bestd=nodedik[i][k]
        node_groups[i]=bestk
    ulist=[x for x in range(num_units) if sol[x]!=k1]
    for i in ulist:
        k=sol[i]
        if nodedik[i][k2]<nodedik[i][k]:
            node_groups[i]=k2
    centersID=[-1 for x in range(num_districts)]
    for x in klist: centersID[x]=x
    update_district_info()
    time_pmp_re_location+=time.time()-t
    return 1,best_obj,klist

def pmp_TB2():
    global node_groups
    global centersID
    global TB_tabu_list
    best_exch=[] #klist,save
    best_obj=biobjective-0.00001
    counts=[[x,0] for x in range(num_districts)]  
    for x in all_solutions:
        for i in range(num_districts):
            if x[1][i]>=0: counts[i][1]+=1
    counts.sort(key=lambda x:-x[1])
    while counts[-1][1]==0: counts.pop()
    list_h=[x[0] for x in counts[:max_num_facility*2/3] ]
    clist=[x for x in centersID if x>=0 and x not in list_h]
    klist=[x[0] for x in counts]
    random.shuffle(clist)
    random.shuffle(klist)
    #print klist
    for k in clist:
        for nk in klist:
            if nk in clist: continue
            obj=0.0
            nklist=clist[:]
            nklist.remove(k)
            nklist.append(nk)
            #for i in range(num_units):
            #    obj+=min(nodedik[i][x] for x in nklist)
            ulist=[x for x in range(num_units) if node_groups[x]==k]
            for i in ulist:
                obj+=min(nodedik[i][x] for x in nklist)
            ulist=[x for x in range(num_units) if node_groups[x]!=k]
            for i in ulist:
                c=node_groups[i]
                obj+=min(nodedik[i][c], nodedik[i][nk])
            if obj<best_obj: 
                best_obj=obj
                best_exch=[k,nk]
        #if best_exch!=[]: break
        #if best_obj-biobjective<-biobjective*0.001: break
    if best_exch==[]: 
        #TB_tabu_list.append(centersID[:])
        #print "tabu len=",len(TB_tabu_list)
        return -1
    #print best_exch,best_obj,biobjective
    k=best_exch[0]
    centersID[k]=-1
    k=best_exch[1]
    centersID[k]=k
    for i in range(num_units):
        for k in NearFacilityList[i]:
            if k in centersID:
                node_groups[i]=k
                break
    update_district_info()
    if spatial_contiguity>=1: repair_fragmented_solution()
    k1=best_exch[0]
    k2=best_exch[1]
    #print int(best_obj),
    #print NearFacilityList[k].index(k2),
    #nearidx=num_units
    #for k in clist:
    #    idx=NearFacilityList[k].index(k2)
    #    if idx<nearidx: nearidx=idx
    #print nearidx
    return 1

    prob = pulp.LpProblem("sub pmp",pulp.LpMinimize)
    centers=dlist[:]
    #for x in dlist: centers+=node_neighbors[x]
    for x in dlist: centers+=NearFacilityList[x][:8]
    centers=list(set(centers))
    centers=[x for x in centers if x in ulist]
    xvariables={}
    costs={}
    for i in ulist:
        for j in centers:
            xvariables["x_" +str(i)+ "_"+ str(j)]=pulp.LpVariable("x_" +str(i)+ "_"+ str(j), 0, 1, pulp.LpBinary)
            costs["x_" +str(i)+ "_"+ str(j)]= nodedik[i][j]
    
    obj=""
    for x in xvariables:
        obj += costs[x]*xvariables[x]
    prob += obj

    #con 1
    s=""
    for k in centers:
        s+=xvariables["x_" +str(k)+"_"+str(k)]
    prob +=s == len(dlist)

    #cons 2
    for i in ulist:
        s=""
        for j in centers:
            s+=xvariables["x_" +str(i)+ "_"+ str(j)]
        prob +=s == 1

    for i in ulist:
        for k in centers:
            s= xvariables["x_" +str(i) + "_"+ str(k) ] - xvariables["x_" +str(k)+"_"+str(k)]
            prob +=s <= 0

    #prob.writeLP("_sub_pmp.lp")
    #maxSeconds=heuristic_time_limit/multi_start_count/2
    #loc_sub_mst_file("_locsub.mst",dlist,ulist)
    gap=mipgap
    timelmt=5
    solver=""
    if mip_solver=='cbc':
        solver=pulp.apis.COIN_CMD(mip=1,msg=solver_message,gapRel = gap,options=['vnd on', 'node hybrid', 'rens on'])
    if mip_solver=='cplex':
        solver=pulp.apis.cplex_api.CPLEX_CMD(mip=1,msg=solver_message, options=['set mip tolerances mipgap '+ str(gap), 'set parallel -1'])
    if mip_solver=='gurobi':
        solver=pulp.apis.GUROBI_CMD(mip=1,msg=solver_message,options=[("MIPGap",gap)])
    solver.setTmpDir()
    solver.actualSolve(prob)

    if prob.status<=0:
        ##noprint "model unsolved..."
        return []
    nlist=[]
    for v in prob.variables():
        if (v.varValue >= 0.90):
            items=v.name.split('_')
            i=int(items[1])
            j=int(items[2])
            if i==j:
                nlist.append(i)
    return nlist

def pop_selection(population):
    population1=copy.deepcopy(population)
    population1.sort(key=lambda x:x[0])
    population2=[] #delete identical solution
    #sol=[best_biobjective_global,best_centersID_global[:],best_solution_global[:],best_objective_global,best_objective_fcost_global,best_overload_global,0]
    #population2.append(copy.deepcopy(sol))
    population2.append(copy.deepcopy(population1[0]))
    sdiff=1
    if location_problem==3:
        sdiff=max_num_facility*5/100
        if sdiff<=5: sdiff=5
        
    for x in population1:
        issimilar=0
        for y in population2:
            rlist=[i for i in range(num_districts) if x[1][i] != y[1][i]]
            if len(rlist)>=sdiff: continue
            else:
                if location_problem>=1:
                    issimilar=1
                    break
            ulist=[i for i in range(num_units) if x[2][i] != y[2][i]]
            #diffpop=sum(nodes[i][3] for i in ulist)
            #if len(ulist)<min(num_units*1.0/num_districts,num_units/30.0) and diffpop*100.0/total_pop < min(3.0,100.0/num_districts): #100.0/num_districts: #<5%
            #print len(ulist),
            if len(ulist)<num_units*(solution_similarity_limit/100.0):
                issimilar=1
                break
        if issimilar==0:
            population2.append(copy.deepcopy(x))
        if len(population2)>=max(multi_start_count*2,10):
            break
    return population2


def update_centers():
    global node_groups
    global centersID
    global time_update_centers
    if location_problem==1 or location_problem==0: return
    t=time.time()
    obj=biobjective
    sol=[-1 for x in range(num_units)]
    centers=[]
    for k in range(num_districts):
        if centersID[k]==-1: continue
        kn,ulist=update_center(k)
        for x in ulist: sol[x]=kn
        centers.append(kn)
        centersID[k]=-1
        #print [k,kn,k in ulist],
    node_groups=sol[:]
    for k in centers:
        centersID[k]=facilityCapacity[k]
    #if location_problem==3:
    #    for i in range(num_units):
    #        for k in NearFacilityList[i]:
    #            if centersID[k]>=0:
    #                node_groups[i]=k
    #                break
    obj=biobjective
    update_district_info()
    update_best_solution()
    #print "updatecenters",biobjective-obj
    time_update_centers+=time.time()-t
    #print obj,obj-biobjective

def update_center(k):
    if location_problem==3: #pmp
        return update_center_pmp(k)
    ulist=[x for x in range(num_units) if node_groups[x]==k]
    if ulist==[]: return k,[]
    best_cost=sum(nodedik[x][k] for x in ulist)
    best_center=k
    if pmp_I_eaquls_J==1:
        for i in ulist:
            cost=sum(nodedik[x][i] for x in ulist)
            if cost<best_cost:
                best_cost=cost
                best_center=i
    if pmp_I_eaquls_J==0:
        for i in range(num_districts):
            cost=sum(nodedik[x][i] for x in ulist)
            if cost<best_cost:
                best_cost=cost
                best_center=i
    return best_center,ulist

    
def print_solution():
    arcpy_print("_______________final solution_______________")
    for i in range(num_units):
        s=""
        for x in nodes[i]:
            s+=str(x)+"\t"
        k=node_groups[i]
        kunit=centersID[k]
        s+=str(nodes[kunit][4])+"\t"
        selected=-1
        if i in facilityCandidate:
            selected=0
            if i in centersID:
                selected=1
        s+=str(selected)
        arcpy_print(s)

def random_centers():
    global centersID
    centers=[]
    while 1:
        if len(centers)==max_num_facility: break
        k=random.randint(0,num_districts-1)
        if k not in centers: centers.append(k)
    for x in range(num_districts):centersID[x]=-1
    for x in centers: centersID[x]=x
    for i in range(num_units):
        for k in NearFacilityList[i]:
            if centersID[k]>=0:
                node_groups[i]=k
                break
    update_district_info()
    return 1
        
def k_medoids_sampling(allulist,cids):
    centers=cids
    num=len(allulist)
    if centers==[]:
        while 1:
            idx=random.randint(0,num-1)
            nid=allulist[idx]
            if nid not in centers:
                centers.append(nid)
            if len(centers)==max_num_facility:
                break
    sol=[-1 for x in allulist]
    loop=0
    distance_obj=MAXNUMBER
    #k-means
    while 1:
        loop+=1
        for i in range(num):
            for k in NearFacilityList[allulist[i]]:
                if k in centers:
                    sol[i]=k
                    break
        obj=0.0
        for k in range(max_num_facility):
            ulist=[allulist[x] for x in range(num) if sol[x]==centers[k]]
            if len(ulist)==0: continue
            clist=ulist
            cid=centers[k]
            mindis=sum(nodedik[x][cid] for x in ulist)
            for i in clist:
                dsum=sum(nodedik[x][i] for x in ulist)
                if dsum<mindis:
                    mindis=dsum
                    cid=i
            centers[k]=cid
            obj+=mindis
        #print loop,distance_obj,sol
        if obj<distance_obj:
            distance_obj=obj
        else:
            break
    return distance_obj,centers


def ils_pmp(numf,multistarts,timelimit,myseed):
    global multi_start_count
    global seed
    global best_objective
    global best_biobjective
    global best_objective_overload
    global best_biobjective_global
    global objective
    global biobjective
    global objective_overload
    global node_groups
    global centersID
    global district_info
    #global node_neighbors
    global move_count
    global region_pool
    global pool_index
    global is_spp_modeling
    global spp_loops
    global adj_pop_loops
    global pop_dis_coeff
    global all_solutions
    global assignment_operators_selected
    global max_num_facility
    global heuristic_time_limit
    global adjust_mutation_rate
    global large_facility_cost
    global initial_solution_method
    global avg_dis_min
    global spatial_contiguity
    global location_tabu_list
    global max_exclusion_list
    max_num_facility=numf
    if location_problem==0:
        max_num_facility=num_districts
    initialize_instance()
    heuristic_time_limit=timelimit
    all_solutions=[]
    multi_start_count=multistarts
    seed=myseed
    if seed<0:
        seed=random.randint(0,100)
    random.seed(seed)
    arcpy_print("seed: "+str(seed))
    region_pool=[]
    t=time.time()
    ga_time=0.0
    best_biobjective_global = MAXNUMBER
    best_biobjective = MAXNUMBER
    district_info = [[0,0.0,0,0,0] for x in range(num_districts)]
    population=[] #all
    pool_index=[]
    if is_spp_modeling>=1:
        pool_index=[[] for x in range(num_districts)]
    t_ga=time.time()
    multi_start=0
    population=generate_initial_solution()
    
    if len(max_exclusion_list)<1:
        max_exclusion_list=[0.0 for x in range(num_districts)]
    not_improved_global=0
    improved_time=time.time()
    loop=-1
    while 1:
        loop+=1
        r=random.random()
        sidx = int(min(multi_start_count,len(population))* pow(r,2)*0.999)
        r=random.random()
        node_groups=population[sidx][2][:]
        centersID=population[sidx][1][:]
        update_district_info()
        current=" current: " +str(int(population[sidx][0]))+" -> "
        not_improved_global+=1
        old_ids=centersID[:]
        best_obj=best_biobjective_global
        cur_obj=biobjective
        op=0
        sta=0
        #check_pmp_sol(0)
        if random.random()>0.5:
            idx,sta=r_r_new_location(-1)
        else:
            assign_ruin_recreate(-1)
            op=1
        if spatial_contiguity>=1: VND_local_search()
        update_centers()
        #if spatial_contiguity>=1: repair_fragmented_solution()
        #    check_pmp_sol(4)
        #check_pmp_sol(1)
        update_district_info()
        update_best_solution()
        update_region_pool_all()
        population.append([biobjective,centersID[:],node_groups[:],objective,objective_fcost,objective_overload,max(0,total_pop-objective_supply)])
        s=""
        if biobjective<best_obj:  #0.1%
            s="*"
            impp=((biobjective-best_obj)/biobjective)*1000 #-0.1%
            not_improved_global+=int(max_loops_solution_not_improved*impp)
            if not_improved_global<0: not_improved_global=0
        elif biobjective<cur_obj-0.000001: 
            s="#"
        else: s="-"
        s=str(op)+" "+s
        bnf=sum(1 for x in best_centersID_global if x>=0)
        s+="Search loop "+str(loop) + " best: " +str(bnf)+" "+str(int(best_biobjective_global))+" "+str(int(best_objective_global))+" "+str(int(best_overload_global))
        nf=sum(1 for x in centersID if x>=0)
        s+=current + str(nf)+" "+str(int(biobjective))+" "+str(int(objective))+" "+str(int(objective_overload))
        f=[x for x in range(num_districts) if max_exclusion_list[x] <=best_biobjective_global+0.00001]
        n=sum(1 for x in range(num_districts) if best_centersID_global[x]>=0 and x in f) 
        s+=" "+ str(n)
        surplus=int(objective_supply-total_pop)
        psize=int(total_pop/nf)
        s+= " Info: " +str(len(population))+" " +str(len(potential_facilities)) +" " + str(surplus)+"/"+str(psize)+" " +str(not_improved_global)+ " "+str(int(time.time()-t_ga))
        arcpy_print(s)
        if not_improved_global >= max_loops_solution_not_improved: break
        population.sort(key=lambda x:x[0])
        population=pop_selection(population)
        #while len(population)>multi_start_count*5:
        #    population.pop()
        all_solutions=population
        if time.time()-t_ga > heuristic_time_limit:  break
    #post procedure
    population.sort(key=lambda x:x[0])
    node_groups=best_solution_global[:] #population[0][2][:]
    centersID=best_centersID_global[:]
    update_district_info()
    ga_time=time.time()-t_ga
    print "Heuristic solution:",biobjective,objective_fcost,objective,objective_overload,ga_time
    t_spp=time.time()
    if is_spp_modeling>=1:
        arcpy_print("SPP modelling..."+str(len(region_pool)) )
        sppmodel(heuristic_time_limit,0.000001)
        #clist=[x for x in centersID if x>=0]
        #for i in range(num_units):
        #    for k in NearFacilityList[i]:
        #        if k in clist: 
        #            node_groups[i]=k
        #            break
        update_district_info()
        update_centers()
        update_best_solution()
        population.append([biobjective,centersID[:],node_groups[:],objective,objective_fcost,objective_overload,0])
        print "spp solution:",biobjective,objective_fcost,objective,objective_overload,time.time()-t_spp, time.time()-t_spp
    print "final solution:",biobjective,objective_fcost,objective,objective_overload
    population.sort(key=lambda x:x[0])
    all_solutions=population
    sta=check_solution_continuality_feasibility(node_groups)
    print "Areal continuality(1 yes, 0 no):",sta
    if spatial_contiguity==1:
        if sta==0:    return "infeasible final solution: continuality"
    ##noprint "time",time.time()-t
    search_stat()
    print "location_tabu_list",len(location_tabu_list)
    return [best_biobjective_global,best_objective_global,best_overload_global,centersID,best_solution_global]

def pmp_sampling_solution(psize,sample_size):
    global node_groups
    global centersID
    global all_solutions
    t=time.time()
    print "psize,sample_size",psize,sample_size
    for idx in range(psize):
        best_biobjective = MAXNUMBER
        random.seed(seed+idx*100)
        ulist=[]
        if sample_size==num_units: ulist=range(num_units)
        else:
            while len(ulist)<sample_size: 
                u=random.randint(0,num_units-1)
                if u not in ulist: ulist.append(u)
        # find k_medoids
        random.shuffle(ulist)
        obj,centers=k_medoids_sampling(ulist,[])
        loop=0
        # #improve by TB method
        # while loop<max_num_facility:
            # loop+=1
            # random.shuffle(centers)
            # random.shuffle(ulist)
            # sta,obj,centers=r_r_reselect_location_homoDP(centers,ulist,ulist)
            # #print "multi_start",multi_start,loop,obj,time.time()-ga_time
            # if sta<=0: break
        #initial solution
        #print centers
        centersID=[-1 for x in range(num_units)]
        for x in centers: centersID[x]=x
        obj=0.0
        for i in range(num_units):
            bestd=MAXNUMBER
            bestk=-1
            for k in centers:
                if nodedik[i][k]<bestd:
                    bestd=nodedik[i][k]
                    bestk=k
            node_groups[i]=bestk
            obj+=bestd
        update_district_info()
        if spatial_contiguity>=1:
            repair_fragmented_solution_large_instance()
        #update_centers()
        all_solutions.append([objective,centersID[:],node_groups[:]])
        update_best_solution()
        print "init. solution",idx,biobjective,objective,objective_overload,time.time()-t
        #print [x for x in centersID if x>=0]
    #return all_solutions
tabu_list_homoDP=[]
def ils_large_homoDP(numf,multistarts):
    global tabu_list_homoDP
    global multi_start_count
    global seed
    global best_objective
    global best_biobjective
    global best_objective_overload
    global best_biobjective_global
    global objective
    global biobjective
    global node_groups
    global centersID
    global district_info
    global region_pool
    global pool_index
    global is_spp_modeling
    global spp_loops
    global all_solutions
    global max_num_facility
    global spatial_contiguity
    max_num_facility=numf
    #initialize_instance()
    find_NearFacilityList(num_units)
    find_near_customer()
    #for x in NearFacilityList: print x[:10]
    node_groups=[-1 for x in range(num_units)]
    centersID=[-1 for x in range(num_districts)]
    all_solutions=[]
    tabu_list_homoDP=[[] for x in range(num_districts)]
    multi_start_count=multistarts

    region_pool=[]
    ga_time=time.time()
    best_biobjective_global = MAXNUMBER
    district_info = [[0,0.0,0,0,0] for x in range(num_districts)]
    population=[] #all
    pool_index=[]
    if is_spp_modeling>=1:
        pool_index=[[] for x in range(num_districts)]
    t_ga=time.time()
    candidate_centers=[]
    centers=[]
    print "step 1: sampling and sub-problem solving ..."
    #sampling and solving sub-problems
    sample_size=max(10*max_num_facility, num_units/10)
    if sample_size>num_units: sample_size=num_units
    num_samplling=multistarts
    if num_samplling<=0:
        num_samplling=3*num_units/sample_size
    #if num_samplling<10*max_num_facility: num_samplling=10*max_num_facility
    #if num_samplling<multi_start_count: num_samplling=multi_start_count
    pmp_sampling_solution(num_samplling,sample_size)
    for x in all_solutions:
        ids=[y for y in x[1] if y>=0]
        candidate_centers+=ids
        print len(candidate_centers),ids
    candidate_centers=list(set(candidate_centers))

    """
    print "step 2: solution for original problen"
    obj,centers=k_medoids_sampling(candidate_centers,centers)
    clusters=[-1 for x in candidate_centers]
    print obj,centers
    num=len(candidate_centers)
    random.shuffle(centers)
    random.shuffle(candidate_centers)
    for i in range(num):
        bestk=-1
        bestd=MAXNUMBER
        u=candidate_centers[i]
        for k in centers:
            if nodedik[u][k] < bestd:
                bestd=nodedik[u][k]
                bestk=k
        clusters[i]=bestk
    """
    print "step 3: solution improvement..."
    all_solutions.sort(key=lambda x:x[0])
    tabu_list_homoDP=[[] for x in range(num_districts)]
    not_improved=0
    loop=0

    if multi_start_count>10: multi_start_count=10
    while 1:
        loop+=1
        r=random.random()
        sidx = int(min(multi_start_count,len(all_solutions))* pow(r,2)*0.999)
        node_groups=all_solutions[sidx][2][:]
        centersID=all_solutions[sidx][1][:]
        update_district_info()
        not_improved+=1
        best_obj=best_biobjective_global-0.000001
        obj_begin=biobjective
        #candidate_centers_selector=candidate_centers
        #if not_improved<max_loops_solution_not_improved*3: 
        #    candidate_centers_selector=candidate_group
        r=random.random()
        sta=1
        if r<0.5:
            print "rr",
            assign_ruin_recreate(-1)
        else:
            print "tb",
            pmp_Whitaker()
        if spatial_contiguity>=1:
            repair_fragmented_solution_large_instance()
        update_region_pool_all()
        update_best_solution()
        if spatial_contiguity>=1:
            pmp_VND_search()#VND_local_search()
            repair_fragmented_solution_large_instance()
        update_region_pool_all()
        update_best_solution()

        if objective<best_obj: 
            print "*",
            delta=max_loops_solution_not_improved*(best_obj-objective)/objective*1000 #100==1%，1000=0.0%
            not_improved=max(0,not_improved-int(delta))
        else:    print "#",
        sol_exist=0
        for x in all_solutions:
            if x[1]==centersID:
                sol_exist=1
                break
        #if sta==-10:
        #    del all_solutions[sidx]
        if sol_exist==0 and sta==1:
            all_solutions.append([biobjective,centersID[:],node_groups[:] ])
            all_solutions.sort(key=lambda x:x[0])
        while len(all_solutions)>2*multi_start_count:all_solutions.pop()
        #if not_improved>=max_loops_solution_not_improved/10:
        #    while 1:
        #        if all_solutions[-1][0]>1.05*all_solutions[0][0]: 
        #            all_solutions.pop()
        #       else: break
        #print best_biobjective_global,best_objective, objective
        print "Search loop",loop,"best",best_biobjective_global,"current",obj_begin,"->",objective, len(all_solutions),not_improved, time.time()-ga_time
        if not_improved>=max_loops_solution_not_improved: break
        #if time.time()-t_ga > heuristic_time_limit:  break
    node_groups=best_solution_global
    centersID=best_centersID_global
    update_district_info()
    if is_spp_modeling>=1:
        sppmodel(heuristic_time_limit,0.000001)
        update_district_info()
        #update_centers()
        update_best_solution()
        print "spp",objective, time.time()-ga_time
    all_solutions.append([biobjective,centersID[:],node_groups[:] ])
    all_solutions.sort(key=lambda x:x[0])
    node_groups=best_solution_global
    centersID=best_centersID_global
    update_district_info()
    ga_time=time.time()-t_ga
    search_stat()
    homoDP_stat()
    print "solution connectivity ..."
    for k in range(num_districts):
        if centersID[k]<0:continue
        ulist=[x for x in range(num_units) if node_groups[x]==k]
        #fu=frag_unit_minority(ulist)
        print "cluster",k,len(ulist),"connectivity:",check_continuality_feasibility(node_groups,k),
        if spatial_contiguity==2:
            ulist=[x for x in range(num_units) if node_groups[x]==k]
            print district_with_multiparts(ulist),
        print 
    return [best_biobjective_global,best_objective_global,best_overload_global,centersID,best_solution_global]
def homoDP_stat():
    attrs=len(nodes_std[0])
    print "-----------------Stats for each attribute--------------------"
    print "Num of attrs", attrs
    print "weights",weight_features,len(weight_features)
    print "Idx----Avg--------Dev---------SSE-------SST------R2-------F-test---------"
    effective_attr=0
    for idx in range(attrs):
        if weight_features[idx]<0.0000001: continue
        effective_attr+=1
        mean=sum(x[idx+4] for x in nodes) / num_units
        sst=sum(weight_features[idx]*pow(x[idx+4]-mean,2) for x in nodes)
        if sst<0.0000001:
            print "Fld"+str(idx), sst
            continue
        dev=pow(sst/num_units,0.5)
        sse=0.0
        means=[0.0 for x in range(num_districts)]
        for k in range(num_districts):
            if centersID[k]<0:continue
            ulist=[x for x in range(num_units) if node_groups[x]==k]
            sumd=sum(nodes[x][idx+4] for x in ulist)
            means[k]=sumd/len(ulist)
        for i in range(num_units):
            k=node_groups[i]
            d=nodes[i][idx+4]-means[k]
            sse+=d*d
        r2=(sst-sse)/sst
        if r2>0.99999: 
            f=10000
        else:
            f=r2/(max_num_facility-1) * (num_units-max_num_facility) / (1-r2)
        print "Fld"+str(idx), mean, dev, sse,sst, r2, f
    SST_all=0.0
    SSE_all=0.0
    for idx in range(attrs):
        if weight_features[idx]<0.0000001: continue
        mean=sum(x[idx] for x in nodes_std) / num_units
        sst=sum(pow(x[idx]-mean,2) for x in nodes_std)
        SST_all+=sst
        #if mean<0.00000000001:  continue
        dev=pow(sst/num_units,0.5)

        sse=0.0
        means=[0.0 for x in range(num_districts)]
        for k in range(num_districts):
            if centersID[k]<0:continue
            ulist=[x for x in range(num_units) if node_groups[x]==k]
            sumd=sum(nodes_std[x][idx] for x in ulist)
            means[k]=sumd/len(ulist)
        for i in range(num_units):
            k=node_groups[i]
            d=nodes_std[i][idx]-means[k]
            sse+=d*d
        SSE_all+=sse
    print SST_all, SSE_all
    r2=(SST_all-SSE_all)/SST_all
    f=r2/(max_num_facility-1) * (num_units-max_num_facility) / (1-r2)
    print "Stat:",  "SSE=",SSE_all,"SST=",SST_all, "R2=",r2,"F-stat=", f

    print "-----------------------mean values of regions-------------------------"
    s="Region"
    klist=[x for x in centersID if x>=0]
    for k in range(max_num_facility):
        s+="\t"+str(klist[k])
    print s

    for idx in range(attrs):
        s=""
        if weight_features[idx]<0.0000001: continue
        s+="Fld"+str(idx)
        for k in range(max_num_facility):
            ulist=[x for x in range(num_units) if node_groups[x]==klist[k]]
            mean=sum(nodes[x][idx+4] for x in ulist) / len(ulist)
            s+="\t"+str(mean)
        print s
    print "-----------------------end of Stat-------------------------"
    #error=0.0
    #for k in range(num_districts):
    #    if centersID[k]<0:continue
    #    ulist=[x for x in range(num_units) if node_groups[x]==k]
    #    for idx in range(attrs):
    #        sumd=sum(nodes[x][idx+4] for x in ulist)
    #        mean=sumd/len(ulist)
    #        d= mean-nodes[k][idx+4]
    #        error+=d*d
    #print "error",error

def search_stat():
    arcpy_print("----------------search statistics----------------------")
    arcpy_print("one unit move, move and time: "+ str(count_op[0])+ ", " +str(time_op[0]) )
    arcpy_print("two unit move, move and time: "+ str(count_op[1])+ ", " +str(time_op[1]) )
    arcpy_print("three unit move, move and time: "+ str(count_op[2])+ ", " +str(time_op[2]) )
    arcpy_print("location swap time: "+ str(time_location[0]) )
    arcpy_print("location drop time: "+ str(time_location[1]) )
    arcpy_print("location add time: "+ str(time_location[2]) )
    arcpy_print("location add-drop time: "+ str(time_location[3]) )
    arcpy_print("location multi-exchange time: "+ str(time_location[4]) )
    arcpy_print("r_r_reselect_location_pmp time: "+ str(time_location[5]) )
    arcpy_print("location TB heur. time: "+ str(time_location[7]) )
    arcpy_print("TB_Whitaker time: "+str(time_Whitaker))
    arcpy_print("location PMP sub_mip time: "+ str(time_location[8]) )
    arcpy_print("location CFLP sub_mip time: "+ str(time_location[9]) )
    arcpy_print("location PMP TB time: "+ str(time_pmp_re_location) )
    arcpy_print("repair time: "+ str(time_repair) )
    arcpy_print("check edge unit time: "+str(time_check_edge_unit))
    arcpy_print("update_centers time: "+ str(time_update_centers) )
    arcpy_print("spp regions: "+ str(len(region_pool)) )
    arcpy_print("spp pooling time: "+ str(time_spp) )
    arcpy_print("connectivity check time: "+ str(time_check))
    arcpy_print("time for ruin and recreate: " + str(time_ruin_recreate))
    if spatial_contiguity==1:
        sta=check_solution_continuality_feasibility(best_solution_global)
        arcpy_print("solution on continuality (0 no, 1 yes) : "+str(sta))
    arcpy_print("----------------end of search statistics----------------")
