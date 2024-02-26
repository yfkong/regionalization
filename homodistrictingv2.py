# -*- coding: utf-8 -*-
## ILS algorithm for regionalization Problems, version 2
## yfkong@henu.edu.cn, Apr.,2024
 
import sys,os,random,time,copy,math
#ArcGIS
has_arcpy=0
try:
	import arcpy
	has_arcpy=1
except:
	has_arcpy=0

#constant
MAXNUMBER=1.0e+20
MINNUMBER=1.0e-10
#instance info
nodes=[]
nodes_std=[] #for homoDP only
weight_features=[] #for homoDP only
num_units=-1
nodedij=[]
nodedik=[]	#weighted cost from i to k, =nodedij*nodes[][3] 
nodendij=[] #network distance
node_neighbors=[]
facility_neighbors=[]
total_pop=0
avg_pop=0
total_supply=0
distance_type=0 #0 Euclidean, 1 Manhattan, 2 Geo 3 Network
all_units_as_candadate_locations=0
facilityCandidate=[] #attribute = 0,1,2...
facilityCapacity=[]
facilityCost=[]
num_facilityCandidate=-1
num_districts=-1 # number of service areas/facilities
avg_dis_min=1.0
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
penalty_on_demand=10000.0
pop_deviation=0.00 #for pdp, 5%
cflp_max_service_radius=10000.0
cflp_service_radius_preference=10000.0
cflp_min_service_coverage_percentage=100.0
envy_service_distance=10000.0
envy_objective_weight_of_travelcost=1 #pmp only
envy_objective_weight=0
envy_service_objective=0 
#0: no envy; 
#1: envy obj with weight; 
#2: envy constraint with CV; 
#3: envy with 0 weight on sum_dis, MELP
#4: envy(abs_dev) obj with weight; 
envy_coefficient_variation=0.5

#current solution
centersID=[]
node_groups=[]
district_info=[] #[[0,0,0.0] for x in range(num_districts)] # solution
objective_overload=0
obj_balance=MAXNUMBER
objective=MAXNUMBER
objective_fcost=0.0
biobjective=MAXNUMBER
objective_supply=0.0
objective_envy=0.0
objective_covered=0
objective_rmax_not_covered=0
objective_pentalty_on_overload=0.0
objective_pentalty_on_rmax=0.0
objective_pentalty_on_covering=0.0

given_solution=0 #reserved
all_solutions=[]

#best solution in each start
best_solution =[] # node_groups[:]
best_centersID=[]
best_biobjective=MAXNUMBER
best_objective=MAXNUMBER
best_objective_overload = MAXNUMBER
best_objective_fcost = 0.0
#global best solution 
#best_centers_global=[]
best_solution_global=[]
best_centersID_global=[]
best_biobjective_global = MAXNUMBER
best_objective_global = MAXNUMBER
best_objective_fcost_global = 0.0
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
	#	 return
	ulist=[x for x in  range(num_units) if node_groups[x]==rid]
	if ulist==[]:
		#print "empty area:",rid,node_groups
		return
	cost1=district_info[rid][2]
	cost2=sum(nodedik[x][rid] for x in ulist)
	if abs(cost1-cost2)>0.001: print rid,cost1,cost2
	obj=district_info[rid][2]+district_info[rid][4]*pop_dis_coeff
	idx=int(obj*100000)
	idx+=sum(x*(x+3) for x in ulist)
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
			feasible=0	#infeas.
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
		flist=frag_unit_minority(ulist,rid)
		if flist==[]: 
			time_check+=time.time()-t
			return 2
	time_check+=time.time()-t
	return 0	#infeasible

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
		k=node_groups[ulist[0]]  
		flist=frag_unit_minority(ulist,k)
		time_check+=time.time()-t
		if flist==[]: return 2
	return 0	#infeasible

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
	#	 p=(4.0*math.sqrt(nu)-4)/nu
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
		district_info[k][0] = len(ulist)
		district_info[k][1] = sum(nodes[x][3] for x in ulist)
		district_info[k][2] = sum(nodedik[x][k] for x in ulist)
		if location_problem==9: district_info[k][2] = sum(nodedij[x][k] for x in ulist)
		district_info[k][3] = facilityCapacity[k] 
		district_info[k][4] = max(0.0,district_info[k][1]-facilityCapacity[k]) # -district_info[k][3]
		if location_problem==3: district_info[k][4]=0 #pmp
		if location_problem==2: #pdp,edp
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
	#	 avg_dis_min=objective/total_pop
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
			nflist=frag_unit_minority(ulist,k)
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
	#print "len(frag_units)",len(frag_units)
	return frag_units	 

def frag_unit_minority(ulist,k):
	final_list=[]
	ulist1=ulist[:]
	ulist1.remove(k)
	ulist2=[k]
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
		if k in x[1]: continue        
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
	#update_centers()
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
	#repair_fragmented_solution()
	time_repair+=time.time()-t
	#print objective,

#repair the fragmented solution
def repair_fragmented_solution():
	if spatial_contiguity==0: return
	#if num_units/max_num_facility>100:
	#	 repair_fragmented_solution_large_instance()
	#	 return
	global node_groups
	global centersID
	global time_repair
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
	#print int(objective),

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

def r_r_pathrelink_AP():
	global node_groups
	for i in range(num_units):
		klist=[x[2][i] for x in all_solutions]
		r=random.random()
		idx=int(r*r*0.999*len(klist))
		#idx=random.randint(0,len(klist)-1)
		k=klist[idx]
		node_groups[i]=k	
	update_district_info()
	if spatial_contiguity>=1:
		repair_fragmented_solution()

def reassign_geo(dlist,nlist):
	global node_groups
	if spatial_contiguity>=1:
		ulist=[x for x in range(num_units) if node_groups[x] in dlist]
		for x in ulist: node_groups[x]=-1
		for x in ulist:
			for k in NearFacilityList[x]:
				if centersID[k]<0: continue
				node_groups[x]=k
				break
		update_district_info()
		return

	ulist=[[x,nodes[x][3]] for x in range(num_units) if node_groups[x] in dlist]
	ulist.sort(key=lambda x:-x[1])
	ulist=[x[0] for x in ulist]
	for x in ulist: node_groups[x]=-1
	for x in nlist: 
		district_info[x][1]=0
		district_info[x][3]=facilityCapacity[x]

	for x in ulist:
		for k in NearFacilityList[x]:
			#if k not in nlist: continue
			if centersID[k]<0: continue
			if district_info[k][1]+nodes[x][3]<=district_info[k][3]:
				node_groups[x]=k
				district_info[k][1]+=nodes[x][3]
				break
	ulist=[x for x in ulist if node_groups[x]<0]
	for x in ulist:
		for k in NearFacilityList[x]:
			if centersID[k]<0: continue	 
			#if k not in nlist: continue
			node_groups[x]=k
			break
	update_district_info()

def reassign(dlist,nlist):
	global node_groups
	if spatial_contiguity>=1:
		return reassign_geo(dlist,nlist)
	ulist=[[x,nodes[x][3]] for x in range(num_units) if node_groups[x] in dlist and node_groups[x] not in nlist]
	ulist.sort(key=lambda x:-x[1])
	ulist=[x[0] for x in ulist]
	for x in ulist: node_groups[x]=-1
	for x in nlist: 
		if x in dlist and x in nlist: continue 
		if x in dlist and x not in nlist: 
			district_info[x][1]=0
			district_info[x][3]=0
		if x not in dlist and x in nlist: 
			district_info[x][1]=0
			district_info[x][3]=facilityCapacity[x]

	for x in ulist:
		for k in NearFacilityList[x]:
			if k not in nlist: continue
			if centersID[k]<0: continue
			if district_info[k][1]+nodes[x][3]<=district_info[k][3]:
				node_groups[x]=k
				district_info[k][1]+=nodes[x][3]
				break
	ulist=[x for x in ulist if node_groups[x]<0]
	for x in ulist:
		for k in NearFacilityList[x]:
			if centersID[k]<0: continue	 
			node_groups[x]=k
			break
	update_district_info()

#r&r method
#assign the removed units to solution
def repair_partial_solution():
	global node_groups
	#if spatial_contiguity!=1:
	#	 repair_partial_solution_nonconnectivity()
	#	 return
	ulist=[x for x in range(num_units) if node_groups[x]==-1] # units not assigned
	for x in ulist:
		d=MAXNUMBER
		k=-1
		for y in range(num_districts):
			if centersID[y]<0: continue
			if nodedik[x][y]<d:
				d=nodedik[x][y]
				k=y
		node_groups[x]=k
	update_district_info()
  
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
	#	 ##noprint "check_solution_continuality_feasibility!!!"
	#	 return improve
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

def region_neighbors(rid):
	if geo_instance==0:
		return region_neighbors2(rid) #for benchmark instance
	ulist=[x for x in range(num_units) if node_groups[x]==rid]
	rlist=[]
	for x in ulist:
		for y in node_neighbors[x]:
			k=node_groups[y]
			if k==rid: continue
			if k not in rlist: rlist.append(k)
	if rid in rlist: rlist.remove(rid)
	random.shuffle(rlist)
	return rlist

def region_neighbors2(rid):
	ulist=[x for x in range(num_units) if node_groups[x]==rid]
	clist=[[x,MAXNUMBER] for x in range(num_districts) if centersID[x]>=0 and x!=rid]
	for i in range(len(clist)):
		k=clist[i][0]
		clist[i][1]=min(nodedij[x][k] for x in ulist)
	clist.sort(key=lambda x:x[1])
	#print clist
	rlist=[x[0] for x in clist if x[1]<=3*search_radius]
	rlist=rlist[:6]
	random.shuffle(rlist)
	k=min(6,len(rlist)-1)
	#print len(rlist),search_radius
	return rlist

#return the neighor regions of unit nid
def neighbor_regions(nid):
	rid=node_groups[nid]
	if spatial_contiguity>=1:
		rlist2=[node_groups[x] for x in node_neighbors[nid] if node_groups[x]!=rid]
		rlist=list(set(rlist2))
		if len(rlist2)>1: random.shuffle(rlist2)
		return rlist
	if spatial_contiguity>=0 and random.random()>0.5: #testing ??? 
		knn=4
		#if knn< num_units/100: knn=num_units/100
		rlist=[]
		for k in NearFacilityList[nid]:
			if k==rid: continue
			if centersID[k]<0: continue
			rlist.append(k)
			if len(rlist)>=knn: break 
		return rlist
	rlist2=[node_groups[x] for x in node_neighbors[nid] if node_groups[x]!=rid]
	rlist=list(set(rlist2))
	if len(rlist2)>1: random.shuffle(rlist2)
	return rlist

#evaluate the possible moves
#return 0 or 1 according to the acceptance rule
#acceptance rule: sa, rrt, oba or others
def isbetter(obj_old,obj_new,bal_old,bal_new):
	#fixed cost is not considered
	penalty=avg_dis_min*pop_dis_coeff
	biobj_new=objective-obj_old+obj_new+bal_new*penalty
	biobj_old=objective+objective_overload*penalty
	if biobj_new<biobj_old-0.000001: #better
		return 1
	if bal_new>bal_old: return 0 #more overload, or not reduce the overload
	if acceptanceRule=="rrt":
		if objective+obj_new-obj_old-best_objective<best_objective*1.0/num_units:
			if random.random()>0.5: return 0
			else: 
				#print ".",
				return 1
	if acceptanceRule=="sa":
		#minobj=min(last_loop_obj,biobjective)
		#a=-(biobj_new-minobj)/minobj*num_units/SA_temperature
		a=-(objective+obj_new-obj_old-best_objective)/best_objective*num_units/SA_temperature
		if a>=0: return 1
		p=math.exp(a)
		if random.random()< p:	return 1
	return 0

def one_unit_move(): #FA
	global node_groups
	global time_op
	global count_op
	if location_problem==3:
		return one_unit_move_pmp_geo()
	if spatial_contiguity>=1:
		return one_unit_move_geo()
	t=time.time()
	ulist=find_edge_units()
	difflist=[]
	for x in ulist:
		rlist=[r for r in NearFacilityList[x] if centersID[r]>=0]
		k=node_groups[x]
		difflist.append([x,nodedik[x][rlist[0]]-nodedik[x][k]])
	difflist.sort(key=lambda x:x[1])
	ulist=[x[0] for x in difflist]
	find_best_solution=0
	for nid in ulist:
		rid=node_groups[nid]
		demand=nodes[nid][3]
		klist=neighbor_regions(nid)
		if len(klist)>1:
			klist2=[x for x in NearFacilityList[nid] if x in klist]
			klist=klist2
		for k in klist:
			#if k==rid: print "k==rid", nid,rid,neighbor_regions(nid)
			if demand+district_info[k][1]>facilityCapacity[k]: continue
			savings=nodedik[nid][rid]-nodedik[nid][k]
			overload_old=district_info[rid][4]+district_info[k][4]
			overload_new= max(0,district_info[rid][1]-nodes[nid][3]-facilityCapacity[rid]) +max(0,district_info[k][1]+nodes[nid][3]-facilityCapacity[k])
			savings*=envy_objective_weight_of_travelcost
			savings+= (overload_old-overload_new)*avg_dis_min*pop_dis_coeff
			if envy_service_objective==1:
				if nodedij[nid][rid]>envy_service_distance:
					savings+= nodes[nid][3]*envy_objective_weight*(nodedij[nid][rid]-envy_service_distance)*(nodedij[nid][rid]-envy_service_distance)
				if nodedij[nid][k]>envy_service_distance:
					savings-= nodes[nid][3]*envy_objective_weight*(nodedij[nid][k]-envy_service_distance)*(nodedij[nid][k]-envy_service_distance)
			if envy_service_objective==4:
				if nodedij[nid][rid]>envy_service_distance:
					savings+= nodes[nid][3]*envy_objective_weight*(nodedij[nid][rid]-envy_service_distance)
				if nodedij[nid][k]>envy_service_distance:
					savings-= nodes[nid][3]*envy_objective_weight*(nodedij[nid][k]-envy_service_distance)
			if savings<=0: continue
			node_groups[nid]=k
			count_op[0]+=1
			#print "OUM",best_savings,biobjective,objective,objective_overload,
			update_district_info()
			#print "->",biobjective,objective,objective_overload
			find_best_solution += 1#update_best_solution()
			break
	if spatial_contiguity>=1: repair_fragmented_solution()
	time_op[0]+=time.time()-t
	return find_best_solution

#move one edge unit to its neighbor region
def one_unit_move_geo():
	#global district_info
	if location_problem==3: #pdp or homoDP
		return one_unit_move_pmp_geo()
	global node_groups
	global time_op
	global count_op
	popa=0.0
	if location_problem==2:
		popa=total_pop*1.0/max_num_facility #location_problem==2
	dev=pop_deviation*popa
	
	t=time.time()
	improve = 0
	ulist=find_edge_units()
	find_best_solution=0
	for n1 in ulist:
		k1 = node_groups[n1]
		for k2 in neighbor_regions(n1):
			obj_new = nodedik[n1][k2]
			obj_old = nodedik[n1][k1]
			old_bal=district_info[k1][4]+district_info[k2][4]
			new_pop1=district_info[k1][1]-nodes[n1][3]
			new_pop2=district_info[k2][1]+nodes[n1][3]
			new_bal=max(0,new_pop1-facilityCapacity[k1]) + max(0,new_pop2-facilityCapacity[k2])
			if location_problem==2: 
				new_bal=0.0
				if new_pop1 >popa+dev: new_bal+=new_pop1-popa-dev
				if new_pop2 >popa+dev: new_bal+=new_pop2-popa-dev
				if new_pop1 <popa-dev: new_bal+=popa-dev-new_pop1
				if new_pop2 <popa-dev: new_bal+=popa-dev-new_pop2
			#if isbetter(obj_old,obj_new,old_bal,new_bal)==0:
			if isbetter(obj_old,obj_new,objective_overload,objective_overload-old_bal+new_bal)==0:
				continue
			sol=node_groups[:]
			sol[n1] = k2
			if spatial_contiguity>=1 and check_continuality_feasibility(sol,k1)==0:
				break
			count_op[0]+=1
			node_groups[n1] = k2
			obj=biobjective
			objb=best_objective
			update_district_info()
			find_best_solution += 1# update_best_solution()
			break
	time_op[0]+=time.time()-t
	return find_best_solution
def one_unit_move_pmp_geo():
	#global district_info
	global node_groups
	global time_op
	global count_op
	t=time.time()
	improve = 0
	ulist=find_edge_units()
	find_best_solution=0
	for n1 in ulist:
		k1 = node_groups[n1]
		for k2 in neighbor_regions(n1):
			if nodedik[n1][k2] >= nodedik[n1][k1]: continue
			sol=node_groups[:]
			sol[n1] = k2
			#if spatial_contiguity>=1 and check_continuality_feasibility(sol,k1)<=0: break
			count_op[0]+=1
			node_groups[n1] = k2
			update_district_info()
			#find_best_solution += update_best_solution()
			#break
	time_op[0]+=time.time()-t
	if spatial_contiguity>=0:
		repair_fragmented_solution()
	return find_best_solution
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
				#	 ulist2=[x for x in range(num_units) if node_groups[x]==k and x!=u]
				#	 if ulist2==[]: continue
				#	 if check_ulist_continuality(ulist2)==0: continue
				node_groups[u] = k2
				improve+=1
				improved+=1
		if improve == 0: break
	time_op[0]+=time.time()-t
	#print "improved",improved,
	update_district_info()
	return improved
def two_unit_move():
	if spatial_contiguity==1 or location_problem==2:
		return two_unit_move_geo()
	global node_groups
	global time_op
	global count_op
	t=time.time()
	find_best_solution=0
	improve = 0
	ulist=find_edge_units()
	if ulist==[]:
		ulist=range(num_units)
	difflist=[]
	for x in ulist:
		rlist=[r for r in NearFacilityList[x] if centersID[r]>=0]
		k=node_groups[x]
		difflist.append([x,nodedij[x][rlist[0]]-nodedij[x][k]])
	if len(difflist)==0:
		print "debug: two_unit_move(): len(difflist)==0"
	difflist.sort(key=lambda x:x[1])
	ulist=[x[0] for x in difflist]
	nu=len(ulist)/2
	fulist=ulist[:nu]
	movelist =[]
	for n1 in ulist:
		k1 = node_groups[n1]
		rlist1=neighbor_regions(n1)
		success_move=0
		ulist2=[x for x in ulist if node_groups[x] in rlist1]
		random.shuffle(ulist2)
		for n2 in ulist2:
			#if n1 == n2: continue
			k2=node_groups[n2]
			#if k2 not in rlist1: continue
			#if district_info[k2][1]+nodes[n1][3]-nodes[n2][3]-facilityCapacity[k2] >district_info[k2][4]: continue
			for k3 in neighbor_regions(n2):
				#if k3!=k1 and district_info[k3][1]+nodes[n2][3]>facilityCapacity[k3]: continue
				#if k3==k1 and district_info[k3][1]-nodes[n1][3]+nodes[n2][3]-facilityCapacity[k3]>district_info[k3][4]: continue
				savings=nodedik[n1][k1]+nodedik[n2][k2]
				savings-=nodedik[n1][k2]+nodedik[n2][k3]
				savings*=envy_objective_weight_of_travelcost

				if envy_service_objective==1:
					if nodedij[n1][k1]>envy_service_distance:
						savings+= nodes[n1][3]*envy_objective_weight*(nodedij[n1][k1]-envy_service_distance)*(nodedij[n1][k1]-envy_service_distance)
					if nodedij[n2][k2]>envy_service_distance:
						savings+= nodes[n2][3]*envy_objective_weight*(nodedij[n2][k2]-envy_service_distance)*(nodedij[n2][k2]-envy_service_distance)
					if nodedij[n1][k2]>envy_service_distance:
						savings-= nodes[n1][3]*envy_objective_weight*(nodedij[n1][k2]-envy_service_distance)*(nodedij[n1][k2]-envy_service_distance)
					if nodedij[n2][k3]>envy_service_distance:
						savings-= nodes[n2][3]*envy_objective_weight*(nodedij[n2][k3]-envy_service_distance)*(nodedij[n2][k3]-envy_service_distance)
				if envy_service_objective==4:
					if nodedij[n1][k1]>envy_service_distance:
						savings+= nodes[n1][3]*envy_objective_weight*(nodedij[n1][k1]-envy_service_distance)
					if nodedij[n2][k2]>envy_service_distance:
						savings+= nodes[n2][3]*envy_objective_weight*(nodedij[n2][k2]-envy_service_distance)
					if nodedij[n1][k2]>envy_service_distance:
						savings-= nodes[n1][3]*envy_objective_weight*(nodedij[n1][k2]-envy_service_distance)
					if nodedij[n2][k3]>envy_service_distance:
						savings-= nodes[n2][3]*envy_objective_weight*(nodedij[n2][k3]-envy_service_distance)

				savings2=0
				if k3!=k1:
					slist=[district_info[k1][1]-nodes[n1][3],district_info[k2][1]+nodes[n1][3]-nodes[n2][3],district_info[k3][1]+nodes[n2][3]]
					savings2=district_info[k1][4]+district_info[k2][4]+district_info[k3][4]
					savings2-=max(0,slist[0]-facilityCapacity[k1])+max(0,slist[1]-facilityCapacity[k2])+max(0,slist[2]-facilityCapacity[k3])
				if k3==k1:
					slist=[district_info[k1][1]-nodes[n1][3]+nodes[n2][3],district_info[k2][1]+nodes[n1][3]-nodes[n2][3]]
					savings2=district_info[k1][4]+district_info[k2][4]
					savings2-=max(0,slist[0]-facilityCapacity[k1])+max(0,slist[1]-facilityCapacity[k2])
				#if district_info[k1][4]==0 and savings<0: continue
				if savings2<0: continue
				if savings2==0 and savings<0: continue
				#if savings<=0: continue
				count_op[1]+=1
				sp=objective_supply
				node_groups[n1] = k2
				node_groups[n2] = k3
				obj=biobjective
				#print "TUM",biobjective,objective,objective_overload,
				sp=objective_supply
				update_district_info()
				#print "->",biobjective,objective,objective_overload,":",savings,savings2,biobjective-obj
				find_best_solution += 1#update_best_solution()
				success_move=1
				break
			if success_move==1: break			 
	if spatial_contiguity>=1: repair_fragmented_solution()
	time_op[1]+=time.time()-t
	return find_best_solution

#for a region
#move out one edge unit to its neighbor region, and
#move in one edge unit from its neighbor region
def two_unit_move_geo():
	#global district_info
	global node_groups
	global time_op
	global count_op
	popa=0.0
	if location_problem==2:
		popa=total_pop*1.0/max_num_facility #location_problem=2,PDP
	dev=popa*pop_deviation

	t=time.time()
	find_best_solution=0
	improve = 0
	ulist=find_edge_units()
	movelist=[]

	difflist=[]
	for x in ulist:
		rlist=[r for r in NearFacilityList[x] if centersID[r]>=0]
		k=node_groups[x]
		difflist.append([x,district_info[k][4]*100000+nodedij[x][rlist[0]]-nodedij[x][k]])
	difflist.sort(key=lambda x:x[1])
	ulist_n1=[x[0] for x in difflist]
	ulist_n1=ulist[:]
	random.shuffle(ulist_n1)
	for n_id1 in ulist_n1:
		#if n_id1 in movelist: continue
		r_id1 = node_groups[n_id1]
		rlist1=neighbor_regions(n_id1)
		success_move=0
		#sol=node_groups[:]
		#sol[n_id1] = -1
		#if spatial_contiguity==1 and check_continuality_feasibility(sol,r_id1)==0:
		#	 movelist.append(n_id1)
		#	 continue

		for n_id2 in ulist:
			if n_id1 == n_id2: continue
			#if n_id1 in movelist: break
			#if n_id2 in movelist: continue
			if node_groups[n_id2] not in rlist1: continue
			new_r_id1=node_groups[n_id2]
			r_id2 = node_groups[n_id2]
			#sol=node_groups[:]
			#sol[n_id2] = -1
			#if spatial_contiguity==1 and check_continuality_feasibility(sol,r_id2)==0:
			#	 movelist.append(n_id2)
			#	 continue
			success_move=0
			for new_r_id2 in neighbor_regions(n_id2):
				obj_new = nodedik[n_id1][new_r_id1] + nodedik[n_id2][new_r_id2]
				obj_old = nodedik[n_id1][r_id1]+nodedik[n_id2][r_id2]
				new_district_info = [x[1] for x in district_info]
				new_district_info[r_id1] -= nodes[n_id1][3]
				new_district_info[r_id2] -= nodes[n_id2][3]
				new_district_info[new_r_id1] += nodes[n_id1][3]
				new_district_info[new_r_id2] += nodes[n_id2][3]
				bal=0.0
				bal=sum(new_district_info[x]-facilityCapacity[x] for x in range(num_districts) if new_district_info[x]>facilityCapacity[x])
				if location_problem==2:
					bal=0.0
					for x in new_district_info:
						if x==0: continue
						if x >popa+dev: bal+=x-popa-dev
						if x <popa-dev: bal+=popa-dev-x
				#yfkong
				#if bal>objective_overload: continue
				if isbetter(obj_old,obj_new,objective_overload,bal)==0:
					continue

				sol=node_groups[:]
				sol[n_id1] = new_r_id1
				sol[n_id2] = new_r_id2
				if spatial_contiguity>=1 and check_continuality_feasibility(sol,r_id1)==0:
					#movelist.append(n_id1)
					success_move=1
					break
				if spatial_contiguity>=1 and check_continuality_feasibility(sol,r_id2)==0:
					#movelist.append(n_id2)
					break	 

				count_op[1]+=1
				node_groups[n_id1] = new_r_id1
				node_groups[n_id2] = new_r_id2
				#movelist.append(n_id1)
				#movelist.append(n_id2)
				obj=objective
				#print [2,obj_old,obj_new,objective_overload,bal,biobjective],
				update_district_info()
				#print biobjective
				#if objective-obj> obj_new-obj_old+0.01:
				#	  print "debug oum!!!"
				find_best_solution += 1#update_best_solution()
				success_move=1
				break
			if success_move==1: break
	time_op[1]+=time.time()-t
	return find_best_solution

# local search
def RRT_local_search():
	#if location_problem==3:
	#	 pmp_VND_search()
	#	 return
	global improved
	#global node_neighbors
	improved=0
	#operators=[0,1,2,3,4]
	operators=assignment_operators_selected[:]
	#if op_random == 1 and random.random()>0.5:
	#	 random.shuffle(operators)
	for op in operators:
		if op == 0:
			one_unit_move()
			#update_region_pool_all()
		if op == 1:
			two_unit_move()
			#update_region_pool_all()
	return

#local search with operator op
def vnd_op_search(op):
	#global node_neighbors
	#for x in node_neighbors:
	#	 random.shuffle(x)
	if op == 0:
		one_unit_move()
	if op == 1:
		two_unit_move()
	#if 2 in assignment_operators_selected: print "VND",op, biobjective
	return

#VND search
def VND_local_search():
	#if location_problem==3:
	#	 pmp_VND_search()
	#	 return
	global improved
	improved=0
	#operators=[0,1,2,3,4]
	operators=assignment_operators_selected[:]	  
	if len(operators)==0: return
	#if op_random == 1:
	#if op_random == 1 and random.random()>0.5:
	#	 random.shuffle(operators)
	obj=biobjective
	while 1:
		vnd_op_search(operators[0])
		if biobjective < obj-0.00001:
			obj=biobjective
			continue
		if len(operators)==1:break
		vnd_op_search(operators[1])
		if biobjective < obj-0.00001:
			obj=biobjective
			continue
		if len(operators)==2:break

		vnd_op_search(operators[2])
		if biobjective < obj-0.00001:
			obj=biobjective
			continue
		if len(operators)==3:break
		vnd_op_search(operators[3])
		if biobjective < obj-0.00001:
			obj=biobjective
			continue
		if len(operators)==4:break
		vnd_op_search(operators[4])
		if biobjective < obj-0.00001:
			obj=biobjective
			continue
		break
	return

#VNs search
def vns_op_search(op):
	#global node_neighbors
	global node_groups
	node_groups=best_solution[:]
	update_district_info()
	obj=biobjective
	shake(op*2+2)
	#assign_ruin_recreate()
	#VND_local_search()
	vnd_op_search(op)
	#update_region_pool_all()
	return obj-biobjective

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
	num_first_feature=4
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
	print
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
		wgt=[x[3] for x in nodes]
		for k in range(num_features):
			if weight_features[k]<0.0000001: continue
			attr=[x[k+4] for x in nodes]
			avg=sum(attr[x]*wgt[x] for x in range(num_units))/area_total
			var=sum(wgt[x]*(attr[x]-avg)*(attr[x]-avg) for x in range (num_units) )
			if var<0.0000001: continue
			var/=area_total
			dev=pow(var,0.5)
			print avg,dev
			for i in range(num_units):
				nodes_std[i][k]=(nodes[i][k+4]-avg)/dev*math.sqrt(weight_features[k])

	nodedik=[[0.0 for x in range(num_units)] for y in range(num_units)]
	#if spatial_contiguity==1:
	nodedij=[[0.0 for x in range(num_units)] for y in range(num_units)]
	#else: nodedij=nodedik
	#nodedij=nodedik
	#print "nodes_std",nodes_std
	if data_standardlize==1:
		print "outlier data:"
		n=0
		'''for j in range(num_features):
			for i in range(num_units):
				x=nodes_std[i][j]
				if x>3 or x<-3: 
					n+=1
					print [nodes[i][0],j, nodes[i][j+4]],
			print
		print "total outlier values (>3, <-3):",n'''
		print "outlier locations:"
		for i in range(num_units):
			s1=sum(abs(x) for x in nodes_std[i])
			s2=max(abs(x) for x in nodes_std[i])
			if s1>7.5 or s2>3.0 : 
				print nodes[i][0],nodes_std[i],nodes[i]
				n+=1
		print
		print "total outlier values (>3, <-3):",n

	print "calculating the matrix ..."
	for i in range(num_units):
		#nodedik[i][i]=0.0
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
			nodedij[i][j]=d
			nodedij[j][i]=d
		if i%10==0: print ".",
		#if i>0 and i%(num_units/100) == 0: print str(i)+"/"+str(num_units)
	print
	#print nodedij
	print "creating node neighbors ..."
	#if spatial_contiguity!=1:
	#nodedij=nodedik
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
	#print node_neighbors
	#create_facility_neighbors()
	#find_NearFacilityList(num_districts)
	#find_near_customer()
	potential_facilities=[x for x in range(num_districts)]
	print "data prepaird!!"
	total_pop=num_units
	return

def read_pmp_instance(pfile):
	global num_units
	global total_pop
	global total_supply
	global nodes
	global node_neighbors
	global nodedij
	global nodedik
	global centersID
	global facilityCandidate
	global facilityCapacity
	global facilityCost
	global num_facilityCandidate
	global num_districts
	global district_info
	global avg_dis_min
	global potential_facilities
	global geo_instance
	geo_instance=0
	node =[0,0,0,1,0,0,0,0,0,0]
	#school_nodes = []
	nodes = []
	f = open(pfile)
	line = f.readline() #I,J
	line=line[0:-1]
	items = line.split(' ')
	if len(items)==1:
		items = line.split('\t')
	num_units=int(items[0])
	num_districts=int(items[0])
	facilityCandidate=[x for x in range(num_districts)]
	facilityCapacity=[num_units for x in range(num_districts)]
	facilityCost=[0.0 for x in range(num_districts)]
	nodes=[node[:] for x in range(num_units) ]
	nodedik=[ [0.0 for x in range(num_districts)] for x in range(num_units) ]
	nodedij=[ [0.0 for x in range(num_districts)] for x in range(num_units) ]
	arcpy_print("M,N: "+str(num_districts)+" "+str(num_units))
	for i in range(num_districts):
		line = f.readline()
		line=line[0:-1]
		items = line.split('\t')
		for j in range(num_units):
			nodedik[j][i]=int(items[j])
			nodedij[j][i]=int(items[j])
	#for i in range(num_districts): print nodedik[i]
	line = f.readline()
	line = f.readline()
	line=line[0:-1]
	items = line.split('\t')
	for i in range(num_units):
		facilityCost[i]=int(items[i])
	f.close()

	centersID=facilityCandidate[:]
	total_pop=sum(x[3] for x in nodes)
	total_supply=sum(facilityCapacity)
	district_info = [[0,0.0,0.0,0.0,0.0] for x in range(num_districts)]
	avg_dis_min=1.0

	create_facility_neighbors()
	find_NearFacilityList(num_districts)
	#print "NearFacilityList:"
	#for i in range(num_districts): print NearFacilityList[i]
	find_near_customer()

	for i in range(num_units):
		node_neighbors.append(NearFacilityList[i][:8])

	potential_facilities=[x for x in range(num_districts)]
	s="total demand: "+ str(total_pop)
	arcpy_print(s)
	s="total supply: "+str(total_supply)
	arcpy_print(s)

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


def create_node_neighbors():
	global node_neighbors
	#rlist=[x for x in range(num_districts)]
	mindij=[[MAXNUMBER for x in range(num_units)] for y in range(num_units)]
	for i in range(num_units):
		for j in range(num_units):
			if j<=i: continue
			dlist=[nodedij[i][x]-nodedij[j][x] for x in range(num_districts)]
			d=sum(x*x for x in dlist)
			mindij[i][j]=d
			mindij[j][i]=d
	node_neighbors = [[]for x in range(num_units)]
	for i in range(num_units):
		dlist=[[x, mindij[i][x]] for x in range(num_units) ]
		dlist.sort(key=lambda x:x[1])
		nn=8
		if nn>num_units: nn=num_units
		nghrs=[dlist[x][0] for x in range(nn)]
		random.shuffle(nghrs) #debug
		node_neighbors[i]=nghrs[:]


def find_nearFacilityFacility():
	global nearFacilityFacility
	nearFacilityFacility=[[] for x in range(num_districts)]
	dkk=[sum(nodedik[x][k]*nodedik[x][k] for x in range(num_units)) for k in range(num_districts)]
	#dkk.sort(key=lambda x:x[1])
	#dlist=[x[0] for x in dkk]
	for k in range(num_districts):
		d=dkk[k]
		dk=[[x,dkk[x]-d] for x in range(num_districts)]
		dk.sort(key=lambda x:x[1])
		del dk[0]
		nearFacilityFacility[k]=[x[0] for x in dk]
	
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
		#max_num_facility=num_districts
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
	#	 nodedik=[[nodedij[y][facilityCandidate[x]]*nodes[y][3] for x in range(num_districts)] for y in range(num_units)]
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
		fdlist=[ [x,nodedij[i][x]] for x in range(num_districts)]
		fdlist.sort(key=lambda x:x[1])
		flist=[x[0] for x in fdlist[0:n]]
		NearFacilityList.append(flist[:])
	print "data samples"
	for i in range(10):
		u=(i*17)%num_units
		print u,NearFacilityList[u][:10]
	if geo_instance==0:
		return
			

TB_tabu_list=[]
def pmp_TB():
	global node_groups
	global centersID
	best_exch=[] #klist,save
	best_obj=biobjective
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
			new_centersID[nk]=facilityCandidate[nk]
			for i in range(num_units):
				j=node_groups[i]
				if j!=k:
					obj+=min(nodedik[i][j],nodedik[i][nk])
				else:
					for j in NearFacilityList[i]:
						if new_centersID[j]>=0:
							obj+=nodedik[i][j]
							break
			if envy_service_objective==1:
				for i in range(num_units):
					for j in NearFacilityList[i]:
						if new_centersID[j]>=0:
							d=nodedij[i][j]
							if d>envy_service_distance:
								obj+= nodes[i][3]*envy_objective_weight*(d-envy_service_distance)*(d-envy_service_distance)
							break
			if envy_service_objective==4:
				for i in range(num_units):
					for j in NearFacilityList[i]:
						if new_centersID[j]>=0:
							d=nodedij[i][j]
							if d>envy_service_distance:
								obj+= nodes[i][3]*envy_objective_weight*(d-envy_service_distance)
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
		return 0
	#print best_exch,best_obj,biobjective
	best_exch.sort(key=lambda x:x[2])
	k=best_exch[0][0]
	centersID[k]=-1
	nk=best_exch[0][1]
	centersID[nk]=facilityCandidate[nk]
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
	#if centersID not in TB_tabu_list: TB_tabu_list.append(centersID[:])
	update_district_info()
	if spatial_contiguity>=1: repair_fragmented_solution()
	#update_centers()
	#print "tb",obj,biobjective,best_exch
	#if centersID not in TB_tabu_list: TB_tabu_list.append(centersID[:])
	return 1

def pmp_find(fi,facility2): #find fr for fi
	if envy_service_objective==1:
		return pmp_envy_find(fi,facility2)
	if envy_service_objective==4:
		return pmp_envy_find(fi,facility2)
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
	#return pmp_TB()
	global node_groups
	global centersID
	global time_Whitaker
	#global TB_tabu_list
	#if centersID in TB_tabu_list: return -1
	t=time.time()
	facility2=[0 for x in range(num_units)] #num_districts
	for i in range(num_units):
		k1=node_groups[i]
		for k2 in NearFacilityList[i]:
			if centersID[k2]>=0 and k2!=k1: 
				facility2[i]=k2
				break
	klist=[x for x in range(num_districts) if centersID[x]<0]
	random.shuffle(klist)
	improved=0
	best_swap=[-1,-1,0.0]
	swap=[]
	for nk in klist:
		k,profit=pmp_find(nk,facility2)
		if profit<=0: continue
		swap.append([k,nk,profit])
		if profit>best_swap[2]:
			best_swap[0]=k
			best_swap[1]=nk
			best_swap[2]=profit
		if len(swap)>=3:
			break
	swap.sort(key=lambda x:-x[2])
	if best_swap[0]>=0:
		#k=best_swap[0]
		#nk=best_swap[1]
		idx=0
		if len(swap)>2:
			idx=random.randint(0,1)
		#idx=0
		k=swap[idx][0]
		nk=swap[idx][1]
		centersID[k]=-1
		centersID[nk]=facilityCandidate[nk]
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
	if spatial_contiguity==1:
		repair_fragmented_solution()
	#if improved==-1: 
	#	 if centersID not in TB_tabu_list: TB_tabu_list.append(centersID[:])
	#print len(swap),
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
	#	 if centersID not in TB_tabu_list: TB_tabu_list.append(centersID[:])
	return improved

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
		if centersID[k]<0: continue
		kn,ulist=update_center(k)
		if kn==-1:
			print "update center error:",[k,sum(1 for x in range(num_units) if node_groups[x]==k)]," ->", [kn,len(ulist)]
			print NearFacilityList[k][:10],ulist
			kn=k
		if location_problem==3 and centersID[kn]>=0 and k!=kn:
			print "update center error: new center is already used:",k,kn
			kn=k
		if len(ulist)==0:
			print "update center error: ulist is empty:",k,kn,ulist
			kn=k
		for x in ulist: sol[x]=kn
		centers.append(kn)
		centersID[k]=-1
	node_groups=sol[:]
	for k in centers:
		centersID[k]=facilityCandidate[k]
	if location_problem==3 and spatial_contiguity==0:
		for i in range(num_units):
			for k in NearFacilityList[i]:
				if centersID[k]>=0:
					node_groups[i]=k
					break
	obj=biobjective
	update_district_info()
	update_best_solution()
	#print "updatecenters",biobjective-obj
	time_update_centers+=time.time()-t

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
def update_center_pmp(k): #need debug for PMP with few cands
	ulist=[x for x in range(num_units) if node_groups[x]==k]
	if ulist==[]: return k,[]
	best_cost=MAXNUMBER
	best_center=-1
	if k not in ulist:
		for i in NearFacilityList[k]:
			if i not in ulist: continue
			cost=sum(nodedik[x][i] for x in ulist)
			if cost<best_cost:
				best_cost=cost
				best_center=i
		return best_center,ulist
		
	for i in NearFacilityList[k][:min(num_districts/10,100)]:
		if i not in ulist: continue
		cost=sum(nodedik[x][i] for x in ulist)
		if cost<best_cost:
			best_cost=cost
			best_center=i
	return best_center,ulist


def location_check(key):
	if -1 in node_groups:
		arcpy_print("debug: "+str(key)+ " unit(s) not assigned! ")
		#return -1
	rlist=list(set(node_groups))
	if -1 in rlist: 
		arcpy_print("debug: "+str(key)+ " demand not assigned"+str(rlist))
		#arcpy_print(str(node_groups))
		rlist.remove(-1)
	if len(rlist)>max_num_facility and adaptive_number_of_facilities==0:
		arcpy_print("debug: "+str(key)+ " too many facilities"+str(rlist))
		#return -1
	for k in range(num_districts):
		if k in rlist and centersID[k]==-1:
			arcpy_print("debug: "+str(key)+ " facilitiy not selected but used")
			#return -1
		if centersID[k]>=0 and k not in rlist:
			arcpy_print("debug: "+str(key)+ " facilitiy selected but not used")
			print k, district_info[k]
			print [x for x in centersID if x>=0]

			#return -1
		uid=centersID[k]
		if spatial_contiguity==1 and uid>-1 and node_groups[uid]!=k:
			arcpy_print("debug: "+str(key)+ " facilitiy unit assigned to other facility"+str(k))
			print k,uid, node_groups[uid]
	'''clist=[]
	for k in range(num_districts):
		if centersID[k]<0: continue
		u=sum (1 for x in range(num_units) if node_groups[x]==k)
		clist.append(u)
	print clist'''
	#return 1

	
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
# multistart wrighted k-medoids
#start_num: number of starts
#loop_num: number of loops 
#connective: constraints on spatial contiguity? 1 for yes and 0 for no 
def multi_w_k_medoids(start_num, loop_num,connective):
	global node_groups
	global avg_pop	
	global avg_dis_min
	avg_pop=total_pop*1.0/max_num_facility
	wk_bestobj=MAXNUMBER
	wk_centers=[]
	wk_sol=[]
	for i in range(start_num):
		centers=k_means()  #change centersxy
		obj,dis,newcenters,sol=weighted_k_medoids(centers,loop_num,connective)
		biobj=obj
		if biobj < wk_bestobj:
			wk_bestobj=biobj
			wk_centers=newcenters[:]
			wk_sol=sol[:]
		avg_dis_min=dis/total_pop
	return wk_bestobj,wk_centers, wk_sol

# return k centers by k-means
def k_means():
	global centersID
	global node_groups
	centers=[]
	centersxy=[[0.0,0.0] for x in range(max_num_facility)]
	sol=[-1 for x in range(num_units)]
	#random centers
	#for i in range(num_units):
	#	 if nodes[i][3]>=avg_pop:
	#		 centers.append(i)
	while 1:
		nid=random.randint(0,num_units-1)
		if nid not in centers:
			centers.append(nid)
		if len(centers)==max_num_facility:
			break
	for k in range(max_num_facility):
		cid=centers[k]
		centersxy[k][0]=nodes[cid][1]
		centersxy[k][1]=nodes[cid][2]
	loop=0
	distance_obj=MAXNUMBER
	#k-means
	while 1:
		loop+=1
		#print	centersxy
		sol = [-1 for x in range(num_units)]
		#random.shuffle(nodelist)
		total_d=0.0
		#assign node to center
		for i in range(num_units):
			cid=-1
			d=MAXNUMBER
			for k in range(max_num_facility):
				dx=nodes[i][1]-centersxy[k][0]
				dy=nodes[i][2]-centersxy[k][1]
				#dxy=pow(dx*dx+dy*dy,0.5)
				dxy=dx*dx+dy*dy
				if dxy<d:
					d=dxy
					cid=k
			sol[i]=cid
		for k in range(max_num_facility):
			klist=[x for x in range(num_units) if sol[x]==k]
			if len(klist)==0:
				continue
			dx=sum(nodes[x][1]*nodes[x][3] for x in klist)
			dy=sum(nodes[x][2]*nodes[x][3] for x in klist)
			dsum=sum(nodes[x][3] for x in klist)
			centersxy[k][0]=dx/dsum
			centersxy[k][1]=dy/dsum
		obj=0.0
		for i in range(num_units):
			k=sol[i]
			dx=nodes[i][1]-centersxy[k][0]
			dy=nodes[i][2]-centersxy[k][1]
			dxy=pow(dx*dx+dy*dy,0.5) * nodes[i][3]
			obj+=dxy
		if obj<distance_obj:
			distance_obj=obj
		else:
			break
		#print loop, distance_obj,obj
	centers=[]
	for k in range(max_num_facility):
		ulist=[x for x in range(num_units) if sol[x]==k]
		kdis=MAXNUMBER
		kcenter=-1
		for x in ulist:
			tmp=0.0
			for y in ulist:
				#tmp+=nodedij[y][x]*nodes[y][3]*(10+random.random())/10
				tmp+=nodedik[y][x]*(9.5+random.random())/10
			if tmp<kdis:
				kdis=tmp
				kcenter=x
		centers.append(kcenter)
	#print "k-means",distance_obj,centers
	centersID=[-1 for x in range(num_districts)]
	for x in centers: centersID[x]=x
	node_groups=[centers[x] for x in sol]
	update_district_info()
	return centers
def random_centers():
	global centersID
	centers=[]
	while 1:
		if len(centers)==max_num_facility: break
		k=random.randint(0,num_districts-1)
		if k not in centers: centers.append(k)
	for x in range(num_districts):centersID[x]=-1
	for x in centers: centersID[x]=facilityCandidate[x]
	for i in range(num_units):
		for k in NearFacilityList[i]:
			if centersID[k]>=0:
				node_groups[i]=k
				break
	update_district_info()
	return 1
		

# return k centers by k-means
def k_medoids():
	global centersID
	global node_groups
	centers=[]
	sol=[-1 for x in range(num_units)]
	sol2=[-1 for x in range(num_units)]
   
	#while 1:
	#	 nid=random.randint(0,num_units-1)
	#	 if nid not in centers:
	#		 centers.append(nid)
	#	 if len(centers)==max_num_facility:
	#		 break
	random_centers()
	centers=[x for x in range(num_districts) if centersID[x] >=0]
	loop=0
	distance_obj=MAXNUMBER
	#k-means
	while 1:
		loop+=1
		for i in range(num_units):
			for k in NearFacilityList[i]:
				if k in centers:
					sol[i]=k
					break
		obj=0.0
		for k in range(max_num_facility):
			ulist=[x for x in range(num_units) if sol[x]==centers[k]]
			if len(ulist)==0: continue
			clist=range(num_districts)
			cid=-1
			mindis=MAXNUMBER
			for i in clist:
				dsum=sum(nodedik[x][i] for x in ulist)
				if dsum<mindis:
					mindis=dsum
					cid=i
			centers[k]=cid
			obj+=mindis
			for i in ulist: sol2[i]=cid
		sol=sol2[:]
		if obj<distance_obj:
			distance_obj=obj
		else:
			break
		#print loop, distance_obj,obj,centers
	#print "k-means",distance_obj,centers
	centersID=[-1 for x in range(num_districts)]
	for x in centers: centersID[x]=facilityCandidate[x]
	for i in range(num_units):
		for k in NearFacilityList[i]:
			if centersID[k]>=0:
				node_groups[i]=k
				break
	update_district_info()
	#update_centers()
	return 1

def k_medoids_sampling(allulist,cids):
	centers=cids
	num=len(allulist)
	if centers==[]:
		while 1:
			k=random.randint(0,num_districts-1)
			if k not in centers:
				centers.append(k)
			if len(centers)==max_num_facility:
				break
	sol=[-1 for x in allulist]
	loop=0
	distance_obj=MAXNUMBER
	#k-means
	while 1:
		loop+=1
		for idx in range(len(allulist)):
			i=allulist[idx]
			for k in NearFacilityList[i]:
				if k in centers:
					sol[idx]=k
					break
		obj=0.0
		for k in centers:
			ulist=[allulist[x] for x in range(num) if sol[x]==k]
			if len(ulist)==0: continue
			mindis=sum(nodedik[x][k] for x in ulist)
			cid=k
			for i in range(num_districts):
				dsum=sum(nodedik[x][i] for x in ulist)
				if dsum<mindis:
					mindis=dsum
					cid=i
			idx=centers.index(k)
			centers[idx]=cid
			obj+=mindis
		#print loop,distance_obj,sol
		if obj<distance_obj:
			distance_obj=obj
		else:
			break
	centers=list(set(centers))
	return distance_obj,centers

#wrighted k-medoids
#start_num: number of starts
#loop_num: number of loops 
#connective: constraints on spatial contiguity? 1 for yes and 0 for no	  
def weighted_k_medoids(centers,loops,connective):
	global node_groups
	global centersID
	weight=[1.0 for x in range(max_num_facility)]
	gaps=[0.0 for x in range(max_num_facility)]
	kw_centers=centers[:]
	sol=[-1 for i in range(num_units)]
	loop=0
	distance_best=MAXNUMBER
	bal_best=MAXNUMBER
	obj_best=MAXNUMBER
	best_centers=[]
	best_sol=[]
	for k in range(max_num_facility):
		weight[k]=1.0
	while loop<loops:
		loop+=1
		sol = [-1 for x in range(num_units)]
		total_d=0.0
		#assign node to center
		for i in range(num_units):
			cid=-1
			d=MAXNUMBER
			for k in range(max_num_facility):
				dxy=nodedij[i][kw_centers[k]]*weight[k]
				if dxy<d:
					d=dxy
					cid=k
			sol[i]=cid
			total_d+=d
		#if connective:
		#	 repair_fragmented_solution()		 

		for k in range(max_num_facility):
			ulist=[x for x in range(num_units) if sol[x]==k]
			kdis=MAXNUMBER
			kcenter=-1
			for x in ulist:
				if x not in facilityCandidate: continue
				tmp=0.0
				for y in ulist:
					tmp+=nodedik[y][x]*(9.5+random.random())/10
				if tmp<kdis:
					kdis=tmp
					kcenter=x
			kw_centers[k]=kcenter
		
		dis=0.0
		for k in range(max_num_facility):
			ulist=[x for x in range(num_units) if sol[x]==k]
			#dis+=sum(nodedij[x][kw_centers[sol[x]]]*nodes[x][3] for x in ulist)
			dis+=sum(nodedik[x][kw_centers[sol[x]]] for x in ulist)
			pop=sum(nodes[x][3] for x in ulist)
			gaps[k]=pop
		gap=0.0
		for pop in gaps:
			if pop<avg_pop*(1-pop_deviation):
				gap+=avg_pop*(1-pop_deviation)-pop
			if pop>avg_pop*(1+pop_deviation):
				gap+=pop-avg_pop*(1+pop_deviation)
		obj=dis + dis/total_pop*gap*pop_dis_coeff
		if obj<obj_best:
			obj_best=obj
			bal_best=gap
			best_centers=kw_centers[:]
			best_sol=sol[:]
			distance_best=dis	
	  
		#update weights
		for k in range(max_num_facility):
			weight[k]*=pow(1-loop*1.0/loops,0.5)+pow(loop*1.0/loops,0.5)*(gaps[k])/avg_pop
		avgw=sum(x for x in weight)/max_num_facility
		for k in range(max_num_facility):
			weight[k]/=avgw
	return obj_best,distance_best,best_centers,best_sol

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
def r_r_perb_center_locations(): #for pmp only, debuging
	global node_groups
	global centersID
	old_centersID=centersID[:]
	while 1:
		u=random.randint(0,num_units-1)
		k=node_groups[u]
		if centersID[k]==-1: continue
		for nk in NearFacilityList[u]:
			if nk==k: continue
			if centersID[nk]<0:
				centersID[k]=-1
				centersID[nk]=facilityCandidate[nk]
				break
		diff=sum(1 for x in range (num_districts) if centersID[x]!=old_centersID[x])
		if diff>=2: break
	for i in range(num_units):
		for k in NearFacilityList[i]:
			if centersID[k]>=0:
				node_groups[i]=k
				break
	update_district_info()
	if spatial_contiguity==1:
		repair_fragmented_solution()


def MIP_location_search():
	global centersID
	old_ids=centersID[:]
	rlist=select_region(-1)
	ulist=[x for x in range(num_units) if node_groups[x] in rlist]
	sol=location_sub_model(rlist,ulist,0.0001)
	if len(sol)==0: return 0
	nlist=list(set(sol))
	nlist.sort()
	if -1 in nlist: nlist.remove(-1)
	for i in range(num_units):
		if sol[i]<0: continue
		node_groups[i]=sol[i]
	centers=set(node_groups)
	centersID=[-1 for x in range(num_districts)]
	for x in centers: centersID[x]=facilityCandidate[x]
	update_district_info()
	#if all_units_as_candadate_locations==1 or num_units==len(facilityCandidate): one_unit_move()
	return 1

def pmp_sampling_solution(psize,sample_size):
	global node_groups
	global centersID
	global all_solutions
	t=time.time()
	print "psize,sample_size",psize,sample_size
	idx=0
	while 1: #for idx in range(psize):
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
		for x in centers: 
			if node_groups[x]!=x:
				node_groups[x]=x
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
		location_check(0)
		if spatial_contiguity>=1:
			repair_fragmented_solution_large_instance()
		update_district_info()
		#update_centers()
		klist=[x for x in range(num_districts) if centersID[x]>=0]
		if len(klist)<max_num_facility: continue
		sta=1
		for k in klist:
			if k not in node_groups:
				sta=0
				break
		'''total_area=sum(x[3] for x in nodes)
		for k in range(num_districts):
			if centersID[k]<0: continue
			ulist=[x for x in range(num_units) if node_groups[x]==k]
			area=sum(nodes[x][3] for x in ulist)
			if area<total_area*spatial_contiguity_minimum_percentage/max_num_facility/10000:
				sta=0
				break				 '''
		if sta==0: continue		   
		location_check(0)

		all_solutions.append([objective,centersID[:],node_groups[:]])
		#update_best_solution()
		idx+=1
		print "init. solution",idx,biobjective,objective,objective_overload,time.time()-t
		#print [x for x in centersID if x>=0]
		if len(all_solutions)>=psize: break
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
		#	 candidate_centers_selector=candidate_group
		location_check(1)
		#sta=1
		assign_ruin_recreate(-1)
		location_check(2)
		update_centers()
		location_check(3)

		#if spatial_contiguity>=1:
		#	 repair_fragmented_solution_large_instance()
		#update_region_pool_all()
		#update_best_solution()
		#update_centers()
		for i in range(6):
			sta=pmp_Whitaker()
			#update_region_pool_all()
			#update_best_solution()
			if sta==0: break
		if spatial_contiguity>=1:
			repair_fragmented_solution_large_instance()
		location_check(4)
		update_region_pool_all()
		update_best_solution()
		location_check(5)
		if best_biobjective_global<best_obj: 
			print "*",
			delta=max_loops_solution_not_improved*(best_obj-best_biobjective_global)/best_biobjective_global*1000 #100==1%，1000=0.0%
			not_improved=max(0,not_improved-int(delta))
		else:	 print "#",

		all_solutions.append([biobjective,centersID[:],node_groups[:] ])
		all_solutions.sort(key=lambda x:x[0])
		while len(all_solutions)>2*multi_start_count:all_solutions.pop()
		klist=list(set(node_groups))
		n=len(klist)
		print "Search loop",loop,"best",best_biobjective_global,"current",obj_begin,"->",n,objective, len(all_solutions),not_improved, time.time()-ga_time
		if not_improved>=max_loops_solution_not_improved: break
		#if time.time()-t_ga > heuristic_time_limit:  break
	node_groups=best_solution_global
	centersID=best_centersID_global
	update_district_info()
	ga_time=time.time()-t_ga
	search_stat()
	homoDP_stat()
	print "solution connectivity ..."
	spatial_contiguity=2
	for k in range(num_districts):
		if centersID[k]<0:continue
		ulist=[x for x in range(num_units) if node_groups[x]==k]
		#fu=frag_unit_minority(ulist)
		sta=check_continuality_feasibility(node_groups,k)
		print "cluster",nodes[k][0],len(ulist),"connectivity:",sta,
		ulist=[x for x in range(num_units) if node_groups[x]==k]
		print district_with_multiparts(ulist)
	return [best_biobjective_global,best_objective_global,best_overload_global,centersID,best_solution_global]

def homoDP_stat():
	if type(nodes[0][4]) == str:
		sst=MAXNUMBER
		cid=-1
		for i in range(num_units):
			var=sum(nodedik[i][x]*nodedik[i][x] for x in range(num_units))
			if var<sst:
				sst=var
				cid=i
			#print [sst,var],
		sse=0.0
		for k in range(num_districts):
			if centersID[k]<0: continue
			for i in range(num_units):
				ulist=[x for x in range(num_units) if node_groups[x]==k]
				sse=sum(nodedik[i][x]*nodedik[i][x] for x in ulist)
		r2=(sst-sse)/sst
		print "sse, sst, r2",sse, sst, r2
		return
	attrs=len(nodes_std[0])
	print "-----------------Stats for each attribute--------------------"
	print "Num of attrs", attrs
	print "weights",weight_features,len(weight_features)
	print "Idx----Avg--------Dev---------SSE-------SST------R2-------F-test---------"
	effective_attr=0
	area_total=sum(x[3] for x in nodes)
	for idx in range(attrs):
		if weight_features[idx]<0.0000001: continue
		effective_attr+=1
		mean=sum(x[3]*x[idx+4] for x in nodes) / area_total
		sst=sum(x[3]*pow(x[idx+4]-mean,2) for x in nodes)
		if sst<0.0000001:
			print "Fld"+str(idx), sst
			continue
		dev=pow(sst/area_total,0.5)
		sse=0.0
		means=[0.0 for x in range(num_districts)]
		for k in range(num_districts):
			if centersID[k]<0:continue
			ulist=[x for x in range(num_units) if node_groups[x]==k]
			sumd=sum(nodes[x][3]*nodes[x][idx+4] for x in ulist)
			means[k]=sumd/sum(nodes[x][3] for x in ulist)
		for i in range(num_units):
			k=node_groups[i]
			d=nodes[i][idx+4]-means[k]
			sse+=nodes[i][3]*d*d
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
		mean=sum(nodes[x][3]*nodes_std[x][idx] for x in range(num_units) ) / area_total
		#sst=sum(pow(x[idx]-mean,2) for x in nodes_std)
		sst=sum(nodes[x][3]*pow(nodes_std[x][idx]-mean,2) for x in range(num_units))
		SST_all+=sst
		#if mean<0.00000000001:	 continue
		dev=pow(sst/num_units,0.5)

		sse=0.0
		means=[0.0 for x in range(num_districts)]
		for k in range(num_districts):
			if centersID[k]<0:continue
			ulist=[x for x in range(num_units) if node_groups[x]==k]
			sumd=sum(nodes[x][3]*nodes_std[x][idx] for x in ulist)
			means[k]=sumd/sum(nodes[x][3] for x in ulist)
		for i in range(num_units):
			k=node_groups[i]
			d=nodes_std[i][idx]-means[k]
			sse+=nodes[i][3]*d*d
		SSE_all+=sse
	print SST_all, SSE_all
	r2=(SST_all-SSE_all)/SST_all
	f=r2/(max_num_facility-1) * (num_units-max_num_facility) / (1-r2)
	print "Stat:",	"SSE=",SSE_all,"SST=",SST_all, "R2=",r2,"F-stat=", f

	print "-----------------------mean values of regions-------------------------"
	s="Region"
	klist=[x for x in centersID if x>=0]
	for k in klist:
		s+="\t"+str(nodes[k][0])
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
	#	 if centersID[k]<0:continue
	#	 ulist=[x for x in range(num_units) if node_groups[x]==k]
	#	 for idx in range(attrs):
	#		 sumd=sum(nodes[x][idx+4] for x in ulist)
	#		 mean=sumd/len(ulist)
	#		 d= mean-nodes[k][idx+4]
	#		 error+=d*d
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
	#if spatial_contiguity==1:
	#	 sta=check_solution_continuality_feasibility(best_solution_global)
	#	 arcpy_print("solution on continuality (0 no, 1 yes) : "+str(sta))
	arcpy_print("----------------end of search statistics----------------")
 
#method
#0 basic 2SFCA: 0, na,na
#1 gravity 2SFCA: 1,1.0,2.0
#2 Gaussian 2SFCA: 2,1.0,na
#3 Kernel Density 2SFCA, 3,1.0,na
#4 OSD
def spatial_accessibility(method,radius,beta):
	#global node_groups
	#global centersID
	#centersID=facilityCandidate
	#node_groups=[-1 for x in range(num_units)]
	t=time.time()
	rlist=[0.0 for x in range(num_districts)]
	for k in range(num_districts):
		dk=0.0
		for i in nearCustomers[k]:
			if nodedij[i][k]>radius: break
			if method==0: dk+=nodes[i][3]
			if method==1: dk+=nodes[i][3]/pow(nodedij[i][k],beta)
			if method==2: dk+=nodes[i][3]* ( math.exp(-nodedij[i][k]*nodedij[i][k]/radius/radius/2)-math.exp(-0.5) )/ ( 1-math.exp(-0.5) )
			if method==3: dk+=nodes[i][3]*(0.75 - 0.75*nodedij[i][k]*nodedij[i][k]/radius/radius ) 
		if dk>0.0000001: rlist[k]=facilityCapacity[k]/dk
	alist=[0.0 for x in range(num_units)]
	for i in range(num_units):
		a=0.0
		for k in NearFacilityList[i]:
			if nodedij[i][k]>radius: break
			if method==0: a+=rlist[k]
			if method==1: a+=rlist[k]/pow(nodedij[i][k],beta)
			if method==2: a+=rlist[k]* ( math.exp(-nodedij[i][k]*nodedij[i][k]/radius/radius/2)-math.exp(-0.5) )/ ( 1-math.exp(-0.5) )
			if method==3: a+=rlist[k]*(0.75 - 0.75*nodedij[i][k]*nodedij[i][k]/radius/radius) 
		alist[i]=a
	print "----------spatial_accessibility-------------"
	print "method, radius, beta", method,radius,beta
	print "time",time.time()-t
	print "Idx x y D ID S r1 r2 A Dis"
	for i in range(num_units):
		#k=NearFacilityList[i][0]
		if nodes[i][3]<=0: continue
		#node_groups[i]=k
		for x in nodes[i]: print x,
		k=node_groups[i]
		print alist[i],nodedij[i][k]
	print "time",time.time()-t
	print "----------stat-------------------"
	pop=[0,0,0,0,0]
	for i in range(num_units):
		a=alist[i]
		if a<0.6: 
			 pop[0]+=nodes[i][3]
		elif a<0.8: 
			 pop[1]+=nodes[i][3]
		elif a<1.2: 
			 pop[2]+=nodes[i][3]
		elif a<1.4: 
			 pop[3]+=nodes[i][3]
		else: 
			 pop[4]+=nodes[i][3]
	print "time",time.time()-t
	for x in pop: print x*10000/total_pop/100.0, 
	print "time",time.time()-t
	update_district_info()
	m=print_equality_measures()
	for x in m: print x,m[x]
def print_equality_measures():
	#print "-----------equality measures---------------"
	maxd=0.0
	mind=MAXNUMBER
	maxdev=0.0
	avgd=objective/total_pop
	absdev=0.0
	stddev=0.0
	Theil=0.0
	schutz=0.0
	upperVar=0.0
	lowerVar=0.0
	for i in range(num_units): 
		k=node_groups[i]
		dis=nodedij[i][k]
		w=nodes[i][3]
		absdev+=w*abs(dis-avgd)
		stddev+=w*(dis-avgd)*(dis-avgd)
		if dis>maxd: maxd=dis
		if dis<mind: mind=dis
		if abs(dis-avgd)>maxdev: maxdev=abs(dis-avgd)
		'''if dis>0:
			a=dis*math.log(dis)-avgd*math.log(avgd)
			Theil+=w*a*a
		else:
			a=avgd*math.log(avgd)
			Theil+=w*a*a'''
		if dis>0:
			a=math.log(dis/avgd)
			Theil+=w*dis*a
		schutz+=w*abs(dis-avgd)
		if dis>avgd:
			upperVar+=w*(dis-avgd)*(dis-avgd)
		else:
			lowerVar+=w*(dis-avgd)*(dis-avgd)
	equality_measures={}
	equality_measures["Mean"]=avgd
	#print "Centre", maxd
	equality_measures["Centre"]= maxd
	#print "Range",maxd-mind
	equality_measures["Range"]=maxd-mind
	#print "MaxDev", maxdev
	equality_measures["MaxDev"]=maxdev
	#print "MeanDev",absdev/total_pop
	equality_measures["MeanDev"]=absdev/total_pop
	#print "StdDev",math.sqrt(stddev/total_pop) #d.num_units
	equality_measures["Varance"]=stddev/total_pop
	equality_measures["StdDev"]=math.sqrt(stddev/total_pop)
	#print "Theil", Theil/avgd/total_pop #d.num_units
	equality_measures["VC"]=equality_measures["StdDev"] /avgd
   
	equality_measures["Theil"]= Theil/avgd/total_pop
	Gini=0.0
	for i in range(num_units): 
		k=node_groups[i]
		d1=nodedij[i][k]
		w1=nodes[i][3]
		for j in range(num_units): 
		   k=node_groups[j]
		   d2=nodedij[j][k]
		   w2=nodes[j][3]
		   Gini+=w1*w2*abs(d1-d2)
	#print "Gini",Gini/total_pop/total_pop/2/avgd #Gini/d.num_units/d.num_units/2/avgd
	equality_measures["Gini"]=Gini/total_pop/total_pop/2/avgd
	equality_measures["Schutz"]=schutz/total_pop/2/avgd
	equality_measures["lowerVar"]=lowerVar
	equality_measures["upperVar"]=upperVar
	return equality_measures

def coverage_stat():
	maxd=0.0
	for i in range(num_units): 
		k=node_groups[i]
		dis=nodedij[i][k]
		if dis>maxd: maxd=dis
	covered=[0 for x in range(20)]
	dlist=[0.5*x for x in range(20)]
	for i in range(num_units):
		k=node_groups[i]
		dis=nodedij[i][k]
		demand=nodes[i][3]
		for j in range(20):
			if dis<=dlist[j]: covered[j]+=demand
	#print dlist
	#print covered
	#print [x*100.0/total_pop for x in covered]
	for i in range(len(dlist)):
		#if covered==total_pop: break
		print dlist[i],(covered[i]*10000/total_pop)/100.0
		if dlist[i]>maxd: break

def PMP_and_minvar_2sfca(num_f,radius,supply_all):
	#pmp_mip_model(num_f,radius,cover_percent,time_limit,mipgap)
	#pmp_mip_model(num_f,20,100,1000,0.000001)
	#Li et al. 2022, CEUS
	facilityList=[x for x in range(num_districts) if centersID[x]>=0]
	print "facility locations:",facilityList
	print "facility idx",[x for x in centersID if x>=0]
	fik=[ [0.0 for x in range(num_f) ] for y in range(num_units)]
	for i in range(num_units):
		for k in range(num_f):
			j=facilityList[k]
			d=nodedij[i][j]
			if d>radius: continue
			fik[i][k]= ( math.exp(-d*d/radius/radius/2)-math.exp(-0.5) )/( 1-math.exp(-0.5))
			#print [i,k,fik[i][k]],
	#print "fik", fik[500:700]

	glist=[100000.0 for x in range(num_f)]
	for k in range(num_f):
		g=sum(nodes[i][3]* fik[i][k] for i in range(num_units))
		glist[k]=1.0/g	  
	#print "glist",glist
	
	dlist=[nodes[x][3] for x in range(num_units)]
	#print "dlist",dlist
	
	#P=FG, H=PDP, Q=-PDA
	pik=[ [0.0 for x in range(num_f) ] for y in range(num_units)]
	for i in range(num_units):
		for k in range(num_f):
			pik[i][k]=fik[i][k]*glist[k]
	#print "pik",pik[:200]

	hlist=[0.0 for k in range(num_f)]
	for k in range(num_f):
		hlist[k]=sum( pik[i][k] * pik[i][k] *dlist[i] for i in range(num_units))  
	#print "hlist",hlist

	qlist=[0.0 for k in range(num_f)]
	a=supply_all*1.0/total_pop
	for k in range(num_f):
		qlist[k]=sum( a*dlist[i]*pik[i][k] for i in range(num_units))
	#print "qlist",qlist

	f = open("pmp_2sfca_"+str(num_f)+"_"+str(radius)+".lp","w")
	f.write("minimize\nobj: ")
	
	s=""
	f.write(" [ ")
	for k in range(num_f):
		s += str(0.5*hlist[k])+ " x_" +str(k)+" ^ 2 + " 
	f.write( s[:-3]+ " ] / 2 ")

	s=""
	for k in range(num_f):
		s+= " - " + str(qlist[k])+ " x_" +str(k) 
	f.write( s +"\n")
	f.write("subject to\n")
	s="_c1: "
	for k in range(num_f):
		s += " x_" +str(k)+" + "
	f.write(s[:-3]+ " <= " + str(supply_all)+"\n")
	f.write("generals\n")
	for k in range(num_f):
		s = " x_" +str(k)+"\n"
		f.write(s)
	f.write("end\n")
	f.close()	 

def sol_stat():
	sol={}
	sol["num_facility"]=sum(1 for x in centersID if x>=0)
	sol["objective"]=biobjective
	sol["TotalCost"]=objective+objective_fcost
	sol["DCost"]=objective
	sol["Fcost"]=objective_fcost
	sol["WeightedEnvyValue"]=objective_envy
	if objective_envy>0:
		sol["EnvyValue"]=objective_envy/envy_objective_weight
	else:
		envy=0.0
		for i in range(num_units):
			k=node_groups[i]
			d=nodedij[i][k]
			if d>envy_service_distance:
				envy+=nodes[i][3]*(d-envy_service_distance)*(d-envy_service_distance)
		sol["EnvyValue"]=envy
	sol["OverloadPenalty"]=objective_pentalty_on_overload
	sol["RMaxPenalty"]=objective_pentalty_on_rmax
	sol["RCoverPenalty"]=objective_pentalty_on_covering
	sol["Penalty"]=penalty_on_demand
	return sol
  
def check_sub(dlist,ulist,centers,exist_facility):
	for k in dlist:
		if centersID[k]<0: print k,"in dlist not used"
		if k not in centers: print k,"in dlist not in centers"
		if k in exist_facility: print k,"in dlist also in exist_facility"
	for i in ulist:
		if node_groups[i] not in dlist: print i, "in ulist not belongs dlist"
	for k in centers:
		if k not in dlist and centersID[k]>=0: print k, "in centers is open"
	for k in exist_facility:
		if centersID[k]<0: print k,"in exist_facility not used"
		if k in dlist: print k,"in exist_facility also in dlist"
	
	


