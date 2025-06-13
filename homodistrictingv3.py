# -*- coding: utf-8 -*-
## ILS algorithm for regionalization Problems, version 2
## yfkong@henu.edu.cn, Apr.,2024
 
import sys,os,random,time,copy,math, keyboard
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
attributes=[]

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
objective_R2=0.0
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
time_check=0.0
time_check_edge_unit=0
time_spp=0.0
time_update_centers=0.0
time_op=[0.0 for x in range(10)]
time_ruin_recreate=[0.0 for x in range(10)]
time_location=[0.0 for x in range(10)]
time_pmp_re_location=0.0
time_Whitaker=0.0
time_repair=0.0
time_Renato_update=0.0
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



stop_searching=0
def on_key_event(e):
	global stop_searching
	#print e.name    
	if e.name=="esc":
		stop_searching=1    
		print "user break...!"
keyboard.on_press_key('esc', on_key_event)
        
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
		ulist1=[y for y in ulist1 if y>=0]
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
	if len(ulist)<2: return []
	final_list=[]
	ulist1=ulist[:]
	if k in ulist1:
		ulist1.remove(k)
		ulist2=[k]
	else:
		ulist2=[ulist1[-1]]
		ulist.pop()
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
	while 1:
		assigned=0
		for i in range(len(frag_units)):
			u=frag_units[i]
			for v in node_neighbors[u]:
				k=node_groups[v]
				if k>=0: 
					assigned=1
					node_groups[u]=k
					frag_units[i]=-1
					break                    
		frag_units=[x for x in frag_units if x>=0]
		if len(frag_units)==0: break
		if assigned==0: break    
    
	'''while len(frag_units)>0:
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
			if nid in frag_units: frag_units.remove(x[0])'''
	for x in frag_units:
		#print len(frag_units)
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
	if num_units/max_num_facility>100:
		 repair_fragmented_solution_large_instance()
		 return
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
	time_repair+=time.time()-t    

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
	global attributes
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
	items=line.split("\t")
	#line=line[:-1]
	attributes=items[num_first_feature:num_first_feature+num_features]
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
	total_pop=sum(x[3] for x in nodes)
	facilityCandidate=range(num_districts)
	facilityCapacity=[total_pop for x in range(num_districts)]
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
		print "attr, mean and stdev:"
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
			print attributes[k],avg,dev
			for i in range(num_units):
				nodes_std[i][k]=(nodes[i][k+4]-avg)/dev*math.sqrt(weight_features[k])

	nodedik=[[0.0 for x in range(num_units)] for y in range(num_units)]
	nodedij=[[0.0 for x in range(num_units)] for y in range(num_units)]
	#else: nodedij=nodedik
	#nodedij=nodedik
	#print "nodes_std",nodes_std
	if data_standardlize==10:
		print "outlier data:"
		n=0
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
		if num_units>1000 and i%100==0: print ".",
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
				if idx1>=0 and idx2>=0:
					break
			if idx1>=0 and idx2>=0:
				node_neighbors[idx1].append(idx2)
		line = f.readline()
	f.close()
	potential_facilities=[x for x in range(num_districts)]
	print "data prepaird!!"
	return

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
		fdlist[i][1]=-1.0
		fdlist.sort(key=lambda x:x[1])
		flist=[x[0] for x in fdlist[0:n]]
		NearFacilityList.append(flist[:])
	for k in range(num_districts):
		u=facilityCandidate[k]
		if NearFacilityList[u][0]!=k: 
			print "near unit error",u,NearFacilityList[u][0]
			
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
	return improved


def Renato_update(Facility2List,gainList,lossList,extraList): 
	global time_Renato_update    
	t=time.time()
	for i in range(num_units):
		k1=node_groups[i]
		for k2 in NearFacilityList[i]:
			if centersID[k2]<0: continue
			if k2==k1: continue
			Facility2List[i]=k2
			break
	for k in range(num_districts):
		gainList[k]=0.0
		if centersID[k]>=0:continue
		for i in range(num_units):
			k1=node_groups[i]
			gainList[k] += max(0, nodedik[i][k1]-nodedik[i][k])
	for k in range(num_districts):
		lossList[k]=0.0
		if centersID[k]<0:continue
		ulist=[x for x in range(num_units) if node_groups[x]==k]
		for i in ulist:
			k2=Facility2List[i]
			lossList[k]+=nodedik[i][k2]-nodedik[i][k]
	for ki in range(num_districts):
		for kr in range(num_districts):
			extraList[ki][kr]=0.0
			if centersID[ki]>=0:continue
			if centersID[kr]<0:continue
			ulist=[x for x in range(num_units) if node_groups[x]==kr]
			for i in ulist:
				k2=Facility2List[i]
				if nodedij[i][ki]<nodedij[i][k2]:
					extraList[ki][kr]+=nodedik[i][k2]-max(nodedik[i][ki],nodedik[i][kr])
	time_Renato_update+= time.time()-t    

#On the implementation of a swap-based local search procedure for the p-median problem
#Renato&Werneck,2007
def pmp_Renato(Facility2List,gainList,lossList,extraList): #bug
	#return pmp_TB()
	global node_groups
	global centersID
	global time_Whitaker
	obj0=biobjective	
	t=time.time()
	kiList=[x for x in range(num_districts) if centersID[x]<0]
	krList=[x for x in range(num_districts) if centersID[x]>=0]
	random.shuffle(kiList)
	random.shuffle(krList)
	improved=0
	bestSwap=[]
	bestprofit=0.0
	bestki=-1
	bestkr=-1
	for kr in krList:
		for ki in kiList:
			profit=gainList[ki]-lossList[kr]+extraList[ki][kr]
			if profit>0.0001: bestSwap.append([ki,kr,profit])
	if bestSwap==[]:return improved
	bestSwap.sort(key=lambda x:-x[2])
	idx=0
	#if len(bestSwap)>3:
	#	idx=random.randint(0,2) #elite improvement, select on from top three best swaps
	bestki=bestSwap[idx][0]
	bestkr=bestSwap[idx][1]
	bestprofit=bestSwap[idx][2]
   
	#print "solution update",bestki,bestkr,bestprofit,gainList[bestki],lossList[bestkr],extraList[bestki][bestkr],facilityCost[bestkr]-facilityCost[bestki],		
	# undo  gainList,lossList,extraList
	alist=[i for i in range(num_units) if node_groups[i]==bestkr or Facility2List[i]==bestkr or nodedij[i][bestki]<nodedij[i][Facility2List[i]] ]
	for i in alist:
		k1=node_groups[i]
		k2=Facility2List[i]
		d1=nodedik[i][k1]
		d2=nodedik[i][k2]
		lossList[k1]-=d2-d1
		for k in range(num_districts):
			if centersID[k]>=0: continue
			d=nodedik[i][k]
			if d<d2:
				gainList[k]-=max(0,d1-d)
				extraList[k][k1]-=d2-max(d,d1)
	#update solution
	centersID[bestkr]=-1
	centersID[bestki]=facilityCandidate[bestki]
	for i in range(num_units):
		for k in NearFacilityList[i]:
			if centersID[k]>=0:
				node_groups[i]=k
				break
		k=node_groups[i]
		if k!=bestkr:
			if nodedij[i][k]>nodedij[i][bestki]: 
				node_groups[i]=bestki
		else:
			k=Facility2List[i]
			if nodedij[i][bestki]<nodedij[i][k]:
				k=bestki
			node_groups[i]=k
	obj=biobjective
	update_district_info()

	#update  Facility2List,gainList,lossList,extraList
	for i in range(num_units):
		k1=node_groups[i]
		for k2 in NearFacilityList[i]:
			if centersID[k2]<0: continue
			if k2==k1: continue
			Facility2List[i]=k2
			break
	for i in alist:
		k1=node_groups[i]
		k2=Facility2List[i]
		d1=nodedik[i][k1]
		d2=nodedik[i][k2]
		lossList[k1]+=d2-d1
		for k in range(num_districts):
			if centersID[k]>=0: continue
			d=nodedik[i][k]
			if d<d2:
				gainList[k]+=max(0,d1-d)
				extraList[k][k1]+=d2-max(d,d1)
	#Renato_update(Facility2List,gainList,lossList,extraList)
	update_district_info()
	#print [bestprofit,obj-biobjective],
	improved=1
	time_Whitaker+=time.time()-t
	#if spatial_contiguity==1:
	#	repair_fragmented_solution()
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
		update_centers()
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
		update_best_solution()        
		all_solutions.append([objective,centersID[:],node_groups[:]])
		#update_best_solution()
		idx+=1
		print "init. solution",idx,biobjective,objective,objective_overload,time.time()-t
		#print [x for x in centersID if x>=0]
		if len(all_solutions)>=psize: break
	#return all_solutions

def pmp_VND_search():
    global node_groups
    global time_op
    global count_op
    improved = 0
    t=time.time()
    if spatial_contiguity==0: return 0
    while 1:
        improve = 0
        ulist=find_edge_units()
        for u in ulist:
            k=node_groups[u]
            for u2 in node_neighbors[u]:
                if node_groups[u2]==k: continue
                k2=node_groups[u2]
                if nodedik[u][k2]>=nodedik[u][k]: continue
                node_groups[u] = k2
                improve+=1
                improved+=1
        if improve == 0: break
    time_op[0]+=time.time()-t
    update_district_info()
    return improved
    
def ils_large_homoDP_fast(numf):
	global best_biobjective_global
	global node_groups
	global centersID
	global district_info
	global region_pool
	global pool_index
	global all_solutions
	global max_num_facility
	max_num_facility=numf
	#initialize_instance()
	find_NearFacilityList(num_units)
	find_near_customer()
	#for x in NearFacilityList: print x[:10]
	node_groups=[-1 for x in range(num_units)]
	centersID=[-1 for x in range(num_districts)]
	all_solutions=[]

	region_pool=[]
	ga_time=time.time()
	best_biobjective_global = MAXNUMBER
	district_info = [[0,0.0,0,0,0] for x in range(num_districts)]
	population=[] #all
	pool_index=[]
	if is_spp_modeling>=1:
		pool_index=[[] for x in range(num_districts)]
	t_ga=time.time()

	print "\nstep 1: generating initial solutions ..."
	#sampling and solving sub-problems
	sample_size=max(100*max_num_facility, num_units/10)
	if sample_size>num_units: sample_size=num_units
	num_samplling=multi_start_count
	if num_samplling<=0:
		num_samplling=3*num_units/sample_size
	#if num_samplling<10*max_num_facility: num_samplling=10*max_num_facility
	#if num_samplling<multi_start_count: num_samplling=multi_start_count
	pmp_sampling_solution(num_samplling,sample_size)

	print "step 2: solution improvement..."
	print "!!!!! Press Esc key to stop searching !!!!!"

	all_solutions.sort(key=lambda x:x[0])
	not_improved=0
	loop=0

	Facility2List=[-1 for x in range(num_units)]
	gainList=[0.0 for x in range(num_districts)] #gian(fi)
	lossList=[0.0 for x in range(num_districts)] #loss(fr)
	extraList=[[0.0 for x in range(num_districts)] for y in range(num_districts)] #extra( fi, fr)

	while 1:
		loop+=1
		r=random.random()
		sidx = int(min(multi_start_count,len(all_solutions))* pow(r,1.2)*0.999)
		node_groups=all_solutions[sidx][2][:]
		centersID=all_solutions[sidx][1][:]
		update_district_info()
		not_improved+=1
		best_obj_begin=best_biobjective_global
		obj_begin=biobjective
		location_check(1)
		#assign_ruin_recreate(-1)
		maxStrength=min(6, max_num_facility*2/3)
		strength=random.randint(3,maxStrength)
		clist=[x for x in range(num_districts) if centersID[x]<0]
		random.shuffle(clist)
		for x in range(strength):
			k=clist[x]
			centersID[k]=facilityCandidate[k]
		clist=[x for x in range(num_districts) if centersID[x]>=0]
		random.shuffle(clist)
		for x in range(strength):
			k=clist[x]
			centersID[k]=-1
		for i in range(num_units):
			for k in NearFacilityList[i]:
				if centersID[k]>=0:
					node_groups[i]=k
					break
		update_district_info()

		location_check(2)
		Renato_update(Facility2List,gainList,lossList,extraList)
		while 1:
			sta=pmp_Renato(Facility2List,gainList,lossList,extraList)#bug
			#sta=pmp_Whitaker()
			update_district_info()
			#print int(biobjective),
			print ".",
			if sta==0: break
			if stop_searching==1: break
		print ""
		'''best_sol=node_groups[:]
		best_cid=centersID[:]
		best_obj=biobjective
		while 1:
			sta=pmp_Renato(Facility2List,gainList,lossList,extraList)#bug
			#sta=pmp_Whitaker()
			update_district_info()
			print int(biobjective),
			if biobjective<	best_obj-0.00001:		
				best_sol=node_groups[:]
				best_cid=centersID[:]
				best_obj=biobjective
			else:
				break    
			if sta==0: break
		node_groups=best_sol[:]
		centersID=best_cid[:]
		update_district_info()'''
		if spatial_contiguity>=1:
			repair_fragmented_solution()
		update_centers()
		update_best_solution()
		pmp_VND_search()            
		if spatial_contiguity>=1:
			repair_fragmented_solution()
		location_check(3)
		update_centers()

		update_region_pool_all()
		update_best_solution()
		all_solutions.append([biobjective,centersID[:],node_groups[:] ])
		#if spatial_contiguity>=1:
		#	repair_fragmented_solution_large_instance()
		#location_check(4)
		#update_region_pool_all()
		#update_best_solution()
		location_check(5)
		if best_biobjective_global<best_obj_begin-0.00001: 
			print "*",
			delta=max_loops_solution_not_improved*(best_obj_begin-best_biobjective_global)/best_biobjective_global*1000 #100==1%，1000=0.0%
			not_improved=max(0,not_improved-int(delta))
		else:	 print "#",

		all_solutions.append([biobjective,centersID[:],node_groups[:] ])
		all_solutions.sort(key=lambda x:x[0])
		while len(all_solutions)>2*multi_start_count:all_solutions.pop()
		klist=list(set(node_groups))
		n=len(klist)
		r2=1-best_biobjective_global/total_pop/sum(weight_features)
		print "Search loop",loop,"best",int(best_biobjective_global), int(r2*10000)/10000.0, "current",int(obj_begin),"->",n,int(biobjective),"stat", len(all_solutions),not_improved,int(time.time()-ga_time)
		if not_improved>=max_loops_solution_not_improved: break
		if stop_searching==1: break
		#if time.time()-t_ga > heuristic_time_limit:  break
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
		sta=check_continuality_feasibility(node_groups,k)
		print "cluster",nodes[k][0],len(ulist),"connectivity:",sta,
		ulist=[x for x in range(num_units) if node_groups[x]==k]
		print district_with_multiparts(ulist)
	return [best_biobjective_global,best_objective_global,best_overload_global,centersID,best_solution_global]

def homoDP_stat():
	global objective_R2
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
	print "Attr----Avg--------Dev---------SSE-------SST------R2-------F-test---------"
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
		print attributes[idx], mean, dev, sse,sst, r2, f
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
	#print SST_all, SSE_all
	r2=(SST_all-SSE_all)/SST_all
	f=r2/(max_num_facility-1) * (num_units-max_num_facility) / (1-r2)
	print "Stat:",	"SSE=",SSE_all,"SST=",SST_all, "R2=",r2,"F-stat=", f
	objective_R2=r2
	print "-----------------------mean values of regions-------------------------"
	s="Region"
	klist=[x for x in centersID if x>=0]
	for k in klist:
		s+="\t"+str(nodes[k][0])
	print s

	for idx in range(attrs):
		s=""
		if weight_features[idx]<0.0000001: continue
		s+=attributes[idx]
		for k in range(max_num_facility):
			ulist=[x for x in range(num_units) if node_groups[x]==klist[k]]
			mean=sum(nodes[x][idx+4] for x in ulist) / len(ulist)
			s+="\t"+str(mean)
		print s
	print "-----------------------end of Stat-------------------------"
	return r2

def search_stat():
	arcpy_print("----------------search statistics----------------------")
	arcpy_print("Unit move time: "+ str(time_op[0]) )
	arcpy_print("Renato update time: "+ str(time_Renato_update) )
	arcpy_print("Renato search time: "+str(time_Whitaker))
	arcpy_print("repair time: "+ str(time_repair) )
	arcpy_print("check edge unit time: "+str(time_check_edge_unit))
	arcpy_print("update_centers time: "+ str(time_update_centers) )
	arcpy_print("spp regions: "+ str(len(region_pool)) )
	arcpy_print("spp pooling time: "+ str(time_spp) )
	arcpy_print("connectivity check time: "+ str(time_check))
	arcpy_print("----------------end of search statistics----------------")
 
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
  

