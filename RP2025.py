import homodistrictingv3 as d
import time,sys,random
#input parameters: filename1 filename2 number_of_facilities pop_size timelimit spp seed max_loops_solution_not_improved
attrs=0 #number of attributes (0 defualt number defined in the data file)
first=4 #number of the column as first attribute 
psize=10 #population size
maxloops=1000 #max_loops_solution_not_improved
timelimit=3600 #timelimit for searching
spp=0 #set partitioning
n=0 # number of regions
if len(sys.argv) >= 2:
    fn=sys.argv[1]
if len(sys.argv) >= 3:
    fn2=sys.argv[2]
if len(sys.argv) >= 4:
    attrs=int(sys.argv[3])
if len(sys.argv) >= 5:
    first=int(sys.argv[4])
if len(sys.argv) >= 6:
    n=int(sys.argv[5])
if len(sys.argv) >= 7:
    psize=int(sys.argv[6])
if len(sys.argv) >= 8:
    maxloops=int(sys.argv[7])

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

d.location_problem=3 # 3 for regionalization
d.spatial_contiguity=2
#0 non-connected, non-spatial clustering 
#1 strictly-connected regions
#2 connected with several parts 
d.spatial_contiguity_minimum_percentage=5.0
#minimun partial regions used when spatial_contiguity=2
#the area of a partial region must >= 5% of average regional area

t0=time.time()
d.read_pmp_climate_instance(fn,fn2,attrs,first) #PMP benchmark

d.is_spp_modeling=0
d.solution_similarity_limit=10 #max(10.0,100.0/n)
d.heuristic_time_limit=timelimit 
d.multi_start_count=psize  #population_size
d.max_loops_solution_not_improved=maxloops #used for search termination

d.ils_large_homoDP_fast(n)

print "=========================Final results========================="
print "objective:",d.biobjective
print "objective_R2",d.objective_R2
print "pool size",len(d.region_pool)
print "total time",time.time()-t0

save_homoDP_sol() 

