#this is a truncated model from the 2014 Freedman paper to estimate the split time of the dingo
#there is no migration in this model

from __future__ import print_function
import msprime
import numpy as np
import random

#setting up demes

N_DNG=(1800+2100)/2   #0
N_NGS=(1800+2100)/2   #1
N_BSJ=(2300+2900)/2   #2
N_CHW=(4700+6100)/2   #3


#converting times to number of generations

gen_time=3

#split_time=12000

T3=13400/gen_time #bottleneck of chinese wolf
T4=14900/gen_time



population_configurations = [
	msprime.PopulationConfiguration(
		initial_size=N_DNG),
	msprime.PopulationConfiguration(
		initial_size=N_NGS),
	msprime.PopulationConfiguration(
		initial_size=N_BSJ),
	msprime.PopulationConfiguration(
		initial_size=N_CHW),
	]

# dated 805-957 years before present (95.4% probability) with median age of 902 years.

ACAD_dingo_age=902
ancient_samples = [
	msprime.Sample(population=0,time=ACAD_dingo_age/gen_time),
	msprime.Sample(population=0,time=ACAD_dingo_age/gen_time),

	msprime.Sample(population=1,time=0),
	msprime.Sample(population=1,time=0),

	msprime.Sample(population=2,time=0),
	msprime.Sample(population=2,time=0),

	msprime.Sample(population=3,time=0),
	msprime.Sample(population=3,time=0),
]
"""
demographic_events=[
 	#joining DIN to NGS
 	msprime.MassMigration(
 		time=T1,source=0,destination=1,proportion=1),
 	#joining DIN with BSJ
 	msprime.MassMigration(
 		time=T2,source=1,destination=2,proportion=1),
 	#bottleneck for CHW
 	msprime.PopulationParametersChange(
 		time=T3,initial_size=13000,population_id=3),
 	#joining CHW with BSJ
 	msprime.MassMigration(
 		time=T4,source=2,destination=3,proportion=1)
]
"""

def split_times(split_time):
	T1=split_time/gen_time
	demographic_events=[
	#joining DIN to NGS
 	msprime.MassMigration(
 		time=T1,source=0,destination=1,proportion=1),
 	#joining DIN with BSJ
 	msprime.MassMigration(
 		time=T2,source=1,destination=2,proportion=1),
 	#bottleneck for CHW
 	msprime.PopulationParametersChange(
 		time=T3,initial_size=13000,population_id=3),
 	#joining CHW with BSJ
 	msprime.MassMigration(
 		time=T4,source=2,destination=3,proportion=1)
	]
	return demographic_events

def double_split_time(t1,t2):
	T1=t1/gen_time
	T2=t2/gen_time
#	print(t1,t2)
	demographic_events=[
	#joining DIN to NGS
 	msprime.MassMigration(
 		time=T1,source=0,destination=1,proportion=1),
 	#joining NGS with BSJ
 	msprime.MassMigration(
 		time=T2,source=1,destination=2,proportion=1),
 	#bottleneck for CHW, we can remove this for now as we are not looking at the CHW genome
 	#msprime.PopulationParametersChange(
 	#	time=T3,initial_size=13000,population_id=3),
 	#joining CHW with BSJ
 	msprime.MassMigration(
 		time=T4,source=2,destination=3,proportion=1)
	]
	return demographic_events

#t1=random.uniform(1,14900)
#t2=random.uniform(t1,14900)
#double_split_time(t1,t2)
#simulation code

def dingo_sim():
	for split_time in range (0,12799,50):
	#	print(" ")
	#	print("Split time is ", split_time)
	#	print(" ")
		fab_list = []
		end=101
		for seed in range(1,end):
			tree_sequence=msprime.simulate(Ne=N_DNG,length=1e8,recombination_rate=5e-7,population_configurations=population_configurations,demographic_events=split_times(split_time),
            	mutation_rate=1e-8, random_seed=seed)
			haplotypes=tree_sequence.haplotypes()
			haplotype_list=[]
			for i in haplotypes:
				haplotype_list.append(i)
			haplotype_list.append(tree_sequence.get_num_mutations())
			sim=haplotype_list
			print_list(sim)
			het_list=het_B_finder(sim)
	#		print("het list is",het_list)
			derived_list=derived_finder(sim,het_list)
			print(split_time,seed,fab_calc(sim,derived_list),len(derived_list[1]))
	#		print("Simulation ", seed, " fab stat is", fab_calc(sim,derived_list))
	#		print(type(fab_avg))
	#		print(type(fab_calc(sim,derived_list)))
	#		fab = fab_calc(sim,derived_list)
			#if fab is not None:
				#fab_list.append(fab)
			#if(seed==end-1):
			#	print("Mean fab is ", np.mean(fab_list))
			#	print("Variance of fab is ", np.var(fab_list))
			#	fab_list=[]

def abc():
	#for seed in range (1,10000):
	count=0
	for _ in range (0,1000000):
		seed=random.randint(1,(2**32)-1)
		#t2=random.uniform(11700,13700)
		t2=12800
		t1=random.uniform(903,t2)
	#	print(seed,t1,t2)
		tree_sequence=msprime.simulate(Ne=N_DNG,length=1e7,recombination_rate=0.33e-8,samples=ancient_samples,population_configurations=population_configurations,demographic_events=double_split_time(t1,t2),
            	mutation_rate=0.66e-8, random_seed=seed)
		haplotypes=tree_sequence.haplotypes()
		haplotype_list=[]
		for i in haplotypes:
			haplotype_list.append(i)
		haplotype_list.append(tree_sequence.get_num_mutations())
		sim=haplotype_list
	#	print_list(sim)
		het_list=het_B_finder(sim)
		derived_list=derived_finder(sim,het_list)
	#	print("het list is",het_list)
		fab_k=fab_calc(sim,derived_list)
	#	print(fab_k,d_calc(fab_k,fab_obs),t1,t2,len(derived_list[1]),seed)
		count=count+1
		print(fab_k,t1,t2,len(derived_list[1]),seed,count)


def test_abc():
	t1=random.uniform(1,14900)
	t2=random.uniform(t1,14900)
	print(t1,t2)
	tree_sequence=msprime.simulate(Ne=N_DNG,length=1000,recombination_rate=5e-7,population_configurations=population_configurations,demographic_events=split_times(t1),
            	mutation_rate=1e-8, random_seed=1,samples=ancient_samples)
	haplotypes=tree_sequence.haplotypes()
	haplotype_list=[]
	for i in haplotypes:
		haplotype_list.append(i)
	haplotype_list.append(tree_sequence.get_num_mutations())
	sim=haplotype_list
	print_list(sim)
	het_list=het_B_finder(sim)
	derived_list=derived_finder(sim,het_list)
	print("het list is",het_list)
	fab_k=fab_calc(sim,derived_list)
	print(fab_k)

#test_abc()





#fab functions go here

#sim=simulate()

def print_list(sim):
	for st in sim:
		print (st)
		print (" ")

#print_list()

def het_B_finder(sim):
	het_B_list=[]
	for base in range(0,int(sim[8])):
		if sim[2][base]!=sim[3][base]:
			het_B_list.append(base)
	return het_B_list
	
#het_list=het_B_finder()
#print("het list is", het_list)

#allele in B that is not in outgroup is the derived allele

def derived_finder(sim,het_list):
	ancestral_allele_list=[]
	ancestral_allele_position_list=[]
	derived_allele_list=[]
	for base in het_list:
		#check if outgroup is homozygous first
		if sim[4][base]==sim[5][base]:
			ancestral_allele_list.append(sim[4][base])
			ancestral_allele_position_list.append(base)
	for base in range(0,len(ancestral_allele_list)):
		if ancestral_allele_list[base]=="0":
			derived_allele_list.append("1")
		elif ancestral_allele_list[base]=="1":
			derived_allele_list.append("0")
		else:
			print ("ERROR IN DERIVED_FINDER FUNCTION")
	derived_allele_list=[derived_allele_list,ancestral_allele_position_list]
	return derived_allele_list

#derived_list=derived_finder()

#print ("derived list is",derived_list)

"""
def fab_calc():
	count=0


	for base in range(0,len(derived_list[1])):
		if sim[0][derived_list[1][base]]==derived_list[0][base]:
			count=count+1
		if sim[1][derived_list[1][base]]==derived_list[0][base]:
			count=count+1
	fab=count/len(derived_list[1])*0.5


	##testing..
allele_1=[]
allele_2=[]
	for base in range(0,len(derived_list[1])):
		allele_1.append(sim[0][derived_list[1][base]])
		allele_2.append(sim[1][derived_list[1][base]])
	allele=[allele_1,allele_2]
	for st in allele:
		print (st)
		print (" ")
	print(derived_list[0])

	for base in range(0,len(derived_list[1])):
		if allele_1[0][base]==derived_list[0][base]:
			count=count+1
		if allele_2[0][base]==derived_list[0][base]:
			count=count+1
	fab=count/len(derived_list[1])*0.5
	return fab
"""

def fab_calc(sim,derived_list):
	if len(derived_list[1])==0:
		return None
	count=0
	allele_1=[]
	allele_2=[]
	for base in range(0,len(derived_list[1])):
		allele_1.append(sim[0][derived_list[1][base]])
		allele_2.append(sim[1][derived_list[1][base]])
	allele=[allele_1,allele_2]

#	for st in allele:
#		print (st)
#		print (" ")
#	print(derived_list[0])
#	print(len(derived_list[1]))
#	print(len(allele[0]))
#	print(len(allele[1]))


	for base in range(0,len(derived_list[1])):
		if allele[0][base]==derived_list[0][base]:
			count=count+1
	#		print("count is ", count)
		if allele[1][base]==derived_list[0][base]:
			count=count+1
	#		print("count is ", count)
#	print("Final count is ", count)
	count=count*1.0
#	print("Did it work??", count)
	fab=count/len(derived_list[1])*0.5
	return fab

#print("fab stat is", fab_calc())



#simulate function returns a string list of haplotypes
"""
def simulate():
	tree_sequence=msprime.simulate(Ne=N_DNG,length=1e8,recombination_rate=1e-9,population_configurations=population_configurations,demographic_events=demographic_events,
            mutation_rate=1e-8, random_seed=50)
	haplotypes=tree_sequence.haplotypes()
	haplotype_list=[]
	for i in haplotypes:
		haplotype_list.append(i)
	haplotype_list.append(tree_sequence.get_num_mutations())
	return haplotype_list
"""

def d_calc(fab_k,fab_obs):
	if fab_k==None:
		return 9000
	else:
		return (fab_obs-fab_k)**2




##run things here

#dingo_sim()
abc()


