import momi

#import os
#os.environ['KMPDUPLICATELIB_OK']='True'

#loading VCF

sfs=momi.Sfs.load("/mnt/students/anthony/momi_work/empirical_data/momi_dem1_filtered_indels.sfs.gz")

#model for inference
Ne=(1800+2100)/2  

freedman_model=momi.DemographicModel(
    N_e=Ne, gen_time=3, muts_per_gen= 1e-8)
    
##Adding leaves

freedman_model.add_leaf("Dingo", N=Ne, t=(805+957)/2)
freedman_model.add_leaf("BorneoVD", N=Ne)
freedman_model.add_leaf("NGSD", N=Ne)
freedman_model.add_leaf("VietnameVD", N=Ne)


#adding parameters
freedman_model.add_time_param("t_split_DNG_NGSD", lower=(805+957)/2)
freedman_model.add_time_param("t_split_Borneo_DNG", lower_constraints=["t_split_DNG_NGSD"])
freedman_model.add_time_param("t_split_Viet_Borneo",lower_constraints=["t_split_Borneo_DNG"] )

#Demographic events
freedman_model.move_lineages("Dingo","NGSD", t="t_split_DNG_NGSD")
freedman_model.move_lineages("BorneoVD", "NGSD", t="t_split_Borneo_DNG")
freedman_model.move_lineages("NGSD", "VietnameVD", t="t_split_Viet_Borneo")

#inference
freedman_model.set_data(sfs)
freedman_model.optimize(method="TNC")
freedman_model.get_params()