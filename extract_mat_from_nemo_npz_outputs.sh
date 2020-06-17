#!/bin/bash

#for i in {1..23}; do python3 nemo_getnumerator.py --input NeMo2_outputs/SUB"$i"_lesion_1mmMNI_nemo_output_chaco_allref.npz --asum NeMo2_outputs/nemo_Asum_weighted_endpoints.npz --output SUB"$i"_output.mat ; done 
#for i in {1..23}; do python3 nemo_getvoxelfraction.py --input NeMo2_outputs/SUB"$i"_lesion_1mmMNI_nemo_output_chaco_allref.npz  --output SUB"$i"_output_fraction.mat ; done 
#for i in {1..23}; do python3 nemo_getvoxelfraction_endpoints.py --input NeMo2_outputs/SUB"$i"_lesion_1mmMNI_nemo_output_chaco_allref.npz --endpointmask NeMo2_outputs/nemo_endpoints_mask.npz --output SUB"$i"_output_fraction_endpoints.mat ; done 
for i in {1..23}; do python3 nemo_getnumerator_allpoints.py --asum ~/Documents/Thesis/colossus_syncs/disconnectivity_maps/nemo_Asum_weighted.npz --input ~/Documents/Thesis/colossus_syncs/nemo_processing/NeMo2_outputs/SUB"$i"_lesion_1mmMNI_nemo_output_chaco_allref.npz  --output ~/Documents/Thesis/colossus_syncs/nemo_processing/SUB"$i"_output_numerator_allpoints.mat ; done 
