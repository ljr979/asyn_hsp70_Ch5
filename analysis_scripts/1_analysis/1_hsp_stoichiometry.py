from src.py4bleaching.py4bleaching import analysis

#This script runs the py4bleaching pipeline on all the trajectories
input_folder = 'python_results/Trajectories/'
output_folder = 'python_results/py4bleaching/Trajectories/coloc/'

#change this according to the model that you'd like to use (from the repo with all the models)
model_name = 'Model_2'



analysis.pipeline(input_folder, output_folder, probability_threshold=0.7, model_name=model_name, x_norm=True)
