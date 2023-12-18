
from enum import unique
import os, re
import threading
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
import math
results_output = f'results/5_molsize_end_mid/'
if not os.path.exists(results_output):
    os.makedirs(results_output)

output_folder=results_output




def plot_violin_count(molecule_counts, output_folder, dicto, order_of_experiment):
#plot distribution of fibril lengths cooc and not coloc
    #for treatment, df in molecule_counts.groupby('treatment'): 
        #fig, ax = plt.subplots()
    
    molecule_counts['treatment_to_plot']=molecule_counts['Treatment'].map(dicto)
    order_of_experiment= order_of_experiment


    ax = sns.violinplot(
        x='treatment_to_plot',
        y= 'last_step_mol_count', 
        data=molecule_counts, 
        scale='width', 
        hue='end_or_middle', 
        palette='Greens', 
        order=order_of_experiment,
        alpha=0.45)
    ax.legend(loc='upper center', ncol=2)
    ax.set_ylim(0,25)
    # ax.annotate(f"median colocalised length = {median_length_colocal_fibrils}nm", xy=(0.1, 0.8), xycoords='figure fraction')
    plt.ylabel(f'# of molecules/foci')
    plt.xlabel('Experimental condition')
    #plt.show()
    
    plt.title(f'# of molecules/foci @ end or middle ')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f'{output_folder}/end_middle_molecules_per_foci.png')
    plt.savefig(f'{output_folder}/end_middle_molecules_per_foci.svg')
    plt.show()


def calculate_percent_summary(molecule_counts, output_folder):
    coloc_dict=[]
    #
    for treatment, df in molecule_counts.groupby('variable1'):
        LIST=[]
        total_hsp_fibs=len(df)
        total_ends=len(df[df['end_or_middle']=='END'])
        total_mids=len(df[df['end_or_middle']=='MIDDLE'])

        percent_end=total_ends/total_hsp_fibs*100

        percent_mid=total_mids/total_hsp_fibs*100

        LIST=treatment, percent_end, percent_mid
        coloc_dict.append(LIST)

    
    coloc_dict=pd.DataFrame(coloc_dict, columns=['treatment', 'percentage_end', 'percentage_middle'])
    coloc_dict.to_csv(f'{current_output}/Summary_percentage_end_middle.csv')
    return coloc_dict



paths = {
    #A8 no ATP
    ('Experiment 75-1', 'A8-noATP-1', 'old'): 'D:/Experiments/Experiment_75/python_results/1/',

    ('Experiment 75-2', 'A8-noATP-2', 'old'): 'D:/Experiments/Experiment_75/python_results/2/',


    #A8 + ATP
    ('Experiment 72-1', 'A8-ATP-1', 'old'): 'D:/Experiments/Experiment_72/python_results/1/',
    #('Experiment 72-2', 'A8-ATP-2'): 'D:/Experiments/Experiment_72/python_results/2/',
    ('Experiment 84-1', 'A8-ATP-3', 'old'): 'D:/Experiments/Experiment_84/python_results/FC1/',


    #JB1+A8+ATP

    #('Experiment 66-2','JB1-A8-1'):'D:',
    # ('Experiment 73-2','JB1-A8-2'):'D:/Experiments/Experiment_73/python_results/2/',
    # ('Experiment 77-2','JB1-A8-3'):'D:/Experiments/Experiment_77/python_results/2/',
    # ('Experiment 77-1','JB1-A8-4'):'D:/Experiments/Experiment_77/python_results/1/',
    # ('Experiment 76-1','JB1-A8-5'):'D:/Experiments/Experiment_76/python_results/1/',
    ('Experiment 84-2', 'JB1-A8-5', 'new'): 'D: / Experiments/Experiment_84/python_results/',

    ('Experiment 86-2', 'JB1-A8-1', 'new'): 'D:/Experiments/Experiment86-FC2-bleaches/python_results/',
    # ('Experiment 87-2','JB1-A8-7'):'D:/Experiments/Experiment87-FC2/python_results/',
    ('Experiment 90-1', 'JB1-A8-2', 'new'): 'D:/Experiments/Experiment_90-FC1/python_results/',
    ('Experiment 90-2', 'JB1-A8-9', 'new'): 'D:/Experiments/Experiment_90-FC2/python_results/',
    # ('Experiment 85-3','JB1-A8-0'):'D:/Experiments/Experiment85-FC3/python_results/FC3/coords_added/',
    # ('Experiment 92-1', 'JB1-A8-A'): 'D:/Experiments/Experiment_92-fibrils_SOD1_control/python_results/',
    # ('Experiment 96-1', 'JB1-A8-B'): 'D:/Experiments/Experiment96-SOD1-control/python_results/',
    #('Experiment 96-2', 'JB1-A8-3'): 'D:/Experiments/Experiment96-FC2-Hsp110-titration/python_results/FC2-1/subunits_per_foci/avg-t0-size-nobleach/',
    ('Experiment 98-1', 'JB1-A8-4', 'new'): 'D:/Experiments/Experiment98-SOD1/python_results/',


    #MIX DARK A8 2UM & LABELLED A8
    ('Experiment 86-1', 'darkA8-647A8-1', 'new'): 'D:/Experiments/Experiment86-FC1/python_results/',

    ('Experiment 87-1', 'darkA8-647A8-2', 'new'): 'D:/Experiments/Experiment87-FC1/python_results/',


    #JB1+A8+1nM hsp110-excessA8
    # ('Experiment 81-1','JB1-A8-110-1nM-1'):'D:/Experiments/Experiment_81/python_results/1/',
    ('Experiment 81-2', 'JB1-A8-110-1nM-1', 'old'): 'D:/Experiments/Experiment_81/python_results/2/',
    # ('Experiment 84-3','JB1-A8-110-1nM-3'):'D:/Experiments/Experiment_84/python_results/FC3/',
    ('Experiment 90-2', 'JB1-A8-110-1nM-2', 'new'): 'D:/Experiments/Experiment_90-FC2/python_results/',
    ('Experiment 100-2', 'JB1-A8-110-5nM-2', 'new'): 'D:/Experiments/Experiment-100_b/python_results/',
    ('Experiment 96-2', 'JB1-A8-110-5nM-3', 'new'): 'D:/Experiments/Experiment96-FC2-Hsp110-titration/python_results/',

    #JB1+A8+0.5uM hsp110- excess A8
    ('Experiment 81-1', 'JB1-A8-110-1nM-1', 'old'): 'D:/Experiments/Experiment_81/python_results/1/',

    ('Experiment 84-2', 'JB1-A8-110-0.5uM-1', 'new'): 'D:/Experiments/Experiment_84/python_results/',

    ('Experiment 86-2', 'JB1-A8-110-0.5uM-2', 'new'): 'D:/Experiments/Experiment86-FC2-bleaches/python_results/',
    # ('Experiment 90-2','JB1-A8-110-0.5uM-3'):'D:/Experiments/Experiment_90-FC2/python_results/',
    ('Experiment 100-1', 'JB1-A8-110-0.5uM-3', 'new'): 'D:/Experiments/Experiment100_hsp110-SOD-compare/python_results/',
    # ('Experiment 87-2', 'JB1-A8-110-0.5uM-4'): 'D:/Experiments/Experiment87-FC2/python_results/',
    ('Experiment 90-1', 'JB1-A8-110-0.5uM-5', 'new'): 'D:/Experiments/Experiment_90-FC1/python_results/',


    #JB1+A8+2uM hsp110- excess A8
    ('Experiment 90-2', 'JB1-A8-110-2uM-1', 'new'): 'D:/Experiments/Experiment_90-FC2/python_results/',
    ('Experiment 96-2', 'JB1-A8-110-2uM-2', 'new'): 'D:/Experiments/Experiment96-FC2-Hsp110-titration/python_results/',
    ('Experiment 90-1', 'JB1-A8-110-2uM-3', 'new'): 'D:/Experiments/Experiment_90-FC1/python_results/',

    #JB1+A8+0.5uM hsp110- NO excess A8
    # ('Experiment 87-2','JB1-110-0.5uM-1'):'D:/Experiments/Experiment87-FC2/python_results/',
    # ('Experiment 90-1','JB1-110-0.5uM-2'):'D:/Experiments/Experiment_90-FC1/python_results/',


    #JB1+A8+2uM hsp110- NO excess A8
    # ('Experiment 90-1','JB1-110-2uM-1'):'D:/Experiments/Experiment_90-FC1/python_results/',



    #jb1_a8+SOD1@ 0.5uM
    ('Experiment 92-1', 'JB1-A8-SOD-0.5uM-1', 'new'): 'D:/Experiments/Experiment_92-fibrils_SOD1_control/python_results/',
    # ('Experiment 96-1', 'JB1-A8-SOD-0.5uM-2'): 'D:/Experiments/Experiment96-SOD1-control/python_results/',
    ('Experiment 98-1', 'JB1-A8-SOD-0.5uM-3', 'new'): 'D:/Experiments/Experiment98-SOD1/python_results/',
    ('Experiment 100-1', 'JB1-A8-SOD-0.5uM-4', 'new'): 'D:/Experiments/Experiment100_hsp110-SOD-compare/python_results/',


}


files=[]
for key, value in paths.items():
    coloc_files =[[f'{root}/{filename}' for filename in files if 'end_vs_middle_collated.csv' in filename] for root, dirs, files in os.walk(f'{value}')]
    
    coloc_files=[item for sublist in coloc_files for item in sublist]
    
    df=pd.DataFrame(coloc_files, columns=['path'])
    df[['Experiment_number', 'Treatment-replicate', 'analysis-version']]=key
    files.append(df)

collated_collated_colocal=pd.concat(files)
old=collated_collated_colocal[collated_collated_colocal['analysis-version']=='old']
new=collated_collated_colocal[collated_collated_colocal['analysis-version']=='new']
#from collate colocalisation




old_collated=[]
paths_read=[]
for  version, df in collated_collated_colocal.groupby('analysis-version'):
     df
     version

     if version=='old':
        for (path,Treatment), row in df.groupby(['path', 'Treatment-replicate']):
            path
            row
            
            Treatment
            path=path.replace('\\','/')
            replicate=Treatment.split('-')[-1]

            treatment_from_file=path.split('/')[-1].split('_')[-4]
            treatment=Treatment.rstrip(Treatment[-1])
            Exp_num=list(row['Experiment_number'])[0]
            if path not in paths_read:
                print(Exp_num)
                coloc_old=pd.read_csv(path)
                paths_read.append(path)
                #fresh_df=pd.DataFrame(coloc['percent_colocalised'])
                coloc_old['Treatment_from_file']=treatment_from_file
                coloc_old['Experiment_number']=Exp_num
                coloc_old['path']=path
                coloc_old['replicate']=replicate
                coloc_old['Treatment']=treatment
                old_collated.append(coloc_old)
            #protein_coloc=list(coloc['proteins_colocalised'])[0]
old_collated=pd.concat(old_collated)
drops=['Unnamed: 0','Unnamed: 0.1','Class', 'Frame','EndX1, EndX2, EndY1, EndY2', 'Contour_ID']
old_collated.drop([col for col in old_collated.columns.tolist() if col in drops], axis=1, inplace=True)
old_collated=old_collated[old_collated['Experiment_number']!='Experiment 72-1']
new_collated=[]

for  version, df in collated_collated_colocal.groupby('analysis-version'):
     df
     version
     if version=='new':

        for (path,Treatment), row in df.groupby(['path', 'Treatment-replicate']):
            path
            


            path=path.replace('\\','/')
            
            replicate_t=Treatment
            replicate=replicate_t.split('-')[-1]

            treatment_from_file=path.split('/')[-1].split('_')[-4]
            treatment=Treatment.rstrip(replicate_t[-1])
            Exp_num=list(row['Experiment_number'])[0]

            row


            coloc_new=pd.read_csv(path)
            paths_read.append(path)
            fresh_df=coloc_new
            
            fresh_df['Treatment_from_file']=treatment_from_file
            fresh_df['Experiment_number']=Exp_num
            fresh_df['path']=path
            fresh_df['replicate']=replicate
            fresh_df['Treatment']=treatment
            
            new_collated.append(fresh_df)

new_collated=pd.concat(new_collated)




new_manual_filtered=[]
old_manual_filtered=[]


t = old_collated[old_collated['Experiment_number']
                      == 'Experiment 75-1']
t_filter = t
old_manual_filtered.append(t_filter)

t = old_collated[old_collated['Experiment_number']
                      == 'Experiment 75-2']
old_manual_filtered.append(t)
#atp------------------------------------------------------
t = old_collated[old_collated['Experiment_number']
                      == 'Experiment 84-1']
t_filter=t[t['variable2']=='a8-atp-flowout']
old_manual_filtered.append(t_filter)


#----------------------------------------------------------
#JB1+A8
t = new_collated[new_collated['Experiment_number']
                      == 'Experiment 98-1']
#FILTER FOR JUST THE A8 AND JB1 HERE ANDDDDD colocalised v non-colocalised ()
t_filter = t[t['Treatment'] == 'JB1-A8-']
t_filter = t_filter[t_filter['variable1'] == 'a8-b1-flowout']

new_manual_filtered.append(t_filter)

#--------


t = new_collated[new_collated['Experiment_number'] == 'Experiment 90-2']
t_filter = t[t['Treatment'] == 'JB1-A8-']
t_filter = t_filter[t_filter['variable1'] == 'a8-b1']
t_filter['replicate'] = 4

new_manual_filtered.append(t_filter)





#------------------------------------
#now onto the next treatment!!!!!!!!!!!!! JB1 + A8 + 110  @ low conc (sub A8)


t = new_collated[new_collated['Experiment_number']
                      == 'Experiment 90-2']
t_filter = t[t['Treatment'] == 'JB1-A8-110-1nM-']
t_filter = t_filter[t_filter['variable1'] == 'a8-10nm-atp-1nmhsp110']
t_filter['Treatment'] = 'JB1-A8-110-5nM-'

new_manual_filtered.append(t_filter)


t = old_collated[old_collated['Experiment_number'] == 'Experiment 81-1']
old_manual_filtered.append(t)


# t = new_collated[new_collated['Experiment_number']
#                       == 'Experiment 100-2']

# t_filter = t[t['Treatment'] == 'JB1-A8-110-5nM-']
# manual_filtered.append(t_filter)
#------------------------------------
#now onto the next treatment!!!!!!!!!!!!! JB1 + A8 + 110  @ middle conc (above A8)

t = new_collated[new_collated['Experiment_number']
                      == 'Experiment 100-1']
t_filter = t[t['Treatment'] == 'JB1-A8-110-0.5uM-']
t_filter = t_filter[t_filter['variable1'] == 'a8-b1-110']

new_manual_filtered.append(t_filter)


t = new_collated[new_collated['Experiment_number']
                      == 'Experiment 90-1']
t_filter = t[t['Treatment'] == 'JB1-A8-110-0.5uM-']
t_filter = t_filter[t_filter['variable1'] == 'a8-b1-0.5umhsp110-1h']

new_manual_filtered.append(t_filter)


t = new_collated[new_collated['Experiment_number']
                      == 'Experiment 86-2']
t_filter = t[t['Treatment'] == 'JB1-A8-110-0.5uM-']
t_filter = t_filter[t_filter['variable1'] == 'a8-jb1-110-0.5um']

new_manual_filtered.append(t_filter)

# t = new_collated[new_collated['Experiment_number']
#                       == 'Experiment 84-2']
# t_filter = t[t['Treatment_from_file'] ==
#              'Experiment84_A8-JB1-110-0.5uM_flowout']
# t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-110-0.5uM-']
# new_manual_filtered.append(t_filter)

#-------------------------------------------
#now onto the next treatment!!!!!!!!!!!!! JB1 + A8 + 110  @ high conc (well above A8)

t = new_collated[new_collated['Experiment_number']
                      == 'Experiment 90-1']
t_filter = t[t['variable1'] == 'a8-b1-2umhsp110-1h']
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-110-2uM-']
new_manual_filtered.append(t_filter)


t = new_collated[new_collated['Experiment_number']
                      == 'Experiment 90-2']
t_filter = t[t['variable1'] == 'a8-10nm-atp-2umhsp110']
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-110-2uM-']
new_manual_filtered.append(t_filter)


#__________________________________________________________
#finally, onto SOD1
t = new_collated[new_collated['Experiment_number']
                      == 'Experiment 100-1']
t_filter = t[t['variable1'] == 'a8-b1-sod1']
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-SOD-0.5uM-']

new_manual_filtered.append(t_filter)


t = new_collated[new_collated['Experiment_number']
                      == 'Experiment 92-1']
t_filter = t[t['variable1'] == 'a8-b1-+sod1-2h']
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-SOD-0.5uM-']

new_manual_filtered.append(t_filter)


t = new_collated[new_collated['Experiment_number']
                      == 'Experiment 98-1']
t_filter = t[t['variable1'] == 'a8-b1-sod1']
t_filter = t_filter[t_filter['Treatment'] == 'JB1-A8-SOD-0.5uM-']


new_manual_filtered.append(t_filter)



#now we want to concat everything in these lists
old_manual_filtered=pd.concat(old_manual_filtered)
new_manual_filtered=pd.concat(new_manual_filtered)

collated=pd.concat([old_manual_filtered, new_manual_filtered])
drop=['Unnamed: 0', 'Contour_ID', 'protein1']
collated.drop([col for col in collated.columns.tolist() if col in drop], axis=1, inplace=True)


collated.to_csv(f'{results_output}molsize_end_mid_collated.csv')
#_______________---------------------------------------------------------


def plot_boxes_count(molecule_counts, output_folder, dicto, order_of_experiment, save=False):
#plot distribution of fibril lengths cooc and not coloc
    #for treatment, df in molecule_counts.groupby('treatment'): 
        #fig, ax = plt.subplots()
    
    molecule_counts['treatment_to_plot']=molecule_counts['Treatment'].map(dicto)
    order_of_experiment= order_of_experiment


    ax = sns.boxplot(
        x='treatment_to_plot',
        y= 'rel_size', 
        data=molecule_counts, 
        hue='end_or_middle', 
        palette='Greens', 
        order=order_of_experiment,
        )
    ax.legend(loc='upper center', ncol=2)
    ax.set_ylim(-0.05,0.8)
    # ax.annotate(f"median colocalised length = {median_length_colocal_fibrils}nm", xy=(0.1, 0.8), xycoords='figure fraction')
    plt.ylabel(f'# of molecules/foci')
    plt.xlabel('Experimental condition')
    #plt.show()
    
    plt.title(f'# of molecules/foci @ end or middle ')
    plt.xticks(rotation=45)
    plt.tight_layout()
    if save==True:
        plt.savefig(f'{output_folder}/box_end_middle_molecules_per_foci.png')
        plt.savefig(f'{output_folder}/box_end_middle_molecules_per_foci.svg')
    plt.show()




dicto={
    'A8-noATP-':'-ATP',
  'A8-ATP-':'+ATP',
    'JB1-A8-110-1nM-': '+110 1nM', 
    'JB1-A8-':'JB1+A8',
       'JB1-A8-110-5nM-':'+110 5nM',
         'JB1-A8-110-0.5uM-':'+110 0.5uM', 
         'JB1-A8-110-2uM-':'+110 2uM',
       'JB1-A8-SOD-0.5uM-':'+SOD 0.5uM'
       }
       


order_of_experiment= ['-ATP',
 '+ATP',
 'JB1+A8',
 '+110 1nM',
 '+110 5nM',
 '+110 0.5uM',
 '+110 2uM',
 '+SOD 0.5uM']

molecule_counts=pd.read_csv(f'{results_output}molsize_end_mid_collated.csv')
output_folder=results_output

plot_violin_count(molecule_counts=molecule_counts, output_folder=output_folder, dicto=dicto, order_of_experiment=order_of_experiment)
plot_boxes_count(molecule_counts, output_folder, dicto, order_of_experiment)
coloc_dict=calculate_percent_summary(molecule_counts, output_folder)







#-----------------------------------------
#now I just want to test whether when I normalise the intensity as I had previously done for molsize experiment, if the terend is the same.

molecule_counts.last_step_mol_count=molecule_counts.last_step_mol_count.mask(molecule_counts.last_step_mol_count.lt(0),1)
normed_store=[]
for experiment, w_df in molecule_counts.groupby('Experiment_number'):
    experiment
    w_df
    maxo=w_df['last_step_mol_count'].max()
    w_df['rel_size']=w_df['last_step_mol_count']/maxo
    normed_store.append(w_df)
nomalised_intensity=pd.concat(normed_store)


nomalised_intensity.to_csv(f'{results_output}normalised_molsize_endmid.csv')



plot_boxes_count(molecule_counts= nomalised_intensity, output_folder=output_folder, dicto=dicto, order_of_experiment=order_of_experiment, save=True)

per_exp_mean=[]
for experiment, df in nomalised_intensity.groupby('Experiment_number'):
    #calculate the ratio between middle and end (subunit/brightness relative to each other)
    for treatment, df1 in df.groupby('Treatment'):
        mean_end=np.mean(df1[df1['end_or_middle']=='END']['rel_size'])
        mean_mid=np.mean(df1[df1['end_or_middle']=='MIDDLE']['rel_size'])
        ratio=mean_mid/mean_end
        per_exp_mean.append([experiment, treatment, mean_end, mean_mid, ratio])

    plot_boxes_count(molecule_counts= df, output_folder=output_folder, dicto=dicto, order_of_experiment=order_of_experiment, save=False)
per_exp_mean=pd.DataFrame(per_exp_mean, columns=['experiment', 'treatment', 'mean_end', 'mean_mid', 'ratio_mid_end'])
#okay that went kind of well but now I want to test whether I can plot the ratio of mid to end for each treatment so it's less complicated.
#descriptivestats saved
descrip=[]
for treatment, df in nomalised_intensity.groupby('Treatment'):
    treatment
    df
    med_end=np.median(df[df['end_or_middle']=='END']['rel_size'])
    med_mid=np.median(df[df['end_or_middle']=='MIDDLE']['rel_size'])
    range_end=np.ptp(df[df['end_or_middle']=='END']['rel_size'])
    range_mid=np.ptp(df[df['end_or_middle']=='MIDDLE']['rel_size'])
    ratio_mid_end_meds=med_mid/med_end
    ratio_mid_end_range=range_mid/range_end
    data=[treatment, med_end, med_mid, range_end, range_mid, ratio_mid_end_meds, ratio_mid_end_range]
    descrip.append(data)



descrip=pd.DataFrame(descrip, columns=['Treatment', 'Median_subunits_end', 'median_subunits_middle','range_subunits_end', 'range_subunits_mid', 'ratio_medians_mid_end', 'ratio_ranges_mid_ends'])
#not sure this tells me anything useful above actually tells me anything? (plot below for visualisation). I think the distribution is more useful?
fig,ax=plt.subplots()
ax=sns.barplot(data=descrip, x='Treatment', y='ratio_medians_mid_end', palette='Blues')
plt.ylim(0.7, 1.1)
plt.xticks(rotation=90)

fig,ax=plt.subplots()
ax=sns.barplot(data=descrip, x='Treatment', y='ratio_ranges_mid_ends', palette='Blues')
plt.ylim(0.8, 1.5)
plt.xticks(rotation=90)

#NOT USEFUL AT ALL. JUST DO KRUSKAL WALLIS ON THE DISTRIBUTION OF THE ORIGINAL DATA! PLOT THIS AS A BARPLOT WITH STANDARD ERROR OR AS A SWARMPLOT MAYBE?

order_of_experiment= ['A8-noATP-',
  'A8-ATP-',
 'JB1-A8-',
 'JB1-A8-110-1nM-',
 'JB1-A8-110-5nM-',
 'JB1-A8-110-0.5uM-',
 'JB1-A8-110-2uM-',
 'JB1-A8-SOD-0.5uM-']
fig, ax= plt.subplots (figsize=(8, 8))
g=sns.barplot(data=nomalised_intensity, x="Treatment", y="rel_size", hue="end_or_middle", 
order=order_of_experiment, 
palette='Blues', alpha=.6, edgecolor="grey")
#g.despine(left=True)
#sns.stripplot(
# x="Treatment", 
# y="Chaperones_per_length", 
# hue="end_or_middle", 
# order=order_of_experiment, 
# data=data, dodge=True, alpha=0.6, ax=g
#)
#g.set_axis_labels(f"Concentration of {protein} titrated onto fibrils", f"# of {protein} bound per unit fibril length (um)")
plt.ylim(0,0.3)
plt.xticks(rotation=90)
plt.ylabel('relative number of subunits/region')
plt.savefig(f'{results_output}number_subunits_per_regio_barplot.svg')
plt.savefig(f'{results_output}number_subunits_per_regio_barplot.png')
#plt.legend(title=f'end or middle', labels=['End', 'Middle'], loc='upper right')



#'Maybe now I'd like to do some stats on the a) distribution of molecules at end v middle for each treatment (kruskal wallis and dunns, on the green boxplots), but also b)the means (anova, blue graph)'


#a) we are going to do kruskal wallis on the distributions of normalised (relative) molecules SIZES i.e. relative number of subunits.
#split this up and turn each timepont into an array, so that we can compare between them easily
from scipy import stats
import scikit_posthocs as sp

import statsmodels.stats.multicomp as mc
stats_output='results/5_molsize_end_mid/'
if not os.path.exists(stats_output):
    os.makedirs(stats_output)
nomalised_intensity=pd.read_csv(f'{stats_output}normalised_molsize_endmid.csv')

def kruskal_treatments(df):
    a8noatp=df[df['Treatment']=='A8-noATP-']
    a8noatp=np.array(a8noatp['rel_size'])

    A8ATP=df[df['Treatment']=='A8-ATP-']
    A8ATP=np.array(A8ATP['rel_size'])

    b1a8=df[df['Treatment']=='JB1-A8-']
    b1a8=np.array(b1a8['rel_size'])

    B1A8110low=df[df['Treatment']=='JB1-A8-110-5nM-']
    B1A8110low=np.array(B1A8110low['rel_size'])

    B1A8110mid=df[df['Treatment']=='JB1-A8-110-0.5uM-']
    B1A8110mid=np.array(B1A8110mid['rel_size'])

    B1A8110high=df[df['Treatment']=='JB1-A8-110-2uM-']
    B1A8110high=np.array(B1A8110high['rel_size'])

    B1A8SOD=df[df['Treatment']=='JB1-A8-SOD-0.5uM-']
    B1A8SOD=np.array(B1A8SOD['rel_size'])


    #this part tells us whether there is a difference in some of the means of these data sets.
    result = stats.kruskal(a8noatp,A8ATP, b1a8,
    B1A8110low,
    B1A8110mid,
    B1A8110high,
    B1A8SOD)
    print(result)
    return result
def dunns_treatments(df, comparisons_list):
    data = [df[df['Treatment'] == 'A8-noATP-']['rel_size'], 
        df[df['Treatment'] == 'A8-ATP-']['rel_size'], 
        df[df['Treatment'] == 'JB1-A8-']['rel_size'],
        df[df['Treatment'] == 'JB1-A8-110-5nM-']['rel_size'],
        df[df['Treatment'] == 'JB1-A8-110-0.5uM-']['rel_size'],
        df[df['Treatment'] == 'JB1-A8-110-2uM-']['rel_size'],
        df[df['Treatment'] == 'JB1-A8-SOD-0.5uM-']['rel_size']]
    p_values = sp.posthoc_dunn(data, p_adjust='holm')


        
    signif=p_values < 0.05
    signif.columns=comparisons_list

    names_rows={
        1:'a8noatp',
        2:'A8ATP',
        3: 'b1a8',
    4:'B1A8110low',
    5:'B1A8110mid',
    6:'B1A8110high',
    7:'B1A8SOD'
    }
    signif.rename(index=names_rows, inplace=True)
    return signif


comparisons_list=['a8noatp','A8ATP', 'b1a8',
'B1A8110low',
'B1A8110mid',
'B1A8110high',
'B1A8SOD']

mids=nomalised_intensity[nomalised_intensity['end_or_middle']=='MIDDLE']
#this test will tell whether BETWEEN TREATMENTS, the distribution of molecule sizes at the MIDDLE of the fibril is different
result=kruskal_treatments(df=mids)
signif=dunns_treatments(df=mids, comparisons_list=comparisons_list)


signif.to_csv(f'{stats_output}stats_kruskal_dunns_mids.csv')






ends=nomalised_intensity[nomalised_intensity['end_or_middle']=='END']
result=kruskal_treatments(df=ends)
signif=dunns_treatments(df=ends, comparisons_list=comparisons_list)


signif.to_csv(f'{stats_output}stats_kruskal_dunns_ends.csv')


#now need to loop through each treatment and compare end and middle



def kruskal_endvmid(df):
    end=df[df['end_or_middle']=='END']
    end=np.array(end['rel_size'])
    
    middle=df[df['end_or_middle']=='MIDDLE']
    middle=np.array(middle['rel_size'])
    #this part tells us whether there is a difference in some of the means of these data sets.
    result = stats.kruskal(end, middle)
    print(result)
    return result

def dunns_endvmid(df):
    data = [df[df['end_or_middle'] == 'END']['rel_size'], 
        df[df['end_or_middle'] == 'MIDDLE']['rel_size']]

    p_values = sp.posthoc_dunn(data, p_adjust='holm')

    signif=p_values < 0.05
        

    return signif


comparisons_list=['end', 'middle']
kruskal_store=[]
Dunns_store=[]
for treatment, df in nomalised_intensity.groupby('Treatment'):
    treatment
    df
    result=kruskal_endvmid(df)
    t=pd.DataFrame(result).T.rename(columns={0:'stat', 1:'pval'})
    names_rows={0:f'{treatment}'}
    t.rename(index=names_rows, inplace=True)
    kruskal_store.append(t)
    if t['pval'].values[0] < 0.05:
        p_values=dunns_endvmid(df=df)
        p_values.columns=['end','middle']
        p_values['Treatment'] = treatment
        names_rows={
            1:'end',
            2:'middle',
        }
        p_values.rename(index=names_rows, inplace=True)

        Dunns_store.append(p_values)


kruskal_end_mid=pd.concat(kruskal_store)
dunns_end_mid=pd.concat(Dunns_store)

kruskal_end_mid.to_csv(f'{stats_output}kruskal_wallis_end_v_middle.csv')

dunns_end_mid.to_csv(f'{stats_output}dunns_posthoc_end_v_middle.csv')



#that was a lot!
#now do an anova on the barplot data!!! this is the means in the blue plot. do this for ends df, and mids df. then repeat for each treatment, end v mid
def anova_treatments(df):

    #NOW TRY AN ANOVA? FOR THE DIF OF THE MEANS
    #df=df[df['Treatment']!='JB1-A8-110-1nM-']
    #df=df[df['Experiment_number']!='Experiment 90-2']
    #dropping replicates I don't want
    #result =df.drop(df[(df['Treatment']=='JB1-A8-SOD-0.5uM-') & (df['Experiment_number']=='Experiment 98-1')].index)
    result=df
    #this part tells us whether there is a difference in some of the means of these data sets.

    stats.f_oneway(result['rel_size'][result['Treatment'] == 'JB1-A8-SOD-0.5uM-'],
                result['rel_size'][result['Treatment']
                                            == 'JB1-A8-110-2uM-'],
                result['rel_size'][result['Treatment']
                                            == 'JB1-A8-110-5nM-'],
                result['rel_size'][result['Treatment']
                                            == 'JB1-A8-110-0.5uM-'],
                result['rel_size'][result['Treatment']
                                            == 'JB1-A8-'],
                result['rel_size'][result['Treatment']
                                            == 'A8-ATP-'],
                result['rel_size'][result['Treatment'] == 'A8-noATP-'])


    #now we have a dataframe where every row is a timepoint, and compares every other timepoint to itself. each column the one being compared moves across. To get a simpler output, we can just say p_values < 0.05, and any entry that is statistically significant will be 'True'.
    comp = mc.MultiComparison(result['rel_size'], result['Treatment'])
    post_hoc_res = comp.tukeyhsd()
    summary=pd.DataFrame(post_hoc_res.summary())
    return summary

def anova_endmid(df, treatment):

    #NOW TRY AN ANOVA? FOR THE DIF OF THE MEANS
    #df=df[df['Treatment']!='JB1-A8-110-1nM-']
    #df=df[df['Experiment_number']!='Experiment 90-2']
    #dropping replicates I don't want
    #result =df.drop(df[(df['Treatment']=='JB1-A8-SOD-0.5uM-') & (df['Experiment_number']=='Experiment 98-1')].index)
    result=df
    #this part tells us whether there is a difference in some of the means of these data sets.

    x=stats.f_oneway(result['rel_size'][result['end_or_middle'] == 'END'],
                result['rel_size'][result['end_or_middle']
                                            == 'MIDDLE'])
    y=pd.DataFrame(x).T.rename(columns={1:'pval'})
    y['Treatment']= treatment
   # if y['pval'].values[0] < 0.05:


    return y


summary_ends=anova_treatments(df=ends)
summary_mids=anova_treatments(df=mids)

anova_store=[]
posthoc=[]
for treatment, df in nomalised_intensity.groupby('Treatment'):
    treatment
    df
    anova=anova_endmid(df, treatment)
    anova_store.append(anova)
    if anova['pval'].values[0] <0.05:
        comp = mc.MultiComparison(df['rel_size'], df['end_or_middle'])
        post_hoc_res = comp.tukeyhsd()
        summary=pd.DataFrame(post_hoc_res.summary())
        summary['Treatment']=treatment
    
    
        posthoc.append(summary)

anova=pd.concat(anova_store)
tukeys=pd.concat(posthoc)

summary_ends.to_csv(f'{stats_output}molsize_treatments_end_anova.csv')
summary_mids.to_csv(f'{stats_output}molsize_treatments_mid_anova.csv')