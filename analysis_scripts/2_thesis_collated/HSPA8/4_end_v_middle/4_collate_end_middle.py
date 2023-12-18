
from enum import unique
import os, re
import threading
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
import math

#output_folder='data/0_colocalisation/'
#defining a dictionary which tells us the experiment number in an easy word as the key, and the PATH to that repo as the value matching the key to call on later

paths = {
#A8 no ATP
('Experiment 75-1','A8-noATP-1', 'old'):'D:/Experiments/Experiment_75/python_results/1/',

('Experiment 75-2','A8-noATP-2', 'old'):'D:/Experiments/Experiment_75/python_results/2/',


#A8 + ATP
('Experiment 72-1','A8-ATP-1', 'old'):'D:/Experiments/Experiment_72/python_results/1/',
('Experiment 72-2','A8-ATP-2', 'old'):'D:/Experiments/Experiment_72/python_results/2/',
('Experiment 84-1','A8-ATP-3','old'):'D:/Experiments/Experiment_84/python_results/FC1/',


#JB1+A8+ATP
#('Experiment 66-2','JB1-A8-1'):'D:',
('Experiment 73-2','JB1-A8-2', 'old'):'D:/Experiments/Experiment_73/python_results/2/',
('Experiment 77-2','JB1-A8-3', 'old'):'D:/Experiments/Experiment_77/python_results/2/',
('Experiment 77-1','JB1-A8-4', 'old'):'D:/Experiments/Experiment_77/python_results/1/',
('Experiment 76-1','JB1-A8-5', 'old'):'D:/Experiments/Experiment_76/python_results/1/',
('Experiment 76-2','JB1-A8-6', 'old'):'D:/Experiments/Experiment_76/python_results/2/',
('Experiment 87-2','JB1-A8-7','new'):'D:/Experiments/Experiment87-FC2/python_results/',
('Experiment 90-1','JB1-A8-8','new'):'D:/Experiments/Experiment_90-FC1/python_results/',('Experiment 90-2','JB1-A8-9','new'):'D:/Experiments/Experiment_90-FC2/python_results/FC2/',
('Experiment 85-3','JB1-A8-0','new'):'D:/Experiments/Experiment85-FC3/python_results/',


#MIX DARK A8 2UM & LABELLED A8
('Experiment 86-1','darkA8-647A8-1','new'):'D:/Experiments/Experiment86-FC1/python_results/',

('Experiment 87-1','darkA8-647A8-2','new'):'D:/Experiments/Experiment87-FC1/python_results/',


#JB1+A8+1nM hsp110-excessA8
('Experiment 81-1','JB1-A8-110-1nM-1', 'old'):'D:/Experiments/Experiment_81/python_results/1/',
('Experiment 81-2','JB1-A8-110-1nM-2', 'old'):'D:/Experiments/Experiment_81/python_results/2/',
('Experiment 84-3','JB1-A8-110-1nM-3','new'):'D:/Experiments/Experiment_84/python_results/FC3/',
('Experiment 90-2','JB1-A8-110-1nM-4','new'):'D:/Experiments/Experiment_90-FC2/python_results/FC2/',
#('Experiment 66-3','JB1-A8-110-1nM-5'):'D:',


#JB1+A8+0.5uM hsp110- excess A8
('Experiment 84-2','JB1-A8-110-0.5uM-1','new'):'D:/Experiments/Experiment_84/python_results/FC2/',
('Experiment 86-2','JB1-A8-110-0.5uM-2','new'):'D:/Experiments/Experiment86-FC2-bleaches/python_results/',
('Experiment 85-3','JB1-A8-110-0.5uM-3','new'):'D:/Experiments/Experiment85-FC3/python_results/',

#JB1+A8+2uM hsp110- excess A8
('Experiment 90-2','JB1-A8-110-2uM-1','new'):'D:/Experiments/Experiment_90-FC2/python_results/FC2/',

#JB1+A8+0.5uM hsp110- NO excess A8
('Experiment 87-2','JB1-110-0.5uM-1','new'):'D:/Experiments/Experiment87-FC2/python_results/',
('Experiment 90-1','JB1-110-0.5uM-2','new'):'D:/Experiments/Experiment_90-FC1/python_results/',
('Experiment 90-1','JB1-110-2uM-1','new'):'D:/Experiments/Experiment_90-FC1/python_results/',

}

files=[]
for key, value in paths.items():
    coloc_files =[[f'{root}/{filename}' for filename in files if 'middle_end_lengths_df' in filename] for root, dirs, files in os.walk(f'{value}')]
    
    coloc_files=[item for sublist in coloc_files for item in sublist ]
 
    df=pd.DataFrame(coloc_files, columns=['path'])
    df[['Experiment_number', 'Treatment-replicate', 'analysis-version']]=key
    files.append(df)

collated_collated_colocal=pd.concat(files)
old=collated_collated_colocal[collated_collated_colocal['analysis-version']=='old']
new=collated_collated_colocal[collated_collated_colocal['analysis-version']=='new']
#from collate colocalisation


new_collated=[]
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

            treatment_from_file=path.split('/')[-3]
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
            

     if version=='new':

        for (path,Treatment), row in df.groupby(['path', 'Treatment-replicate']):
            path
            


            path=path.replace('\\','/')

            replicate_t=Treatment
            replicate=replicate_t.split('-')[-1]

            treatment_from_file=path.split('/')[-3]
            treatment=Treatment.rstrip(replicate_t[-1])
            Exp_num=list(row['Experiment_number'])[0]

            row
            if path not in paths_read:
                print(Exp_num)
                coloc_new=pd.read_csv(path)
                paths_read.append(path)
                fresh_df=coloc_new
                protein=fresh_df['protein'].tolist()[0]
                # if protein=='JB1':

                #         print(Exp_num)
                #         if Exp_num == 'Experiment 85-3':
                #             print('shit')
                if 'DNAJB1_middle_end_lengths_df' not in path:

                    fresh_df['Treatment_from_file']=treatment_from_file
                    fresh_df['Experiment_number']=Exp_num
                    fresh_df['path']=path
                    fresh_df['replicate']=replicate
                    fresh_df['Treatment']=treatment
                    new_collated.append(fresh_df)

new_collated=pd.concat(new_collated)
old_collated=pd.concat(old_collated)


end_mid=pd.concat([new_collated, old_collated])
end_mid.drop([col for col in end_mid.columns.tolist() if 'Class' in col], axis=1, inplace=True)
end_mid.drop([col for col in end_mid.columns.tolist() if 'Unnamed: 0' in col], axis=1, inplace=True)
end_mid.drop([col for col in end_mid.columns.tolist() if 'Unnamed: 0.1' in col], axis=1, inplace=True)
end_mid.drop([col for col in end_mid.columns.tolist() if 'Unnamed: 0.1.1' in col], axis=1, inplace=True)
end_mid.drop([col for col in end_mid.columns.tolist() if 'Frame' in col], axis=1, inplace=True)
g=pd.Series(end_mid['path'].unique()).value_counts()

#the lines bellow (commented out, uncomment when you need them) allow us to check for proteins, and then see which experiments they came from so I can trace back in my lab book at repositories to see whether it's been labelled the correct protein, and rerun if need be. 
# end_mid['protein'].unique()
# t=end_mid[end_mid['protein']=='JB1']
# t['Experiment_number'].unique()
# t['path'].unique()
# q=end_mid[end_mid['protein']=='HSPA8']
# q['Experiment_number'].unique()
#s=end_mid[end_mid['Experiment_number']=='Experiment 85-3']

#end_mid is now the concatinated 'middle_end_lengths_df' from each experiment. 
#now need to run the rest of the 'df for plotting' script here to make sure they're all the same. Had to get them from this point and not earlier, because the earlier steps in the script need information about the zoom of the microscope which can change day to day, and that is accounted for already in these dataframes

# also need to filter concentrations the same way that I did in mols/subunit analysis.
treatments_to_filter=list(end_mid['concentration'].unique()) 


filter_list=[
 'A8-B1-0.5HSP110-t1h',
 'A8-JB1-excess',
 'A8-dark1uM-A8-647-10nM',
 'Experiment86-A8-B1-excess',
 '+hsp110-2h',
 'Experiment87-A8-B1-excess',
 'A8-JB1-110-0.5uM-flowout',
 'A8-JB1-110-1nM-flowout',

 'Experiment90-A8-B1-0.5hsp110',
 'Experiment90-A8-B1-2uMhsp110',

 'Experiment90-A8-B1-1nMhsp110',
 'Experiment90-A8-B1',

 '10nM',
 'darkJB1',
 'noATP-10nM',
 'mixed',
 'HSPA8-JB1-110',
 'excess_A8']

end_mid_filter=end_mid[end_mid['concentration'].isin(filter_list)]
#now need to match up the concentration label with the treatment it is supposed to be
fix_treatment={
 'A8-JB1-excess':'JB1-A8-', 
 'Experiment86-A8-B1-excess':'JB1-A8-',
 ' Experiment87-A8-B1-excess':'JB1-A8-',
 'darkJB1':'JB1-A8-',
 'A8-B1-0.5HSP110-t1h':'JB1-A8-110-0.5uM-', 
 'A8-JB1-110-0.5uM-flowout':'JB1-A8-110-0.5uM-', 
 'A8-dark1uM-A8-647-10nM':'darkA8-647A8-', 
 '+hsp110-2h':'JB1-110-0.5uM-',
 '+hsp110-2h':'JB1-110-0.5uM-',
 'A8-JB1-110-1nM-flowout':'JB1-A8-110-1nM-',
'Experiment90-A8-B1-1nMhsp110':'JB1-A8-110-1nM-',
 'HSPA8-JB1-110':'JB1-A8-110-1nM-',
 '10nM':'A8-ATP-',
 'excess_A8':'A8-ATP-',
 'noATP-10nM':'A8-noATP-'}
end_mid_filter['fix_label_treatment']=end_mid_filter['concentration'].map(fix_treatment)




fixup_ninetytwo={ 'Experiment90-A8-B1-1nMhsp110':'JB1-A8-110-1nM-',
 'Experiment90-A8-B1-2uMhsp110':'JB1-A8-110-2uM-',
 'Experiment90-A8-B1':'JB1-A8-'}
fixup_ninetyone={ 'Experiment90-A8-B1-0.5hsp110':'JB1-110-0.5uM',
 'Experiment90-A8-B1-2uMhsp110':'JB1-110-2uM',}
checklist=['Experiment90-2','Experiment90-1']
new=[]
for (Experiment, Treatment), row in end_mid_filter.groupby(['Experiment_number', 'Treatment']):

    Experiment
    Treatment
    row
    replicate=row['replicate'].unique()
    if Experiment not in checklist:
        new.append(row)
    if Experiment == 'Experiment 90-2':
        row['fix_label_treatment']=row['concentration'].map(fixup_ninetytwo)
        new.append(row)
    if Experiment == 'Experiment 90-1':
        row['fix_label_treatment']=row['concentration'].map(fixup_ninetytwo)
        new.append(row)
    
        
new=pd.concat(new)
new.drop([col for col in new.columns.tolist() if 'Treatment' in col], axis=1, inplace=True)
new=new.rename(columns={'fix_label_treatment':'Treatment'})

df_melt=pd.melt(new, id_vars=['Contour_ID',
 'Pos.',
 'X',
 'Y',
 'Length',
 'distance',
 'Point_X',
 'Point_Y',
 'concentration',
 'protein',
 'image_number',
 'ID_length_conc_image',
 'new_ID_hspX_hspY_num',
 'EndX1, EndX2, EndY1, EndY2',
 'EndX1',
 'EndX2',
 'EndY1',
 'EndY2',
 'dist_1',
 'dist_2',
 'Where',
 'end_length_pixels',
 'end_length_um',
 'middle_length_um',
 'number_of_chaps_END',
 'number_of_chaps_MIDDLE',
 'Experiment_number',
 'path',
 'replicate',
 'Treatment'], value_vars=['chaps_per_end_length','chaps_per_middle_length'], var_name="end_or_middle", value_name='Chaperones_per_length')




data_output='data/4_end_middle/'
if not os.path.exists(data_output):
         os.makedirs(data_output)
results_output='results/4_end_middle/'
if not os.path.exists(results_output):
         os.makedirs(results_output)
new.to_csv(f'{data_output}end_middle_filtered.csv')



data=df_melt[df_melt['Chaperones_per_length']>0]

palette='Blues_r'
protein='HSPA8'
legendtitle='Location'
plot_title=f'Location of {protein} on fibril'
#order_of_experiment=[
# 'darkA8-647A8-',
#  'A8-noATP-',
#  'A8-ATP-',
#  'JB1-A8-',
#  'JB1-A8-110-1nM-',
#  'JB1-A8-110-0.5uM-',
#  'JB1-A8-110-2uM-',
#  'JB1-110-0.5uM-',
#  'JB1-110-2uM-',
#  ]
output_folder=results_output
g=sns.barplot(data=data, x="Treatment", y="Chaperones_per_length", hue="end_or_middle", 
#order=order_of_experiment, 
palette=palette, alpha=.6, edgecolor="grey")
#g.despine(left=True)
sns.stripplot(x="Treatment", 
y="Chaperones_per_length", 
hue="end_or_middle", #order=order_of_experiment, 
data=data, palette='Blues', dodge=True, alpha=0.1, ax=g
)
#g.set_axis_labels(f"Concentration of {protein} titrated onto fibrils", f"# of {protein} bound per unit fibril length (um)")
plt.xticks(rotation=90)
plt.ylim(0,10)
plt.legend(title=f'{legendtitle}', loc='upper right')
plt.title(f'{plot_title}')
plt.tight_layout()
plt.savefig(f'{output_folder}End_vs_middle_{protein}.svg')
plt.savefig(f'{output_folder}End_vs_middle_{protein}.png')


def plot_ends(condition, conditions_list, output_folder, protein, palette, legendtitle, plot_title, treat):
        
    g=sns.barplot(data=condition, x="Treatment", y="Chaperones_per_length",         #hue="end_or_middle",
        order=conditions_list, palette=palette, alpha=.6, edgecolor="grey")
    
    #g.despine(left=True)
    sns.stripplot(x="Treatment", 
    y="Chaperones_per_length", 
    #hue="end_or_middle", 
    # order=order_of_experiment, 
    data=condition, order=conditions_list, palette=palette, dodge=True, alpha=0.1, ax=g
    )
    h,l = g.get_legend_handles_labels()
    plt.legend(h[0:2],l[0:2], title=f'{legendtitle}',bbox_to_anchor=(1.05, 1),loc='upper right')

    #g.set_axis_labels(f"Concentration of {protein} titrated onto fibrils", f"# of {protein} bound per unit fibril length (um)")
    plt.xticks(rotation=90)
    plt.ylabel('# HSPA8 foci per region')
    plt.ylim(0,10)

    plt.title(f'{plot_title}')
    plt.tight_layout()
    plt.savefig(f'{output_folder}{treat}_End_vs_middle_{protein}.svg')
    plt.savefig(f'{output_folder}{treat}_End_vs_middle_{protein}.png')

#plotting separately
#this stuff remains the same for all treatments
palette='Blues_r'
protein='HSPA8'
legendtitle='Location'
plot_title=f'Location of {protein} on fibril'
output_folder='results/4_end_middle/'
change_labels={'chaps_per_end_length':'End', 'chaps_per_middle_length':'Middle'}
data['end_or_middle']=data['end_or_middle'].map(change_labels)
#this part changes depending on treatment
conditions_list=['A8-noATP-']
treat='A8_noATP'
condition=data[data['Treatment'].isin(conditions_list)]

plot_ends(condition, conditions_list, output_folder, protein, palette, legendtitle, plot_title, treat)





conditions_list=['A8-noATP-','A8-ATP-']
treat='+A8_ATP'
condition=data[data['Treatment'].isin(conditions_list)]

plot_ends(condition, conditions_list, output_folder, protein, palette, legendtitle, plot_title, treat)






#+jb1
conditions_list=['A8-noATP-','A8-ATP-','JB1-A8-']
treat='b1_a8_atp'
condition=data[data['Treatment'].isin(conditions_list)]


plot_ends(condition, conditions_list, output_folder, protein, palette, legendtitle, plot_title, treat)



#+110 1nM
conditions_list=['A8-noATP-','A8-ATP-','JB1-A8-','JB1-A8-110-1nM-']
treat='b1_a8__110_1nm_atp'
condition=data[data['Treatment'].isin(conditions_list)]



plot_ends(condition, conditions_list, output_folder, protein, palette, legendtitle, plot_title, treat)






#+110 500nM
conditions_list=['A8-noATP-','A8-ATP-','JB1-A8-','JB1-A8-110-1nM-', 'JB1-A8-110-0.5uM-']
treat='b1_a8__110_0.5nm_atp'
condition=data[data['Treatment'].isin(conditions_list)]



plot_ends(condition, conditions_list, output_folder, protein, palette, legendtitle, plot_title, treat)



#PLOT FOR EMBO
palette = {'A8-ATP-': '#f46d43', 'JB1-A8-': '#fee999', 'JB1-A8-110-1nM-': '#e3f399', 'JB1-A8-110-0.5uM-': '#9dd569'}
#palette = 'Oranges'
protein = 'HSPA8'
legendtitle = 'Location'
plot_title = f'Location of {protein} on fibril'
output_folder = 'D:/presentations_conferences/2023/EMBO/figs_from_fibrils_summary_repo/'

change_labels = {'chaps_per_end_length': 'End',
                 'chaps_per_middle_length': 'Middle'}
data['end_or_middle'] = data['end_or_middle'].map(change_labels)

conditions_list=['A8-ATP-','JB1-A8-','JB1-A8-110-1nM-', 'JB1-A8-110-0.5uM-']
treat='b1_a8__110_0.5nm_atp_for_poster'
condition=data[data['Treatment'].isin(conditions_list)]



ax=sns.catplot(
    data=condition, x="Treatment", y="Chaperones_per_length", col="end_or_middle",
    kind="bar", height=4, aspect=.6, palette=palette
plt.xticks(rotation=90)
plt.ylabel('# HSPA8 foci per region')
plt.ylim(0, 4)

plt.title(f'{plot_title}')
#plt.tight_layout()
plt.savefig(f'{output_folder}{treat}_End_vs_middle_{protein}.svg')
plt.savefig(f'{output_folder}{treat}_End_vs_middle_{protein}.png')
#plot_ends(condition, conditions_list, output_folder, protein, palette, legendtitle, plot_title, treat)










def make_summary(data, summary_df):
    
    for timepoint, df in data.groupby('plot_labels'): 
        timepoint
        mean_end=df['number_of_chaps_END'].mean()
        sem_end=df['number_of_chaps_END'].sem()
        mean_MIDDLE=df['number_of_chaps_MIDDLE'].mean()
        sem_MIDDLE=df['number_of_chaps_MIDDLE'].sem()
        summary_df[timepoint]=[mean_end,sem_end,mean_MIDDLE,sem_MIDDLE]

    summary_df=pd.DataFrame(summary_df.items(), columns=['timepoint', 'meansemENDmeansemMIDDLE'])
    summary_df[['mean_end', 'sem_END','mean_middle', 'sem_MIDDLE']]=pd.DataFrame(summary_df.meansemENDmeansemMIDDLE.tolist(), index= summary_df.index)
    return summary_df



end_mid_all_summary=[]

for key, value in paths.items():
    Experiment=key[0]
    treatment=key[1]
    top_path=f'{value}'
    python_results=f'{top_path}python_results/'
    replicates=[folder for folder in os.listdir(f'{python_results}')]
    for replicate in replicates: 
        
        coloc= pd.read_csv(f'{python_results}{replicate}/end_vs_middle/step_two/middle_end_lengths_df.csv')
        coloc['replicate']=replicate

        coloc['Experiment']=Experiment
        coloc['treatment']=treatment
        end_mid_all_summary.append(coloc)
end_mid_all_summary=pd.concat(end_mid_all_summary)
end_mid_summary=end_mid_all_summary[end_mid_all_summary['concentration'].isin(filter_treatments)]

t=end_mid_summary[end_mid_summary['replicate']=='3_110']
t.replace('DARKJB1-flowA8', 'darkb1-flowA8-flow110', inplace=True)
end_mid_summary=end_mid_summary[end_mid_summary['replicate']!='3_110']
end_mid_summary=end_mid_summary.append(t)

plot_label_dict={'HSPA8_titration':'A8',
 'DARKJB1-flowA8': 'J+A-flow',
 'A8-noATP':'-ATP',
 'A8_titration-darkJB1':'J+A-flow',
 'premix-darkjb1-A8':'J+A-mix',
 'premix-darkjb1-A8-flowonhsp110-100min':'J+A-flow110-100min',
 'premix-hsp110-A8-JB1':'J+A+110-mix',
 'darkb1-flowA8-flow110': 'J+A-flow+110-60min'}

end_mid_summary['plot_labels']=end_mid_summary['treatment'].map(plot_label_dict)


summary_df={}
data=end_mid_summary
end_mid_summary=make_summary(data, summary_df)
end_mid_summary.to_csv(f'data/4_end_middle/all_end_middle_summary.csv')








def plotting_chaps_per_um_end_middle(data, palette, protein, legendtitle, plot_title, summary):
    melt_summary=pd.melt(summary, id_vars=['timepoint'], value_vars=['sem_END','sem_MIDDLE'], var_name="end_or_middle", value_name='sem')
    #g=sns.catplot(data=data, kind='bar', x="concentration", y="Chaperones_per_length", hue="end_or_middle", 
    #order=['0-0', '0min', '20min', '40min', '60min', '80min'], palette=palette, 
           #alpha=.6, 
          # )

    g=sns.barplot(x='plot_labels', y="Chaperones_per_length", data=data, ci=None, palette=palette, hue='end_or_middle', alpha=0.45, edgecolor='black')
    plt.ylabel(f"# of {protein} bound per unit fibril length (um)")
    plt.xlabel(f"Time after hsp110 flowed onto fibrils")

    #sns.stripplot(x='concentration', y="Chaperones_per_length", data=data, palette='Purples', hue='end_or_middle',dodge=True, alpha=0.7, size=4, edgecolor='Darkgrey', linewidth=1, ax=g)
    x_coords = [p.get_x() + 0.5 * p.get_width() for p in g.patches]
    y_coords = [p.get_height() for p in g.patches]
    g.errorbar(x=x_coords, y=y_coords, yerr=melt_summary["sem"], fmt="none",  elinewidth = 1, capsize = 4, color = 'black')
    #g.despine(left=True)


    plt.xticks(rotation=45)
    plt.ylim(0,2)
    plt.legend(title=f'{legendtitle}', labels=['End', 'Middle'], loc='upper center')
    plt.title(f'{plot_title}')
    plt.tight_layout()
    plt.savefig(f'{output_folder}End_vs_middle.png')


output_folder='results/4_end_middle/'

two_colours=['#d20f46', '#d39ca2']
sns.set_palette(sns.color_palette(two_colours))
palette=two_colours
protein= 'HSPA8'
legendtitle='Location of molecule'
plot_title='HSPA8 foci per um fibril length'
data=end_mid
summary=end_mid_summary

plotting_chaps_per_um_end_middle(data, palette, protein, legendtitle, plot_title, summary)





output='results/2_num_foci'
ax = sns.violinplot(x='plot_labels', y='foci_per_pixel', data=end_mid,scale='width', palette='coolwarm_r')
plt.title(f'foci per length unit fibril (pixel)')
plt.ylabel('foci (per pixel)')

plt.xlabel(f'treatment')
plt.xticks(rotation=45)
plt.savefig(f'{output}/foci_per_length_unit_colocalised.png')

