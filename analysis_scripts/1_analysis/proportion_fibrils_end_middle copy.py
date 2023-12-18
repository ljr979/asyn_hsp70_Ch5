import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

input_folder='python_results/FC2/coords_added/'
output_folder='python_results/FC2/coords_added/end_vs_middle/'


#this script aims to find the proportion of the regions we have defined on each fibril that are bound by a hspa8 molecule
#THIS IS DIFFERENT TO % OF HPSA8 BOUND TO FIBRILS, becaause instead of finding ALL the hspa8 molecules that are bound to a fibril, and saying how many of THEM are at the end or the middle? IT is instead asking- of the FIBRILS that are bound, what percentage of alllllll the ends that exist in this entire treatment, are bound compared to all of the middles

#look for end_vs_middle.csv files
if not os.path.exists(output_folder):
        os.makedirs(output_folder)

yikes_=[[f'{root}/{name}' for name in files if 'end_vs_middle.csv' in name]for root, dirs, files in os.walk(f'{input_folder}/')]
yikes_=[item for sublist in yikes_ for item in sublist if 'jb1' not in item]



def gather (yikes_):
    end_vs_middle=[]
    for item in yikes_:
            df=pd.read_csv(f'{item}')
            item=item.replace('//', '/')
            item=item.replace('\\', '/')
            
            treatment=item.split('/')[-1].split('_')[-4]
            df['treatment']=treatment
            end_vs_middle.append(df)

    end_vs_middle=pd.concat(end_vs_middle)
    return end_vs_middle

def find_proportion_end_middle(df):

    #find the total sum of the LENGTH column (in pixels)
    num_fibs=len(df['ID_length_conc_image'])
    total_lengths=sum(df['Length'])
    #tfind total ENDS lengths in pixels by multiplying end pixels (4 or 6 pixels) by the number of fibrils 
    both_ends=4
    all_ends=num_fibs*both_ends
    #total middle pixels is the total length minus total ends
    all_mids=total_lengths-all_ends
    ##find total NUMBER OF CHAPERONES at an end and total number at a middle pixel
    mid=len(df[df['Where']=='MIDDLE'])
    end=mid=len(df[df['Where']=='END'])

    #then DIVIDE THE NUMBER of chaperones at end lengths by the number of total end pixels
    end_perc=end/all_ends*100
    #divide number of chaperones at mid lengths by number of mid pixels
    mid_perc=mid/all_mids*100
    alls=treatment, total_lengths, end_perc, mid_perc

    return alls





end_vs_middle=gather(yikes_)


proportion=[]
for treatment, df in end_vs_middle.groupby('treatment'):
    treatment
    df
    t=find_proportion_end_middle(df=df)
    proportion.append(t)


proportion=pd.DataFrame(proportion, columns=['treatment', 'total_length_fibs', 'end_perc', 'middle_perc'])

order=['Experiment87-A8-B1-excess', '+hsp110-2h']

conc_plot=pd.melt(proportion, id_vars=['treatment' ], value_vars=['end_perc','middle_perc'], var_name="end or middle", value_name='percentage of total fibril pixels bound at each location')
f=sns.barplot(data=conc_plot, x='treatment', y='percentage of total fibril pixels bound at each location', order=order, hue='end or middle', palette='crest', alpha=0.5, edgecolor='grey')
plt.title('Percentage of colocalised fibril pixels @ end vs middle ')
plt.savefig(f'{output_folder}/percent_pixels_endmiddle.svg')
plt.savefig(f'{output_folder}/percent_pixels_endmiddle.png')
plt.show()
