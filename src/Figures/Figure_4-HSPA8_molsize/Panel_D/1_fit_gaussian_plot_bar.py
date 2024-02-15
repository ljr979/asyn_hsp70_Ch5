import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
from lmfit import models

def fit_single_gauss(df, treatment, mu_1, sigma_1, amplitude_1):
    filt_df = df[df['Treatment'] == treatment]
    bins = np.arange(0, 0.4, 0.025)
    inds = np.digitize(filt_df['rel_intens'].astype(float), bins)
    xdata, ydata = np.unique(inds, return_counts=True)
    ydata = ydata[1:-1] #### trim off outside range bins at the end
    xdata = [np.mean(bins[x : x + 2]) for x in range(len(bins)- 1)]  ##### convert bin edges to bin centres, therefore end up with one less bin
    sns.lineplot(xdata, ydata)
    model_1 = models.GaussianModel(prefix='m1_')
    model = model_1 
    model_1.set_param_hint('m1_center', vary=True)

    params_1 = model_1.make_params(center = mu_1, sigma = sigma_1, amplitude = amplitude_1, min = 0)
    params = params_1.update(params_1)
    #params = params.update(params_3)
    output = model.fit((ydata/np.max(ydata)), params, x=xdata)
    fig = output.plot(data_kws={'markersize': 3})
    paramaters = {name:output.params[name].value for name in output.params.keys()}
    fitx = np.arange(-0.2, 1.2, 0.025)
    fit1 = model_1.eval(x = fitx, center = paramaters['m1_center'], amplitude = abs(paramaters['m1_amplitude']), sigma = paramaters['m1_sigma'])


    sns.lineplot(fitx, fit1)

    #sns.lineplot(fitx, fit3)
    plt.show()
    # Calculate area under the curve for each gaussian
    aoc_m1 = paramaters['m1_amplitude']


    sum_aoc = aoc_m1 
    aoc_m1_percent_of_total = (aoc_m1/sum_aoc)*100

    #aoc_m3_percent_of_total = (aoc_m3/sum_aoc)*100
    list_of_gaus_proportion = [aoc_m1_percent_of_total]
    labels_of_gaus_proportion = ['m1']
    proportion_df = pd.DataFrame([labels_of_gaus_proportion, list_of_gaus_proportion])
    proportion_df.columns = proportion_df.iloc[0]
    proportion_df = proportion_df.drop(0)
    proportion_df['treatment'] = treatment
    proportion_df.to_csv(f'{output_folder}/single_gaussian_proportions_for_{treatment}.csv')
    return proportion_df, paramaters

def three_fit_gauss_dif_constrained_nativespont(df, treatment, mu_1, sigma_1, amplitude_1, mu_2, sigma_2, amplitude_2, mu_3, sigma_3, amplitude_3):
    filt_df = df[df['Treatment'] == treatment]
    bins = np.arange(0.025, 0.9, 0.025)
    inds = np.digitize(filt_df['rel_intens'].astype(float), bins)
    xdata, ydata = np.unique(inds, return_counts=True)
    ydata = ydata[1:-1] #### trim off outside range bins at the end
    xdata = [np.mean(bins[x : x + 2]) for x in range(len(bins)- 1)]  ##### convert bin edges to bin centres, therefore end up with one less bin
    sns.lineplot(xdata, ydata)
    model_1 = models.GaussianModel(prefix='m1_')
    model_2 = models.GaussianModel(prefix='m2_')
    model_3 = models.GaussianModel(prefix='m3_')
    model = model_1 + model_2 + model_3
    model_1.set_param_hint('m1_center', vary=True)
    model_1.set_param_hint('m1_sigma', vary=False)
    params_1 = model_1.make_params(center = mu_1, sigma = sigma_1, amplitude = amplitude_1, min = 0)
    params_2 = model_2.make_params(center = mu_2, sigma = sigma_2, amplitude = amplitude_2, min = 0)
    params_3 = model_3.make_params(center = mu_3, sigma = sigma_3, amplitude = amplitude_3, min = 0)



    params = params_1.update(params_2)
    params = params.update(params_3)
    output = model.fit((ydata/np.max(ydata)), params, x=xdata)
    fig = output.plot(data_kws={'markersize': 3})
    paramaters = {name:output.params[name].value for name in output.params.keys()}
    fitx = np.arange(-0.2, 1.2, 0.025)
    fit1 = model_1.eval(x = fitx, center = paramaters['m1_center'], amplitude = abs(paramaters['m1_amplitude']), sigma = paramaters['m1_sigma'])
    fit2 = model_2.eval(x = fitx, center = paramaters['m2_center'], amplitude = abs(paramaters['m2_amplitude']), sigma = paramaters['m2_sigma'], fwhm = paramaters['m2_fwhm'])
    fit3 = model_3.eval(x = fitx, center = paramaters['m3_center'], amplitude = abs(paramaters['m3_amplitude']), sigma = paramaters['m3_sigma'] )
    sns.lineplot(fitx, fit1, color='red')
    sns.lineplot(fitx, fit2, color='green')
    sns.lineplot(fitx, fit3, color= 'blue')
    plt.title(f'{treatment}')
    plt.show()
    # Calculate area under the curve for each gaussian
    aoc_m1 = paramaters['m1_amplitude']
    aoc_m2 = paramaters['m2_amplitude']
    aoc_m3 = paramaters['m3_amplitude']
    # aoc_m1 = (paramaters['m1_amplitude']*paramaters['m1_sigma'])/0.3989
    # aoc_m2 = (paramaters['m2_amplitude']*paramaters['m2_sigma'])/0.3989
    # aoc_m3 = (paramaters['m3_amplitude']*paramaters['m3_sigma'])/0.3989
    sum_aoc = aoc_m1 + aoc_m2 + aoc_m3
    aoc_m1_percent_of_total = (aoc_m1/sum_aoc)*100
    aoc_m2_percent_of_total = (aoc_m2/sum_aoc)*100
    aoc_m3_percent_of_total = (aoc_m3/sum_aoc)*100
    list_of_gaus_proportion = [aoc_m1_percent_of_total, aoc_m2_percent_of_total, aoc_m3_percent_of_total]
    labels_of_gaus_proportion = ['m1', 'm2', 'm3']
    proportion_df = pd.DataFrame([labels_of_gaus_proportion, list_of_gaus_proportion])
    proportion_df.columns = proportion_df.iloc[0]
    proportion_df = proportion_df.drop(0)
    proportion_df['treatment'] = treatment
    proportion_df.to_csv(f'{output_folder}/three_gaussian_proportions_for_{treatment}.csv')
    return proportion_df, parameters


if __name__ == "__main__":
    input_folder = 'data/Figures/Figure_4/Panel_B/'
    output_folder = 'data/Figures/Figure_4/Panel_D/'
    #read in normalised molsize data
    compiled_df = pd.read_csv(f'{input_folder}1_normed_mol_per_foci.csv', header="infer")

    collated = []
    #the below way fits a single skewed gaussian distribution for the A8 only treatment, which is used to guide the multiple gaussian populations in the later fitting steps
    single_gauss_set, parameters = fit_single_gauss(compiled_df, 'A8-ATP-', 0.1, 0.009, 1)
    #append to list
    collated.append(single_gauss_set)
    #define parameters for fitting the A8 curve in subsequent treatments, based on initial fit
    mu_1 = parameters['m1_center']
    sigma_1 = parameters['m1_sigma']
    amplitude_1 = parameters['m1_amplitude']

    #dictionary guiding the parameters for the 2nd and third gaussian distributions within each treatment.
    para_dict = {

        'JB1-A8-':[0.3, 0.05, 1, 0.5, .05, 1],
        'JB1-A8-110-0.5uM-':[0.3, 0.05, 1, 0.5, .05, 1],
        'JB1-A8-SOD-0.5uM-':[0.3, 0.05,1, 0.5, .05, 1]
    }
    #treatments to loop over
    treats = ['JB1-A8-', 'JB1-A8-110-0.5uM-', 'JB1-A8-SOD-0.5uM-']
    #loop through and fit each of the datasets
    for treat in treats:
        treat
        three_gaussian_kj_skew_con_nat, p = three_fit_gauss_dif_constrained_nativespont(compiled_df, treat, mu_1, sigma_1, amplitude_1, para_dict[treat][0],para_dict[treat][1],para_dict[treat][2],para_dict[treat][3],para_dict[treat][4],para_dict[treat][5])
        collated.append(three_gaussian_kj_skew_con_nat)
    collated = pd.concat(collated)
    #save the collated percentage of area under the curve that is assigned to each population
    collated.to_csv(f'{output_folder}gaussians_collated.csv')


    colors = plt.cm.BuPu(np.linspace(0.4, 1, 4))
    collated.plot(
        x='treatment', 
        y=['m1', 'm2', 'm3'], 
        kind='bar',
        stacked=True,
        figsize=(10, 6),
        color=colors, 
        title=f"A8 molecules in each population (%)",
        linewidth=2, 
        edgecolor='0',
        alpha=0.8
        )
    plt.legend(loc="upper left", ncol=3, title='Population')
    plt.xlabel("Treatment")
    plt.ylabel("Molecules in each peak(%)")
    plt.show()













#this also is not working? wtf
gaussian_kj_skew_con_nat = fit_gauss_dif_constrained_nativespont(compiled_df, 'JB1-A8-110-1nM-',mu_1, sigma_1, amplitude_1, 0.2, .05, 1)
collated.append(gaussian_kj_skew_con_nat)


gaussian_kj_skew_con_nat = three_fit_gauss_dif_constrained_nativespont(compiled_df, 'JB1-A8-110-1nM-', mu_1, sigma_1, amplitude_1,  0.2, 0.05 ,2, 0.5, .05, 1)

gaussian_kj_skew_con_nat = fit_single_gauss2(compiled_df, 'JB1-A8-110-1nM-', mu_1, sigma_1, amplitude_1)

