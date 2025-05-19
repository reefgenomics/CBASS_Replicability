"""
Script to see how good a predictor one set of CBASS samples are of the other.
"""

from collections import defaultdict
import pandas as pd
import itertools
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import statistics
from math import comb
import time
import random
import scipy.stats

class gsPilot:
    # NB We tried median and it doesn't work. I have left the code in place but will not develop it further and do 
    # not intend for it to be run.
    # method: can be 'range' or 'median'

    # split: can be True or False. If True then the 40 set is split into first and second 20
    
    # experiment: can be 'predictive' or 'barshis'. If 'predictive', we ask whether the each of the
    # colonies is below or above the "threshold" (which is separate for each of the replicates) in both replicates.
    # If it is, it scores a 1, if it is not it scores a 0. E.g. we calculate the threshold for R1 and for R2, then for each
    # colony we look to see if the R1 is above or below the R1 threshold, if it is below, then it will score a 1 if the R2
    # representative of that colony is below the R2 threshold. If the R2 is above the R2 threshold (i.e. they disagree), then a 0 is scored.
    # If the experiement is 'barshis', then there is no threshold considered and rather we are splitting the points of each
    # replicate run arbitrarily into the top half and bottom half (we only work with even sample numbers). Then we ask the question 
    # what proportion of the colonies agree with the top or bottom half assignment.
    # I.e. for colony 1, if in R1 it is a top half and in R2 it is a bottom half (i.e. they disagree), then a 0 is scored. If they agree, a 1 is scored.
    # NB Originally I was reporting this as the average proportion of correct for each sub sample, but then Chris asked for the actual
    # number of time that there was perfect agreement i.e. for a subsample of 10 samples, how many times the top 5 matched exactly in R1 and R2.
    # These two methods of reporting are now proportion_of_top_half_agreement and proportion_of_absolute_agreement respectively.
    # Chris then asked to see the proportion_of_absolute_agreement for not only perfect match cases (i.e. all top 5 agree) but
    # also 0.8, 0.6, 0.4, 0.2, and 0. This is implemented using the plot parameter (see below). If 1_point then I only show the
    # cases where there were perfect matches of the top halfs in both replicates. If 'line_series' then it is for the other levels of
    # agreement. I.e. agreement of 0.8 would mean that 4 of the top 5 agree in R1 and R2 (for the 10 sample subset example).

    # reporting (only for barshis): can be 'proportion_of_absolute_agreement' or 'proportion_of_top_half_agreement'
    # proportion_of_absolute_agreement reports the proportion of the number of times the subsets showed agreement at or above the given level
    # proportion_of_top_half_agreement reports the average proportions of each subset that showed agreement. (see code calculation for exact definitions)

    # plots = '1_point' or 'line_series': only makes difference if experiment is 'barshis'
    # if 1_point then only the values for exact matches I.e. value of 1 are plotted, else a series of lines are plotted for 1, 0.8, 0.6, 0.4, 0.2, 0
    def __init__(self, method="range", split=True, experiment='predictive', reporting="proportion_of_absolute_agreement", plots='1_point'):
        self.method = method
        self.reporting = reporting
        self.plots = plots

        df = pd.read_csv("pilot_data.csv")
        self.df40 = df[["40_genet", "40_ED50_1", "40_ED50_2"]]
        self.df40.columns = ["genet", "ED50_1", "ED50_2"]
        self.df40.set_index("genet", drop=True, inplace=True)
        
        self.df16 = df[["16_genet", "16_ED50_1", "16_ED50_2"]].dropna()
        self.df16.columns = ["genet", "ED50_1", "ED50_2"]
        self.df16.set_index("genet", drop=True, inplace=True)
        self.split = split
        self.experiment = experiment
        if self.split:
            self.df40_1 =self.df40.iloc[:20,]
            self.df40_2 =self.df40.iloc[20:,]
            self.data_sets = [self.df40_1, self.df40_2, self.df16]
            self.data_set_names = ["40 1-20", "40 21-40", "16"]
            fig, ax_arr = plt.subplots(nrows=1, ncols=2)
            ax = ax_arr[0]
            x = self.df40_1["ED50_1"]
            y = self.df40_1["ED50_2"]
            ax.scatter(x, y)
            ax.set_title("40 1-20")
            ax.set_xlim(34,38.5)
            ax.set_ylim(34,38.5)
            m, b = np.polyfit(x, y, 1)
            ax.plot(x, m*x + b)
            slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
            rsquared = r_value**2

            ax = ax_arr[1]
            x = self.df40_2["ED50_1"]
            y = self.df40_2["ED50_2"]
            ax.scatter(x, y)
            ax.set_title("40 21-40")
            ax.set_xlim(34,38.5)
            ax.set_ylim(34,38.5)
            m, b = np.polyfit(x, y, 1)
            slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
            rsquared = r_value**2
            ax.plot(x, m*x + b)

            

            plt.tight_layout()
            plt.savefig("split_data_scatter.png")
            plt.close()
        else:
            fig, ax = plt.subplots(nrows=1, ncols=1)
            
            x = self.df40["ED50_1"]
            y = self.df40["ED50_2"]
            ax.scatter(x, y)
            ax.set_title("40")
            ax.set_xlim(34,38.5)
            ax.set_ylim(34,38.5)
            m, b = np.polyfit(x, y, 1)
            ax.plot(x, m*x + b)
            slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(self.df40["ED50_1"], self.df40["ED50_2"])
            rsquared = r_value**2
            plt.tight_layout()
            plt.savefig("non_split_data_scatter.png")
            plt.close()


            self.data_sets = [self.df40, self.df16]
            self.data_set_names = ["40", "16"]


    def _compute(self):
        # For both the df16 dataset and the df40 dataset
        data_set_results = []
        fig, ax_arr = plt.subplots(nrows=1, ncols=len(self.data_sets))
        for i, ds in enumerate(self.data_sets):
            # For each number of samples starting at all samples and going down to 5
            # samples
            

            num_samples_to_num_subsets = defaultdict(int)
            if self.experiment == "predictive":
                num_samples = list(range(len(ds.index), 3, -1))
                num_sample_result_dd = defaultdict(list)
            elif self.experiment == "barshis":
                if len(ds.index) % 2 == 0:
                    # Then we have an even number and we can start with all samples
                    num_samples = list(range(len(ds.index), 3, -2))
                else:
                    # we have an odd number of samples and so we start with one less
                    num_samples = list(range(len(ds.index)-1, 3, -2))
            
                if self.reporting == "proportion_of_top_half_agreement":
                    num_sample_result_dd = defaultdict(list)
                elif self.reporting == "proportion_of_absolute_agreement":
                    num_sample_result_dd = defaultdict(list)
                    results_df = pd.DataFrame(columns=[0, .2, .4, .6, .8, 1], index=num_samples)

            for sample_num in num_samples:
                # For each of the combinations of samples
                absolute_agreement_count = 0
                subset_len = comb(len(ds.index), sample_num)
                print(f"Sample number {sample_num}: {subset_len} subsets")
                start = time.time()
                if subset_len > 1000:
                    indices_as_list = list(ds.index.values)
                    subsets = []
                    while len(subsets) < 1000:
                        new_set = set(random.sample(indices_as_list, sample_num))
                        if new_set not in subsets:
                            subsets.append(new_set)
                        # subsets = random.sample(list(itertools.combinations(ds.index.values, sample_num)), 1000)
                else:
                    subsets = itertools.combinations(ds.index.values, sample_num)
                
                for subset in subsets:
                    num_samples_to_num_subsets[sample_num] += 1
                    # Calculate the dividing threshold
                    # Currently the dividng threshold is  (max-min)/2
                    if self.experiment == "predictive":
                        ser = ds.loc[subset,]
                        min_val_R1 = ser["ED50_1"].min()
                        min_val_R2 = ser["ED50_2"].min()
                        if self.method == "range":
                            threshold_R1 = ((ser["ED50_1"].max() - min_val_R1) / 2) + min_val_R1
                            threshold_R2 = ((ser["ED50_2"].max() - min_val_R2) / 2) + min_val_R2

                            results = [1 if (val_1 > threshold_R1) == (val_2 > threshold_R2) else 0 for val_1, val_2 in zip(ser["ED50_1"], ser["ED50_2"])]
                        elif self.method == "median":
                            threshold_R1 = ser["ED50_1"].median()
                            threshold_R2 = ser["ED50_2"].median()
                            results = [1 if (val_1 >= threshold_R1) == (val_2 >= threshold_R2) else 0 for val_1, val_2 in zip(ser["ED50_1"], ser["ED50_2"])]
                        # Now for each each of the pairs of values
                        # we simply want the proportion of the pairs are found on the same side of the threshold
                        num_sample_result_dd[sample_num].append(sum(results) / len(results))
                    elif self.experiment == "barshis":
                        ser = ds.loc[subset,]
                        ordered_ser_R1 = [_[0] for _ in sorted(ser["ED50_1"].items(), key=lambda x: x[1], reverse=True)]
                        ordered_ser_R2 = [_[0] for _ in sorted(ser["ED50_2"].items(), key=lambda x: x[1], reverse=True)]
                        # Get the top n/2 for each, make a set and see how many values are found in common.
                        # You get a proportion from the number found in common
                        half_n = int(len(subset)/2)
                        set_R1 = set(ordered_ser_R1[:half_n])
                        set_R2 = set(ordered_ser_R2[:half_n])
                        if self.reporting == "proportion_of_top_half_agreement":
                            prop = len(set_R1.intersection(set_R2))/half_n
                            num_sample_result_dd[sample_num].append(prop)
                        elif self.reporting == "proportion_of_absolute_agreement":
                            prop = len(set_R1.intersection(set_R2))/half_n
                            num_sample_result_dd[sample_num].append(prop)
                    else:
                        raise RuntimeError("unknown experiment type")
                if self.reporting == "proportion_of_absolute_agreement" and self.experiment == 'barshis':
                    # here we populate one row of the df
                    for proportion in results_df.columns:
                        # Value should be number of proportion => the proportion
                        results_df.at[sample_num, proportion] = len([_ for _ in num_sample_result_dd[sample_num] if _ >= proportion])/len(num_sample_result_dd[sample_num])    
                print(f"done: took {time.time() - start}s")
            # At this point we're have the raw data. Now we want to plot these up as means as a product of number of samples
            if self.reporting == "proportion_of_absolute_agreement" and self.experiment == "barshis":
                if self.plots == "line_series":
                    results_df.plot.line(ax=ax_arr[i])
                    # wins_props = [num_sample_result_dd[_] for _ in num_samples]
                    # ax_arr[i].scatter(x=num_samples, y=wins_props, c="black")
                    for num_sample in num_samples:
                        ax_arr[i].text(s=num_samples_to_num_subsets[num_sample], x=num_sample, y=1.1, ha="center", va="center", rotation="vertical", fontsize=8)
                else:
                    ax_arr[i].scatter(x=num_samples, y=results_df[1], c="black")
                    for num_sample in num_samples:
                        ax_arr[i].text(s=num_samples_to_num_subsets[num_sample], x=num_sample, y=results_df.at[num_sample, 1] + 0.1, ha="center", va="center", rotation="vertical", fontsize=8)
            else:
                data_lists = [num_sample_result_dd[samp_num] for samp_num in num_samples]
                sample_averages = [sum(_)/len(_) for _ in data_lists]
                stdevs = [statistics.pstdev(_) if len(_) > 1 else 0 for _ in data_lists]
                ax_arr[i].errorbar(x=num_samples, y=sample_averages, yerr = stdevs, fmt='o', c="black")
            ax_arr[i].set_xlim(num_samples[0] +1, num_samples[-1] - 1)
            ax_arr[i].set_ylim(0,1.2)
            if self.experiment == "predictive":
                ax_arr[i].set_ylabel("proportion correct predictions")
            elif self.experiment == "barshis":
                if self.reporting == "proportion_of_top_half_agreement":
                    ax_arr[i].set_ylabel("proportion top half agreement")
                elif self.reporting == "proportion_of_absolute_agreement":
                    ax_arr[i].set_ylabel("proportion of absolute agreement")
            ax_arr[i].set_xlabel("number of samples")
            ax_arr[i].set_title(self.data_set_names[i])
        if self.experiment == 'barshis':
            fig.suptitle(f"predictions_split_{self.split}_experiment_{self.experiment}_reporting_{self.reporting}")
            plt.tight_layout()
            plt.savefig(f"predictions_split_{self.split}_experiment_{self.experiment}_reporting_{self.reporting}.png")
        elif self.experiment == 'predictive':
            fig.suptitle(f"predictions_split_{self.split}_method_{self.method}_experiment_{self.experiment}")
            plt.tight_layout()
            plt.savefig(f"predictions_split_{self.split}_method_{self.method}_experiment_{self.experiment}.png")
        

gsPilot()._compute()

