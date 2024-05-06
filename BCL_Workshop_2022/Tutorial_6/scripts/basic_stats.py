#! /usr/bin/env python
import os
import matplotlib as mpl
import numpy as np
import math
import random
from scipy import stats
from sklearn.metrics import mean_squared_error
if os.environ.get('DISPLAY', '') == '':
    mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser

####################################
# Janky old stats script
#####################################

def main():
    
    # Read in command-line arguments
    args = cmdlineparse()
    if args.has_id_labels:
        pred, exp = np.loadtxt(args.input, dtype="float", delimiter=",", unpack=True, usecols=(1,2))
    else:
        pred, exp = np.loadtxt(args.input, dtype="float", delimiter=",", unpack=True, usecols=(0,1))
        
    # Make cross-validation splits
    pred_exp_full = [(x,y) for x,y in zip(pred, exp)]
    if args.const_random_seed:
        random.shuffle(pred_exp_full, int(args.const_random_seed))
    else:
        random.shuffle(pred_exp_full)
    if int(args.cv) > int(0):
        n_pairs_per_chunk = int(len(pred_exp_full) / int(args.cv))
        pred_exp_full = [pred_exp_full[i * n_pairs_per_chunk:(i + 1) * n_pairs_per_chunk] for i in range((len(pred_exp_full) + n_pairs_per_chunk - 1) / n_pairs_per_chunk )]
    
    # Open output file
    ofile = open(args.output, "w+")
    ofile.write("CASFv2016 score function metrics \n")
    
    if args.R:
        R, R_p = CalcR(pred, exp, args)
        ofile.write("Correlation coefficient: \n")
        ofile.write("R: " + str(round(R,3)) + "\n")
        
        # Compute metric for all CV splits
        if int(args.cv) > int(0):
            ofile.write("R CV: ")
            for cv_i in range(0, int(args.cv)):
                R, R_p = CalcR(np.transpose(pred_exp_full[cv_i])[0], np.transpose(pred_exp_full[cv_i])[1], args)
                ofile.write(str(round(R,3)) + ",")
            ofile.write("\n")
                
        ofile.write("p-value: " + str(round(R_p,3)) + "\n")
    
    if args.rho:
        rho, rho_p = CalcRho(pred, exp)
        ofile.write("Ranked correlation coefficient: \n")
        ofile.write("rho: " + str(round(rho,3)) + "\n")
        
        # Compute metric for all CV splits
        if int(args.cv) > int(0):
            ofile.write("rho CV: ")
            for cv_i in range(0, int(args.cv)):
                rho, rho_p = CalcRho(np.transpose(pred_exp_full[cv_i])[0], np.transpose(pred_exp_full[cv_i])[1])
                ofile.write(str(round(rho,3)) + ",")
            ofile.write("\n")
        
        ofile.write("p-value: " + str(round(rho_p,3)) + "\n")

    if args.tau:
        tau, tau_p = CalcKendallTau(pred, exp)
        ofile.write("Kendall's tau: \n")
        ofile.write("tau: " + str(round(tau,3)) + "\n")
        
        # Compute metric for all CV splits
        if int(args.cv) > int(0):
            ofile.write("Kendall's tau CV: ")
            for cv_i in range(0, int(args.cv)):
                tau, tau_p = CalcKendallTau(np.transpose(pred_exp_full[cv_i])[0], np.transpose(pred_exp_full[cv_i])[1])
                ofile.write(str(round(tau,3)) + ",")
            ofile.write("\n")
        
        ofile.write("p-value: " + str(round(tau_p,3)) + "\n")
    
    if args.mae:
        mae, mad, nmae = CalcMAE(pred, exp)
        ofile.write("Mean Absolute Error: \n")
        ofile.write("MAE: " + str(round(mae,3)) + "\n")
        ofile.write("Mean Absolute Deviation: \n")
        ofile.write("MAD: " + str(round(mad,3)) + "\n")
        ofile.write("Normalized Mean Absolute Error: \n")
        ofile.write("NMAE: " + str(round(nmae,3)) + "\n")
        
        # Compute metric for all CV splits
        if int(args.cv) > int(0):
            mae_l = []
            mad_l = []
            nmae_l = []
            for cv_i in range(0, int(args.cv)):
                mae, mad, nmae = CalcMAE(np.transpose(pred_exp_full[cv_i])[0], np.transpose(pred_exp_full[cv_i])[1])
                mae_l.append(mae)
                mad_l.append(mad)
                nmae_l.append(nmae)
            ofile.write("MAE CV: ")
            for cv_i in range(0, int(args.cv)):
                ofile.write(str(round(mae_l[cv_i],3)) + ",")
            ofile.write("\n")
            ofile.write("MAD CV: ")
            for cv_i in range(0, int(args.cv)):
                ofile.write(str(round(mad_l[cv_i],3)) + ",")
            ofile.write("\n")
            ofile.write("NMAE CV: ")
            for cv_i in range(0, int(args.cv)):
                ofile.write(str(round(nmae_l[cv_i],3)) + ",")
            ofile.write("\n")
        
    if args.rmse:
        RMSE = CalcRMSE(pred, exp)
        ofile.write("Root Mean Squared Error: \n")
        ofile.write("RMSE: " + str(round(RMSE,3)) + "\n")
        
        # Compute metric for all CV splits
        if int(args.cv) > int(0):
            ofile.write("Root Mean Squared Error CV: ")
            for cv_i in range(0, int(args.cv)):
                rmse = CalcRMSE(np.transpose(pred_exp_full[cv_i])[0], np.transpose(pred_exp_full[cv_i])[1])
                ofile.write(str(round(rmse,3)) + ",")
            ofile.write("\n")
        
    return 0

def cmdlineparse():
    parser = ArgumentParser(description="command line arguments")
    parser.add_argument("-input", dest="input", required=True, 
                        help="Three-column CSV file where the far left column is the group id, followed by the predicted Ki (or Kd) values and then the actual values. \n\
                        Note that this is just the format of the bcl.exe model:Test output file.", metavar="<input>")
    parser.add_argument("-output", dest="output", required=True, help="Name of output file", metavar="<output file>")
    parser.add_argument("-R", dest="R",required=False,default=False,action="store_true",help="Enable to compute the Pearson correlation (scoring power metric)")
    parser.add_argument("-rho", dest="rho",required=False,default=False,action="store_true",help="Enable to compute Spearman rank correlation (ranking power metric)")
    parser.add_argument("-mae", dest="mae",required=False,default=False,action="store_true",help="Enable to compute the mean absolute error")
    parser.add_argument("-rmse", dest="rmse",required=False,default=False,action="store_true",help="Enable to compute the root-mean-squared error")
    parser.add_argument("-tau",dest="tau",required=False,default=False,action="store_true",help="Enable to compute Kendall's tau")
    parser.add_argument("-cv", dest="cv", required=False, default=0, help="Perform N-fold random-split cross-validation of the input set in addition to computing statistics across the whole set")
    parser.add_argument("-const_random_seed", dest="const_random_seed", required=False, 
                        help="Pass a constant random seed to use for randomization in cross-validation splits")
    parser.add_argument("-has_id_labels", dest="has_id_labels",required=False,default=False,action="store_true",help="Get the correct columns for the calculation")
    parser.add_argument("-plot", dest="plot", required=False, default=False, action="store_true", help="Plot results")
    parser.add_argument("-x_range", dest="x_range", required=False,help="Set x-axis range", nargs=2, metavar="<x range>")
    parser.add_argument("-y_range", dest="y_range", required=False,help="Set y-axis range", nargs=2, metavar="<y range>")
    parser.add_argument("-x_label", dest="x_label", required=False,default="Experimental", help="Set x-axis label", metavar="<x label>")
    parser.add_argument("-y_label", dest="y_label", required=False,default="Predicted", help="Set y-axis label", metavar="<y label>")
    parser.add_argument("-font_size", dest="font_size", required=False,default=8, help="Font size", metavar="<font size>")
    args = parser.parse_args()
    return args

def CalcR(predicted, experimental, args):

    # Scoring power part 1: Compute Pearson Correlation Coefficient
    # Initialize component piece 0 values for calculation of R
    numerator = float(0.0)
    denominator_pred = float(0.0)
    denominator_exp = float(0.0)

    # Compute means of each dataset
    pred_ave = predicted.mean()
    exp_ave = experimental.mean()

    # Determine all component pieces
    for row in range(0, len(predicted)):
        numerator += (predicted[row] - pred_ave) * (experimental[row] - exp_ave)
        denominator_pred += np.square((predicted[row] - pred_ave))
        denominator_exp += np.square((experimental[row] - exp_ave))
    denominator_pred = np.sqrt(denominator_pred)
    denominator_exp = np.sqrt(denominator_exp)

    # Final calculation of R
    R = float(numerator / (denominator_pred * denominator_exp))

    # Just use scipy to compute R
    R_scipy, p_scipy = stats.pearsonr(predicted, experimental)

    # Scoring power part 2: Compute linear regression and standard deviation
    # Make heavy use of scipy for this part - get the regression
    m, b, r_value, p_value, std_err = stats.linregress(predicted, experimental)

    # Compute the standard deviation in linear regression
    SD = float(0.0)
    for row in range(0, len(predicted)):
        SD += np.square(experimental[row] - (m * predicted[row] + b))
    SD = float(np.sqrt(SD / (len(predicted) - 1)))

    # Plot the curve
    if args.plot:
        ofile = open(str(args.output)+".png", "w+")

        # Set global plot options
        font = {'family': 'sans-serif',
                'weight': 'normal',
                'size': 8}
        plt.rc('font', **font)
        plt.rcParams["axes.edgecolor"] = "black"
        plt.rcParams["axes.linewidth"] = 1

        # Plot data
        #plt.xlabel("Experimental")
        #plt.ylabel("Predicted")

        # Plot the data
        sns.set(color_codes=True)
        sns.set_style("whitegrid")

        points = plt.scatter(experimental, predicted, marker='o', color='lime', edgecolors='black')
        ax = sns.regplot(x=experimental, y=predicted, scatter=False, color=".1", line_kws={'lw': 1.5})
        ax.grid(True)
        ax.set_facecolor('white')
        ax.axis('on')
        ax.patch.set_edgecolor('black')
        ax.patch.set_linewidth('1')
        # ax.set_xticks(np.arange(experimental))
        # ax.set_yticks(np.arange(predicted))

        # Set axes options
        ax.set_xlabel(args.x_label, fontsize=args.font_size)
        ax.set_ylabel(args.y_label, fontsize=args.font_size)
        if args.x_range:
            ax.set_xlim([float(args.x_range[0]), float(args.x_range[1])])
        if args.y_range:
            ax.set_ylim([float(args.y_range[0]), float(args.y_range[1])])

        # Save the figure to output and be done
        plt.savefig(str(args.output)+".png", bbox_inches='tight', dpi=600)
        ofile.close()

    return R_scipy, p_scipy

def CalcRho(predicted, experimental):
    # Ranking power part 1: Compute Spearman's rank correlation coefficient
    # Just use SciPy

    rho_scipy, rho_p_scipy = stats.spearmanr(predicted, experimental)
    return rho_scipy, rho_p_scipy

def CalcMAE(predicted, experimental):
    # Compute the MAE
    ae = np.zeros(len(predicted))
    for row in range(0, len(predicted)):
        ae[row] = np.sqrt(np.square(predicted[row] - experimental[row]))
    mae = np.mean(ae)

    # Compute the MAD
    ad = []
    mean = float(np.mean(experimental))
    for row in range(0, len(experimental)):
        ad.append(np.sqrt(np.square(experimental[row] - mean)))
    mad = np.mean(ad)
    return mae, mad, mae / mad

def CalcRMSE(predicted, experimental):
    rmse = math.sqrt(mean_squared_error(experimental, predicted))
    return rmse

def CalcKendallTau(predicted, experimental):
    tau, p = stats.kendalltau(predicted, experimental)
    return tau, p
    
# Execute main() only if run directly
if __name__ == '__main__':
    main()
