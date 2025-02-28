## This is tu run the TF analysis. It will ask as parameter the path to the yml file with the comparisons and the path to save the results
# as for run_DDS. It will also ask for pvalue, threshold, interest, and species='human', but these are optional.
# Here we will check if the pavalue is > 0 and <=1, if the threshold is >=0, and if the interest is in ['stat', 'log2FoldChange']
# but only if they are specifies, if no, no argument will be passed to the class.
import argparse
import TF
import pandas as pd
import os


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run DEA')
    parser.add_argument('--dds_files', type=str, help='YAML file with the names of the comparisons and the path to the DDS file')
    parser.add_argument('--save_name', type=str, help='Path to save the results')
    parser.add_argument('--pvalue', type=float, help='Pvalue to filter the results', default=0.05)
    parser.add_argument('--threshold', type=float, help='Threshold to filter the results', default=None)
    parser.add_argument('--interest', type=str, help='Interest to filter the results', default='stat')
    parser.add_argument('--species', type=str, help='Species to filter the results', default='human')

    args = parser.parse_args()
    if args.pvalue is not None:
        if args.pvalue<=0 or args.pvalue>1:
            raise ValueError('Pvalue must be >0 and <=1')
    else:
        args.pvalue = 0.05
    if args.threshold is not None:
        if args.threshold<0:
            raise ValueError('Threshold must be >=0')
    else:
        args.threshold = 0
    if args.interest is not None:
        if args.interest not in ['stat', 'log2FoldChange']:
            raise ValueError('Interest must be in ["stat", "log2FoldChange"]')
    else:
        args.interest = 'stat'
    tf = TF.TranFact(dds_files_path=args.dds_files, pvalue=args.pvalue, threshold=args.threshold, interest=args.interest, species=args.species)
    tf.save_tf_df_of_combined_comparisons(args.save_name)
    print('Done')
