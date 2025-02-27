import DEA
import argparse
### Make this to get as parameters the experiments to do the DEA, the metadata of 
# the experiments and the path to save the results, the comparison groups. Check that the comparisons are onle "young", "old" or "middle"
# and that the metadata has the columns

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run DEA')
    parser.add_argument('--experiment', type=str, help='Experiment name')
    parser.add_argument('--metadata', type=str, help='Metadata file')
    parser.add_argument('--save_folder', type=str, help='Path to save the results')
    parser
    args = parser.parse_args()
    # Create DEA object
    dea = DEA.DEA()
    # Run DEA
    dea.get_DEA()