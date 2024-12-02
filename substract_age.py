import os
import pandas as pd
from scipy import stats
from multiprocessing import Pool

def remove_trend(y, x):
    slope, _, _, _, _ = stats.linregress(x, y)
    return y - slope * x

def substract_age(file_path):
    try:
        data = pd.read_pickle(file_path)
        x = data['Age']
        data.iloc[:, :-2] = data.iloc[:, :-2].apply(remove_trend, x=x)
        return data

    except ValueError as e:
        print(f'Error processing {file_path}: {e}')
        return None


def process_and_save(file_path):
    processed_df = substract_age(file_path)
    if processed_df is not None:
        output_path = os.path.join(substracted_age, os.path.basename(file_path))
        processed_df.to_pickle(output_path)
        print(f'Processed and saved: {output_path}')

if __name__ == '__main__':
    initial_datasets = '/tank/projects/vpalagina_hobotnica/hobotnica/data_imputed_with_meta_for_regression_pkl'
    substracted_age = '/tank/projects/vpalagina_hobotnica/hobotnica/data_imputed_with_meta_substracted_age'

    # Get a list of all files in the dataset directory
    files = [os.path.join(initial_datasets, filename) for filename in os.listdir(initial_datasets)]

    # Create a multiprocessing Pool
    with Pool(processes=os.cpu_count()) as pool:
        # Map the files to the processing function
        pool.map(process_and_save, files)
