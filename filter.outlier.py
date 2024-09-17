import pandas as pd
import argparse
import os
import numpy as np

# Define the function to filter outliers based on log-transformed p-values
def filter_outliers(input_file, output_dir):
    # Ensure the input file is located within the specified output directory
    input_file_path = os.path.join(output_dir, input_file)

    # Check if the file exists
    if not os.path.exists(input_file_path):
        print(f"Error: Input file {input_file_path} not found.")
        return

    # Read the file, assuming it's tab-delimited
    df = pd.read_csv(input_file_path, delimiter='\t', header=0)

    # Check if the correct columns are present
    if df.shape[1] < 4:
        print(f"Error: Expected at least 4 columns in the input file.")
        return

    # Assuming the 4th column is 'lrt-pvalue'
    p_values = df.iloc[:, 3]  # Column index starts at 0, so column 3 is the 4th column

    # Apply -log10 transformation to p-values to work in log scale
    log_p_values = -np.log10(p_values)

    # IQR method on the log-transformed values to find high outliers
    Q1 = log_p_values.quantile(0.25)
    Q3 = log_p_values.quantile(0.75)
    IQR = Q3 - Q1

    print(f"Log-transformed p-values (first 5):\n{log_p_values.head()}")
    print(f"Q1: {Q1}, Q3: {Q3}, and IQR: {IQR}")

    # Define the upper bound for normal values in the log scale
    upper_bound_log = Q3 + 1.5 * IQR

    # Filter out high outliers in log scale (which correspond to significant k-mers)
    high_outliers_df = df[log_p_values > upper_bound_log]

    # Keep normal and low values in log scale
    filtered_df = df[log_p_values <= upper_bound_log]

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Strip leading/trailing whitespace from all columns before saving
    filtered_df = filtered_df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    high_outliers_df = high_outliers_df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    
    # Save the filtered data (normal + low values) in the output directory
    filtered_file = os.path.join(output_dir, 'filtered_kmers_normal_and_low.txt')
    high_outliers_file = os.path.join(output_dir, 'high_outliers_kmers.txt')
    
    filtered_df.to_csv(filtered_file, sep='\t', index=False)
    high_outliers_df.to_csv(high_outliers_file, sep='\t', index=False)

    print(f"Saved {len(high_outliers_df)} high outliers to {high_outliers_file}")
    print(f"Saved {len(filtered_df)} normal/low values to {filtered_file}")

# Define the main function to handle argument parsing
def main():
    parser = argparse.ArgumentParser(description="Filter high outliers in the log-transformed 4th column and save them separately.")
    
    # Define required arguments
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Input file name located in the output directory')
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='Directory to read input and save output files')

    # Parse the arguments
    args = parser.parse_args()

    # Call the function to filter outliers
    filter_outliers(args.input_file, args.output_dir)

if __name__ == "__main__":
    main()
