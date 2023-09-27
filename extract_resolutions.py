import argparse
import pandas as pd
import numpy as np
import os
import re

# Calculate the average resolution while handling empty arrays
def calculate_average_resolution(resolutions):
    if isinstance(resolutions, list) and len(resolutions) == 0:
        return np.nan  # Handle empty list by returning NaN
    else:
        return np.mean(resolutions)

# Calculate the median resolution while handling empty arrays
def calculate_median_resolution(resolutions):
    if isinstance(resolutions, list) and len(resolutions) == 0:
        return np.nan  # Handle empty list by returning NaN
    else:
        return np.median(resolutions)



def get_sample_name(file_name):
# Define the pattern
    pattern = r'^(.*?)_Run\d+_updated_coords_(.*?)\.bedpe$'

    # Example string that matches the pattern

    # Use re.match to extract values
    match = re.match(pattern, file_name)

    if match:
        prefix = match.group(1)  # Extract the prefix
        sample_id = match.group(2)  # Extract the sample_id
        #print(f"Prefix: {prefix}")
        #print(f"Sample ID: {sample_id}")
    else:
        print("No match found.",file_name )
        sample_id = ""
    return sample_id

def modify_arr(arr):
    if(len(arr) ==1):
        a = str(arr[0])
        modify_arr = a.replace('kb','')
        modify_arr = modify_arr.replace('kb','')
        modify_arr = modify_arr.replace('Mb','000')
        #print("modified:",modify_arr,"original:",a)
        return int(modify_arr)
        
    else:
        modify_arr = []
        for a in arr:
            a = str(a)
            a=a.replace('kb','')
            a=a.replace('Kb','')
            a= a.replace('Mb','000')
            modify_arr.append(int(a))
        return modify_arr    

def split_sample_name(df, column_name):
    # Initialize empty lists for 'Sample_name' and 'Replicate'
    sample_name_only = []
    sample_rep = []

    # Loop through the input column in the DataFrame
    for input_string in df[column_name]:
        # Split the string by "_Rep" if it exists, or use the whole string
        parts = input_string.split("_Rep")
        sample_id = parts[0]
        replicate = "rep" + parts[1] if len(parts) > 1 else "rep1"
        
        # Append the extracted values to the respective lists
        sample_name_only.append(sample_id)
        sample_rep.append(replicate)

    # Create new columns 'Sample_name' and 'Replicate' in the DataFrame
    df['Sample_name'] = sample_name_only
    df['Replicate'] = sample_rep

    return df

def add_run_to_validation_column(df):
    # Check if 'Validation Run' column exists
    if 'Validation Run' not in df.columns:
        df['Validation Run'] = 'Validation Run'  # Create the column if it doesn't exist

    # Concatenate 'Validation Run' with 'Run' for each row
    df['Validation Run'] += df['Run'].astype(str)

    return df


def extract_resolutions(start_dir):
    dataset = []

    # Check if the directory exists
    if os.path.exists(start_dir) and os.path.isdir(start_dir):
        # List the contents of the directory
        for i in range(1,10):
            print("Processing Run",i)
            contents_protean = os.listdir(start_dir+"Run"+str(i)+"/Protean_samples")
            contents_arima = os.listdir(start_dir+"Run"+str(i)+"/Arima_samples")

            matching_files = [item for item in contents_arima if item.endswith(".bedpe")]
            # Print the list of matching files
            for file_name in contents_arima:
                if file_name.endswith(".bedpe"):
                    df=pd.read_csv(start_dir+"Run"+str(i)+"/Arima_samples/"+file_name,sep="\t")
                    #print(file_name, df.shape)
                    resolutions = []
                    for j in range(df.shape[0]):
                        resolutions.append(df[df.columns[8]][j])
                        
                    dataset.append([i,file_name, "Arima",resolutions[:]])
                
            matching_files = [item for item in contents_protean if item.endswith(".bedpe")]
            # Print the list of matching files
            for file_name in matching_files:
                df=pd.read_csv(start_dir+"Run"+str(i)+"/Protean_samples/"+file_name,sep="\t")
                #print(file_name, df.shape)
                resolutions = []
                for j in range(df.shape[0]):
                    resolutions.append(df[df.columns[8]][j])
                        
                dataset.append([i,file_name, "Protean",resolutions[:]])

            #print("\n\n\n")

    else:
        print(f"The directory {start_dir} does not exist.")
    final_dataset = pd.DataFrame(dataset)
    final_dataset.columns = ["Run","File","Source","Resolutions"]
    return final_dataset


def main():
    """ This script filters a GTF file for all genes that have at least one transcript with a biotype.
    The complete list of biotypes is passed and the desired biotypes are marked by a boolean Y or N.
    The new gtf file attributes exon_id, gene_id and transcript_id are modified by removing the versions to match Cellranger IDs.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-dir",
        "-i",
        dest="input_dir",
        default=None,
        required=True,
        help="input directory containing the manually curated bedpe files",
    )
    parser.add_argument(
        "--outputname",
        "-o",
        dest="output_name",
        default=None,
        help="output csv file name"
    )
    args = parser.parse_args()
    # Extract the resolutions
    print("Extracting resolutions")
    dataset = extract_resolutions(args.input_dir)
    dataset['Resolutions'] = dataset['Resolutions'].apply(modify_arr)
    
    # Apply the functions to create the 'average_resolution' and 'median_resolution' columns
    dataset['average_resolution'] = dataset['Resolutions'].apply(calculate_average_resolution)
    dataset['median_resolution'] = dataset['Resolutions'].apply(calculate_median_resolution)

    # Extract the sample name from the file name    
    dataset['Sample_name'] = dataset['File'].apply(get_sample_name)

    dataset = split_sample_name(dataset, 'Sample_name')
    dataset = add_run_to_validation_column(dataset)

    # save the dataset into a csv file
    print("Saving the dataset")
    dataset.to_csv(args.output_name,index=False)


if __name__ == "__main__":
    main()