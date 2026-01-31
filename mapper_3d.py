import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import gget
import requests
import os
from Bio.PDB.DSSP import make_dssp_dict
import time


def download_af_model(uniprot_id, output_dir):
    # Variables for API
    file_name = f"AF-{uniprot_id}-F1-model_v6.cif"
    url = f"https://alphafold.ebi.ac.uk/files/{file_name}"
    output_path = os.path.join(output_dir, file_name)

    # Check if the file exists
    if os.path.exists(output_path):
        print("File already exists")
    else:

        print(f"Attempting to download: {file_name}...")

        # Request
        response = requests.get(url, stream=True)
        
        # Status Code 200 = OK
        if response.status_code == 200:
            with open(output_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            print(f"Success! Saved to: {output_path}")
        else:
            print(f"Error {response.status_code}: Could not find model for ID {uniprot_id}")


def run_dssp(cif_file, output_dir):
    # Conf
    URL_PDB = "https://pdb-redo.eu/dssp/do"
    input_file_path = os.path.join(output_dir, cif_file) 
    output_file = f"{cif_file}.dssp"
    output_path = os.path.join(output_dir, output_file)                  

    if os.path.exists(output_path):
        print("DSSP file already exists")
    else:
        # Prepare the payload
        form_data = {
            'format': 'dssp'
        }

        # Open the file and send the POST request
        print(f"Sending {cif_file} to DSSP server...")

        with open(input_file_path, 'rb') as f:
            # 'files' matches the file upload (-F "data=<filename")
            # The key 'data' corresponds to the field name the server expects
            files_payload = {'data': f}

            try:
                response = requests.post(URL_PDB, data=form_data, files=files_payload)

                # Status Code 200 = OK
                if response.status_code == 200:
                    # Write the result to the output file
                    with open(output_path, 'wb') as out_f:
                        out_f.write(response.content)
                    print(f"Success! Output saved to: {output_path}")
                else:
                    print(f"Error {response.status_code}: {response.text}")
            except Exception as e:
                print(f"An error occurred: {e}")

def process_dssp(file_path):
    # make_dssp_dict parses the file into a dictionary
    # Key = (chain, res_id), Value = (aa, sec_struct, acc)
    parsed_data = []
    if os.path.exists(file_path):
        dssp_dict, keys = make_dssp_dict(file_path)

        for key in keys:
            chain_id = key[0]
            res_id   = key[1][1] # BioPython residue ID tuple

            # dssp_dict returns a tuple.
            # Index 0 = Amino Acid
            # Index 1 = Secondary Structure
            # Index 2 = ACC (Accessibility)
            data = dssp_dict[key]

            parsed_data.append({
                'chain': chain_id,
                'residue': res_id,
                'AA': data[0],
                'SS': data[1],
                'ACC': data[2]
            })
    else:
        parsed_data.append({
                'chain': np.nan,
                'residue': "",
                'AA': "",
                'SS': "",
                'ACC': ""
            })
        

    return parsed_data


def mapper_3d(data, raw):
    """
    Docstring for mapper_3d
    """

    df = pd.read_csv(data)
    aa_pos = df['HGVSp'].str.split(":", expand=True)[1].str.split("p.", regex=False, expand=True)[1].str[3:-3]
    df["aa_pos"] = aa_pos
    df = df[df["CSQ"] == "missense_variant"]

    ensp_to_uniprot = pd.read_csv("idmapping_ensp_ids.tsv", sep="\t")
    df_final = df.merge(ensp_to_uniprot, on='ensp', how='left')

    df_final.dropna(subset=['Entry'])
    df_final.rename(columns={'Entry': 'uniprot_id'}, inplace=True)
    
    """uniprot_ids = [gget.info(id)["uniprot_id"].item() for id in df["ensp"]]
    df["uniprot_id"] = uniprot_ids"""

    #[download_af_model(id, output_dir=raw) for id in df_final["uniprot_id"]]

    for id in df_final["uniprot_id"]:
        try:
            run_dssp(f'AF-{id}-F1-model_v6.cif', output_dir=raw)
        except:
            print(f"AlphaFold model for {id} not found")
        
    ss_list = []

    for i in range(len(df_final)):
        id = df_final.iloc[i]["uniprot_id"]
        var_aa_pos = int(df_final.iloc[i]["aa_pos"])


        info_3d = process_dssp(f"./threed_data/AF-{id}-F1-model_v6.cif.dssp")
        # Handling vor variations that may me mapped to non-canonical forms
        try:
            ss_list.append(pd.DataFrame(info_3d[var_aa_pos - 1], index=[0]))
        except IndexError:
            print(f'IndexError: {var_aa_pos + 1} is out of range of {id} model (range: {len(info_3d)})')
            oof_row = []
            oof_row.append({
                'chain': np.nan,
                'residue': "oor",
                'AA': "oor",
                'SS': "oor",
                'ACC': "oor"
            })
            ss_list.append(pd.DataFrame(oof_row[0], index=[0]))


    
    ss_df = pd.concat(ss_list)

    ss_df = ss_df.reset_index()
    print(ss_df)
    final_df = pd.concat([df_final, ss_df], axis=1)

    return final_df



def init_argparser():
    '''
    Initializes an argument parser
    '''

    parser = argparse.ArgumentParser(
        usage='%(prog)s',
        description='Retrieving 3D data for genebass variants', 
    )

    parser.add_argument('-mapping', '--mapping_result')
    parser.add_argument('-out_dir', '--output_directory')
    parser.add_argument('-raw_dir', '--rawdata_directory')

    return parser

def main() -> None:
    
    """
    Main function for [program's name].
    """

    parser = init_argparser()
    args = parser.parse_args()

    print("Initializing")


    threed_df = mapper_3d(f"./output/{args.mapping_result}", 
                            raw=f"./{args.rawdata_directory}")
    threed_df.to_csv(Path(f"{args.output_directory}","phenotype_mapped_ccrs_3d.csv"))


if __name__ == '__main__':
    main()