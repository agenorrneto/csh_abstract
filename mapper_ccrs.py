import pandas as pd
import numpy as np
import argparse
from pathlib import Path

dtype_mapping = {"chrom": str, "start": int, "end": int, "gene": str, 
 "transcript": str, "exon": str, "ranges": str, "varflag": str, 
 "syn_density": float, "cpg": float, "cov_score": float, 
 "resid": float, "resid_pctile": float, "weighted_pct": float}

MANE = pd.read_csv("MANE.GRCh38.v1.4.summary.txt", sep ="\t", header=0)
ccrs_hg38 = pd.read_csv("sort_weightedresiduals-cpg-synonymous-novariant.txt", sep="\t", low_memory=False, dtype=dtype_mapping)

def mapper(assoc, mane, phe_name, ccrs_hg38=ccrs_hg38):
    '''
    Docstring for mapper
    
    :param assoc: Description
    :param mane: Description
    :param ccrs: Description
    '''
    assoc = pd.read_csv(assoc)

    assoc['ensp'] = assoc['HGVSp'].str.split(":", expand=True)[0]

    assoc_annot = pd.merge(
                assoc,
                mane[['symbol', 'Ensembl_prot', 'HGNC_ID']],
                right_on='Ensembl_prot',
                left_on= 'ensp',
                how='inner'
                )

    # Extracting coordinates from assoc df
    genomic_pos = assoc_annot['Variant ID'].str.split('-', expand=True)[1].astype(int)
    assoc_annot['genomic_pos'] = genomic_pos

    # Same for CCRs
    ccrs_hg38['start'] = ccrs_hg38['start'].astype(int)
    ccrs_hg38['end'] = ccrs_hg38['end'].astype(int)
    ccrs_hg38['c_start'] = np.minimum(ccrs_hg38['start'], ccrs_hg38['end'])
    ccrs_hg38['c_end'] = np.maximum(ccrs_hg38['start'], ccrs_hg38['end'])

    # Merging assoc variants and CCRs coordinates datasets
    merged = pd.merge(
        assoc_annot,
        ccrs_hg38[['chrom','gene', 'transcript', 'exon', 'c_start', 'c_end', "weighted_pct", "ranges", "varflag"]],
        left_on='symbol', # Validation question: are these the best keys for merging?
        right_on='gene',
        how='inner'
    )

    # Checking for intersection
    subset_mask = (merged['c_start'] <= merged['genomic_pos']) & (merged['genomic_pos'] <= merged['c_end'])

    # Assigning status
    merged['status'] = np.nan
    merged.loc[subset_mask, 'status'] = 'variant within'

    # Filtering for only relevant results
    results = merged.dropna(subset=['status'])

    
    # Descriptive statistics
    significant_results = results[(results['P-Value'] < 8e-9) & (results["varflag"] == "VARFALSE")]
    significant_results["trait"] = phe_name

    print(f"-------{phe_name}------")
    print(f"Variants: {len(significant_results)}")
    print(f"CCRspct mean {significant_results.weighted_pct.mean()}")
    print(f"CCRspct standart deviation {significant_results.weighted_pct.std()}")
    print(f"CCRspct min {min(significant_results.weighted_pct)}")
    print(f"CCRspct min {max(significant_results.weighted_pct)}")


    return significant_results


def init_argparser():
    '''
    Initializes an argument parser
    '''

    parser = argparse.ArgumentParser(
        usage='%(prog)s',
        description='Mapping variant association into CCRs', 
    )

    parser.add_argument('-assoc', '--assoc_result')
    parser.add_argument('-out_dir', '--output_directory', nargs=1)
    parser.add_argument('-phe', '--phenotype_list')

    return parser

def main() -> None:
    
    """
    Main function for [program's name].
    """

    parser = init_argparser()
    args = parser.parse_args()

    print("Initializing")
    phe_list = pd.read_csv(args.phenotype_list)

    for phe_name, phe_path in zip(phe_list["phe_name"], phe_list["file_name"]):
        mapped_df = mapper(f"./{args.assoc_result}/{phe_path}", 
                            mane=MANE, ccrs_hg38=ccrs_hg38, phe_name=phe_name)
        mapped_df.to_csv(Path(f"{args.output_directory[0]}",f"{phe_name}_mapped_ccrs.csv"))


if __name__ == '__main__':
    main()