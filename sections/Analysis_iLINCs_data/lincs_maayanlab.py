#!/usr/bin/env python

"""

"""

import maayanlab_bioinformatics
from cmapPy.pandasGEXpress.parse import parse
import pandas as pd
import numpy as np
import argparse
import scipy.stats as stats
import os

# Names of the files we need, must be in the same folder as this script
siginfo_file = "siginfo_beta.txt"
instinfo_file = "instinfo_beta.txt"
gene_metadata_file = "geneinfo_beta.txt"
gctx_file = "level3_beta_all_n3026460x12328.gctx"


def add_parser_args(parser):
    parser.add_argument(
        "-c", "--cell_name",
        required=True,
        type=str)
    parser.add_argument(
        "-p", "--perturbagen",
        required=True,
        type=str)
    parser.add_argument(
        "-t", "--timepoint",
        required=True,
        type=str)
    parser.add_argument(
        "-d", "--dose",
        required=True,
        type=str)
    parser.add_argument(
        "-a", "--alpha",
        required=False,
        type=float)
    parser.add_argument(
        "-o", "--output",
        required=False,
        type=str)


def write_tsv(data, fname):
    """
    Given a dataframe, write it to a tsv
    """
    print(f"writing to {fname}")
    data.to_csv(fname, sep="\t", index=False, encoding="utf-8")


def read_file(fname):
    """
    Reads in a text file as a tsv and into a dataframe
    """
    return pd.read_csv(fname, sep="\t", encoding="utf-8", low_memory=False)


def identify_signature_of_interest(fields, metadata):
    """
    Given a set of parameters, find the rows that match the parameters
    in the dataframe
    """
    sig_of_int = metadata[((metadata['cell_iname'] == fields.cell_name) &
                          (metadata['cmap_name'] == fields.perturbagen) &
                          (metadata['pert_itime'] == fields.timepoint) &
                          (metadata['pert_idose'] == fields.dose))]
    return sig_of_int


def main():
    """

    """
    # 2. This defines the parameters for the signature of interest
    parser = argparse.ArgumentParser()
    add_parser_args(parser)
    fields = parser.parse_args()

    # Reading in the geneinfo_beta text file as a dataframe
    gene_info = read_file(gene_metadata_file)

    # 3. Extract and store the row names, column names, and data matrix
    # from the level 3 GCTX file
    # NOTE: we don't extract and store the data matrix as it's too big.
    # We will instead read in slices of the GCTX file whenever we need
    rnames = parse(gctx_file, row_meta_only=True)
    rnames = rnames.reset_index()
    cnames = parse(gctx_file, col_meta_only=True)
    cnames = cnames.reset_index()

    # 4. Store the batch IDs for all signatures in the signature metadata file
    signature_metadata = read_file(siginfo_file)
    # Get only the signature id and store it in batch_id
    signature_metadata['batch_id'] = signature_metadata['sig_id']
    # The batch id is the text that occurs before the colon
    # so split by : and get the first element
    signature_metadata['batch_id'] = signature_metadata[
            'batch_id'].str.split(":").str[0]

    # 5. Store the batch IDs for all perturbation instances in the
    # instance metadata file
    instance_metadata = read_file(instinfo_file)
    instance_metadata['batch_id'] = instance_metadata['sample_id']
    # The batch id for the instance file is the first three underscore
    # delimited terms so we just split by _ and get the first three elements
    instance_metadata['batch_id'] = instance_metadata[
            'batch_id'].str.split("_").str[:3]
    # Column is a list of the three elements so need to join back together
    instance_metadata['batch_id'] = instance_metadata[
            'batch_id'].apply(lambda x: "_".join(x))

    # 6. Identify the signature of interest
    print("These parameters define the signature of interest")
    for f in vars(fields):
        print(f"{f}: {getattr(fields, f)}")

    # If an output folder is specified, use that
    # else we create an output folder based on the inputs
    if fields.output is None:
        t = fields.timepoint.split(" ")[0]
        d = fields.dose.split(" ")[0]
        output = f"{fields.cell_name}_{fields.perturbagen}_{t}_{d}"
    else:
        # An output folder was specified on input so save to that
        output = fields.output

    # Create the output directory if it doesn't already exist
    os.makedirs(output, exist_ok=True)
    print(f"Files will be saved to the {output}/ folder")

    sig_of_int = identify_signature_of_interest(fields, signature_metadata)
    print(f"There are {len(sig_of_int.index)} signature(s) of interest")

    print("===================================================")

    # If multiple signatures from different batches match,
    # perform the rest of this protocol for each signature individually
    for index, soi in sig_of_int.iterrows():
        print("Signature of Interest with:")
        print(f"signature_id={soi['sig_id']}")
        print(f"batch_id={soi['batch_id']}")

        sig_batch_id = soi['batch_id']

        # 7a. Use the distil_ids column to identify the instance IDs
        # corresponding to the signature of interest
        distil_ids = soi['distil_ids'].split("|")
        # This is a pipe-delimited list of IDs to need to split by |
        # to get the instance IDs separated
        # print("These are the distil_ids:")
        # print(", ".join(distil_ids))

        # 7b. Slice those instances from the level 3 GCTX file and store
        # as the treatment matrix
        treatment_idx = list(cnames[cnames['cid'].isin(distil_ids)].index)
        treatment = parse(gctx_file, cidx=treatment_idx)
        treatment = treatment.data_df
        # print("This is the treatment data matrix")
        # print(treatment)

        # 8a. Filter the instance metadata by the batch ID corresponding
        # to the signature of interest
        other_instances = instance_metadata[
                instance_metadata['batch_id'] == sig_batch_id]
        oi_size = len(other_instances.index)

        # 8b. Remove the treatment instances
        # Filter out any of the treatment instance IDs from this dataframe
        control_instances = other_instances[
                ~other_instances['sample_id'].isin(distil_ids)]

        ci_size = len(control_instances.index)
        removed = oi_size - ci_size
        print(
            f"{removed} instance(s) were removed to create the control matrix")

        # 8c. Slice the instances from the level 3 GCTX file and store
        # as the control matrix
        control_idx = list(control_instances.index)
        control = parse(gctx_file, cidx=control_idx)
        control = control.data_df
        # print("This is the control data matrix")
        # print(control)

        # 9. Run the Characteristic Direction (CD) method on the treatment
        # data and the control data
        cd = maayanlab_bioinformatics.dge.characteristic_direction(
                control, treatment)

        # Move the index to a column so it can be merged with gene_info later
        cd['gene_id'] = cd.index.astype(int)

        # 10a. Compute the 2-tailed z-score and p-values for each coefficient
        # on the CD results
        mean = cd['CD-coefficient'].mean()
        std_dev = cd['CD-coefficient'].std()
        print(f"mean={mean}, standard deviation={std_dev}")
        cd['z_score'] = (cd['CD-coefficient'] - mean) / std_dev
        cd['p_value'] = 2 * (1 - stats.norm.cdf(np.abs(cd['z_score'])))

        # 10b. Map the row IDs in step 3 to each coefficient and it's p-value
        # using the gene metadata file, the gene_id column provides the index
        # of each gene
        cd = cd.merge(gene_info, how="left", on=['gene_id'])

        # Just saving the whole file with the CD coefficients and don't
        # need the blow code to separate the up vs down regulated genes
        file_stem = f"{soi['sig_id']}_{soi['batch_id']}"
        write_tsv(cd, f"{output}/{file_stem}_characteristic_direction.tsv")

        # alpha = fields.alpha

        # # 11. Identify all genes that correspond to a CD coefficient greater
        # # than 0 and p-value < alpha
        # up_regulated = cd[((cd['CD-coefficient'] > 0) &
        #                  (cd['p_value'] < alpha))]
        # print("Up-regulated genes")
        # print(up_regulated['gene_symbol'])

        # # 12. Identify all genes that correspond to a CD coefficient less
        # # than 0 and p-value < alpha
        # down_regulated = cd[((cd['CD-coefficient'] < 0) &
        #                    (cd['p_value'] < alpha))]
        # print("Down-regulated genes")
        # print(down_regulated['gene_symbol'])

        # stem = f"{soi['sig_id']}_{soi['batch_id']}"
        # # Write these genes to tsv files
        # write_tsv(up_regulated, f"{stem}_up_regulated.tsv")
        # write_tsv(down_regulated, f"{stem}_down_regulated.tsv")

        print("===================================================")


if __name__ == '__main__':
    main()
