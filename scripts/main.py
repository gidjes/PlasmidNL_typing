#!/usr/bin/env python

import os
import pandas as pd
from functools import partial
from multiprocessing import Pool, set_start_method
import argparse

from . import type_plasmid, typing_functions


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run PlasmidNL typing pipeline (in parallel)"
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input directory containing plasmid FASTA files",
    )

    parser.add_argument(
        "-c", "--card", action="store_true", help="Run CARD resistance typing"
    )

    parser.add_argument(
        "-sa",
        "--skip_amrfinder",
        action="store_true",
        help="Ignore resistance genes from AMRFinderPlus",
    )

    parser.add_argument(
        "-mge",
        "--custom_mge_scheme",
        action="store_true",
        help="Use a custom mge-cluster scheme",
    )

    parser.add_argument(
        "-n",
        "--jobs",
        type=int,
        default=1,
        help="Number of parallel workers (default: 1)",
    )

    parser.add_argument(
        "-rf",
        "--rerun-failed",
        action="store_true",
        help="Rerun only sequences listed in failed_sequences.csv",
    )

    return parser.parse_args()


def main():
    print(
        """                                                                                                            
        PlasmidNL Typing pipeline        
                        .   ...::-=-.     
                        :====-::=-:--.    
                    ..  .======-::::::.   
                    :.. .-=========:::    
                    :==-.  ..-==::=-:=:.  
                    :==-.  :.:::::--::.   
                .=+=: :=::::===-::.       
                -=--:====::::=::--:.      
                .-::==:---====:====::     
            .:-:-==::--:-==--::=-:.       
            :::----==-=-:::=::.:..        
            ..::-:::::==::===:            
        .---:::-=::::=-:::=-::            
        :=-.:::..  .. :-::==:::           
        . :::           :=::. .:::.       
                :..::...::.  :: .:. ..    
                ..  ::.. .. :::. ::-:.    
                .......::.. .::  :..:.    
    """
    )
    args = parse_args()
    in_dir = args.input
    n_jobs = args.jobs
    run_card = args.card
    skip_amrfinder = args.skip_amrfinder
    rerun_failed = args.rerun_failed
    custom_mge_scheme = args.custom_mge_scheme

    # -------------------------------------------------
    # Determine inputs
    # -------------------------------------------------
    all_inputs = os.listdir(in_dir)
    all_inputs = [input.split(".fa")[0] for input in all_inputs]
    if rerun_failed:
        failed_file = "failed_sequences.csv"

        if not os.path.exists(failed_file):
            raise FileNotFoundError(
                f"{failed_file} not found. Cannot rerun failed plasmids."
            )

        failed_df = pd.read_csv(failed_file, sep=";")

        # assuming the sequence/sample name column is called "input"
        inputs = failed_df["Plasmid"].tolist()

        print(f"Rerunning {len(inputs)} failed plasmids")

    else:
        inputs = all_inputs

    # -------------------------------------------------
    # Run pipeline
    # -------------------------------------------------
    with Pool(processes=n_jobs) as pool:
        pool.map(
            partial(
                type_plasmid.process_plasmid,
                in_dir=in_dir,
                ignore_amr=skip_amrfinder,
                card=run_card,
                custom_scheme=custom_mge_scheme,
            ),
            inputs,
        )

    # -------------------------------------------------
    # Process results
    # -------------------------------------------------
    all_final_dfs = [typing_functions.safe_final_dfs(input) for input in all_inputs]
    final_df = pd.concat(all_final_dfs)
    final_df = typing_functions.determine_upper_group(
        final_df, "replicon", "replicon_family"
    )
    final_df = typing_functions.determine_upper_group(final_df, "amr", "amr_classes")
    final_df = typing_functions.determine_upper_group(
        final_df, "metal", "metal_classes"
    )
    final_df = typing_functions.correct_names(final_df)
    final_df.to_csv("PlasmidNL_report.csv", sep=";", index=False)
    failed_df = final_df.loc[final_df["replicon"] == "FAILED"]

    if not failed_df.empty:
        failed_df.to_csv("failed_sequences.csv", sep=";", index=False)
        n_failed = len(failed_df)
        print(
            f"{n_failed} plasmids failed in the pipeline. Check log files to troubleshoot."
        )
    elif os.path.isfile("failed_sequences.csv"):
        os.remove("failed_sequences.csv")
        print("No more failed sequences in data.")

    print(
        f"{len(final_df)} plasmids processed.\nThank you for using the PlasmidNL typing pipeline"
    )


if __name__ == "__main__":
    # Set multiprocessing parameters
    set_start_method("spawn")
    main()
