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
        "--input", required=True, help="Input directory containing plasmid FASTA files"
    )

    parser.add_argument(
        "--card", action="store_true", help="Run CARD resistance typing"
    )

    parser.add_argument(
        "--skip_amrfinder",
        action="store_true",
        help="Ignore resistance genes from AMRFinderPlus",
    )

    parser.add_argument(
        "--jobs", type=int, default=1, help="Number of parallel workers (default: 1)"
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

    inputs = os.listdir(in_dir)
    inputs = [input.split(".fa")[0] for input in inputs]
    with Pool(processes=n_jobs) as pool:
        pool.map(
            partial(
                type_plasmid.process_plasmid,
                in_dir="input",
                ignore_amr=skip_amrfinder,
                card=run_card,
            ),
            inputs,
        )
        pool.close()
        pool.join()

    all_final_dfs = [typing_functions.safe_final_dfs(input) for input in inputs]
    final_df = pd.concat(all_final_dfs)
    final_df.to_csv("PlasmidNL_report.csv", sep=";", index=False)
    failed_df = final_df.loc[final_df["replicon"] == "FAILED"]
    if not failed_df.empty:
        failed_df.to_csv("failed_sequences.csv", sep=";", index=False)
        n_failed = len(failed_df)
        print(
            f"{n_failed} plasmids failed in the pipeline. Check log files to troubleshoot."
        )

    print(
        f"{len(final_df)} plasmids processed.\nThank you for using the PlasmidNL typing pipeline"
    )


if __name__ == "__main__":
    # Set multiprocessing parameters
    set_start_method("spawn")
    main()
