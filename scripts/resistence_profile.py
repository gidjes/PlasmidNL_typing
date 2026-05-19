import os
import pandas as pd
import time
from . import typing_functions as tf


def safe_read_csv(filepath, sep=";", retries=4, delay=60):
    for i in range(retries):
        try:
            return pd.read_csv(filepath, sep=sep)
        except Exception as e:
            print(f"Retry {i+1}/{retries} reading {filepath} — {e}")
            time.sleep(delay)
    raise FileNotFoundError(f"Failed to read {filepath} after {retries} attempts.")


def set_up_amrfinder(plasmid, in_dir, log):
    outdir = os.path.join("output", plasmid, "amrfinder")
    result = os.path.join(outdir, f"{plasmid}_amrfinder.txt")
    fasta = os.path.join(in_dir, f"{plasmid}.fasta")
    os.makedirs(outdir, exist_ok=True)

    if os.path.exists(result):
        log.write("[amrfinder] already completed\n")
        return

    if not os.path.exists(fasta):
        raise FileNotFoundError(f"Missing FASTA: {fasta}")

    ncbi_cmd = f"amrfinder -n {fasta} -o {result} --plus"

    tf.run_cmd(ncbi_cmd, log, "amrfinder", plasmid)


def process_amrfinder(plasmid, ignore_amr=False):
    data = safe_read_csv(
        f"output/{plasmid}/amrfinder/{plasmid}_amrfinder.txt", sep="\t"
    )
    if data.empty:
        amrfinder_df = pd.DataFrame(
            columns=[
                "gene",
                "Identity",
                "Coverage",
                "Contig",
                "Start",
                "Stop",
                "seq",
                "db",
                "type",
            ]
        )
        return amrfinder_df
    else:
        amrfinder_df = data[
            [
                "Element symbol",
                "Subtype",
                "% Coverage of reference",
                "% Identity to reference",
                "Contig id",
                "Start",
                "Stop",
            ]
        ]
        rename_dict = {
            "% Identity to reference": "Identity",
            "% Coverage of reference": "Coverage",
            "Subtype": "type",
            "Element symbol": "gene",
            "Contig id": "Contig",
        }
        amrfinder_df = amrfinder_df.rename(columns=rename_dict)

        if ignore_amr:
            amrfinder_df.loc[amrfinder_df["type"] != "AMR"]

        amrfinder_df.loc[:, "seq"] = plasmid
        amrfinder_df.loc[:, "db"] = "amrfinder"
        amrfinder_df = amrfinder_df.loc[
            (amrfinder_df["Identity"] >= 95) & (amrfinder_df["Coverage"] >= 100)
        ]
        amrfinder_df.sort_values(
            by=["Coverage", "Identity"], ascending=[False, False], inplace=True
        )
        amrfinder_df.drop_duplicates(subset=["gene"], inplace=True, keep="first")

        return amrfinder_df


def set_up_card(plasmid, in_dir, log):
    outdir = os.path.join("output", plasmid, "card")
    result = os.path.join(outdir, f"{plasmid}_card.txt")
    fasta = os.path.join(in_dir, f"{plasmid}.fasta")
    os.makedirs(outdir, exist_ok=True)

    if os.path.exists(result):
        log.write("[card] already completed\n")
        return

    if not os.path.exists(fasta):
        raise FileNotFoundError(f"Missing FASTA: {fasta}")

    card_cmd = f"conda run -n CARD rgi main -i {fasta} -o {outdir}/{plasmid}_card -d plasmid --clean"

    tf.run_cmd(card_cmd, log, "card", plasmid)


def process_card(plasmid):
    data = safe_read_csv(f"output/{plasmid}/card/{plasmid}_card.txt", sep="\t")
    if data.empty:
        card_df = pd.DataFrame(
            columns=[
                "gene",
                "Identity",
                "Coverage",
                "Contig",
                "Start",
                "Stop",
                "seq",
                "db",
                "type",
            ]
        )
        return card_df
    else:
        card_df = data[
            [
                "Best_Hit_ARO",
                "Best_Identities",
                "Percentage Length of Reference Sequence",
                "Contig",
                "Start",
                "Stop",
            ]
        ]
        rename_dict = {
            "Best_Identities": "Identity",
            "Percentage Length of Reference Sequence": "Coverage",
            "Best_Hit_ARO": "gene",
        }
        card_df.rename(columns=rename_dict, inplace=True)
        card_df["seq"] = plasmid
        card_df["Contig"] = ["_".join(x.split("_")[:-1]) for x in card_df["Contig"]]
        card_df["db"] = "CARD"
        card_df["type"] = "AMR"
        card_df = card_df.loc[
            (card_df["Identity"] >= 95) & (card_df["Coverage"] >= 100)
        ]
        card_df.sort_values(
            by=["Coverage", "Identity"], ascending=[False, False], inplace=True
        )
        card_df.drop_duplicates(subset=["gene"], inplace=True, keep="first")
        return card_df


def set_up_resfinder(plasmid, in_dir, log, db_dir: str = "data/databases/"):
    outdir = os.path.join("output", plasmid, "resfinder")
    result = os.path.join(outdir, f"ResFinder_results.txt")
    fasta = os.path.join(in_dir, f"{plasmid}.fasta")

    if os.path.exists(result):
        log.write("[resfinder] already completed\n")
        return

    if not os.path.exists(fasta):
        raise FileNotFoundError(f"Missing FASTA: {fasta}")

    resfinder_cmd = f"python -m resfinder -ifa {fasta} -o {outdir} -db_res {db_dir}/resfinder_db -acq"
    tf.run_cmd(resfinder_cmd, log, "resfinder", plasmid)


def process_resfinder(plasmid):
    all_data = safe_read_csv(
        f"output/{plasmid}/resfinder/ResFinder_results_tab.txt", sep="\t"
    )
    if all_data.empty:
        resfinder_df = pd.DataFrame(
            columns=[
                "gene",
                "Identity",
                "Coverage",
                "Contig",
                "Start",
                "Stop",
                "seq",
                "db",
                "type",
            ]
        )
        return resfinder_df
    else:
        resfinder_df = all_data[
            ["Resistance gene", "Identity", "Coverage", "Contig", "Position in contig"]
        ]
        resfinder_df.loc[:, ["Start", "Stop"]] = resfinder_df[
            "Position in contig"
        ].str.split("..", expand=True)
        resfinder_df = resfinder_df.drop("Position in contig", axis=1)
        resfinder_df.loc[:, ["Start", "Stop"]] = resfinder_df[["Start", "Stop"]].apply(
            pd.to_numeric
        )
        rename_dict = {"Resistance gene": "gene"}
        resfinder_df.rename(columns=rename_dict, inplace=True)
        resfinder_df["seq"] = plasmid
        resfinder_df["db"] = "resfinder"
        resfinder_df["type"] = "AMR"
        resfinder_df = resfinder_df.loc[
            (resfinder_df["Identity"] >= 95) & (resfinder_df["Coverage"] >= 100)
        ]
        resfinder_df.sort_values(
            by=["Coverage", "Identity"], ascending=[False, False], inplace=True
        )
        resfinder_df.drop_duplicates(subset=["gene"], inplace=True, keep="first")
        return resfinder_df


def choose_better(group, i, j, DB_PRIORITY):
    """Return the index of the better interval between i and j."""
    row_i = group.loc[i]
    row_j = group.loc[j]

    if row_i["Identity"] > row_j["Identity"]:
        return i
    if row_i["Identity"] < row_j["Identity"]:
        return j

    if row_i["Coverage"] > row_j["Coverage"]:
        return i
    if row_i["Coverage"] < row_j["Coverage"]:
        return j

    pri_i = DB_PRIORITY.get(str(row_i["db"]).lower(), 0)
    pri_j = DB_PRIORITY.get(str(row_j["db"]).lower(), 0)

    return i if pri_i >= pri_j else j


def process_group(group):
    # --- Preprocessing ---
    group = group.sort_values(by="Start").reset_index(drop=True)
    group["Coverage"] = group["Coverage"].clip(upper=100)

    # Database priority
    DB_PRIORITY = {"resfinder": 3, "amrfinder": 2, "card": 1}

    # --- Step 1: Sort and build intervals ---
    active = []  # currently overlapping intervals
    keep = []  # indices of intervals we’ll keep

    for i, row in group.iterrows():
        # Remove inactive intervals (those that ended before this one starts)
        active = [a for a in active if group.loc[a, "Stop"] >= row["Start"]]

        # Compare new interval to currently active overlaps
        overlaps = [
            a
            for a in active
            if not (
                group.loc[a, "Stop"] < row["Start"]
                or group.loc[a, "Start"] > row["Stop"]
            )
        ]

        # If no overlaps, keep it
        if not overlaps:
            active.append(i)
            keep.append(i)
            continue

        # Otherwise, determine the "best" among all overlapping intervals
        competing = overlaps + [i]
        best = competing[0]
        for c in competing[1:]:
            best = choose_better(group, best, c, DB_PRIORITY)

        # Keep only the best one and discard the others
        active = [best]
        keep = [k for k in keep if k in active] + [best]

    # Return the filtered group
    return group.loc[sorted(set(keep))].reset_index(drop=True)


def remove_overlapping(df: pd.DataFrame):
    # Sort the DataFrame by the 'start' column
    df = df.sort_values(by="Start").reset_index(drop=True)
    # Apply the process_group function to each group
    df_filtered = (
        df.groupby("Contig", group_keys=False, sort=False)
        .apply(process_group, include_groups=False)
        .reset_index(drop=True)
    )
    return df_filtered


def combine_finder_reports(plasmid, resfinder_df, amrfinder_df, card_df=pd.DataFrame()):

    dfs = [resfinder_df, amrfinder_df, card_df]
    dfs = [df for df in dfs if not df.empty]
    combined_report = pd.concat(dfs, ignore_index=True)
    combined_report_filtered = remove_overlapping(combined_report)
    combined_report_filtered.to_csv(
        f"output/{plasmid}/resistance_report.csv",
        index=False,
        sep=";",
    )


def get_carbapenemases(df: pd.DataFrame, gene_col: str, log) -> list[str]:
    """
    Extract carbapenemase genes from a given DataFrame column.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing gene information.
    gene_col : str
        Name of the column in `df` that contains gene data.

    Returns
    -------
    List[str]
        List of comma-separated carbapenemase alleles for each row.
    """
    log.write("Highlighting Carbapenemases\n")
    # Load reference carbapenemase alleles
    carba_alleles_df = pd.read_csv("data/bldb_carbapenemases.csv", sep=";")
    carba_alleles = set(
        carba_alleles_df["NAME"].dropna().tolist()
        + carba_alleles_df["blaNAME"].dropna().tolist()
    )

    # Ensure gene column is safe to process
    genes = df[gene_col].fillna("-").astype(str).str.split(",")

    # Extract matching carbapenemase alleles
    return [
        ",".join([g for g in gene_list if g in carba_alleles]) for gene_list in genes
    ]
