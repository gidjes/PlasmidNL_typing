import subprocess
import pandas as pd
import os
import datetime
import json


# ------------------------------
# 1. General Functions
# ------------------------------
def run_cmd(cmd, log, step_name, plasmid):
    """
    Run shell command with logging and error handling.
    """

    timestamp = datetime.datetime.now().isoformat()

    log.write(f"\n[{timestamp}] START {step_name}\n")
    log.write(f"CMD: {cmd}\n")
    log.flush()

    try:
        subprocess.run(
            cmd,
            shell=True,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )

    except subprocess.CalledProcessError:
        log.write(f"[{timestamp}] ERROR in {step_name}\n")
        log.flush()
        raise

    log.write(f"[{timestamp}] END {step_name}\n")
    log.flush()


def summarise_plasmid_data(data: str, log) -> pd.DataFrame:
    log.write("combining outputs...\n")
    # Acquire (and flatten) replicon output
    replicon_df = pd.read_csv(
        f"output/{data}/plasmidfinder/{data}_replicon.csv", sep=";"
    )
    replicons = replicon_df["plasmid"].values.tolist()
    replicon = ",".join(replicons)

    # Gather mobility data
    mobility_df = pd.read_csv(f"output/{data}/mob_suite/{data}_result.txt", sep="\t")
    mobility = mobility_df["predicted_mobility"].values[0]

    # Get cluster and coordinates
    if os.path.exists(f"output/{data}/mge_cluster/{data}_prediction.csv"):
        mge_cluster = pd.read_csv(
            filepath_or_buffer=f"output/{data}/mge_cluster/{data}_prediction.csv",
            sep=",",
        )
        cluster = str(mge_cluster["Standard_Cluster"].values[0])
        tsne1D = mge_cluster["tsne1D"].values[0]
        tsne2D = mge_cluster["tsne2D"].values[0]
    elif os.path.exists(f"output/{data}/mge_cluster/{data}_prediction.rtab"):
        mge_cluster = pd.read_csv(
            filepath_or_buffer=f"output/{data}/mge_cluster/{data}_prediction.rtab",
            sep="\t",
        )
        unitigs = mge_cluster[f"{data}"].sum()
        if unitigs == 0:
            cluster = "-"
            tsne1D = "-"
            tsne2D = "-"
        else:
            print(f"Mge-cluster failed for isolates {data}. Please repeat...")
            exit()
    else:
        print(f"Mge-cluster failed for isolates {data}. Please repeat...")
        exit()

    # Get gene data, split by type and flatten
    resistence = pd.read_csv(f"output/{data}/resistance_report.csv", sep=";")
    amr_genes = resistence["gene"].loc[resistence["type"] == "AMR"].values.tolist()
    amr = ",".join(amr_genes)

    metal_genes = resistence["gene"].loc[resistence["type"] == "METAL"].values.tolist()
    metal = ",".join(metal_genes)

    virulance_genes = (
        resistence["gene"].loc[resistence["type"] == "VIRULENCE"].values.tolist()
    )
    virulence = ",".join(virulance_genes)

    biocide_genes = (
        resistence["gene"].loc[resistence["type"] == "BIOCIDE"].values.tolist()
    )
    biocide = ",".join(biocide_genes)

    heat_genes = resistence["gene"].loc[resistence["type"] == "HEAT"].values.tolist()
    heat = ",".join(heat_genes)

    acid_genes = resistence["gene"].loc[resistence["type"] == "ACID"].values.tolist()
    acid = ",".join(acid_genes)

    # Make DataFrame
    final_df = pd.DataFrame(
        {
            "Plasmid": data,
            "replicon": replicon,
            "mobility": mobility,
            "mge_cluster": cluster,
            "tnse1D": tsne1D,
            "tnse2D": tsne2D,
            "amr": amr,
            "virulence": virulence,
            "metal": metal,
            "biocide": biocide,
            "heat": heat,
            "acid": acid,
        },
        index=[0],
    )

    return final_df


def safe_final_dfs(plasmid: str) -> pd.DataFrame:
    if os.path.isfile(f"output/{plasmid}/final_report.csv"):
        final_report = pd.read_csv(f"output/{plasmid}/final_report.csv", sep=";")
    else:
        final_report = pd.DataFrame(
            {
                "Plasmid": plasmid,
                "replicon": "FAILED",
                "mobility": "FAILED",
                "mge_cluster": "FAILED",
                "tnse1D": "FAILED",
                "tnse2D": "FAILED",
                "amr": "FAILED",
                "virulence": "FAILED",
                "metal": "FAILED",
                "biocide": "FAILED",
                "heat": "FAILED",
                "acid": "FAILED",
            },
            index=[0],
        )
    return final_report


# ------------------------------
# 2. PlasmidFinder Functions
# ------------------------------
def run_plasmidfinder(plasmid, in_dir, log):
    outdir = os.path.join("output", plasmid, "plasmidfinder")
    result = os.path.join(outdir, f"{plasmid}_replicon.csv")
    fasta = os.path.join(in_dir, f"{plasmid}.fasta")

    if os.path.exists(result):
        log.write("[plasmidfinder] already completed\n")
        return

    if not os.path.exists(fasta):
        raise FileNotFoundError(f"Missing FASTA: {fasta}")

    os.makedirs(outdir, exist_ok=True)

    cmd = f"plasmidfinder.py " f"-i {fasta} " f"-o {outdir}"

    run_cmd(cmd, log, "plasmidfinder", plasmid)


def get_replicon_data(data, isolate_name):
    replicon_list = []
    if data == "No hit found":
        replicon_list = [
            {
                "isolate_key": isolate_name,
                "contig_name": "-",
                "plasmid": "No hit found",
                "identity": "-",
                "template_length": "-",
                "coverage": "-",
                "position": "-",
                "note": "-",
            }
        ]
    else:
        for plasmid_info in data:
            plasmid_data = data[plasmid_info]
            contig_name = plasmid_data["contig_name"]
            plasmid_name = plasmid_data["plasmid"]
            identity = plasmid_data["identity"]
            template_length = plasmid_data["template_length"]
            coverage = plasmid_data["coverage"]
            position = plasmid_data["positions_in_contig"]
            note = plasmid_data["note"]
            # Store the extracted data in a dictionary
            replicon_list += [
                {
                    "isolate_key": isolate_name,
                    "contig_name": contig_name,
                    "plasmid": plasmid_name,
                    "identity": identity,
                    "template_length": template_length,
                    "coverage": coverage,
                    "position": position,
                    "note": note,
                }
            ]

    return replicon_list


def process_group(group):
    group = group.sort_values(by="start").reset_index(drop=True)
    indices_to_drop = set()

    for i in range(1, len(group)):
        if group.loc[i, "start"] <= group.loc[i - 1, "end"]:
            if group.loc[i, "identity"] > group.loc[i - 1, "identity"]:
                indices_to_drop.add(group.index[i - 1])
            else:
                indices_to_drop.add(group.index[i])

    return group.drop(indices_to_drop)


def remove_overlapping(replicon_df: pd.DataFrame):
    # Split the 'position' column into 'start' and 'end' columns
    replicon_df[["start", "end"]] = [
        [int(x.split("..")[0]), int(x.split("..")[1])]
        for x in replicon_df["position"].values
    ]
    # Sort the DataFrame by the 'start' column
    replicon_df = replicon_df.sort_values(by="start").reset_index(drop=True)

    # Apply the process_group function to each group
    df_filtered = (
        replicon_df.groupby("contig_name", group_keys=False, sort=False)
        .apply(process_group, include_groups=False)
        .reset_index(drop=True)
    )
    return df_filtered


def parse_output(isolate_name):
    # List of input data.json files
    file_name = f"output/{isolate_name}/plasmidfinder/data.json"

    with open(file_name, "r") as file:
        data = json.load(file)
        # Extract information from the JSON structure
        # filename = data["plasmidfinder"]["user_input"]["filename(s)"][0]

        replicons = []
        for classification in data["plasmidfinder"]["results"].keys():
            for group in data["plasmidfinder"]["results"][classification].keys():
                replicons += get_replicon_data(
                    data["plasmidfinder"]["results"][classification][group],
                    isolate_name,
                )

        # Remove all No hits
        final_replicon_dict = [x for x in replicons if x["plasmid"] != "No hit found"]

        # Repopulate if empty
        if len(final_replicon_dict) == 0:
            final_replicon_dict.append(
                {
                    "isolate_key": isolate_name,
                    "contig_name": "-",
                    "plasmid": "No hit found",
                    "identity": "-",
                    "template_length": "-",
                    "coverage": "-",
                    "position": "-",
                    "note": "-",
                }
            )

        replicon_df = pd.DataFrame(final_replicon_dict)
        if len(replicon_df["plasmid"].values) >= 2:
            replicon_df = remove_overlapping(replicon_df)
        # Create an output file for this input file
        output_file_name = (
            f"output/{isolate_name}/plasmidfinder/{isolate_name}_replicon.csv"
        )
        replicon_df.to_csv(output_file_name, sep=";", index=False)


# ------------------------------
# 3. MOB-suite Function
# ------------------------------
def run_mob_suite(plasmid, in_dir, log):

    out_file = os.path.join("output", plasmid, "mob_suite", f"{plasmid}_result.txt")
    fasta = os.path.join(in_dir, f"{plasmid}.fasta")

    if os.path.exists(out_file):
        log.write("[mob_suite] already completed\n")
        return

    if not os.path.exists(fasta):
        raise FileNotFoundError(f"Missing FASTA: {fasta}")

    os.makedirs(f"output/{plasmid}/mob_suite", exist_ok=True)

    cmd = (
        f"conda run -n mob-suite mob_typer "
        f"--infile {fasta} "
        f"--out_file {out_file}"
    )

    run_cmd(cmd, log, "mob_suite", plasmid)


# ------------------------------
# 4. mge-cluster Function
# ------------------------------
def run_mge_cluster(plasmid, in_dir, log, custom_scheme):

    outdir = os.path.join("output", plasmid, "mge_cluster")

    rtab = os.path.join(outdir, f"{plasmid}_prediction.rtab")
    csv = os.path.join(outdir, f"{plasmid}_prediction.csv")

    if os.path.exists(rtab) and os.path.exists(csv):
        log.write("[mge_cluster] already completed\n")
        return

    if os.path.exists(rtab):
        rtab_df = pd.read_csv(rtab, sep="\t")
        if rtab_df[plasmid].sum() == 0:
            log.write("[mge_cluster] already completed (empty)\n")
            return

    if custom_scheme:
        use_custom = "1"
    else:
        use_custom = "0"

    os.makedirs(outdir, exist_ok=True)

    cmd = (
        f"conda run -n mge-cluster "
        f"./scripts/run_mge.sh {plasmid} '{in_dir}' {use_custom}"
    )

    run_cmd(cmd, log, "mge_cluster", plasmid)


# ------------------------------
# 4. Length and GC Function
# ------------------------------
def get_seqlength_and_gc(sample, log, in_dir="input"):
    """
    Calculate total sequence length and GC content across all contigs
    in a FASTA file.

    Parameters
    ----------
    row : pandas.Series
        A row from a DataFrame. Should contain a column 'Plasmid'.
    fasta_dir : str
        Directory containing FASTA files.

    Returns
    -------
    tuple
        (total_length, gc_content_percent)
    """

    log.write("Getting length + GC%\n")
    fasta_path = os.path.join(in_dir, f"{sample}.fasta")

    total_length = 0
    gc_count = 0

    with open(fasta_path) as fasta:
        for line in fasta:
            if not line.startswith(">"):
                seq_line = line.strip().upper()

                total_length += len(seq_line)
                gc_count += seq_line.count("G") + seq_line.count("C")

    if total_length == 0:
        gc_content = 0.0
    else:
        gc_content = (gc_count / total_length) * 100

    return total_length, gc_content


# ------------------------------
# 5. Final output file cleaning
# ------------------------------
def determine_upper_group(
    df_in: pd.DataFrame, col_items: str, col_groups: str
) -> pd.DataFrame:
    df = df_in.copy()

    # Create lookup dictionary
    if col_items == "replicon":
        family_df = pd.read_csv("data/replicon_classifications.csv", sep=";")
        items = col_items
        classification = col_groups
    elif col_items in ["amr", "metal"]:
        family_df = pd.read_csv("data/gene_db.csv", sep=";")
        items = "Gene"
        classification = "Class"
    else:
        print("invalid column detected, skipping step...")
        return df

    group_map = dict(zip(family_df[items], family_df[classification]))

    # Match the families
    def map_groups(value):
        # Handle NaN / None / empty strings
        if pd.isna(value) or not str(value).strip():
            return ""

        groups = {
            group_map[v.strip()]
            for v in str(value).split(",")
            if v.strip() in group_map
        }

        return ",".join(sorted(groups))

    df[col_groups] = df[col_items].apply(map_groups)

    return df


def correct_names(df_in: pd.DataFrame) -> pd.DataFrame:
    df = df_in.copy()

    # Create lookup dictionary
    family_df = pd.read_csv("data/gene_db.csv", sep=";")
    correction_map = dict(zip(family_df["Gene"], family_df["Gene_corrected"]))

    def correct_values(value):
        # Handle empty values
        if pd.isna(value) or not str(value).strip():
            return ""

        corrected = []

        for v in str(value).split(","):
            v = v.strip()

            # Replace if mapping exists
            new_value = correction_map.get(v)

            # Skip missing mappings if desired
            if new_value is not None and str(new_value).strip():
                corrected.append(new_value)

        # Remove duplicates while preserving order
        corrected = list(dict.fromkeys(corrected))
        corrected.sort()

        return ",".join(corrected)

    df["amr"] = df["amr"].apply(correct_values)
    return df


if __name__ == "__main__":
    print("Helper functions to type plasmids")
