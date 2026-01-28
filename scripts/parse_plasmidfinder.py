# This script simplifies the plasmidfinder output into a more simple (human-readable) file.
import json
import pandas as pd
import sys


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
        replicon_df.groupby("contig_name", group_keys=False, as_index=False, sort=False)
        .apply(process_group)
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


if __name__ == "__main__":
    isolate_name = sys.argv[1]
    parse_output(isolate_name)
