import typing_functions as tf
import resistence_profile as rp
import os
import sys


def process_plasmid(
    plasmid: str, in_dir: str, ignore_amr: bool = False, card: bool = True
):
    os.makedirs("logs", exist_ok=True)
    log_path = os.path.join("logs", f"{plasmid}.log")

    with open(log_path, "a") as log:

        log.write(f"\n===== START {plasmid} =====\n")
        log.flush()

        try:
            # First run all the tools
            tf.run_plasmidfinder(plasmid, in_dir, log)
            tf.run_mob_suite(plasmid, in_dir, log)
            tf.run_mge_cluster(plasmid, in_dir, log)
            rp.set_up_resfinder(plasmid, in_dir, log)
            rp.set_up_amrfinder(plasmid, in_dir, log)
            if card:
                rp.set_up_card(plasmid, in_dir, log)

            # Clean up output
            tf.parse_output(plasmid)
            res_df = rp.process_resfinder(plasmid)
            amr_df = rp.process_amrfinder(plasmid, ignore_amr)
            if card:
                card_df = rp.process_card(plasmid)
                rp.combine_finder_reports(plasmid, res_df, amr_df, card_df)
            else:
                rp.combine_finder_reports(plasmid, res_df, amr_df)

            # Combine outputs
            final_df = tf.summarise_plasmid_data(plasmid, log)

            # Add last pieces
            final_df[["bp_length", "GC%"]] = tf.get_seqlength_and_gc(
                plasmid, log, in_dir
            )
            final_df["carba_alleles"] = rp.get_carbapenemases(final_df, "amr", log)

            final_df.to_csv(f"output/{plasmid}/final_report.csv", sep=";", index=False)

        except Exception as e:
            log.write(f"\nFATAL ERROR: {e}\n")
            log.write("PIPELINE FAILED FOR THIS SAMPLE\n")
            log.flush()

            return False

        log.write(f"===== END {plasmid} =====\n")
        log.flush()

    return True


def main():
    sample = sys.argv[1]
    in_dir = sys.argv[2]
    print(process_plasmid(sample, in_dir, True))


if __name__ == "__main__":
    main()
