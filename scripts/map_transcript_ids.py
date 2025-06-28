#!/usr/bin/env python3

import argparse
import os

def extract_transcript_id(attributes):
    for attr in attributes.strip().split(";"):
        if "transcript_id" in attr:
            return attr.strip().split(" ")[-1].strip('"')
    return None

def build_gtf_dict(gtf_file):
    gtf_dict = {}
    with open(gtf_file) as fh:
        for line in fh:
            if line.strip() == "" or line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            key = "\t".join(fields[:8])  # All but attributes
            tid = extract_transcript_id(fields[8])
            if tid:
                gtf_dict[key] = tid
    return gtf_dict

def main():
    parser = argparse.ArgumentParser(description="Map to-be-masked away transcript IDs between GTF files.")
    parser.add_argument("--lst", required=True, help="List of old (bad) transcript IDs")
    parser.add_argument("--original", required=True, help="Original GTF file (bad IDs)")
    parser.add_argument("--modified", required=True, help="Modified GTF file (new IDs)")
    args = parser.parse_args()

    # Output file: <original_lst_name>_mapped.lst
    output_file = os.path.splitext(args.lst)[0] + "_mapped.lst"

    # Read list of bad transcript IDs
    with open(args.lst) as f:
        bad_ids = set(line.strip() for line in f if line.strip())

    # Build GTF dictionaries
    orig_dict = build_gtf_dict(args.original)
    mod_dict = build_gtf_dict(args.modified)

    # Map original bad transcript IDs to keys
    id_to_key = {v: k for k, v in orig_dict.items() if v in bad_ids}

    # Write matched new transcript IDs
    with open(output_file, "w") as out:
        for old_id, key in id_to_key.items():
            new_id = mod_dict.get(key)
            if new_id:
                out.write(new_id + "\n")

if __name__ == "__main__":
    main()
