#!/usr/bin/env python3

import argparse
import re
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def clean_sequence(seq):
    return re.sub(r"[^ATCGatcgNn]", "N", seq)

def load_genome(fasta_file):
    genome = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        record.seq = Seq(clean_sequence(str(record.seq)))
        genome[record.id] = record
    return genome

def parse_gtf(gtf_file):
    tx2seq = {}
    tx2str = {}
    cds = defaultdict(lambda: defaultdict(list))
    introns = defaultdict(lambda: defaultdict(int))
    with open(gtf_file, "r") as gtf_handle:
        for line in gtf_handle:
            if re.match(r"\S+\t[^\t]+\tCDS\t\d+\t\d+\t\S+\t\S+\t\d\t.*transcript_id (\S+)", line):
                seq_id, st, en, stx, fr, tx_id = re.match(
                    r"(\S+)\t[^\t]+\tCDS\t(\d+)\t(\d+)\t\S+\t(\S+)\t(\d)\t.*transcript_id (\S+)", line).groups()
                tx_id = re.sub(r'\"(\S+)\"', r'\1', tx_id)
                tx_id = re.sub(r';', r'', tx_id)
                cds[seq_id][tx_id].append({
                    'start': int(st), 'end': int(en), 'strand': stx, 'frame': int(fr)
                })
                tx2seq[tx_id] = seq_id
                tx2str[tx_id] = stx
            elif '\tintron\t' in line:
                match = re.search(r'transcript_id "?([^";]+)"?', line)
                if match:
                    seq_id = line.split('\t')[0]
                    tx_id = match.group(1)
                    introns[seq_id][tx_id] += 1
    return cds, tx2seq, tx2str, introns

def extract_coding_sequences(cds, genome):
    codingseq = {}
    for seqid in cds:
        for tx in cds[seqid]:
            txseq = Seq("")
            txstrand = cds[seqid][tx][0]['strand']
            entries = cds[seqid][tx]
            for i in range(0, len(entries)):
                cds_line = entries[i]
                if i == 0 and txstrand == '+' and cds_line['frame'] != 0:
                    txseq += Seq('N' * (3 - cds_line['frame']))
                frag = genome[seqid].seq[cds_line['start'] - 1:cds_line['end']]
                txseq += frag
                if i == (len(entries) - 1):
                    if txstrand == '+' and len(txseq) % 3 != 0:
                        txseq += Seq('N' * (3 - len(txseq) % 3))
                    elif txstrand == '-' and cds_line['frame'] != 0:
                        txseq = Seq('N' * (3 - cds_line['frame'])) + txseq
                        if len(txseq) % 3 != 0:
                            txseq = Seq('N' * (3 - len(txseq) % 3)) + txseq
            if txstrand == '-':
                txseq = txseq.reverse_complement()
            codingseq[tx] = SeqRecord(txseq, id=tx, description="")
    return codingseq

def determine_stop_inclusion(codingseq, stop_codons, max_check=1000):
    count_stop_included = 0
    count_stop_excluded = 0
    for tx in list(codingseq.keys())[:max_check]:
        seq = str(codingseq[tx].seq).upper().rstrip("N")
        if len(seq) < 3:
            continue
        last_codon = seq[-3:]
        if last_codon in stop_codons:
            count_stop_included += 1
        else:
            count_stop_excluded += 1
    return count_stop_included >= count_stop_excluded

def check_valid_codons(codingseq, start_codons, stop_codons, introns, cds, tx2seq, stop_included=True):
    valid_tx = set()
    bad_summary = []
    for tx, rec in codingseq.items():
        seq = str(rec.seq).upper().strip("N")
        if len(seq) < 6:
            bad_summary.append((tx, "short_seq", ""))
            continue
        if len(seq) % 3 != 0:
            bad_summary.append((tx, "bad_frame", ""))
            continue
        start = seq[:3]
        stop = seq[-3:] if stop_included else seq[-6:-3]
        has_start = start in start_codons
        has_stop = stop in stop_codons

        seqid = cds[tx2seq[tx]][tx][0]['start']  # any CDS line
        n_introns = introns.get(tx2seq[tx], {}).get(tx, 0)
        n_cds = len(cds[tx2seq[tx]][tx])

        if n_introns >= n_cds:
            bad_summary.append((tx, "too_many_introns", f"{n_introns} introns, {n_cds} CDS"))
            continue

        if has_start and has_stop:
            valid_tx.add(tx)
        else:
            bad_summary.append((tx, start if not has_start else "", stop if not has_stop else ""))
    return valid_tx, bad_summary

def filter_gtf(gtf_file, output_file, valid_tx):
    with open(gtf_file) as inp, open(output_file, "w") as out:
        for line in inp:
            if line.startswith("#"):
                out.write(line)
                continue
            match = re.search(r'transcript_id "?([^";]+)"?', line)
            if match:
                tx_id = match.group(1)
                if tx_id in valid_tx:
                    out.write(line)

def write_coding_sequences(codingseq, filename="all_codingseqs.fasta"):
    with open(filename, "w") as out:
        SeqIO.write(codingseq.values(), out, "fasta")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genome", required=True, help="FASTA genome file")
    parser.add_argument("-t", "--gtf", required=True, help="Input GTF file")
    parser.add_argument("-o", "--output", required=True, help="Filtered GTF output")
    parser.add_argument("--debug", action="store_true", help="Print bad transcripts")
    args = parser.parse_args()

    start_codons = {"ATG", "CTG", "GTG", "TTG"}
    stop_codons = {"TAA", "TAG", "TGA"}

    genome = load_genome(args.genome)
    cds, tx2seq, tx2str, introns = parse_gtf(args.gtf)
    codingseq = extract_coding_sequences(cds, genome)

    if args.debug:
        write_coding_sequences(codingseq)

    stop_included = determine_stop_inclusion(codingseq, stop_codons)

    valid_tx, bad_summary = check_valid_codons(codingseq, start_codons, stop_codons, introns, cds, tx2seq, stop_included)

    filter_gtf(args.gtf, args.output, valid_tx)

    print(f"Valid transcripts: {len(valid_tx)}")
    print(f"Invalid transcripts: {len(bad_summary)}")
    if args.debug:
        for txid, issue1, issue2 in bad_summary:
            issues = []
            if issue1:
                issues.append(issue1)
            if issue2:
                issues.append(issue2)
            print(f"transcript:{txid}\t" + "\t".join(issues))

if __name__ == "__main__":
    main()
