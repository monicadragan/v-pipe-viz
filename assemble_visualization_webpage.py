import numpy as np
import os
import pandas as pd
import re
import sys
import json
from pathlib import Path

import vcf
from BCBio import GFF
from Bio import SeqIO


def convert_vcf(fname):
    """Convert VCF to JSON."""
    output = []

    with open(fname) as fd:
        vcf_reader = vcf.Reader(fd)

        for record in vcf_reader:
            output.append({
                'position': record.POS,
                'variant': [v.sequence for v in record.ALT],
                'pvalue': record.INFO['Pval']
            })

    return json.dumps(output)


def parse_gff(fname):
    """Convert GFF to map."""
    features = []

    with open(fname) as fd:
        for record in GFF.parse(fd):
            for feature in record.features:
                features.append({
                    'id': record.id,
                    'type': feature.type,
                    'name': feature.id,
                    'start': int(feature.location.start),
                    'end': int(feature.location.end)
                })
    return features

def arrange_gff_data(features):
  """Add row number to each feature.""" 
  features.sort(key=lambda f: f['start'])

  rows = []
  for feature in features:
    if not rows:
      feature["row_cnt"] = 0
      rows.append([feature])
    else:
      found = False
      for idx, row in enumerate(rows):
        if row[-1]["end"] <= feature["start"]:
          feature["row_cnt"] = idx
          row.append(feature)
          found = True
          break
      if not found:
        feature["row_cnt"] = len(rows)
        rows.append([feature])

  return [item for row in rows for item in row]

def get_gff_data(gff_dir):
  """Returns a map with filename key and gff json data."""
  gff_map = {}
  for path in os.listdir(gff_dir):
    full_path = os.path.join(gff_dir, path)
    filename = os.path.splitext(path)[0]
    gff_map[filename] = arrange_gff_data(parse_gff(full_path))
  return gff_map

def merge_coverage(csv_files, genome_length):
  """ Parses the coverage CSV files and returns a coverage array."""
  coverage_data = []
  for file in csv_files:
    csv = pd.read_csv(file, sep='\t', index_col=0, header=None) 
    coverage_data.extend([{'start':row[1], "end":row[2], "coverage":row[3]} for row in csv.values])

  coverage = [0] * genome_length
  for data in coverage_data:
    for idx in range(data["start"], data["end"] + 1):
      coverage[idx] += data["coverage"]
  return coverage

def main(sample_name, coverage_files, vcf_file, gff_dir, consensus_fasta, html_file_in, html_file_out):

    consensus = next(SeqIO.parse(consensus_fasta, "fasta")).seq.upper()
    vcf_json = convert_vcf(vcf_file)
    gff_map = get_gff_data(gff_dir)
    coverage = merge_coverage(coverage_files, len(consensus))
 
    embed_code = f"""
        var sample_name = \"{sample_name}\"
        var consensus = \"{consensus}\"
        var coverage = {coverage}
        var vcfData = {vcf_json}
        var gffData = {gff_map}
    """
    
    # assemble webpage
    with open(html_file_in) as fd:
        raw_html = fd.read()

    # TODO: make this more robust
    mod_html = raw_html.replace('{EXTERNAL_SNAKEMAKE_CODE_MARKER}', embed_code)

    with open(html_file_out, 'w') as fd:
        fd.write(mod_html)

_CONSENSUS = "/Users/mdragan/ivan_vpipe_branch//V-pipe/samples/SRR10903401/20200102/references/ref_majority.fasta"
_VCF_FILE = "/Users/mdragan/ivan_vpipe_branch//V-pipe/samples/SRR10903402/20200102/variants/SNVs/REGION_1/snv/SNVs_0.010000_final.vcf"
_GFF_DIRECTORY = "/Users/mdragan/Downloads/gffs/"
_SNV_DIR = "/Users/mdragan/ivan_vpipe_branch/V-pipe/samples/SRR10903401/20200102/variants/SNVs"

coverage_files = [os.path.join(_SNV_DIR, f, "coverage.txt") for f in os.listdir(_SNV_DIR) if re.match("REGION_[1-9]+[0-9]*", f)]

main("Sample", coverage_files, _VCF_FILE, _GFF_DIRECTORY, _CONSENSUS, "coverage.html", "out.html")


'''
if __name__ == '__main__':
    main(
        sys.argv[1], sys.argv[2],
        sys.argv[3], sys.argv[4]
    )
'''
