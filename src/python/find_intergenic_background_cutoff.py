#! /usr/local/packages/Python-2.6.4/bin/python

from sys import *
from collections import defaultdict
import optparse
import re

###############################################################################
# command line parameters

usage = """find_intergenic_background_cutoff.py [options] zcontig_length_file gff3_file wig_file*

This script produces a depth-of-coverage cut-off intended for transcript
finding, defined as a particular quantile of the distribution of
depths-of-coverage over what we hope are dependably intergenic regions.
These are positions meeting the following criteria:
  No feature covers the position
  The nearest flanking features both point away from the position
  The distance to those flanking features is neither too short nor too long

The contig length file should be tab-delimited, with no header and two columns: contig ID and length
"""

parser = optparse.OptionParser(usage=usage)
parser.add_option('-q', '--quantile', type='float', default=0.7,
    help='quantile (0-1) of coverage to output (default 0.7)')
parser.add_option('-n', '--min_interbutt', type='int', default=50,
    help='minimum distance from nearest flanking feature (default 50)')
parser.add_option('-x', '--max_interbutt', type='int', default=1000,
    help='maximum distance from nearest flanking feature (default 1000)')
parser.add_option('-g', '--gff3_file', help='Path to a GFF3 file')
parser.add_option('-c', '--contig_length_file', help='Path to a contig lengths file')
parser.add_option('-w', '--wig_file', help='Path to a WIG file')
parser.add_option('-W', '--second_wig', help='Path to a second WIG file that is a pair to the file specified in --wig_file')

(options, args) = parser.parse_args()

quantile = options.quantile
min_interbutt = options.min_interbutt
max_interbutt = options.max_interbutt

contig_length_file = open(options.contig_length_file)
gff3_file = open(options.gff3_file)

wig_files = []
wig1 = open(options.wig_file)
wig_files.append(wig1)
if options.second_wig:
    wig2 = open(options.second_wig)
    wig_files.append(wig2)

###############################################################################
# read contig length file

contig_length = {}
for line in contig_length_file:
  contig, length = line[:-1].split('\t')
  contig_length[contig] = int(length)
contig_length_file.close()

###############################################################################
# read gff3 file

contig_direction_position_genic = {}

for line in gff3_file:
  if line.startswith('#'):
    continue

  contig, unk1, span_type, start, end, unk2, direction, unk3, att_val_pairs = line[:-1].split('\t')

  start, end = map(int, (start, end))
  direction = intern(direction)

  try:
    direction_position_genic = contig_direction_position_genic[contig]
  except:
    direction_position_genic = contig_direction_position_genic[contig] = {
        '+' : [False] * contig_length[contig],
        '-' : [False] * contig_length[contig]
        }

  position_genic = direction_position_genic[direction]

  for position in range(start - 1, end):
    position_genic[position] = True
gff3_file.close()

###############################################################################
# calculate interbutts: for each position on each contig, is the position
# outside of any annotated gene, are the nearest flanking genes both oriented
# away from the current position, and if so, what is the distance to the
# nearest of the two flanking genes

def calculate_distance_from_most_recent_gene(position_genic):
  distance_from_most_recent_gene = []
  distance = 0
  for position, genic in enumerate(position_genic):
    if genic:
      distance = 0
    else:
      distance += 1
    distance_from_most_recent_gene.append(distance)
  return distance_from_most_recent_gene

contig_position_interbutt = {}
for contig, direction_position_genic in contig_direction_position_genic.iteritems():
  distance_from_plus_left = calculate_distance_from_most_recent_gene(direction_position_genic['+'])
  distance_from_minus_left = calculate_distance_from_most_recent_gene(direction_position_genic['-'])
  distance_from_plus_right = list(reversed(calculate_distance_from_most_recent_gene(reversed(direction_position_genic['+']))))
  distance_from_minus_right = list(reversed(calculate_distance_from_most_recent_gene(reversed(direction_position_genic['-']))))

  position_interbutt = contig_position_interbutt[contig] = []
  for position in range(contig_length[contig]):
    position_interbutt.append(
        distance_from_minus_left[position] < distance_from_plus_left[position]
        and distance_from_plus_right[position] < distance_from_minus_right[position] 
        and min(distance_from_minus_left[position], distance_from_plus_right)
        )

###############################################################################
# read WIG files

#One of two possible header types
variableStep_header_line_pat = re.compile('^variableStep chrom=(.*)$')
fixedStep_header_line_pat = re.compile(r'^fixedStep chrom=(\S+) start=(\d+) step=(\d+)')

def read_wig(file):
  contig_position_count = defaultdict(lambda: defaultdict(lambda: 0))
  contig = None
  in_fixed = in_variable = False
  for line in file:
    #Either match the fixed step header or the variable step header
    match = fixedStep_header_line_pat.match(line)
    if match:
      contig, start, step = match.groups()
      position = int(start)
      step = int(step)
      in_fixed = True
      in_variable = False
    else:
      match = variableStep_header_line_pat.match(line)
      if match:
        contig, = match.groups()
        in_variable = True
        in_fixed = False
      elif in_fixed:
        assert contig is not None
        count = int(line.strip())
        contig_position_count[contig][position] = count
        position += step
      else:
        assert in_variable
        assert contig is not None
        position, count = map(float, line[:-1].split('\t'))
        contig_position_count[contig][position] = count

  return contig_position_count

contig_position_counts = []
for wig_file in wig_files:
  contig_position_counts.append(read_wig(wig_file))
  wig_file.close()

###############################################################################
# for each WIG file, collect read counts for each position meeting the
# interbutt criteria

target_zone_countss = []

for contig_position_count in contig_position_counts:
  target_zone_counts = []
  for contig, position_interbutt in contig_position_interbutt.iteritems():
    position_count = contig_position_count[contig]

    for position, interbutt in enumerate(position_interbutt):
      count = position_count[position]
      if interbutt is not False and min_interbutt < interbutt < max_interbutt:
        target_zone_counts.append(count)
  target_zone_countss.append(target_zone_counts)

###############################################################################
# calculate a cutoff as the requested quantile for each set of counts, and
# average the individual cutoffs as the final output cutoff

def get_quantile(data, fraction=0.5):
  data = list(sorted(data))
  position = fraction * (len(data) - 1)

  if position % 1 == 0.0:
    return data[int(position)]
  else:
    position = int(position)
    return 0.5 * (data[position] + data[position + 1])

cutoffs = [get_quantile(counts, quantile) for counts in target_zone_countss]
print sum(cutoffs) / float(len(cutoffs))

