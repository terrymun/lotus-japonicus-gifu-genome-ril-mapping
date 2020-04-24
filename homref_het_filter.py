#!/bin/env python
import os, glob, re, sys, csv
from optparse import OptionParser

# Get script path
currentDir = os.path.dirname(os.path.realpath(__file__))

# Parser
parser = OptionParser()
parser.add_option('-l', '--log', dest="logfile", help="Write output log to file", type="string")
parser.add_option('-t', '--homhet', dest="homhetThreshold", help="Upper threshold for the proportion of hom_ref and het, compared to all genotypes scanned", type=float)
parser.add_option('-m', '--missing', dest="missingThreshold", help="Upper threshold for the proportion of missing genotypes, compared to all genotypes scanned", type=float)

# Default arguments
parser.set_defaults(logfile=os.path.join(currentDir, 'log.csv'), homhetThreshold=0.75, missingThreshold=0.25)

# Parse arguments
(options, args) = parser.parse_args()

# Cap homhet ratio to 1
if options.homhetThreshold > 1:
	options.homhetThreshold = 1

# Glob files
files = args

# Log file
logfile = open(options.logfile, 'w')
_l = csv.writer(logfile, delimiter=",", lineterminator="\n")
_l.writerow(['Contig', 'Sites', 'WrongParentGenotype', 'Filtered', 'MissingGenotype', 'OK'])

# Iterate through each file
for file in files:

	# Filename
	filename = re.sub('.vcf', '', os.path.basename(file))

	# Open file for reading
	with \
	open(file, 'r') as f, \
	open(os.path.dirname(os.path.abspath(file))+'/re'+filename+'.vcf', 'w') as w:

		# CSV reader
		_f = csv.reader(f, delimiter="\t")

		# Counts
		counts = {'sites': 0, 'wrongParentGenotype': 0, 'filtered': 0, 'missing': 0, 'ok': 0}

		# Iterate throguh each line
		for l in _f:

			# Write headers to new file
			if re.match(r'^#', l[0]):

				w.write("\t".join(l)+"\n")

			# Only parse lines that not are commented out
			else:

				# Column 9 belongs to the Gifu Illumina reads
				GifuGenotype = l[9].split(':')[0]

				# Column 10 belongs to the second parent Illumina reads
				Parent2Genotype = l[10].split(':')[0]

				# Check column 11 and onwards
				genotypes = l[11:]

				# Track hom_ref, hom_var, het
				hom_ref = 0
				hom_var = 0
				het = 0
				missing = 0

				# Iterate through each genotype call
				for g in genotypes:
					_g = g.split(':')[0]

					if _g == '0/0':
						hom_ref += 1
					elif _g == '0/1':
						het += 1
					elif _g == '1/1':
						hom_var += 1
					elif _g == './.':
					  missing += 1
					else:
					  print('WARNING: Unknown genotype detected for '+l[0])

				# Do not write lines if anything other than
				# - hom_ref is called on GifuGenotype
				# - hom_alt is called on Parent2Genotype
				if GifuGenotype != '0/0' or Parent2Genotype != '1/1':
					counts['wrongParentGenotype'] += 1

				# Do not write lines if the line has too many missing genotypes
				elif missing >= len(genotypes)*options.missingThreshold:
				  counts['missing'] += 1

				# Do not write lines if the line has only hom_ref and het
				elif hom_ref + het >= len(genotypes)*options.homhetThreshold:
					counts['filtered'] += 1

				# Write lines finally
				else:
					counts['ok'] += 1
					w.write("\t".join(l)+"\n")

				# Update all counts
				counts['sites'] += 1

				if counts['sites'] % 10000 == 0:
					print('Sites processed: '+str(counts['sites']))

		# Print summary
		print('# Summary for '+file)
		print('- Total sites: '+str(counts['sites']))
		print('- Wrong Gifu genotype called: '+str(counts['wrongParentGenotype']))
		print('- Filtered (too much hom_ref / het): '+str(counts['filtered']))
		print('- Missing (too many missing calls): '+str(counts['missing']))
		print('- Accepted sites: '+str(counts['ok']))

		# Write to log file
		contig = re.search(r'_(tig\d+)\.vcf$', file)
		if contig:
			_l.writerow([contig.group(1), counts['sites'], counts['wrongParentGenotype'], counts['filtered'], counts['missing'], counts['ok']])
		else:
			print('!!! Unable to extract contig name from file handler')

