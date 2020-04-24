#!/usr/env python
import os, csv, re, sys, collections, argparse, random, math, numpy, itertools

# Current directory
currentDir = os.path.dirname(os.path.abspath(__file__))

# Argument parser
parser = argparse.ArgumentParser(description='Smoothening and tranposing genotype call data produced')
parser.add_argument('-i', '--input',
	dest='inputFile', 
	action='store',
	default=False,
	required=True,
	help='Relative or absolute path to the input file (merged.012). This will be internally resolved to an absolute path.')
parser.add_argument('-c', '--contigLengths',
	dest='contigLengthFile',
	action='store',
	default=False,
	required=True,
	help='The file path to a tab separated file containing contig IDs in the first column and lengths in the next column.')
parser.add_argument('--passLimit',
	dest='passLimit',
	action='store',
	default=10,
	type=int,
	help='Maximum number of passes when smoothing out genotyping calls.')
parser.add_argument('--cherryPick',
	dest='cherryPick',
	action='store',
	default=1,
	type=float,
	help='Cherry pick genotype call positions. Value should be between 0 and 1 (default), which represents the proportion of random call positions to pick after setting aside n number of calls for the start, middle, and end of the contig (n is 3x the value of the --window option). Value of 0 means only a maximum of 12 calls are chosen. Value of 1 means all calls are chosen.')
parser.add_argument('--popType',
	dest='populationType',
	action='store',
	default=False,
	required=True,
	help='The type of mapping population used. Possibe values are DH and RIL<n>, where n is the number of generations.')
parser.add_argument('--popName',
	dest='populationName',
	action='store',
	default=False,
	required=True,
	help='Population name. Only characters within the set [a-zA-z0-9] are accepted.')
parser.add_argument('--distFun',
	dest='distanceFunction',
	action='store',
	default='kosambi',
	help='The distance function used. Possible values: kosambi and haldane. Default: kosambi')
parser.add_argument('--pvalue',
	dest='pvalue',
	action='store',
	default=0.000001,
	type=float,
	help='The threshold to be used for clustering the markers into LGs. A reasonable choice of p_value is 0.000001 (default). Alternatively, the user can turn off this feature by setting this argument to any number larger than 1. If the user does so, MSTMap assumes that all markers belong to one single linakge group.')
parser.add_argument('--noMapDist',
	dest='noMapDist',
	action='store',
	default=15,
	type=int,
	help='Defines an isolated marker group that is this number of centimorgans away from the rest of the markers. Default: 15.')
parser.add_argument('--noMapSize',
	dest='noMapSize',
	action='store',
	default=1,
	type=int,
	help='Defines an isolated marker group that belongs to a set of this size or smaller. A reasonable choice of this parameter is 1 (default) or 2. Disable bad marker detection alltogether by setting this to 0.')
parser.add_argument('--missingThreshold',
	dest='missingThreshold',
	action='store',
	default=1.00,
	type=float,
	help='Eliminate markers that have more than this proportion of missing observation.')
parser.add_argument('--estimation',
	dest='estimation',
	action='store',
	default='yes',
	help='A binary flag which can be set to yes (default) or no. If set to yes, then MSTMap will try to estimate missing data before clustering the markers into linkage groups.')
parser.add_argument('--detectBadData',
	dest='detectBadData',
	action='store',
	default='yes',
	help='A binary flag which can be set to yes (default) or no. If set to "yes", then MSTMap will try to detect bad data during the map construction process. Those suspicious genotype data will be printed to the console for user inspection. The error detection feature can be turned off by setting this option to "no".')
parser.add_argument('--objFun',
	dest='objectiveFunction',
	action='store',
	default='COUNT',
	help='The objective function to be used. Possible choices are "COUNT" (default) and "ML". COUNT refers to the commonly used sum of recombination events objective function and ML refers to the commonly used maximum likelihood objective function.')
parser.add_argument('--skipEnds',
	dest='skipEnds',
	action='store',
	default=5000,
	type=int,
	help='Number of bases at the start and end of a contig where genotype calls will be ignored')
parser.add_argument('--minSitesPerContig',
	dest='minSitesPerContig',
	action='store',
	default=5,
	type=int,
	help='Minimum number of sites per contig before random sampling (based on --cherryPick parameter) kicks in')
parser.add_argument('-o', '--out',
	dest='out',
	action='store',
	default=None,
	type=str,
	help='Option alternate output folder. By default output is written to the same folder where the input file is located')

args = parser.parse_args()

# Parse argument
mergedFile = os.path.abspath(args.inputFile)
mergedDir = os.path.dirname(mergedFile)

if args.out is not None:
	outDir = os.path.abspath(args.out)
	if not os.path.exists(outDir):
		os.makedirs(outDir)
else:
	outDir = mergedDir


# Cache
contigConsensusCache = []


# Read files
with \
open(mergedFile, 'r') as callsFile, \
open(os.path.join(mergedDir, 'refiltered.merged.012.pos'), 'r') as posFile, \
open(os.path.join(outDir, 'refiltered.merged.012.smoothed'), 'w') as outFile:

	# Get the number of positions per contig
	contigs = collections.OrderedDict()
	posFileReader = csv.reader(posFile, delimiter="\t")

	# Iterate through all lines in the position file
	for p in posFileReader:
		
		# Assign variables
		contig = p[0]
		position = p[1]

		# Store in dict
		if contig not in contigs:
			contigs[contig] = [p[1]]
		else:
			contigs[contig].append(p[1])

	# Collect number of positions in each contig into a list
	contigPosCounts = [0]
	for contig, pos in contigs.items():
		contigPosCounts.append(len(pos))

	print('Total number of positions detected: '+str(sum(contigPosCounts)))

	# Read each line in the calls file
	callsFileReader = csv.reader(callsFile, delimiter="\t")
	callsFileWriter = csv.writer(outFile, delimiter="\t", lineterminator="\n")
	countRIL = 0
	countCall = 0

	# Iterate through each line (RIL individual)
	# We skip the first two lines because they belong to parental lines and should have uniform genotype calls
	p1 = next(callsFileReader)
	p2 = next(callsFileReader)
	rep = ['A', 'X', 'B', '-']
	callsFileWriter.writerow([rep[int(x)] for x in p1])
	callsFileWriter.writerow([rep[int(x)] for x in p2])
	countRIL += 2
	countCall += len(p1) + len(p2) - 2

	for c in callsFileReader:

		# Store new calls
		_c = [c[0]]
		__c = c[1:]

		# Split row into a list of genotype calls per contig
		# Skip the first column because it's the index of the RIL individual
		for s in range(1, len(contigPosCounts)):

			_thisPos = sum(contigPosCounts[0:s])
			_nextPos = _thisPos + contigPosCounts[s]
			#_thisPos = contigPosCounts[s]
			#_nextPos = _thisPos + contigPosCounts[s+1]
			
			# Genotype calls originating from a specific contig
			_calls = __c[_thisPos:_nextPos]

			# Now we check for contiguous heterozygous calls
			# Do this iteratively with a limit
			passCount = 0
			callString = ''.join([rep[int(x)] for x in _calls])

			while re.search(r"([AB])([X\-]+)\1", callString):

				passCount += 1

				callString = re.sub(r"([AB])([X\-]+)\1",lambda m: "{0}{1}{0}".format(m.group(1),m.group(1)*len(m.group(2))),callString)

				# Break out of while loop when recursion exceeds 10
				if passCount > args.passLimit:
					break

			_c += list(callString)

		# Write to smoothened output
		callsFileWriter.writerow(_c)

		# Track counts
		countRIL += 1
		countCall += len(c) - 1

		# Early break
		#if countRIL > 10:
		#	break

		print('Processing RIL #'+str(countRIL))

	# Summarise
	print("RILs processed: "+str(countRIL))
	print("Genotype calls processed: "+str(countCall))

callsFile.close()
posFile.close()
outFile.close()

# Open files for tranposing
with \
open(os.path.join(outDir, 'refiltered.merged.012.smoothed'), 'r') as inFile, \
open(os.path.join(outDir, 'refiltered.merged.012.transposed'), 'w') as transposedFile, \
open(os.path.join(outDir, 'refiltered.merged.012.tmp'), 'w') as tempFile, \
open(os.path.join(mergedDir, 'refiltered.merged.012.pos'), 'r') as posFile, \
open(os.path.join(mergedDir, 'refiltered.merged.012.indv'), 'r') as indvFile, \
open(os.path.join(outDir, 'contigConsensus.csv'), 'w') as contigConsensusFile, \
open(os.path.join(outDir, 'contigConsensus.tmp'), 'w') as contigConsensusTempFile, \
open(os.path.join(outDir, 'contigConsensus.txt'), 'w') as contigConsensusOutFile, \
open(os.path.join(outDir, 'tig00000000.txt'), 'w') as contigSampleFile, \
open(args.contigLengthFile, 'r') as contigLengthFile:

	# File handlers
	transposedFileWriter = csv.writer(transposedFile, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE, escapechar=None)
	tempFileWriter = csv.writer(tempFile, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE, escapechar=None)
	posFileReader = list(csv.reader(posFile, delimiter="\t"))
	indvFileReader = list(csv.reader(indvFile, delimiter="\t"))
	contigConsensusWriter = csv.writer(contigConsensusFile, delimiter="\t", lineterminator="\n")
	contigConsensusTempWriter = csv.writer(contigConsensusTempFile, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE, escapechar=None)
	contigConsensusOutWriter = csv.writer(contigConsensusOutFile, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE, escapechar=None)
	contigLengthFileReader = list(csv.reader(contigLengthFile, delimiter="\t"))

	# Write header
	tempFileWriter.writerow(['locus_name'] + [re.sub('-', '', x[0]) for x in indvFileReader])
	contigConsensusTempWriter.writerow(['locus_name'] + [re.sub('-', '', x[0]) for x in indvFileReader])
	contigConsensusWriter.writerow(['Contig'] + [re.sub('-', '', x[0]) for x in indvFileReader[2:]] * 2 + ['AverageConfidenceScore', 'ConsensusRate'])

	# Write metadata to transposed file
	transposedFileWriter.writerow(['population_type '+args.populationType])
	transposedFileWriter.writerow(['population_name '+args.populationName])
	transposedFileWriter.writerow(['distance_function '+args.distanceFunction])
	transposedFileWriter.writerow(['cut_off_p_value '+str('{:.16f}'.format(args.pvalue))])
	transposedFileWriter.writerow(['no_map_dist '+str(args.noMapDist)])
	transposedFileWriter.writerow(['no_map_size '+str(args.noMapSize)])
	transposedFileWriter.writerow(['missing_threshold '+str(args.missingThreshold)])
	transposedFileWriter.writerow(['estimation_before_clustering '+str(args.estimation)])
	transposedFileWriter.writerow(['detect_bad_data '+str(args.detectBadData)])
	transposedFileWriter.writerow(['objective_function '+str(args.objectiveFunction)])

	# Write metadata to contig out file
	contigConsensusOutWriter.writerow(['population_type '+args.populationType])
	contigConsensusOutWriter.writerow(['population_name '+args.populationName])
	contigConsensusOutWriter.writerow(['distance_function '+args.distanceFunction])
	contigConsensusOutWriter.writerow(['cut_off_p_value '+str('{:.16f}'.format(args.pvalue))])
	contigConsensusOutWriter.writerow(['no_map_dist '+str(args.noMapDist)])
	contigConsensusOutWriter.writerow(['no_map_size '+str(args.noMapSize)])
	contigConsensusOutWriter.writerow(['missing_threshold '+str(args.missingThreshold)])
	contigConsensusOutWriter.writerow(['estimation_before_clustering '+str(args.estimation)])
	contigConsensusOutWriter.writerow(['detect_bad_data '+str(args.detectBadData)])
	contigConsensusOutWriter.writerow(['objective_function '+str(args.objectiveFunction)])

	# Store contig lengths as dict
	contigLengths = {}
	for contig in contigLengthFileReader:
		if contig:
			contigLengths[contig[0]] = int(contig[1])

	# Store content in memory
	content = []

	# Tranpose in memory
	a = zip(*csv.reader(inFile, delimiter="\t"))

	# Skip first row
	_a = list(a)[1:]

	# Track number of sites removed
	sitesRemoved = 0
	siteCount = 0
	sitesIgnored = 0
	sitesPruned = 0
	sitesUnique = []

	# Bind row values
	for i in range(0, len(_a)):
		_a[i] =  ['.'.join(posFileReader[i])] + list(_a[i])

	# Iterate through each rows by contig SNP counts
	for s in range(1, len(contigPosCounts)):

		# Early Break
		#if s > 10:
		#	break

		# Remember that position counts for each contig is stored on a per contig basis
		# Which means that when slicing the array, we need cumulative counts!
		_start = sum(contigPosCounts[0:s])
		_end = _start + contigPosCounts[s]

		_aList = _a[_start:_end]
		_contig = _aList[0][0].split('.')[0]
		print('Performing sequential pairwise comparison for contig '+_contig+' (#'+str(s+1)+') containing '+str(contigPosCounts[s])+' sites')
		#print('Start '+str(_start))
		#print('End '+str(_end))

		# Purify _aList if skipEnds is turned on
		if args.skipEnds > 0:

			# Store pruned list
			_aListPruned = []
			_pos = []

			# Iterate through list to store positions
			for t in _aList:
				pos = int(t[0].split('.')[1])
				_pos.append(pos)

			# Iterate through list again, this time to store items that are not in ends
			_sitesPruned = 0
			
			if _contig in contigLengths:
				_contigLength = contigLengths[_contig]
			else:
				print('Warning: Contig not found in dictionary of contig lengths: '+_contig)
				continue

			for u in _aList:
				pos = int(u[0].split('.')[1])
				if pos < args.skipEnds:
					print('Pruned position '+str(pos)+', close to start')
					sitesPruned += 1
					_sitesPruned += 1
				elif pos > max(_pos) - args.skipEnds:
					print('Pruned position '+str(pos)+', close to end')
					sitesPruned += 1
					_sitesPruned += 1
				else:
					_aListPruned.append(u)

			print('Pruned '+str(_sitesPruned)+' positions')

			_aList = _aListPruned

		# Make a copy of the list
		_bList = _aList[:]
		
		# Perform sequential pairwise comparisons
		identicalIndex = []
		for t in range(0, len(_aList)-1):

			# Ignore first column
			if _aList[t][1:] == _aList[t+1][1:]:
				#print('Identical sites: '+str(t)+' and '+str(t+1))
				identicalIndex.append(t+1)

		_aList = [v for i,v in enumerate(_aList) if i not in frozenset(identicalIndex)] 

		sitesRemoved += len(identicalIndex)
		print('Removed '+str(len(identicalIndex))+' sites')

		# If cherry picking is enabled
		if args.cherryPick < 1:
			
			if(len(_aList)*args.cherryPick < args.minSitesPerContig):
				__aList = _aList
				print('Insufficient sample size, using all')
			else:
				__aList = random.sample(_aList, int(len(_aList)*args.cherryPick))
				__aList.sort(key=lambda x: int(x[0].split('.')[1]))
				print('Picked '+str(len(__aList))+' positions out of '+str(len(_aList)))

		else:
			__aList = _aList

		# Write unique entries
		_aListUnique = []
		for u in __aList:
			if u[0] not in sitesUnique:

				# If not too many calls are missing
				if u[3:].count('X') >= len(u[3:])*args.missingThreshold:

					sitesIgnored += 1

				else:

					tempFileWriter.writerow(u)

					# Store position in unique list
					sitesUnique.append(u[0])
					siteCount += 1


		# Generate consensus genotype call per contig
		# In other words, we collapse all call positions in each contig into just one
		# _bList is a copy of _aList, which is not trimmed or pruned otherwise
		contigConsensus = []
		for i, site in enumerate(_bList):

			# If first iteration, create data structure
			if i == 0:
				contigConsensus = [{ 'calls': { 'A': 0, 'B': 0, 'X': 0, '-': 0 }, 'consensus': None, 'total': 0, 'confidence': 0 } for j in range(len(site[3:]))]

			# Iterate through all columns in each item
			# Exclude first 3 columns because they are contig information and parental genotype calls
			for j, call in enumerate(site[3:]):

				# Store counts
				contigConsensus[j]['calls'][call] += 1

				# Store total count
				contigConsensus[j]['total'] += 1

		# Iterate through contigConsensus list again, to determine consensus call per contig
		for pos in contigConsensus:

			calls = pos['calls']
			consensusCall = max(calls, key=lambda key: calls[key])

			# Set consensus genotype
			pos['consensus'] = consensusCall
			pos['confidence'] = float(pos['calls'][consensusCall]) / float(pos['total'])

		averageConfidence = numpy.mean([c['confidence'] for c in contigConsensus], dtype=np.float64)

		# Write to contig consensus temp file (for MSTmap)
		# Only pick genotype call positions that are identical
		contigConsensusCalls = [c['consensus'] for c in contigConsensus]
		sitesMatchingConsensus = 0

		# We create a regex pattern out of the contigCalls, since heterozygous genotypes (X) can match to both A and B
		contigCalls = re.sub('X', '[AB]' ,''.join(contigConsensusCalls))
		contigCallPattern = re.compile(contigCalls)

		if _contig == 'tig00000000':
			contigSampleFile.write(''.join(contigConsensusCalls)+'\n\n')

		if 'X' in contigConsensusCalls:
			print('Ambiguous contig consensus for '+_contig+': '+contigCalls)

		for k, site in enumerate(_bList):
			siteCalls = ''.join(site[3:])
			siteCongruent = []
			if contigCallPattern.match(siteCalls):
				siteCongruent.append(site)
				sitesMatchingConsensus += 1

			if _contig == 'tig00000000':
				contigSampleFile.write(siteCalls+'\n')

		# If cherry picking is enabled
		if args.cherryPick < 1:
			
			if(len(siteCongruent)*args.cherryPick < args.minSitesPerContig):
				_siteCongruent = siteCongruent
			else:
				_siteCongruent = random.sample(siteCongruent, int(len(siteCongruent)*args.cherryPick))
				_siteCongruent.sort(key=lambda x: int(x[0].split('.')[1]))

		else:
			_siteCongruent = siteCongruent

		# Iterate through all congruent sites
		for s in _siteCongruent:
			contigConsensusTempWriter.writerow(site)

		# Write to contig consensus file (for checking)
		consensusRate = float(sitesMatchingConsensus) / float(len(_bList))
		if _contig == 'tig00000000':
			contigSampleFile.write('\nTotal: ' + str(len(_bList)))
			contigSampleFile.write('\nSites matching consensus: ' + str(sitesMatchingConsensus))
		contigConsensusWriter.writerow([_contig] + contigConsensusCalls + [c['confidence'] for c in contigConsensus] + [averageConfidence] + [consensusRate])
		contigConsensusCache.append({ 'contig': _contig, 'consensus': contigConsensusCalls[2:] })

	# Summarise
	print('Total number of sites removed: '+str(sitesRemoved))
	print('Total number of sites ignored: '+str(sitesIgnored))
	print('Total number of sites pruned: '+str(sitesPruned))
	print('Total number of remaining sites: '+str(siteCount))

	# Write
	transposedFileWriter.writerow(['number_of_loci '+str(siteCount)])
	transposedFileWriter.writerow(['number_of_individual '+str(countRIL)])
	contigConsensusOutWriter.writerow(['number_of_loci '+str(siteCount)])
	contigConsensusOutWriter.writerow(['number_of_individual '+str(countRIL)])

	# Status update
	print('Writing to temp files complete, now concatenating...')

	# Close handlers
	transposedFile.close()
	tempFile.close()
	contigConsensusFile.close()
	contigConsensusOutFile.close()
	contigConsensusTempFile.close()

# Open files for concatenation
with \
open(os.path.join(outDir, 'refiltered.merged.012.transposed'), 'a') as transposedFile, \
open(os.path.join(outDir, 'refiltered.merged.012.tmp'), 'r') as tempFile, \
open(os.path.join(outDir, 'contigConsensus.txt'), 'a') as contigConsensusOutFile, \
open(os.path.join(outDir, 'contigConsensus.tmp'), 'r') as contigConsensusTempFile, \
open(os.path.join(outDir, 'contigConsensusMatrix.txt'), 'w') as contigConsensusMatrixFile:

	# Append tempFile to end of transposedFile
	print('Writing transposed file...')
	transposedFile.write('\n\n')
	for l in tempFile:
		transposedFile.write(l)
	print('Writing transposed file completed.\n')

	# Append tempFile to end of contigConsensusFile
	contigConsensusOutFile.write('\n\n')
	for m in contigConsensusTempFile:
		contigConsensusOutFile.write(m)

	matrixWriter = csv.writer(contigConsensusMatrixFile, lineterminator='\n')
	matrixWriter.writerow(['Contig1', 'Contig2', 'SimilarityScore', 'IdentityScore'])
	pairCount = 0
	print('Performing pairwise comparison to compute similarity matrix across contig consensus...')
	for pair in list(itertools.combinations(contigConsensusCache, 2)):
		c1 = pair[0]
		c2 = pair[1]

		# Compute similarity score
		cs1 = c1['consensus'][2:]
		cs2 = c2['consensus'][2:]

		similarity = 0
		identity = 0
		for i in range(0, len(cs1)):

			if cs1[i] in ['A','B', '-']:
				if cs1[i] == cs2[i]:
					similarity += 1
					identity += 1

			elif cs1[i] == 'X':
				if cs2[i] in ['A','B']:
					similarity += 0.5

				if cs2[i] == 'X':
					identity += 1

		similarityScore = float(similarity) / float(len(cs1))
		identityScore = float(identity) / float(len(cs1))

		matrixWriter.writerow([c1['contig'], c2['contig'], similarityScore, identityScore])

		pairCount += 1
		if pairCount % 100000 == 0:
			print('Compared pairs: '+str(pairCount))

	# Status update:
	print('Number of unique pairwise comparisons: '+str(pairCount))
	print('Writing completed :)')

# Remove temp files
os.remove(os.path.join(outDir, 'refiltered.merged.012.tmp'))
os.remove(os.path.join(outDir, 'contigConsensus.tmp'))

