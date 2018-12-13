import sys
import pandas
import numpy
from Bio import SeqIO

region_filename = sys.argv[1]
peaks_filename = sys.argv[2]
THRESHOLD = int(sys.argv[3])


# Read data frame from output of 'depth-high_coverage.py' and 'terminal-find_peaks.py'

region_file = pandas.read_csv(region_filename, names=('reference_name', 'beg', 'end'), delimiter='\t')
peaks_file = pandas.read_csv(peaks_filename, names=('reference_name', 'position'), delimiter='\t')


# Cut

reference_names = region_file.reference_name.drop_duplicates().values

for reference_name in reference_names:
	
	# Extract data frame for each reference.
	regions = region_file[region_file.reference_name.map(str) == str(reference_name)]
	peaks = peaks_file[peaks_file.reference_name.map(str) == str(reference_name)]
	
	# Cut each region.
	for index, region in regions.iterrows():
		
		# Append the beginning of the region to the list.
		positions = [region.beg]
		
		# Append the peaks in the region to the List.
		for index2, peak in peaks.iterrows():
			if region.beg < peak.position < region.end:
				positions.append(peak.position)
		
		# Append the end of the region to the list.
		positions.append(region.end + 1)
		
		# Print the region list.
		# chrom chromStat chromEnd name score strand thickStart thickEnd itemRgb
		for i in range(len(positions) - 1):
			if positions[i + 1] - positions[i] >= THRESHOLD:
				print(str(reference_name) + '\t' + str(positions[i]) + '\t' + str(positions[i + 1] - 1))

