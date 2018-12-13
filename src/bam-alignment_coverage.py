import sys
import pysam

bam_filename = sys.argv[1]


# Set target_id, reference_length and an empty array from reference_name

def reference_initialize(bam_file, reference_name):
	target_id = bam_file.get_tid(reference_name)
	reference_length = bam_file.header['SQ'][target_id]['LN']
	coverage_count = [0 for i in range(reference_length)]
	
	return target_id, reference_length, coverage_count


# Output coverage_count.

def print_count(reference_name, reference_length, coverage_count):
	print('variableStep chrom=' + reference_name)
	for i in range(reference_length):
		print(str(i + 1) + '\t' + str(coverage_count[i])) # 1-origin


# main

bam_file = pysam.AlignmentFile(bam_filename)

for read in bam_file:
	if 'reference_name' not in locals():
		# Initialize the coverage count.
		reference_name = read.reference_name
		target_id, reference_length, coverage_count = reference_initialize(bam_file, reference_name)
	
	#elif read.reference_id == -1:
	#	# End of alignment.
	#	print_count(reference_name, reference_length, coverage_count)
	#	break
	
	elif reference_name != read.reference_name:
		# Print the coverage count.
		print_count(reference_name, reference_length, coverage_count)
		
		# Initialize the coverage count.
		reference_name = read.reference_name
		target_id, reference_length, coverage_count = reference_initialize(bam_file, reference_name)
	
	# Calculate the alignment terminals from pos and CIGAR.
	beg = read.reference_start   # 0-origin
	end = read.reference_end - 1 # 0-origin
	
	# Increase the coverage count in the positions of alignment coverages.
	for i in range(beg, end):
		coverage_count[i] += 1
print_count(reference_name, reference_length, coverage_count)



