import sys
import pysam

bam_filename = sys.argv[1]


# Set target_id, reference_length and an empty array from reference_name

def reference_initialize(bam_file, reference_name):
	target_id = bam_file.get_tid(reference_name)
	reference_length = bam_file.header['SQ'][target_id]['LN']
	terminal_count = [0 for i in range(reference_length)]
	
	return target_id, reference_length, terminal_count


# Output terminal_count.

def print_count(reference_name, reference_length, terminal_count):
	for i in range(reference_length):
		print(reference_name + '\t' + str(i + 1) + '\t' + str(terminal_count[i])) # 1-origin


# main

bam_file = pysam.AlignmentFile(bam_filename)

for read in bam_file:
	
	if 'reference_name' not in locals():
		# Initialize the terminal count.
		reference_name = read.reference_name
		target_id, reference_length, terminal_count = reference_initialize(bam_file, reference_name)
	
	elif read.reference_id == -1:
		# End of alignment.
		break
	
	elif reference_name != read.reference_name:
		# Print the terminal count.
		print_count(reference_name, reference_length, terminal_count)
		
		# Initialize the terminal count.
		reference_name = read.reference_name
		target_id, reference_length, terminal_count = reference_initialize(bam_file, reference_name)
	
	# Calculate the alignment terminals from pos and CIGAR.
	beg = read.reference_start   # 0-origin
	end = read.reference_end - 1 # 0-origin
	
	# Increase the terminal count in the positions of alignment terminals.
	terminal_count[beg] += 1
	terminal_count[end] += 1

# Print the terimal count of the last reference.
print_count(reference_name, reference_length, terminal_count)


