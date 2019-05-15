import sys
import wiggelen

wig_fname = sys.argv[1]
threshold = int(sys.argv[2])	# the lowest coverage of high-coverage region.

# Read data frame from output of 'samtools depth'.

wig_file = wiggelen.walk(open(wig_fname))


# Set the beginning state 'LOW'.

state = 'LOW'

# Search and print regions.
for reference, position, value in wig_file:
	
	if state == 'LOW':
		if value >= threshold:
			reference_name = reference
			beg = position
			end = position
			state = 'HIGH'
	
	elif state == 'HIGH':
		
		# The next reference has started.
		if reference_name != reference:
			print(str(reference_name) + '\t' + str(beg) + '\t' + str(end))
			state = 'LOW'
		
		# THe high-coverage region is continuing.
		elif value >= threshold:
			end = position
		
		# The coverage is too small to permit.
		elif value < threshold:
			print(str(reference_name) + '\t' + str(beg) + '\t' + str(end))
			state = 'LOW'
# Print the last region.
if(reference_name != None):
	print(str(reference_name) + '\t' + str(beg) + '\t' + str(end))

