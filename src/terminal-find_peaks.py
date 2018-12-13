import sys
import pandas
import scipy.stats as st
import scipy.signal as sg
import matplotlib.pyplot as plt

terminal_filename = sys.argv[1]
ORDER = int(sys.argv[2]) # scipy.signal.argrelmax(..., oreder=ORDER)
THRESHOLD = int(sys.argv[3])

PRINTPEAKS = True
SHOWPLOT = False



# Read data frame from output of 'bam-alignment_terminal.py'

terminal_file = pandas.read_csv(terminal_filename, names=('reference_name', 'position', 'terminal_count'), delimiter='\t')

reference_names = terminal_file.reference_name.drop_duplicates()


# Find peaks for each references.

for reference in reference_names:
	
	print('Find peaks for reference "' + str(reference) + '"...', file=sys.stderr)
	
	# Extract data frame for each references.
	terminal = terminal_file[terminal_file.reference_name == reference]
	count = terminal.terminal_count.values
	
	# Find peaks with scipy.signal.argrelmax().
	# I don't know why '[0]' is needed...
	peaks = sg.argrelmax(count, order=ORDER)[0]
	
	# Extract enough high peaks.
	peaks = peaks[count[peaks] >= THRESHOLD]
	
	# If PRINTPEAKS is True, print peaks to stdout.
	if PRINTPEAKS:
		for peak in peaks:
			print(str(reference) + '\t' + str(peak))
	
	# If SHOWPLOT is True, show plot with matplotlib.
	if SHOWPLOT:
		plt.bar(terminal.position, terminal.terminal_count)
		plt.plot(terminal.position.values[peaks], terminal.terminal_count.values[peaks], 'o')
		
		plt.title(str(reference))
		plt.xlabel('Reference Position')
		plt.ylabel('Terminal Count')
		
		plt.savefig(str(reference).replace('/', '__') + '.png')
		plt.close()

