import os
import os.path
from cStringIO import StringIO

def main():

	for i in xrange(100):

		index = str(i + 1)

		for j in xrange(4 - len(str(i + 1))):
			index = "0" + index

		inFile = "real-100taxa/COG" + index + ".newick"

		if not os.path.isfile(inFile):
			continue

		os.system("./Jane/jane-cli.sh -p 100 -i 100 -c 0 2 3 1 0 treeFiles/COG" + index + ".tree > janeOut100/COG" + index + ".txt")

if __name__ == "__main__": 
	main()