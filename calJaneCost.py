import os
import os.path
from cStringIO import StringIO

def main():

	for i in xrange(100):

		index = str(i + 1)

		for j in xrange(4 - len(str(i + 1))):
			index = "0" + index

		f = "janeOut100/COG" + index + ".txt"

		if not os.path.isfile(f):
			continue

		inFile = open(f, 'r')

		ans = 0
		d = 2
		t = 3
		l = 1

		for line in inFile:
			if line.startswith("Duplication: "):
				num = int(line[13:])
				ans += num * d
				print num
			elif line.startswith("Host Switch: "):
				num = int(line[12:])
				ans += num * t
				print num
			elif line.startswith("Loss: "):
				num = int(line[6:])
				ans += num * l
				print num

		print ans

		inFile.close()
		
		outFile = open(f, 'a')
		outFile.write(str(ans))
		outFile.close()

		print "Done ", index



if __name__ == "__main__": 
	main()