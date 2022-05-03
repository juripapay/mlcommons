# Remove tabs from python file.

inputFile = open('stemdl_light.py', "r")
exportFile = open('stemdl_light_1.py', "w")
for line in inputFile:
   new_line = line.replace('\t', '    ')
   exportFile.write(new_line)

inputFile.close()
exportFile.close()

