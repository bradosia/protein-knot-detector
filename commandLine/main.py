import subprocess

spooky_path = 'C:/github/protein-knot-detector/commandLine/protein-knot-detector'
cmd = [spooky_path, '-l']
process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
output = process.communicate()[0]
print "Output:", output
process.wait()
print('Done')