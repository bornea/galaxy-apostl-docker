import os 
import sys
import time

input_list = open(sys.argv[1], 'r')
prey_input = open(sys.argv[2], 'r')
stamped_app = r"shiny_bubble" + str(time.strftime('_%d_%m_%Y_%H_%M')) 
cmd = r"cp -r /srv/shiny-server/shiny_bubble /srv/shiny-server/" + str(stamped_app) 
os.system(cmd)

if sys.argv[3] != 'None':
	crapome = open(sys.argv[3], 'r')
	crap_file = open('/srv/shiny-server/'+ str(stamped_app) + '/craptest.txt', 'w')
	glob_manip = open('/srv/shiny-server/shiny_bubble/global.R', 'r')
	glob_write = open('/srv/shiny-server/'+ str(stamped_app) + '/global.R', 'w')
	for code_line in glob_manip:
		if r"main.data <- as.data.frame\(merge_files" in code_line:
			glob_write.write(r"main.data <- as.data.frame(merge_files(\"test_list.txt\", \"preytest.txt\", \"craptest.txt\"))")
		else:
			glob_write.write(code_line)
	for line in crapome:
		crap_file.write(line)

input_file = open('/srv/shiny-server/'+ str(stamped_app) + '/test_list.txt', 'w')
for line in input_list:
	input_file.write(line)
prey_file = open('/srv/shiny-server/'+ str(stamped_app) + '/preytest.txt', 'w')
for line in prey_input:
	prey_file.write(line)




#cmd1 = r"touch '/srv/shiny-server/" + str(stamped_app) + r"/restart.txt"
#os.system(cmd1)

with open("shiny.txt", "wt") as x:
	x.write("<html><body> open <a href=\"http://localhost:3838/"+ str(stamped_app) + "\">Shiny Bubblebeam</a> in your browser to view shiny app. If there are issues with the sizing within galaxy you can right click and open in a new tab or window.</body></html>")

os.rename('shiny.txt', str(sys.argv[4]))
