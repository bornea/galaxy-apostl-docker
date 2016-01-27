#######################################################################################
# Python-code: Shiny Bubblebeam wrapper
# Author: Adam L Borne
# Contributers: Paul A Stewart, Brent Kuenzi
#######################################################################################
# This program runs the R script that generates a bubble plot in shiny. Generates
# a unique app for each run of the tool for galaxy integration.
#######################################################################################
# Copyright (C)  Adam Borne.
# Permission is granted to copy, distribute and/or modify this document
# under the terms of the GNU Free Documentation License, Version 1.3
# or any later version published by the Free Software Foundation;
# with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.
# A copy of the license is included in the section entitled "GNU
# Free Documentation License".
#######################################################################################
## REQUIRED INPUT ##

# 1) list_file: SaintExpress output file.
# 2) prey_file: Prey file listing gene name, sequence legnth, and gene id.
# 3) crapome: Crapome file can be created at http://crapome.org. (default = "None")
#######################################################################################
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
    x.write("<html><body> open <a href=\"http://proteomicspipeline.moffitt.usf.edu:3838/" +
            str(stamped_app) + r"\">Shiny Bubblebeam</a> in your browser to view shiny app.</body></html>")

os.rename('shiny.txt', str(sys.argv[4]))
