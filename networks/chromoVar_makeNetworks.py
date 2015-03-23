from optparse import OptionParser
import os
import math
from time import gmtime, strftime
import re
import gzip
import random
import numpy 
'''
Author:Oana Ursu
'''

def main():
    parser=OptionParser()
    
    parser.add_option('--peak2peak',dest='peak2peak',default='')
    parser.add_option('--mode',dest='mode',help='peak2peak, network')
    opts,args=parser.parse_args()
   

    if opts.mode=='peak2peak':
        #read table into a 






















def qsub_a_command(cmd,shell_script_name,split_string=',',memory_number='20G'):
    f=open(shell_script_name,'w')
    cmds=cmd.split(split_string)
    for i in range(len(cmds)):
        f.write(cmds[i]+'\n') #just write the command
    f.close()
    #make runnable                                                                                                                                                      
    os.system('chmod 711 '+shell_script_name)
    #Qsub the script                                                                                                                                                    
    os.system("qsub -l mem_free="+memory_number+" -l h_vmem="+memory_number+" -l h_rt=90:00:00 -o "+shell_script_name+'.o'+' -e '+shell_script_name+'.e'+' '+shell_script_name)
        
            

main()
