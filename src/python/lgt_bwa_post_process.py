# This is component program to use the lgtanalysis lirary.
# AUTHOR: Ashwinkumar Ganesan.
# DATE  : 02-18-2011

from lgtanalysis import LGTNode
import sys
import os
import argparse
import glob
import re

RefFile1List = []
RefFile2List = []
RefFile2UMUMCreated = []
RefFile2UMMCreated = []

class ParseListFile(argparse.Action):
    def __call__(self,parser,namespace,values,option_string=None):
        if(values != ''):
            FILE = open(values,"r")
            retlist = []
            for line in FILE:
                line.rstrip('\n')
                retlist.append(line)
            if(self.dest == 'ref_list1'):
                RefFile1List.extend(retlist)
            if(self.dest == 'ref_list2'):
                RefFile2UMUMCreated.append(0)
                RefFile2UMMCreated.append(0)
                RefFile2List.extend(retlist)
                setattr(namespace, self.dest,retlist);

class ParseList(argparse.Action):
    def __call__(self,parser,namespace,values,option_string=None):
       if(values != ''):
           retlist = values.split(',')
           if(self.dest == 'ref1'):
               RefFile1List.extend(retlist)
           if(self.dest == 'ref2'):
               RefFile2UMUMCreated.append(0)
               RefFile2UMMCreated.append(0)
               RefFile2List.extend(retlist)
               setattr(namespace, self.dest, retlist)

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--output')
parser.add_argument('-rl1', '--ref_list1',action=ParseListFile)
parser.add_argument('-rl2', '--ref_list2',action=ParseListFile)
parser.add_argument('-r1', '--ref1',action=ParseList)
parser.add_argument('-r2', '--ref2',action=ParseList)
parser.add_argument('-i', '--input')
args = parser.parse_args();

#FileLocation = args.input
regex = re.compile('\.lite')
FileLocation = regex.sub('',args.input)
print FileLocation
#FileLocation = sys.argv[1]
FilePath = FileLocation.split("/")
FilePathString=""
for i in range(1,len(FilePath)-1):
  FilePathString = FilePathString + "/" + FilePath[i]

FileName = FilePath[len(FilePath)-1]
#RefFileLocation1 = sys.argv[2]
#RefFileLocation1 = 
#print RefFileLocation1

#RefFileLocation2 = sys.argv[3]
#print RefFileLocation2

SingleEnd = glob.glob(FilePathString + "/*" + FileName + "_single_read.txt")
print SingleEnd

#RefFilePtr1 = open(RefFileLocation1,"r")
#RefFilePtr2 = open(RefFileLocation2,"r")

#RefFile2UMUMCreated = []
#RefFile2UMMCreated = []

print RefFile1List
print RefFile2List
#for RefFile1 in RefFilePtr1:
#   RefFile1List.append(RefFile1)
 
#for RefFile2 in RefFilePtr2:
#   RefFile2List.append(RefFile2)
#   RefFile2UMUMCreated.append(0)
#   RefFile2UMMCreated.append(0)

CheckSingleEnd = all((os.access(f, os.F_OK) for f in SingleEnd))
if(len(SingleEnd) == 0 and CheckSingleEnd):
   for RefFile1 in RefFile1List:
      print "File 1: "+RefFile1
      # De-construct the reference file name.
      RefFileNameList1 = RefFile1.split("/")
      RefFileName1 = RefFileNameList1[len(RefFileNameList1) - 1]
      RefFileName1 = os.path.splitext(RefFileName1)[0]
      j = 0
      for RefFile2 in RefFile2List:
         print "File 2: "+RefFile2
         # De-construct the reference file name.
         RefFileNameList2 = RefFile2.split("/")
         RefFileName2 = RefFileNameList2[len(RefFileNameList2) - 1]
         RefFileName2 = os.path.splitext(RefFileName2)[0]
         BacteriaFile = FilePathString + "/" + RefFileName1 + "_" + FilePath[len(FilePath) - 1] + ".sam" 
         print BacteriaFile
      
         SymbionFile = FilePathString + "/" + RefFileName2 + "_" + FilePath[len(FilePath) - 1] + ".sam" 
         print SymbionFile
      
         x=LGTNode(BacteriaFile,SymbionFile)
         RefFile2UMUMCreated[j] = x.Find_UM_UM_Mates(0,1,"\t", RefFile2UMUMCreated[j])
         RefFile2UMMCreated[j] = x.Find_MM_UM_Mates(0,1,"\t", RefFile2UMMCreated[j])
         j += 1

         print("File Processing Done\n")
else:
   print("\nIt is a single read SRA file.. hence no post processing!!!")

# End of Script.
