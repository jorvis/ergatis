# This is component program to use the lgtanalysis lirary.
# AUTHOR: Ashwinkumar Ganesan.
# DATE  : 02-18-2011

from lgtanalysis import LGTNode
import sys
import os

FileLocation = sys.argv[1]
FilePath = FileLocation.split("/")
print FilePath

RefFileLocation1 = sys.argv[2]
print RefFileLocation1

RefFileLocation2 = sys.argv[3]
print RefFileLocation2

CheckSingleFastq = FileLocation + "_1.fastq"
print CheckSingleFastq

RefFilePtr1 = open(RefFileLocation,"r")
RefFilePtr2 = open(RefFileLocation,"r")

RefFile1List = []
RefFile2List = []

RefFile2UMUMCreated = []
RefFile2UMMCreated = []

for RefFile1 in RefFilePtr1:
   RefFile1List.append(RefFile1)
 
for RefFile2 in RefFilePtr2:
   RefFile2List.append(RefFile2)
   RefFile2UMUMCreated.append(0)
   RefFile2UMMCreated.append(0)

if(os.access(CheckSingleFastq, os.F_OK) != 0):
   FilePathString=""
   for i in range(1,len(FilePath)-1):
      FilePathString = FilePathString + "/" + FilePath[i]

   for RefFile1 in RefFile1List:
      # De-construct the reference file name.
      RefFileNameList1 = RefFile1.split("/")
      RefFileName1 = RefFileNameList1[len(RefFileNameList1) - 1]

      j = 0
      for RefFile2 in RefFile2List:
         # De-construct the reference file name.
         RefFileNameList2 = RefFile2.split("/")
         RefFileName2 = RefFileNameList2[len(RefFileNameList2) - 1]

         BacteriaFile = FilePathString + "/" + RefFileName1 + "_" + FilePath[len(FilePath) - 1] + ".sam" 
         print BacteriaFile
      
         SymbionFile = FilePathString + "/" + RefFileName2 + "_" + FilePath[len(FilePath) - 1] + ".sam" 
         print SymbionFile
      
         x=LGTNode(BacteriaFile,HumanFile)
         RefFile2UMUMCreated[j] = x.Find_UM_UM_Mates(0,1,"\t", RefFile2UMUMCreated[j])
         RefFile2UMMCreated[j] = x.Find_MM_UM_Mates(0,1,"\t", RefFile2UMMCreated[j])
         j += 1

         print("File Processing Done\n")
else:
   print("\nIt is a single read SRA file.. hence no post processing!!!")

# End of Script.
