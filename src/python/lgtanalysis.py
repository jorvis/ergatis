# This is a program to analyse SAM files for Lateral Genome Transfers.
# Ashwinkumar Ganesan.
# Date: 01/25/2010

# This module is access os functions, like check if a file exists or not.
import os

# Class to work with SAM Files.
from seqfileparse import SeqFileAnalyse

# Class for LGT.
# It contains the function to analyse and create statistics for Lateral Genome Transfer.
# The two SAM Files are created using the same comparison genome, but reference files for the alignment are different.
# SAMFile1 and SAMFile2 are the location of the two sam files which are used to analysing lateral genome transfers.
class LGTNode:
   def __init__(self, SAMFile1, SAMFile2):
      self.FirstFileLocation = SAMFile1
      self.SecondFileLocation = SAMFile2
      
      # 0 index position is for combinations where First File Mate 1 is mapped and Second File Mate 2 is mapped.
      # 1 index position is for combinations where First File Mate 2 is mapped and Second File Mate 1 is mapped.
      self.UM_M_UM_M_MapCounter = [0,0]
      self.M_M_UM_M_MapCounter = [0,0]
 
   # This a function to convert a number into binary form.
   # The number of bits in the number is upto 16.
   # The function returns the binary value as a list of bit values.
   def Int2Binary(self, Number):
      BinaryValue = []
      while(Number > 0):
         Rem = Number % 2
         BinaryValue.append(str(Rem))
         Number = Number / 2
         
      while(len(BinaryValue) < 16):
         BinaryValue.append('0')
            
      return BinaryValue

   # This function is to index a file.
   # The index works only on CSV Type of files.
   # Index only the first read in the pair.
   def IndexFile(self, FileName, IndexColumn, DeLimiter):
      FilePointer = open(FileName, "r")
      SamAnalyse = SeqFileAnalyse(FileName)
      ReadDict = {}

      # Read all the contents of the file and put in an array.
      while 1:
         Position1 = FilePointer.tell()
         FileLine1 = FilePointer.readline()
         ParsedFileLine1 = SamAnalyse.ParseLine(FileLine1, DeLimiter)
     
         # Exit the loop if all the lines in the file are read.
         if not FileLine1:
            break
         pass
         
         # Read the other Read-Pair.
         FileLine2 = FilePointer.readline()
          
         # Add the value to Dictionary.
         Index = ParsedFileLine1[IndexColumn]
         ReadDict[Index] = Position1

      FilePointer.close()
      return ReadDict

      
   # This is a function to find the combinations of mapping amongst sam files for those reads
   # where one mate is mapped and the other mate is unmapped.

   # Since the Reads are in a sequence in the SAM file, once a read is found then, searching can continue in from the same location.
   # U - Unmapped.
   # M - Mapped.

   # Paramters: 1. SamFile1FilterOrNot - Flag option to decide, if the SAM File should consider a filtered output or unfiltered output.
   #            0 - to filter, 1 - Not filter.
   #            e.g - This can be used for Human Genome (req. unfiltered SAM File) to Bacterial compared SAM File which should be filtered. 
   # CreateAlways is a flag - so that the Sam file parsing is done everytime.
   def Find_UM_UM_Mates(self,SamFile1FilterOrNot, SamFile2FilterOrNot, DeLimiter, CreateAlways):
      self.UM_M_UM_M_MapCounter = [0,0]
      Location = self.FirstFileLocation + "_UM_M_UM_M.sam"
      MappedMatesFile = open(Location,"w")

      FirstFileLocation = self.FirstFileLocation + "_UM_M.sam"
      SamAnalysis1 = SeqFileAnalyse(self.FirstFileLocation)
      # "False" means the file does not exist. Hence create the file.
      if(os.access(FirstFileLocation, os.F_OK) == 0):
         SamAnalysis1.OpenFile()
         SamAnalysis1.ParseOut(DeLimiter,[0,1,1,1],0.8,SamFile1FilterOrNot)
         SamAnalysis1.PrintStat()
         SamAnalysis1.CloseFile()

      SecondFileLocation = self.SecondFileLocation + "_UM_M.sam"
      SamAnalysis2 = SeqFileAnalyse(self.SecondFileLocation)
      # "False" means the file does not exist.
      if(os.access(SecondFileLocation, os.F_OK) == 0 and (CreateAlways == 0)):
         SamAnalysis2.OpenFile()
         SamAnalysis2.ParseOut(DeLimiter,[1,1,1,0],0.8,SamFile2FilterOrNot)
         SamAnalysis2.PrintStat()
         SamAnalysis2.CloseFile()

      # Index the Files.
      FirstFileDict = self.IndexFile(FirstFileLocation, 0, DeLimiter)
      SecondFileDict = self.IndexFile(SecondFileLocation, 0, DeLimiter)

      # Once the Files are created, verify the values.
      FirstMappedFile = open(FirstFileLocation,'r')
      SecondMappedFile = open(SecondFileLocation,'r')
     
      # The DeLimiter has to be the same accross the files. 
      for Key in FirstFileDict:
         if((Key in SecondFileDict) == 1):
            # Get the location of the read in the first file.
            FirstFileReadLocation = FirstFileDict[Key]
            FirstMappedFile.seek(FirstFileReadLocation, 0)
            
            # Read the First File pair.
            FirstFileMate1 = FirstMappedFile.readline()
            ParsedFirstFileMate1 = SamAnalysis1.ParseLine(FirstFileMate1, DeLimiter)
            FirstFileMate2 = FirstMappedFile.readline()


            # Get the location of the read in the second file.
            SecondFileReadLocation = SecondFileDict[Key]
            SecondMappedFile.seek(SecondFileReadLocation, 0)

            # Read Second File Pair.
            SecondFileMate1 = SecondMappedFile.readline()
            ParsedSecondFileMate1 = SamAnalysis2.ParseLine(SecondFileMate1, DeLimiter)
            SecondFileMate2 = SecondMappedFile.readline()
           
            # Using Binary values to check between singletons.
            FirstFileMate1BinVal = self.Int2Binary(int(ParsedFirstFileMate1[1]))
            SecondFileMate1BinVal = self.Int2Binary(int(ParsedSecondFileMate1[1]))

            if(FirstFileMate1BinVal[2] == '1' and SecondFileMate1BinVal[3] == '1'):
               print(str(FirstFileMate1BinVal) + " " + str(SecondFileMate1BinVal))
               self.UM_M_UM_M_MapCounter[0] = self.UM_M_UM_M_MapCounter[0] + 1
               MappedMatesFile.write(FirstFileMate1)
               MappedMatesFile.write(FirstFileMate2)
               MappedMatesFile.write(SecondFileMate1)
               MappedMatesFile.write(SecondFileMate2)
            elif(FirstFileMate1BinVal[3] == '1' and SecondFileMate1BinVal[2] == '1'):
               print(str(FirstFileMate1BinVal) + " " + str(SecondFileMate1BinVal))
               self.UM_M_UM_M_MapCounter[1] = self.UM_M_UM_M_MapCounter[1] + 1
               MappedMatesFile.write(FirstFileMate1)
               MappedMatesFile.write(FirstFileMate2)
               MappedMatesFile.write(SecondFileMate1)
               MappedMatesFile.write(SecondFileMate2)

      MappedMatesFile.close()
      print str(self.UM_M_UM_M_MapCounter[0]) + " " + str(self.UM_M_UM_M_MapCounter[1])

      CounterFile = self.FirstFileLocation + "_UM_M_UM_M_Numbers.txt"
      CounterFilePointer = open(CounterFile,"w")
      WriteString = str(self.UM_M_UM_M_MapCounter[0]) + "," + str(self.UM_M_UM_M_MapCounter[1])
      CounterFilePointer.write(WriteString)
      CounterFilePointer.close()

      # Close the other files.
      FirstMappedFile.close()
      SecondMappedFile.close()

      return 1


   # This function looks for the combinations where the first file has a read pair, where both are mapped.

   def Find_MM_UM_Mates(self,SamFile1FilterOrNot, SamFile2FilterOrNot, DeLimiter, CreateAlways):
      self.M_M_UM_M_MapCounter = [0,0]
      Location = self.FirstFileLocation + "_M_M_UM_M.sam"
      MappedMatesFile = open(Location,"w")

      FirstFileLocation = self.FirstFileLocation + "_M_M.sam"
      SamAnalysis1 = SeqFileAnalyse(self.FirstFileLocation)
      # "False" means the file does not exist. Hence create the file.
      if(os.access(FirstFileLocation, os.F_OK) == 0):
         SamAnalysis1.OpenFile()
         SamAnalysis1.ParseOut(DeLimiter,[1,1,1,1],0.8,SamFile1FilterOrNot)
         SamAnalysis1.PrintStat()
         SamAnalysis1.CloseFile()

      SecondFileLocation = self.SecondFileLocation + "_UM_M.sam"
      SamAnalysis2 = SeqFileAnalyse(self.SecondFileLocation)
      # "False" means the file does not exist.
      if(os.access(SecondFileLocation, os.F_OK) == 0 and (CreateAlways == 0)):
         SamAnalysis2.OpenFile()
         SamAnalysis2.ParseOut(DeLimiter,[1,1,1,1],0.8,SamFile2FilterOrNot)
         SamAnalysis2.PrintStat()
         SamAnalysis2.CloseFile()

      # Index the Files.
      FirstFileDict = self.IndexFile(FirstFileLocation, 0, DeLimiter)
      SecondFileDict = self.IndexFile(SecondFileLocation, 0, DeLimiter)

      # Once the Files are created, verify the values.
      FirstMappedFile = open(FirstFileLocation,'r')
      SecondMappedFile = open(SecondFileLocation,'r')
     
      # The DeLimiter has to be the same accross the files. 
      for Key in FirstFileDict:
         if((Key in SecondFileDict) == 1):
            # Get the location of the read in the first file.
            FirstFileReadLocation = FirstFileDict[Key]
            FirstMappedFile.seek(FirstFileReadLocation, 0)
            
            # Read the First File pair.
            FirstFileMate1 = FirstMappedFile.readline()
            FirstFileMate2 = FirstMappedFile.readline()

            # Get the location of the read in the second file.
            SecondFileReadLocation = SecondFileDict[Key]
            SecondMappedFile.seek(SecondFileReadLocation, 0)

            # Read Second File Pair.
            SecondFileMate1 = SecondMappedFile.readline()
            SecondFileMate2 = SecondMappedFile.readline()
           
            self.M_M_UM_M_MapCounter[0] = self.M_M_UM_M_MapCounter[0] + 1
            MappedMatesFile.write(FirstFileMate1)
            MappedMatesFile.write(FirstFileMate2)
            MappedMatesFile.write(SecondFileMate1)
            MappedMatesFile.write(SecondFileMate2)

      MappedMatesFile.close()
      print str(self.M_M_UM_M_MapCounter[0]) + " " + str(self.M_M_UM_M_MapCounter[1])

      CounterFile = self.FirstFileLocation + "_M_M_UM_M_Numbers.txt"
      CounterFilePointer = open(CounterFile,"w")
      WriteString = str(self.M_M_UM_M_MapCounter[0]) + "," + str(self.M_M_UM_M_MapCounter[1])
      CounterFilePointer.write(WriteString)
      CounterFilePointer.close()

      # Close the other Files. 
      FirstMappedFile.close()
      SecondMappedFile.close()

      return 1

# End Of Class.
# End Of Program      
