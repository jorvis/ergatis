# This is progam that has function to work with a sequence file.
# Author: Ashwinkumar Ganesan.
# Date: 1/26/2010

# This is the class for analysing a sequence file.
# The input file is a SAM File. (Sequence Alignment Mapping File).
class SeqFileAnalyse:
   def __init__(self,FileLocation):
      self.FilePath = FileLocation
      FilePathContent = FileLocation.split('/')
      self.FileName = FilePathContent[len(FilePathContent) - 1]
       
      # Getting the File Directory.
      self.FileDirectory = FilePathContent[0]
      for iterator in range(1,len(FilePathContent)-1):
         self.FileDirectory = self.FileDirectory + "/" + FilePathContent[iterator]

      self.FileDescriptor = ''
      
      # This is a counter for different mapped combinations.
      self.ReadCounter = [0,0,0,0]

   def OpenFile(self):
      self.FileDescriptor = open(self.FilePath)

   def CloseFile(self):
      self.FileDescriptor.close()

   # This function is to Parse the Line from the File.
   # The DeLimiter should be in single quotes.
   def ParseLine(self, FileLine, DeLimiter):
      ParsedLine = FileLine.split(DeLimiter)
      NoOfEmptyValues = ParsedLine.count('')
         
      # This is to remove all the values that are null values in the list.
      for Iterator in range(NoOfEmptyValues):
         ParsedLine.remove('')
      
      return ParsedLine

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


   # This function is to check if the given input read is a singleton or not.
   # 0 means the Sequence is not a singleton.
   # 1 means the Sequence is a singleton.
   def ReadType(self, ParsedLine):
      # This is for SQ information at the start of the file.
      if(len(ParsedLine) < 11):
         return -1

      # Check the bitvalues in the Flag Integer.
      BinaryValList = self.Int2Binary(int(ParsedLine[1]))
      # BinaryValList.reverse()

      # Check for Cigar Value, Mapped Value and Pair Mapping.
      # *,* means both are unmapped.
      if(BinaryValList[2] == '1' and BinaryValList[3] == '1'):
         return 0
      # *,!* is a singleton for unmapped, mapped combination.
      elif(BinaryValList[2] == '1' and BinaryValList[3] == '0'):
         return 1
      # This is special case, as it means the next read is the singleton.
      # !*,* is a singleton for mapped, unmapped combination. (This means that the current read is properly mapped, but its pair is not).
      elif(BinaryValList[2] == '0' and BinaryValList[3] == '1'):
         return 2
      # !*,!* is a singleton for mapped, mapped combination.
      elif(BinaryValList[2] == '0' and BinaryValList[3] == '0'):
         return 3
      
      print ParsedLine
      print BinaryValList

   # This function parses a file with a specific delimiter and writes the contents to a new file.
   # Parameters: 1. DeLimiter - Is the character or set of characters based on which the File needs to be parsed.
   #             2. FileWrite - Is a flag if the file needs to be written to a file. It is a list. [<UM,UM>,<UM,M>,<M,UM>,<M,M>]
   #                            UM - Umapped Read and M - Mapped Read.
   #             3. ReadPercentage - Defines the minimum quality of the read. It takes a value from 0 to 1.
   #             4. FilterOrNot - Flag Value to check if filtering should be done or not. 0 - means filter and 1 - means do not filter.
   def ParseOut(self, DeLimiter, FileWrite, ReadPercentage, FilterOrNot):
      # Initialize Counts and open files for writing.
      # Initializations for singletons.
      self.ReadCounter = [0,0,0,0]

      # Unmapped - Unmapped
      if(FileWrite[0] == 1):
         UM_UM_FileLocation = self.FileDirectory + "/" + self.FileName + "_UM_UM.sam"
         UM_UM_File = open(UM_UM_FileLocation,'w')

      # Unmapped - Mapped
      if(FileWrite[1] == 1 or FileWrite[2] == 1):
         UM_M_FileLocation = self.FileDirectory + "/" + self.FileName + "_UM_M.sam"
         UM_M_File = open(UM_M_FileLocation,'w')

      # Mapped - Mapped
      if(FileWrite[3] == 1):
         M_M_FileLocation = self.FileDirectory + "/" + self.FileName + "_M_M.sam"
         M_M_File = open(M_M_FileLocation,'w')

      # The seek function is to start the parsing from the start of the file.
      # The file pointer may not be at the start if the function were called multiple times.
      self.FileDescriptor.seek(0,0)

      # Set Filter options.
      # If Filter option is disabled, then allow all reads to written to a file.
      if(FilterOrNot == 1):
         FilterRead = 1
         FilterReadMate1 = 1
         FilterReadMate2 = 1

      while 1:
         FileLine = self.FileDescriptor.readline()
         ParsedLine = self.ParseLine(FileLine, DeLimiter)
         if not FileLine:
            break
         pass
       
         # Find all combinations. 
         TypeCheck = self.ReadType(ParsedLine)
         if(TypeCheck != -1): # Not a Read if the value is -1.
            ReadPairLine = self.FileDescriptor.readline()
           
            # Filtering is being for only reads which have been mapped.
            # In case both reads are mapped, the quality of both reads, needs to be checked.
            # Do Filtering, if allowed.
            if(FilterOrNot == 0):
               if(TypeCheck == 1): # Check Mate 2.
                  ParsedPairLine = self.ParseLine(ReadPairLine, DeLimiter)
                  FilterRead = self.ReadFilter(ParsedPairLine, ReadPercentage)
               elif(TypeCheck == 2): # Check Mate 1.
                  FilterRead = self.ReadFilter(ParsedLine, ReadPercentage)
               elif(TypeCheck == 3): # Check both the pair.
                  ParsedPairLine = self.ParseLine(ReadPairLine, DeLimiter)
                  FilterReadMate1 = self.ReadFilter(ParsedLine, ReadPercentage)
                  FilterReadMate2 = self.ReadFilter(ParsedPairLine, ReadPercentage)

            if(FileWrite[TypeCheck] == 1):
               if(TypeCheck == 0):
                  self.ReadCounter[TypeCheck] = self.ReadCounter[TypeCheck] + 1  # Increment the counter according to the type of reading.
                  UM_UM_File.write(FileLine)
                  UM_UM_File.write(ReadPairLine)
               elif((TypeCheck == 1 or TypeCheck == 2) and FilterRead != 0):
                  self.ReadCounter[TypeCheck] = self.ReadCounter[TypeCheck] + 1  # Increment the counter according to the type of reading.
                  UM_M_File.write(FileLine)
                  UM_M_File.write(ReadPairLine)
               elif(TypeCheck == 3 and FilterReadMate1 != 0 and FilterReadMate2 != 0):
                  self.ReadCounter[TypeCheck] = self.ReadCounter[TypeCheck] + 1  # Increment the counter according to the type of reading.
                  M_M_File.write(FileLine)
                  M_M_File.write(ReadPairLine)
            
      # File Close.   
      if(FileWrite[0] == 1):
         UM_UM_File.close()
      elif(FileWrite[1] == 1 or FileWrite[2] == 1):
         UM_M_File.close()
      elif(FileWrite[3] == 1):
         M_M_File.close()

   # This function filters out low-quality reads.
   # The function looks for the number of character which are same in the read and the ratio of that to the total number of characters in the read.
   # If the ratio is higher than the ReadQuality value, then the read can be filtered. 

   # A better filter condition (later on) can the number of consecutive characters which are same.
   def ReadFilter(self, ParsedLine, ReadPercentage):
      # Convert the sequence string to a list and then count the elements.
      Seq = ParsedLine[9]
      SeqInList = list(Seq)
      SeqInListLen = len(SeqInList)
      
      if((float(SeqInList.count('T'))/float(SeqInListLen)) > ReadPercentage):
         return 0
      elif((float(SeqInList.count('G'))/float(SeqInListLen)) > ReadPercentage):
         return 0
      elif((float(SeqInList.count('C'))/float(SeqInListLen)) > ReadPercentage):
         return 0
      elif((float(SeqInList.count('A'))/float(SeqInListLen)) > ReadPercentage):
         return 0

      return 1

   # This is the function print the number of the type of pairs found in the SAM file.
   def PrintStat(self):
      print "Number of pairs with both Mates Unmapped (UM - UM): " + str(self.ReadCounter[0])
      print "Number of pairs with First Mate Unmapped and Second Mate Mapped (UM - M): " + str(self.ReadCounter[1])
      print "Number of pairs with First Mate Mapped and Second Mate Unmapped (M - UM): " + str(self.ReadCounter[2])
      print "Number of pairs with both Mates Mapped (M - M): " + str(self.ReadCounter[3]) + "\n"

      CounterFileName = self.FileDirectory + "/" + self.FileName + "_mapping_numbers.txt"
      CounterFilePointer = open(CounterFileName,"w")
      WriteString = str(self.ReadCounter[0]) + "," + str(self.ReadCounter[1]) + "," + str(self.ReadCounter[2]) + "," + str(self.ReadCounter[3])
      CounterFilePointer.write(WriteString)
      CounterFilePointer.close()

# End of Program.
