CPP = g++ 
#CPP = g++ -pg
FLAGS = -O9 -Wall 
#FLAGS = -g -O0 -Wall -DNDEBUG 
OBJS = 	./src/basics/io.o \
				./src/basics/buffer.o\
				./src/basics/BitSequenceRG.o\
				./src/coders/Coder.o \
				./src/coders/THuff.o \
				./src/coders/Huffman.o \
				./src/cqual/CQualPBlock.o\
				./src/cqual/CQualRBlock.o\
				./src/cqual/CQualLL.o\
				./src/cqual/CQualBlock.o\
				./src/crps/CRPSPreMF.o\
				./src/crps/CRPSPreMFAll.o\
				./src/crps/CRePoSe.o\
				./src/cfields/CFields.o\
				./src/csam/CSAM.o\


BIN = CompressSAM DecompressSAM CompressQual DecompressQual CompressSeq DecompressSeq GetIntervalSeq GetIntervalSeqSample CountReadsSample GetIntervalSAM GetIntervalSAMSample GetIntervalSSN
%.o: %.cpp
	@echo " [C++] Compiling $<"
	@$(CPP) $(FLAGS) -c $< -o $@ 

all: clean $(OBJS) $(BIN)


IntervalTest:
	$(CPP) $(FLAGS) $(OBJS) -o IntervalTest ./tests/IntervalTest.cpp -lboost_iostreams

CompressQual:
	$(CPP) $(FLAGS) $(OBJS) -o CompressQual ./tests/CompressQual.cpp -lboost_iostreams

DecompressQual:
	$(CPP) $(FLAGS) $(OBJS) -o DecompressQual ./tests/DecompressQual.cpp -lboost_iostreams

CompressSeq:
	$(CPP) $(FLAGS) $(OBJS) -o CompressSeq ./tests/CompressSeq.cpp -lboost_iostreams

DecompressSeq:
	$(CPP) $(FLAGS) $(OBJS) -o DecompressSeq ./tests/DecompressSeq.cpp -lboost_iostreams

CompressSAM:
	$(CPP) $(FLAGS) $(OBJS) -o CompressSAM ./tests/CompressSAM.cpp -lboost_iostreams

DecompressSAM:
	$(CPP) $(FLAGS) $(OBJS) -o DecompressSAM ./tests/DecompressSAM.cpp -lboost_iostreams

GetIntervalSeq:	 
	$(CPP) $(FLAGS) $(OBJS) -o GetIntervalSeq ./tests/GetIntervalSeq.cpp -lboost_iostreams

GetIntervalSeqSample:  
	$(CPP) $(FLAGS) $(OBJS) -o GetIntervalSeqSample ./tests/GetIntervalSeqSample.cpp -lboost_iostreams

CountReadsSample:
	$(CPP) $(FLAGS) $(OBJS) -o CountReadsSample ./tests/CountReadsSample.cpp -lboost_iostreams

GetIntervalSAM:  
	$(CPP) $(FLAGS) $(OBJS) -o GetIntervalSAM ./tests/GetIntervalSAM.cpp -lboost_iostreams

GetIntervalSAMSample:  
	$(CPP) $(FLAGS) $(OBJS) -o GetIntervalSAMSample ./tests/GetIntervalSAMSample.cpp -lboost_iostreams

GetIntervalSSN:
	$(CPP) $(FLAGS) $(OBJS) -o GetIntervalSSN ./tests/GetIntervalSSN.cpp -lboost_iostreams

doc:
	@echo " [DOC] Generating documentation"
	@doxygen


clean:
	@echo " [CLN] Removing object files"
	@rm -f $(OBJS) $(BIN) *~

