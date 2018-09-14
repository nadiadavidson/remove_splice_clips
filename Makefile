SEQ=../SeqLib
INC=-I$(SEQ) -I$(SEQ)/htslib
LD=$(SEQ)/bin/libseqlib.a $(SEQ)/bin/libbwa.a $(SEQ)/bin/libfml.a $(SEQ)/bin/libhts.a

all: remove_splice_clips fix_contam_clips

remove_splice_clips: remove_splice_clips.c++
	g++ -O3 $(INC) -o $@ $< $(LD) -lz -lbz2 -llzma

fix_contam_clips: fix_contam_clips.c++
	g++ -O3 $(INC) -o $@ $< $(LD) -lz -lbz2 -llzma

clean:                   
	-rm *~ remove_splice_clips

