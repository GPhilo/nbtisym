all: parser.cc parser.hh scanner.cc

parser.cc parser.hh: parser.yy
	bison --defines=parser.hh -oparser.cc parser.yy

scanner.cc: scanner.ll
	flex -oscanner.cc scanner.ll

.PHONY: 
	clean

clean:
	-rm *.o parser.hh parser.cc scanner.cc location.hh position.hh stack.hh calc++.out

