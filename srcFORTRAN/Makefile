.PHONY: clean all

all:
	(cd lapack-3.10.0;	make)
	(cd srcNMCT;        make)
    
clean:
	-rm -rf *.exe tl.bin

cleanall:
	(cd lapack-3.10.0;	make clean)
	(cd srcNMCT;	    make clean)  
