all: make_gcell make_aspect

make_gcell:
	make -f makefile_g

make_aspect:
	make -f makefile_a

clean_gcell:
	make -f makefile_g clean

clean_aspect:
	make -f makefile_a clean

clean: clean_gcell clean_aspect

.PHONY: all make_gcell make_aspect clean_gcell clean_aspect clean