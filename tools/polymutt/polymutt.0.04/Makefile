TARGET = bin/polymutt
SUBDIRS=core src

# Build in all subdirectories.
.PHONY: $(SUBDIRS) $(TARGET) all clean help 


all: $(TARGET)
	(cd core; make)
	(cd src; make)
	@mkdir -p bin
	@cp src/polymutt $(TARGET)
	@echo "Compiling finished successfully.\n"
	@echo "Type ./$(TARGET) to invoke the command and see the command line options.\n"

noopenmp: $(TARGET)
	(cd core; make)
	(cd src; make -f Makefile.noopenmp)
	@mkdir -p bin
	@cp src/polymutt $(TARGET)
	@echo "Compiling finished successfully.\n"
	@echo "Type ./$(TARGET) to invoke the command and see the command line options.\n"

clean: 
	for i in $(SUBDIRS); do cd $$i && make clean && cd ..; done;
	@rm -f $(TARGET)

help : 
	@echo "Makefile help"
	@echo "-------------"
	@echo "Type...           To..."
	@echo "make              Compile everything "
	@echo "make help         Display this help screen"
	@echo "make clean        Delete temporary files"

test:
	(cd example; ./run.sh)
