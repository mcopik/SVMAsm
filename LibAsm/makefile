CC=nasm
LD=ld
SOURCES=foo.asm
OBJ=foo.o
AFLAGS=-f elf32 -g -F dwarf
LDFLAGS=-shared
OUTPUT=../bin/libasm.so
all:
	$(CC) $(SOURCES) $(AFLAGS) -o $(OBJ)
	$(LD) $(LDFLAGS) $(OBJ) -o $(OUTPUT)