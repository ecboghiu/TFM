# The compiler:
CC= gcc

# Flags:
# see https://stackoverflow.com/questions/3375697/useful-gcc-flags-for-c
CFLAGS= @compilerFlags
# C files:
SOURCES=  main.c stat.c adj.c syncperc.c nuc.c

# Name of output:
TARGET= $(name)

$(TARGET):
	$(CC) $(SOURCES) $(CFLAGS) -o $(TARGET).out

clean:
	rm -f $(TARGET)

run:
	./$(name).out
	#nohup ./$(name).out > $(name).txt &

valgrind:
	valgrind --leak-check=yes ./$(TARGET)
