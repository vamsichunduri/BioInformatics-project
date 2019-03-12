# PROGRAM:    Project1
# PROGRAMMER: David Bierc
# LOGON ID:   z1737661
# DATE DUE:   9/18/2018
#

# Compiler variables
CCFLAGS = -ansi -Wall -std=c++11

# Rule to link object code files to create second executable file
Project_1: Project_1.o
	g++ $(CCFLAGS) -o Project_1 Project_1.o

# Rules to compile source code files to second object code
Project_1.o: Project_1.cpp
	g++ $(CCFLAGS) -c Project_1.cpp

#Execution command
run:
	./Project_1
# Pseudo-target to remove object code and executable files
clean:
	-rm *.o Project_1

