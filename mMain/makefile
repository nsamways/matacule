CC = g++

OBJS = universe.o environmentArtifact.o chemistry.o agent.o bacteria.o genome.o genotype.o operon.o protein.o parameter.o

CFLAGS = -O3 -funroll-all-loops 

# CFLAGS = -g -Wall -O0 

RM = rm -f
CP = cp -r
program: ${OBJS}
	${CC} ${CFLAGS} -o matacule ${OBJS} main.cpp
	@echo
	@echo -- Done --
	@echo

universe.o: universe.cpp universe.h environmentArtifact.o chemistry.o protein.o parameter.o bacteria.o global.h 
	${CC} ${CFLAGS} -c universe.cpp

environmentArtifact.o: environmentArtifact.cpp environmentArtifact.h
	${CC} ${CFLAGS} -c environmentArtifact.cpp

chemistry.o: chemistry.cpp chemistry.h
	${CC} ${CFLAGS} -c chemistry.cpp

agent.o: agent.cpp agent.h
	${CC} ${CFLAGS} -c agent.cpp

bacteria.o: bacteria.cpp bacteria.h
	${CC} ${CFLAGS} -c bacteria.cpp

genome.o: genome.cpp genome.h
	${CC} ${CFLAGS} -c genome.cpp

genotype.o: genotype.cpp genotype.h
	${CC} ${CFLAGS} -c genotype.cpp

operon.o: operon.cpp operon.h
	${CC} ${CFLAGS} -c operon.cpp

protein.o: protein.cpp protein.h
	${CC} ${CFLAGS} -c protein.cpp

parameter.o: parameter.cpp parameter.h
	${CC} ${CFLAGS} -c parameter.cpp

clean:
	@echo -- cleaning objects --
	-$(RM) $(OBJS) *~
	-$(RM) matacule 
	@echo done

install:
	@echo -- installing --
	-$(CP) matacule ~/bin
	@echo done

uninstall:
	@echo -- uninstalling --
	-$(RM) ~/bin/matacule
