CC=g++
ROOTCFLAGS=`root-config --cflags`
ROOTGLIBS=`root-config --glibs`
CFLAGS=-std=c++11 -g -Wall $(ROOTCFLAGS)

INCLDIR=./include
SRCDIR=./src
OBJDIR=./objs

CPPFLAGS=-I $(INCLDIR)
LDFLAGS=$(ROOTGLIBS)

SRC=$(wildcard $(SRCDIR)/*.cpp)
OBJS=$(SRC:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

DICT_PAGES=$(INCLDIR)/DataStructs.h $(INCLDIR)/LinkDef_sps.h
DICT=$(SRCDIR)/sps_dict.cxx
LIB=$(OBJDIR)/sps_dict.o

EXE=bin/wkill

.PHONY: all clean

all: $(EXE)

$(EXE): $(LIB) $(OBJS)
	$(CC) $^ -o $@ $(LDFLAGS)

$(LIB): $(DICT)
	$(CC) $(CFLAGS) $(CPPFLAGS) -I ./ -o $@ -c $^
	mv $(SRCDIR)/*.pcm ./bin/

$(DICT): $(DICT_PAGES)
	rootcling -f $@ $^

clean:
	$(RM) $(OBJS) $(EXE) $(LIB) $(DICT) ./bin/*.pcm

VPATH=$(SRCDIR)
$(OBJDIR)/%.o: %.cpp
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ -c $^

