CC=g++

OBJDIR=./obj
SRCDIR=./
CURRENTDIR=$(shell pwd)

CFLAGS=-c -std=c++0x -O3 -Wall
LINKLAGS=-O3

#debug = true
ifdef debug
	CFLAGS +=-g
	LINKFLAGS += -flto
endif

INCDIR=\
	-I/usr/local/include \
	-I/usr/include \
	-I/usr/X11/include \
	-I./include \

LIBDIR=-L/usr/X11R6 -L/usr/local/lib -L./lib

LIBS=

OS_NAME:=$(shell uname -s)

ifeq ($(OS_NAME),Darwin)
	LIBS += \
		-framework Cocoa \
		-framework OpenGL \
		-framework IOKit \
		-framework CoreVideo
endif
ifeq ($(OS_NAME),Linux)
	LIBS += \
		-lGL \
		-lXi \
		-lglfw \
		-lGLEW
endif
ifeq ($(OS),Windows_NT)
	LIBS += \
		-lglfw3 \
		-lglew32 \
		-lopengl32 \
		-lgdi32		
endif
	
SOURCES=$(wildcard $(SRCDIR)/*cpp)
OBJECTS=$(addprefix $(OBJDIR)/,$(notdir $(SOURCES:.cpp=.o)))
EXECUTABLE=mapMaker.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LINKFLAGS) $(OBJECTS) -o $@ $(LIBS) $(LIBDIR)
	
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) $(INCDIR) $< -o $@ 
	
clean:
	rm $(OBJDIR)/*.o $(EXECUTABLE)
