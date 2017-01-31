CC = g++
CFLAGS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE
ORIG_CXXFLAGS = -Wall -I.

ifeq ($(BUILD), Debug)
CXXFLAGS = $(ORIG_CXXFLAGS) -D_DEBUG -g
else
CXXFLAGS = $(ORIG_CXXFLAGS) -DNDEBUG -O2
endif
ifeq ($(LINK), Static)
LDFLAGS = -static-libgcc
else
LDFLAGS = -shared-libgcc
endif

OBJ_FILES = $(SRC:%=%.o)
%.o:	%.c
	$(CC) $(CFLAGS) $(CXXFLAGS) -c $<

$(APP):	$(OBJ_FILES)

all:	$(APP)

clean:
	rm -rf core
	rm -rf *.o
	rm -rf $(APP)
