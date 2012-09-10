# Set compiler's options depending on the architecture
ifeq ($(ARCHI),i386)
  CC     = gcc
  CFLAGS = -O3 -g -c -Wuninitialized -Wunused -Winline -Wshadow \
           -fexpensive-optimizations -funroll-loops
# CFLAGS   = -g -c
  LDFLAGS= -g -static-libgcc
endif

# working dirs
EXEDIR = $(HOME)/bin/$(ARCHI)
SRCDIR = sources
OBJDIR = objects/$(ARCHI)
ARCDIR = archives
DIRDIR = objects $(EXEDIR) $(OBJDIR) $(ARCDIR)
INCDIR = 
LDLDIR = 
VPATH  = $(SRCDIR)

# objects list
src    = $(wildcard $(SRCDIR)/*.c)
header = $(wildcard $(SRCDIR)/*.h)
objs   = $(patsubst $(SRCDIR)%,$(OBJDIR)%,$(src:.c=.o))
prog   = mmg3d5

#.SILENT:

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(INCDIR) $(CFLAGS) $< -o $@

$(EXEDIR)/$(prog):$(DIRDIR) $(objs)
	echo "#define COMPIL " '"' `date` '"' > $(SRCDIR)/compil.date
	$(CC) $(CFLAGS) -c $(INCDIR) $(SRCDIR)/mmg3d.c -o $(OBJDIR)/mmg3d.o
	echo "$(CC) $(LDFLAGS) $(LDLDIR) $(objs) -o $@ -lm "
	$(CC) $(LDFLAGS) $(LDLDIR) $(objs) -o $@ -lm 

$(objs):$(header)

$(DIRDIR):
	@[ -d $@ ] || mkdir $@

clean:
	-rm $(objs) $(EXEDIR)/$(prog)

tar:$(DIRDIR)
	tar czfL $(ARCDIR)/$(prog).`date +"%Y.%m.%d"`.tgz sources makefile

target: $(EXEDIR)/$(prog)
