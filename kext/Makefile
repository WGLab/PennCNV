SHELL = /bin/sh
.SUFFIXES: .c .o

#note for Windows: replace khmm.so with khmm.dll below
#note for Mac OS X: replace -shared with -dynamiclib, replace khmm.so with khmm.dylib
#note for Solaris: replace gcc with cc below

CC = gcc
LD = gcc
CFLAGS = `perl -MExtUtils::Embed -e ccopts` -fPIC
LDFLAGS = `perl -MExtUtils::Embed -e ldopts`
ARCH = `perl -MConfig -e 'print $$Config{archname}'`
VER = `perl -MConfig -e 'print $$Config{version}'`

khmm.so: khmm_wrap.o khmm.o kc.o khmmDev.o
	$(LD) -shared -o $@ $^ $(LDFLAGS)
	mkdir -p $(VER)
	mkdir -p $(VER)/$(ARCH)/
	mkdir -p $(VER)/$(ARCH)/auto/
	mv $@ $(VER)/$(ARCH)/

clean:
	rm -f *.o *.so
