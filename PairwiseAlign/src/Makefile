LDFLAGS = -lm
CFLAGS = -g -Wall -Wextra -std=c99
CC = gcc

# objects common to linear and affine algorithms
OBJCOMMON = aliCost.o aliGetSeq.o mem.o
# Smith-Waterman objects
OBJSW = aliCalcSW.o aliOutSW.o
# Altschul-Erickson objects
OBJAE = aliCalcAE.o aliOutAE.o

# first target is the default for "make", aliDnaLinear is the
# first one you should work on, you can add *Affine* targets
# from "all:" after implementing Altschul-Erickson
etu: aliDnaLinear

all: aliDnaLinear aliDnaAffine aliDnaAffineLongGaps aliProtAffine printStrings

aliDnaLinear: aliDnaLinear.o $(OBJCOMMON) $(OBJSW)
	gcc  $^  -o $@ $(LDFLAGS)

aliDnaAffine: aliDnaAffine.o $(OBJCOMMON) $(OBJAE)
	gcc  $^  -o $@ $(LDFLAGS)

aliDnaAffineLongGaps: aliDnaAffineLongGaps.o $(OBJCOMMON) $(OBJAE)
	gcc  $^  -o $@ $(LDFLAGS)

aliProtAffine: aliProtAffine.o $(OBJCOMMON) $(OBJAE)
	gcc  $^  -o $@ $(LDFLAGS)

printStrings: printStrings.o
	gcc  $^  -o $@ $(LDFLAGS)


%.o: %.c
	$(CC) $(CFLAGS) -MMD -c -o $@ $<

# deps made thanks to -MMD above, just include them
-include $(patsubst %.c,%.d,$(wildcard *.c))

clean:
	rm -f *.o *.d *~

