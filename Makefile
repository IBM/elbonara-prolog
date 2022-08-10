




all: $(patsubst %.in,%.out,$(wildcard *.in))

%.out: %.in elbo.pl
	swipl $< > $@

