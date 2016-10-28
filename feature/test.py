import khmer

# make a new ktable, L=6
ktable = khmer.new_ktable(6)

# count all k-mers in the given string
ktable.consume("ATGAGAGACACAGGGAGAGACCCAATTAGAGAATTGGACC")

# run through all entries. if they have nonzero presence, print.
for i in range(0, ktable.n_entries()):
   n = ktable.get(i)
   if n:
      print ktable.reverse_hash(i), "is present", n, "times."