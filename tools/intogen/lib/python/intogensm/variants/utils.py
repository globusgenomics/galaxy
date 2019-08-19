__CB = {
	"A" : "T",
	"T" : "A",
	"G" : "C",
	"C" : "G"
}

def neutral_head_length(ref, alt):
	i = 0
	while i < len(ref) and ref[i] == "N" and i < len(alt) and alt[i] == "N":
		i += 1
	return i

def neutral_tail_length(ref, alt):
	i = len(ref) - 1
	j = len(alt) - 1
	while i >= 0 and ref[i] == "N" and j >= 0 and alt[j] == "N":
		i -= 1
		j -= 1
	return len(ref) - i - 1

def remove_neutral_caps(start, ref, alt):
	head_len = neutral_head_length(ref, alt)
	ins_correction = 1 if len(ref) < len(alt) else 0
	start += max(0, head_len - ins_correction)
	alt = alt[head_len:]
	ref = ref[head_len:]

	tail_len = neutral_tail_length(ref, alt)
	if tail_len > 0:
		alt = alt[:-tail_len]
		ref = ref[:-tail_len]

	return start, ref, alt

def prefix_length(ref, alt):
	i = 0
	while i < len(ref) and i < len(alt) and ref[i] == alt[i]:
		i += 1
	return i

def suffix_length(ref, alt):
	i = len(ref) - 1
	j = len(alt) - 1
	while i >= 0 and j >= 0 and ref[i] == alt[j]:
		i -= 1
		j -= 1
	return len(ref) - i - 1

def remove_common(start, ref, alt):
	prefix_len = prefix_length(ref, alt)
	ins_correction = 1 if len(ref) < len(alt) else 0
	start += max(0, prefix_len - ins_correction)
	alt = alt[prefix_len:]
	ref = ref[prefix_len:]

	suffix_len = suffix_length(ref, alt)
	if suffix_len > 0:
		alt = alt[:-suffix_len]
		ref = ref[:-suffix_len]

	return start, ref, alt

def complementary_sequence(seq):
	return "".join([__CB[base] if base in __CB else base for base in seq.upper()])

def var_to_tab(var, force_strand=None):

	start, ref, alt = remove_neutral_caps(var.start, var.ref, var.alt)

	start, ref, alt = remove_common(start, ref, alt)

	ref_len = len(ref)
	alt_len = len(alt)

	if alt_len > ref_len: # Insertion
		end = start
		start = start + 1
	elif alt_len < ref_len: # Deletion
		end = start + ref_len - 1
	else: # Substitution
		end = start + alt_len - 1

	if ref_len == 0:
		ref = "-"
	if alt_len == 0:
		alt = "-"

	if force_strand is not None and force_strand != var.strand:
		ref = complementary_sequence(ref)
		alt = complementary_sequence(alt)

	return start, end, ref, alt