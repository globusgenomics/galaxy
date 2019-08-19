def ascendant(ontology, rel="is_a"):
	asc_tree = {}
	for term in ontology.get_stanzas("term"):
		term_id = term.get_id()
		if term.contains_tag(rel):
			for isa in term.get_tag(rel):
				parent = isa.content
				if term_id not in asc_tree:
					asc_tree[term_id] = [parent]
				else:
					asc_tree[term_id] += [parent]
	return asc_tree

def descendant(ontology, rel="is_a"):
	des_tree = {}
	for term in ontology.get_stanzas("term"):
		term_id = term.get_id()
		if term.contains_tag(rel):
			for isa in term.get_tag(rel):
				parent = isa.content
				if parent not in des_tree:
					des_tree[parent] = [term_id]
				else:
					des_tree[parent] += [term_id]
	return des_tree

def all(ontology, rel="is_a"):
	asc_tree = {}
	des_tree = {}
	for term in ontology.get_stanzas("term"):
		term_id = term.get_id()
		if term.contains_tag(rel):
			for isa in term.get_tag(rel):
				parent = isa.content
				if term_id not in asc_tree:
					asc_tree[term_id] = [parent]
				else:
					asc_tree[term_id] += [parent]

				if parent not in des_tree:
					des_tree[parent] = [term_id]
				else:
					des_tree[parent] += [term_id]
	return (asc_tree, des_tree)
