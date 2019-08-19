#include "GeneticCode.h"

const char * codons[64] = {
	"AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG",
	"AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC",
	"CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA",
	"GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT",
	"GTA", "GTC", "GTG", "GTT", "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG",
	"TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT" };

const char * codons2aminoAcids [64] = {
	"Lys", "Asn", "Lys", "Asn", "Thr", "Thr", "Thr", "Thr", "Arg", "Ser", "Arg",
	"Ser", "Ile", "Ile", "Met", "Ile", "Gln", "His", "Gln", "His", "Pro", "Pro",
	"Pro", "Pro", "Arg", "Arg", "Arg", "Arg", "Leu", "Leu", "Leu", "Leu", "Glu",
	"Asp", "Glu", "Asp", "Ala", "Ala", "Ala", "Ala", "Gly", "Gly", "Gly", "Gly",
	"Val", "Val", "Val", "Val", "stop", "Tyr", "stop", "Tyr", "Ser", "Ser", "Ser",
	"Ser", "stop", "Cys", "Trp", "Cys", "Leu", "Phe", "Leu", "Phe" };

const char * codons2letters [64] = {
	"K", "N", "K", "N", "T", "T", "T", "T", "R", "S", "R", "S", "I", "I", "M",
	"I", "Q", "H", "Q", "H", "P", "P", "P", "P", "R", "R", "R", "R", "L", "L",
	"L", "L", "E", "D", "E", "D", "A", "A", "A", "A", "G", "G", "G", "G", "V",
	"V", "V", "V", "*", "Y", "*", "Y", "S", "S", "S", "S", "*", "C", "W", "C",
	"L", "F", "L", "F" };

const char * aminoAcids [21] = {
	"Ter", "Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu",
	"Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr" };

const char * aminoAcidLetters [21] = {
	"*", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q",
	"R", "S", "T", "V", "W", "Y" };


