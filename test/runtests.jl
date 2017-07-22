using AlignmentStatistics
using Base.Test
using FastaIO

# write your own tests here
aaix = read_AAindex1() #test func.
@test size(aaix) == (534,22) #change the 534 accordingly if more indexes are contained

daaix = AAindex1_to_Dict(aaix,1) #test func.
@test daaix["A"] == 4.35
@test daaix["F"] == 4.66

fasta = readfasta("interleukin-4-aln.fa")
labs = get_sequence_labels(fasta) #test func.
@test labs[1] == "P07750"
@test labs[length(labs)] == "O77762"

run(`touch test.fa`)
run(`rm -f test.fa`)
rfasta = rectangularize_alignment(fasta) #test func.
@test rfasta[2,1] == "M"
@test rfasta[2,155] == "S"
@test rfasta[3,155] == "-"
export_fasta("test.fa",labs,rfasta) #test func.
fasta2 = readfasta("test.fa")
@test fasta2 == fasta

sfasta = split_fasta_sequences(fasta)
x = sequence_composition(sfasta[1][1:10]) #test func.
@test x["G"] == 1
@test x["V"] == 2
@test x["L"] == 2

@test PDB_AA_to_1_letter_code["ALA"] == "A" #test Dict
@test PDB_AA_to_1_letter_code["GLN"] == "Q"

seqs = readfasta("interleukin-4.fa")
seqs_only = split_fasta_sequences(seqs) 
nseqs = sequences_to_numerical_properties(daaix, seqs_only) #test func.
@test nseqs[1][1] == 4.52

x = sequence_lengths(seqs) #test func.
@test length(x) == 5
@test x[135] == 1

sis = symbols_in_sequences(sfasta[1], rfasta) #test func.
@test sis[2][1] == 5
@test sis[2][length(sis[2])] == 4
@test sis[3][length(sis[2])] == 0.8

mc = majority_consensus(rfasta) #test func.
@test mc[2][1] == 5
@test mc[2][length(mc[2])] == 4
@test mc[3][length(mc[2])] == 0.8

lfasta = readfasta("interleukin-4-labs.fa")
labs2 = extract_label_element(get_sequence_labels(lfasta),".",2) #test func.
slabs2 = sorted_label_frequencies(labs2) #test func.
@test slabs2[1][1] == "alpha"
@test slabs2[2][2] == 2







