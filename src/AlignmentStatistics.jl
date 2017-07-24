module AlignmentStatistics

# package code goes here
using FastaIO, StatsBase, HypothesisTests, MultipleTesting,
DataFrames, Distributions, GaussDCA

export AAindex1_to_Dict,
read_AAindex1,
export_fasta,
get_sequence_labels,
split_fasta_sequences,
split_sequences,
clean_sequences,
sequence_composition,
sequences_to_numerical_properties,
PDB_AA_to_1_letter_code,
sequence_lengths,
symbols_in_sequences,
majority_consensus,
extract_label_element,
sorted_label_frequencies,
rectangularize_alignment
# binomial_CIs,
# compute_distances_dca_scores_table,
# dca_with_reference_sequence,
# Dirichlet_K3,
# Fisher_test_sequence_sets,
# PDB_to_single_chain_fasta,
# random_rows_from_msa,
# reference_sequence_column_labels,
# sequence_composition_difference,
          
PDB_AA_to_1_letter_code = Dict(
"ALA" => "A",
"ARG" => "R",
"ASN" => "N",
"ASP" => "D",
"CYS" => "C",
"GLN" => "Q",
"GLU" => "E",
"GLY" => "G",
"HIS" => "H",
"ILE" => "I",
"LEU" => "L",
"LYS" => "K",
"MET" => "M",
"PHE" => "F",
"PRO" => "P",
"SER" => "S",
"THR" => "T",
"TRP" => "W",
"TYR" => "Y",
"VAL" => "V"
)

"""
Input: fasta data read with FastaIO function readfasta

Output: sequences split up into single characters.
    
"""
function split_fasta_sequences(fasta::Array{Any,1})
    return map(i -> convert(Array{String,1},split(fasta[i][2],"")), 1:length(fasta))
end

"""
Can apply the following operations to sequences read with FastaIO.readfasta:

1. remove last char of each sequence (often "*"),

2. remove sequences that contain any of a set of given substrings, e.g. "X".

3. remove sequences that have a length different from a required length


Input:

- fasta array as read with FastaIO.readfasta

- optionally: should last character be removed of each sequence (default: remove_last_char=false); removing the last character is done first (if requested)

- optionally: remove sequences with any of the given strings (default: remove_seqs_with=["X", "\#", "*"])

- optionally: remove sequences with a length different from a required length (default: required_length = 0, i.e. option switched off)

    - Consider also using the rectangularize_alignment command as a way to produce alignment with homogeneous lengths
    
Output:

- Cleaned-up fasta (works with a copy of the original data)

"""
function clean_sequences(fasta::Array{Any,1};
                         remove_last_char::Bool=false,
                         remove_seqs_with::Array{String,1}=["X","#","*"],
                         required_length::Int64=0)
    n_seqs=length(fasta)
    new_fasta = copy(fasta)
    clean_seqs = fill(true, n_seqs)
    n_removal_triggers = length(remove_seqs_with)
    if remove_last_char
        for i in 1:n_seqs
            new_fasta[i] = (new_fasta[i][1], new_fasta[i][2][1:(end-1)])
        end
    end
    for i in 1:n_seqs
        for j in 1:n_removal_triggers
            if contains(new_fasta[i][2],remove_seqs_with[j]) 
                clean_seqs[i] = false
                break
            end
        end
    end
    if required_length != 0
        for i in 1:n_seqs
            for j in 1:n_removal_triggers
                if (clean_seqs[i] == true) &
                   (length(new_fasta[i][2])!=required_length)
                    clean_seqs[i] = false
                    break
                end
            end
        end
    end
    return new_fasta[clean_seqs]
end

# """
# Input:

# - fasta_msa: MSA as read with readfasta of FastaIO

# - n_select: number of rows to be selected

# Output:

# - MSA in fasta format consisting of n_select rows of fasta_msa
# """
# function random_rows_from_msa(fasta_msa::Array{Any,1}, n_select::Int64)
#     n_rows = length(fasta_msa)
#     if n_select >= n_rows
#         error("requested more random rows than available")
#     end
#     fasta_msa[sample(1:n_rows, n_select, replace=false)]
# end

"""
Input: array of strings (=sequences)

Output: sequences split up into single characters.

    Hint: if you have read the sequences from a fasta file with readfasta (package FastaIO), use split_fasta_sequences instead of split_sequences.

"""
function split_sequences(seqs::Any)
    return map(i -> convert(Array{String,1},split(seqs[i],"")), 1:length(seqs))
end
    
# """
#     Input: fasta data (in general multiple sequences) read with FastaIO function readfasta.

#     Output: vector of sequences lengths.
# """
# function fasta_lengths(fasta::Array{Any,1})
#     return map(i -> length(fasta[i][2]), 1:length(fasta))
# end
       
"""
Takes an alignment, as read by FastaIO.readall, and appends gap symbols '-' to shorter sequences to make all sequences the same length. Returns this rectangularized sequence array.
"""
function rectangularize_alignment(raw_ali::Array{Any,1})
    seqs = map(x -> convert(Array{String,1},split(x[2],"")), raw_ali[1:end])
    maxlen = maximum(map(x -> length(x), seqs))
    
    #Fill everything with gaps to make the array rectangular (otherwise we have dangling ends):
    Seqs = fill("-",length(seqs),maxlen);

    for i in 1:length(seqs)
        L = length(seqs[i])
        for j in 1:L
            Seqs[i,j] = seqs[i][j]
        end
    end

    return Seqs
end


"""
Returns sequence labels from sequences read with FastaIO.readall (without the leading '>').
"""
function get_sequence_labels(raw_ali::Array{Any,1})
    return convert(Array{String,1},map(x -> x[1], raw_ali[1:end]))
end


"""
Input: rectangularized sequence array.

Output: array of three components:

(1) consensus sequence where at each position the majority symbol in the respective alignment column is given;

(2) absolute frequency of majority symbols;

(3) relative frequency of majority symbols (absolute / n_rows).
"""
function majority_consensus(seqs::Array{String,2})

    n_rows, n_cols = size(seqs)

    maxelem = fill("-", n_cols)
    maxcount = zeros(Int64, n_cols)
                     
    for i in 1:n_cols
        coltab = countmap(seqs[:,i])
        vals = collect(values(coltab))
        maxval = maximum(vals)
        maxcount[i] = maxval
        maxelem[i] = collect(keys(coltab))[vals .== maxval][1]
    end

    return maxelem, maxcount, (convert(Array{Float64,1},maxcount) ./ n_rows)
end

"""
Input: sequence (ASCIIArray)

Output: dictionary of sequence composition
"""
function sequence_composition(seq::Array{String,1})
    return countmap(seq)
end

"""
Input:

(1) String array of sequence labels (eg from FASTA header)

(2) delimiter character

(3) number of element to be extracted

Output: array of extracted label elements (Strings)
"""
function extract_label_element(labels::Array{String,1}, delim::String, elno::Int64)
    return map(x -> convert(String,x[elno]), map(x -> split(x,delim), labels))
end


"""
# Input: array of labels.

# Output: 2D array with labels in first row and frequencies in second, sorted according to decreasing frequency.
# """
function sorted_label_frequencies(labels::Array{String, 1})
    label_counts = countmap(labels)
    p = sortperm(collect(values(label_counts)),rev=true)
    return permute!(collect(keys(label_counts)), p), permute!(collect(values(label_counts)), p)
end 


"""
Input:

(1) symbols: sequence of symbols (String) of the same length as each of the sequences in seqs.

(2) seqs: rectangular set of sequences (2D String array)

Output:

(1) absolute frequency of symbols in the respective columns of seqs.

(2) relative frequency of symbols in respective columns of seqs (absolute / n_rows).

"""
function symbols_in_sequences(symbols::Array{String,1}, seqs::Array{String,2})

    n_rows, n_cols = size(seqs)

    symbol_count = zeros(Int64, n_cols)
                     
    for i in 1:n_cols
        symbol_count[i] = count(x -> x == symbols[i], seqs[:,i])
    end

    return symbols, symbol_count, (convert(Array{Float64,1},symbol_count) ./ n_rows)

end

# """
# Input: two sequence sets, seqsA and seqsB, and a sequence of reference symbols, ref_syms.

# Optionally:
# - set parameter 'FDR' true to correct for multiple testing by control of false discovery rate (method of Benjamini and Hochberg); default value is false.
# - choose 'tail' of Fisher test between :left, :both, and :right

# Output: p-values from Fisher's exact test carried out for each column;
#     we test with 2x2 contingency table:
#     n_ref_sym_i(seqsAi), n_rows_A - n_ref_sym_i(seqsAi)
#     n_ref_sym_i(seqsBi), n_rows_B - n_ref_sym_i(seqsBi)    
# """
# function Fisher_test_sequence_sets(
#     seqsA::Array{String,2}, 
#     seqsB::Array{String,2},
#     ref_syms::Array{String,1};
#     FDR::Bool=false,
#     tail=:both
# )
#     rowsA, colsA = size(seqsA)
#     rowsB, colsB = size(seqsB)
#     len_ref = length(ref_syms)
    
#     if colsA != colsB || colsA != len_ref || colsB != len_ref
#         error("Fisher_test_sequence_sets: inconsistent lengths of sequences.")
#     end
    
#     pvals = zeros(len_ref)
    
#     #go through all columns and carry out test
#     for i in 1:len_ref
#         x = ref_syms[i]
#         nAx = count(y -> y == x, seqsA[:,i])
#         nBx = count(y -> y == x, seqsB[:,i])
#         if nAx + nBx != 0 && (nAx != rowsA || nBx != rowsB)
#             pvals[i] = pvalue(FisherExactTest(nAx, rowsA-nAx, nBx, rowsB-nBx), tail=tail)
#         else
#             pvals[i] = 1.0
#         end
#     end

#     if FDR==true #if user wishes FDR correction for multiple testing
#         return adjust(pvals, BenjaminiHochberg())
#     else #no correction
#         return pvals
#     end

# end

# """
#     Computes column labels of a multiple sequence alignment based on a reference sequence

#     Input:

#     - label of reference sequence (FASTA header without "")

#     - array of all labels (1D ASCII string array)

#     - sequence array (2D ASCII array)

#     Output: Labels for each column:

#     - columns in which the reference sequence has a letter are labelled with the number of this letter in the sequence (without gaps)

#     - columns in which the reference sequence has a gap are labelled with the number of the closest left neighbor letter and the number of gap symbols to the current column
# """

# function reference_sequence_column_labels(ref_label::String,
#                                           labels::Array{String,1},
#                                           seqs::Array{String,2})
#     ref_seq = vec(seqs[labels .== ref_label,:])
#     ali_len = length(ref_seq)
#     res_num = 0
#     gap_num = 1
#     ref_col = fill(" ", ali_len)
#     for i in 1:ali_len
#         if ref_seq[i] == "-"
#             ref_col[i] = join([string(res_num),".",string(gap_num)],"")
#             gap_num += 1
#         else
#             gap_num = 1
#             res_num += 1
#             ref_col[i] = string(res_num)
#         end 
#     end 
#     return ref_col
# end

"""
     Input:
     - filename for fasta export
     - label vector (1D String array)
     - sequence vector (2D String array)

     Output:
     writes a fasta file with given labels as headers
"""
function export_fasta(
                                filename::String,
                                labels::Array{String,1},
                                seqs::Array{String,2}
                                )
    n_seqs = length(labels)
    to_export = Array{Any}(n_seqs)
    for i in 1:n_seqs
        to_export[i] = (labels[i], join(collect(seqs[i,:]),""))
    end
    writefasta(filename, to_export)
end

"""
     Input:
     - filename for fasta export
     - label vector (1D String array)
     - sequence vector (array of String arrays)

     Output:
     writes a fasta file with given labels as headers
"""
function export_fasta(
                                filename::String,
                                labels::Array{String,1},
                                seqs::Array{String,1}
                                )
    n_seqs = length(labels)
    to_export = Array{Any}(n_seqs)
    for i in 1:n_seqs
        to_export[i] = (labels[i], seqs[i])
    end
    writefasta(filename, to_export)
end

# """
#     Input:
#     - array of successes x (Int64)
#     - array of trials n (Int64)

#     Output: three arrays
#     - fractions of successes
#     - lower boundaries of 95% confidence intervals (CIs)
#     - upper boundaries of 95% confidence intervals

#     If no x = n or x = 0, the CIs have a width of zero
# """
# function fractions_and_binomial_CIs(x::Array{Int64,1}, n::Array{Int64,1})

#     len = length(x)
#     if (len != length(n))
#         error("lengths of x and n do not match")
#     elseif sum(x .> n) > 0
#         error("x > n")
#     end
    
#     fraction = zeros(len)
#     lowerCI = zeros(len)
#     upperCI = zeros(len)
    
#     for j in 1:len
#         if ((!isnan(n[j])) & (n[j] != 0))
#             fraction[j] = x[j]/n[j]
#             if x[j] == n[j]
#                 (lowerCI[j],upperCI[j]) = (1.0,1.0)
#             elseif x[j] == 0
#                 (lowerCI[j],upperCI[j]) = (0.0,0.0)                
#             else 
#                 (lowerCI[j],upperCI[j]) = ci(BinomialTest(x[j],n[j]))
#             end
#         else
#             fraction[j] = NA
#         end 
#     end

#     return fraction, lowerCI, upperCI
# end

"""
   Input:
            
     - translation dictionary (AA -> numerical property)

     - sequences (array of arrays of strings)

   Output:
            
     - sequences translated to array of array of numerical properties
            
"""
function sequences_to_numerical_properties(trans::Dict, seqs::Array{Array{String,1},1})
    nseqs = length(seqs)
    props = Array{Any}(nseqs)
    for i in 1:nseqs
        props[i] = map(aa -> trans[aa], seqs[i])
    end
    return convert(Array{Array{Float64,1}},props)
end

"""
  Input:
          
     - AAindex1 data frame (see function read_AAindex1())

     - index to be used

    
  Output:

     - dictionary relating AAs with numbers (e.g. "A" => 3.95, ...)
"""                            
function AAindex1_to_Dict(AAindex1::DataFrame, ix::Int64)
    nAAs = 20 
    offset = 2
    AAs = map(x->string(x), names(AAindex1)[3:22])
    trans =
         Dict(AAs[i] => convert(Float64,AAindex1[ix, i+offset]) for i in 1:nAAs)
    return trans
end

"""
   Input: none (reads file AAindex1.csv from data directory of
         package AlignmentStatistics)

   Output: DataFrame of AAindex1
       
"""
function read_AAindex1()
    readtable(Pkg.dir("AlignmentStatistics","data","AAindex1.csv"))
end

# """
# Input:
         
#  - Three element Float64 array alpha, e.g. probabilities or counts of an AA to have charge -1, 0, 1.
#    All elements should be positive.

#  - Optional: minimum probability pmin, maximum probability pmax, number of steps np.

# Output:

#  - Array of probability densities computed for a Dirichlet distribution with concentration parameter
#      vector alpha, evaluated at (p_m1 = pmin, ..., pmax; p_0 = pmin, ..., pmax; p_1 = 1.0-p_m1-p_0).
#      This can e.g. correspond to the probabilities for the three charge states -1 (m1), 0 and 1.
         
# """
# function Dirichlet_K3(alpha::Array{Float64,1}; pmin::Float64=0.01, pmax::Float64=0.99, np::Int64=100)
#     if length(alpha) != 3
#         error("needs alpha of length 3")
#     elseif sum(alpha .<= 0.0) != 0
#         error("all elements of alpha must be > 0")
#     elseif (pmin <= 0.) | (pmax <= 0.) 
#         error("pmin and pmax must be positive")
#     elseif pmax <= pmin
#         error("pmax must be > pmin")
#     end
    
#     p_m1 = linspace(pmin,pmax,np)
#     p_0 = linspace(pmin,pmax,np)
#     pDir = zeros(np,np)
#     for i in 1:np 
#         pm1 = p_m1[i]
#         for j in 1:np
#             p0 = p_0[j]
#             if pm1 + p0 >= 1.
#                 break
#             end
#             p1 = 1.0 - pm1 - p0
#             pDir[i,j] = pdf(Dirichlet(alpha),[pm1;p0;p1])
#         end
#     end
#     return collect(p_m1), collect(p_0), pDir
# end

# AAs = ["A";"C";"D";"E";"F";"G";"H";"I";"K";"L";"M";"N";"P";"Q";"R";"S";"T";"V";"W";"Y"]
# DNAs = ["A";"C";"G";"T"]
# RNAs = ["A";"C";"G";"U"]

# """
# Input:
#          - string s (either String or Array{String,1})
#          - alphabet

# Output:
#          DataFrame with columns:
#          - alphabet: elements of the alphabet
#          - abs: absolute frequency of each element of the alphabet in s
#          - rel: relative frequency of each element of the alphabet in s

#          Note that the relative frequency does add up to 1.
#          The absolute frequency is lower than the number of characters in s if s
#              contains characters that are not elements of the alphabet.
# """
# function sequence_composition(s::Any, alphabet::Array{String,1})
#     if typeof(s)==String
#         a = split(s,"")
#     elseif typeof(s) != Array{String,1}
#         error("input has to be String or Array{String,1}")
#     end 
#     n = length(alphabet)
#     c = zeros(Int64,n)
#     i = 1
#     result = DataFrame(alphabet = alphabet , abs = zeros(Int64,n), rel = zeros(n))
#     for letter in alphabet
#         c[i] = count(x -> x==letter, a)
#         i+=1
#     end
#     result[:abs] = c
#     result[:rel] = c/sum(c)
#     result
# end

# """
# Input: two DataFrames as returned by function sequence_composition, one DataFrame
#     for each sequence s1 and s2

# Output: \sum_{i \in alphabet} |f_i(s_1)-f_i(s_2)| with relative frequencies f_i

# """
# function sequence_composition_difference(c1::DataFrame, c2::DataFrame)
#     if c1[:alphabet] != c2[:alphabet]
#         error("different alphabets")
#     end
#     sum(abs(c1[:rel] .- c2[:rel]))
# end 

"""
Sequence length distribution (as table) in a fasta array, read with FastaIO.readfasta.

Input:

- fasta array

Output:

- table of sequence lengths, generated with countmap
             
"""
function sequence_lengths(fasta::Array{Any,1})
    countmap(map(x->length(fasta[x][2]), 1:length(fasta)))
end

# """
# Input:

#          - pdb_input: pdb file name
         
#          - fasta_out_name: fasta output file name

#          - fasta header: header of fasta output file

# Output:

#          - fasta file with sequence extracted from pdb input and interpreted as single chain
# """
# function PDB_to_single_chain_fasta(
#                                    pdb_input::String,
#                                    fasta_out_name::String,
#                                    fasta_header::String
#                                    )
#     pdb = read(pdb_input,PDB)
#     CAs = collectatoms(pdb, calphaselector)
#     seq = (map(i -> PDB_AA_to_1_letter_code[CAs[i].res_name], 1:length(CAs)))'
#     export_fasta(fasta_out_name, [fasta_header], seq)
# end

# """
# Given a PDB of a protein (interpreted as single chain) and a MSA of sequences related to this protein, the function computes a table of distances and Direct Coupling Analysis (DCA) scores for all Calpha pairs.

# Input:

#     - pdb_file: name of PDB file

#     - msa_file: name of MSA file

#     - clean_up_msa=true (default): apply some clean-up on the MSA, e.g. give all columns the same lengths, remove rows (=sequences) with non-AA symbols. To switch this off, use clean_up_msa=false.

#     - remove_last_char=false (default): sometimes fasta blocks in MSA are finished with a "*" symbol. If remove_last_char=true, the last symbol is removed from all fasta blocks (one per block).

#     - dca_type="gDCA_FRN" (default): use package GaussDCA with Frobenius norm as score. Alternative are currently: "gDCA_DI" (direct information), and a file name from which an externally computed DCA output can be read (text array of tuples (i,j,score)) [the latter is currently not available].

#     - contact_threshold=8 (default): contact threshold in Angstrom
    
# Output:

#     - data frame with columns i, j, dij (=Calpha-Calpha distance), score (=DCA score), contact, true_positive_rate
# """
# function compute_distances_dca_scores_table(
#                 pdb_file::String, 
#                 msa_file::String;
#                 clean_up_msa::Bool=true,
#                 remove_last_char::Bool=false, #remove last char of each fasta block
#                 dca_type::String="gDCA_FRN",
#                 contact_threshold::Float64=8.0
#             )
    
#     #read pdb file and extract a fasta file -> write fasta
#     pdb_fasta_file = pdb_file*".fa"
#     PDB_to_single_chain_fasta(pdb_file, pdb_fasta_file, pdb_file)

#     #if requested: clean up MSA and write file with cleaned-up 
#     if clean_up_msa
#         seq = readfasta(msa_file)
#         clean_seq = clean_sequences(seq, remove_last_char=remove_last_char)
#         rect_seq = rectangularize_alignment(clean_seq)
#         msa_file = msa_file*"_clean.fa"
#         export_fasta(msa_file, get_sequence_labels(clean_seq), rect_seq)
#     end

#     #generate MSA from MSA file and pdb-fasta, with pdb-fasta as last entry
#     msa_out = (splitdir(pdb_file)[end])*(splitdir(msa_file)[end])
#     run(pipeline(`mafft --add $pdb_fasta_file $msa_file`, stdout=msa_out))

#     #compute DCA based on current MSA
#     dca_out = 0
#     if dca_type == "gDCA_FRN"
#         dca_out = gDCA(msa_out)
#     elseif dca_type == "gDCA_DI"
#         dca_out = gDCA(msa_out, pseudocount=0.2, score=:DI)
#     else # the DCA has been computed beforehand; dca_type is interpreted as name
#          # of file that contains tuples (i,j,score)
#         error("external DCA currently not available")
#         external = readdlm(dca_type,';')
#         n = length(external)
#         dca_out = Array{Any,1}(n)
#         for ix in 1:n
#             s = split(external[ix],['(',',',')'])
#             dca_out[ix]= (parse(s[2]),parse(s[3]),parse(s[4]))
#         end
#     end

#     #put DCA results in data frame for later join
#     n_dca = length(dca_out)
#     dca_df = DataFrame(
#         i_dca = map(i->dca_out[i][1],1:n_dca),
#         j_dca = map(i->dca_out[i][2],1:n_dca),
#         score = map(i->dca_out[i][3],1:n_dca)
#     )
    
#     #prepare computation of Calpha-Calpha distances
#     struct = read(pdb_file,PDB)
#     struct_ca = collectatoms(struct,calphaselector)
    
#     #find mapping of pdb-fasta to MSA columns necessary to map distances to DCA scores
#     struct_msa = readfasta(msa_out)
#     rect_struct_msa = rectangularize_alignment(struct_msa)

#          #the structure-sequence is the last sequence in the MSA
#     seq_struct_ali = collect(rect_struct_msa[end,:])
#     seq_struct = seq_struct_ali[find(x->x!="-", seq_struct_ali)]
#     indices = zeros(Int64,length(seq_struct))
#     j = 1
#     for i in 1:length(seq_struct_ali)
#         if (seq_struct_ali[i]!="-")
#             indices[j]=i
#             j += 1
#         end
#     end
#     #indices contains now the indices of the columns in which the pdb-fasta residues are
    
    
#     #generate data frame of Calpha-Calpha distances
#     n_ca = length(struct_ca)
#     n_pairs = round(Int64,n_ca*(n_ca-1)/2)

#     Dij = DataFrame()
#     Dij[:ix] = zeros(Int64,n_pairs)
#     Dij[:jx] = zeros(Int64,n_pairs)
#     Dij[:i_dca] = zeros(Int64,n_pairs)
#     Dij[:j_dca] = zeros(Int64,n_pairs)
#     Dij[:dij]=zeros(n_pairs)
#     k=1
#     for i in 1:(n_ca-1)
#         for j in (i+1):n_ca
#             dr = struct_ca[i].coords-struct_ca[j].coords
#             Dij[k,:dij] = sqrt(dot(dr,dr))
#             Dij[k,:ix] = i
#             Dij[k,:jx] = j
#             Dij[k,:i_dca] = indices[i]
#             Dij[k,:j_dca] = indices[j] 
#             k += 1
#         end
#     end
    
#     #join data frame with distances and DCA scores
#     Dij = join(Dij, dca_df, on=[:i_dca,:j_dca])

#     #sort according to decreasing DCA scores
#     Dij = sort(Dij, cols = [:score], rev=true)

#     #add contact columns
#     Dij[:contact] = (Dij[:dij] .<= contact_threshold)

#     #compute true positive rate up to distance pair i, j
#     Dij[:true_positive_rate] = cumsum(Dij[:contact]) ./ (1:length(Dij[:contact]))

#     #output
#     Dij
# end

# """
# Given an MSA of sequences of a protein and a reference sequence in that MSA given by the fasta header (without initial ">") the function computes DCA scores and provides the corresponding pair indices of the reference sequence in the output data frame.

# Input:

#     - msa_file: name of MSA file

#     - reference_header: fasta header of reference sequence (without initial ">")

#     - clean_up_msa=true (default): apply some clean-up on the MSA, e.g. give all columns the same lengths, remove rows (=sequences) with non-AA symbols. To switch this off, use clean_up_msa=false.

#     - remove_last_char=false (default): sometimes fasta blocks in MSA are finished with a "*" symbol. If remove_last_char=true, the last symbol is removed from all fasta blocks (one per block).

#     - dca_type="gDCA_FRN" (default): use package GaussDCA with Frobenius norm as score. Alternative are currently: "gDCA_DI" (direct information), and a file name from which an externally computed DCA output can be read (text array of tuples (i,j,score)) [the latter is currently not available].

#     - data frame with columns iref, jref, imsa, jmsa, score (=DCA score), 
# """
# function dca_with_reference_sequence(
#                 reference_label::String,                     
#                 msa_file::String;
#                 clean_up_msa::Bool=true,
#                 remove_last_char::Bool=false, #remove last char of each fasta block
#                 dca_type::String="gDCA_FRN"
#             )
    
#     #if requested: clean up MSA and write file with cleaned-up 
#     if clean_up_msa
#         seq = readfasta(msa_file)
#         clean_seq = clean_sequences(seq, remove_last_char=remove_last_char)
#         rect_seq = rectangularize_alignment(clean_seq)
#         msa_file = msa_file*"_clean.fa"
#         export_fasta(msa_file, get_sequence_labels(clean_seq), rect_seq)
#     end

#     #compute DCA based on MSA
#     dca_out = 0
#     if dca_type == "gDCA_FRN"
#         dca_out = gDCA(msa_file)
#     elseif dca_type == "gDCA_DI"
#         dca_out = gDCA(msa_file, pseudocount=0.2, score=:DI)
#     else # the DCA has been computed beforehand; dca_type is interpreted as name
#          # of file that contains tuples (i,j,score)
#         error("external DCA currently not available")
#         external = readdlm(dca_type,';')
#         n = length(external)
#         dca_out = Array{Any,1}(n)
#         for ix in 1:n
#             s = split(external[ix],['(',',',')'])
#             dca_out[ix]= (parse(s[2]),parse(s[3]),parse(s[4]))
#         end
#     end

#     #put DCA results in data frame for later join
#     n_dca = length(dca_out)
#     dca_df = DataFrame(
#         imsa = map(i->dca_out[i][1],1:n_dca),
#         jmsa = map(i->dca_out[i][2],1:n_dca),
#         score = map(i->dca_out[i][3],1:n_dca)
#     )
    
#     #rectangularize MSA and compute reference
#     msa = readfasta(msa_file)
#     rect_msa = rectangularize_alignment(msa)
#     labels = get_sequence_labels(msa)
#     ref_cols = reference_sequence_column_labels(reference_label, labels, rect_msa)

#     dca_df[:iref] = ref_cols[dca_df[:imsa]]
#     dca_df[:jref] = ref_cols[dca_df[:jmsa]]    
  
#     dca_df

# end

end # module










