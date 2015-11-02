module AlignmentStatistics

# package code goes here
using FastaIO, StatsBase, HypothesisTests, PValueAdjust

export rectangularize_ali,
          get_seq_labels,
          majority_consensus,
          sequence_composition,
          extract_label_element,
          sorted_label_frequencies,
          symbols_in_seqs,
          Fisher_test_sequence_sets
          
"""
Takes an alignment, as read by FastaIO.readall, and appends gap symbols '-' to shorter sequences to make all sequences the same length. Returns this rectangularized sequence array.
"""
function rectangularize_ali(raw_ali::Array{Any,1})
    seqs = map(x -> convert(Array{AbstractString,1},split(x[2],"")), raw_ali[1:end])
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
function get_seq_labels(raw_ali::Array{Any,1})
    return map(x -> x[1], raw_ali[1:end])
end


"""
Input: rectangularized sequence array.

Output: array of three components:

(1) consensus sequence where at each position the majority symbol in the respective alignment column is given;

(2) absolute frequency of majority symbols;

(3) relative frequency of majority symbols (absolute / n_rows).
"""
function majority_consensus(seqs::Array{ASCIIString,2})

    n_rows, n_cols = size(seqs)

    maxelem = fill("-", n_cols)
    maxcount = zeros(Int64, n_cols)
                     
    for i in 1:n_cols
        coltab = countmap(seqs[:,i])
        maxval = maximum(collect(values(coltab)))
        maxcount[i] = maxval
        maxelem[i] = coltab.keys[coltab.vals .== maxval][1]
    end

    return maxelem, maxcount, (convert(Array{Float64,1},maxcount) ./ n_rows)
end

"""
Input: sequence (ASCIIArray)

Output: dictionary of sequence composition
"""
function sequence_composition(seq::Array{ASCIIString,1})
    return countmap(seq)
end

"""
Input:

(1) ASCIIString array of sequence labels (eg from FASTA header)

(2) delimiter character

(3) number of element to be extracted

Output: array of extracted label elements (ASCIIStrings)
"""
function extract_label_element(labels::Array{Any,1}, delim::ASCIIString, elno::Int64)
    return map(x -> x[elno], map(x -> split(x,delim), labels))
end


"""
Input: array of labels.

Output: 2D array with labels in first row and frequencies in second, sorted according to decreasing frequency.
"""
function sorted_label_frequencies(labels::Array{Any, 1})
    label_counts = countmap(labels)
    p = sortperm(collect(values(label_counts)),rev=true)
    return permute!(collect(keys(label_counts)), p), permute!(collect(values(label_counts)), p)
end 


"""
Input:

(1) symbols: sequence of symbols (ASCIIString) of the same length as each of the sequences in seqs.

(2) seqs: rectangular set of sequences (2D ASCIIString array)

Output:

(1) absolute frequency of symbols in the respective columns of seqs.

(2) relative frequency of symbols in respective columns of seqs (absolute / n_rows).

"""
function symbols_in_seqs(symbols::Array{ASCIIString,1}, seqs::Array{ASCIIString,2})

    n_rows, n_cols = size(seqs)

    symbol_count = zeros(Int64, n_cols)
                     
    for i in 1:n_cols
        symbol_count[i] = count(x -> x == symbols[i], seqs[:,i])
    end

    return symbols, symbol_count, (convert(Array{Float64,1},symbol_count) ./ n_rows)

end

"""
Input: two sequence sets, seqsA and seqsB, and a sequence of reference symbols, ref_syms.

Optionally: set last parameter (FDR) true to correct for multiple testing by control of false discovery rate (method of Benjamini and Hochberg).

Output: p-values from Fisher's exact test carried out for each column;
    we test with 2x2 contingency table:
    n_ref_sym_i(seqsAi), n_rows_A - n_ref_sym_i(seqsAi)
    n_ref_sym_i(seqsBi), n_rows_B - n_ref_sym_i(seqsBi)    
"""
function Fisher_test_sequence_sets(
    seqsA::Array{ASCIIString,2}, 
    seqsB::Array{ASCIIString,2},
    ref_syms::Array{ASCIIString,1},
    FDR::Bool=false 
)
    rowsA, colsA = size(seqsA)
    rowsB, colsB = size(seqsB)
    len_ref = length(ref_syms)
    
    if colsA != colsB || colsA != len_ref || colsB != len_ref
        error("Fisher_test_sequence_sets: inconsistent lengths of sequences.")
    end
    
    pvals = zeros(len_ref)
    
    #go through all columns and carry out test
    for i in 1:len_ref
        x = ref_syms[i]
        nAx = count(y -> y == x, seqsA[:,i])
        nBx = count(y -> y == x, seqsB[:,i])
        if nAx + nBx != 0 && (nAx != rowsA || nBx != rowsB)
            pvals[i] = pvalue(FisherExactTest(nAx, rowsA-nAx, nBx, rowsB-nBx))
        else
            pvals[i] = 1.0
        end
    end

    if FDR==true #if user wishes FDR correction for multiple testing
        return PValueAdjust.padjust(pvals, BenjaminiHochberg)
    else #no correction
        return pvals
    end

end

end # module










