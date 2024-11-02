"""
Pauli Operations
----------------
This module contains several useful operations on Pauli strings. The identity along with the three Pauli matrices are
interchangeably referred to as (1, 2, 3, 4) or (-, X, Y, Z).
"""

# These arrays are used to find products of Pauli matrices.
const rules = Int8[1 2 3 4;
    2 1 4 3;
    3 4 1 2;
    4 3 2 1]
const sign_rules = Complex{Int8}[1 1 1 1;
    1 1 1im -1im;
    1 -1im 1 1im;
    1 1im -1im 1]

function pauliprod(string1::Vector{Int8},
    string2::Vector{Int8})::Tuple{Vector{Int8},Complex{Int8},Bool}

    length(string1) == length(string2) || throw(DimensionMismatch("Strings have different lengths ($(length(string1)) and $(length(string2)))"))
    index = CartesianIndex.(string1, string2)
    result = rules[index]
    sign = prod(sign_rules[index])
    imag(sign) == 0 ? (return result, sign, true) : (return result, sign, false)
end

function pauliprod(string1::SubArray,
    string2::Vector{Int8})::Tuple{Vector{Int8},Complex{Int8},Bool}

    length(string1) == length(string2) || throw(DimensionMismatch("Strings have different lengths ($(length(string1)) and $(length(string2)))"))
    index = CartesianIndex.(string1, string2)
    result = rules[index]
    sign = prod(sign_rules[index])
    imag(sign) == 0 ? (return result, sign, true) : (return result, sign, false)
end

function pauliprod(string1::Vector{Int8},
    string2::SubArray)::Tuple{Vector{Int8},Complex{Int8},Bool}

    length(string1) == length(string2) || throw(DimensionMismatch("Strings have different lengths ($(length(string1)) and $(length(string2)))"))
    index = CartesianIndex.(string1, string2)
    result = rules[index]
    sign = prod(sign_rules[index])
    imag(sign) == 0 ? (return result, sign, true) : (return result, sign, false)
end

function pauliprod(string1::SubArray,
    string2::SubArray)::Tuple{Vector{Int8},Complex{Int8},Bool}

    length(string1) == length(string2) || throw(DimensionMismatch("Strings have different lengths ($(length(string1)) and $(length(string2)))"))
    index = CartesianIndex.(string1, string2)
    result = rules[index]
    sign = prod(sign_rules[index])
    imag(sign) == 0 ? (return result, sign, true) : (return result, sign, false)
end

function _paulisum(sentence1::Dict{Vector{Int8},Float64},
    sentence2::Dict{Vector{Int8},Float64})::Dict{Vector{Int8},Float64}

    result = copy(sentence1)
    for (key, value) in sentence2
        result[key] = get(result, key, 0.0) + value
    end
    return result
end

function paulisum(sentences::Dict{Vector{Int8},Float64}...;
    tol=0.0)::Dict{Vector{Int8},Float64}

    result = sentences[1]
    for sentence in sentences[2:end]
        result = _paulisum(result, sentence)
    end
    for (key, value) in result
        abs(value) <= tol && pop!(result, key)
    end
    return result
end

function _paulidiff(sentence1::Dict{Vector{Int8},Float64},
    sentence2::Dict{Vector{Int8},Float64})::Dict{Vector{Int8},Float64}

    result = copy(sentence1)
    for (key, value) in sentence2
        result[key] = get(result, key, 0.0) - value
    end
    return result
end

function paulidiff(sentences::Dict{Vector{Int8},Float64}...;
    tol=0.0)::Dict{Vector{Int8},Float64}

    result = sentences[1]
    for sentence in sentences[2:end]
        result = _paulidiff(result, sentence)
    end
    for (key, value) in result
        abs(value) <= tol && pop!(result, key)
    end
    return result
end

function _pauliprod(sentence1::Dict{Vector{Int8},Number},
    sentence2::Dict{Vector{Int8},Number})::Dict{Vector{Int8},Number}

    result = Dict{Vector{Int8},Number}()
    for (key1, value1) in sentence1
        for (key2, value2) in sentence2
            string, sign, _ = pauliprod(key1, key2)
            result[string] = get(result, string, 0) + sign * value1 * value2
        end
    end
    return result
end

function pauliprod(sentences::Dict{Vector{Int8},Number}...;
    tol=0.0)::Dict{Vector{Int8},Number}

    result = sentences[1]
    for sentence in sentences[2:end]
        result = _pauliprod(result, sentence)
    end
    for (key, value) in result
        abs(value) <= tol && pop!(result, key)
    end
    return result
end

function _conjugate(sentence::Dict{Vector{Int8},Float64},
    generator::Vector{Int8},
    angle::Float64;
    tol::Float64=0.0)::Dict{Vector{Int8},Float64}

    # trash = Matrix{Int8}(undef, length(generator), 0)
    # for (key, value) in sentence
    #     if !(key in eachcol(trash))
    #         string, sign, c = pauliprod(generator, key)
    #         trash = [trash string]
    #         if !c
    #             sentence[key] = cos(2 * angle) * value + imag(sign) * sin(2 * angle) * get(sentence, string, 0.0)
    #             sentence[string] = cos(2 * angle) * get(sentence, string, 0.0) - imag(sign) * sin(2 * angle) * value
    #             abs(sentence[string]) <= tol && pop!(sentence, string)
    #         end
    #         abs(sentence[key]) <= tol && pop!(sentence, key)
    #     end
    # end
    result = Dict{Vector{Int8},Float64}()
    for (key, value) in sentence
        string, sign, c = pauliprod(generator, key)
        if c
            result[key] = get(result, key, 0.0) + value
        else
            result[key] = get(result, key, 0.0) + value * cos(2 * angle)
            result[string] = get(result, string, 0.0) - imag(sign) * value * sin(2 * angle)
            abs(result[key]) <= tol && pop!(result, key)
            abs(result[string]) <= tol && pop!(result, string)
        end
    end
    return result
end

function _conjugate(sentence::Dict{Vector{Int8},Float64},
    generator::SubArray,
    angle::Float64;
    tol::Float64=0.0)::Dict{Vector{Int8},Float64}

    result = Dict{Vector{Int8},Float64}()
    for (key, value) in sentence
        string, sign, c = pauliprod(generator, key)
        if c
            result[key] = get(result, key, 0.0) + value
        else
            result[key] = get(result, key, 0.0) + value * cos(2 * angle)
            result[string] = get(result, string, 0.0) - imag(sign) * value * sin(2 * angle)
            abs(result[key]) <= tol && pop!(result, key)
            abs(result[string]) <= tol && pop!(result, string)
        end
    end
    return result
end

function conjugate(sentence::Dict{Vector{Int8},Float64},
    generators::Matrix{Int8},
    angles::Vector{Float64};
    tol=0.0)

    size(generators, 2) == length(angles) || throw(DimensionMismatch("Generators and angles need to have equal size ($(size(generators, 2)) and $(length(angles)))"))

    result = sentence
    for i in eachindex(angles)[end:-1:1]
    #     chunks = Iterators.partition(result, max(1, length(result) ÷ Threads.nthreads()))
    #     tasks = map(chunks) do chunk
    #         Threads.@spawn _conjugate(Dict(chunk), view(generators, :, i), angles[i], tol=tol)
    #     end
    #     result = paulisum(fetch.(tasks)..., tol=tol)
        result = _conjugate(result, generators[:, i], angles[i], tol=tol)
    end
    return result
end

function conjugate(sentence::Dict{Vector{Int8},Float64},
    generators::SubArray,
    angles::Vector{Float64};
    tol=0.0)

    size(generators, 2) == length(angles) || throw(DimensionMismatch("Generators and angles need to have equal size ($(size(generators, 2)) and $(length(angles)))"))

    result = sentence
    for i in eachindex(angles)[end:-1:1]
    #     chunks = Iterators.partition(result, max(1, length(result) ÷ Threads.nthreads()))
    #     tasks = map(chunks) do chunk
    #         Threads.@spawn _conjugate(Dict(chunk), view(generators, :, i), angles[i], tol=tol)
    #     end
    #     result = paulisum(fetch.(tasks)..., tol=tol)
        result = _conjugate(result, generators[:, i], angles[i], tol=tol)
    end
    return result
end


