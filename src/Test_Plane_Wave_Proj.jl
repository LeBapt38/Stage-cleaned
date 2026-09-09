## Testing the Julia implementation by comparing with the mathematica results
include("Plane_Wave_Proj.jl")

const ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
const DATA_DIR = joinpath(ROOT_DIR, "data")
const COEFF_PATH = joinpath(DATA_DIR, "SDPB_d=6_ci=0_lmax=16_Nmax=20_minmax=max_improv={TH, subTH, SPC}.txt")
const TEST_J_PATH = joinpath(DATA_DIR, "testJ.txt")

const TEST_D = 6
const TEST_NMAX = 20
const ERROR_COUNT_THRESHOLD = 1e-2
const ERROR_REPORT_THRESHOLD = 0.1
const MAX_REPORTED_ERRORS = 5

#%% Small Mathematica-input parser

function _strip_mathematica_precision_marks(s::AbstractString)
    # Examples handled here include 1.23`50, 1.23` and 1.23`50*^-20.
    return replace(s, r"(?<=[0-9.])`(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)?" => "")
end

function _replace_mathematica_numeric_heads(s::AbstractString)
    out = s
    # FullForm may print complex numbers as Complex[re, im].  Do this before
    # converting generic List[...] brackets.
    complex_pat = r"\bComplex\s*\[\s*([^,\[\]]+)\s*,\s*([^\[\]]+?)\s*\]"
    while occursin(complex_pat, out)
        out = replace(out, complex_pat => s"complex(\1, \2)")
    end

    rational_pat = r"\bRational\s*\[\s*([^,\[\]]+)\s*,\s*([^\[\]]+?)\s*\]"
    while occursin(rational_pat, out)
        out = replace(out, rational_pat => s"((\1) / (\2))")
    end
    return out
end

function _mathematica_to_julia_syntax(s::AbstractString)
    out = strip(s)
    out = _strip_mathematica_precision_marks(out)
    out = replace(out, "*^" => "e")
    out = _replace_mathematica_numeric_heads(out)
    out = replace(out, r"\b(?:List|Rule)\s*\[" => "Any[")
    out = replace(out, "\\[ImaginaryI]" => "im")
    out = replace(out, r"\bI\b" => "im")
    out = replace(out, "{" => "Any[")
    out = replace(out, "}" => "]")
    # Mathematica prints products as `3 I`; Julia needs `3*im` unless the
    # numeric-literal coefficient form is already adjacent.
    out = replace(out, r"(?<=[0-9\]])\s+im\b" => "*im")
    return out
end

function parse_mathematica_value(s::AbstractString)
    julia_syntax = _mathematica_to_julia_syntax(s)
    try
        return eval(Meta.parse(julia_syntax))
    catch err
        preview = first(julia_syntax, min(lastindex(julia_syntax), 500))
        throw(ErrorException("could not parse Mathematica record after conversion. Converted preview: $preview\nOriginal error: $err"))
    end
end

function _to_complex_number(x)
    x isa Complex && return x
    x isa Real && return complex(x)
    throw(ArgumentError("expected a real or complex number, got $(typeof(x)): $x"))
end

_is_real_or_complex_number(x) = x isa Real || x isa Complex

function _is_projection_leaf(row)
    row isa AbstractVector || return false
    # Current testJ format:
    # {Re[J], Im[J], Re[s], Im[s], Re[value], Im[value]}.
    return length(row) == 6 && all(_is_real_or_complex_number, row)
end

function _normalise_projection_row(row)
    if _is_projection_leaf(row)
        return (complex(row[1], row[2]), complex(row[3], row[4]), complex(row[5], row[6]))
    else
        throw(ArgumentError("unsupported Mathematica projection row format: $row"))
    end
end

function foreach_projection_row(data, f::Function)
    if _is_projection_leaf(data)
        f(_normalise_projection_row(data))
    elseif data isa AbstractVector
        for item in data
            foreach_projection_row(item, f)
        end
    else
        throw(ArgumentError("unsupported Mathematica projection data element: $data"))
    end
    return nothing
end

"""
    each_mathematica_record(path, f)

Stream top-level records from a Mathematica list without materialising the full
reference file as Julia data.  This is meant for exports of the form
`{{...}, {...}, ...}`.  A line-by-line fallback is used if no outer list records
are detected.
"""
function each_mathematica_record(path::AbstractString, f::Function)
    found_records = false
    open(path, "r") do io
        depth = 0
        in_record = false
        record = IOBuffer()
        while !eof(io)
            ch = read(io, Char)
            if ch == '{'
                depth += 1
                if depth == 2 && !in_record
                    in_record = true
                    truncate(record, 0)
                    seekstart(record)
                end
            end

            if in_record
                write(record, ch)
            end

            if ch == '}'
                if in_record && depth == 2
                    found_records = true
                    f(String(take!(record)))
                    in_record = false
                end
                depth -= 1
            end
        end
    end

    if !found_records
        open(path, "r") do io
            square_depth = 0
            token = IOBuffer()
            in_record = false
            record = IOBuffer()
            while !eof(io)
                ch = read(io, Char)

                if ch == '['
                    head = String(take!(token))
                    starts_record = head in ("List", "Rule")
                    square_depth += 1
                    if starts_record && square_depth == 2 && !in_record
                        found_records = true
                        in_record = true
                        truncate(record, 0)
                        seekstart(record)
                        write(record, head)
                        write(record, '[')
                        continue
                    end
                end

                if in_record
                    write(record, ch)
                end

                if ch == ']'
                    if in_record && square_depth == 2
                        f(String(take!(record)))
                        in_record = false
                    end
                    square_depth -= 1
                end

                if isletter(ch)
                    write(token, ch)
                elseif !isspace(ch)
                    truncate(token, 0)
                    seekstart(token)
                end
            end
        end
    end

    if !found_records
        for line in eachline(path)
            stripped = strip(line)
            isempty(stripped) && continue
            stripped in ("{", "}") && continue
            stripped = rstrip(stripped, [','])
            f(stripped)
        end
    end
    return nothing
end

#%% Coefficients

function _parse_mathematica_real_bigfloat(s::AbstractString)
    out = strip(_strip_mathematica_precision_marks(s))
    out = replace(out, "*^" => "e")
    out = replace(out, r"\s+" => "")
    if occursin('/', out)
        num, den = split(out, '/'; limit=2)
        return parse(BigFloat, num) / parse(BigFloat, den)
    end
    return parse(BigFloat, out)
end

function _parse_mathematica_scalar(s::AbstractString)
    # Preserve the high-precision SDPB coefficients instead of first converting
    # them through Julia Float64 literals.  Complex coefficients are not expected
    # for these alpha rules, but keep a fallback for possible future data files.
    if occursin("I", s) || occursin("ImaginaryI", s)
        return _to_complex_number(parse_mathematica_value(s))
    end
    return complex(_parse_mathematica_real_bigfloat(s))
end

function load_alpha_coefficients(path::AbstractString)
    # The Mathematica code uses mysdpbout[[2]], whose file begins as
    #   List[..., List[Rule[\[Alpha]0x0x0x0, value], ...]]
    # Accept both that full-form `Rule[symbol, value]` syntax and the infix
    # `symbol -> value` syntax Mathematica may emit in other exports.
    txt = read(path, String)
    coeffs = Dict{AlphaIndex, Complex{BigFloat}}()

    alpha_name = raw"(?:\\\[Alpha\]|α|Alpha|\\alpha|alpha)"
    int_pat = raw"(-?\d+)"
    # Coefficients in the SDPB file are scalar Mathematica numbers.  Stop the
    # value before the closing bracket/comma of Rule/List; precision marks and
    # *^ exponents are converted later by `_parse_mathematica_scalar`.
    value_pat = raw"([^,\]]+)"

    function store_match!(m, offset::Integer)
        idx = AlphaIndex(parse(Int, m.captures[offset]),
                         parse(Int, m.captures[offset + 1]),
                         parse(Int, m.captures[offset + 2]),
                         parse(Int, m.captures[offset + 3]))
        coeffs[idx] = _parse_mathematica_scalar(strip(m.captures[offset + 4]))
    end

    full_form_rule = Regex("Rule\\[\\s*" * alpha_name * int_pat * "x" * int_pat * "x" * int_pat * "x" * int_pat * "\\s*,\\s*" * value_pat * "\\]")
    infix_rule = Regex(alpha_name * int_pat * "x" * int_pat * "x" * int_pat * "x" * int_pat * "\\s*(?:->|:>)\\s*" * value_pat)

    for m in eachmatch(full_form_rule, txt)
        store_match!(m, 1)
    end
    for m in eachmatch(infix_rule, txt)
        store_match!(m, 1)
    end

    isempty(coeffs) && error("No alpha replacement rules were found in $path")
    return coeffs
end

#%% Projection reference rows

function load_projection_reference_rows(path::AbstractString)
    reference_rows = Tuple{Any, Any, Any}[]

    each_mathematica_record(path, function (raw_record)
        parsed_record = try
            parse_mathematica_value(raw_record)
        catch err
            preview = first(raw_record, min(lastindex(raw_record), 500))
            throw(ErrorException("failed to parse record in $path. Raw preview: $preview\nOriginal error: $err"))
        end

        try
            foreach_projection_row(parsed_record, function (normalised_row)
                push!(reference_rows, normalised_row)
            end)
        catch err
            preview = first(raw_record, min(lastindex(raw_record), 500))
            throw(ErrorException("failed to normalise rows in $path near row $(length(reference_rows) + 1). Raw preview: $preview\nOriginal error: $err"))
        end
    end)

    return reference_rows
end

#%% Error accounting and reporting

relative_error(mathematica_value, julia_value) = abs((mathematica_value - julia_value) / mathematica_value)

function compare_projection_reference_file(label::AbstractString, path::AbstractString, evaluator::Function;
                                           count_threshold=ERROR_COUNT_THRESHOLD,
                                           report_threshold=ERROR_REPORT_THRESHOLD,
                                           max_reported=MAX_REPORTED_ERRORS)
    reference_rows = load_projection_reference_rows(path)

    # Parsing is small compared with evaluating the high-precision expression;
    # collect rows first so the expensive independent evaluations can use all threads.
    evaluations = Vector{Any}(undef, length(reference_rows))
    Threads.@threads for row_index in eachindex(reference_rows)
        J, s, mathematica_value = reference_rows[row_index]
        evaluations[row_index] = try
            evaluator(J, s)
        catch err
            CapturedException(err, catch_backtrace())
        end
    end

    nbad_count = 0
    first_large_errors = NamedTuple[]

    for row_index in eachindex(reference_rows)
        J, s, mathematica_value = reference_rows[row_index]
        julia_value = evaluations[row_index]
        if julia_value isa CapturedException
            throw(ErrorException("failed to evaluate Julia projection for row $row_index in $path, J=$J, s=$s\n$(julia_value)"))
        end
        err = relative_error(mathematica_value, julia_value)

        if err > count_threshold
            nbad_count += 1
        end
        if err > report_threshold && length(first_large_errors) < max_reported
            push!(first_large_errors, (row=row_index,
                                       J=J,
                                       s=s,
                                       mathematica=mathematica_value,
                                       julia=julia_value,
                                       error=err))
        end
    end

    nrows = length(reference_rows)

    println()
    println(label)
    println("  rows compared: ", nrows)
    println("  number of errors above ", count_threshold, ": ", nbad_count)
    println("  first ", max_reported, " errors above ", report_threshold, ":")
    if isempty(first_large_errors)
        println("    none")
    else
        for item in first_large_errors
            println("    row = ", item.row)
            println("      J            = ", item.J)
            println("      s            = ", item.s)
            println("      mathematica  = ", item.mathematica)
            println("      julia        = ", item.julia)
            println("      error        = ", item.error)
        end
    end

    return (rows=nrows, errors_above_count_threshold=nbad_count, first_large_errors=first_large_errors)
end

function run_plane_wave_J_projection_comparison()
    coeffs = load_alpha_coefficients(COEFF_PATH)
    println("Loaded ", length(coeffs), " alpha coefficients from ", COEFF_PATH)
    model = make_projection_model(TEST_D, coeffs; Nmax=TEST_NMAX,
        numeric_type=:big)
    println("Precomputed ", length(model.terms), " shared amplitude/discontinuity terms")

    projection_summary = compare_projection_reference_file(
        "Froissart-Gribov / plane-wave projection exprJ[J, model, s]",
        TEST_J_PATH,
        (J, s) -> exprJ(J, model, s; p=8),
    )

    return (projection=projection_summary,)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_plane_wave_J_projection_comparison()
end