# ============================================================================
#  Generic Abaqus/CalculiX/Gmsh `.inp` reader used by the pre-processing step
#  (Post2D / Post3D).
#
#  Supported keywords (case-insensitive):
#    *NODE [NSET=]              2D (x,y) or 3D (x,y,z) coordinates
#    *ELEMENT TYPE= [ELSET=]    any element type, multi-line connectivity
#    *ELSET ELSET= [GENERATE] [INSTANCE=]   ids and/or names of other elsets
#    *NSET  NSET=  [GENERATE] [ELSET=] [INSTANCE=]
#    *SURFACE NAME= [TYPE=ELEMENT|NODE]     "elset_or_id, S1" lines
#    *SOLID/SHELL/MEMBRANE SECTION ELSET= MATERIAL=
#    *PART / *END PART, *ASSEMBLY, *INSTANCE (translation + rotation),
#    *END INSTANCE, *INCLUDE INPUT=
#  Everything else (materials, steps, outputs, ...) is skipped silently.
#
#  "Physical groups" are the named sets of the file: ELSET, NSET and SURFACE.
#  Gmsh writes every physical group as an *ELSET (and optionally *NSET),
#  Abaqus/CAE writes *Elset, *Nset and *Surface.
# ============================================================================

"""
    InpMesh

Mesh and named groups read from an Abaqus-style `.inp` file by [`read_inp`](@ref).

Fields
- `filepath`: absolute path of the file
- `nodes`: node id => `[x, y, z]`
- `elements`: element id => node ids (full connectivity, as in the file)
- `element_types`: element id => element type string (e.g. `"CPS4"`)
- `element_dims`: element id => topological dimension (1, 2 or 3)
- `elsets`: name => sorted element ids
- `nsets`: name => sorted node ids
- `surfaces`: name => `(element id, face number)` pairs (`face = 0`: the element itself)
- `sections`: elset name => material name (from `*... SECTION`)

Print the mesh (`println(mesh)`) to see all groups available in the file.
"""
struct InpMesh
    filepath::String
    nodes::Dict{Int,Vector{Float64}}
    elements::Dict{Int,Vector{Int}}
    element_types::Dict{Int,String}
    element_dims::Dict{Int,Int}
    elsets::Dict{String,Vector{Int}}
    nsets::Dict{String,Vector{Int}}
    surfaces::Dict{String,Vector{Tuple{Int,Int}}}
    sections::Dict{String,String}
end

# ----------------------------------------------------------------------------
#  Element type information
# ----------------------------------------------------------------------------

const _INP_1D_RE = r"^(T2D|T3D|B2|B3|PIPE|DC1D|R2D|RB2D|RAX|SAX|GK2D|FRAME)"
const _INP_3D_RE = r"^([A-Z]{0,3}C3D|SC)\d"
const _INP_2D_RE = r"^(CPS|CPE|CPEG|CAX|CGAX|DC2D|DCAX|DCC2D|DCCAX|AC2D|ACAX|S\d|STRI|DS\d|M3D|SFM3D|R3D|M2D|CPEG|SFMAX|MAX|MGAX)"
const _INP_NNODE_RE = r"(?:C3D|CPS|CPE|CAX|DC2D|DC3D|DCAX|AC2D|AC3D|T2D|T3D)(\d+)"

"""
    inp_element_dim(type) -> Int

Topological dimension of an Abaqus element type, 0 if unknown/unsupported
(springs, masses, connectors, ...).
"""
function inp_element_dim(type::AbstractString)
    t = uppercase(strip(type))
    occursin(_INP_1D_RE, t) && return 1
    occursin(_INP_3D_RE, t) && return 3
    occursin(_INP_2D_RE, t) && return 2
    return 0
end

# expected number of nodes of an element type, or 0 if it cannot be told from the name
function _inp_expected_nnodes(type::AbstractString)
    m = match(_INP_NNODE_RE, uppercase(type))
    m === nothing && return 0
    return parse(Int, m.captures[1])
end

# number of corner (vertex) nodes from dimension and number of nodes.
# Corner nodes always come first in Abaqus/Gmsh connectivity.
function _inp_ncorners(dim::Int, n::Int)
    if dim == 1
        n >= 2 && return 2
    elseif dim == 2
        n in (3, 6, 7) && return 3
        n in (4, 8, 9) && return 4
    elseif dim == 3
        n in (4, 10, 11) && return 4
        n in (5, 13, 14) && return 5
        n in (6, 15, 18) && return 6
        n in (8, 20, 27) && return 8
    end
    return 0
end

# Local faces (corner indices, cyclic order, Abaqus face numbering S1, S2, ...)
const _INP_FACES = Dict{Tuple{Int,Int},Vector{Vector{Int}}}(
    (2, 3) => [[1, 2], [2, 3], [3, 1]],
    (2, 4) => [[1, 2], [2, 3], [3, 4], [4, 1]],
    (3, 4) => [[1, 2, 3], [1, 4, 2], [2, 4, 3], [3, 4, 1]],
    (3, 5) => [[1, 2, 3, 4], [1, 5, 2], [2, 5, 3], [3, 5, 4], [4, 5, 1]],
    (3, 6) => [[1, 2, 3], [4, 6, 5], [1, 4, 5, 2], [2, 5, 6, 3], [3, 6, 4, 1]],
    (3, 8) => [[1, 2, 3, 4], [5, 8, 7, 6], [1, 5, 6, 2], [2, 6, 7, 3], [3, 7, 8, 4], [4, 8, 5, 1]],
)

# ----------------------------------------------------------------------------
#  Low level line handling
# ----------------------------------------------------------------------------

struct _InpLine
    text::String
    file::String
    lineno::Int
end

function _inp_collect_lines!(out::Vector{_InpLine}, filepath::String, depth::Int=0)
    depth > 20 && error("INP: *INCLUDE nested too deeply (recursive include?) in $filepath")
    isfile(filepath) || error("INP file does not exist: $filepath")
    for (lineno, raw) in enumerate(eachline(filepath))
        line = strip(raw)
        isempty(line) && continue
        startswith(line, "**") && continue
        if startswith(line, "*")
            kw, params = _inp_keyword(line)
            if kw == "INCLUDE"
                inc = get(params, "INPUT", "")
                isempty(inc) && error("INP: *INCLUDE without INPUT= ($filepath:$lineno)")
                incpath = isabspath(inc) ? inc : joinpath(dirname(filepath), inc)
                _inp_collect_lines!(out, incpath, depth + 1)
                continue
            end
        end
        push!(out, _InpLine(String(line), filepath, lineno))
    end
    return out
end

# "*Element, type=C3D8R, elset=Part-1" -> ("ELEMENT", Dict("TYPE"=>"C3D8R", "ELSET"=>"Part-1"))
function _inp_keyword(line::AbstractString)
    parts = split(line[2:end], ',')
    kw = uppercase(join(split(strip(parts[1])), " "))
    params = Dict{String,String}()
    for p in parts[2:end]
        p = strip(p)
        isempty(p) && continue
        if occursin('=', p)
            k, v = split(p, '='; limit=2)
            params[uppercase(strip(k))] = strip(strip(v), '"')
        else
            params[uppercase(p)] = ""
        end
    end
    return kw, params
end

_inp_tokens(line::AbstractString) = [String(strip(t)) for t in split(line, ',') if !isempty(strip(t))]

function _inp_int(tok::AbstractString, l::_InpLine)
    v = tryparse(Int, tok)
    v === nothing && (v = tryparse(Float64, tok); v !== nothing && isinteger(v) && (v = Int(v)))
    v isa Int || error("INP: expected an integer, got \"$tok\" ($(l.file):$(l.lineno))")
    return v
end

function _inp_float(tok::AbstractString, l::_InpLine)
    v = tryparse(Float64, replace(tok, r"[dD]" => "e"))
    v === nothing && error("INP: expected a number, got \"$tok\" ($(l.file):$(l.lineno))")
    return v
end

# ----------------------------------------------------------------------------
#  Parsing scopes (root / part / instance / assembly)
# ----------------------------------------------------------------------------

mutable struct _InpScope
    nodes::Dict{Int,Vector{Float64}}
    elements::Dict{Int,Vector{Int}}
    etypes::Dict{Int,String}
    elsets::Dict{String,Vector{Int}}
    nsets::Dict{String,Vector{Int}}
    surfaces::Dict{String,Vector{Tuple{Int,Int}}}
    raw_surfaces::Vector{Tuple{String,String,Vector{Tuple{String,String}}}} # (name, type, [(target, face)])
    sections::Dict{String,String}
end
_InpScope() = _InpScope(Dict(), Dict(), Dict(), Dict(), Dict(), Dict(), [], Dict())

mutable struct _InpInstance
    name::String
    part::String
    translation::Vector{Float64}
    rotation::Vector{Float64}   # empty or [ax, ay, az, bx, by, bz, angle_deg]
    scope::_InpScope            # mesh written directly inside *Instance (independent instance)
end

# case-insensitive lookup of a set name
function _inp_find(d::Dict{String,V}, name::AbstractString) where {V}
    haskey(d, name) && return d[name]
    u = uppercase(name)
    for (k, v) in d
        uppercase(k) == u && return v
    end
    return nothing
end

function _inp_addset!(d::Dict{String,Vector{Int}}, name::AbstractString, ids)
    key = name
    for k in keys(d)
        if uppercase(k) == uppercase(name)
            key = k
            break
        end
    end
    append!(get!(d, key, Int[]), ids)
    return d
end

# ----------------------------------------------------------------------------
#  Public reader
# ----------------------------------------------------------------------------

"""
    read_inp(filepath; verbose=true) -> InpMesh

Read an Abaqus / CalculiX / Gmsh `.inp` file: nodes, elements of any type, and all
named groups (`*ELSET`, `*NSET`, `*SURFACE`, section materials). Keywords are
case-insensitive, comments (`**`), continuation lines, `GENERATE`, `*INCLUDE` and
Abaqus `*PART`/`*INSTANCE` assemblies are supported.

```julia
mesh = read_inp("model.inp")
println(mesh)               # lists all physical groups of the file
geo  = Post2D(mesh; boundary=["Wall"])
```
"""
function read_inp(filepath::AbstractString; verbose::Bool=true)
    isempty(filepath) && error("File path cannot be empty!")
    path = abspath(filepath)
    lines = _inp_collect_lines!(_InpLine[], path)

    root = _InpScope()
    parts = Dict{String,_InpScope}()
    part_order = String[]
    instances = _InpInstance[]
    assembly = _InpScope()
    assembly_inst_sets = Tuple{Symbol,String,String,Vector{Int}}[] # (:elset/:nset, name, instance, ids)
    assembly_inst_surfs = Tuple{String,String,Vector{Tuple{String,String}},String}[]

    scope = root
    context = :root           # :root, :part, :assembly, :instance
    current_instance = nothing
    i = 1
    n = length(lines)
    while i <= n
        l = lines[i]
        if !startswith(l.text, "*")
            i += 1      # data line of an unsupported keyword
            continue
        end
        kw, params = _inp_keyword(l.text)
        # gather data lines of this keyword
        j = i + 1
        while j <= n && !startswith(lines[j].text, "*")
            j += 1
        end
        data = @view lines[i+1:j-1]

        if kw == "PART"
            name = get(params, "NAME", "PART-$(length(parts)+1)")
            scope = _InpScope()
            parts[name] = scope
            push!(part_order, name)
            context = :part
        elseif kw == "END PART"
            scope = root
            context = :root
        elseif kw == "ASSEMBLY"
            scope = assembly
            context = :assembly
        elseif kw == "END ASSEMBLY"
            scope = root
            context = :root
        elseif kw == "INSTANCE"
            name = get(params, "NAME", "INSTANCE-$(length(instances)+1)")
            part = get(params, "PART", "")
            tr = zeros(3)
            rot = Float64[]
            if length(data) >= 1
                v = [_inp_float(t, data[1]) for t in _inp_tokens(data[1].text)]
                tr[1:min(3, length(v))] .= v[1:min(3, length(v))]
            end
            if length(data) >= 2
                rot = [_inp_float(t, data[2]) for t in _inp_tokens(data[2].text)]
                length(rot) == 7 || error("INP: instance rotation needs 7 values ($(data[2].file):$(data[2].lineno))")
            end
            current_instance = _InpInstance(name, part, tr, rot, _InpScope())
            push!(instances, current_instance)
            scope = current_instance.scope
            context = :instance
        elseif kw == "END INSTANCE"
            scope = assembly
            context = :assembly
            current_instance = nothing
        elseif kw == "NODE"
            nset = get(params, "NSET", "")
            ids = Int[]
            for d in data
                tok = _inp_tokens(d.text)
                length(tok) >= 3 || error("INP: node line needs id and at least 2 coordinates ($(d.file):$(d.lineno)): $(d.text)")
                id = _inp_int(tok[1], d)
                xyz = zeros(3)
                for k in 1:min(3, length(tok) - 1)
                    xyz[k] = _inp_float(tok[k+1], d)
                end
                scope.nodes[id] = xyz
                push!(ids, id)
            end
            isempty(nset) || _inp_addset!(scope.nsets, nset, ids)
        elseif kw == "ELEMENT"
            etype = uppercase(get(params, "TYPE", ""))
            elset = get(params, "ELSET", "")
            expected = _inp_expected_nnodes(etype)
            ids = Int[]
            buf = String[]
            first_line = nothing
            for (k, d) in enumerate(data)
                first_line === nothing && (first_line = d)
                append!(buf, _inp_tokens(d.text))
                # expected node count known from the type: use it (trailing commas are then irrelevant),
                # otherwise a trailing comma marks a continuation line
                continued = expected > 0 ? (length(buf) - 1 < expected) : endswith(d.text, ",")
                if continued && k < length(data)
                    continue
                end
                length(buf) >= 2 || error("INP: element line without nodes ($(d.file):$(d.lineno))")
                id = _inp_int(buf[1], first_line)
                scope.elements[id] = [_inp_int(t, first_line) for t in buf[2:end]]
                scope.etypes[id] = etype
                push!(ids, id)
                empty!(buf)
                first_line = nothing
            end
            isempty(elset) || _inp_addset!(scope.elsets, elset, ids)
        elseif kw == "ELSET" || kw == "NSET"
            isel = kw == "ELSET"
            name = get(params, kw, "")
            isempty(name) && error("INP: *$kw without $kw= ($(l.file):$(l.lineno))")
            # the set may refer to an instance from the assembly level
            target = scope
            inst = get(params, "INSTANCE", "")
            ids = Int[]
            if !isel && haskey(params, "ELSET")
                src = _inp_find(scope.elsets, params["ELSET"])
                src === nothing && error("INP: *NSET refers to unknown elset \"$(params["ELSET"])\" ($(l.file):$(l.lineno))")
                for e in src
                    append!(ids, scope.elements[e])
                end
            end
            for d in data
                tok = _inp_tokens(d.text)
                if haskey(params, "GENERATE")
                    a = _inp_int(tok[1], d)
                    b = length(tok) >= 2 ? _inp_int(tok[2], d) : a
                    s = length(tok) >= 3 ? _inp_int(tok[3], d) : 1
                    append!(ids, a:s:b)
                else
                    for t in tok
                        v = tryparse(Int, t)
                        if v !== nothing
                            push!(ids, v)
                        else
                            ref = _inp_find(isel ? scope.elsets : scope.nsets, t)
                            ref === nothing && error("INP: *$kw $name refers to unknown set \"$t\" ($(d.file):$(d.lineno))")
                            append!(ids, ref)
                        end
                    end
                end
            end
            if context == :assembly && !isempty(inst)
                push!(assembly_inst_sets, (isel ? :elset : :nset, name, inst, ids))
            else
                _inp_addset!(isel ? target.elsets : target.nsets, name, ids)
            end
        elseif kw == "SURFACE"
            name = get(params, "NAME", "")
            isempty(name) && error("INP: *SURFACE without NAME= ($(l.file):$(l.lineno))")
            stype = uppercase(get(params, "TYPE", "ELEMENT"))
            entries = Tuple{String,String}[]
            for d in data
                tok = _inp_tokens(d.text)
                isempty(tok) && continue
                push!(entries, (tok[1], length(tok) >= 2 ? uppercase(tok[2]) : ""))
            end
            push!(scope.raw_surfaces, (name, stype, entries))
        elseif endswith(kw, " SECTION")
            elset = get(params, "ELSET", "")
            mat = get(params, "MATERIAL", "")
            if !isempty(elset) && !isempty(mat)
                scope.sections[elset] = mat
            end
        end
        i = j
    end

    # ---------------- flatten everything into one global mesh ----------------
    nodes = Dict{Int,Vector{Float64}}()
    elements = Dict{Int,Vector{Int}}()
    etypes = Dict{Int,String}()
    elsets = Dict{String,Vector{Int}}()
    nsets = Dict{String,Vector{Int}}()
    surfaces = Dict{String,Vector{Tuple{Int,Int}}}()
    sections = Dict{String,String}()

    _resolve_surfaces!(root)
    _merge_scope!(nodes, elements, etypes, elsets, nsets, surfaces, sections, root, 0, 0, "", identity)

    if isempty(instances) && !isempty(parts)
        # parts without an assembly: use every part once, untransformed
        for name in part_order
            push!(instances, _InpInstance(name, name, zeros(3), Float64[], _InpScope()))
        end
    end
    if !isempty(instances)
        node_off = isempty(nodes) ? 0 : maximum(keys(nodes))
        elem_off = isempty(elements) ? 0 : maximum(keys(elements))
        inst_offsets = Dict{String,Tuple{Int,Int}}()
        single = length(instances) == 1
        for inst in instances
            src = haskey(parts, inst.part) ? parts[inst.part] : _inp_find(parts, inst.part)
            if src === nothing || (isempty(src.nodes) && !isempty(inst.scope.nodes))
                src = inst.scope
            end
            _resolve_surfaces!(src)
            tf = _inp_instance_transform(inst)
            prefix = single ? "" : inst.name * "."
            _merge_scope!(nodes, elements, etypes, elsets, nsets, surfaces, sections, src,
                          node_off, elem_off, prefix, tf)
            if single
                # also keep the qualified name, as Abaqus/CAE refers to it that way
                _merge_scope!(nothing, nothing, nothing, elsets, nsets, surfaces, nothing, src,
                              node_off, elem_off, inst.name * ".", tf)
            end
            inst_offsets[uppercase(inst.name)] = (node_off, elem_off)
            node_off = isempty(nodes) ? node_off : max(node_off, maximum(keys(nodes)))
            elem_off = isempty(elements) ? elem_off : max(elem_off, maximum(keys(elements)))
        end
        for (kind, name, inst, ids) in assembly_inst_sets
            off = get(inst_offsets, uppercase(inst), nothing)
            off === nothing && error("INP: set \"$name\" refers to unknown instance \"$inst\"")
            o = kind == :elset ? off[2] : off[1]
            _inp_addset!(kind == :elset ? elsets : nsets, name, ids .+ o)
        end
        # assembly level surfaces may reference instance-qualified elsets
        assembly.elsets = merge(elsets, assembly.elsets)
        assembly.elements = elements
        _resolve_surfaces!(assembly)
        _merge_scope!(nothing, nothing, nothing, elsets, nsets, surfaces, sections, assembly, 0, 0, "", identity;
                      skip_elsets=true)
    else
        _resolve_surfaces!(assembly)
        _merge_scope!(nodes, elements, etypes, elsets, nsets, surfaces, sections, assembly, 0, 0, "", identity)
    end

    for v in values(elsets); unique!(sort!(v)); end
    for v in values(nsets); unique!(sort!(v)); end
    for v in values(surfaces); unique!(sort!(v)); end

    # element dimensions + sanity checks
    dims = Dict{Int,Int}()
    unknown = Set{String}()
    for (id, conn) in elements
        d = inp_element_dim(etypes[id])
        d == 0 && push!(unknown, isempty(etypes[id]) ? "<no TYPE>" : etypes[id])
        dims[id] = d
        for nid in conn
            haskey(nodes, nid) || error("INP: element $id ($(etypes[id])) uses undefined node $nid")
        end
    end
    if !isempty(unknown) && verbose
        @warn "INP: element types not used for PD discretisation (unknown/unsupported): $(join(sort!(collect(unknown)), ", "))"
    end
    isempty(nodes) && error("INP: no nodes found in $path")
    isempty(elements) && error("INP: no elements found in $path")

    mesh = InpMesh(path, nodes, elements, etypes, dims, elsets, nsets, surfaces, sections)
    verbose && println("Read INP file: $path ($(length(nodes)) nodes, $(length(elements)) elements, ",
                       "$(length(elsets)) elsets, $(length(nsets)) nsets, $(length(surfaces)) surfaces)")
    return mesh
end

function _inp_instance_transform(inst::_InpInstance)
    tr = inst.translation
    rot = inst.rotation
    if isempty(rot) || rot[7] == 0
        all(iszero, tr) && return identity
        return x -> x .+ tr
    end
    a = rot[1:3]
    b = rot[4:6]
    k = b .- a
    k ./= norm(k)
    θ = deg2rad(rot[7])
    # Abaqus: translation first, then rotation about the axis a->b (a, b given after translation)
    return function (x)
        y = x .+ tr .- a
        y = y .* cos(θ) .+ cross(k, y) .* sin(θ) .+ k .* dot(k, y) .* (1 - cos(θ))
        return y .+ a
    end
end

function _resolve_surfaces!(s::_InpScope)
    for (name, stype, entries) in s.raw_surfaces
        if stype == "NODE"
            ids = Int[]
            for (t, _) in entries
                v = tryparse(Int, t)
                if v !== nothing
                    push!(ids, v)
                else
                    ref = _inp_find(s.nsets, t)
                    ref === nothing && error("INP: *SURFACE $name refers to unknown nset \"$t\"")
                    append!(ids, ref)
                end
            end
            _inp_addset!(s.nsets, name, ids)
        else
            faces = get!(s.surfaces, name, Tuple{Int,Int}[])
            for (t, f) in entries
                m = match(r"^S(\d+)$", f)
                face = m === nothing ? 0 : parse(Int, m.captures[1])
                v = tryparse(Int, t)
                ids = v !== nothing ? [v] : _inp_find(s.elsets, t)
                ids === nothing && error("INP: *SURFACE $name refers to unknown elset \"$t\"")
                for e in ids
                    push!(faces, (e, face))
                end
            end
        end
    end
    empty!(s.raw_surfaces)
    return s
end

function _merge_scope!(nodes, elements, etypes, elsets, nsets, surfaces, sections, s::_InpScope,
                       noff::Int, eoff::Int, prefix::String, tf; skip_elsets::Bool=false)
    if nodes !== nothing
        for (id, x) in s.nodes
            haskey(nodes, id + noff) && error("INP: duplicate node id $(id + noff)")
            nodes[id+noff] = tf(x)
        end
    end
    if elements !== nothing
        for (id, conn) in s.elements
            haskey(elements, id + eoff) && error("INP: duplicate element id $(id + eoff)")
            elements[id+eoff] = conn .+ noff
            etypes[id+eoff] = s.etypes[id]
        end
    end
    if !skip_elsets
        for (k, v) in s.elsets
            _inp_addset!(elsets, prefix * k, v .+ eoff)
        end
    end
    for (k, v) in s.nsets
        _inp_addset!(nsets, prefix * k, v .+ noff)
    end
    for (k, v) in s.surfaces
        append!(get!(surfaces, prefix * k, Tuple{Int,Int}[]), [(e + eoff, f) for (e, f) in v])
    end
    if sections !== nothing
        for (k, v) in s.sections
            sections[prefix*k] = v
        end
    end
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", m::InpMesh)
    dimcount = Dict{Int,Int}()
    for d in values(m.element_dims)
        dimcount[d] = get(dimcount, d, 0) + 1
    end
    println(io, "InpMesh: ", m.filepath)
    println(io, "  nodes: ", length(m.nodes), ", elements: ", length(m.elements),
            " (", join(["$(v) × $(k)D" for (k, v) in sort!(collect(dimcount))], ", "), ")")
    if !isempty(m.elsets)
        println(io, "  element sets (ELSET):")
        for name in sort!(collect(keys(m.elsets)))
            ids = m.elsets[name]
            ds = sort!(unique([get(m.element_dims, e, -1) for e in ids]))
            mat = _inp_find(m.sections, name)
            println(io, "    ", rpad(name, 24), lpad(length(ids), 8), " elements, dim ", join(ds, "/"),
                    mat === nothing ? "" : ", material = $mat")
        end
    end
    if !isempty(m.nsets)
        println(io, "  node sets (NSET):")
        for name in sort!(collect(keys(m.nsets)))
            println(io, "    ", rpad(name, 24), lpad(length(m.nsets[name]), 8), " nodes")
        end
    end
    if !isempty(m.surfaces)
        println(io, "  surfaces (SURFACE):")
        for name in sort!(collect(keys(m.surfaces)))
            println(io, "    ", rpad(name, 24), lpad(length(m.surfaces[name]), 8), " faces")
        end
    end
end
Base.show(io::IO, m::InpMesh) = print(io, "InpMesh(", length(m.nodes), " nodes, ", length(m.elements), " elements)")

# ----------------------------------------------------------------------------
#  Geometry helpers shared by Post2D / Post3D
# ----------------------------------------------------------------------------

_corners(m::InpMesh, e::Int) = (c = m.elements[e]; c[1:_inp_ncorners(m.element_dims[e], length(c))])

function _group_names_message(m::InpMesh)
    names = String[]
    append!(names, sort!(collect(keys(m.elsets))))
    append!(names, sort!(collect(keys(m.nsets))))
    append!(names, sort!(collect(keys(m.surfaces))))
    return isempty(names) ? "(the file defines no named groups)" : join(unique(names), ", ")
end

"""
    domain_elements(mesh, dim; domain=nothing) -> Vector{Int}

Sorted ids of the elements that become PD points: all elements of dimension `dim`,
or only those in the element sets listed in `domain`. The PD point `i` corresponds
to element `domain_elements(...)[i]`.
"""
function domain_elements(m::InpMesh, dim::Int; domain=nothing)
    if domain === nothing
        ids = [e for (e, d) in m.element_dims if d == dim]
    else
        names = domain isa AbstractString || domain isa Symbol ? [string(domain)] : string.(collect(domain))
        ids = Int[]
        for name in names
            set = _inp_find(m.elsets, name)
            set === nothing && error("Domain group \"$name\" is not an element set. Available groups: $(_group_names_message(m))")
            sel = [e for e in set if get(m.element_dims, e, 0) == dim]
            isempty(sel) && error("Domain group \"$name\" contains no $(dim)D elements.")
            append!(ids, sel)
        end
    end
    unique!(sort!(ids))
    if isempty(ids)
        found = sort!(unique(collect(values(m.element_dims))))
        error("No $(dim)D elements found in $(m.filepath) (element dimensions present: $(join(found, ", "))). " *
              "Use Post2D for 2D meshes and Post3D for 3D meshes.")
    end
    for e in ids
        _inp_ncorners(dim, length(m.elements[e])) == 0 &&
            error("Element $e ($(m.element_types[e]), $(length(m.elements[e])) nodes) is not supported for the PD discretisation.")
    end
    return ids
end

# all faces (corner node ids, cyclic) of a domain element
function _element_faces(m::InpMesh, e::Int)
    c = _corners(m, e)
    tbl = get(_INP_FACES, (m.element_dims[e], length(c)), nothing)
    tbl === nothing && error("Element $e ($(m.element_types[e])) has no face table.")
    return [c[f] for f in tbl]
end

# face key -> list of (element, local face index) for all domain elements
function _face_map(m::InpMesh, dom::Vector{Int})
    fm = Dict{Vector{Int},Vector{Tuple{Int,Int}}}()
    for e in dom
        for (k, f) in enumerate(_element_faces(m, e))
            push!(get!(fm, sort(f), Tuple{Int,Int}[]), (e, k))
        end
    end
    return fm
end

"""
Resolve a list of boundary group names into `(element, face corner nodes)` pairs of
domain elements. A group may be
- an ELSET of lower-dimensional elements (Gmsh physical lines/surfaces),
- an ELSET of domain elements (their free faces are used),
- a SURFACE (`elset, S1` element faces),
- an NSET (free faces whose corner nodes all lie in the set),
- `:all` (every free face of the domain).
"""
# Topology of the domain, built once and reused for every boundary group
struct _FaceCtx
    domset::Set{Int}
    fm::Dict{Vector{Int},Vector{Tuple{Int,Int}}}          # sorted face nodes -> owners
    node2dom::Dict{Int,Vector{Int}}                        # corner node -> domain elements
    free::Vector{Tuple{Vector{Int},Int,Int}}               # free faces: (sorted nodes, element, local face)
end

function _face_ctx(m::InpMesh, dom::Vector{Int})
    fm = _face_map(m, dom)
    node2dom = Dict{Int,Vector{Int}}()
    for e in dom, nid in _corners(m, e)
        push!(get!(node2dom, nid, Int[]), e)
    end
    free = [(key, owners[1][1], owners[1][2]) for (key, owners) in fm if length(owners) == 1]
    return _FaceCtx(Set(dom), fm, node2dom, free)
end

function _boundary_faces(m::InpMesh, dom::Vector{Int}, groups, dim::Int; warn::Bool=true,
                         ctx::_FaceCtx=_face_ctx(m, dom))
    domset = ctx.domset
    fm = ctx.fm
    node2dom = ctx.node2dom
    result = Dict{Tuple{Int,Vector{Int}},Nothing}()   # ordered by insertion is not needed
    out = Tuple{Int,Vector{Int}}[]
    add!(e, f) = (key = (e, sort(f)); haskey(result, key) || (result[key] = nothing; push!(out, (e, f))))
    free(f) = length(fm[sort(f)]) == 1

    # lower-dimensional element: add every domain element that contains all its corners
    function add_lower!(be::Int, name)
        c = _corners(m, be)
        isempty(c) && return 0
        cands = get(node2dom, c[1], Int[])
        cnt = 0
        for e in cands
            ce = _corners(m, e)
            if all(in(ce), c)
                for f in _element_faces(m, e)
                    if sort(f) == sort(c)
                        add!(e, f); cnt += 1
                    end
                end
            end
        end
        return cnt
    end

    if groups === :all || groups == "all" || groups == [:all]
        for (_, e, k) in ctx.free
            add!(e, _element_faces(m, e)[k])
        end
        return sort!(out; by=x -> (x[1], sort(x[2])))
    end

    names = groups isa AbstractString || groups isa Symbol ? [string(groups)] : string.(collect(groups))
    for name in names
        nadded = length(out)
        if (set = _inp_find(m.elsets, name)) !== nothing
            for be in set
                d = get(m.element_dims, be, 0)
                if d == dim - 1
                    add_lower!(be, name)
                elseif d == dim && be in domset
                    for f in _element_faces(m, be)
                        free(f) && add!(be, f)
                    end
                end
            end
        elseif (surf = _inp_find(m.surfaces, name)) !== nothing
            for (e, k) in surf
                d = get(m.element_dims, e, 0)
                if k == 0 && d == dim - 1
                    add_lower!(e, name)
                elseif k > 0 && e in domset
                    faces = _element_faces(m, e)
                    k <= length(faces) || error("Surface \"$name\": element $e has no face S$k")
                    add!(e, faces[k])
                end
            end
        elseif (ns = _inp_find(m.nsets, name)) !== nothing
            nodeset = Set(ns)
            for (key, e, k) in ctx.free
                all(in(nodeset), key) || continue
                add!(e, _element_faces(m, e)[k])
            end
        else
            error("Boundary group \"$name\" not found in $(m.filepath). Available groups: $(_group_names_message(m))")
        end
        warn && length(out) == nadded && @warn "Boundary group \"$name\" does not touch any $(dim)D domain element."
    end
    return sort!(out; by=x -> (x[1], sort(x[2])))   # deterministic order
end

function _element_centroid_measure(m::InpMesh, e::Int, dim::Int)
    c = [m.nodes[nid] for nid in _corners(m, e)]
    n = length(c)
    x = sum(c) ./ n
    if dim == 2
        a = 0.0
        for j in 1:n
            k = mod(j, n) + 1
            a += c[j][1] * c[k][2] - c[k][1] * c[j][2]
        end
        return x, 0.5 * abs(a)
    end
    v = if n == 4
        tetvol(c[1], c[2], c[3], c[4])
    elseif n == 8
        tetvol(c[1], c[2], c[4], c[5]) + tetvol(c[2], c[3], c[4], c[7]) +
        tetvol(c[2], c[5], c[6], c[7]) + tetvol(c[4], c[5], c[7], c[8]) +
        tetvol(c[2], c[4], c[5], c[7])
    elseif n == 6
        tetvol(c[1], c[2], c[3], c[6]) + tetvol(c[1], c[2], c[6], c[5]) + tetvol(c[1], c[5], c[6], c[4])
    elseif n == 5
        tetvol(c[1], c[2], c[3], c[5]) + tetvol(c[1], c[3], c[4], c[5])
    else
        error("Only tetrahedral, pyramid, wedge or hexahedral elements are supported (element $e).")
    end
    return x, v
end

# pos (3×N) and measure (area or volume) of the PD points
function _pd_points(m::InpMesh, dom::Vector{Int}, dim::Int)
    pos = zeros(3, length(dom))
    meas = zeros(length(dom))
    for (i, e) in enumerate(dom)
        x, v = _element_centroid_measure(m, e, dim)
        pos[:, i] .= x
        meas[i] = v
        v > 0 || @warn "Element $e has zero $(dim == 2 ? "area" : "volume") (degenerate element)."
    end
    return pos, meas
end

"""
    pd_sets(mesh::InpMesh, geo::Union{Post2D,Post3D}; domain=nothing) -> Dict{Symbol,Vector{Int}}

PD point indices of every named group of the mesh, ready for `point_set!`:
- ELSET of domain elements (e.g. one physical surface/volume per material): the PD
  points of these elements,
- ELSET of lower-dimensional elements, SURFACE, NSET (e.g. a physical line/surface on
  the boundary): the PD points whose element has a face on that group.

`domain` must be the same value that was given to `Post2D`/`Post3D`.

```julia
mesh = read_inp("probe.inp")
geo  = Post2D(mesh; boundary=["Wall"])
sets = pd_sets(mesh, geo)
body = Body(BBTMMaterial(), geo.pos, h .* geo.area)
point_set!(body, :Cu, sets[:Cu])
material!(body, :Cu; ...)
```
"""
function pd_sets(m::InpMesh, geo; domain=nothing)
    dim = geo isa Post2D ? 2 : geo isa Post3D ? 3 : error("pd_sets: geo must be a Post2D or Post3D")
    dom = domain_elements(m, dim; domain=domain)
    length(dom) == size(geo.pos, 2) ||
        error("pd_sets: the mesh/domain does not match geo ($(length(dom)) elements vs $(size(geo.pos, 2)) PD points). " *
              "Pass the same `domain` that was used to build geo.")
    pos, _ = _pd_points(m, dom, dim)
    isapprox(pos, geo.pos; atol=1e-10 * max(1.0, maximum(abs, pos))) ||
        error("pd_sets: PD points of geo do not match this mesh (was geo built from another file or domain?).")
    pd_of = Dict(e => i for (i, e) in enumerate(dom))
    ctx = _face_ctx(m, dom)      # face topology built once for all groups

    sets = Dict{Symbol,Vector{Int}}()
    for name in keys(m.elsets)
        ids = [pd_of[e] for e in m.elsets[name] if haskey(pd_of, e)]
        if isempty(ids)
            faces = _boundary_faces_quiet(m, dom, name, dim, ctx)
            ids = [pd_of[e] for (e, _) in faces]
        end
        isempty(ids) || (sets[Symbol(name)] = unique!(sort!(ids)))
    end
    for dict in (m.surfaces, m.nsets)
        for name in keys(dict)
            haskey(sets, Symbol(name)) && continue
            ids = [pd_of[e] for (e, _) in _boundary_faces_quiet(m, dom, name, dim, ctx)]
            isempty(ids) || (sets[Symbol(name)] = unique!(sort!(ids)))
        end
    end
    return sets
end

_boundary_faces_quiet(m, dom, name, dim, ctx) = _boundary_faces(m, dom, [name], dim; warn=false, ctx=ctx)

"""
    section_materials(mesh::InpMesh) -> Dict{String,String}

Element set => material name, from the `*SOLID SECTION` (etc.) keywords of the file.
"""
section_materials(m::InpMesh) = copy(m.sections)
